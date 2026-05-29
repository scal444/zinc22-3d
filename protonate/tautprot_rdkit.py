#!/usr/bin/env python3
"""RDKit + Dimorphite-DL replacement for tautprot_cxcalc.sh.

The legacy ChemAxon wrapper emitted:

    SMILES name tautomer-dist protomer-dist score

RDKit and Dimorphite-DL do not produce population percentages, so this script
emits deterministic ranking scores in those columns. The downstream coalescer
only uses the last extra column for sorting before assigning protomer ids.
"""

import argparse
import os
import sys

from dimorphite_dl import protonate_smiles
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


DEFAULT_TAUT_PROT_CUTOFF = float(os.environ.get("TAUT_PROT_CUTOFF", 10))
DEFAULT_TAUTOMER_LIMIT = float(os.environ.get("TAUTOMER_LIMIT", 20))
DEFAULT_PROTOMER_LIMIT = float(os.environ.get("PROTOMER_LIMIT", 20))
DEFAULT_MAX_TAUTOMERS = int(os.environ.get("RDKIT_MAX_TAUTOMERS", 1))
DEFAULT_MAX_PROTOMERS = int(os.environ.get("RDKIT_MAX_PROTOMERS", 1))
DEFAULT_TAUTOMER_MODE = os.environ.get("RDKIT_TAUTOMER_MODE", "canonical")
DEFAULT_DIMORPHITE_PRECISION = float(os.environ.get("DIMORPHITE_PRECISION", 1.0))


def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def pseudo_percent(rank):
    """Convert a zero-based rank into a ChemAxon-like sortable percentage."""
    return max(1.0, 100.0 - (rank * 10.0))


def ranked_tautomers(mol, mode, max_tautomers):
    enumerator = rdMolStandardize.TautomerEnumerator()
    if mode == "canonical":
        tautomer = enumerator.Canonicalize(mol)
        return [(Chem.MolToSmiles(tautomer, isomericSmiles=True), 100.0)]

    canonical = Chem.MolToSmiles(enumerator.Canonicalize(mol), isomericSmiles=True)
    seen = set()
    tautomers = []
    for tautomer in enumerator.Enumerate(mol):
        smiles = Chem.MolToSmiles(tautomer, isomericSmiles=True)
        if smiles in seen:
            continue
        seen.add(smiles)
        score = enumerator.ScoreTautomer(tautomer)
        tautomers.append((smiles, score))

    tautomers.sort(key=lambda item: (item[0] != canonical, -item[1], item[0]))
    ranked = []
    for rank, (smiles, _score) in enumerate(tautomers[:max_tautomers]):
        ranked.append((smiles, pseudo_percent(rank)))
    return ranked


def ranked_protomers(smiles, ph, precision, max_protomers):
    try:
        variants = protonate_smiles(
            smiles,
            ph_min=ph,
            ph_max=ph,
            precision=precision,
            max_variants=max_protomers,
            validate_output=True,
        )
    except Exception:
        variants = []

    seen = set()
    protomers = []
    for variant in variants:
        tokens = variant.split()
        canonical = canonical_smiles(tokens[0]) if tokens else None
        if canonical is None or canonical in seen:
            continue
        seen.add(canonical)
        protomers.append(canonical)

    fallback = canonical_smiles(smiles)
    if fallback is not None and fallback not in seen:
        protomers.append(fallback)

    protomers.sort()
    return [
        (protomer, pseudo_percent(rank))
        for rank, protomer in enumerate(protomers[:max_protomers])
    ]


def generate_records(
    smiles,
    name,
    ph,
    tautomer_mode,
    max_tautomers,
    max_protomers,
    tautomer_limit,
    protomer_limit,
    score_cutoff,
    precision,
):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")

    seen = set()
    records = []
    for tautomer_smiles, tautomer_score in ranked_tautomers(mol, tautomer_mode, max_tautomers):
        if tautomer_score < tautomer_limit:
            continue
        for protomer_smiles, protomer_score in ranked_protomers(
            tautomer_smiles,
            ph,
            precision,
            max_protomers,
        ):
            if protomer_score < protomer_limit:
                continue
            score = (tautomer_score * protomer_score) / 100.0
            if score < score_cutoff or protomer_smiles in seen:
                continue
            seen.add(protomer_smiles)
            records.append((protomer_smiles, name, tautomer_score, protomer_score, score))

    records.sort(key=lambda item: (-item[-1], -item[2], -item[3], item[0]))
    return records


def input_rows(handle):
    for idx, line in enumerate(handle, start=1):
        line = line.strip()
        if not line:
            continue
        tokens = line.split()
        smiles = tokens[0]
        name = tokens[1] if len(tokens) > 1 else "mol{}".format(idx)
        yield smiles, name


def main(argv=None, stdout=sys.stdout, stderr=sys.stderr, stdin=None):
    if stdin is None:
        stdin = sys.stdin
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--help", action="help", help="show this help message and exit")
    parser.add_argument("input", nargs="?", default="-")
    parser.add_argument("-H", "--pH", dest="ph", type=float, default=7.4)
    parser.add_argument("-c", "--cutoff", type=float, default=DEFAULT_TAUT_PROT_CUTOFF)
    parser.add_argument("-t", "--tautomer-limit", type=float, default=DEFAULT_TAUTOMER_LIMIT)
    parser.add_argument("-p", "--protomer-limit", type=float, default=DEFAULT_PROTOMER_LIMIT)
    parser.add_argument("-h", "--header", action="store_true")
    parser.add_argument("--tautomer-mode", choices=("enumerate", "canonical"), default=DEFAULT_TAUTOMER_MODE)
    parser.add_argument("--max-tautomers", type=int, default=DEFAULT_MAX_TAUTOMERS)
    parser.add_argument("--max-protomers", type=int, default=DEFAULT_MAX_PROTOMERS)
    parser.add_argument("--dimorphite-precision", type=float, default=DEFAULT_DIMORPHITE_PRECISION)
    args = parser.parse_args(argv)

    handle = stdin if args.input == "-" else open(args.input)
    close_handle = handle is not stdin
    try:
        if args.header:
            print("smiles name tautomer-dist protomer-dist score", file=stdout)
        for smiles, name in input_rows(handle):
            try:
                records = generate_records(
                    smiles=smiles,
                    name=name,
                    ph=args.ph,
                    tautomer_mode=args.tautomer_mode,
                    max_tautomers=args.max_tautomers,
                    max_protomers=args.max_protomers,
                    tautomer_limit=args.tautomer_limit,
                    protomer_limit=args.protomer_limit,
                    score_cutoff=args.cutoff,
                    precision=args.dimorphite_precision,
                )
            except Exception as exc:
                print("tautprot_rdkit skipped {}: {}".format(name, exc), file=stderr)
                continue
            for record in records:
                print("{} {} {:.1f} {:.1f} {:.1f}".format(*record), file=stdout)
    finally:
        if close_handle:
            handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
