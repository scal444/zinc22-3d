#!/usr/bin/env python3
"""RDKit conformer generation for the ZINC22 3D ligand build."""

import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


CONF_ENERGY_PROP = "rdkit_mmff94s_energy"
CONF_FORCEFIELD_PROP = "rdkit_forcefield"

DEFAULT_CONF_BUDGET_BASE = 600
DEFAULT_ENERGY_WINDOW = 12.0
DEFAULT_RMSD_THRESHOLD = 0.5
DEFAULT_TIMEOUT = 120
DEFAULT_SEED = 0xF00D


def _env_int(name, default):
    value = os.environ.get(name, "").strip()
    if not value:
        return default
    return int(value, 0)


def _env_float(name, default):
    value = os.environ.get(name, "").strip()
    if not value:
        return default
    return float(value)


def _conformer_budget(num_rotatable_terminal_h):
    num_confs = _env_int("RDKIT_CONF_BUDGET_BASE", DEFAULT_CONF_BUDGET_BASE)
    if num_rotatable_terminal_h > 5:
        raise ValueError("Refusing to build conformations with >5 rotatable hydrogens")
    if num_rotatable_terminal_h >= 4:
        return num_confs // 30
    if num_rotatable_terminal_h >= 2:
        return num_confs // 3
    return num_confs


def _etkdg_params(seed, timeout):
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    params.pruneRmsThresh = -1.0
    params.useRandomCoords = False
    if hasattr(params, "timeout"):
        params.timeout = int(timeout)
    if hasattr(params, "clearConfs"):
        params.clearConfs = True
    return params


def _optimize_conformers(mol):
    try:
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ValueError("MMFF94s parameters unavailable")
        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            mmffVariant="MMFF94s",
            maxIters=1000,
            ignoreInterfragInteractions=False,
        )
        forcefield = "MMFF94s"
    except Exception:
        # OMEGA used MMFF94Smod. RDKit MMFF94s is the closest available analog,
        # but molecules outside MMFF typing fall back to UFF rather than being
        # dropped from the build.
        results = AllChem.UFFOptimizeMoleculeConfs(
            mol,
            maxIters=1000,
            ignoreInterfragInteractions=False,
        )
        forcefield = "UFF"
    energies = [energy for _not_converged, energy in results]
    return energies, forcefield


def _copy_conformers(mol, conf_ids, energies, forcefield):
    result = Chem.Mol(mol)
    result.RemoveAllConformers()
    for conf_id in conf_ids:
        conf = Chem.Conformer(mol.GetConformer(conf_id))
        conf.SetId(len(result.GetConformers()))
        conf.SetDoubleProp(CONF_ENERGY_PROP, float(energies[conf_id]))
        conf.SetProp(CONF_FORCEFIELD_PROP, forcefield)
        result.AddConformer(conf, assignId=False)
    return result


def _copy_single_conformer(mol, conf_id, energies, forcefield):
    return _copy_conformers(mol, [conf_id], energies, forcefield)


def _rmsd_by_atom_order(mol, prb_id, ref_id):
    atom_map = [(idx, idx) for idx in range(mol.GetNumAtoms())]
    rmsd, _transform = rdMolAlign.GetAlignmentTransform(
        mol,
        mol,
        prbCid=prb_id,
        refCid=ref_id,
        atomMap=atom_map,
    )
    return rmsd


def _embed_conformers(mol, num_confs, seed):
    timeout = _env_int("RDKIT_CONF_TIMEOUT", DEFAULT_TIMEOUT)
    params = _etkdg_params(seed, timeout)
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(num_confs), params=params))
    if not cids:
        raise ValueError("RDKit ETKDG generated zero conformers")
    return cids


def generate_seed_conformation(mol, seed=DEFAULT_SEED):
    """Return a one-conformer RDKit Mol suitable for the AMSOL seed mol2."""
    seed = _env_int("RDKIT_CONF_SEED", seed)
    mol = Chem.Mol(mol)
    cids = _embed_conformers(mol, 1, seed)
    energies, forcefield = _optimize_conformers(mol)
    energy_by_id = {cid: energies[i] for i, cid in enumerate(cids)}
    return _copy_single_conformer(mol, cids[0], energy_by_id, forcefield)


def generate_conformations(mol, num_rotatable_terminal_h, seed=DEFAULT_SEED):
    """Return (seed_mol_1conf, ensemble_mol_Nconfs).

    mol: an RDKit Mol (sanitized, 2D or no coords), explicit Hs added.
    num_rotatable_terminal_h: int from hydrogens.count_hydrogens().

    The ensemble matches OMEGA's production budget gates: ETKDGv3 embeds up to
    the requested budget, MMFF94s minimizes, a 12 kcal/mol window filters, then
    energy-sorted greedy all-atom GetBestRMS pruning keeps conformers >=0.5 A
    from previously accepted conformers. MMFF94s is not OMEGA's MMFF94Smod, so
    absolute and relative energies are expected to differ from legacy builds.
    """
    seed = _env_int("RDKIT_CONF_SEED", seed)
    num_confs = _conformer_budget(num_rotatable_terminal_h)
    energy_window = _env_float("RDKIT_CONF_ENERGY_WINDOW", DEFAULT_ENERGY_WINDOW)
    rmsd_threshold = _env_float("RDKIT_CONF_RMSD", DEFAULT_RMSD_THRESHOLD)

    mol = Chem.Mol(mol)
    cids = _embed_conformers(mol, num_confs, seed)
    energies, forcefield = _optimize_conformers(mol)
    energy_by_id = {cid: energies[i] for i, cid in enumerate(cids)}

    min_energy = min(energy_by_id.values())
    energy_sorted_cids = sorted(
        (cid for cid in cids if energy_by_id[cid] - min_energy <= energy_window),
        key=lambda cid: energy_by_id[cid],
    )

    kept = []
    for cid in energy_sorted_cids:
        if all(_rmsd_by_atom_order(mol, cid, k) >= rmsd_threshold for k in kept):
            kept.append(cid)
        if len(kept) >= num_confs:
            break

    if not kept:
        raise ValueError("RDKit conformer filtering removed every conformer")

    seed_mol = _copy_single_conformer(mol, kept[0], energy_by_id, forcefield)
    ensemble = _copy_conformers(mol, kept, energy_by_id, forcefield)
    return seed_mol, ensemble
