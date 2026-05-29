#!/usr/bin/env python3
"""Conformer generation for the ZINC22 3D ligand build.

The default backend is RDKit. The nvmolkit backend is intentionally selected
only through CONFORMER_BACKEND=nvmolkit so the RDKit port remains the baseline.
"""

import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


CONF_ENERGY_PROP = "rdkit_mmff94s_energy"
CONF_FORCEFIELD_PROP = "rdkit_forcefield"
CONF_BACKEND_PROP = "conformer_backend"

DEFAULT_CONF_BUDGET_BASE = 600
DEFAULT_ENERGY_WINDOW = 12.0
DEFAULT_RMSD_THRESHOLD = 0.5
DEFAULT_TIMEOUT = 120
DEFAULT_SEED = 0xF00D
DEFAULT_BACKEND = "rdkit"


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


def _env_bool(name, default=False):
    value = os.environ.get(name, "").strip().lower()
    if not value:
        return default
    return value in {"1", "true", "yes", "on"}


def _env_int_list(name, default):
    value = os.environ.get(name, "").strip()
    if not value:
        return list(default)
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def selected_backend():
    backend = os.environ.get("CONFORMER_BACKEND", DEFAULT_BACKEND).strip().lower()
    if backend in {"", "rdkit"}:
        return "rdkit"
    if backend in {"nv", "nvidia", "nvmolkit"}:
        return "nvmolkit"
    raise ValueError("Unsupported CONFORMER_BACKEND={!r}".format(backend))


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


def _nvmolkit_params(seed):
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    params.pruneRmsThresh = -1.0
    params.useRandomCoords = True
    if hasattr(params, "clearConfs"):
        params.clearConfs = True
    return params


def _nvmolkit_hardware_options():
    try:
        from nvmolkit.types import HardwareOptions
    except ImportError as exc:
        raise ImportError(
            "CONFORMER_BACKEND=nvmolkit requires the nvmolkit package"
        ) from exc

    kwargs = {
        "preprocessingThreads": _env_int("NVMOLKIT_PREPROCESSING_THREADS", -1),
        "batchSize": _env_int("NVMOLKIT_BATCH_SIZE", -1),
        "batchesPerGpu": _env_int("NVMOLKIT_BATCHES_PER_GPU", 4),
        "gpuIds": _env_int_list("NVMOLKIT_GPU_IDS", []),
    }
    try:
        return HardwareOptions(**kwargs)
    except TypeError:
        options = HardwareOptions()
        for key, value in kwargs.items():
            if hasattr(options, key):
                setattr(options, key, value)
        return options


def _optimize_conformers_rdkit(mol):
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


def _normalize_nvmolkit_energies(raw_energies):
    if hasattr(raw_energies, "numpy"):
        raw_energies = raw_energies.numpy()
    if hasattr(raw_energies, "tolist"):
        raw_energies = raw_energies.tolist()
    if not isinstance(raw_energies, (list, tuple)):
        raw_energies = [raw_energies]
    if raw_energies and isinstance(raw_energies[0], (list, tuple)):
        raw_energies = raw_energies[0]
    return [float(energy) for energy in raw_energies]


def _optimize_conformers_nvmolkit(mol):
    if not _env_bool("NVMOLKIT_OPTIMIZE", True):
        return _optimize_conformers_rdkit(mol)
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        return _optimize_conformers_rdkit(mol)
    try:
        from nvmolkit.mmffOptimization import MMFFOptimizeMoleculesConfs
    except ImportError as exc:
        raise ImportError(
            "CONFORMER_BACKEND=nvmolkit requires nvmolkit.mmffOptimization"
        ) from exc

    raw_energies = MMFFOptimizeMoleculesConfs(
        molecules=[mol],
        maxIters=_env_int("NVMOLKIT_MMFF_MAX_ITERS", 1000),
        ignoreInterfragInteractions=False,
        hardwareOptions=_nvmolkit_hardware_options(),
    )
    return _normalize_nvmolkit_energies(raw_energies), "nvMolKit-MMFF94"


def _optimize_conformers(mol, backend):
    if backend == "nvmolkit":
        return _optimize_conformers_nvmolkit(mol)
    return _optimize_conformers_rdkit(mol)


def _copy_conformers(mol, conf_ids, energies, forcefield, backend):
    result = Chem.Mol(mol)
    result.RemoveAllConformers()
    for conf_id in conf_ids:
        conf = Chem.Conformer(mol.GetConformer(conf_id))
        conf.SetId(len(result.GetConformers()))
        conf.SetDoubleProp(CONF_ENERGY_PROP, float(energies[conf_id]))
        conf.SetProp(CONF_FORCEFIELD_PROP, forcefield)
        conf.SetProp(CONF_BACKEND_PROP, backend)
        result.AddConformer(conf, assignId=False)
    return result


def _copy_single_conformer(mol, conf_id, energies, forcefield, backend):
    return _copy_conformers(mol, [conf_id], energies, forcefield, backend)


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


def _embed_conformers_rdkit(mol, num_confs, seed):
    timeout = _env_int("RDKIT_CONF_TIMEOUT", DEFAULT_TIMEOUT)
    params = _etkdg_params(seed, timeout)
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(num_confs), params=params))
    if not cids:
        raise ValueError("RDKit ETKDG generated zero conformers")
    return cids


def _embed_conformers_nvmolkit(mol, num_confs, seed):
    try:
        from nvmolkit.embedMolecules import EmbedMolecules
    except ImportError as exc:
        raise ImportError(
            "CONFORMER_BACKEND=nvmolkit requires the nvmolkit package"
        ) from exc

    mol.RemoveAllConformers()
    EmbedMolecules(
        molecules=[mol],
        params=_nvmolkit_params(seed),
        confsPerMolecule=int(num_confs),
        maxIterations=_env_int("NVMOLKIT_MAX_ITERATIONS", -1),
        hardwareOptions=_nvmolkit_hardware_options(),
    )
    cids = [conf.GetId() for conf in mol.GetConformers()]
    if not cids:
        raise ValueError("nvMolKit ETKDG generated zero conformers")
    return cids


def _embed_conformers(mol, num_confs, seed, backend):
    if backend == "nvmolkit":
        return _embed_conformers_nvmolkit(mol, num_confs, seed)
    return _embed_conformers_rdkit(mol, num_confs, seed)


def generate_seed_conformation(mol, seed=DEFAULT_SEED):
    """Return a one-conformer RDKit Mol suitable for the AMSOL seed mol2."""
    backend = selected_backend()
    seed = _env_int("RDKIT_CONF_SEED", seed)
    mol = Chem.Mol(mol)
    cids = _embed_conformers(mol, 1, seed, backend)
    energies, forcefield = _optimize_conformers(mol, backend)
    energy_by_id = {cid: energies[i] for i, cid in enumerate(cids)}
    return _copy_single_conformer(mol, cids[0], energy_by_id, forcefield, backend)


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
    backend = selected_backend()
    seed = _env_int("RDKIT_CONF_SEED", seed)
    num_confs = _conformer_budget(num_rotatable_terminal_h)
    energy_window = _env_float("RDKIT_CONF_ENERGY_WINDOW", DEFAULT_ENERGY_WINDOW)
    rmsd_threshold = _env_float("RDKIT_CONF_RMSD", DEFAULT_RMSD_THRESHOLD)

    mol = Chem.Mol(mol)
    cids = _embed_conformers(mol, num_confs, seed, backend)
    energies, forcefield = _optimize_conformers(mol, backend)
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
        raise ValueError("conformer filtering removed every conformer")

    seed_mol = _copy_single_conformer(mol, kept[0], energy_by_id, forcefield, backend)
    ensemble = _copy_conformers(mol, kept, energy_by_id, forcefield, backend)
    return seed_mol, ensemble
