import pytest
from rdkit import Chem

from generate.rdkit_conformers import (
    CONF_BACKEND_PROP,
    CONF_ENERGY_PROP,
    CONF_FORCEFIELD_PROP,
    conformer_batch_size,
    generate_conformation_batch,
    generate_conformations,
    generate_seed_conformation,
    selected_backend,
    selected_seed_backend,
)


def mol_from_smiles(smiles):
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def test_refuses_more_than_five_rotatable_terminal_hydrogens():
    with pytest.raises(ValueError, match=">5 rotatable hydrogens"):
        generate_conformations(mol_from_smiles("CCO"), 6)


def test_generates_energy_sorted_pruned_conformers(monkeypatch):
    monkeypatch.setenv("RDKIT_CONF_BUDGET_BASE", "20")

    _seed, ensemble = generate_conformations(mol_from_smiles("CCCCCC"), 0)

    assert 1 <= ensemble.GetNumConformers() <= 20
    energies = [
        ensemble.GetConformer(i).GetDoubleProp(CONF_ENERGY_PROP)
        for i in range(ensemble.GetNumConformers())
    ]
    assert energies == sorted(energies)
    assert ensemble.GetConformer(0).GetProp(CONF_FORCEFIELD_PROP) in {"MMFF94s", "UFF"}
    assert ensemble.GetConformer(0).GetProp(CONF_BACKEND_PROP) == "rdkit"


def test_rejects_unknown_conformer_backend(monkeypatch):
    monkeypatch.setenv("CONFORMER_BACKEND", "not-a-backend")

    with pytest.raises(ValueError, match="Unsupported CONFORMER_BACKEND"):
        selected_backend()


def test_seed_backend_defaults_to_rdkit_when_ensemble_backend_is_nvmolkit(monkeypatch):
    monkeypatch.setenv("CONFORMER_BACKEND", "nvmolkit")

    assert selected_seed_backend() == "rdkit"


def test_seed_backend_uses_seed_override(monkeypatch):
    monkeypatch.setenv("SEED_CONFORMER_BACKEND", "nvmolkit")

    assert selected_seed_backend() == "nvmolkit"


def test_rdkit_batch_api_uses_existing_single_molecule_path(monkeypatch):
    monkeypatch.setenv("RDKIT_CONF_BUDGET_BASE", "5")

    results, failures = generate_conformation_batch([
        ("a", mol_from_smiles("CCCC"), 0),
        ("b", mol_from_smiles("CCO"), 0),
    ])

    assert failures == {}
    assert set(results) == {"a", "b"}
    assert 1 <= results["a"][1].GetNumConformers() <= 5
    assert results["a"][1].GetConformer(0).GetProp(CONF_BACKEND_PROP) == "rdkit"
    assert generate_conformation_batch.last_stats["single_molecule_calls"] == 2


def test_conformer_batch_size_defaults_to_one(monkeypatch):
    monkeypatch.delenv("NVMOLKIT_CONFORMER_BATCH_SIZE", raising=False)

    assert conformer_batch_size() == 1
