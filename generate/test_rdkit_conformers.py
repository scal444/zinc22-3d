import pytest
from rdkit import Chem

from generate.rdkit_conformers import (
    CONF_ENERGY_PROP,
    CONF_FORCEFIELD_PROP,
    generate_conformations,
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
