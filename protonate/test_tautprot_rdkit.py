import io

from protonate.tautprot_rdkit import generate_records, main


def test_generate_records_returns_legacy_columns():
    records = generate_records(
        smiles="c1ccccc1O",
        name="phenol",
        ph=7.4,
        tautomer_mode="canonical",
        max_tautomers=5,
        max_protomers=5,
        tautomer_limit=20,
        protomer_limit=20,
        score_cutoff=10,
        precision=1.0,
    )

    assert records
    assert len(records[0]) == 5
    assert records[0][1] == "phenol"
    assert records[0][-1] >= 10


def test_main_writes_chem_axon_compatible_rows():
    stdout = io.StringIO()
    stdin = io.StringIO("CC(=O)O acetate\n")

    rc = main(
        [
            "-H",
            "7.4",
            "--tautomer-mode",
            "canonical",
            "--max-protomers",
            "4",
            "-",
        ],
        stdout=stdout,
        stderr=io.StringIO(),
        stdin=stdin,
    )

    assert rc == 0
    fields = stdout.getvalue().splitlines()[0].split()
    assert len(fields) == 5
    assert fields[1] == "acetate"
