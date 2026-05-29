# ZINC22 3D Building Pipeline

This branch builds ZINC-style 3D ligand outputs from SMILES using a
license-free toolchain:

- RDKit for seed structures, conformer generation, and tautomer handling.
- Dimorphite-DL for protonation.
- OpenBabel for file-format conversion.
- AMSOL 7.1 for solvation.

OpenEye OMEGA, ChemAxon/JChem, and CORINA are not required on this branch.

## Repository Layout

Some legacy scripts expect the repository to be available as
`$DOCKBASE/ligand`. Either clone the repository directly at that path or create
a symlink:

```bash
export DOCKBASE=$HOME/experiments/zinc
git clone git@github.com:scal444/zinc22-3d.git "$DOCKBASE/zinc22-3d"
ln -sfn "$DOCKBASE/zinc22-3d" "$DOCKBASE/ligand"
cd "$DOCKBASE/zinc22-3d"
git switch rdkit_port
```

## Conda Environment

Use the compact RDKit environment for new installs:

```bash
mamba env create -f environment-rdkit.yml
conda activate zinc22-rdkit
```

If `mamba` is unavailable, `conda env create -f environment-rdkit.yml` also
works, just more slowly.

The older `environment.yml` is retained for compatibility with historical
builds; new users should start from `environment-rdkit.yml`.

## AMSOL

AMSOL is not a Python package and is not vendored in this repository. Download
AMSOL 7.1 from:

```text
https://comp.chem.umn.edu/amsol/
```

Build it on the target machine and point `AMSOLEXE` at the executable. One
working layout is:

```bash
mkdir -p "$DOCKBASE/third_party/amsol/amsol7.1"
# place or build the executable here:
export AMSOLEXE="$DOCKBASE/third_party/amsol/amsol7.1/amsol7.1.exe"
chmod +x "$AMSOLEXE"
```

## Runtime Environment

After activating the conda environment, set:

```bash
export DOCKBASE=${DOCKBASE:-$HOME/experiments/zinc}
export AMSOLEXE="$DOCKBASE/third_party/amsol/amsol7.1/amsol7.1.exe"
export OBABELBASE="$CONDA_PREFIX"
export CSH="$CONDA_PREFIX/bin/tcsh"
export PYTHONPATH="$DOCKBASE/ligand/mol2db2_py3_strain:$DOCKBASE/ligand/strain"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
```

Conformer defaults:

```bash
export RDKIT_CONF_BUDGET_BASE=${RDKIT_CONF_BUDGET_BASE:-600}
export RDKIT_CONF_ENERGY_WINDOW=${RDKIT_CONF_ENERGY_WINDOW:-12}
export RDKIT_CONF_RMSD=${RDKIT_CONF_RMSD:-0.5}
export RDKIT_CONF_TIMEOUT=${RDKIT_CONF_TIMEOUT:-120}
export RDKIT_CONF_SEED=${RDKIT_CONF_SEED:-0xf00d}
```

## Smoke Test

```bash
python -m pytest -q
bash -n generate/build_database_ligand_strain_noH_btingle.sh

head -n 1 validation/zinc22_random_10.smi > /tmp/zinc22-one.smi
bash generate/build_database_ligand_strain_noH_btingle.sh \
  -H 7.4 --no-db \
  -d /tmp/zinc22-smoke \
  /tmp/zinc22-one.smi
```

A successful smoke test should report `outputs built: 1` and create
`/tmp/zinc22-smoke/working/output.tar.gz`.

## Production Submission

For batch submission, provide the variables checked by
`submit/submit-all.bash`, especially:

```bash
export INPUT_FILE=/path/to/input.smi
export OUTPUT_DEST=/path/to/output-root
export SOFT_HOME=/path/to/software-root
```

Then run:

```bash
bash submit/submit-all.bash
```
