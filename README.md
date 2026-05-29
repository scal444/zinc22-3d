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
git switch nvmolkit_port
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

### nvMolKit GPU Environment

On `nvmolkit_port`, use the GPU environment instead:

```bash
mamba env create -f environment-nvmolkit.yml
conda activate zinc22-nvmolkit
```

This spec intentionally pins CUDA PyTorch builds:

```text
pytorch=2.10.0=*cuda129*
libtorch=2.10.0=*cuda129*
triton=3.5.1
cuda-version=12.9
nvmolkit=0.5.0
```

Avoid an unconstrained PyTorch install for this branch. Conda can otherwise
solve to a CPU-only PyTorch build or to a CUDA build that does not support newer
GPUs such as RTX 50-series cards.

Verify the GPU stack before benchmarking:

```bash
nvidia-smi
python - <<'PY'
import torch
import nvmolkit
from nvmolkit.embedMolecules import EmbedMolecules
print("torch", torch.__version__)
print("torch cuda", torch.version.cuda)
print("cuda available", torch.cuda.is_available())
print("device count", torch.cuda.device_count())
if torch.cuda.is_available():
    print("device 0", torch.cuda.get_device_name(0))
print("nvmolkit", getattr(nvmolkit, "__version__", "unknown"))
print("EmbedMolecules import ok")
PY
```

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

For nvMolKit conformer generation:

```bash
export CONFORMER_BACKEND=nvmolkit
export SEED_CONFORMER_BACKEND=${SEED_CONFORMER_BACKEND:-rdkit}
export NVMOLKIT_CONFORMER_BATCH_SIZE=${NVMOLKIT_CONFORMER_BATCH_SIZE:-1}
```

The AMSOL seed conformer defaults to RDKit even when the ensemble conformer
backend is nvMolKit. That avoids GPU launch overhead for the one-conformer seed
task.

`NVMOLKIT_CONFORMER_BATCH_SIZE` is experimental. Values above 1 group molecules
with the same conformer budget into one nvMolKit call inside each worker
process. Assess it with longer benchmarks because short runs are sensitive to
GPU warmup, MPS state, and tranche composition.

## CUDA MPS for Parallel nvMolKit Runs

If several worker processes share one GPU, enable CUDA MPS before launching the
workers. Without MPS, multi-process nvMolKit runs can serialize enough GPU work
that the speedup is much smaller.

```bash
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$USER
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-mps-$USER-log
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
nvidia-cuda-mps-control -d
```

Run the benchmark or production workers with those two environment variables
still set. Stop the temporary MPS daemon after the run:

```bash
printf 'quit\n' | nvidia-cuda-mps-control
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
