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
git switch nvmolkit_profiling
```

## Conda Environment

Use the compact RDKit environment for new installs:

```bash
mamba env create -f environment-rdkit.yml
conda activate zinc22-rdkit
source ./zinc22-env.sh
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
source ./zinc22-env.sh
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
and build it from the University of Minnesota distribution site. The AMSOL
compile script expects an old `g77` command, so create a local shim that calls
the conda Fortran compiler with legacy flags.

The conda environment files include `gfortran_linux-64`. If the environment was
created before that dependency was added, install it into the active env:

```bash
mamba install -c conda-forge gfortran_linux-64
```

Then download and build AMSOL:

```bash
export DOCKBASE=${DOCKBASE:-$HOME/experiments/zinc}
mkdir -p "$DOCKBASE/third_party/amsol/bin"
cd "$DOCKBASE/third_party/amsol"

curl -L -o amsol7.1.tar.xz \
  https://comp.chem.umn.edu/sds/amsol/amsol7.1.tar.xz
tar -xf amsol7.1.tar.xz

GFORTRAN="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gfortran"
test -x "$GFORTRAN"

cat > "$DOCKBASE/third_party/amsol/bin/g77" <<EOF
#!/bin/sh
exec "$GFORTRAN" -std=legacy -fdec -fallow-argument-mismatch "\$@"
EOF
chmod +x "$DOCKBASE/third_party/amsol/bin/g77"

cd "$DOCKBASE/third_party/amsol/amsol7.1"
# Modern gfortran rejects AMSOL's old NAME= OPEN specifier. FILE= is the
# equivalent standard spelling.
sed -i "s/OPEN(19,NAME='fort.19')/OPEN(19,FILE='fort.19')/" new/amsol.f
sed -i "s/OPEN(20,NAME='fort.20')/OPEN(20,FILE='fort.20')/" new/amsol.f

printf 'man\nlinux\namsol7.1.exe\nsn\n' \
  | env PATH="$DOCKBASE/third_party/amsol/bin:$PATH" tcsh -f ./amsol.compile

export AMSOLEXE="$DOCKBASE/third_party/amsol/amsol7.1/amsol7.1.exe"
test -x "$AMSOLEXE"
```

The build emits many legacy Fortran warnings. The important success marker is
`AMSOL Compiled Successfully` followed by a present executable at `$AMSOLEXE`.

Do not commit the downloaded AMSOL source or binary to this repository.

## Runtime Environment

After activating the conda environment, source the runtime helper:

```bash
source ./zinc22-env.sh
```

It sets the paths expected by the legacy shell/Python pipeline:

```bash
DOCKBASE
AMSOLEXE
OBABELBASE
CSH
SHELL
PYTHONPATH
LD_LIBRARY_PATH
```

It also sets the conformer defaults, which can be overridden before running a
benchmark:

```bash
CONFORMER_BACKEND=rdkit
SEED_CONFORMER_BACKEND=rdkit
RDKIT_CONF_BUDGET_BASE=600
RDKIT_CONF_ENERGY_WINDOW=12
RDKIT_CONF_RMSD=0.5
RDKIT_CONF_TIMEOUT=120
RDKIT_CONF_SEED=0xf00d
```

For nvMolKit conformer generation, set:

```bash
export CONFORMER_BACKEND=nvmolkit
export SEED_CONFORMER_BACKEND=rdkit
```

The AMSOL seed conformer defaults to RDKit even when the ensemble conformer
backend is nvMolKit. That avoids GPU launch overhead for the one-conformer seed
task.

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

## NVTX Profiling

Set `ZINC_NVTX=1` to emit optional NVTX ranges around the seed, solvation,
conformer, strain, DB2, conversion, and archive-write stages. With the variable
unset, the annotations are no-ops.

Example single-process profile:

```bash
source ./zinc22-env.sh
export ZINC_NVTX=1
nsys profile --trace=cuda,nvtx,osrt -o /tmp/zinc22-profile \
  bash generate/build_database_ligand_strain_noH_btingle.sh \
    -H 7.4 --no-db \
    -d /tmp/zinc22-profile-run \
    /path/to/input.smi
```

## Benchmarks

Benchmark inputs are plain whitespace-delimited SMILES files:

```text
SMILES ZINC_ID
```

The committed `validation/zinc22_random_10.smi` file is useful for smoke tests,
short profiling runs, and checking that the full path still produces output. It
is not large enough for stable throughput estimates.

For throughput numbers, use a single ZINC22 tranche rather than a chemically
mixed random set. Runtime and final DB2 set counts depend strongly on heavy atom
count, logP bin, charge, and terminal-H expansion. The reference run in
`validation/budget600_500_perf_summary.md` used 500 neutral molecules from
public ZINC22 `H24P100-N-oaa`, selected from IDs with official DB2 headers and
matching tranche SMILES.

The public ZINC file server URLs used for that reference run have returned HTTP
403 in later checks from this machine. Do not use a nonexistent
`validation/bench_inputs/...` path as a benchmark default. Use the committed
validation input for a runnable pipeline check, or point `INPUT` at a local
tranche file in `SMILES ZINC_ID` format when you have one.

If a local tranche file is already in `SMILES ZINC_ID` format, sample it:

```bash
mkdir -p validation/bench_inputs
grep -v '^\s*$' /path/to/local-tranche.smi \
  | shuf -n 500 \
  > validation/bench_inputs/local_tranche_500.smi
```

If the file has extra columns, keep only the first two:

```bash
gzip -dc /path/to/local-tranche.smi.gz \
  | awk '{print $1, $2}' \
  | shuf -n 500 \
  > validation/bench_inputs/local_tranche_500.smi
```

### Single-Worker Benchmark

Use this for quick RDKit versus nvMolKit checks and single-process Nsight
profiles:

```bash
source ./zinc22-env.sh

INPUT=${INPUT:-validation/zinc22_random_10.smi}
STAMP=$(date +%Y%m%d-%H%M%S)
WORK=/tmp/zinc22-bench-single-$STAMP

test -s "$INPUT"
test -x "$AMSOLEXE"

export CONFORMER_BACKEND=${CONFORMER_BACKEND:-rdkit}
export SEED_CONFORMER_BACKEND=rdkit
export RDKIT_CONF_BUDGET_BASE=20
export RDKIT_CONF_TIMEOUT=20

/usr/bin/time -v bash generate/build_database_ligand_strain_noH_btingle.sh \
  -H 7.4 --no-db \
  -d "$WORK" \
  "$INPUT" \
  2>&1 | tee "$WORK.log"
```

The log should include `outputs built`, failure counts, and an `elapsed times`
block with seed, solvation, conformer, strain, DB2, conversion, and total times.
If the log prints `found`, the work directory already had an `output.tar.gz` and
the run is exercising restart behavior. Delete the work directory or choose a
fresh `WORK` before timing a clean build.

For the ZINC-equivalent conformer budget used in the validation summary, run the
same command with:

```bash
export RDKIT_CONF_BUDGET_BASE=600
export RDKIT_CONF_TIMEOUT=120
```

### Four-Worker nvMolKit/MPS Benchmark

Use this for GPU utilization experiments. It splits one input file into four
line-balanced shards and runs four independent full-pipeline workers. Use a
larger tranche sample for throughput numbers; the 10-molecule validation input
is only a short profiling workload.

```bash
source ./zinc22-env.sh

INPUT=${INPUT:-validation/zinc22_random_10.smi}
STAMP=$(date +%Y%m%d-%H%M%S)
WORK=/tmp/zinc22-bench-mps4-$STAMP

test -s "$INPUT"
test -x "$AMSOLEXE"

export CONFORMER_BACKEND=nvmolkit
export SEED_CONFORMER_BACKEND=rdkit
export RDKIT_CONF_BUDGET_BASE=${RDKIT_CONF_BUDGET_BASE:-600}
export RDKIT_CONF_TIMEOUT=${RDKIT_CONF_TIMEOUT:-120}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NVMOLKIT_PREPROCESSING_THREADS=1

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$USER-$STAMP
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-mps-$USER-$STAMP-log
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
nvidia-cuda-mps-control -d

mkdir -p "$WORK/splits" "$WORK/logs" "$WORK/jobs"
split -n l/4 -d -a 2 "$INPUT" "$WORK/splits/input_"

pids=()
for shard in "$WORK"/splits/input_*; do
  name=$(basename "$shard")
  bash generate/build_database_ligand_strain_noH_btingle.sh \
    -H 7.4 --no-db \
    -d "$WORK/jobs/$name" \
    "$shard" \
    > "$WORK/logs/$name.log" 2>&1 &
  pids+=($!)
done

failed=0
for pid in "${pids[@]}"; do
  wait "$pid" || failed=1
done

printf 'quit\n' | nvidia-cuda-mps-control
grep -H "outputs built:\|conformer failures:\|strain failures:\|db2 failures:\|archive conversion failures:\|confs:" "$WORK"/logs/*.log
if [ "$failed" -ne 0 ]; then
  echo "one or more workers failed" >&2
fi
```

For an Nsight Systems version of this run, wrap the four-worker shell block in
`nsys profile --trace=cuda,nvtx,osrt --trace-fork-before-exec=true --wait=all`.
Some Nsight versions require `sudo` for GPU metrics. If using `sudo`, pass the
conda/runtime environment explicitly and delete the selected `WORK` path before
each profile.

## Smoke Test

```bash
source ./zinc22-env.sh

export RDKIT_CONF_BUDGET_BASE=20
export RDKIT_CONF_TIMEOUT=20

python -m pytest -q
bash -n generate/build_database_ligand_strain_noH_btingle.sh

test -x "$AMSOLEXE"
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
