# Source this after activating the conda environment.
#
# Example:
#   conda activate zinc22-nvmolkit
#   source ./zinc22-env.sh

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "CONDA_PREFIX is not set. Activate the conda environment first." >&2
    return 1 2>/dev/null || exit 1
fi

_ZINC22_ENV_FILE="${BASH_SOURCE:-$0}"
_ZINC22_REPO="$(cd "$(dirname "$_ZINC22_ENV_FILE")" && pwd -P)"

export DOCKBASE="${DOCKBASE:-$(dirname "$_ZINC22_REPO")}"

if [ ! -e "$DOCKBASE/ligand" ]; then
    ln -sfn "$_ZINC22_REPO" "$DOCKBASE/ligand" 2>/dev/null || true
fi

export AMSOLEXE="${AMSOLEXE:-$DOCKBASE/third_party/amsol/amsol7.1/amsol7.1.exe}"
export OBABELBASE="$CONDA_PREFIX"
export CSH="$CONDA_PREFIX/bin/tcsh"
export SHELL="${SHELL:-/bin/bash}"
export PYTHONPATH="$DOCKBASE/ligand/mol2db2_py3_strain:$DOCKBASE/ligand/strain${PYTHONPATH:+:$PYTHONPATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$CONDA_PREFIX/bin:$PATH"

export CONFORMER_BACKEND="${CONFORMER_BACKEND:-rdkit}"
export SEED_CONFORMER_BACKEND="${SEED_CONFORMER_BACKEND:-rdkit}"
export RDKIT_CONF_BUDGET_BASE="${RDKIT_CONF_BUDGET_BASE:-600}"
export RDKIT_CONF_ENERGY_WINDOW="${RDKIT_CONF_ENERGY_WINDOW:-12}"
export RDKIT_CONF_RMSD="${RDKIT_CONF_RMSD:-0.5}"
export RDKIT_CONF_TIMEOUT="${RDKIT_CONF_TIMEOUT:-120}"
export RDKIT_CONF_SEED="${RDKIT_CONF_SEED:-0xf00d}"
