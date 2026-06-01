# Compatibility wrapper for the license-free RDKit/nvMolKit branches.
# Source this after activating the conda environment.

_ENV_SH_FILE="${BASH_SOURCE:-$0}"
_ENV_SH_DIR="$(cd "$(dirname "$_ENV_SH_FILE")" && pwd -P)"
. "$_ENV_SH_DIR/zinc22-env.sh"

export BINDIR="$DOCKBASE/ligand/submit"
export SOFT_HOME="${SOFT_HOME:-$DOCKBASE/third_party}"
export SUBMIT_MODE="${SUBMIT_MODE:-TEST_LOCAL}"
export BUILD_MOL2="${BUILD_MOL2:-false}"
export SKIP_DELETE="${SKIP_DELETE:-false}"
