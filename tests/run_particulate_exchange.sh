#!/usr/bin/env bash
# Run the particulate-exchange unit tests against the actual repository modules.
# Usage: bash tests/run_particulate_exchange.sh
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT

fc="${FC:-gfortran}"
cd "$build_dir"

"$fc" -std=f2008 -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    -c "$repo_root/src/utils/precision_types.F90"
"$fc" -std=f2008 -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    -c "$repo_root/src/bio/sediment_exchange.F90"
"$fc" -std=f2008 -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    precision_types.o sediment_exchange.o \
    "$repo_root/tests/test_particulate_exchange.F90" -o test_particulate_exchange

./test_particulate_exchange
