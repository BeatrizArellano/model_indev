#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT

(
    cd "$work_dir"
    "${FC:-gfortran}" -std=f2008 -Wall -Wextra "$script_dir/test_restart_csv_format.F90" -o restart_format_test
    ./restart_format_test
    python3 - <<'PY'
import csv
import math

with open('restart_format_test.csv', newline='') as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 6
for row in rows:
    for field in ('depth', 'value'):
        assert 'E' in row[field].upper(), (field, row[field])
        assert math.isfinite(float(row[field])), row
assert float(rows[1]['value']) > 0.0  # minimum positive normal
assert float(rows[2]['value']) > 0.0  # minimum positive subnormal
assert float(rows[3]['value']) == 0.0
assert float(rows[5]['value']) < 0.0
print('Python CSV parsing checks passed')
PY
)
