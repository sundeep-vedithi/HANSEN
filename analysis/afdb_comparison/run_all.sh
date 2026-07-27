#!/bin/bash
# Full pipeline, in dependency order. 03_analyse.py rewrites summary.json, so the
# two scripts that append to it (03b, 05) must run after it and before 04.
set -euo pipefail
cd "$(dirname "$0")"

# Machine-specific input locations, if present (not in the public release).
[ -f ./local_paths.sh ] && . ./local_paths.sh
PY=python3
VPY=./.venv/bin/python     # venv carries matplotlib / openpyxl / python-docx / tmtools

$PY  01_build_manifest.py
$PY  02_compare.py
$PY  02b_crossmethod.py
$VPY 03_analyse.py
$PY  03b_worked_example.py
$VPY 03c_profiles.py
$VPY 05_validate.py
$VPY 04_report.py
echo "done"
