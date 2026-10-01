#!/usr/bin/env bash
# Regenerate every number, table and figure in the manuscript, then validate.
#
# This is the command sequence Appendix F of the paper refers to.  It consumes
# the adopted chain read-only -- repair_v8 since 2026-10-01 (issue #19; set in
# scripts/chain.py, override with CYGOB2_CHAIN) -- and rewrites only generated
# products:
# manuscript/numbers.tex, manuscript/tables_generated.tex, the paper figures and
# the WP12 tables.  It does NOT re-run WP1-WP9; those are the frozen chain.
#
# Usage, from the repository root:
#   bash scripts/run_manuscript_chain.sh
#
# Requires the `cygob2-gaia` environment
# (provenance/environment_cygob2-gaia.yml, or the smaller
# environment_cygob2-gaia_from-history.yml) and, for the final step only, a
# LaTeX toolchain providing latexmk or tectonic.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts"
PYTHON_BIN="${PYTHON_BIN:-python}"

step() { printf '\n=== %s\n' "$1"; }

step "0  WP12 preregistration (idempotent: refuses to overwrite an existing one)"
if [ -f provenance/wp12_revision_prereg.json ]; then
  echo "    already written -- a preregistration is written once, not re-run"
else
  "$PYTHON_BIN" scripts/wp12_prereg.py
fi

step "1  authorize and audit the manuscript inputs"
"$PYTHON_BIN" scripts/wp10_inputs.py

step "2  WP12.1  residual-gate landscape and the Appendix A branch table"
"$PYTHON_BIN" scripts/wp12_gate_landscape.py

step "3  WP12.2  subgroup closure, closing slopes, mixed-slope ledger"
"$PYTHON_BIN" scripts/wp12_closure_slopes.py

step "4  WP12.3/12.4  C4 scan, C3 subtype mappings, scenario score"
"$PYTHON_BIN" scripts/wp12_scenario_score.py

step "5  WP12.5  neighbouring-association budget (gate R6 fails by design)"
"$PYTHON_BIN" scripts/wp12_neighbour_budget.py

step "6  generated LaTeX tables"
"$PYTHON_BIN" scripts/wp12_tables.py

step "7  generated numbers -> manuscript/numbers.tex"
"$PYTHON_BIN" scripts/wp10_numbers.py

step "8  paper figures"
"$PYTHON_BIN" scripts/wp12_figures.py

step "9  the WP10 validation gate (all seven checks must pass)"
"$PYTHON_BIN" scripts/wp10_validate.py

step "10 file inventory with checksums"
"$PYTHON_BIN" audit.py

step "11 compile"
if command -v latexmk >/dev/null 2>&1; then
  ( cd manuscript && latexmk -pdf main.tex )
elif command -v tectonic >/dev/null 2>&1; then
  ( cd manuscript && tectonic --keep-logs --keep-intermediates main.tex )
  echo "    compiled with tectonic (latexmk was not available)"
else
  echo "    no supported LaTeX engine found (need latexmk or tectonic)" >&2
  exit 1
fi

if [ -f manuscript/main.pdf ]; then
  echo
  echo "PDF SHA-256:"
  shasum -a 256 manuscript/main.pdf
else
  echo "    compilation returned without producing manuscript/main.pdf" >&2
  exit 1
fi
