#!/usr/bin/env python3
"""Re-base the WP13 brief's disclosed prior knowledge (section 2) on repair_v9.

Same computation as scripts/wp13_brief_rebase.py (its rows_for is imported,
unchanged), run on the repair_v8 and repair_v9 products, so the method is
checked against the repair_v8 numbers the brief already prints.  Read-only:
writes one table and one record.

Outputs:
  tables/wp13_brief_rebase_repair_v9.csv
  provenance/wp13_brief_rebase_repair_v9.json

Run:  PYTHONPATH=scripts python3 scripts/wp13_brief_rebase_v9.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd

import wp5_common as w
from wp13_brief_rebase import CHAINS, rows_for

OUT = w.TABLES / "wp13_brief_rebase_repair_v9.csv"
OUT_EXEC = w.PROVENANCE / "wp13_brief_rebase_repair_v9.json"
V9 = {"scan": "tables/wp7_age_sensitivity_repair_v9.csv",
      "ledger": "tables/wp7_ledger_repair_v9.csv",
      "wp5": "data/processed/wp5_imf_normalization_repair_v9.parquet"}


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; nothing is overwritten")
    chains = {"repair_v8": CHAINS["repair_v8"], "repair_v9": V9}
    rows = []
    for chain, paths in chains.items():
        rows += rows_for(chain, paths)
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_brief_rebase_v9.py",
        "purpose": "re-base tasks/wp13_pooled_vs_resolved_ablation_brief.md section 2 on repair_v9",
        "inputs": {p: w.sha256(w.ROOT / p) for c in chains.values() for p in c.values()},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
