#!/usr/bin/env python3
"""Re-base the WP13 brief's disclosed prior knowledge (section 2) on repair_v8.

The brief's section 2.1 interpolates the WP7 common-age scan and section 2.2
quotes the WP5 counts-based age ranges; both were computed on repair_v7.  This
script recomputes them from the repair_v8 products, and recomputes the
repair_v7 rows by the same code so the method is checked against the numbers
the brief already prints.  Read-only: it writes one table and one record.

Outputs:
  tables/wp13_brief_rebase_repair_v8.csv
  provenance/wp13_brief_rebase_repair_v8.json

Run:
  PYTHONPATH=scripts python3 scripts/wp13_brief_rebase.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

OUT = w.TABLES / "wp13_brief_rebase_repair_v8.csv"
OUT_EXEC = w.PROVENANCE / "wp13_brief_rebase_repair_v8.json"

CHAINS = {
    "repair_v7": {"scan": "tables/wp7_age_sensitivity.csv",
                  "ledger": "tables/wp7_ledger.csv",
                  "wp5": "data/processed/wp5_imf_normalization_repair_v7.parquet"},
    "repair_v8": {"scan": "tables/wp7_age_sensitivity_repair_v8.csv",
                  "ledger": "tables/wp7_ledger_repair_v8.csv",
                  "wp5": "data/processed/wp5_imf_normalization_repair_v8.parquet"},
}


def rows_for(chain: str, paths: dict) -> list[dict]:
    scan = pd.read_csv(w.ROOT / paths["scan"])
    ages = scan.assumed_age_Myr.to_numpy(float)
    n_of = scan.N_SN_mean.to_numpy(float)
    p_of = scan.P_last_SN_within_100kyr.to_numpy(float)
    ledger = pd.read_csv(w.ROOT / paths["ledger"])
    base = ledger[ledger.scope.eq("association") & ledger.family.eq("PARSEC")
                  & ledger.R_V.eq(3.1) & ledger.alpha.eq(2.3)
                  & ledger.sf_duration_Myr.eq(0) & ledger.explodability.eq("all_explode")]
    n_m1 = float(base.N_SN_mean.iloc[0])
    p_m1 = float(base.P_last_SN_within_100kyr.iloc[0])
    wp5 = pd.read_parquet(w.ROOT / paths["wp5"])
    cell = wp5[wp5.family.eq("PARSEC") & wp5.R_V.eq(3.1) & wp5.alpha.eq(2.3)].set_index("subgroup")
    k = cell.k_median
    age = cell.truth_age_posterior_mean_Myr
    k_mean_age = float((k * age).sum() / k.sum())
    rising = ages >= 3.0                       # N is 0 below 3 Myr; invert on the rising part
    age_for_n = float(np.interp(n_m1, n_of[rising], ages[rising]))
    out = [
        {"chain": chain, "row": "resolved baseline (M1)",
         "age_Myr": "/".join(f"{age[s]:.2f}" for s in ("CygOB2-A", "CygOB2-B", "CygOB2-C")),
         "N_death": n_m1, "P_last_lt_100kyr": p_m1},
        {"chain": chain, "row": "common age forced to 4.0", "age_Myr": "4.000",
         "N_death": float(np.interp(4.0, ages, n_of)),
         "P_last_lt_100kyr": float(np.interp(4.0, ages, p_of))},
        {"chain": chain, "row": "common age = k-weighted mean of the counts-based ages",
         "age_Myr": f"{k_mean_age:.3f}", "N_death": float(np.interp(k_mean_age, ages, n_of)),
         "P_last_lt_100kyr": float(np.interp(k_mean_age, ages, p_of))},
        {"chain": chain, "row": "common age that reproduces the resolved count",
         "age_Myr": f"{age_for_n:.3f}", "N_death": n_m1,
         "P_last_lt_100kyr": float(np.interp(age_for_n, ages, p_of))},
        {"chain": chain, "row": "scan P(last<100kyr) range over 3.25-6 Myr",
         "age_Myr": "3.25-6.00", "N_death": np.nan,
         "P_last_lt_100kyr_min": float(p_of[ages >= 3.25].min()),
         "P_last_lt_100kyr_max": float(p_of[ages >= 3.25].max())},
    ]
    ranges = wp5.groupby("subgroup").truth_age_posterior_mean_Myr.agg(["min", "max"])
    for subgroup, r in ranges.iterrows():
        out.append({"chain": chain, "row": f"counts-based age range, 18 cells, {subgroup}",
                    "age_min_Myr": float(r["min"]), "age_max_Myr": float(r["max"])})
    return out


def main() -> None:
    rows = []
    for chain, paths in CHAINS.items():
        rows += rows_for(chain, paths)
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_brief_rebase.py",
        "purpose": "re-base tasks/wp13_pooled_vs_resolved_ablation_brief.md section 2 on repair_v8",
        "inputs": {p: w.sha256(w.ROOT / p) for c in CHAINS.values() for p in c.values()},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
