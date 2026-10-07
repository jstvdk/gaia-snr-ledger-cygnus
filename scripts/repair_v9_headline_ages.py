#!/usr/bin/env python3
"""repair_v9 -- the headline WP4 age table (owner decisions of 2026-10-07).

A and C carry their adopted Stage 1 rows unchanged: the Phase A' test-c
spectroscopic-HRD posteriors (eps 0.05).  B was not measurable (4
spectroscopic stars; its M1 photometric fit swings 1.4-5.7 Myr with R_V), so
by decision 3 it BORROWS the A + C combined posterior: the product of the two
test-c likelihood curves of the same family and R_V, i.e. the sum of their
ln L curves on the shared native age grid.  Every B row is flagged
"not measured for B".

Test c has no distance or f_bin variants, so every (f_bin, dmu) key of a
(subgroup, family, R_V) carries the same curve -- exactly as the Stage 1
table already does for A and C (fallback suffix no_dmu_refit).

The Stage 1 table data/processed/wp4_age_posteriors_repair_v9.parquet is read,
never written.  Only the adopted 'ums' rows are carried: downstream code
selects exactly one 'ums' row per (subgroup, family, R_V, f_bin=0.4, dmu=0)
and wp4_repair_common.age_posterior_nodes silently falls back to a median
over subgroups if it finds anything other than one -- this script refuses to
write a table that would trigger that.

Outputs:
  data/processed/wp4_age_posteriors_repair_v9_headline.parquet
  provenance/repair_v9_headline_ages_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/repair_v9_headline_ages.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp4v9_common as V
import wp5_common as w

STAGE1 = w.PROC / "wp4_age_posteriors_repair_v9.parquet"
OUT = w.PROC / "wp4_age_posteriors_repair_v9_headline.parquet"
OUT_JSON = w.PROVENANCE / "repair_v9_headline_ages_execution.json"
EPS = 0.05                                   # test-c outlier fraction (Stage 1 fallback)
KEY = ["subgroup", "family", "R_V", "f_bin", "dmu"]
B_FLAG = "borrowed_A_plus_C:not_measured_for_B"


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; nothing is overwritten")
    V.check_method_unchanged()
    stage1 = pd.read_parquet(STAGE1)
    curves = json.loads(V.frozen("spectroscopic_curves").read_text())["loglike_curves"]
    tests = pd.read_csv(V.frozen("spectroscopic_tests"))
    native = {f: np.array(v) for f, v in json.loads(
        (w.PROVENANCE / "wp4v9_fit_execution.json").read_text())["native_ages"].items()}

    def spec(sg: str, fam: str, rv: float) -> tuple[np.ndarray, int]:
        row = tests[tests.test.eq("c") & np.isclose(tests.eps, EPS) & tests.subgroup.eq(sg)
                    & tests.family.eq(fam) & np.isclose(tests.R_V, rv)]
        if len(row) != 1:
            raise RuntimeError(f"no unique test-c row for {sg}/{fam}/{rv}")
        curve = np.array(curves[f"c|{EPS}|{sg}|{fam}|{rv}"], dtype=float)
        if curve.shape != native[fam].shape:
            raise RuntimeError("test-c curve is not on the native age grid")
        return curve, int(row.n_stars.iloc[0])

    adopted = stage1[stage1.indicator.eq("ums")].copy()
    if adopted.duplicated(KEY).any():
        raise RuntimeError("Stage 1 table has more than one 'ums' row per key")

    rows, combined = [], {}
    for _, r in adopted.iterrows():
        if r.subgroup != "CygOB2-B":
            if r.window != "spectroscopic_hrd" or r.error_model != "test_c":
                raise RuntimeError(f"{r.subgroup} adopted row is not the test-c fallback")
            rows.append(r.to_dict())
            continue
        fam, rv = r.family, float(r.R_V)
        ca, na = spec("CygOB2-A", fam, rv)
        cc, nc = spec("CygOB2-C", fam, rv)
        total = ca + cc
        s = V.posterior(native[fam], total, na + nc)
        out = r.to_dict()
        out.update({k: s[k] for k in ("age_map", "age_lo68", "age_hi68", "age_lo90",
                                      "age_hi90", "age_mean", "grid_railed",
                                      "measurable", "exclusion_reason")})
        out.update({"n_stars": na + nc, "window": "spectroscopic_hrd_A_plus_C",
                    "error_model": "test_c_product_A_C",
                    "fallback": B_FLAG + ("" if r.dmu == 0.0 and r.f_bin == 0.4
                                          else ":no_dmu_refit"),
                    "adopted": True, "branch": "", "photometry_only": False})
        rows.append(out)
        combined[f"{fam}|{rv}"] = {
            "age_map": s["age_map"], "age_lo68": s["age_lo68"], "age_hi68": s["age_hi68"],
            "age_median": s["age_median"], "n_stars_A_plus_C": na + nc,
            "bimodal": bool(V.bimodal(native[fam], total)),
            "lnL_curve": [float(x) for x in total]}

    table = pd.DataFrame(rows)[list(stage1.columns)]
    table = table.sort_values(KEY).reset_index(drop=True)

    # one adopted row per key, and exactly one per node-rule key
    if table.duplicated(KEY).any():
        raise RuntimeError("duplicate keys in the headline table")
    node_keys = table[np.isclose(table.f_bin, 0.4) & np.isclose(table.dmu, 0.0)]
    expected = len(V.SUBGROUPS) * len(V.FAMILIES) * len(V.R_V_BRANCHES)
    if len(node_keys) != expected or node_keys.duplicated(["subgroup", "family", "R_V"]).any():
        raise RuntimeError("node-rule keys are not exactly one per subgroup x family x R_V")
    if not table.measurable.all():
        raise RuntimeError("an adopted headline age is not measurable")

    table.to_parquet(OUT, index=False)
    decision = table[np.isclose(table.R_V, 3.1) & np.isclose(table.f_bin, 0.4)
                     & np.isclose(table.dmu, 0.0)]
    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_headline_ages.py",
        "decision_record": "provenance/decisions_2026_10_07.json (decisions 6 and 8)",
        "rule": ("A and C: the Stage 1 adopted 'ums' row (test c, eps 0.05), copied. "
                 "B: product of the A and C test-c likelihood curves (sum of ln L) of the "
                 "same family and R_V, summarised with wp4v9_common.posterior; flagged "
                 f"'{B_FLAG}'."),
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (
            STAGE1, V.frozen("spectroscopic_curves"), V.frozen("spectroscopic_tests"),
            w.PROVENANCE / "wp4v9_fit_execution.json")},
        "b_combined_curves": combined,
        "decision_cell_R_V_3p1": {
            f"{r.subgroup}|{r.family}": {"age_map": r.age_map, "age_lo68": r.age_lo68,
                                         "age_hi68": r.age_hi68, "n_stars": int(r.n_stars),
                                         "window": r.window}
            for r in decision.itertuples()},
        "rows": int(len(table)),
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    for r in decision.itertuples():
        print(f"{r.subgroup} {r.family:6s} {r.age_map:.2f} ({r.age_lo68:.2f}-{r.age_hi68:.2f}) "
              f"n={r.n_stars} {r.window}")
    print(f"wrote {OUT.relative_to(w.ROOT)} ({len(table)} rows)")


if __name__ == "__main__":
    main()
