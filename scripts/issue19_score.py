#!/usr/bin/env python3
"""Issue #19 -- score P1-P6 and tabulate every headline number before/after.

Scores the predictions of provenance/issue19_repair_v8_prereg.json exactly as
their ``pass_if`` states, against the repair_v7 artifacts hashed there and the
repair_v8 products of scripts/run_repair_v8_chain.sh.  Failed predictions are
recorded as failed.  Nothing here feeds back into any product.

Outputs:
  provenance/issue19_repair_v8_outcome.json
  tables/issue19_before_after.csv

Run:
  PYTHONPATH=scripts python3 scripts/issue19_score.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

PREREG = w.PROVENANCE / "issue19_repair_v8_prereg.json"
INTEGRITY = w.PROVENANCE / "issue19_repair_v8_integrity.json"
OUT = w.PROVENANCE / "issue19_repair_v8_outcome.json"
TABLE = w.TABLES / "issue19_before_after.csv"

KEY = ["subgroup", "family", "R_V", "alpha"]
SUB = ["CygOB2-A", "CygOB2-B", "CygOB2-C"]


def proc(name: str) -> pd.DataFrame:
    return pd.read_parquet(w.PROC / name)


def js(rel: str) -> dict:
    return json.loads((w.ROOT / rel).read_text())


def norm(v: str) -> pd.DataFrame:
    return proc(f"wp5_imf_normalization_{v}.parquet").set_index(KEY)


def ledger(v: str) -> pd.DataFrame:
    return pd.read_csv(w.TABLES / ("wp7_ledger.csv" if v == "repair_v7"
                                   else f"wp7_ledger_{v}.csv"))


def closure(v: str) -> pd.DataFrame:
    return pd.read_csv(w.TABLES / f"wp6_closure_{v}.csv")


def tagged(rel: str, v: str) -> str:
    if v == "repair_v7":
        return rel
    stem, dot, ext = rel.rpartition(".")
    return f"{stem}_{v}.{ext}"


# ------------------------------------------------------------- predictions
def p1() -> dict:
    v7, v8 = norm("repair_v7"), norm("repair_v8")
    loss = {}
    for rv in w.R_V_BRANCHES:
        idx = ("CygOB2-A", "PARSEC", rv, 2.3)
        loss[f"rv{rv:.1f}"] = round(float(
            v7.loc[idx, "membership_weighted_calibration_sources"]
            - v8.loc[idx, "membership_weighted_calibration_sources"]), 4)
    ok = all(7.5 <= x <= 8.5 for x in loss.values())
    return {"passed": ok, "measured": {"loss_by_R_V": loss}}


def p2() -> dict:
    v7, v8 = norm("repair_v7"), norm("repair_v8")
    k = (v8.k_median / v7.k_median - 1).rename("rel").reset_index()
    cells = k[k.subgroup.isin(["CygOB2-B", "CygOB2-C"]) | k.family.eq("MIST")]
    k_worst = cells.loc[cells.rel.abs().idxmax()]
    k_fail = cells[cells.rel.abs() > 0.01]

    l7, l8 = ledger("repair_v7"), ledger("repair_v8")
    cols = ["scope", "subgroup", "family", "R_V", "alpha", "sf_duration_Myr",
            "explodability"]
    j = l7.merge(l8, on=cols, suffixes=("_v7", "_v8"))
    j = j[j.explodability.eq("all_explode")]
    sel = j[(j.scope.eq("subgroup") & j.subgroup.isin(["CygOB2-B", "CygOB2-C"]))
            | j.family.eq("MIST")].copy()
    sel["abs"] = (sel.N_SN_mean_v8 - sel.N_SN_mean_v7).abs()
    sel["rel"] = sel["abs"] / sel.N_SN_mean_v7.replace(0.0, np.nan)
    sel["ok"] = np.where(sel.N_SN_mean_v7 < 1.0, sel["abs"] <= 0.01,
                         sel["rel"] <= 0.01)
    bad = sel[~sel.ok]
    material = sel[sel.N_SN_mean_v7 >= 1.0]
    worst = material.loc[material.rel.idxmax()]
    return {
        "passed": bool(k_fail.empty and bad.empty),
        "measured": {
            "k_cells": int(len(cells)),
            "k_cells_outside_1pct": k_fail[KEY + ["rel"]].round(5).to_dict("records"),
            "k_worst": {**{c: k_worst[c] for c in KEY}, "rel": round(float(k_worst.rel), 5)},
            "ledger_rows": int(len(sel)),
            "ledger_rows_outside_tolerance": int(len(bad)),
            "ledger_rows_outside_tolerance_examples": bad[
                cols + ["N_SN_mean_v7", "N_SN_mean_v8"]].head(12).round(4).to_dict("records"),
            "ledger_worst_relative_material_row": {
                **{c: worst[c] for c in cols},
                "N_SN_mean_v7": round(float(worst.N_SN_mean_v7), 4),
                "N_SN_mean_v8": round(float(worst.N_SN_mean_v8), 4),
                "rel": round(float(worst.rel), 5),
            },
        },
    }


def baseline_assoc(v: str) -> float:
    l = ledger(v)
    row = l[l.scope.eq("association") & l.family.eq("PARSEC") & l.R_V.eq(3.1)
            & l.alpha.eq(2.3) & l.sf_duration_Myr.eq(0.0)
            & l.explodability.eq("all_explode")]
    return float(row.N_SN_mean.iloc[0])


def p3() -> dict:
    value = baseline_assoc("repair_v8")
    return {"passed": bool(8.20 <= value <= 8.43),
            "measured": {"N_death_baseline_repair_v8": round(value, 4),
                         "repair_v7": round(baseline_assoc("repair_v7"), 4)}}


def p4() -> dict:
    c7 = closure("repair_v7").set_index(KEY)
    c8 = closure("repair_v8").set_index(KEY)
    rows = {}
    for rv in w.R_V_BRANCHES:
        for alpha in w.IMF_SLOPES:
            idx = ("CygOB2-A", "PARSEC", rv, alpha)
            rows[f"rv{rv:.1f}_a{alpha:.1f}"] = round(float(
                c8.loc[idx, "closure_ratio"] / c7.loc[idx, "closure_ratio"] - 1), 4)
    ok = all(0.08 <= x <= 0.20 for x in rows.values())
    return {"passed": ok, "measured": {"rise_by_cell": rows,
                                       "median_rise_not_gating": round(float(np.median(list(rows.values()))), 4),
                                       "cells_in_range": int(sum(0.08 <= x <= 0.20 for x in rows.values()))}}


def p5() -> dict:
    rec = js("provenance/wp12_closure_slopes_execution_repair_v8.json")
    value = float(rec["closing_slope_grid_median_by_subgroup"]["association"])
    old = js("provenance/wp12_closure_slopes_execution.json")
    return {"passed": bool(2.10 <= value <= 2.30),
            "measured": {"association_closing_slope_repair_v8": value,
                         "repair_v7": old["closing_slope_grid_median_by_subgroup"]["association"]}}


def p6() -> dict:
    post = proc("wp4_age_posteriors_repair_v5.parquet")
    kept = post[post.measurable.astype(bool) & ~post.grid_railed.astype(bool)]
    lo, hi = round(float(kept.age_map.min()), 2), round(float(kept.age_map.max()), 2)
    pms = int(kept.indicator.eq("pms").sum())
    return {"passed": bool(lo == 2.00 and hi == 4.01 and pms == 0),
            "measured": {"retained_rows": int(len(kept)), "min_Myr": lo,
                         "max_Myr": hi, "pms_rows": pms}}


# ----------------------------------------------------------- before/after
def before_after() -> pd.DataFrame:
    rows = []

    def add(quantity, before, after, source):
        rows.append({"quantity": quantity, "repair_v7": before, "repair_v8": after,
                     "change": (after - before) if isinstance(before, (int, float))
                     and isinstance(after, (int, float)) else None,
                     "source": source})

    # ledger
    for scope, sub in [("association", None)] + [("subgroup", s) for s in SUB]:
        vals = []
        for v in ("repair_v7", "repair_v8"):
            l = ledger(v)
            sel = l[l.scope.eq(scope) & l.family.eq("PARSEC") & l.R_V.eq(3.1)
                    & l.alpha.eq(2.3) & l.sf_duration_Myr.eq(0.0)
                    & l.explodability.eq("all_explode")]
            if sub:
                sel = sel[sel.subgroup.eq(sub)]
            vals.append(round(float(sel.N_SN_mean.iloc[0]), 3))
        add(f"N_death baseline, {sub or 'association'}", *vals, "wp7_ledger")
    for label, fn in [
        ("P(>=1 event) baseline", lambda r: r.P_at_least_one),
        ("P(last < 100 kyr) baseline", lambda r: r.P_last_SN_within_100kyr),
    ]:
        vals = []
        for v in ("repair_v7", "repair_v8"):
            l = ledger(v)
            r = l[l.scope.eq("association") & l.family.eq("PARSEC") & l.R_V.eq(3.1)
                  & l.alpha.eq(2.3) & l.sf_duration_Myr.eq(0.0)
                  & l.explodability.eq("all_explode")].iloc[0]
            vals.append(round(float(fn(r)), 4))
        add(label, *vals, "wp7_ledger")
    for label, cond in [("36 headline branches", lambda a: a != 2.6),
                        ("all 54 branches", lambda a: True)]:
        for stat in ("min", "max", "median"):
            vals = []
            for v in ("repair_v7", "repair_v8"):
                l = ledger(v)
                a = l[l.scope.eq("association") & l.explodability.eq("all_explode")]
                a = a[a.alpha.map(cond)]
                vals.append(round(float(getattr(a.N_SN_mean, stat)()), 2))
            add(f"N_death {stat}, {label}", *vals, "wp7_ledger")
    # WP5
    n7, n8 = norm("repair_v7"), norm("repair_v8")
    for s in SUB:
        idx = (s, "PARSEC", 3.1, 2.3)
        add(f"k_median baseline, {s}", round(float(n7.loc[idx, 'k_median']), 1),
            round(float(n8.loc[idx, 'k_median']), 1), "wp5_imf_normalization")
        add(f"age posterior mean baseline, {s} (Myr)",
            round(float(n7.loc[idx, 'truth_age_posterior_mean_Myr']), 3),
            round(float(n8.loc[idx, 'truth_age_posterior_mean_Myr']), 3),
            "wp5_imf_normalization")
    add("WP5 cells passing residual gate (of 54)",
        int(n7.residual_gate_pass.sum()), int(n8.residual_gate_pass.sum()),
        "wp5_imf_normalization")
    a7 = proc("wp5_association_mass_repair_v7.parquet")
    a8 = proc("wp5_association_mass_repair_v8.parquet")
    base = lambda a: a[a.family.eq("PARSEC") & a.R_V.eq(3.1) & a.alpha.eq(2.3)].iloc[0]
    add("association stellar mass baseline (Msun)",
        round(float(base(a7).multiplicity_adjusted_mass_median_Msun), 0),
        round(float(base(a8).multiplicity_adjusted_mass_median_Msun), 0),
        "wp5_association_mass")
    g7 = js("provenance/wp5_imf_fit_execution_repair_v7.json")["gate"]
    g8 = js("provenance/wp5_imf_fit_execution_repair_v8.json")["gate"]
    add("WP5 gate G3 pass", str(g7["G3_pass"]), str(g8["G3_pass"]), "wp5_imf_fit_execution")
    # WP6
    c7, c8 = closure("repair_v7"), closure("repair_v8")
    for label, frame_fn in [("closure grid median at alpha 2.3 (18 cells)",
                             lambda c: c[c.alpha.eq(2.3)].closure_ratio.median())]:
        add(label, round(float(frame_fn(c7)), 3), round(float(frame_fn(c8)), 3), "wp6_closure")
    for s in SUB:
        sel = lambda c: c[c.subgroup.eq(s) & c.family.eq("PARSEC") & c.R_V.eq(3.1)
                          & c.alpha.eq(2.3)].iloc[0]
        add(f"closure ratio baseline, {s}", round(float(sel(c7).closure_ratio), 3),
            round(float(sel(c8).closure_ratio), 3), "wp6_closure")
    cs7 = js("provenance/wp12_closure_slopes_execution.json")
    cs8 = js("provenance/wp12_closure_slopes_execution_repair_v8.json")
    for s in SUB + ["association"]:
        add(f"closing slope grid median, {s}",
            round(float(cs7["closing_slope_grid_median_by_subgroup"][s]), 3),
            round(float(cs8["closing_slope_grid_median_by_subgroup"][s]), 3),
            "wp12_closure_slopes_execution")
    add("star-weighted association closure, grid median at 2.3",
        cs7["association_closure_grid_median_at_2p3"],
        cs8["association_closure_grid_median_at_2p3"], "wp12_closure_slopes_execution")
    m7 = pd.read_csv(w.TABLES / "wp6_massive_census.csv")
    m8 = pd.read_csv(w.TABLES / "wp6_massive_census_repair_v8.csv")
    for s in SUB:
        sel = lambda m: m[m.subgroup.eq(s) & m.family.eq("PARSEC") & m.R_V.eq(3.1)].iloc[0]
        add(f"census above 8 Msun (probabilistic), {s}",
            round(float(sel(m7).observed_above_8_probabilistic), 2),
            round(float(sel(m8).observed_above_8_probabilistic), 2), "wp6_massive_census")
        add(f"census above 8 Msun (thresholded), {s}",
            int(sel(m7).observed_above_8_thresholded),
            int(sel(m8).observed_above_8_thresholded), "wp6_massive_census")
    add("living ledger above 8 Msun",
        js("provenance/wp6_ledger_execution.json")["total_living_above_8_Msun"],
        js("provenance/wp6_ledger_execution_repair_v8.json")["total_living_above_8_Msun"],
        "wp6_ledger_execution")
    # WP9 / WP12
    v7 = js("provenance/wp9_verdict_execution.json")["verdict"]
    v8 = js("provenance/wp9_verdict_execution_repair_v8.json")["verdict"]
    for k in ("P_verdict_min", "P_verdict_median", "P_verdict_max"):
        add(k, v7[k], v8[k], "wp9_verdict_execution")
    add("WP9 verdict", v7["outcome"], v8["outcome"], "wp9_verdict_execution")
    gl7 = js("provenance/wp12_gate_landscape_execution.json")
    gl8 = js("provenance/wp12_gate_landscape_execution_repair_v8.json")
    add("headline cells passing (of 36)", gl7["headline_cells"]["passing"],
        gl8["headline_cells"]["passing"], "wp12_gate_landscape_execution")
    add("headline combinations passing all subgroups (of 12)",
        gl7["combinations"]["headline_passing_all_subgroups"],
        gl8["combinations"]["headline_passing_all_subgroups"], "wp12_gate_landscape_execution")
    add("headline branches on all-passing cells (of 36)",
        gl7["headline_branches"]["built_entirely_on_passing_cells"],
        gl8["headline_branches"]["built_entirely_on_passing_cells"],
        "wp12_gate_landscape_execution")
    for a in ("alpha_2", "alpha_2.3"):
        add(f"branches with P_verdict > 0.5, {a} (of 18)",
            gl7["alpha_split"]["all_headline"][a]["above_0p5"],
            gl8["alpha_split"]["all_headline"][a]["above_0p5"],
            "wp12_gate_landscape_execution")
    return pd.DataFrame(rows)


def main() -> None:
    prereg = json.loads(PREREG.read_text())
    integrity = json.loads(INTEGRITY.read_text()) if INTEGRITY.exists() else {}
    scorers = {"P1": p1, "P2": p2, "P3": p3, "P4": p4, "P5": p5, "P6": p6}
    scored = []
    for spec in prereg["predictions"]:
        result = scorers[spec["id"]]()
        scored.append({
            "id": spec["id"],
            "statement": spec["statement"],
            "pass_if": spec["pass_if"],
            "outcome": "PASS" if result["passed"] else "FAIL",
            "measured": result["measured"],
        })
    table = before_after()
    table.to_csv(TABLE, index=False)
    checks = integrity.get("checks", {})
    adopted = all(checks.get(i, {}).get("pass") for i in ("I1", "I2", "I3", "I4"))
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue19_score.py",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "preregistration": str(PREREG.relative_to(w.ROOT)),
        "preregistration_sha256": w.sha256(PREREG),
        "preregistration_created_utc": prereg["created_utc"],
        "predictions": scored,
        "predictions_passed": sum(p["outcome"] == "PASS" for p in scored),
        "predictions_total": len(scored),
        "integrity_checks": {i: checks.get(i, {}).get("pass") for i in ("I1", "I2", "I3", "I4")},
        "adoption_rule": prereg["adoption_rule"],
        "adopted": adopted,
        "before_after_table": str(TABLE.relative_to(w.ROOT)),
        "before_after": table.to_dict("records"),
    }
    w.write_json(OUT, record)
    for p in scored:
        print(f"{p['id']}: {p['outcome']}  {json.dumps(p['measured'])[:300]}")
    print(f"integrity: {record['integrity_checks']}  -> adopted: {adopted}")
    print(table.to_string(index=False))
    print(f"\nwrote {OUT.relative_to(w.ROOT)} and {TABLE.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
