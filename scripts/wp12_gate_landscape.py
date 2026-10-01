#!/usr/bin/env python3
"""WP12.1 -- the residual-gate landscape, and the branch/gate table (Appendix A).

Two facts the manuscript currently blurs together:

  * the BASELINE branch passes the WP5 residual gate in all three subgroups;
  * 40 of 54 subgroup fit CELLS pass it.

Neither implies that the 36-branch headline set is validated.  The headline set
was selected on alpha, not on gate status, and it contains cells that failed.
This script measures exactly how many, and produces the strict
all-subgroup-pass subset the brief asks for, plus the Appendix A table.

It also records that the failure BREAKDOWN differs between repair_v6 (the
accepted gate record) and repair_v7 (the chain the ledger actually ran on).
Both give 40/54.  repair_v7 is used here and in the manuscript because it is
the chain that produced every downstream number.

Outputs:
  tables/wp12_wp5_gate_map.csv          54 cells with their gate statistics
  tables/wp12_combination_gate.csv      18 family-R_V-alpha combinations
  tables/wp12_branch_gate_table.csv     Appendix A, 54 all-explode branches
  provenance/wp12_gate_landscape_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp12_gate_landscape.py
"""
from __future__ import annotations

import pandas as pd

import chain as C
import wp5_common as w
import wp12_common as W


def breakdown(frame: pd.DataFrame) -> dict:
    fail = frame[~frame.residual_gate_pass]
    return {
        "cells_total": int(len(frame)),
        "cells_passing": int(frame.residual_gate_pass.sum()),
        "cells_failing": int(len(fail)),
        "failures_by_subgroup": {
            str(k): int(v) for k, v in fail.subgroup.value_counts().sort_index().items()
        },
        "failures_by_family": {
            str(k): int(v) for k, v in fail.family.value_counts().sort_index().items()
        },
        "failures_by_R_V": {
            f"{k:g}": int(v) for k, v in fail.R_V.value_counts().sort_index().items()
        },
        "failures_by_alpha": {
            f"{k:g}": int(v) for k, v in fail.alpha.value_counts().sort_index().items()
        },
    }


def main() -> None:
    cells = W.gate_map()
    combos = W.combination_pass()

    # The repair_v6 comparison.  Same headline count, different failing cells.
    v6 = pd.read_parquet(w.PROC / "wp5_imf_normalization_repair_v6.parquet")
    v6_breakdown = breakdown(v6)
    v7_breakdown = breakdown(cells)

    headline_cells = cells[cells.in_headline_set]
    headline_combos = combos[combos.in_headline_set]
    n_all_pass_combos = int(combos.all_subgroup_pass.sum())
    n_headline_all_pass_combos = int(headline_combos.all_subgroup_pass.sum())
    # Each combination spawns three ledger branches, one per star-formation
    # duration.  The WP5 gate does not depend on the duration.
    n_sf = 3
    all_pass_headline_branches = n_headline_all_pass_combos * n_sf

    # ------------------------------------------------ Appendix A branch table
    assoc = W.association_all_explode()
    closure = W.closure()
    verdict = W.verdict()
    sets = W.branch_sets()

    closure_wide = closure.pivot_table(
        index=["family", "R_V", "alpha"], columns="subgroup",
        values="closure_ratio",
    ).rename(columns={s: f"closure_{s[-1]}" for s in w.SUBGROUPS}).reset_index()
    # The association-level closure ratio.  WP6 defines closure_ratio as
    # OBSERVED / PREDICTED (wp6_closure_test.py:237), so the star-weighted
    # aggregate is summed observed over summed predicted -- not the mean of the
    # three subgroup ratios, which is a different quantity.
    assoc_closure = closure.groupby(["family", "R_V", "alpha"], as_index=False).apply(
        lambda g: pd.Series(
            {"closure_association":
             g.observed_living.sum() / g.predicted_observed_living.sum()}
        ),
        include_groups=False,
    )
    closure_wide = closure_wide.merge(assoc_closure, on=["family", "R_V", "alpha"])

    gate_wide = cells.pivot_table(
        index=["family", "R_V", "alpha"], columns="subgroup",
        values="residual_gate_pass", aggfunc="first",
    ).rename(columns={s: f"wp5_gate_{s[-1]}" for s in w.SUBGROUPS}).reset_index()
    gate_wide = gate_wide.merge(
        combos[["family", "R_V", "alpha", "all_subgroup_pass", "cells_passing"]],
        on=["family", "R_V", "alpha"],
    )

    keys = ["family", "R_V", "alpha", "sf_duration_Myr"]
    table = (
        assoc[keys + ["N_SN_mean", "N_SN_median", "N_SN_p16", "N_SN_p84",
                      "P_at_least_one", "P_last_SN_within_100kyr",
                      "t_last_median_Myr"]]
        .merge(gate_wide, on=["family", "R_V", "alpha"])
        .merge(closure_wide, on=["family", "R_V", "alpha"])
        .merge(
            sets[keys + ["min_turnoff_Msun", "min_dead_progenitor_Msun",
                         "fraction_of_SNe_below_52Msun", "N_SN_bh_cut_30",
                         "N_SN_bh_cut_40"]],
            on=keys,
        )
        .merge(
            verdict[keys + ["C1_age", "C1_age_permissive", "C3_stripped_fraction",
                            "C4_in_situ", "P_verdict", "P_verdict_permissive",
                            "in_headline_set"]],
            on=keys,
        )
    )
    table = table.rename(columns={
        "P_verdict": "scenario_score",
        "P_verdict_permissive": "scenario_score_permissive",
    })
    table["headline_and_all_subgroup_pass"] = (
        table.in_headline_set & table.all_subgroup_pass
    )

    def note(row: pd.Series) -> str:
        parts = []
        if not row.in_headline_set:
            parts.append("sensitivity only (alpha=2.6, excluded from headline set)")
        if not row.all_subgroup_pass:
            failed = [
                s[-1] for s in w.SUBGROUPS if not row[f"wp5_gate_{s[-1]}"]
            ]
            parts.append(f"WP5 residual gate FAILS in {', '.join(failed)}")
        if row.min_turnoff_Msun >= 119.9:
            parts.append("turnoff at the 120 Msun IMF ceiling")
        return "; ".join(parts) if parts else "baseline-quality: all cells pass"

    table["notes"] = table.apply(note, axis=1)
    is_base = (
        table.family.eq(W.BASE["family"]) & table.R_V.eq(W.BASE["R_V"])
        & table.alpha.eq(W.BASE["alpha"])
        & table.sf_duration_Myr.eq(W.BASE["sf_duration_Myr"])
    )
    table.loc[is_base, "notes"] = "BASELINE BRANCH; " + table.loc[is_base, "notes"]
    table = table.sort_values(keys).reset_index(drop=True)

    # ----------------------------------------------- the alpha split, strictly
    head = table[table.in_headline_set]
    strict = head[head.all_subgroup_pass]
    split = {}
    for label, frame in (("all_headline", head), ("all_subgroup_pass_only", strict)):
        entry = {}
        for alpha in W.HEADLINE_ALPHAS:
            arm = frame[frame.alpha.eq(alpha)]
            entry[f"alpha_{alpha:g}"] = {
                "branches": int(len(arm)),
                "above_0p5": int((arm.scenario_score > W.SUPPORT_THRESHOLD).sum()),
                "score_min": round(float(arm.scenario_score.min()), 4) if len(arm) else None,
                "score_max": round(float(arm.scenario_score.max()), 4) if len(arm) else None,
                "N_SN_min": round(float(arm.N_SN_mean.min()), 3) if len(arm) else None,
                "N_SN_max": round(float(arm.N_SN_mean.max()), 3) if len(arm) else None,
            }
        split[label] = entry

    # ------------------------------------------------------------ predictions
    r1_measured = {
        "all_pass_headline_branches": all_pass_headline_branches,
        "headline_branches": int(len(head)),
        "headline_combinations_passing_all_subgroups": n_headline_all_pass_combos,
        "headline_combinations": int(len(headline_combos)),
        "note": (
            "NOT BLIND -- this value was computed before the preregistration "
            "was written and is disclosed as prior knowledge there.  Recorded "
            "as an expectation, not as a test."
        ),
    }
    r1 = W.score_prediction("R1", all_pass_headline_branches < 18, r1_measured)

    a20 = strict[strict.alpha.eq(2.0)]
    a23 = strict[strict.alpha.eq(2.3)]
    r2_pass = bool(
        (a20.scenario_score > W.SUPPORT_THRESHOLD).sum() >= 1
        and (a23.scenario_score > W.SUPPORT_THRESHOLD).sum() == 0
    )
    r2 = W.score_prediction("R2", r2_pass, {
        "alpha_2p0_all_pass_branches": int(len(a20)),
        "alpha_2p0_above_0p5": int((a20.scenario_score > W.SUPPORT_THRESHOLD).sum()),
        "alpha_2p3_all_pass_branches": int(len(a23)),
        "alpha_2p3_above_0p5": int((a23.scenario_score > W.SUPPORT_THRESHOLD).sum()),
    })

    # -------------------------------------------------------------- write out
    out_cells = C.tag(w.TABLES / "wp12_wp5_gate_map.csv")
    out_combo = C.tag(w.TABLES / "wp12_combination_gate.csv")
    out_table = C.tag(w.TABLES / "wp12_branch_gate_table.csv")
    cells.to_csv(out_cells, index=False)
    combos.to_csv(out_combo, index=False)
    table.to_csv(out_table, index=False)

    payload = {
        "item": "WP12.1 -- residual-gate landscape and Appendix A branch table",
        "wp5_version_used": C.V["wp5"],
        "chain": C.CHAIN,
        "breakdown_key_note": (
            "the key repair_v7_breakdown is kept for its readers "
            "(wp12_tables, wp10_numbers); it holds the breakdown of "
            "wp5_version_used, which on a later chain is that chain's"
        ),
        "why_repair_v7": (
            "the ledger, closure test and verdict all consumed the repair_v7 "
            "normalization.  The accepted GATE RECORD is repair_v6 and the two "
            "agree on 40/54, but they disagree on which 14 cells fail, so the "
            "failure map quoted in the paper must be repair_v7's."
        ),
        "repair_v7_breakdown": v7_breakdown,
        "repair_v6_breakdown_for_comparison": v6_breakdown,
        "breakdowns_agree_on_total": (
            v6_breakdown["cells_passing"] == v7_breakdown["cells_passing"]
        ),
        "combinations": {
            "total": int(len(combos)),
            "passing_all_subgroups": n_all_pass_combos,
            "headline_total": int(len(headline_combos)),
            "headline_passing_all_subgroups": n_headline_all_pass_combos,
        },
        "headline_cells": {
            "total": int(len(headline_cells)),
            "passing": int(headline_cells.residual_gate_pass.sum()),
        },
        "headline_branches": {
            "total": int(len(head)),
            "built_entirely_on_passing_cells": all_pass_headline_branches,
        },
        "alpha_split": split,
        "predictions": [r1, r2],
    }
    rec = W.record(
        "scripts/wp12_gate_landscape.py", payload,
        {"cells": out_cells, "combos": out_combo, "table": out_table},
    )
    w.write_json(C.tag(w.PROVENANCE / "wp12_gate_landscape_execution.json"), rec)

    print("WP12.1 residual-gate landscape")
    print(f"  repair_v7 cells passing        : {v7_breakdown['cells_passing']}/54")
    print(f"    failures by subgroup         : {v7_breakdown['failures_by_subgroup']}")
    print(f"    failures by family           : {v7_breakdown['failures_by_family']}")
    print(f"    failures by R_V              : {v7_breakdown['failures_by_R_V']}")
    print(f"    failures by alpha            : {v7_breakdown['failures_by_alpha']}")
    print(f"  repair_v6 for comparison       : {v6_breakdown['failures_by_subgroup']}")
    print(f"  combinations all-subgroup pass : {n_all_pass_combos}/18")
    print(f"  headline cells passing         : "
          f"{int(headline_cells.residual_gate_pass.sum())}/{len(headline_cells)}")
    print(f"  headline combos all-pass       : {n_headline_all_pass_combos}/12")
    print(f"  headline branches all-pass     : {all_pass_headline_branches}/36")
    for label, entry in split.items():
        print(f"  {label}:")
        for k, v in entry.items():
            print(f"    {k}: {v['above_0p5']}/{v['branches']} above 0.5 "
                  f"(score {v['score_min']}--{v['score_max']})")
    for p in (r1, r2):
        print(f"  {p['id']} {p['outcome']}")
    print("wrote tables/wp12_wp5_gate_map.csv, wp12_combination_gate.csv, "
          "wp12_branch_gate_table.csv")


if __name__ == "__main__":
    main()
