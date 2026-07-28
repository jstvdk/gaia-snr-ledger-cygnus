#!/usr/bin/env python3
"""Step 5 closure for the #1c gated plan: repair_v4 gate record and reports.

Writes the versioned WP5 gate record and, when the gate passes, the
superseding completion report.  Nothing is overwritten: the frozen
``wp5_completion_report.md`` keeps its blocked-gate text and gains a
supersession banner pointing at the new report, the same pattern used for the
superseded task briefs.

``downstream_wp6_authorized`` reflects the WP5 products only.  Issue #3
(bright-mass completeness plateaus near 0.8, not 1.0) is carried explicitly as
a blocking precondition for WP6 step 2(a) so authorization here cannot be
misread as clearance to run WP6's closure test unchanged.

Run:
  PYTHONPATH=scripts python3 scripts/wp5_repair_v4_finalize.py \
      --wp5-version repair_v4 --upstream-version repair_v3
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

SUPERSESSION_BANNER = """> **SUPERSEDED 2026-07-27.** This report describes the frozen pre-repair WP5
> run (0/54 branches passing) and is retained as a historical record. The
> current WP5 verdict is in
> [wp5_completion_report_{version}.md](wp5_completion_report_{version}.md);
> the repair history is in [PROJECT_TRACE.md](PROJECT_TRACE.md) §8.

"""


def both_readings(version: str, upstream: str) -> dict:
    """Gate G3's no-regression clause under both readings (CUTS section 14.7).

    The strict per-branch reading chosen on 2026-07-27 stays binding and is
    what ``accepted`` is computed from.  Section 14.7 additionally requires
    that the refined reading -- in which a cell indeterminate in BOTH versions
    counts as neither a regression nor an improvement -- be reported alongside
    it, so that a divergence between them can never pass unnoticed.  Produced
    by scripts/wp5_verdict_stability.py, which must be run before this script.
    """
    path = w.PROVENANCE / "wp5_verdict_stability_execution.json"
    if not path.exists():
        return {
            "available": False,
            "note": (
                "run scripts/wp5_verdict_stability.py to populate the refined "
                "reading required by CUTS_AND_THRESHOLDS.md section 14.7"
            ),
        }
    record = json.loads(path.read_text(encoding="utf-8"))
    blocks = [
        block
        for block in record.get("gate_G3_no_regression_clause", [])
        if block["comparison"].endswith(f"-> {version}")
    ]
    if not blocks:
        return {
            "available": False,
            "note": (
                f"wp5_verdict_stability_execution.json carries no comparison "
                f"ending at {version}; re-run it including that version"
            ),
        }
    block = blocks[-1]
    strict = block["strict_reading"]
    refined = block["refined_reading_section_14_7"]
    return {
        "available": True,
        "source": "provenance/wp5_verdict_stability_execution.json",
        "comparison": block["comparison"],
        "strict_reading_passes": strict["passes"],
        "strict_reading_regressions": strict["count"],
        "refined_reading_passes": refined["passes"],
        "refined_reading_regressions": refined["count"],
        "discounted_as_indeterminate_in_both": refined[
            "discounted_as_indeterminate_in_both"
        ],
        "readings_agree": block["readings_agree"],
        "binding": (
            "the STRICT reading remains binding; the refined reading is "
            "reported for transparency and, where the two disagree, WP6 "
            "authorization must not be flipped on the refined reading alone "
            "without that being stated explicitly"
        ),
        "surviving_regressions": refined["regressions"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wp5-version", default="repair_v4")
    parser.add_argument("--upstream-version", default="repair_v3")
    args = parser.parse_args()
    version = args.wp5_version
    upstream = args.upstream_version

    fit_record = json.loads(
        (w.PROVENANCE / f"wp5_imf_fit_execution_{version}.json").read_text(
            encoding="utf-8"
        )
    )
    gate = fit_record["gate"]
    normalization = pd.read_parquet(
        w.PROC / f"wp5_imf_normalization_{version}.parquet"
    )
    association = pd.read_parquet(
        w.PROC / f"wp5_association_mass_{version}.parquet"
    )
    baseline = normalization[
        normalization.family.eq("PARSEC")
        & normalization.R_V.eq(3.1)
        & normalization.alpha.eq(2.3)
    ]
    baseline_assoc = association[
        association.family.eq("PARSEC")
        & association.R_V.eq(3.1)
        & association.alpha.eq(2.3)
    ].iloc[0]
    accepted = bool(gate["G3_pass"])

    unmet = [
        name
        for name, ok in [
            ("baseline_all_subgroups_residual_gate", gate["baseline_all_subgroups_residual_gate"]),
            ("no_A_or_C_branch_regression", gate["no_A_or_C_branch_regression"]),
            ("baseline_mass_within_factor_two", gate["baseline_mass_within_factor_two"]),
        ]
        if not ok
    ]
    blocking_reason = (
        None
        if accepted
        else (
            "G3 sub-criteria not met: "
            + ", ".join(unmet)
            + ".  Adopted reading of 'A and C must not regress' is the STRICT "
            "per-branch one (science decision, 2026-07-27).  The two flagged "
            "CygOB2-C cells (MIST, R_V=3.1, alpha=2.0 and 2.3) are "
            "single-truth-age-node cells whose estimator is provably identical "
            "to repair_v3, and their flips are measured Monte-Carlo noise on "
            "the six-bin rank trend statistic, not an effect of the age "
            "marginalization — see "
            "provenance/wp5_trend_stability_check_execution.json.  Under the "
            "strict reading this still blocks acceptance, and the trend "
            "diagnostic itself is the next work item (issue #11)."
        )
    )

    passing = int(normalization.residual_gate_pass.sum())
    failures = normalization[~normalization.residual_gate_pass]
    failure_corners = (
        failures.groupby(["R_V", "alpha"]).size().sort_values(ascending=False)
        if len(failures)
        else pd.Series(dtype=int)
    )

    record = {
        "wp5_version": version,
        "upstream_repair_version": upstream,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_repair_v4_finalize.py",
        "issue": "#1c step 5 closure",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "cause_of_the_blocking_failure": (
            "The WP5 injection truth model generated synthetic photometry at a "
            "single age (the upper-MS MAP).  At CygOB2-B's MAP age the "
            "PMS/Henyey isochrone fold sits on the 3.17-4.0 Msun bin boundary, "
            "so a truth-age error displaced the model caustic and produced a "
            "sharp, subgroup-specific, mass-localized residual that survived "
            "every instrumental re-run.  Evidence: "
            "provenance/wp5_bump_shape_diagnostic_execution.json (no single "
            "IMF slope explains it; the fold location), "
            "provenance/wp5_age_scan_execution.json (the residual is strictly "
            "monotone in truth age and vanishes toward older ages), "
            "provenance/wp5_stable_label_refit_execution.json (labels "
            "exonerated), provenance/wp5_lower_edge_scan_execution.json "
            "(completeness ramp excluded)."
        ),
        "fix": (
            "Truth-side age marginalization by a joint age-k fit with the WP4 "
            "age posterior as prior, applied identically to all three "
            "subgroups.  No gate threshold was moved; no per-subgroup age was "
            "chosen; with a single truth-age node the estimator reduces to the "
            "unmodified wp5_fit_imf.fit_one bit-for-bit."
        ),
        "baseline_definition": "PARSEC, R_V=3.1, alpha=2.3",
        "baseline_subgroups": [
            {
                "subgroup": row.subgroup,
                "poisson_chi_square_p": float(row.poisson_chi_square_p),
                "residual_trend_p": float(row.residual_trend_p),
                "max_abs_pearson_residual": float(row.max_abs_pearson_residual),
                "k_median": float(row.k_median),
                "truth_age_posterior_mean_Myr": float(
                    row.truth_age_posterior_mean_Myr
                ),
                "residual_gate_pass": bool(row.residual_gate_pass),
            }
            for row in baseline.itertuples()
        ],
        "branch_grid": {
            "passing": passing,
            "total": int(len(normalization)),
            "policy": "CUTS_AND_THRESHOLDS.md section 13 (issue #5)",
            "failure_concentration": {
                f"R_V={rv},alpha={alpha}": int(count)
                for (rv, alpha), count in failure_corners.items()
            },
        },
        "association_mass_median_Msun": float(
            baseline_assoc.multiplicity_adjusted_mass_median_Msun
        ),
        "association_mass_lo68_Msun": float(
            baseline_assoc.multiplicity_adjusted_mass_lo68_Msun
        ),
        "association_mass_hi68_Msun": float(
            baseline_assoc.multiplicity_adjusted_mass_hi68_Msun
        ),
        "within_factor_two_literature": bool(
            baseline_assoc.within_factor_two_literature
        ),
        "gate_G3": gate,
        "gate_G3_reading": (
            "STRICT per-branch reading of 'A and C must not regress' (science "
            "decision, 2026-07-27).  The alternative baseline-only reading "
            "would pass, and is deliberately not adopted."
        ),
        "gate_G3_both_readings": both_readings(version, upstream),
        "accepted": accepted,
        "blocking_reason": blocking_reason,
        "model_change_verdict": (
            "The fix itself works and is not in question: the baseline branch "
            "passes for all three subgroups, CygOB2-B's worst residual falls "
            "3.62 -> 2.51 with k moving only -3.6%, and the branch grid "
            "improves 26 -> 29 of 54.  What blocks acceptance is a measurement "
            "property of the trend diagnostic, not the age marginalization."
        ),
        "downstream_wp6_authorized": accepted,
        "blocking_preconditions_for_wp6": [
            {
                "issue": 3,
                "statement": (
                    "Bright-mass completeness plateaus near 0.8, not 1.0.  WP6 "
                    "step 2(a) must divide by the injection response or it will "
                    "manufacture a spurious ~20% massive-star deficit and bias "
                    "N_SN high."
                ),
                "status": "OPEN — must be handled before WP6 runs",
            }
        ],
        "carried_limitations": [
            {
                "issue": "1d",
                "statement": (
                    "CygOB2-B's A_V is set almost entirely by broadband "
                    "photometry (no spectroscopic anchors near B), so the "
                    "Teff/A_V degeneracy the anchors were meant to break is "
                    "reintroduced.  An extinction-scale error moves B's stars "
                    "across the isochrone fold coherently."
                ),
            },
            {
                "issue": 9,
                "statement": (
                    "The fitted truth-age posterior for CygOB2-B sits above its "
                    "upper-MS MAP, in the direction of the documented "
                    "upper-MS/PMS indicator disagreement.  The joint fit uses "
                    "the same counts to weight the age nodes and to score the "
                    "gate, and the chi-square dof is not reduced for the "
                    "effectively fitted age."
                ),
            },
            {
                "issue": 4,
                "statement": (
                    "No absolute 95% completeness edge exists on this field; "
                    "the corrected_no_absolute95_edge fallback is in force."
                ),
            },
        ],
        "preserved_history": {
            "frozen_wp5_products_overwritten": False,
            "repair_v1_v2_v3_products_overwritten": False,
            "upstream_wp3_wp4_rerun": False,
            "note": (
                "repair_v4 is a WP5-only version; WP3 extinction and WP4 ages "
                "and masses are the repair_v3 products, consumed unchanged."
            ),
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp5_imf_normalization_{version}.parquet",
                w.PROC / f"wp5_mass_function_bins_{version}.parquet",
                w.PROC / f"wp5_association_mass_{version}.parquet",
                w.PROC / f"wp5_completeness_curves_{version}.parquet",
                w.PROVENANCE / f"wp5_imf_fit_execution_{version}.json",
                w.PROVENANCE / f"wp5_injections_agenodes_execution_{version}.json",
                w.PROVENANCE / "wp5_age_scan_execution.json",
                w.PROVENANCE / "wp5_joint_fit_baseline_check_execution.json",
                w.PROVENANCE / "wp5_stable_label_refit_execution.json",
                w.PROVENANCE / "wp2_label_stability_execution.json",
            ]
        },
    }
    gate_path = w.PROVENANCE / f"wp5_{version}_gate.json"
    w.write_json(gate_path, record)

    rows = "\n".join(
        f"| {entry['subgroup']} | {entry['poisson_chi_square_p']:.4f} | "
        f"{entry['residual_trend_p']:.3f} | "
        f"{entry['max_abs_pearson_residual']:.2f} | "
        f"{entry['truth_age_posterior_mean_Myr']:.3f} | "
        f"{'yes' if entry['residual_gate_pass'] else 'no'} |"
        for entry in record["baseline_subgroups"]
    )
    verdict = "ACCEPTED" if accepted else "BLOCKED"
    completion = f"""# WP5 completion report ({version})

**Verdict: {verdict}.**

This report supersedes [wp5_completion_report.md](wp5_completion_report.md)
(frozen pre-repair run, 0/54 branches passing). Gate record:
[provenance/wp5_{version}_gate.json](provenance/wp5_{version}_gate.json).

## Cause of the blocking failure (issue #1c)

{record['cause_of_the_blocking_failure']}

## Fix

{record['fix']}

## Baseline gate (PARSEC, R_V = 3.1, α = 2.3)

| subgroup | χ² p | trend p | max abs residual | fitted truth age (Myr) | pass |
|---|---:|---:|---:|---:|:--:|
{rows}

## Why this verdict

{record['model_change_verdict']}

{('**Blocking reason.** ' + blocking_reason) if blocking_reason else '**All G3 sub-criteria met.**'}

Branch grid: **{passing}/{len(normalization)}** subgroup-branch fits pass,
under the retention policy in `CUTS_AND_THRESHOLDS.md` §13.
Association mass {record['association_mass_median_Msun']:.0f}
[{record['association_mass_lo68_Msun']:.0f},
{record['association_mass_hi68_Msun']:.0f}] M☉, within a factor two of the
16,500 M☉ literature scale: {record['within_factor_two_literature']}.

## Downstream authorization

`downstream_wp6_authorized = {str(accepted).lower()}`.

**WP6 may not start until open issue #3 is handled**: bright-mass completeness
plateaus near 0.8, not 1.0, so WP6 step 2(a) must divide by the injection
response or it will manufacture a spurious ~20% massive-star deficit and bias
N_SN high.

## Carried limitations

{chr(10).join(f"- **Issue #{item['issue']}** — {item['statement']}" for item in record['carried_limitations'])}

## Reproduction

```bash
bash scripts/run_repair_v4_chain.sh                      # nodes + joint fit
PYTHONPATH=scripts python3 scripts/wp5_report.py --wp5-version {version}
PYTHONPATH=scripts python3 scripts/wp5_repair_v4_finalize.py --wp5-version {version}
```

Full evidence chain: `reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md` and the
provenance records listed in the gate JSON.
"""
    completion_path = w.ROOT / f"wp5_completion_report_{version}.md"
    completion_path.write_text(completion, encoding="utf-8")

    frozen = w.ROOT / "wp5_completion_report.md"
    text = frozen.read_text(encoding="utf-8")
    if accepted and "SUPERSEDED" not in text:
        frozen.write_text(
            SUPERSESSION_BANNER.format(version=version) + text, encoding="utf-8"
        )

    print(json.dumps({
        "accepted": accepted,
        "downstream_wp6_authorized": accepted,
        "branch_grid": record["branch_grid"],
        "baseline_subgroups": record["baseline_subgroups"],
    }, indent=2))
    print(f"wrote {gate_path.relative_to(w.ROOT)}")
    print(f"wrote {completion_path.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
