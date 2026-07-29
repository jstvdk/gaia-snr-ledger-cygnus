#!/usr/bin/env python3
"""Run and score the repair_v7 go/no-go discriminator.

Pre-registered in provenance/wp5_fbin_discriminator_prereg.json.  Read that
first: the model, the two measurements, the 2% threshold and the decision rule
were all written down before this ran.

Three nodes (one per subgroup, at the highest-weight truth age of the reporting
branch), each injected twice on the frozen MASS_GRID from a fresh
default_rng(SEED):

  ctl  frozen constant F_BINARY = 0.40
  trt  extended f_bin(M) with the below-8 Msun rise switched on

The per-star binary threshold consumes the same single rng.random(n) draw in
both arms, so every other realization is bit-identical and realization noise
cancels in the difference.  That is why three nodes suffice.

Outputs:
  data/processed/wp5_fbindisc_{ctl,trt}_{sg}_age{age}_response.parquet
  tables/wp5_fbin_discriminator.csv
  provenance/wp5_fbin_discriminator_execution.json

Run (about 10 min):
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_fbin_discriminator.py
"""
from __future__ import annotations

import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn

import wp5_common as w
import wp5_injections_repair as R
import wp5_joint_age_fit as J
from wp3_repair_common import (
    ANCHOR_PRIOR_MODE,
    REPAIR_VERSION,
    AnchorMap,
    load_template_library,
)
from wp5_fbin_discriminator_prereg import (
    BASELINE_FAMILY,
    BASELINE_RV,
    CALIBRATION_WINDOW,
    DECISION_THRESHOLD,
    UPSCATTER_WINDOW,
    WEIGHT_ALPHA,
    extended_binary_fraction,
)

WP5_VERSION = "repair_v6"
SN_THRESHOLD_MSUN = 8.0
ARMS = ("ctl", "trt")


def arm_path(arm: str, subgroup: str, age: float) -> Path:
    label = subgroup.replace("CygOB2-", "")
    return w.PROC / (
        f"wp5_fbindisc_{arm}_{label}_age{age:.3f}".replace(".", "p")
        + "_response.parquet"
    )


def recovery_by_mass(response: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """R(recovered | true mass) -- the fraction with a finite recovered mass."""
    masses = np.sort(response["true_primary_mass"].unique())
    out = np.zeros(len(masses))
    for index, value in enumerate(masses):
        rows = response[response["true_primary_mass"].eq(value)]
        out[index] = float(np.isfinite(rows["recovered_mass"]).mean())
    return masses, out


def above_threshold_by_mass(
    response: pd.DataFrame, threshold: float
) -> tuple[np.ndarray, np.ndarray]:
    """R(estimated above threshold | true mass), over the posterior draws."""
    draw_columns = sorted(
        c for c in response.columns if c.startswith("recovered_mass_draw_")
    )
    masses = np.sort(response["true_primary_mass"].unique())
    out = np.zeros(len(masses))
    for index, value in enumerate(masses):
        rows = response[response["true_primary_mass"].eq(value)]
        active = [c for c in draw_columns if rows[c].notna().any()]
        if not active or not len(rows):
            continue
        draws = rows[active].to_numpy(float)
        finite = np.isfinite(draws)
        out[index] = np.sum(finite & (draws >= threshold)) / (
            len(active) * len(rows)
        )
    return masses, out


def imf_weighted_mean(
    masses: np.ndarray, values: np.ndarray, window: tuple[float, float]
) -> float:
    inside = (masses >= window[0]) & (masses <= window[1])
    if not inside.any():
        return float("nan")
    weights = masses[inside] ** (-WEIGHT_ALPHA)
    return float(np.sum(weights * values[inside]) / np.sum(weights))


def main() -> None:
    prereg = json.loads(
        (w.PROVENANCE / "wp5_fbin_discriminator_prereg.json").read_text(
            encoding="utf-8"
        )
    )
    if (REPAIR_VERSION, ANCHOR_PRIOR_MODE) != ("repair_v5", "kriging"):
        raise RuntimeError(
            "the discriminator holds accepted repair_v6 fixed, whose upstream "
            f"is repair_v5/kriging; got {(REPAIR_VERSION, ANCHOR_PRIOR_MODE)!r}"
        )

    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet"
    )
    native = J.native_isochrone_ages(BASELINE_FAMILY)

    plan = []
    for subgroup in w.SUBGROUPS:
        nodes = J.truth_age_nodes(
            age_posterior, subgroup, BASELINE_FAMILY, BASELINE_RV, native,
            snap=not J.uses_age_interpolation(WP5_VERSION),
        )
        age = max(nodes, key=nodes.get)
        plan.append(
            {"subgroup": subgroup, "truth_age_Myr": float(age),
             "prior_weight": float(nodes[age])}
        )
    print(f"discriminator plan: {len(plan)} nodes x {len(ARMS)} arms, "
          f"{len(w.MASS_GRID)} masses each\n")

    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = R.build_donor_pool(classifier)
    donor_pool = R.augment_donor_pool(donor_pool)
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("discriminator QMC validation failed")
    posterior = np.load(w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    repair_provenance = json.loads(
        (w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8")
    )
    branch_sigma = float(
        repair_provenance["configuration"][
            "template_branch_uncertainty_calibration"
        ][f"rv{BASELINE_RV:.1f}"]["adopted_template_branch_sigma_mag"]
    )
    _, template_magnitudes, template_weights = load_template_library()

    rows = []
    bit_checks = []
    for entry in plan:
        subgroup, age = entry["subgroup"], entry["truth_age_Myr"]
        per_arm = {}
        for arm in ARMS:
            out = arm_path(arm, subgroup, age)
            started = time.time()
            if out.exists():
                response = pd.read_parquet(out)
                elapsed = 0.0
            else:
                _, response, summary = R.inject_curve(
                    subgroup, BASELINE_FAMILY, BASELINE_RV, classifier,
                    donor_pool, donor_model, normal_points,
                    np.random.default_rng(w.SEED),
                    posterior["source_id"].astype("int64"),
                    posterior["probability"], AnchorMap.from_frozen_wp3(),
                    template_magnitudes, template_weights, branch_sigma,
                    age_posterior,
                    truth_age_override=age,
                    interpolate_truth_age=True,
                    truth_binary_fraction=(
                        None if arm == "ctl" else extended_binary_fraction
                    ),
                )
                response.to_parquet(out, index=False)
                elapsed = time.time() - started
                print(f"  {subgroup} {arm}: f_bin realized "
                      f"{summary['binary_fraction_realized']:.3f}, "
                      f"{elapsed:.0f}s", flush=True)
            per_arm[arm] = response

        # Bit-preservation: the control arm differs from the accepted repair_v6
        # node in nothing at all, so it must reproduce it exactly.
        accepted = J.node_response_path(
            subgroup, BASELINE_FAMILY, BASELINE_RV, age, WP5_VERSION
        )
        match = None
        if accepted.exists():
            match = bool(w.sha256(arm_path("ctl", subgroup, age)) == w.sha256(accepted))
        bit_checks.append(
            {"subgroup": subgroup, "accepted_artifact": str(
                accepted.relative_to(w.ROOT)) if accepted.exists() else None,
             "control_reproduces_it": match}
        )

        measured = {}
        for arm in ARMS:
            masses, recovered = recovery_by_mass(per_arm[arm])
            _, above = above_threshold_by_mass(per_arm[arm], SN_THRESHOLD_MSUN)
            measured[arm] = {
                "D1": imf_weighted_mean(masses, recovered, CALIBRATION_WINDOW),
                "D2": imf_weighted_mean(masses, above, UPSCATTER_WINDOW),
            }
        rows.append(
            {
                "subgroup": subgroup,
                "truth_age_Myr": round(age, 3),
                "prior_weight": round(entry["prior_weight"], 4),
                "D1_ctl": measured["ctl"]["D1"], "D1_trt": measured["trt"]["D1"],
                "D2_ctl": measured["ctl"]["D2"], "D2_trt": measured["trt"]["D2"],
                "D1_shift_fraction": (
                    measured["trt"]["D1"] / measured["ctl"]["D1"] - 1.0
                ),
                "D2_shift_fraction": (
                    measured["trt"]["D2"] / measured["ctl"]["D2"] - 1.0
                ),
            }
        )

    table = pd.DataFrame(rows)
    out_csv = w.TABLES / "wp5_fbin_discriminator.csv"
    table.to_csv(out_csv, index=False)

    d1 = float(table.D1_shift_fraction.mean())
    d2 = float(table.D2_shift_fraction.mean())
    d1_spread = float(table.D1_shift_fraction.std())
    d2_spread = float(table.D2_shift_fraction.std())
    justified = bool(abs(d1) > DECISION_THRESHOLD or abs(d2) > DECISION_THRESHOLD)

    g1_pass = bool(d1 > 0 and d2 > 0)
    g2_pass = bool(abs(d2) > abs(d1))

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_fbin_discriminator.py",
        "status": "SUCCESS",
        "work_package": "repair_v7 go/no-go discriminator",
        "prereg": "provenance/wp5_fbin_discriminator_prereg.json",
        "prereg_created_utc": prereg["created_utc"],
        "instrument_not_result": (
            "a go/no-go test.  It authorizes or declines a repair_v7; it "
            "changes no number in the accepted chain either way."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__,
        },
        "branch": f"{BASELINE_FAMILY} R_V={BASELINE_RV}",
        "nodes": plan,
        "qmc_validation": validation,
        "bit_preservation": {
            "claim": (
                "the control arm differs from the accepted repair_v6 node "
                "response in nothing, so it must reproduce it byte for byte"
            ),
            "checks": bit_checks,
            "all_match": bool(
                all(c["control_reproduces_it"] for c in bit_checks
                    if c["control_reproduces_it"] is not None)
            ),
        },
        "per_node": rows,
        "results": {
            "D1_mean_shift_fraction": round(d1, 5),
            "D1_node_spread": round(d1_spread, 5),
            "D2_mean_shift_fraction": round(d2, 5),
            "D2_node_spread": round(d2_spread, 5),
            "threshold": DECISION_THRESHOLD,
        },
        "predictions": {
            "G1": {
                "statement": prereg["predictions"][0]["statement"],
                "pass": g1_pass,
                "evidence": f"D1 {d1:+.4f}, D2 {d2:+.4f}",
            },
            "G2": {
                "statement": prereg["predictions"][1]["statement"],
                "pass": g2_pass,
                "evidence": f"|D2| {abs(d2):.4f} vs |D1| {abs(d1):.4f}",
            },
        },
        "decision": {
            "repair_v7_justified": justified,
            "rule_applied": prereg["decision_rule"]["rule"],
            "reasoning": (
                f"D1 = {d1:+.2%} and D2 = {d2:+.2%} against a pre-declared 2% "
                + ("threshold; at least one exceeds it, so repair_v7 is "
                   "JUSTIFIED and should be scheduled."
                   if justified else
                   "threshold; both stay below it, so the effect is recorded "
                   "as a carried systematic and the accepted chain stands "
                   "unchanged.")
            ),
            "conservatism_note": (
                "the model holds f_bin = 0.40 at 2 Msun where Duchene & Kraus "
                "2013 suggest ~50%, so the measured shifts are LOWER BOUNDS.  "
                + ("This strengthens a positive decision."
                   if justified else
                   "A negative decision is therefore weaker than it looks and "
                   "should be revisited if the low-mass anchor is ever raised.")
            ),
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet",
                w.PROVENANCE / "wp5_fbin_discriminator_prereg.json",
            ]
        },
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp5_fbin_discriminator_execution.json", record)

    print("\nrepair_v7 discriminator — scored against the prereg\n")
    print(f"  {'subgroup':12s} {'D1 ctl':>8s} {'D1 trt':>8s} {'shift':>8s} "
          f"{'D2 ctl':>8s} {'D2 trt':>8s} {'shift':>8s}")
    for row in rows:
        print(f"  {row['subgroup']:12s} {row['D1_ctl']:8.4f} {row['D1_trt']:8.4f} "
              f"{row['D1_shift_fraction']:+8.2%} {row['D2_ctl']:8.4f} "
              f"{row['D2_trt']:8.4f} {row['D2_shift_fraction']:+8.2%}")
    print(f"\n  D1 (k / calibration window) {d1:+.2%}  +- {d1_spread:.2%}")
    print(f"  D2 (up-scatter 4-8 Msun)    {d2:+.2%}  +- {d2_spread:.2%}")
    print(f"  threshold {DECISION_THRESHOLD:.0%}")
    print(f"\n  G1 (both positive): {'PASS' if g1_pass else 'FAIL'}")
    print(f"  G2 (D2 > D1):       {'PASS' if g2_pass else 'FAIL'}")
    print(f"  bit-preservation:   "
          f"{'PASS' if record['bit_preservation']['all_match'] else 'FAIL'}")
    print(f"\n  DECISION: repair_v7 "
          f"{'JUSTIFIED' if justified else 'NOT justified'}")
    print("\nwrote provenance/wp5_fbin_discriminator_execution.json")


if __name__ == "__main__":
    main()
