#!/usr/bin/env python3
"""Step 3c, part 1: generate the truth-age node injections for repair_v4.

For every (subgroup, family, R_V) branch this generates one injection response
per truth-age node of the uniform rule in ``wp5_joint_age_fit.truth_age_nodes``
-- the recovery side's own nine-node WP4 posterior discretization snapped to
native isochrone ages.  Each node of a branch is generated with a fresh
``default_rng(wp5_common.SEED)``, so the nodes of a branch share their donor,
binary and extinction realization and differ only in truth age; that is the
same paired recipe the gate-G2 scan used, and its CygOB2-B baseline nodes are
reused rather than recomputed.

Nothing upstream of WP5 changes: repair_v4 consumes the repair_v3 WP3/WP4
products unchanged, exactly as repair_v2 consumed repair_v1.

Outputs (new files only):
  data/processed/wp5_agenode_{sg}_{family}_rv{rv}_age{age}_{version}_{response,curve}.parquet
  provenance/wp5_injections_agenodes_execution_{version}.json

Run (about 3 hours for the full grid):
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_injections_agenodes.py \
      --output-version repair_v4
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
import sklearn

import wp5_common as w
import wp5_injections_repair as R
from wp5_fbin_discriminator_prereg import extended_binary_fraction
import wp5_joint_age_fit as J
from wp3_repair_common import ANCHOR_PRIOR_MODE, REPAIR_VERSION, AnchorMap, load_template_library


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-version", default="repair_v4")
    parser.add_argument(
        "--baseline-only",
        action="store_true",
        help="restrict to the PARSEC R_V=3.1 baseline branch",
    )
    parser.add_argument(
        "--fbin-model",
        choices=["constant", "extended"],
        default="constant",
        help=(
            "truth-side binary fraction.  'constant' is the frozen F_BINARY at "
            "every mass (repair_v1..v6).  'extended' switches on the "
            "mass-dependent rise below and above 8 Msun that the pre-registered "
            "discriminator justified (repair_v7); see "
            "provenance/wp5_fbin_discriminator_execution.json."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report the node plan and exit without injecting",
    )
    args = parser.parse_args()
    version = args.output_version
    supported = {("repair_v3", "variogram"), ("repair_v5", "kriging")}
    if (REPAIR_VERSION, ANCHOR_PRIOR_MODE) not in supported:
        raise RuntimeError(
            "upstream/prior-mode combination not supported; expected one of "
            f"{sorted(supported)}, got {(REPAIR_VERSION, ANCHOR_PRIOR_MODE)!r}"
        )

    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet"
    )
    families = ["PARSEC"] if args.baseline_only else w.FAMILIES
    rv_branches = [3.1] if args.baseline_only else w.R_V_BRANCHES
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    interpolate = J.uses_age_interpolation(version)
    # repair_v7 only: the truth-side binary fraction becomes mass dependent.
    # None reproduces repair_v1..v6 bit for bit (verified by V4).
    fbin_model = (
        extended_binary_fraction if args.fbin_model == "extended" else None
    )
    plan = []
    for family in families:
        for rv in rv_branches:
            for subgroup in w.SUBGROUPS:
                nodes = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family],
                    snap=not interpolate,
                )
                if len({f"{age:.3f}" for age in nodes} ) != len(nodes):
                    raise RuntimeError(
                        f"{subgroup}/{family}/R_V={rv}: truth-age nodes collide at "
                        "the 3-decimal precision used by the node file stem"
                    )
                for age, weight in nodes.items():
                    reuse = (
                        None
                        if interpolate
                        else J.reusable_scan_snapshot(subgroup, family, rv, age)
                    )
                    target = J.node_response_path(subgroup, family, rv, age, version)
                    plan.append(
                        {
                            "subgroup": subgroup,
                            "family": family,
                            "R_V": rv,
                            "truth_age_Myr": age,
                            "prior_weight": weight,
                            "reuse_scan_snapshot": (
                                str(reuse[0].relative_to(w.ROOT)) if reuse else None
                            ),
                            "already_generated": target.exists(),
                        }
                    )
    to_run = [
        entry
        for entry in plan
        if entry["reuse_scan_snapshot"] is None and not entry["already_generated"]
    ]
    print(f"node plan: {len(plan)} nodes total, {len(to_run)} require injection")
    for entry in plan:
        status = (
            "reuse-scan"
            if entry["reuse_scan_snapshot"]
            else ("present" if entry["already_generated"] else "INJECT")
        )
        print(
            f"  {entry['family']:6s} rv{entry['R_V']} {entry['subgroup'][-1]} "
            f"age={entry['truth_age_Myr']:.3f} w={entry['prior_weight']:.3f}  {status}"
        )
    if args.dry_run:
        return

    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = R.build_donor_pool(classifier)
    donor_pool = R.augment_donor_pool(donor_pool)
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("repair injection QMC validation failed")
    posterior = np.load(w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    posterior_ids = posterior["source_id"].astype("int64")
    posterior_cube = posterior["probability"]
    anchor_map = AnchorMap.from_frozen_wp3()
    _, template_magnitudes, template_weights = load_template_library()
    repair_provenance = json.loads(
        (w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8")
    )
    branch_sigma = {
        rv: float(
            repair_provenance["configuration"][
                "template_branch_uncertainty_calibration"
            ][f"rv{rv:.1f}"]["adopted_template_branch_sigma_mag"]
        )
        for rv in w.R_V_BRANCHES
    }

    generated = []
    for index, entry in enumerate(to_run, start=1):
        subgroup = entry["subgroup"]
        family = entry["family"]
        rv = entry["R_V"]
        age = entry["truth_age_Myr"]
        started = time.time()
        print(
            f"=== [{index}/{len(to_run)}] {subgroup} {family} R_V={rv} "
            f"truth age {age:.3f} Myr ===",
            flush=True,
        )
        curve, response, summary = R.inject_curve(
            subgroup,
            family,
            rv,
            classifier,
            donor_pool,
            donor_model,
            normal_points,
            np.random.default_rng(w.SEED),
            posterior_ids,
            posterior_cube,
            anchor_map,
            template_magnitudes,
            template_weights,
            branch_sigma[rv],
            age_posterior,
            truth_age_override=age,
            interpolate_truth_age=interpolate,
            truth_binary_fraction=fbin_model,
        )
        response_path = J.node_response_path(subgroup, family, rv, age, version)
        curve_path = J.node_curve_path(subgroup, family, rv, age, version)
        response.to_parquet(response_path, index=False)
        curve.to_parquet(curve_path, index=False)
        generated.append(
            {
                **entry,
                "native_isochrone_age_Myr": float(summary["age_isochrone_Myr"]),
                "membership_pass": int(summary["membership_pass"]),
                "mass_recovered": int(summary["mass_recovered"]),
                "elapsed_seconds": round(time.time() - started, 1),
                "outputs": {
                    str(response_path.relative_to(w.ROOT)): w.sha256(response_path),
                    str(curve_path.relative_to(w.ROOT)): w.sha256(curve_path),
                },
            }
        )

    inputs = [
        w.PROC / "wp1_gaia_narrow.parquet",
        w.PROC / "wp1_2mass_join.parquet",
        w.PROC / "wp2_members.parquet",
        w.TABLES / "wp2_subgroup_labels.parquet",
        w.PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz",
        w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet",
        w.PROVENANCE / "wp2_membership_manifest.json",
        w.PROVENANCE / "wp3_repair_execution.json",
    ]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_injections_agenodes.py",
        "status": "SUCCESS",
        "issue": "#1c step 3c — truth-age node injections for the joint age–k fit",
        "upstream_repair_version": REPAIR_VERSION,
        "output_version": version,
        "anchor_prior_mode": ANCHOR_PRIOR_MODE,
        "upstream_unchanged": (
            "repair_v4 is a WP5-only version: WP3 extinction and WP4 ages and "
            "masses are consumed from repair_v3 unchanged"
        ),
        "seed": w.SEED,
        "truth_binary_fraction_model": (
            "constant F_BINARY at every mass (repair_v1..v6)"
            if args.fbin_model == "constant"
            else "mass-dependent, extended below and above 8 Msun "
                 "(repair_v7); see provenance/wp5_fbin_discriminator_prereg.json"
        ),
        "seed_recipe": (
            "fresh default_rng(SEED) per node, so every node of a branch shares "
            "its donor, binary and extinction realization and differs only in "
            "truth age"
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "truth_age_interpolation": interpolate,
        "node_rule": (
            (
                "wp4_repair_common.age_posterior_nodes kept unsnapped, each "
                "with weight 1/N; the truth isochrone is interpolated between "
                "native ages by wp5_common.load_isochrone_between_ages, the "
                "same bracketing and linear blend the recovery side already "
                "uses.  Issue #13: the rule is Lipschitz-1 in the WP4 "
                "posterior instead of discontinuous.  Identical rule for every "
                "subgroup and branch; no new parameter."
            )
            if interpolate
            else (
                "wp4_repair_common.age_posterior_nodes snapped to native "
                "isochrone ages; prior weight = summed node count.  Identical "
                "rule for every subgroup and branch; no new parameter."
            )
        ),
        "qmc_validation": validation,
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in inputs},
        "node_plan": plan,
        "generated": generated,
        "frozen_wp5_outputs_overwritten": False,
    }
    w.write_json(
        w.PROVENANCE / f"wp5_injections_agenodes_execution_{version}.json", record
    )
    print(f"\ngenerated {len(generated)} node responses")
    print(f"wrote provenance/wp5_injections_agenodes_execution_{version}.json")


if __name__ == "__main__":
    main()
