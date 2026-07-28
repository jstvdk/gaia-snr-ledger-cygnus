#!/usr/bin/env python3
"""Issue #15: paired injections measuring the multiplicity bias on the census.

Pre-registered in provenance/wp6_multiplicity_prereg.json.  Read that first:
the predictions M1/M2/M3 and the decision rule were written before this ran.

WHY THIS IS A PAIRED RUN
------------------------
The naive experiment -- inject with the mass-dependent f_bin and compare against
the closure ratio already on disk -- is contaminated.  The response on disk for
M >= 8 comes from two files generated on different mass grids (the frozen
0.5-18 MASS_GRID and the 20-115 extension), so a run restricted to M >= 8 draws
a different RNG stream even with the same seed.  Any shift would then mix the
multiplicity change with a change of Monte Carlo realization.

So each node is injected TWICE on the identical M >= 8 grid, from a fresh
default_rng(SEED) each time:

  control    truth_binary_fraction = None  -> the frozen constant 0.40
  treatment  truth_binary_fraction = f(M)  -> the measured mass-dependent rate

Because the per-star binary threshold consumes the same single
rng.random(n_injected) draw in both arms, every donor, extinction, photometric
and QMC realization is bit-identical between them.  The two arms differ in
exactly one thing: which stars got a companion.  The difference in closure
ratio is therefore the multiplicity effect and nothing else.

The recovery side is unchanged in BOTH arms and keeps assuming 0.40.  That is
the point: nature makes binaries at the true rate while the estimator assumes
0.40, and that mismatch IS the bias under test.

Outputs:
  data/processed/wp6_mult{ctl,trt}_{sg}_{family}_rv{rv}_age{age}_response.parquet
  provenance/wp6_multiplicity_injections_execution.json

Run (about 3 h; --dry-run to see the plan, --preflight for the V4 check only):
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_multiplicity_injections.py
"""
from __future__ import annotations

import argparse
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
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, WP6_MASS_EXTENSION, turnoff_mass
from wp6_massive_injections import response_path as extension_path
from wp6_multiplicity_prereg import truth_binary_fraction

WP5_VERSION = "repair_v6"
SN_THRESHOLD_MSUN = 8.0
ARMS = ("ctl", "trt")


def closure_grid(family: str, age: float) -> np.ndarray:
    """The closure window's mass grid: the frozen grid above 8 Msun joined to
    the WP6 extension, truncated at this node's turnoff.  Identical to the set
    of masses the accepted closure test integrates over, so the diagnostic and
    the test see the same window."""
    cap = min(turnoff_mass(family, age), IMF_UPPER_LIMIT)
    base = w.MASS_GRID[w.MASS_GRID >= SN_THRESHOLD_MSUN]
    extension = WP6_MASS_EXTENSION[WP6_MASS_EXTENSION <= cap]
    return np.unique(np.concatenate([base[base <= cap], extension]))


def arm_path(arm: str, subgroup: str, family: str, rv: float, age: float) -> Path:
    label = subgroup.replace("CygOB2-", "")
    rv_tag = f"{rv:.1f}".replace(".", "p")
    age_tag = f"{age:.3f}".replace(".", "p")
    return (
        w.PROC
        / f"wp6_mult{arm}_{label}_{family}_rv{rv_tag}_age{age_tag}_response.parquet"
    )


def preflight(context: dict) -> dict:
    """V4: the default truth_binary_fraction=None must be bit-preserving.

    Regenerates one already-published WP6 extension node with the current code
    and compares the parquet hash to the stored one.  If adding the kwarg had
    perturbed the RNG stream, this fails and nothing else runs.
    """
    plan = json.loads(
        (w.PROVENANCE / "wp6_massive_injections_execution.json").read_text(
            encoding="utf-8"
        )
    )
    entry = plan["generated"][0]
    subgroup, family = entry["subgroup"], entry["family"]
    rv, age = float(entry["R_V"]), float(entry["truth_age_Myr"])
    stored = extension_path(subgroup, family, rv, age)
    expected = w.sha256(stored)

    _, response, _ = R.inject_curve(
        subgroup, family, rv, context["classifier"], context["donor_pool"],
        context["donor_model"], context["normal_points"], np.random.default_rng(w.SEED),
        context["posterior_ids"], context["posterior_cube"], context["anchor_map"],
        context["template_magnitudes"], context["template_weights"],
        context["branch_sigma"][rv], context["age_posterior"],
        truth_age_override=age, interpolate_truth_age=True,
        mass_grid=np.array(entry["extension_masses"], dtype=float),
        truth_binary_fraction=None,
    )
    scratch = w.PROC / "wp6_multiplicity_preflight.parquet"
    response.to_parquet(scratch, index=False)
    actual = w.sha256(scratch)
    scratch.unlink()
    return {
        "node": f"{subgroup} {family} R_V={rv} age={age:.3f}",
        "artifact": str(stored.relative_to(w.ROOT)),
        "expected_sha256": expected,
        "actual_sha256": actual,
        "V4_pass": bool(actual == expected),
        "claim": (
            "adding truth_binary_fraction to inject_curve leaves the default "
            "path bit-identical"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    if (REPAIR_VERSION, ANCHOR_PRIOR_MODE) != ("repair_v5", "kriging"):
        raise RuntimeError(
            "the multiplicity diagnostic holds the accepted repair_v6 WP5 "
            f"fixed, whose upstream is repair_v5/kriging; got "
            f"{(REPAIR_VERSION, ANCHOR_PRIOR_MODE)!r}"
        )

    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet"
    )
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    plan = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                nodes = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family],
                    snap=not J.uses_age_interpolation(WP5_VERSION),
                )
                for age, weight in nodes.items():
                    masses = closure_grid(family, age)
                    plan.append(
                        {
                            "subgroup": subgroup, "family": family,
                            "R_V": float(rv), "truth_age_Myr": float(age),
                            "prior_weight": float(weight),
                            "turnoff_Msun": round(
                                min(turnoff_mass(family, age), IMF_UPPER_LIMIT), 2
                            ),
                            "n_masses": int(len(masses)),
                            "masses": [float(m) for m in masses],
                            "arms_done": [
                                arm for arm in ARMS
                                if arm_path(arm, subgroup, family, rv, age).exists()
                            ],
                        }
                    )
    to_run = [entry for entry in plan if len(entry["arms_done"]) < len(ARMS)]
    injections = sum(
        entry["n_masses"] * w.N_INJECT_PER_MASS * (len(ARMS) - len(entry["arms_done"]))
        for entry in to_run
    )
    print(
        f"multiplicity plan: {len(plan)} nodes x {len(ARMS)} arms, "
        f"{len(to_run)} nodes outstanding, {injections:,} injections"
    )
    print("  truth-side f_bin: " + ", ".join(
        f"{m:g}->{truth_binary_fraction(m):.3f}" for m in [8.0, 10.0, 12.0, 16.0, 40.0]
    ))
    if args.dry_run:
        for entry in plan[:6]:
            print(f"  {entry['family']:6s} rv{entry['R_V']} {entry['subgroup'][-1]} "
                  f"age={entry['truth_age_Myr']:.3f} "
                  f"turnoff={entry['turnoff_Msun']:.1f} "
                  f"masses={entry['n_masses']}")
        return

    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = R.build_donor_pool(classifier)
    donor_pool = R.augment_donor_pool(donor_pool)
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("multiplicity injection QMC validation failed")
    posterior = np.load(w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    repair_provenance = json.loads(
        (w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8")
    )
    context = {
        "classifier": classifier,
        "donor_pool": donor_pool,
        "donor_model": donor_model,
        "normal_points": normal_points,
        "posterior_ids": posterior["source_id"].astype("int64"),
        "posterior_cube": posterior["probability"],
        "anchor_map": AnchorMap.from_frozen_wp3(),
        "template_magnitudes": load_template_library()[1],
        "template_weights": load_template_library()[2],
        "age_posterior": age_posterior,
        "branch_sigma": {
            rv: float(
                repair_provenance["configuration"][
                    "template_branch_uncertainty_calibration"
                ][f"rv{rv:.1f}"]["adopted_template_branch_sigma_mag"]
            )
            for rv in w.R_V_BRANCHES
        },
    }

    print("\nV4 preflight — is the default path still bit-identical?", flush=True)
    v4 = preflight(context)
    print(f"  {v4['node']}")
    print(f"  V4 {'PASS' if v4['V4_pass'] else 'FAIL'} — {v4['artifact']}")
    if not v4["V4_pass"]:
        raise RuntimeError(
            "V4 FAILED: adding truth_binary_fraction perturbed the default RNG "
            "stream, so the control arm is not comparable to the published "
            "response.  Refusing to run the diagnostic."
        )
    if args.preflight:
        return

    generated = []
    for index, entry in enumerate(to_run, start=1):
        subgroup, family = entry["subgroup"], entry["family"]
        rv, age = entry["R_V"], entry["truth_age_Myr"]
        masses = np.array(entry["masses"], dtype=float)
        print(
            f"=== [{index}/{len(to_run)}] {subgroup} {family} R_V={rv} "
            f"age {age:.3f} Myr, {len(masses)} masses "
            f"to {masses.max():.0f} Msun ===",
            flush=True,
        )
        node = {**entry, "arms": {}}
        for arm in ARMS:
            out = arm_path(arm, subgroup, family, rv, age)
            if out.exists():
                node["arms"][arm] = {"reused": True}
                continue
            started = time.time()
            _, response, summary = R.inject_curve(
                subgroup, family, rv, classifier, donor_pool, donor_model,
                normal_points, np.random.default_rng(w.SEED),
                context["posterior_ids"], context["posterior_cube"],
                context["anchor_map"], context["template_magnitudes"],
                context["template_weights"], context["branch_sigma"][rv],
                age_posterior,
                truth_age_override=age,
                interpolate_truth_age=True,
                mass_grid=masses,
                truth_binary_fraction=(
                    None if arm == "ctl" else truth_binary_fraction
                ),
            )
            response.to_parquet(out, index=False)
            node["arms"][arm] = {
                "reused": False,
                "binary_fraction_realized": round(
                    float(summary["binary_fraction_realized"]), 4
                ),
                "membership_pass": int(summary["membership_pass"]),
                "mass_recovered": int(summary["mass_recovered"]),
                "elapsed_seconds": round(time.time() - started, 1),
            }
            print(
                f"    {arm}: f_bin realized "
                f"{node['arms'][arm]['binary_fraction_realized']:.3f}, "
                f"{node['arms'][arm]['mass_recovered']} masses recovered, "
                f"{node['arms'][arm]['elapsed_seconds']:.0f}s",
                flush=True,
            )
        generated.append(node)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_multiplicity_injections.py",
        "status": "SUCCESS",
        "work_package": "WP6 issue #15 diagnostic",
        "prereg": "provenance/wp6_multiplicity_prereg.json",
        "wp5_version_held_fixed": WP5_VERSION,
        "upstream_repair_version": REPAIR_VERSION,
        "anchor_prior_mode": ANCHOR_PRIOR_MODE,
        "design": (
            "paired.  Each node is injected twice on the identical M >= 8 grid "
            "from a fresh default_rng(SEED): control keeps the frozen constant "
            "F_BINARY = 0.40, treatment applies the measured mass-dependent "
            "f_bin(M).  The per-star threshold consumes the same single "
            "rng.random(n_injected) draw in both arms, so donor, extinction, "
            "photometric and QMC realizations are bit-identical and the arms "
            "differ in exactly one thing: which stars got a companion."
        ),
        "why_paired": (
            "the published M >= 8 response is assembled from two files "
            "generated on different mass grids, so a grid-restricted re-run "
            "draws a different RNG stream even at the same seed.  Comparing "
            "the treatment directly against it would mix the multiplicity "
            "change with a change of Monte Carlo realization."
        ),
        "recovery_side_unchanged": (
            "both arms recover under the unchanged estimator, which keeps "
            "assuming 0.40.  That mismatch IS the bias under test."
        ),
        "seed": w.SEED,
        "f_bin_model": {
            "anchors_Msun": [8.0, 16.0],
            "anchors_f_bin": [0.40, 0.70],
            "interpolation": "linear in log mass, flat outside",
            "grid": {
                f"{m:g}": round(float(truth_binary_fraction(m)), 4)
                for m in [8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 40.0, 100.0]
            },
        },
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__,
        },
        "qmc_validation": validation,
        "V4_preflight": v4,
        "node_plan": plan,
        "generated": generated,
        "scope_limit": (
            "this is a DIAGNOSTIC.  It measures the size of the multiplicity "
            "effect on the WP6 closure ratio while holding accepted WP5 "
            "repair_v6 fixed.  Adopting a mass-dependent f_bin would also "
            "perturb the WP5 response inside 2-8 Msun through down-scatter "
            "from above, which requires a full repair_v7 chain re-run.  Masses "
            "below 8 Msun are deliberately not re-injected here."
        ),
    }
    w.write_json(
        w.PROVENANCE / "wp6_multiplicity_injections_execution.json", record
    )
    print(f"\ngenerated {len(generated)} paired nodes")
    print("wrote provenance/wp6_multiplicity_injections_execution.json")


if __name__ == "__main__":
    main()
