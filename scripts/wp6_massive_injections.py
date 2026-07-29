#!/usr/bin/env python3
"""WP6 step 0a: inject the response above WP5's frozen 18 Msun ceiling.

Design and predicted downstream effects were recorded BEFORE this ran, in
provenance/wp6_mass_extension_decision.json.  Read that first.

Each repair_v6 truth-age node gets a companion injection on the extension grid,
truncated at that node's isochrone turnoff.  The seed recipe is the node
recipe -- fresh default_rng(SEED) per node -- so an extension shares its
parent node's donor, binary and extinction realizations.

wp5_common.MASS_GRID is NOT touched.  Outputs use their own file prefix and no
WP5 script globs it, so every accepted repair_v6 artifact is untouched; that is
verified afterwards by scripts/wp6_verify_wp5_untouched.py.

Outputs:
  data/processed/wp6_massext_{sg}_{family}_rv{rv}_age{age}_{response,curve}.parquet
  provenance/wp6_massive_injections_execution.json

Run (about 1 h):
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_massive_injections.py
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
from wp3_repair_common import ANCHOR_PRIOR_MODE, REPAIR_VERSION, AnchorMap, load_template_library
from wp5_fbin_discriminator_prereg import extended_binary_fraction
from wp6_mass_extension_decision import (
    IMF_UPPER_LIMIT,
    WP6_MASS_EXTENSION,
    turnoff_mass,
)

WP5_VERSION = "repair_v6"


def stem(subgroup: str, family: str, rv: float, age: float) -> str:
    label = subgroup.replace("CygOB2-", "")
    return (
        f"wp6_massext_{label}_{family}_rv{rv:.1f}".replace(".", "p")
        .replace(f"rv{rv:.1f}".replace(".", "p"), f"rv{rv:.1f}".replace(".", "p"))
        + f"_age{age:.3f}".replace(".", "p")
    )


def response_path(
    subgroup: str, family: str, rv: float, age: float, version: str = "repair_v6"
) -> Path:
    label = subgroup.replace("CygOB2-", "")
    rv_tag = f"{rv:.1f}".replace(".", "p")
    age_tag = f"{age:.3f}".replace(".", "p")
    # repair_v6 keeps the original unsuffixed name so every accepted artifact
    # and every hash already recorded against it stays valid.
    tag = "" if version == "repair_v6" else f"_{version}"
    return w.PROC / (
        f"wp6_massext_{label}_{family}_rv{rv_tag}_age{age_tag}{tag}_response.parquet"
    )


def curve_path(
    subgroup: str, family: str, rv: float, age: float, version: str = "repair_v6"
) -> Path:
    return Path(str(response_path(subgroup, family, rv, age, version)).replace(
        "_response.parquet", "_curve.parquet"
    ))


def extension_masses(family: str, age: float) -> np.ndarray:
    cap = min(turnoff_mass(family, age), IMF_UPPER_LIMIT)
    return WP6_MASS_EXTENSION[WP6_MASS_EXTENSION <= cap]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--wp5-version", default=WP5_VERSION)
    parser.add_argument(
        "--fbin-model", choices=["constant", "extended"], default="constant",
        help="truth-side binary fraction; 'extended' is repair_v7",
    )
    args = parser.parse_args()
    version = args.wp5_version
    fbin_model = (
        extended_binary_fraction if args.fbin_model == "extended" else None
    )
    if (REPAIR_VERSION, ANCHOR_PRIOR_MODE) != ("repair_v5", "kriging"):
        raise RuntimeError(
            "WP6 consumes the accepted repair_v6 WP5, whose upstream is "
            f"repair_v5/kriging; got {(REPAIR_VERSION, ANCHOR_PRIOR_MODE)!r}"
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
                    snap=not J.uses_age_interpolation(version),
                )
                for age, weight in nodes.items():
                    masses = extension_masses(family, age)
                    plan.append(
                        {
                            "subgroup": subgroup, "family": family, "R_V": float(rv),
                            "truth_age_Myr": float(age), "prior_weight": float(weight),
                            "turnoff_Msun": round(
                                min(turnoff_mass(family, age), IMF_UPPER_LIMIT), 2
                            ),
                            "extension_masses": [float(m) for m in masses],
                            "already_generated": response_path(
                                subgroup, family, rv, age, version
                            ).exists(),
                        }
                    )
    to_run = [entry for entry in plan if not entry["already_generated"]]
    print(f"extension plan: {len(plan)} nodes, {len(to_run)} require injection")
    if args.dry_run:
        for entry in plan[:6]:
            print(f"  {entry['family']:6s} rv{entry['R_V']} {entry['subgroup'][-1]} "
                  f"age={entry['truth_age_Myr']:.3f} turnoff={entry['turnoff_Msun']:.1f} "
                  f"masses={len(entry['extension_masses'])}")
        return

    classifier = w.reconstruct_wp2_classifier()
    donor_pool = R.augment_donor_pool(R.build_donor_pool(classifier)[0])
    donor_model = R.build_donor_pool(classifier)[1]
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("WP6 extension QMC validation failed")
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
        subgroup, family = entry["subgroup"], entry["family"]
        rv, age = entry["R_V"], entry["truth_age_Myr"]
        masses = np.array(entry["extension_masses"], dtype=float)
        started = time.time()
        print(
            f"=== [{index}/{len(to_run)}] {subgroup} {family} R_V={rv} "
            f"age {age:.3f} Myr, {len(masses)} masses to "
            f"{masses.max():.0f} Msun (turnoff {entry['turnoff_Msun']:.1f}) ===",
            flush=True,
        )
        curve, response, summary = R.inject_curve(
            subgroup, family, rv, classifier, donor_pool, donor_model,
            normal_points, np.random.default_rng(w.SEED), posterior_ids,
            posterior_cube, anchor_map, template_magnitudes, template_weights,
            branch_sigma[rv], age_posterior,
            truth_age_override=age,
            interpolate_truth_age=True,
            mass_grid=masses,
            truth_binary_fraction=fbin_model,
        )
        response.to_parquet(
            response_path(subgroup, family, rv, age, version), index=False
        )
        curve.to_parquet(
            curve_path(subgroup, family, rv, age, version), index=False
        )
        generated.append(
            {
                **entry,
                "membership_pass": int(summary["membership_pass"]),
                "mass_recovered": int(summary["mass_recovered"]),
                "elapsed_seconds": round(time.time() - started, 1),
            }
        )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_massive_injections.py",
        "status": "SUCCESS",
        "work_package": "WP6 step 0a",
        "decision_record": "provenance/wp6_mass_extension_decision.json",
        "wp5_version_consumed": version,
        "truth_binary_fraction_model": args.fbin_model,
        "upstream_repair_version": REPAIR_VERSION,
        "anchor_prior_mode": ANCHOR_PRIOR_MODE,
        "frozen_and_not_touched": (
            "wp5_common.MASS_GRID and every repair_v1..v6 artifact; the "
            "extension uses its own file prefix and its own mass grid, passed "
            "through the new inject_curve(mass_grid=...) argument whose default "
            "is the frozen grid"
        ),
        "seed": w.SEED,
        "seed_recipe": (
            "fresh default_rng(SEED) per node, matching the repair_v6 node "
            "recipe, so each extension shares its parent node's donor, binary "
            "and extinction realization"
        ),
        "extension_grid_Msun": [float(m) for m in WP6_MASS_EXTENSION],
        "imf_upper_limit_Msun": IMF_UPPER_LIMIT,
        "per_node_turnoff_cap": (
            "each node's grid is truncated at the interpolated isochrone "
            "turnoff at its truth age (issue #14); above it the PMS/MS locus "
            "does not exist and the star is dead by construction"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__,
        },
        "qmc_validation": validation,
        "node_plan": plan,
        "generated": generated,
        "verification_obligations": [
            "V1 — run scripts/wp6_verify_wp5_untouched.py",
            "V3 — no extension mass exceeds its node's turnoff (enforced here)",
        ],
    }
    name = (
        "wp6_massive_injections_execution.json" if version == "repair_v6"
        else f"wp6_massive_injections_execution_{version}.json"
    )
    w.write_json(w.PROVENANCE / name, record)
    print(f"\ngenerated {len(generated)} extension responses")
    print("wrote provenance/wp6_massive_injections_execution.json")


if __name__ == "__main__":
    main()
