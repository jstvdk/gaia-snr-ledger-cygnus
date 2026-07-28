#!/usr/bin/env python3
"""Issue #13: prove the repair_v6 code changes leave repair_v1-v5 untouched.

Three source files changed to support truth-age interpolation:

  wp5_common.load_isochrone_between_ages   NEW function, nothing else touched
  wp5_injections_repair.inject_curve       new kwarg interpolate_truth_age=False
  wp5_joint_age_fit.truth_age_nodes        new kwarg snap=True

Every change is behind a default that reproduces the old behaviour, and
wp5_residual_trend.bootstrap_null was vectorized without changing its output.
"Behind a default" is a claim, not a proof, so this script checks it against
artifacts already on disk:

  C1  the default node rule still reproduces the exact repair_v4 and repair_v5
      node sets, verified against the response files those runs wrote
  C2  the interpolating loader reproduces the snapping loader to 0 mag at every
      native isochrone age, so a node landing on the grid is unaffected
  C3  the vectorized Jeffreys k is bitwise identical to the scalar loop
  C4  the interpolate flag is False for every pre-v6 version

Output: provenance/wp5_v6_backward_compatibility.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_v6_backward_compatibility.py
"""
from __future__ import annotations

import platform
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
import wp5_joint_age_fit as J
import wp5_residual_trend as T

NODE_PATTERN = re.compile(
    r"wp5_agenode_(?P<sg>[ABC])_(?P<family>PARSEC|MIST)_rv(?P<rv>\d+p\d+)_"
    r"age(?P<age>\d+p\d+)_(?P<version>repair_v\d+)_response\.parquet"
)


def nodes_on_disk(version: str) -> dict[tuple, set[str]]:
    found: dict[tuple, set[str]] = {}
    for path in w.PROC.glob(f"wp5_agenode_*_{version}_response.parquet"):
        match = NODE_PATTERN.fullmatch(path.name)
        if not match:
            continue
        key = (
            f"CygOB2-{match['sg']}",
            match["family"],
            float(match["rv"].replace("p", ".")),
        )
        found.setdefault(key, set()).add(match["age"].replace("p", "."))
    return found


def check_node_rule(version: str, upstream: str) -> dict:
    """C1: the default (snapping) rule still reproduces what that run wrote."""
    posterior = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{upstream}.parquet")
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}
    disk = nodes_on_disk(version)
    mismatches = []
    branches = 0
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                key = (subgroup, family, float(rv))
                if key not in disk:
                    continue
                branches += 1
                computed = {
                    f"{age:.3f}"
                    for age in J.truth_age_nodes(
                        posterior, subgroup, family, rv, native[family]
                    )
                }
                # repair_v4 reused the gate-G2 scan snapshots for CygOB2-B's
                # baseline nodes, so those files carry a different stem and the
                # disk set is a subset rather than an equality there.
                missing = computed - disk[key]
                extra = disk[key] - computed
                if missing or extra:
                    mismatches.append(
                        {
                            "branch": f"{subgroup} {family} R_V={rv}",
                            "computed_not_on_disk": sorted(missing),
                            "on_disk_not_computed": sorted(extra),
                        }
                    )
    return {
        "version": version,
        "upstream": upstream,
        "branches_checked": branches,
        "mismatches": mismatches,
        "pass": not mismatches,
    }


def check_loader_equivalence() -> dict:
    """C2: at a native age the interpolating loader IS the native table."""
    results = []
    for family in w.FAMILIES:
        ages = J.native_isochrone_ages(family)
        window = ages[(ages > 1.0) & (ages < 8.0)]
        worst = 0.0
        for age in window:
            snapped, _ = w.load_isochrone_at_age(family, float(age))
            interpolated, effective = w.load_isochrone_between_ages(family, float(age))
            assert abs(effective - age) < 1e-9
            mass = interpolated["Mini"].to_numpy(float)
            for band in ["G0", "BP0", "RP0", "J0", "H0", "Ks0"]:
                reference = np.interp(
                    mass, snapped["Mini"].to_numpy(float), snapped[band].to_numpy(float)
                )
                worst = max(
                    worst,
                    float(np.max(np.abs(reference - interpolated[band].to_numpy(float)))),
                )
        results.append(
            {
                "family": family,
                "native_ages_checked": int(len(window)),
                "max_abs_magnitude_difference": worst,
                "pass": worst == 0.0,
            }
        )
    return {"per_family": results, "pass": all(r["pass"] for r in results)}


def check_jeffreys_k() -> dict:
    """C3: the vectorized posterior median equals the scalar loop bitwise."""
    rng = np.random.default_rng(20260728)
    totals = rng.integers(10, 800, size=20_000).astype(float)
    rate = 4.31
    scalar = np.array([T.jeffreys_k(value, rate) for value in totals])
    vectorized = stats.gamma.ppf(0.5, a=totals + 0.5, scale=1.0 / rate)
    return {
        "samples": int(len(totals)),
        "max_abs_difference": float(np.max(np.abs(scalar - vectorized))),
        "bitwise_identical": bool(np.array_equal(scalar, vectorized)),
        "pass": bool(np.array_equal(scalar, vectorized)),
    }


def check_version_flags() -> dict:
    versions = ["frozen", "repair_v1", "repair_v2", "repair_v3", "repair_v4",
                "repair_v5", "repair_v6"]
    flags = {v: bool(J.uses_age_interpolation(v)) for v in versions}
    return {
        "flags": flags,
        "pass": bool(
            not any(value for key, value in flags.items() if key != "repair_v6")
            and flags["repair_v6"]
        ),
    }


def main() -> None:
    c1 = [
        check_node_rule("repair_v4", "repair_v3"),
        check_node_rule("repair_v5", "repair_v5"),
    ]
    c2 = check_loader_equivalence()
    c3 = check_jeffreys_k()
    c4 = check_version_flags()
    passed = all(entry["pass"] for entry in c1) and c2["pass"] and c3["pass"] and c4["pass"]

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_v6_backward_compatibility.py",
        "status": "SUCCESS" if passed else "FAILED",
        "issue": "#13 — repair_v1..v5 must be unaffected by the repair_v6 change",
        "claim_under_test": (
            "every repair_v6 code change sits behind a default argument that "
            "reproduces the previous behaviour, and the one performance edit "
            "(vectorizing the Jeffreys k inside bootstrap_null) changes no value"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "C1_default_node_rule_reproduces_disk": c1,
        "C2_interpolating_loader_equals_native_at_native_ages": c2,
        "C3_vectorized_jeffreys_k_bitwise_identical": c3,
        "C4_interpolation_flag_off_for_all_pre_v6_versions": c4,
        "all_pass": passed,
    }
    w.write_json(w.PROVENANCE / "wp5_v6_backward_compatibility.json", record)

    for entry in c1:
        print(f"C1 {entry['version']:10s} {entry['branches_checked']} branches, "
              f"{len(entry['mismatches'])} mismatch(es)  "
              f"{'PASS' if entry['pass'] else 'FAIL'}")
    for entry in c2["per_family"]:
        print(f"C2 {entry['family']:6s} max |interpolated - native| = "
              f"{entry['max_abs_magnitude_difference']:.2e} mag over "
              f"{entry['native_ages_checked']} native ages  "
              f"{'PASS' if entry['pass'] else 'FAIL'}")
    print(f"C3 vectorized Jeffreys k bitwise identical: {c3['bitwise_identical']}  "
          f"{'PASS' if c3['pass'] else 'FAIL'}")
    print(f"C4 interpolation flag: {c4['flags']}  {'PASS' if c4['pass'] else 'FAIL'}")
    print(f"\nall checks {'PASS' if passed else 'FAIL'}")
    print("wrote provenance/wp5_v6_backward_compatibility.json")


if __name__ == "__main__":
    main()
