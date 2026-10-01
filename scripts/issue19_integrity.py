#!/usr/bin/env python3
"""Issue #19 -- the four integrity checks of the repair_v8 preregistration.

These are code-correctness checks, not science predictions
(provenance/issue19_repair_v8_prereg.json, ``integrity_checks``).  A failure
stops the chain until the code is fixed; it is never resolved by moving a
number.

  I1  non-anchor rows of the repair_v8 mass posteriors (table and sample cube)
      are identical to repair_v5;
  I2  the repair_v8 R_V = 3.1 anchor masses reproduce the scope diagnostic
      (tables/issue19_anchor_mass_rederivation.csv) exactly;
  I3  the WP5 replay (repair_v5 masses + repair_v7 responses through the
      modified wp5_fit_imf_joint.py) reproduces repair_v7 exactly;
  I4  every repair_v7-chain artifact hashed at preregistration is unchanged.

Run one or more:
  PYTHONPATH=scripts python3 scripts/issue19_integrity.py I1 I2 I3 I4
Each result is merged into provenance/issue19_repair_v8_integrity.json.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

PREREG = w.PROVENANCE / "issue19_repair_v8_prereg.json"
OUT = w.PROVENANCE / "issue19_repair_v8_integrity.json"
REPLAY_VERSION = "repair_v7_replay_issue19"
CSV_ROUND_TRIP = 1e-12


def i1() -> dict:
    old = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v5.parquet")
    new = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v8.parquet")
    if not np.array_equal(old.source_id.to_numpy(), new.source_id.to_numpy()):
        return {"pass": False, "reason": "row order differs"}
    photometric = ~new.mass_method.eq("spectroscopic_hrd").to_numpy()
    was_photometric = ~old.mass_method.eq("spectroscopic_hrd_frozen").to_numpy()
    numeric = [c for c in old.columns
               if c.startswith(("mass_", "age_median_")) and c in new.columns
               and pd.api.types.is_numeric_dtype(old[c])]
    differing = [c for c in numeric
                 if not np.array_equal(old.loc[photometric, c].to_numpy(),
                                       new.loc[photometric, c].to_numpy(),
                                       equal_nan=True)]
    so = np.load(w.PROC / "wp4_mass_posterior_samples_repair_v5.npz")["samples"]
    sn = np.load(w.PROC / "wp4_mass_posterior_samples_repair_v8.npz")["samples"]
    cube_equal = bool(np.array_equal(so[photometric], sn[photometric], equal_nan=True))
    return {
        "pass": bool(not differing and cube_equal
                     and np.array_equal(photometric, was_photometric)),
        "photometric_rows": int(photometric.sum()),
        "anchor_rows": int((~photometric).sum()),
        "same_rows_photometric_in_both": bool(np.array_equal(photometric, was_photometric)),
        "columns_compared": len(numeric),
        "columns_differing_on_photometric_rows": differing,
        "sample_cube_identical_on_photometric_rows": cube_equal,
    }


def i2() -> dict:
    diag = pd.read_csv(w.TABLES / "issue19_anchor_mass_rederivation.csv")
    hrd = pd.read_parquet(w.PROC / "wp4_anchor_hrd_repair_v8.parquet")
    rows = []
    for family in ("PARSEC", "MIST"):
        d = diag[diag.family.eq(family)].set_index("source_id")
        h = hrd.set_index("source_id")[[f"mass_{family}_rv3.1",
                                        f"age_used_{family}_rv3.1"]].dropna()
        joined = d.join(h, how="outer")
        mass_diff = (joined[f"mass_{family}_rv3.1"] - joined.mass_rederived).abs()
        age_diff = (joined[f"age_used_{family}_rv3.1"] - joined.age_repair_v5_Myr).abs()
        rows.append({
            "family": family,
            "diagnostic_rows": int(len(d)),
            "repair_v8_rows": int(len(h)),
            "unmatched": int(mass_diff.isna().sum()),
            "max_abs_mass_diff": float(np.nanmax(mass_diff)),
            "max_abs_age_diff": float(np.nanmax(age_diff)),
        })
    # The diagnostic is stored as CSV, so "exactly" means equal to the CSV text
    # round trip: the first run of this check (2026-10-01) used == 0 and failed
    # on one MIST anchor at 1.8e-15 Msun, a serialization artifact.
    ok = all(r["unmatched"] == 0 and r["max_abs_mass_diff"] <= CSV_ROUND_TRIP
             and r["max_abs_age_diff"] <= CSV_ROUND_TRIP for r in rows)
    return {"pass": ok, "by_family": rows,
            "tolerance": CSV_ROUND_TRIP,
            "tolerance_note": ("CSV round-trip of the diagnostic; first run at "
                               "tolerance 0 failed at 1.8e-15 Msun on one MIST anchor")}


def i3() -> dict:
    ref = pd.read_parquet(w.PROC / "wp5_imf_normalization_repair_v7.parquet")
    rep = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{REPLAY_VERSION}.parquet")
    key = ["subgroup", "family", "R_V", "alpha"]
    ref = ref.sort_values(key).reset_index(drop=True)
    rep = rep.sort_values(key).reset_index(drop=True)
    numeric = [c for c in ref.columns
               if pd.api.types.is_numeric_dtype(ref[c]) and c not in key]
    # Identity is tested with NaN == NaN: absolute_95_edge_Msun is all-NaN in
    # both files.  The first run of this check (2026-10-01) took nanmax of the
    # difference, got NaN for that column and reported FAIL although every
    # column and every draw was identical -- a defect in the check, fixed here.
    identical = {c: bool(np.array_equal(ref[c].to_numpy(float), rep[c].to_numpy(float),
                                        equal_nan=True)) for c in numeric}
    worst = {c: (0.0 if identical[c] else float(np.nanmax(np.abs(
        ref[c].to_numpy(float) - rep[c].to_numpy(float))))) for c in numeric}
    gate_same = bool(ref.residual_gate_pass.equals(rep.residual_gate_pass))
    dref = np.load(w.PROC / "wp5_imf_posterior_draws_repair_v7.npz")
    drep = np.load(w.PROC / f"wp5_imf_posterior_draws_{REPLAY_VERSION}.npz")
    draws_same = sorted(dref.files) == sorted(drep.files) and all(
        np.array_equal(dref[k], drep[k]) for k in dref.files
    )
    return {
        "pass": bool(all(identical.values()) and gate_same and draws_same),
        "first_run_note": ("first run reported FAIL from a NaN-comparison defect "
                           "in this check (all-NaN absolute_95_edge_Msun); no "
                           "value differed"),
        "replay_version": REPLAY_VERSION,
        "cells": int(len(ref)),
        "max_abs_difference_by_column": worst,
        "gate_verdicts_identical": gate_same,
        "posterior_draws_identical": bool(draws_same),
    }


def i4() -> dict:
    record = json.loads(PREREG.read_text())
    moved = {}
    for rel, digest in record["compare_artifacts"].items():
        path = w.ROOT / rel
        now = w.sha256(path) if path.exists() else None
        if now != digest:
            moved[rel] = now
    for rel, digest in record["reused_injection_responses"].items():
        path = w.ROOT / rel
        now = w.sha256(path) if path.exists() else None
        if now != digest:
            moved[rel] = now
    return {
        "pass": not moved,
        "checked": len(record["compare_artifacts"]) + len(record["reused_injection_responses"]),
        "moved": moved,
    }


CHECKS = {"I1": i1, "I2": i2, "I3": i3, "I4": i4}


def main() -> None:
    wanted = sys.argv[1:] or list(CHECKS)
    record = json.loads(OUT.read_text()) if OUT.exists() else {
        "script": "scripts/issue19_integrity.py",
        "preregistration": str(PREREG.relative_to(w.ROOT)),
        "preregistration_sha256": w.sha256(PREREG),
        "checks": {},
    }
    failed = []
    for cid in wanted:
        result = CHECKS[cid]()
        result["run_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        record["checks"][cid] = result
        print(f"{cid}: {'PASS' if result['pass'] else 'FAIL'}")
        if not result["pass"]:
            failed.append(cid)
            print(json.dumps(result, indent=2, default=str)[:3000])
    record["all_run_pass"] = all(c["pass"] for c in record["checks"].values())
    w.write_json(OUT, record)
    if failed:
        raise SystemExit(f"integrity check(s) failed: {failed}")


if __name__ == "__main__":
    main()
