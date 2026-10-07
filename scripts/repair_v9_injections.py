#!/usr/bin/env python3
"""repair_v9 -- parallel, resumable driver for the WP5 truth-age node
injections and the WP6 mass-extension injections.

The science is NOT re-implemented here.  Every node is one call of the
unchanged ``wp5_injections_repair.inject_curve`` with a fresh
``default_rng(wp5_common.SEED)`` -- exactly what ``wp5_injections_agenodes.py``
(kind wp5) and ``wp6_massive_injections.py`` (kind wp6ext) do one node after
another.  Because every node restarts the same seed, a node's output does not
depend on which process runs it or in what order, so running nodes in parallel
changes nothing.  Integrity check I1 (provenance/repair_v9_prereg.json)
verifies that on a regenerated subset of the repair_v7 nodes.

What this driver adds, and only this:
  * --age-version: the WP4 age table the node plan AND the recovery side read
    (the originals read wp4_age_posteriors_<WP_REPAIR_VERSION>, which also
    fixes the WP3 extinction and so cannot move to repair_v9);
  * a process pool (--workers);
  * atomic writes (.tmp then rename) and a per-node JSON-lines log, so a
    killed run resumes where it stopped and never leaves a half-written node;
  * it never reuses the gate-G2 scan snapshots that the snapped node rule of
    wp5_injections_agenodes.py falls back to: every node is injected.

Provenance is written in the originals' formats, under their names, so the
WP5 fit and the WP6 closure test read it unchanged:
  provenance/wp5_injections_agenodes_execution_<version>.json     (kind wp5)
  provenance/wp6_massive_injections_execution_<version>.json      (kind wp6ext)

Run (from the repository root):
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging OMP_NUM_THREADS=1 \\
  PYTHONPATH=scripts python3 scripts/repair_v9_injections.py --kind wp5 \\
      --age-version repair_v9_headline --output-version repair_v9 --workers 28
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from multiprocessing import get_context

import numpy as np
import pandas as pd
import scipy
import sklearn

import wp5_common as w
import wp5_injections_repair as R
import wp5_joint_age_fit as J
from wp3_repair_common import ANCHOR_PRIOR_MODE, REPAIR_VERSION, AnchorMap, load_template_library
from wp5_fbin_discriminator_prereg import extended_binary_fraction
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass
from wp6_massive_injections import (
    curve_path as ext_curve_path,
    extension_masses,
    response_path as ext_response_path,
)

LOGS = w.PROC / "repair_v9_logs"
_G: dict = {}          # per-worker heavy state, built once by _init


def node_paths(kind: str, entry: dict, version: str):
    sg, fam, rv, age = entry["subgroup"], entry["family"], entry["R_V"], entry["truth_age_Myr"]
    if kind == "wp5":
        return (J.node_response_path(sg, fam, rv, age, version),
                J.node_curve_path(sg, fam, rv, age, version))
    return (ext_response_path(sg, fam, rv, age, version),
            ext_curve_path(sg, fam, rv, age, version))


def build_plan(kind: str, age_posterior: pd.DataFrame, version: str, only: set | None) -> list[dict]:
    interpolate = J.uses_age_interpolation(version)
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}
    plan = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            if only and f"{family}:{rv}" not in only:
                continue
            for subgroup in w.SUBGROUPS:
                nodes = J.truth_age_nodes(age_posterior, subgroup, family, rv,
                                          native[family], snap=not interpolate)
                if len({f"{a:.3f}" for a in nodes}) != len(nodes):
                    raise RuntimeError(f"{subgroup}/{family}/R_V={rv}: nodes collide at 3 decimals")
                for age, weight in nodes.items():
                    entry = {"subgroup": subgroup, "family": family, "R_V": float(rv),
                             "truth_age_Myr": float(age), "prior_weight": float(weight)}
                    if kind == "wp6ext":
                        masses = extension_masses(family, age)
                        entry["turnoff_Msun"] = round(min(turnoff_mass(family, age), IMF_UPPER_LIMIT), 2)
                        entry["extension_masses"] = [float(m) for m in masses]
                        if not len(masses):
                            continue
                    entry["already_generated"] = all(p.exists() for p in node_paths(kind, entry, version))
                    plan.append(entry)
    return plan


def _init(age_version: str) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = R.build_donor_pool(classifier)
    donor_pool = R.augment_donor_pool(donor_pool)
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("injection QMC validation failed")
    posterior = np.load(w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    provenance = json.loads((w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8"))
    _, template_magnitudes, template_weights = load_template_library()
    _G.update(
        classifier=classifier, donor_pool=donor_pool, donor_model=donor_model,
        normal_points=normal_points, validation=validation,
        posterior_ids=posterior["source_id"].astype("int64"),
        posterior_cube=posterior["probability"], anchor_map=AnchorMap.from_frozen_wp3(),
        template_magnitudes=template_magnitudes, template_weights=template_weights,
        branch_sigma={rv: float(provenance["configuration"]["template_branch_uncertainty_calibration"]
                                [f"rv{rv:.1f}"]["adopted_template_branch_sigma_mag"])
                      for rv in w.R_V_BRANCHES},
        age_posterior=pd.read_parquet(w.PROC / f"wp4_age_posteriors_{age_version}.parquet"),
    )


def _atomic_parquet(frame: pd.DataFrame, path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    frame.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def _run_node(job: tuple) -> dict:
    kind, version, fbin, entry = job
    started = time.time()
    sg, fam, rv, age = entry["subgroup"], entry["family"], entry["R_V"], entry["truth_age_Myr"]
    interpolate = True if kind == "wp6ext" else J.uses_age_interpolation(version)
    curve, response, summary = R.inject_curve(
        sg, fam, rv, _G["classifier"], _G["donor_pool"], _G["donor_model"],
        _G["normal_points"], np.random.default_rng(w.SEED), _G["posterior_ids"],
        _G["posterior_cube"], _G["anchor_map"], _G["template_magnitudes"],
        _G["template_weights"], _G["branch_sigma"][rv], _G["age_posterior"],
        truth_age_override=age,
        interpolate_truth_age=interpolate,
        mass_grid=(np.array(entry["extension_masses"], dtype=float) if kind == "wp6ext" else None),
        truth_binary_fraction=(extended_binary_fraction if fbin == "extended" else None),
    )
    response_path, curve_path = node_paths(kind, entry, version)
    _atomic_parquet(response, response_path)
    _atomic_parquet(curve, curve_path)
    return {**entry, "already_generated": False,
            "native_isochrone_age_Myr": float(summary["age_isochrone_Myr"]),
            "membership_pass": int(summary["membership_pass"]),
            "mass_recovered": int(summary["mass_recovered"]),
            "elapsed_seconds": round(time.time() - started, 1),
            "worker_pid": os.getpid(),
            "outputs": {str(response_path.relative_to(w.ROOT)): w.sha256(response_path),
                        str(curve_path.relative_to(w.ROOT)): w.sha256(curve_path)}}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["wp5", "wp6ext"], required=True)
    ap.add_argument("--age-version", required=True)
    ap.add_argument("--output-version", required=True)
    ap.add_argument("--fbin-model", choices=["constant", "extended"], default="extended")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 4))
    ap.add_argument("--only", nargs="*", default=None,
                    help="restrict to FAMILY:R_V branches, e.g. PARSEC:3.1 (I1 subset)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    kind, version = args.kind, args.output_version
    if (REPAIR_VERSION, ANCHOR_PRIOR_MODE) != ("repair_v5", "kriging"):
        raise RuntimeError(f"expected repair_v5/kriging upstream, got {(REPAIR_VERSION, ANCHOR_PRIOR_MODE)!r}")
    record_path = (w.PROVENANCE / f"wp5_injections_agenodes_execution_{version}.json" if kind == "wp5"
                   else w.PROVENANCE / f"wp6_massive_injections_execution_{version}.json")
    if record_path.exists():
        raise SystemExit(f"{record_path.relative_to(w.ROOT)} exists; this injection set is complete")

    age_path = w.PROC / f"wp4_age_posteriors_{args.age_version}.parquet"
    age_posterior = pd.read_parquet(age_path)
    plan = build_plan(kind, age_posterior, version, set(args.only) if args.only else None)
    to_run = [e for e in plan if not e["already_generated"]]
    print(f"{kind} {version}: {len(plan)} nodes, {len(to_run)} to inject, "
          f"{args.workers} workers, ages from {age_path.name}", flush=True)
    if args.dry_run:
        for e in plan:
            print(f"  {e['family']:6s} rv{e['R_V']} {e['subgroup'][-1]} age={e['truth_age_Myr']:.3f} "
                  f"w={e['prior_weight']:.3f} {'present' if e['already_generated'] else 'INJECT'}")
        return

    LOGS.mkdir(exist_ok=True)
    log_path = LOGS / f"{kind}_{version}_nodes.jsonl"
    started = time.time()
    if to_run:
        jobs = [(kind, version, args.fbin_model, e) for e in to_run]
        ctx = get_context("fork")
        with ctx.Pool(min(args.workers, len(jobs)), initializer=_init,
                      initargs=(args.age_version,)) as pool:
            for i, rec in enumerate(pool.imap_unordered(_run_node, jobs), start=1):
                with log_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(rec) + "\n")
                print(f"  [{i}/{len(jobs)}] {rec['family']} rv{rec['R_V']} {rec['subgroup'][-1]} "
                      f"age {rec['truth_age_Myr']:.3f} {rec['elapsed_seconds']:.0f}s", flush=True)

    # every planned node must now exist; the record lists every generated node,
    # including those an interrupted earlier invocation wrote
    missing = [e for e in plan if not all(p.exists() for p in node_paths(kind, e, version))]
    if missing:
        raise RuntimeError(f"{len(missing)} planned nodes missing after the run")
    generated = {}
    if log_path.exists():
        for line in log_path.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            generated[(rec["subgroup"], rec["family"], rec["R_V"], f"{rec['truth_age_Myr']:.6f}")] = rec
    _init(args.age_version)                     # QMC validation figure for the record
    inputs = [w.PROC / "wp1_gaia_narrow.parquet", w.PROC / "wp1_2mass_join.parquet",
              w.PROC / "wp2_members.parquet", w.TABLES / "wp2_subgroup_labels.parquet",
              w.PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet",
              w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz",
              age_path, w.PROVENANCE / "wp2_membership_manifest.json",
              w.PROVENANCE / "wp3_repair_execution.json"]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/repair_v9_injections.py",
        "driver_of": ("scripts/wp5_injections_agenodes.py" if kind == "wp5"
                      else "scripts/wp6_massive_injections.py"),
        "status": "SUCCESS",
        "preregistration": "provenance/repair_v9_prereg.json",
        "upstream_repair_version": REPAIR_VERSION,
        "age_version": args.age_version,
        "output_version": version,
        "anchor_prior_mode": ANCHOR_PRIOR_MODE,
        "seed": w.SEED,
        "seed_recipe": ("fresh default_rng(SEED) per node, so every node of a branch shares its "
                        "donor, binary and extinction realization and differs only in truth age; "
                        "node results do not depend on worker or order"),
        "parallel_workers": args.workers,
        "wall_seconds": round(time.time() - started, 1),
        "truth_age_interpolation": True if kind == "wp6ext" else J.uses_age_interpolation(version),
        "truth_binary_fraction_model": args.fbin_model,
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "numpy": np.__version__, "pandas": pd.__version__,
                        "scipy": scipy.__version__, "sklearn": sklearn.__version__},
        "qmc_validation": _G["validation"],
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in inputs},
        "node_plan": plan,
        "generated": list(generated.values()),
        "frozen_outputs_overwritten": False,
    }
    if kind == "wp6ext":
        record.update({"work_package": "WP6 step 0a", "wp5_version_consumed": version,
                       "decision_record": "provenance/wp6_mass_extension_decision.json",
                       "imf_upper_limit_Msun": IMF_UPPER_LIMIT})
    w.write_json(record_path, record)
    print(f"wrote {record_path.relative_to(w.ROOT)} ({len(generated)} generated, "
          f"{time.time() - started:.0f}s)")


if __name__ == "__main__":
    main()
