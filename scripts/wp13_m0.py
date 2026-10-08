#!/usr/bin/env python3
"""WP13 -- driver for M0, the pooled one-population model (prereg A5).

M0 is the repair_v9 chain re-run with ONE label.  The science is not
re-implemented: every engine script runs unchanged, in this process (runpy),
after two in-memory patches that touch no file on disk:

  * pooled labels -- ``tables/wp2_subgroup_labels.parquet`` is read back with
    CygOB2-A/B/C replaced by wp13_rules.POOLED_LABEL, and the shared
    ``wp4_common.SUBGROUPS`` list (the same object as ``wp5_common.SUBGROUPS``)
    becomes [POOLED_LABEL].  The 61 unlabelled members have no row there, so
    they stay excluded, as in M1;
  * truth-age interpolation -- ``wp5_joint_age_fit.AGE_INTERPOLATED_VERSIONS``
    gains this driver's versions, so M0 uses repair_v9's interpolated node rule.

The injection driver forks its workers, so both patches reach them.  With
pooling off (the ``identity`` stages) only the interpolation patch applies, and
W13-I1 checks the code path reproduces repair_v9.

Stages (each refuses to overwrite its outputs):
  ages              pooled test-c age table + W13-I4    -> wp4_age_posteriors_wp13_m0.parquet
  anchors, masses   WP4 anchors and mass posteriors at the pooled age
  inject            WP5 truth-age nodes, pooled footprint (9 nodes x 6 cells)
  fit               WP5 joint fit                       -> wp5_imf_posterior_draws_wp13_m0.npz
  identity_*        the same with pooling off (W13-I1); identity_inject is PARSEC R_V 3.1
                    only and identity_fit reuses the repair_v9 node responses
  plan              dry run: the node plan only, nothing written

Run:  PYTHONPATH=scripts python3 scripts/wp13_m0.py <stage> [--workers N]
"""
from __future__ import annotations

import os

os.environ.setdefault("WP_REPAIR_VERSION", "repair_v5")      # WP3 extinction, as in repair_v9
os.environ.setdefault("WP3_ANCHOR_PRIOR_MODE", "kriging")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import runpy
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp13_rules as R
import wp4_common
import wp4v9_common as V
import wp5_common as w
import wp5_joint_age_fit as J

M0 = "wp13_m0"
IDENTITY = "wp13_identity"
PREREG = w.PROVENANCE / "wp13_prereg.json"
HEADLINE = "repair_v9_headline"
LABELLED = ("CygOB2-A", "CygOB2-B", "CygOB2-C")
AGES_OUT = w.PROC / f"wp4_age_posteriors_{M0}.parquet"
AGES_JSON = w.PROVENANCE / f"{M0}_ages_execution.json"
KEY = ["subgroup", "family", "R_V", "f_bin", "dmu"]
POOLED_FLAG = "pooled_A_B_C:wp13_m0"

_read_parquet = pd.read_parquet


def _pooled_read_parquet(path, *args, **kwargs):
    frame = _read_parquet(path, *args, **kwargs)
    if str(path).endswith("wp2_subgroup_labels.parquet"):
        frame = frame.copy()
        frame.loc[frame["subgroup"].isin(LABELLED), "subgroup"] = R.POOLED_LABEL
    return frame


def patch(pooled: bool) -> None:
    J.AGE_INTERPOLATED_VERSIONS = frozenset(J.AGE_INTERPOLATED_VERSIONS | {M0, IDENTITY})
    if pooled:
        pd.read_parquet = _pooled_read_parquet
        wp4_common.SUBGROUPS[:] = [R.POOLED_LABEL]
        if w.SUBGROUPS != [R.POOLED_LABEL]:
            raise RuntimeError("wp5_common.SUBGROUPS is not the shared wp4_common list")


def run_script(script: str, *argv: str) -> None:
    print(f"--- {script} {' '.join(argv)}", flush=True)
    saved = sys.argv
    sys.argv = [script, *argv]
    try:
        runpy.run_path(str(w.ROOT / "scripts" / script), run_name="__main__")
    finally:
        sys.argv = saved


def require_prereg() -> None:
    if not PREREG.exists():
        raise SystemExit("provenance/wp13_prereg.json missing -- pre-register and commit first")


# ----------------------------------------------------------------- ages
def stage_ages() -> None:
    """Pooled test-c curve = sum of the stored A, B and C curves (exact: the
    test-c ln L is a membership-weighted sum over stars)."""
    if AGES_OUT.exists():
        raise SystemExit(f"{AGES_OUT.relative_to(w.ROOT)} exists; nothing is overwritten")
    V.check_method_unchanged()
    curves = json.loads(V.frozen("spectroscopic_curves").read_text())["loglike_curves"]
    tests = pd.read_csv(V.frozen("spectroscopic_tests"))
    native = {f: np.array(v) for f, v in json.loads(
        (w.PROVENANCE / "wp4v9_fit_execution.json").read_text())["native_ages"].items()}
    head_record = json.loads((w.PROVENANCE / "repair_v9_headline_ages_execution.json").read_text())
    head = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{HEADLINE}.parquet")

    def spec(sg: str, fam: str, rv: float) -> tuple[np.ndarray, int]:
        row = tests[tests.test.eq("c") & np.isclose(tests.eps, R.EPS) & tests.subgroup.eq(sg)
                    & tests.family.eq(fam) & np.isclose(tests.R_V, rv)]
        if len(row) != 1:
            raise RuntimeError(f"no unique test-c row for {sg}/{fam}/{rv}")
        curve = np.array(curves[f"c|{R.EPS}|{sg}|{fam}|{rv}"], dtype=float)
        if curve.shape != native[fam].shape:
            raise RuntimeError("test-c curve is not on the native age grid")
        return curve, int(row.n_stars.iloc[0])

    rows, record, i4 = [], {}, {}
    template = head[head.subgroup.eq("CygOB2-A")]
    for (fam, rv), block in template.groupby(["family", "R_V"]):
        parts = {sg: spec(sg, fam, float(rv)) for sg in LABELLED}
        a_plus_c = parts["CygOB2-A"][0] + parts["CygOB2-C"][0]
        stored = np.array(head_record["b_combined_curves"][f"{fam}|{float(rv)}"]["lnL_curve"])
        i4[f"{fam}|{float(rv)}"] = float(np.max(np.abs(a_plus_c - stored)))
        total = sum(c for c, _ in parts.values())
        n = sum(k for _, k in parts.values())
        s = V.posterior(native[fam], total, n)
        for _, r in block.iterrows():
            out = r.to_dict()
            out.update({k: s[k] for k in ("age_map", "age_lo68", "age_hi68", "age_lo90", "age_hi90",
                                          "age_mean", "grid_railed", "measurable", "exclusion_reason")})
            out.update({"subgroup": R.POOLED_LABEL, "n_stars": n,
                        "window": "spectroscopic_hrd_A_plus_B_plus_C",
                        "error_model": "test_c_sum_A_B_C",
                        "fallback": POOLED_FLAG + ("" if r.dmu == 0.0 and r.f_bin == 0.4
                                                   else ":no_dmu_refit"),
                        "adopted": True, "branch": "", "photometry_only": False})
            rows.append(out)
        record[f"{fam}|{float(rv)}"] = {
            "n_stars": {sg: k for sg, (_, k) in parts.items()}, "age_map": s["age_map"],
            "age_lo68": s["age_lo68"], "age_hi68": s["age_hi68"], "age_median": s["age_median"],
            "bimodal": bool(V.bimodal(native[fam], total)), "lnL_curve": [float(x) for x in total]}

    table = pd.DataFrame(rows)[list(head.columns)].sort_values(KEY).reset_index(drop=True)
    if table.duplicated(KEY).any():
        raise RuntimeError("duplicate keys in the pooled age table")
    node_keys = table[np.isclose(table.f_bin, 0.4) & np.isclose(table.dmu, 0.0)]
    if len(node_keys) != len(w.FAMILIES) * len(w.R_V_BRANCHES):
        raise RuntimeError("node-rule keys are not exactly one per family x R_V")
    if not table.measurable.all():
        raise RuntimeError("a pooled age is not measurable")
    i4_pass = max(i4.values()) <= 1e-9 and bool((table.n_stars == 106).all())
    table.to_parquet(AGES_OUT, index=False)
    w.write_json(AGES_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_m0.py ages",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "rule": ("per family x R_V, the sum of the stored A, B and C test-c ln L curves (eps 0.05), "
                 "summarised with wp4v9_common.posterior; every (f_bin, dmu) key carries the same curve"),
        "W13-I4": {"max_abs_diff_A_plus_C_vs_stored_B_curve": i4, "n_stars_all_106": bool(
            (table.n_stars == 106).all()), "pass": i4_pass},
        "pooled": record,
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (
            V.frozen("spectroscopic_curves"), V.frozen("spectroscopic_tests"),
            w.PROVENANCE / "wp4v9_fit_execution.json",
            w.PROVENANCE / "repair_v9_headline_ages_execution.json",
            w.PROC / f"wp4_age_posteriors_{HEADLINE}.parquet")},
        "outputs": {str(AGES_OUT.relative_to(w.ROOT)): w.sha256(AGES_OUT)},
    })
    print(f"W13-I4 {'PASS' if i4_pass else 'FAIL'}; wrote {AGES_OUT.relative_to(w.ROOT)}")
    if not i4_pass:
        raise SystemExit("W13-I4 failed: stop and report")


# ----------------------------------------------------------------- engines
def stage_anchors(age_version: str, version: str) -> None:
    run_script("wp4_anchors_hrd.py", "--age-version", age_version, "--extinction-version", "repair_v5",
               "--output-version", version, "--preregistration", str(PREREG.relative_to(w.ROOT)))


def stage_masses(age_version: str, version: str) -> None:
    run_script("wp4_mass_posteriors_repair.py", "--age-version", age_version,
               "--anchor-version", version, "--output-version", version)


def stage_inject(age_version: str, version: str, workers: int, only: list[str] | None,
                 dry: bool = False) -> None:
    argv = ["--kind", "wp5", "--age-version", age_version, "--output-version", version,
            "--workers", str(workers)]
    if only:
        argv += ["--only", *only]
    if dry:
        argv.append("--dry-run")
    run_script("repair_v9_injections.py", *argv)


def stage_fit(mass_version: str, age_version: str, response_version: str, version: str) -> None:
    run_script("wp5_fit_imf_joint.py", "--upstream-version", "repair_v5", "--mass-version", mass_version,
               "--age-version", age_version, "--response-version", response_version,
               "--wp5-version", version, "--compare-version", "none_wp13")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["ages", "anchors", "masses", "inject", "fit", "plan",
                                      "identity_anchors", "identity_masses", "identity_inject",
                                      "identity_fit", "identity_plan"])
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()
    if args.stage not in ("plan", "identity_plan"):
        require_prereg()
    pooled = not args.stage.startswith("identity")
    if pooled and args.stage not in ("ages", "plan") and not AGES_OUT.exists():
        raise SystemExit("run the 'ages' stage first: pooled engines need the pooled age table")
    patch(pooled)
    if args.stage == "ages":
        stage_ages()
    elif args.stage == "anchors":
        stage_anchors(M0, M0)
    elif args.stage == "masses":
        stage_masses(M0, M0)
    elif args.stage == "inject":
        stage_inject(M0, M0, args.workers, None)
    elif args.stage == "fit":
        stage_fit(M0, M0, M0, M0)
    elif args.stage == "plan":
        # Before the pooled table exists this plans against the headline table:
        # wp4_repair_common.age_posterior_nodes finds no pooled row and silently
        # falls back to the median of the subgroup rows.  That shows the node
        # layout and that pooling reaches the injector; it reads no M0 number.
        # The real stages never hit the fallback (stage_ages enforces one row
        # per key, and the pooled engines refuse to run without that table).
        stage_inject(M0 if AGES_OUT.exists() else HEADLINE, M0, args.workers, None, dry=True)
    elif args.stage == "identity_anchors":
        stage_anchors(HEADLINE, IDENTITY)
    elif args.stage == "identity_masses":
        stage_masses(HEADLINE, IDENTITY)
    elif args.stage == "identity_inject":
        stage_inject(HEADLINE, IDENTITY, args.workers, ["PARSEC:3.1"])
    elif args.stage == "identity_fit":
        stage_fit(IDENTITY, HEADLINE, "repair_v9", IDENTITY)
    elif args.stage == "identity_plan":
        stage_inject(HEADLINE, IDENTITY, args.workers, ["PARSEC:3.1"], dry=True)


if __name__ == "__main__":
    main()
