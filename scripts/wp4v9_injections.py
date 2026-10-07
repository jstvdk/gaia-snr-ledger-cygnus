#!/usr/bin/env python3
"""repair_v9 Stage 1 -- injection validation (provenance/wp4v9_age_prereg.json, "injections").

Synthetic subgroups at 2.5 / 3.2 / 4.0 / 5.0 Myr (native-snapped), R_V 3.1,
f_bin 0.4, built from the real stars' photometric errors and A_V posteriors
(wp4v9_common.synthetic), fitted with M1 and with the old WP4 model.
100 realisations for M1 in the bright window (GA1); 20 for everything else.

Outputs:
  tables/wp4v9_injection_validation.csv
  provenance/wp4v9_injection_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp4v9_injections.py [--workers 14]
"""
from __future__ import annotations

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp4v9_common as V
import wp5_common as w

OUT = w.TABLES / "wp4v9_injection_validation.csv"
OUT_EXEC = w.PROVENANCE / "wp4v9_injection_execution.json"
SEED = 20261006
TRUE_AGES = (2.5, 3.2, 4.0, 5.0)
R_PRIMARY, R_OTHER = 100, 20
RV, F_BIN = 3.1, 0.4

_STATE = {}


def _init():
    ext = pd.read_parquet(V.frozen("extinction_v5"))
    labels = pd.read_parquet(V.frozen("subgroup_labels"))
    npz = np.load(V.frozen("av_posterior_v5"))
    stars = V.load_stars(ext, labels, npz, RV)
    _STATE["stars"] = {sg: stars[stars.subgroup.eq(sg)].copy() for sg in V.SUBGROUPS}
    for sg in V.SUBGROUPS:
        _STATE["stars"][sg].attrs.update(stars.attrs)
    _STATE["parts"] = {}
    for fam in V.FAMILIES:
        iso = pd.read_parquet(V.frozen(f"isochrones_{fam.lower()}"))
        _STATE["parts"][fam] = V.particles_by_age(iso, F_BIN)


def _fit(stars, parts, window, model):
    ages, plist = parts
    if model == "M1":
        tot = np.array([stars.P.values @ V.m1_loglike(stars, p, window) for p in plist])
    else:
        tot = np.array([stars.P.values @ V.gauss_loglike(stars, p, window, correlated=False,
                                                         av_col=f"av_rv{RV:.1f}",
                                                         sig_col=f"av_err_rv{RV:.1f}") for p in plist])
    return V.posterior(ages, tot, len(stars))


def task(args):
    fam, sg, true_age, r, windows = args
    stars_all = _STATE["stars"][sg]
    ages, plist = _STATE["parts"][fam]
    j = int(np.argmin(np.abs(ages - true_age)))
    rng = np.random.default_rng([SEED, V.FAMILIES.index(fam), V.SUBGROUPS.index(sg), j, r])
    rows = []
    for window, models in windows:
        n_target = len(V.in_window(stars_all, window))
        syn = V.synthetic(stars_all, plist[j], window, n_target, rng)
        for model in models:
            p = _fit(syn, (ages, plist), window, model)
            rows.append({"family": fam, "subgroup": sg, "true_age": float(ages[j]),
                         "realisation": r, "window": window, "model": model,
                         "n_stars": n_target, **p,
                         "covered68": bool(p["age_lo68"] <= ages[j] <= p["age_hi68"]),
                         "bias_median": p["age_median"] - float(ages[j])})
    return rows


PARTS = w.PROC / "wp4v9_injection_parts"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    args = ap.parse_args()
    t0 = time.time()
    method = V.check_method_unchanged()
    PARTS.mkdir(parents=True, exist_ok=True)
    # one part file per family x subgroup, written when that block completes;
    # a restart skips finished blocks (same seeds, so results are identical)
    for fam in V.FAMILIES:
        for sg in V.SUBGROUPS:
            part = PARTS / f"{fam}_{sg}.csv"
            if part.exists():
                print(f"skip {fam} {sg} (done)", flush=True)
                continue
            jobs = []
            for a in TRUE_AGES:
                for r in range(R_PRIMARY):
                    windows = [("bright", ["M1"] + (["old"] if r < R_OTHER else []))]
                    if r < R_OTHER:
                        windows += [("faint", ["M1", "old"]), ("full", ["M1", "old"])]
                    jobs.append((fam, sg, a, r, windows))
            rows = []
            with ProcessPoolExecutor(max_workers=args.workers, initializer=_init) as ex:
                for k, out in enumerate(ex.map(task, jobs, chunksize=1)):
                    rows += out
                    if k % 50 == 0:
                        print(f"{fam} {sg}: {k}/{len(jobs)} jobs, {time.time() - t0:.0f} s", flush=True)
            tmp = part.with_suffix(".tmp")
            pd.DataFrame(rows).to_csv(tmp, index=False)
            tmp.rename(part)
            print(f"block {fam} {sg} written, {time.time() - t0:.0f} s", flush=True)
    parts = [PARTS / f"{fam}_{sg}.csv" for fam in V.FAMILIES for sg in V.SUBGROUPS]
    table = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    table.to_csv(OUT, index=False)
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp4v9_injections.py",
        "preregistration": {"path": "provenance/wp4v9_age_prereg.json", "sha256": w.sha256(V.PREREG_PATH)},
        "method_check": method,
        "seed": SEED, "true_ages_requested": TRUE_AGES, "R_primary": R_PRIMARY, "R_other": R_OTHER,
        "parts": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in parts},
        "note": "run in resumable blocks; a 2026-10-06 run was killed before writing anything",
        "runtime_s_this_invocation": round(time.time() - t0, 1),
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    print(f"done, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
