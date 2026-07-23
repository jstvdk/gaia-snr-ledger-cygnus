#!/usr/bin/env python3
"""WP4 age fit over the full branch grid.

For every combination of
    subgroup   in {CygOB2-A, B, C}
    family     in {PARSEC, MIST}                  (Class E - never averaged)
    R_V        in {3.0, 3.1, 3.5}                 (Class E - via de-reddened data)
    f_bin      in {0.3, 0.4, 0.5}                 (Class D/E - binary fraction)
    indicator  in {ums (upper main seq), pms (turn-on/lower CMD)}
we fit the de-reddened CMD by the forward-model likelihood of wp4_common,
weighting each star by its WP2 membership probability, and record the age
posterior (MAP + 68/90% credible intervals).

The distance is fixed at 1.6245 kpc; its +/-0.045 kpc (=> Delta mu = +/-0.060
mag) is propagated by refitting at mu +/- sigma_mu for the baseline
(R_V=3.1, f_bin=0.4) and recorded as the distance systematic.

Outputs:
    data/processed/wp4_age_posteriors.parquet   (one row per branch combination)
    data/processed/wp4_posterior_curves.npz     (fine posterior curves, baseline)
    provenance/wp4_fit_ages_execution.json
"""
from __future__ import annotations

import json
import hashlib
import datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w

BASE_RV = 3.1
BASE_FBIN = 0.4


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run():
    m = w.load_members()
    iso = w.load_isochrones()

    rows = []
    curves = {}   # keyed for the baseline branch, for the figures

    for rv in w.R_V_BRANCHES:
        fit = w.branch_photometry(m, rv)
        for family in ["PARSEC", "MIST"]:
            isodf = iso[family]
            for f_bin in w.F_BIN_BRANCHES:
                for sub in w.SUBGROUPS:
                    fsub = fit[fit.subgroup == sub]
                    for ind in ["ums", "pms"]:
                        # distance branch: baseline gets +/- shifts too
                        dmus = [0.0]
                        if abs(rv - BASE_RV) < 1e-9 and abs(f_bin - BASE_FBIN) < 1e-9:
                            dmus = [0.0, +w.DIST_MODULUS_ERR, -w.DIST_MODULUS_ERR]
                        for dmu in dmus:
                            ages, ll, n = w.fit_age_grid(
                                fsub, isodf, family, f_bin, ind, dmu=dmu)
                            post = w.posterior_from_loglike(ages, ll)
                            measurable = (n >= w.N_MIN_INDICATOR) and np.isfinite(post["map"])
                            rows.append(dict(
                                subgroup=sub, family=family, R_V=rv, f_bin=f_bin,
                                indicator=ind, dmu=round(dmu, 4), n_stars=n,
                                measurable=bool(measurable),
                                age_map=post["map"], age_lo68=post["lo68"],
                                age_hi68=post["hi68"], age_lo90=post["lo90"],
                                age_hi90=post["hi90"], age_mean=post["mean"],
                            ))
                            # keep fine curves only for the baseline, dmu=0
                            if (abs(rv - BASE_RV) < 1e-9 and abs(f_bin - BASE_FBIN) < 1e-9
                                    and dmu == 0.0 and post["ages_fine"] is not None):
                                key = f"{sub}|{family}|{ind}"
                                curves[key + "|age"] = post["ages_fine"]
                                curves[key + "|post"] = post["post"]
                                curves[key + "|grid_age"] = ages
                                curves[key + "|grid_ll"] = ll

    df = pd.DataFrame(rows)
    out_parquet = w.PROC / "wp4_age_posteriors.parquet"
    df.to_parquet(out_parquet, index=False)

    out_npz = w.PROC / "wp4_posterior_curves.npz"
    np.savez_compressed(out_npz, **curves)

    # ---- console summary of the baseline (R_V=3.1, f_bin=0.4, dmu=0) ----
    base = df[(df.R_V == BASE_RV) & (df.f_bin == BASE_FBIN) & (df.dmu == 0.0)]
    print("=== Baseline ages (R_V=3.1, f_bin=0.4) ===")
    for _, r in base.sort_values(["indicator", "subgroup", "family"]).iterrows():
        flag = "" if r.measurable else "  [NOT MEASURABLE]"
        print(f"  {r.indicator:3s} {r.subgroup} {r.family:6s} "
              f"MAP={r.age_map:5.2f} Myr  68%=[{r.age_lo68:.2f},{r.age_hi68:.2f}] "
              f"n={r.n_stars}{flag}")

    # ---- provenance execution log ----
    exec_log = {
        "script": "scripts/wp4_fit_ages.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "wp3_extinction": sha256(w.PROC / "wp3_extinction.parquet"),
            "wp2_subgroup_labels": sha256(w.TABLES / "wp2_subgroup_labels.parquet"),
            "isochrones_parsec": sha256(w.PROC / "wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha256(w.PROC / "wp3_isochrones_mist.parquet"),
        },
        "method": ("forward-model CMD likelihood; per-age synthetic population "
                   "(IMF-weighted single + unresolved-binary q~U[0.1,1]); "
                   "membership-P-weighted total loglike; posterior with uniform "
                   "prior in log-age over 1-10 Myr grid"),
        "config": {
            "IMF_SLOPE": w.IMF_SLOPE, "F_BIN_BRANCHES": w.F_BIN_BRANCHES,
            "Q_MIN": w.Q_MIN, "SIGMA_INT": w.SIGMA_INT, "MAG_FLOOR": w.MAG_FLOOR,
            "UMS_FAINT_EDGE": w.UMS_FAINT_EDGE, "PMS_BRIGHT_EDGE": w.PMS_BRIGHT_EDGE,
            "N_MIN_INDICATOR": w.N_MIN_INDICATOR,
            "R_V_BRANCHES": w.R_V_BRANCHES,
            "distance_kpc": w.D_KPC, "distance_kpc_err": w.D_KPC_ERR,
            "dist_modulus": w.DIST_MODULUS, "dist_modulus_err": w.DIST_MODULUS_ERR,
        },
        "n_branch_rows": len(df),
        "outputs": {
            "data/processed/wp4_age_posteriors.parquet": None,   # filled below
            "data/processed/wp4_posterior_curves.npz": None,
        },
    }
    exec_log["outputs"]["data/processed/wp4_age_posteriors.parquet"] = sha256(out_parquet)
    exec_log["outputs"]["data/processed/wp4_posterior_curves.npz"] = sha256(out_npz)
    with open(w.ROOT / "provenance" / "wp4_fit_ages_execution.json", "w") as f:
        json.dump(exec_log, f, indent=2)
    print(f"\nWrote {out_parquet.name} ({len(df)} rows), {out_npz.name}, execution log.")
    return df


if __name__ == "__main__":
    run()
