#!/usr/bin/env python3
"""repair_v9 Stage 1 -- real-data fits of the WP4 age redesign
(provenance/wp4v9_age_prereg.json).  Run only after wp4v9_injections.py.

Per R_V x family x f_bin x subgroup:
  M1 in every window (bright = primary, bright_m0.5, bright_m1.5, faint, full);
  M1 bright on the non-anchor stars only (for the joint fallback);
  the old WP4 model (bright, faint, full) and the correlated Gaussian (bright, full);
at the decision cell also the distance refits (dmu +/- 0.060, M1 bright, all
stars and non-anchors), and at f_bin 0.4 the x2 widened-posterior sensitivity;
the faint-end study M4 at the decision cell.

Outputs:
  tables/wp4v9_real_fits.csv
  tables/wp4v9_faint_end_study.csv
  provenance/wp4v9_fit_execution.json    (includes every fit's lnL curve)

Run:
  PYTHONPATH=scripts python3 scripts/wp4v9_fit.py [--workers 14]
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

OUT = w.TABLES / "wp4v9_real_fits.csv"
OUT_FAINT = w.TABLES / "wp4v9_faint_end_study.csv"
OUT_EXEC = w.PROVENANCE / "wp4v9_fit_execution.json"
DMU = 0.060
OUTLIER_EPS, OUTLIER_DENSITY = 0.05, 0.2

_S = {}


def _init():
    ext = pd.read_parquet(V.frozen("extinction_v5"))
    labels = pd.read_parquet(V.frozen("subgroup_labels"))
    npz = np.load(V.frozen("av_posterior_v5"))
    _S["stars"] = {rv: V.load_stars(ext, labels, npz, rv) for rv in V.R_V_BRANCHES}
    _S["iso"] = {f: pd.read_parquet(V.frozen(f"isochrones_{f.lower()}")) for f in V.FAMILIES}
    _S["parts"] = {}


def _parts(fam, f_bin, drop_pms=False):
    key = (fam, f_bin, drop_pms)
    if key not in _S["parts"]:
        _S["parts"][key] = V.particles_by_age(_S["iso"][fam], f_bin, drop_pms=drop_pms)
    return _S["parts"][key]


def _sub(rv, sg, nonanchor=False, p_min=None, widen=False):
    st = _S["stars"][rv]
    s = st[st.subgroup.eq(sg)]
    if nonanchor:
        s = s[~s.is_anchor]
    if p_min is not None:
        s = s[s.P > p_min]
    s = s.copy()
    s.attrs.update(st.attrs)
    return V.widened(s) if widen else s


def _fit(stars, parts, window, model, dmu=0.0, eps=0.0):
    ages, plist = parts
    s = V.in_window(stars, window) if model != "old" else stars
    if model == "old":
        rv = stars.attrs["R_V"]
        kG = stars.attrs["kG"]
        mg_wp3 = stars["g_obs"] - kG * stars[f"av_rv{rv:.1f}"]
        lo, hi = V.WINDOWS[window]
        s = stars[(mg_wp3 > lo) & (mg_wp3 <= hi)]
    s = s.copy()
    s.attrs.update(stars.attrs)
    if model == "M1":
        tot = np.array([s.P.values @ V.m1_loglike(s, p, window, dmu, eps, OUTLIER_DENSITY) for p in plist])
    elif model == "old":
        tot = np.array([s.P.values @ V.gauss_loglike(s, p, window, False, dmu,
                                                     av_col=f"av_rv{rv:.1f}", sig_col=f"av_err_rv{rv:.1f}")
                        for p in plist])
    else:
        tot = np.array([s.P.values @ V.gauss_loglike(s, p, window, True, dmu) for p in plist])
    return V.posterior(ages, tot, len(s)), tot.tolist(), len(s)


def task(cell):
    rv, fam, f_bin, sg = cell
    parts = _parts(fam, f_bin)
    rows = []

    def add(label, stars, window, model, dmu=0.0, eps=0.0, prts=None, extra=None):
        p, curve, n = _fit(stars, prts or parts, window, model, dmu, eps)
        rows.append({"subgroup": sg, "family": fam, "R_V": rv, "f_bin": f_bin, "dmu": dmu,
                     "fit": label, "window": window, "error_model": model, "n_stars": n,
                     **p, **(extra or {}), "_curve": curve})

    allst = _sub(rv, sg)
    for window in V.WINDOWS:
        add("M1", allst, window, "M1")
    add("M1_nonanchor", _sub(rv, sg, nonanchor=True), "bright", "M1")
    for window in ("bright", "faint", "full"):
        add("old", allst, window, "old")
    for window in ("bright", "full"):
        add("corr", allst, window, "corr")
    if f_bin == 0.4:
        add("M1_widened_x2", _sub(rv, sg, widen=True), "bright", "M1")
    if rv == 3.1 and f_bin == 0.4:
        for d in (+DMU, -DMU):
            add("M1", allst, "bright", "M1", dmu=d)
            add("M1_nonanchor", _sub(rv, sg, nonanchor=True), "bright", "M1", dmu=d)
        # M4 faint-end study
        add("M4_base", allst, "faint", "M1")
        if fam == "MIST":
            add("M4_i_no_pms", allst, "faint", "M1", prts=_parts(fam, f_bin, drop_pms=True))
        for fb in (0.6, 0.7):
            add(f"M4_ii_fbin{fb}", allst, "faint", "M1", prts=_parts(fam, fb))
        add("M4_iii_P_gt_0.9", _sub(rv, sg, p_min=0.9), "faint", "M1")
        add("M4_iv_outliers", allst, "faint", "M1", eps=OUTLIER_EPS)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    args = ap.parse_args()
    if not (w.TABLES / "wp4v9_injection_validation.csv").exists():
        raise SystemExit("run wp4v9_injections.py first (brief 3.3: validation before real data)")
    t0 = time.time()
    method = V.check_method_unchanged()
    cells = [(rv, fam, fb, sg) for rv in V.R_V_BRANCHES for fam in V.FAMILIES
             for fb in V.F_BINS for sg in V.SUBGROUPS]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init) as ex:
        for k, out in enumerate(ex.map(task, cells, chunksize=1)):
            rows += out
            print(f"{k + 1}/{len(cells)} cells, {time.time() - t0:.0f} s", flush=True)
    curves = {"|".join(str(r[k]) for k in ("fit", "window", "error_model", "subgroup", "family",
                                           "R_V", "f_bin", "dmu")): r.pop("_curve") for r in rows}
    table = pd.DataFrame(rows)
    table[~table.fit.str.startswith("M4")].to_csv(OUT, index=False)
    table[table.fit.str.startswith("M4")].to_csv(OUT_FAINT, index=False)
    ages = {f: V.particles_by_age(pd.read_parquet(V.frozen(f"isochrones_{f.lower()}")), 0.4)[0].tolist()
            for f in V.FAMILIES}
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp4v9_fit.py",
        "preregistration": {"path": "provenance/wp4v9_age_prereg.json", "sha256": w.sha256(V.PREREG_PATH)},
        "method_check": method,
        "inputs": V.prereg()["consumed_inputs"],
        "native_ages": ages,
        "loglike_curves": curves,
        "runtime_s": round(time.time() - t0, 1),
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (OUT, OUT_FAINT)},
    })
    print(f"done, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
