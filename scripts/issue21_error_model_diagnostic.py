#!/usr/bin/env python3
"""Issue #21 (found 2026-10-06 in issue #20 Phase A') -- the WP4 upper-MS age
depends on how the per-star A_V uncertainty is propagated.  Post-hoc diagnostic.

wp4_common.branch_photometry adds (k_G av_err)^2 to sigma(M_G0)^2 and
((k_BP - k_RP) av_err)^2 to sigma(colour)^2 as INDEPENDENT errors.  An A_V
error moves a star along the reddening vector, so the two are fully
correlated.  With repair_v5 av_err of ~0.6-0.9 mag for photometric stars the
difference is large.  This script refits each subgroup's upper-MS age on
exactly WP4's star sample and extinction with WP4's own particles
(wp4_common.build_model_particles), changing ONLY the error model, and repeats
it in a bright window (M_G0 <= -1.0) where the turnoff information lives.

Nothing is adopted from it.  It is not pre-registered: it diagnoses why Phase
A' test a and its control disagree with WP4.

Outputs:
  tables/issue21_error_model_diagnostic.csv
  provenance/issue21_error_model_diagnostic_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue21_error_model_diagnostic.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.special import logsumexp

import issue20_common as I20
import wp4_common as W4
import wp5_common as w
from wp3_extinction_law import band_coefficients

OUT = w.TABLES / "issue21_error_model_diagnostic.csv"
OUT_EXEC = w.PROVENANCE / "issue21_error_model_diagnostic_execution.json"
INPUTS = ["data/processed/wp3_extinction_repair_v5.parquet", "tables/wp2_subgroup_labels.parquet",
          "data/processed/wp3_isochrones_parsec.parquet", "data/processed/wp3_isochrones_mist.parquet",
          "data/processed/wp4_age_posteriors_repair_v5.parquet"]
RV, F_BIN = 3.1, 0.4


def loglike(stars, colour, mg, weight, edge, correlated, kg, kc):
    """WP4's per-star mixture likelihood, window-renormalised, with the A_V term
    either independent (WP4) or fully correlated along the reddening vector."""
    inside = mg <= edge
    cm, gm = colour[inside], mg[inside]
    lw = np.log(weight[inside] / weight[inside].sum())
    phot_g = np.sqrt(stars.sig_g0 ** 2 + W4.SIGMA_INT ** 2).to_numpy()
    phot_c = np.sqrt(stars.sig_c0 ** 2 + W4.SIGMA_INT ** 2).to_numpy()
    sa = stars.av_err.to_numpy()
    out = np.empty(len(stars))
    for i, (g, c) in enumerate(zip(stars.MG0.to_numpy(), stars.colour.to_numpy())):
        c11 = phot_g[i] ** 2 + (kg * sa[i]) ** 2
        c22 = phot_c[i] ** 2 + (kc * sa[i]) ** 2
        c12 = kg * kc * sa[i] ** 2 if correlated else 0.0
        det = c11 * c22 - c12 ** 2
        rg, rc = g - gm, c - cm
        q = (c22 * rg ** 2 - 2 * c12 * rg * rc + c11 * rc ** 2) / det
        out[i] = logsumexp(lw - 0.5 * q) - 0.5 * np.log(4 * np.pi ** 2 * det)
    return out


def main() -> None:
    ext = pd.read_parquet(w.ROOT / INPUTS[0]).drop(columns=["subgroup"], errors="ignore")
    labels = pd.read_parquet(w.ROOT / INPUTS[1])
    m = ext.merge(labels, on="source_id")
    stored = pd.read_parquet(w.ROOT / INPUTS[4])
    k = band_coefficients(RV)
    kg, kc = k["G"], k["BP"] - k["RP"]
    rows = []
    for fam in w.FAMILIES:
        iso = pd.read_parquet(w.ROOT / INPUTS[2 if fam == "PARSEC" else 3])
        ages = np.sort(iso.age_Myr.unique())
        parts = [W4.build_model_particles(iso[np.isclose(iso.age_Myr, a)], F_BIN) for a in ages]
        for sg in w.SUBGROUPS:
            sub = m[m.subgroup.eq(sg)]
            fit = W4.branch_photometry(sub, RV)
            f = sub.set_index("source_id").loc[fit.source_id]
            fit = fit.assign(
                av_err=np.nan_to_num(f[f"av_err_rv{RV:.1f}"].to_numpy(), nan=0.0),
                sig_g0=np.sqrt(np.nan_to_num(f.G_err.to_numpy(), nan=0.02) ** 2 + W4.MAG_FLOOR ** 2),
                sig_c0=np.sqrt(np.nan_to_num(f.BP_err.to_numpy(), nan=0.02) ** 2
                               + np.nan_to_num(f.RP_err.to_numpy(), nan=0.02) ** 2 + W4.MAG_FLOOR ** 2))
            for edge in (W4.UMS_FAINT_EDGE, -1.0):
                stars = fit[(fit.MG0 <= edge) & np.isfinite(fit.colour) & np.isfinite(fit.MG0)]
                for correlated in (False, True):
                    total = np.array([stars.P.to_numpy() @ loglike(stars, *p, edge, correlated, kg, kc)
                                      for p in parts])
                    s = I20.posterior_summary(ages, total)
                    rows.append({"family": fam, "subgroup": sg, "window_edge": edge,
                                 "av_error_model": "correlated" if correlated else "independent (WP4)",
                                 "n_stars": len(stars), "av_err_median": float(np.median(stars.av_err)),
                                 **{x: s[x] for x in ("map", "lo68", "hi68", "railed")}})
                    print(rows[-1], flush=True)
            ref = stored[stored.subgroup.eq(sg) & stored.family.eq(fam) & stored.R_V.eq(RV)
                         & stored.f_bin.eq(F_BIN) & stored.indicator.eq("ums") & stored.dmu.eq(0.0)]
            rows.append({"family": fam, "subgroup": sg, "window_edge": W4.UMS_FAINT_EDGE,
                         "av_error_model": "stored repair_v5 WP4 MAP", "map": float(ref.age_map.iloc[0])})
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue21_error_model_diagnostic.py",
        "status": "post-hoc diagnostic, not pre-registered, nothing adopted",
        "inputs": {p: w.sha256(w.ROOT / p) for p in INPUTS},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })


if __name__ == "__main__":
    main()
