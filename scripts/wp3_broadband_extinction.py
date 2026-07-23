#!/usr/bin/env python3
"""WP3 step: broadband multiband extinction A_V for non-anchor members.

For each member we fit observed G/BP/RP/J/H/Ks against reddened synthetic
photometry from the WP3 isochrone families at the group distance modulus:

    m_X = M_X0(template) + mu + k_X(R_V) * A_V

The template library is the union of PARSEC and MIST isochrone points, ages
1-10 Myr, solar Z (the same families WP4 will use -- consistency requirement).
For every template, A_V is the closed-form WLS optimum over the available bands;
the star's A_V estimate is the chi^2-posterior-weighted mean over templates, with
the weighted std as the statistical error. Runs independently for R_V=3.0/3.1/3.5.

The 2MASS J/H/Ks bands break the Teff/A_V degeneracy (extinction is small in the
IR), so stars WITHOUT a 2MASS match get a genuinely weaker constraint; their
errors are additionally inflated (plan WP3 step 2) rather than dropping them.

Outputs: data/processed/wp3_broadband_extinction.parquet
Provenance: provenance/wp3_broadband_extinction_execution.json
"""
import hashlib, json
import numpy as np, pandas as pd
from wp3_common import BANDS, DIST_MODULUS
from wp3_extinction_law import band_coefficients, R_V_BRANCHES

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

NO2MASS_ERR_INFLATE = 1.8   # extra multiplicative error inflation for optical-only fits
AV_MIN, AV_MAX = -1.0, 30.0  # allowed A_V search span (negatives kept as noise diagnostic)

def load_templates():
    cols = ["logAge", "Mini", "logTe", "family"] + [f"{b}0" for b in BANDS]
    par = pd.read_parquet("data/processed/wp3_isochrones_parsec.parquet")
    mist = pd.read_parquet("data/processed/wp3_isochrones_mist.parquet")
    par = par[[c for c in cols if c in par.columns]]
    mist = mist[[c for c in cols if c in mist.columns]]
    t = pd.concat([par, mist], ignore_index=True).dropna(subset=[f"{b}0" for b in BANDS])
    M0 = np.column_stack([t[f"{b}0"].values for b in BANDS])   # (T, B)
    return t.reset_index(drop=True), M0

def fit_star(y, w, M0, k):
    """y = m_X - mu (B,), w weights (B,), M0 (T,B), k (B,).
    Returns posterior-weighted A_V, its std, chi2_min, best template index."""
    wk = w * k                       # (B,)
    den = np.sum(wk * k)             # scalar
    if den <= 0:
        return np.nan, np.nan, np.nan, -1
    num = (y * wk)[None, :].sum(axis=1) - (M0 * wk[None, :]).sum(axis=1)   # (T,)
    av = num / den
    av = np.clip(av, AV_MIN, AV_MAX)
    resid = y[None, :] - M0 - av[:, None] * k[None, :]                     # (T,B)
    chi2 = (w[None, :] * resid ** 2).sum(axis=1)                           # (T,)
    chi2 -= np.nanmin(chi2)
    post = np.exp(-0.5 * chi2)
    post /= post.sum()
    av_mean = float(np.sum(post * av))
    av_std = float(np.sqrt(np.sum(post * (av - av_mean) ** 2)))
    best = int(np.argmin(chi2))
    return av_mean, av_std, float(np.min((y[None, :] - M0 - av[:, None] * k[None, :]) ** 2 @ w)), best

def main():
    ph = pd.read_parquet("data/processed/wp3_member_photometry.parquet")
    anc = pd.read_parquet("data/processed/wp3_anchor_extinction.parquet")
    anchor_ids = set(anc.source_id.dropna().astype("int64"))

    tmpl, M0 = load_templates()
    kX = {f"{rv:.1f}": np.array([band_coefficients(rv)[b] for b in BANDS]) for rv in R_V_BRANCHES}
    mu = DIST_MODULUS

    recs = []
    for _, r in ph.iterrows():
        sid = int(r.source_id)
        is_anchor = sid in anchor_ids
        mags = np.array([r.get(b, np.nan) for b in BANDS], float)
        errs = np.array([r.get(f"{b}_err", np.nan) for b in BANDS], float)
        avail = np.isfinite(mags) & np.isfinite(errs) & (errs > 0)
        has_2mass = bool(r.has_2mass)
        rec = {"source_id": sid, "is_anchor": is_anchor, "has_2mass": has_2mass,
               "n_bands": int(avail.sum()), "l_deg": r.l_deg, "b_deg": r.b_deg}
        errs_eff = errs.copy()
        if not has_2mass:
            errs_eff = errs_eff * NO2MASS_ERR_INFLATE
        y = mags - mu
        for rv, k in kX.items():
            if avail.sum() < 2:
                rec[f"av_rv{rv}"] = np.nan; rec[f"av_err_rv{rv}"] = np.nan
                continue
            w = np.where(avail, 1.0 / errs_eff ** 2, 0.0)
            av, avstd, chi2, best = fit_star(y, w, M0, k)
            # add a floor error and the 2MASS-missing inflation on the reported error
            av_err = np.sqrt(avstd ** 2 + 0.05 ** 2)
            if not has_2mass:
                av_err = np.sqrt(av_err ** 2 + 0.5 ** 2)   # explicit optical-only penalty
            rec[f"av_rv{rv}"] = av
            rec[f"av_err_rv{rv}"] = av_err
            if rv == "3.1":
                rec["best_template_family"] = tmpl.family.iloc[best] if best >= 0 else None
                rec["best_template_logTe"] = float(tmpl.logTe.iloc[best]) if (best >= 0 and "logTe" in tmpl) else np.nan
                rec["chi2_min"] = chi2
        recs.append(rec)

    out = pd.DataFrame(recs)
    non_anchor = out[~out.is_anchor]
    p = "data/processed/wp3_broadband_extinction.parquet"
    out.to_parquet(p, index=False)

    av = non_anchor["av_rv3.1"].dropna()
    prov = {
        "script": "scripts/wp3_broadband_extinction.py",
        "inputs": {
            "wp3_member_photometry": sha("data/processed/wp3_member_photometry.parquet"),
            "wp3_anchor_extinction": sha("data/processed/wp3_anchor_extinction.parquet"),
            "isochrones_parsec": sha("data/processed/wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha("data/processed/wp3_isochrones_mist.parquet"),
        },
        "model": "m_X = M_X0(template) + mu + k_X A_V; per-template WLS A_V; chi2 posterior over PARSEC+MIST templates",
        "distance_modulus_used": mu,
        "n_templates": int(len(tmpl)),
        "no2mass_error_inflation": NO2MASS_ERR_INFLATE,
        "n_members": int(len(out)),
        "n_non_anchor": int(len(non_anchor)),
        "n_non_anchor_no2mass": int((~non_anchor.has_2mass).sum()),
        "av_rv3.1_non_anchor_summary": {
            "median": float(av.median()), "mean": float(av.mean()),
            "p05": float(av.quantile(0.05)), "p95": float(av.quantile(0.95)),
            "min": float(av.min()), "max": float(av.max()),
            "n_negative": int((av < 0).sum()),
            "n_negative_beyond_2sigma": int((non_anchor["av_rv3.1"] + 2 * non_anchor["av_err_rv3.1"] < 0).sum()),
        },
        "output": {p: sha(p), "rows": int(len(out))},
    }
    json.dump(prov, open("provenance/wp3_broadband_extinction_execution.json", "w"), indent=2)
    s = prov["av_rv3.1_non_anchor_summary"]
    print(f"non-anchor members fit: {len(non_anchor)}  (no-2MASS {prov['n_non_anchor_no2mass']})")
    print(f"A_V(R_V=3.1): median={s['median']:.2f} mean={s['mean']:.2f} "
          f"[{s['p05']:.1f},{s['p95']:.1f}] range[{s['min']:.1f},{s['max']:.1f}]")
    print(f"negatives: {s['n_negative']}  beyond-2sigma-negative: {s['n_negative_beyond_2sigma']}")

if __name__ == "__main__":
    main()
