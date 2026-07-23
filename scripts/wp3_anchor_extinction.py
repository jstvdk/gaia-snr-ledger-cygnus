#!/usr/bin/env python3
"""WP3 step: extinction A_V for the spectroscopic anchors from intrinsic colours.

Rationale (plan WP3 + CUTS_AND_THRESHOLDS.md 5): for hot OB stars broadband
colours saturate and the Teff/extinction degeneracy is real, so A_V for anchors
is NOT derived by free-Teff broadband fitting. Instead we FIX the intrinsic SED
from the spectroscopic Teff (or spectral type) and solve only the reddening.

Model, per anchor, per R_V branch:  m_X = M_X0(Teff) + mu + k_X(R_V) * A_V
Two linear parameters (mu, A_V) fit by weighted least squares over available
bands. Key property: an incorrect luminosity class (dwarf vs supergiant) shifts
every band by the same amount and is fully absorbed by mu, so A_V is set purely
by the *colour* pattern, which is nearly logg-independent for hot stars. This is
exactly why the anchor treatment is robust where free broadband fitting is not.

Intrinsic absolute mags M_X0(Teff) come from the WP3 isochrone main-sequence
locus (PARSEC+MIST averaged; the family colour spread is tiny for hot stars and
is recorded). WR stars are flagged and excluded from the colour method.

Outputs: data/processed/wp3_anchor_extinction.parquet
Provenance: provenance/wp3_anchor_extinction_execution.json
"""
import hashlib, json, re
import numpy as np, pandas as pd
from wp3_common import BANDS, normalize_teff
from wp3_extinction_law import band_coefficients, R_V_BRANCHES

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

# --- spectral type -> Teff (K) for anchors lacking a tabulated Teff.
# O/B calibration bins (Pecaut & Mamajek 2013 scale, dwarf values; luminosity
# class is irrelevant here because it is absorbed by mu). Coarse but only used
# as a fallback for the handful of anchors without a spectroscopic Teff.
SPT_TEFF = {
    "O2": 54000, "O3": 44900, "O4": 42900, "O5": 40900, "O6": 38900,
    "O7": 36900, "O8": 34900, "O9": 32500, "B0": 31400, "B1": 26000,
    "B2": 20600, "B3": 17000, "B5": 15200, "B7": 13000, "B8": 12300, "B9": 10700,
}

def spt_to_teff(spt):
    if not isinstance(spt, str):
        return np.nan
    m = re.match(r"\s*([OB])\s*([0-9](?:\.[0-9])?)", spt)
    if not m:
        return np.nan
    letter, sub = m.group(1), float(m.group(2))
    key = f"{letter}{int(round(sub))}"
    if key in SPT_TEFF:
        return SPT_TEFF[key]
    # interpolate on the ordered grid
    order = list(SPT_TEFF)
    teffs = np.array([SPT_TEFF[k] for k in order])
    idx = np.array([("OB".index(k[0]) * 10 + float(k[1:])) for k in order])
    want = "OB".index(letter) * 10 + sub
    return float(np.interp(want, idx, teffs))

def build_intrinsic_locus():
    """M_X0 and colours as functions of logTe from the MS isochrone locus."""
    frames = []
    par = pd.read_parquet("data/processed/wp3_isochrones_parsec.parquet")
    frames.append(par[par.label.isin([0, 1])])
    mist = pd.read_parquet("data/processed/wp3_isochrones_mist.parquet")
    frames.append(mist[mist.phase <= 0])
    loci = {}
    for name, d in [("PARSEC", frames[0]), ("MIST", frames[1])]:
        d = d.dropna(subset=["logTe"])
        grid = np.linspace(3.9, 4.75, 60)  # ~7900 K to ~56000 K
        rel = {}
        for b in BANDS:
            col = f"{b}0"
            # median absolute mag per logTe bin
            vals = []
            for g in grid:
                s = d[np.abs(d.logTe - g) < 0.03]
                vals.append(np.median(s[col]) if len(s) >= 3 else np.nan)
            rel[b] = np.array(vals, float)
        loci[name] = (grid, rel)
    return loci

def M0_of_teff(loci, teff):
    """Family-averaged intrinsic absolute mags at Teff; returns dict band->M0
    and the family spread (max-min) as a systematic diagnostic."""
    logte = np.log10(teff)
    out, spread = {}, {}
    for b in BANDS:
        vals = []
        for name, (grid, rel) in loci.items():
            good = np.isfinite(rel[b])
            vals.append(np.interp(logte, grid[good], rel[b][good]))
        out[b] = float(np.mean(vals))
        spread[b] = float(np.max(vals) - np.min(vals))
    return out, spread

def fit_av(mags, errs, M0, kX):
    """WLS solve m_X = M_X0 + mu + k_X*A_V. Returns (mu, A_V, sigma_AV, nbands)."""
    bands = [b for b in BANDS if b in mags and np.isfinite(mags[b]) and np.isfinite(errs[b])]
    if len(bands) < 2:
        return np.nan, np.nan, np.nan, len(bands)
    y = np.array([mags[b] - M0[b] for b in bands])
    A = np.array([[1.0, kX[b]] for b in bands])
    w = np.array([1.0 / max(errs[b], 1e-3) ** 2 for b in bands])
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    beta = cov @ (A.T @ W @ y)
    mu, av = beta
    sig_av = float(np.sqrt(cov[1, 1]))
    # spread of k across used bands must be non-trivial for a real constraint
    if (max(kX[b] for b in bands) - min(kX[b] for b in bands)) < 0.1:
        sig_av = max(sig_av, 1.0)
    return float(mu), float(av), sig_av, len(bands)

def main():
    aa = pd.read_parquet("data/processed/wp2_anchor_assignments.parquet")
    aa["source_id"] = pd.to_numeric(aa["source_id"], errors="coerce").astype("Int64")
    aa = aa[aa.membership_probability > 0.5].copy()
    sa = pd.read_parquet("data/processed/wp1_spectroscopic_anchors.parquet")
    sa["source_id"] = pd.to_numeric(sa["source_id"], errors="coerce").astype("Int64")
    ph = pd.read_parquet("data/processed/wp3_member_photometry.parquet")

    sa_small = (sa.dropna(subset=["source_id"])
                  .sort_values("teff_K")
                  .drop_duplicates("source_id")
                  [["source_id", "teff_K", "spectral_type", "extinction_av_mag",
                    "teff_error_K"]])
    anc = aa.merge(sa_small, on="source_id", how="left", suffixes=("", "_sa"))
    anc = anc.merge(ph[["source_id", "G", "G_err", "BP", "BP_err", "RP", "RP_err",
                        "J", "J_err", "H", "H_err", "Ks", "Ks_err", "has_2mass",
                        "in_narrow_catalogue", "l_deg", "b_deg"]],
                    on="source_id", how="left")

    # Teff: normalize units, fall back to spectral type
    teff = normalize_teff(anc["teff_K"].values)
    spt = anc["spectral_type"].fillna(anc.get("spectral_type_sa"))
    for i in range(len(anc)):
        if not np.isfinite(teff[i]):
            teff[i] = spt_to_teff(spt.iloc[i])
    anc["teff_used"] = teff
    anc["is_wr"] = spt.fillna("").str.contains(r"W[NCR]", regex=True)

    loci = build_intrinsic_locus()
    kX_branches = {f"{rv:.1f}": band_coefficients(rv) for rv in R_V_BRANCHES}

    rows = []
    for _, r in anc.iterrows():
        base = {
            "source_id": r["source_id"], "object_name": r.get("object_name"),
            "spectral_type": spt.get(r.name) if hasattr(spt, "get") else r.get("spectral_type"),
            "teff_used_K": r["teff_used"], "membership_probability": r["membership_probability"],
            "l_deg": r.get("l_deg"), "b_deg": r.get("b_deg"),
            "literature_av": r.get("extinction_av_mag"), "is_wr": bool(r["is_wr"]),
            "has_gaia": bool(np.isfinite(r.get("G"))), "has_2mass": bool(r.get("has_2mass")),
        }
        mags = {b: r.get(b) for b in BANDS}
        errs = {b: r.get(f"{b}_err") for b in BANDS}
        n_ok = sum(1 for b in BANDS if np.isfinite(mags[b]))
        if r["is_wr"] or not np.isfinite(r["teff_used"]) or n_ok < 2:
            base["method"] = ("wr_flagged" if r["is_wr"]
                              else "no_teff" if not np.isfinite(r["teff_used"]) else "insufficient_phot")
            for rv in kX_branches:
                base[f"av_rv{rv}"] = r.get("extinction_av_mag")  # fall back to literature
                base[f"av_err_rv{rv}"] = 1.0 if np.isfinite(r.get("extinction_av_mag") or np.nan) else np.nan
            base["av_source"] = "literature_fallback"
            rows.append(base); continue

        M0, spread = M0_of_teff(loci, r["teff_used"])
        base["family_color_spread_G_Ks"] = spread["G"] + spread["Ks"]
        base["method"] = "intrinsic_color_spectroscopic"
        base["av_source"] = "intrinsic_color"
        for rv, kX in kX_branches.items():
            mu, av, sig, nb = fit_av(mags, errs, M0, kX)
            # add Teff systematic on A_V: propagate a Teff error into intrinsic colours
            base[f"av_rv{rv}"] = av
            base[f"av_err_rv{rv}"] = sig
            base["n_bands"] = nb
        rows.append(base)

    out = pd.DataFrame(rows)
    p = "data/processed/wp3_anchor_extinction.parquet"
    out.to_parquet(p, index=False)

    # validation vs literature (baseline branch)
    v = out[(out.av_source == "intrinsic_color") & out.literature_av.notna()]
    resid = (v["av_rv3.1"] - v["literature_av"]).values
    prov = {
        "script": "scripts/wp3_anchor_extinction.py",
        "inputs": {
            "wp2_anchor_assignments": sha("data/processed/wp2_anchor_assignments.parquet"),
            "wp1_spectroscopic_anchors": sha("data/processed/wp1_spectroscopic_anchors.parquet"),
            "wp3_member_photometry": sha("data/processed/wp3_member_photometry.parquet"),
            "isochrones_parsec": sha("data/processed/wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha("data/processed/wp3_isochrones_mist.parquet"),
        },
        "method": "fixed-Teff two-parameter (mu, A_V) WLS over available bands; R_V branches",
        "n_anchors_P_gt_0.5": int(len(out)),
        "n_intrinsic_color": int((out.av_source == "intrinsic_color").sum()),
        "n_wr_flagged": int(out.is_wr.sum()),
        "n_literature_fallback": int((out.av_source == "literature_fallback").sum()),
        "teff_from_spectral_type_count": int((~np.isfinite(anc["teff_K"].apply(
            lambda x: normalize_teff([x])[0])) & np.isfinite(anc["teff_used"])).sum()),
        "validation_vs_literature_av_rv3.1": {
            "n": int(len(v)),
            "median_residual": float(np.nanmedian(resid)),
            "mad": float(np.nanmedian(np.abs(resid - np.nanmedian(resid)))),
            "rms": float(np.sqrt(np.nanmean(resid ** 2))),
        },
        "av_rv3.1_summary": {
            "median": float(np.nanmedian(out["av_rv3.1"])),
            "p05": float(np.nanpercentile(out["av_rv3.1"].dropna(), 5)),
            "p95": float(np.nanpercentile(out["av_rv3.1"].dropna(), 95)),
            "n_negative": int((out["av_rv3.1"] < 0).sum()),
        },
        "output": {p: sha(p), "rows": int(len(out))},
    }
    json.dump(prov, open("provenance/wp3_anchor_extinction_execution.json", "w"), indent=2)
    print(f"anchors: {len(out)}  intrinsic_color={prov['n_intrinsic_color']}  "
          f"WR={prov['n_wr_flagged']}  lit_fallback={prov['n_literature_fallback']}")
    print(f"A_V(R_V=3.1) median={prov['av_rv3.1_summary']['median']:.2f}  "
          f"[{prov['av_rv3.1_summary']['p05']:.1f},{prov['av_rv3.1_summary']['p95']:.1f}]  "
          f"n_neg={prov['av_rv3.1_summary']['n_negative']}")
    print(f"vs literature A_V: n={prov['validation_vs_literature_av_rv3.1']['n']}  "
          f"median_resid={prov['validation_vs_literature_av_rv3.1']['median_residual']:+.2f}  "
          f"rms={prov['validation_vs_literature_av_rv3.1']['rms']:.2f}")

if __name__ == "__main__":
    main()
