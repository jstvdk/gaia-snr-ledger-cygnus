#!/usr/bin/env python3
"""WP3 step: assemble the per-star extinction catalogue and de-reddened photometry.

Combines the two extinction estimators:
  * spectroscopic anchors  -> intrinsic-colour A_V (wp3_anchor_extinction)
  * all other members      -> broadband multiband A_V (wp3_broadband_extinction)
for each R_V branch (3.0 / 3.1 / 3.5). Produces de-reddened absolute magnitudes
and colours using the single-population distance posterior (1.6245 kpc); the
brightest, well-measured members use their individual parallax instead.

Outputs:
  data/processed/wp3_extinction.parquet   (all branches, one row per member)
  tables/wp3_extinction_rv{30,31,35}.cat  (per-branch ECSV catalogues)
Provenance: provenance/wp3_build_catalogue_execution.json
"""
import hashlib, json
import numpy as np, pandas as pd
from astropy.table import Table
from wp3_common import (BANDS, DIST_MODULUS, DIST_MODULUS_ERR, D_KPC, D_KPC_ERR,
                        BRIGHT_G_MAX, PARALLAX_SNR_MIN)
from wp3_extinction_law import band_coefficients, R_V_BRANCHES

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

def main():
    ph = pd.read_parquet("data/processed/wp3_member_photometry.parquet")
    anc = pd.read_parquet("data/processed/wp3_anchor_extinction.parquet")
    bb = pd.read_parquet("data/processed/wp3_broadband_extinction.parquet")
    anc = anc.dropna(subset=["source_id"]).copy()
    anc["source_id"] = anc["source_id"].astype("int64")
    anchor_ids = set(anc.source_id)

    df = ph.copy()
    # per-member A_V per branch and method
    branches = [f"{rv:.1f}" for rv in R_V_BRANCHES]
    anc_idx = anc.set_index("source_id")
    bb_idx = bb.set_index("source_id")

    method, avsrc = [], []
    av = {b: [] for b in branches}
    averr = {b: [] for b in branches}
    for sid in df.source_id.astype("int64"):
        if sid in anchor_ids and anc_idx.loc[sid, "av_source"] == "intrinsic_color":
            src = "anchor_intrinsic_color"; row = anc_idx.loc[sid]
        elif sid in anchor_ids and pd.notna(anc_idx.loc[sid].get("av_rv3.1")):
            src = "anchor_literature"; row = anc_idx.loc[sid]
        else:
            src = "broadband"; row = bb_idx.loc[sid]
        method.append(row.get("method", "broadband") if src != "broadband" else "broadband_multiband")
        avsrc.append(src)
        for b in branches:
            av[b].append(float(row.get(f"av_rv{b}", np.nan)))
            averr[b].append(float(row.get(f"av_err_rv{b}", np.nan)))
    df["av_method"] = method
    df["av_source"] = avsrc
    for b in branches:
        df[f"av_rv{b}"] = av[b]
        df[f"av_err_rv{b}"] = averr[b]

    # distance choice: group posterior, override with individual parallax for
    # the brightest, well-measured members.
    plx = df["parallax_corrected"].values
    plx_err = df["parallax_error"].values
    snr = np.where(plx_err > 0, plx / plx_err, 0.0)
    use_indiv = (df["G"].values < BRIGHT_G_MAX) & (snr > PARALLAX_SNR_MIN) & (plx > 0)
    mu_ind = np.where(plx > 0, 5.0 * np.log10(1000.0 / np.clip(plx, 1e-6, None)) - 5.0, np.nan)
    mu = np.where(use_indiv, mu_ind, DIST_MODULUS)
    df["distance_basis"] = np.where(use_indiv, "individual_parallax", "group_posterior")
    df["dist_modulus_used"] = mu

    # de-reddened absolute magnitudes + colours per branch
    kcoef = {b: band_coefficients(float(b)) for b in branches}
    for b in branches:
        k = kcoef[b]
        AV = df[f"av_rv{b}"].values
        for band in BANDS:
            A_band = k[band] * AV
            df[f"{band}0_abs_rv{b}"] = df[band].values - mu - A_band
            if band == "G":
                df[f"A_G_rv{b}"] = A_band
        df[f"BPRP0_rv{b}"] = (df["BP"].values - k["BP"] * AV) - (df["RP"].values - k["RP"] * AV)
        df[f"GKs0_rv{b}"] = (df["G"].values - k["G"] * AV) - (df["Ks"].values - k["Ks"] * AV)

    out_p = "data/processed/wp3_extinction.parquet"
    df.to_parquet(out_p, index=False)

    # per-branch ECSV catalogues (wp3_extinction.cat, one per R_V branch)
    cat_paths = {}
    for b in branches:
        cols = ["source_id", "ra", "dec", "l_deg", "b_deg", "membership_probability",
                "av_method", "av_source", "has_2mass", "n_bands", "distance_basis",
                "dist_modulus_used", f"av_rv{b}", f"av_err_rv{b}", f"A_G_rv{b}",
                "G", "BP", "RP", "J", "H", "Ks",
                f"G0_abs_rv{b}", f"BPRP0_rv{b}", f"GKs0_rv{b}"]
        cols = [c for c in cols if c in df.columns]
        sub = df[cols].rename(columns={f"av_rv{b}": "A_V", f"av_err_rv{b}": "A_V_err",
                                       f"A_G_rv{b}": "A_G", f"G0_abs_rv{b}": "M_G0",
                                       f"BPRP0_rv{b}": "BP_RP_0", f"GKs0_rv{b}": "G_Ks_0"})
        t = Table.from_pandas(sub)
        t.meta["R_V"] = float(b)
        t.meta["distance_kpc"] = D_KPC
        t.meta["distance_kpc_err"] = D_KPC_ERR
        t.meta["extinction_law"] = "CCM89 + O'Donnell 1994, R_V-dependent"
        t.meta["WP"] = "WP3"
        rvtag = b.replace(".", "")
        cp = f"tables/wp3_extinction_rv{rvtag}.cat"
        t.write(cp, format="ascii.ecsv", overwrite=True)
        cat_paths[b] = cp

    # summary
    summ = {}
    for b in branches:
        a = df[f"av_rv{b}"].dropna()
        summ[b] = {"median": float(a.median()), "p05": float(a.quantile(0.05)),
                   "p95": float(a.quantile(0.95)), "min": float(a.min()), "max": float(a.max()),
                   "n_negative": int((a < 0).sum()),
                   "n_negative_beyond_2sigma": int((df[f"av_rv{b}"] + 2 * df[f"av_err_rv{b}"] < 0).sum())}
    prov = {
        "script": "scripts/wp3_build_catalogue.py",
        "inputs": {k: sha(f"data/processed/wp3_{k}.parquet") for k in
                   ["member_photometry", "anchor_extinction", "broadband_extinction"]},
        "distance": {"kpc": D_KPC, "kpc_err": D_KPC_ERR, "modulus": DIST_MODULUS,
                     "modulus_err": DIST_MODULUS_ERR,
                     "bright_override": f"G<{BRIGHT_G_MAX} and parallax SNR>{PARALLAX_SNR_MIN}"},
        "n_members": int(len(df)),
        "n_distance_individual": int(use_indiv.sum()),
        "av_source_counts": df.av_source.value_counts().to_dict(),
        "av_summary_per_branch": summ,
        "outputs": {out_p: sha(out_p), **{cat_paths[b]: sha(cat_paths[b]) for b in branches}},
    }
    json.dump(prov, open("provenance/wp3_build_catalogue_execution.json", "w"), indent=2)
    print(f"members: {len(df)}  distance_individual: {int(use_indiv.sum())}")
    print("A_V source counts:", prov["av_source_counts"])
    for b in branches:
        s = summ[b]
        print(f"  R_V={b}: A_V median={s['median']:.2f} [{s['p05']:.1f},{s['p95']:.1f}] "
              f"max={s['max']:.1f} n_neg={s['n_negative']} (beyond2sig={s['n_negative_beyond_2sigma']})")
    print("wrote:", list(cat_paths.values()))

if __name__ == "__main__":
    main()
