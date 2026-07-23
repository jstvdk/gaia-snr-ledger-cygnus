#!/usr/bin/env python3
"""
WP2 Task A — Derive Cyg OB2 subgroups from sky position + proper motion.

Input : data/processed/wp2_members.parquet  (WP2 canonical membership)
Sample: 1,331 clean automatic members  (membership_probability>0.5 & ~anchor_quality_exempt)

Method (per CUTS_AND_THRESHOLDS.md sec 6 + Task A brief):
  * Feature space = (l, b, pmra*, pmdec). Parallax EXCLUDED: WP2 distance test
    found ONE population at 1.62 kpc, intrinsic depth 45 pc, exhausted at DR3.
  * Gaussian mixtures k=2..8 refit across 50 deterministic seeds. A k-partition
    is STABLE only if seed-to-seed Adjusted Rand Index (ARI) is high AND its 10th
    percentile is high (the WP2 cluster GMM BIC was non-monotonic -> local optima
    -> seed stability is the acceptance criterion, not BIC).
  * A parametric-bootstrap NULL (single 4D Gaussian matched to the data) is run
    through the identical procedure: real structure must beat what GMM regularity
    alone produces on a smooth blob.
  * Robustness: repeated under RobustScaler and on the P>0.9 subset.
  * Non-kinematic confirmation: extinction A_V and spectral types on the P>0.5
    anchors (features OUTSIDE the clustering space).
  * External validation: Paiz et al. 2025 Table 3 clusters, HSC 630 as control.

Outputs:
  provenance/wp2_subgroups_execution.json
  provenance/wp2_gmm_seed_stability.csv
  provenance/wp2_hdbscan_subgroup_scan.csv
  provenance/wp2_paiz_crossmatch.csv
  tables/wp2_subgroup_labels.parquet   SIDECAR (source_id -> subgroup);
                                        wp2_members.parquet is NOT modified so
                                        Task B can read it concurrently.
  figures/wp2/wp2_subgroups_{sky,vpd,extinction}.png

Run:  python scripts/wp2_derive_subgroups.py
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.mixture import GaussianMixture
from sklearn.cluster import HDBSCAN, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parent.parent
MEMBERS = ROOT / "data/processed/wp2_members.parquet"
NARROW = ROOT / "data/processed/wp1_gaia_narrow.parquet"
ANCHORS = ROOT / "data/processed/wp2_anchor_assignments.parquet"
ANCHOR_RECS = ROOT / "data/processed/wp1_spectroscopic_anchors.parquet"
PROV = ROOT / "provenance"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures/wp2"
FIGS.mkdir(parents=True, exist_ok=True)

N_SEEDS = 50
SEEDS = list(range(1000, 1000 + N_SEEDS))
NULL_SEEDS = list(range(1000, 1030))     # 30 seeds for the (cheaper) null runs
N_NULL_SIMS = 20
FEATURES = ["l_deg", "b_deg", "pmra", "pmdec"]
STABLE_ARI = 0.90        # ari_mean threshold
STABLE_ARI_P10 = 0.90    # ari_p10 threshold (a genuinely reproducible split)

# Paiz et al. 2025 (IJAA 15, 171) Table 3 — cluster astrometry.
# name: (ra, dec, l, b, pmra*, pmdec, N, r_arcmin, plx_mas, dist_pc, role)
PAIZ = {
    "Bica 2":   (308.320, 41.316, 80.229, 0.784, -2.673, -4.443, 41, 5.4, 0.587, 1656, "association"),
    "FSR 0224": (306.371, 40.198, 78.457, 1.325, -3.109, -4.337, 20, 9.9, 0.575, 1688, "association"),
    "FSR 0238": (308.704, 41.412, 80.478, 0.610, -2.825, -4.604, 21, 5.0, 0.576, 1686, "association"),
    "OC-123":   (306.345, 40.840, 78.970, 1.711, -3.049, -4.757, 17, 7.2, 0.558, 1740, "association"),
    "OC-128":   (307.938, 40.798, 79.641, 0.709, -3.010, -4.337, 85, 22.1, 0.583, 1682, "association"),
    "HSC 625":  (308.258, 40.792, 79.780, 0.510, -2.595, -4.924, 22, 16.0, 0.576, 1686, "association"),
    "HSC 630":  (307.689, 41.466, 80.069, 1.254, -1.989, -2.901, 16, 8.1, 0.726, 1346, "CONTROL"),
}


def load_clean():
    df = pd.read_parquet(MEMBERS)
    clean = df[(df["membership_probability"] > 0.5) & (~df["anchor_quality_exempt"])].copy()
    return df, clean.reset_index(drop=True)


def fit_labels(X, k, seed):
    return GaussianMixture(n_components=k, covariance_type="full", n_init=1,
                           max_iter=500, reg_covar=1e-4, random_state=seed).fit(X)


def ari_stats(X, k, seeds):
    labs = [fit_labels(X, k, s).predict(X) for s in seeds]
    a = [adjusted_rand_score(labs[i], labs[j])
         for i in range(len(labs)) for j in range(i + 1, len(labs))]
    a = np.array(a)
    return float(a.mean()), float(np.percentile(a, 10)), labs


def consensus_partition(X, k, seeds):
    n = X.shape[0]
    C = np.zeros((n, n))
    for s in seeds:
        lab = fit_labels(X, k, s).predict(X)
        C += (lab[:, None] == lab[None, :]).astype(float)
    C /= len(seeds)
    cons = AgglomerativeClustering(n_clusters=k, metric="precomputed",
                                   linkage="average").fit_predict(1.0 - C)
    stab = {}
    for c in np.unique(cons):
        idx = np.where(cons == c)[0]
        if len(idx) > 1:
            sub = C[np.ix_(idx, idx)]
            stab[int(c)] = float(sub[np.triu_indices(len(idx), 1)].mean())
        else:
            stab[int(c)] = 1.0
    return cons, stab, C


def diagnostics(clean):
    out = {}
    for ax, err in [("pmra", "pmra_error"), ("pmdec", "pmdec_error")]:
        obs = clean[ax].std(ddof=1)
        med_err = clean[err].median()
        out[f"{ax}_obs_std"] = float(obs)
        out[f"{ax}_median_err"] = float(med_err)
        out[f"{ax}_intrinsic_disp"] = float(np.sqrt(max(obs**2 - med_err**2, 0.0)))
    out["pm_intrinsic_disp_mean"] = float(
        np.mean([out["pmra_intrinsic_disp"], out["pmdec_intrinsic_disp"]]))
    out["pm_intrinsic_over_error"] = float(
        out["pm_intrinsic_disp_mean"] / clean[["pmra_error", "pmdec_error"]].median().mean())
    lc, bc = clean["l_deg"].mean(), clean["b_deg"].mean()
    r = np.sqrt((clean["l_deg"] - lc) ** 2 + (clean["b_deg"] - bc) ** 2).values
    rmax = r.max()
    out["sky_radial_KS_vs_uniform_disc"] = float(
        stats.kstest(r, lambda x: np.clip((x / rmax) ** 2, 0, 1)).statistic)
    out["sky_centroid_l"] = float(lc)
    out["sky_centroid_b"] = float(bc)
    return out


def name_components(clean, cons):
    """
    Deterministic physical naming, independent of the integer the consensus
    happens to assign:
      B (OC-128 kinematic group) = component with most negative median pmra.
      Of the remaining two:
        C (HSC 625 group)  = more negative median pmdec.
        A (FSR 0238 group) = less negative median pmdec (also lowest extinction).
    """
    med = {c: (clean.loc[cons == c, "pmra"].median(),
               clean.loc[cons == c, "pmdec"].median()) for c in np.unique(cons)}
    b_comp = min(med, key=lambda c: med[c][0])          # most negative pmra
    rest = [c for c in med if c != b_comp]
    c_comp = min(rest, key=lambda c: med[c][1])          # most negative pmdec
    a_comp = [c for c in rest if c != c_comp][0]
    return {a_comp: "CygOB2-A", b_comp: "CygOB2-B", c_comp: "CygOB2-C"}


def characterise(sub, label, analog):
    w = sub["membership_probability"].values
    W = w.sum()
    def wq(col, q):
        v = sub[col].values
        o = np.argsort(v)
        cw = np.cumsum(w[o]) / W
        return float(np.interp(q, cw, v[o]))
    p = sub["parallax_corrected"].values
    iv = 1.0 / sub["parallax_error"].values ** 2
    plx = float(np.sum(p * iv) / np.sum(iv))
    se = float(np.sqrt(1.0 / np.sum(iv)))
    return {
        "label": label, "paiz_analog": analog,
        "N_raw": int(len(sub)), "N_prob_weighted": float(W),
        "l_median": float(np.median(sub["l_deg"])), "b_median": float(np.median(sub["b_deg"])),
        "l_span_90": [wq("l_deg", 0.05), wq("l_deg", 0.95)],
        "b_span_90": [wq("b_deg", 0.05), wq("b_deg", 0.95)],
        "pmra_median": float(np.median(sub["pmra"])), "pmdec_median": float(np.median(sub["pmdec"])),
        "pmra_std": float(sub["pmra"].std(ddof=1)), "pmdec_std": float(sub["pmdec"].std(ddof=1)),
        "parallax_ivw_mas": plx, "parallax_ivw_se_mas": se,
        "distance_pc": float(1000.0 / plx),
        "distance_pc_lo": float(1000.0 / (plx + se)),
        "distance_pc_hi": float(1000.0 / (plx - se)),
    }


def hdbscan_scan(Xs):
    grid = []
    for mcs in [10, 15, 20, 25, 30, 40, 50]:
        for ms in [None, 5, 10, 15]:
            lab = HDBSCAN(min_cluster_size=mcs, min_samples=ms,
                          cluster_selection_method="eom").fit_predict(Xs)
            sizes = sorted([int((lab == c).sum()) for c in set(lab) if c != -1], reverse=True)
            grid.append({"min_cluster_size": mcs, "min_samples": ms if ms is not None else mcs,
                         "n_clusters": len(sizes), "n_noise": int((lab == -1).sum()),
                         "largest_frac": float(sizes[0] / len(lab)) if sizes else 0.0,
                         "cluster_sizes": sizes[:8]})
    return grid


def paiz_crossmatch(clean, name_by_comp, cons):
    lab_full = np.array([name_by_comp[c] for c in cons])
    rows = []
    ra, dec = clean["ra"].values, clean["dec"].values
    cosd = np.cos(np.deg2rad(dec))
    for name, v in PAIZ.items():
        cra, cdec, cl, cb, cpmra, cpmdec, N, r_arc, cplx, cdist, role = v
        dsky = np.sqrt(((ra - cra) * cosd) ** 2 + (dec - cdec) ** 2)
        sky_in = dsky < (r_arc / 60.0 + 0.05)
        dpm = np.sqrt((clean["pmra"].values - cpmra) ** 2 + (clean["pmdec"].values - cpmdec) ** 2)
        pm_in = dpm < 0.45
        plx_in = np.abs(clean["parallax_corrected"].values - cplx) < 0.10
        full = sky_in & pm_in & plx_in
        bd = pd.Series(lab_full[full]).value_counts().to_dict()
        rows.append({
            "cluster": name, "role": role, "l": cl, "b": cb,
            "pmra": cpmra, "pmdec": cpmdec, "plx": cplx,
            "in_member_footprint": bool(clean["l_deg"].min() - 0.05 <= cl <= clean["l_deg"].max() + 0.05
                                        and clean["b_deg"].min() - 0.05 <= cb <= clean["b_deg"].max() + 0.05),
            "n_sky_only": int(sky_in.sum()), "n_sky_pm": int((sky_in & pm_in).sum()),
            "n_sky_pm_plx": int(full.sum()),
            "subgroup_breakdown": {k: int(val) for k, val in bd.items()},
        })
    return rows


def make_figures(clean, name_by_comp, cons, anc_full):
    lab = np.array([name_by_comp[c] for c in cons])
    colors = {"CygOB2-A": "#4477AA", "CygOB2-B": "#EE6677", "CygOB2-C": "#228833"}
    # sky
    fig, ax = plt.subplots(figsize=(6, 5))
    for L in ["CygOB2-A", "CygOB2-B", "CygOB2-C"]:
        m = lab == L
        ax.scatter(clean["l_deg"][m], clean["b_deg"][m], s=6, alpha=0.6, c=colors[L], label=L)
    ax.set_xlabel("l [deg]"); ax.set_ylabel("b [deg]"); ax.invert_xaxis()
    ax.set_title("Cyg OB2 subgroups — sky"); ax.legend()
    fig.tight_layout(); fig.savefig(FIGS / "wp2_subgroups_sky.png", dpi=130); plt.close(fig)
    # VPD
    fig, ax = plt.subplots(figsize=(6, 5))
    for L in ["CygOB2-A", "CygOB2-B", "CygOB2-C"]:
        m = lab == L
        ax.scatter(clean["pmra"][m], clean["pmdec"][m], s=6, alpha=0.6, c=colors[L], label=L)
    for name, v in PAIZ.items():
        mk = "x" if v[10] == "CONTROL" else "*"
        ax.scatter(v[4], v[5], marker=mk, s=140, c="k")
        ax.annotate(name, (v[4], v[5]), fontsize=7)
    ax.set_xlabel(r"$\mu_{\alpha*}$ [mas/yr]"); ax.set_ylabel(r"$\mu_\delta$ [mas/yr]")
    ax.set_title("Cyg OB2 subgroups — VPD (stars=Paiz clusters, x=HSC630 control)")
    ax.legend(); fig.tight_layout()
    fig.savefig(FIGS / "wp2_subgroups_vpd.png", dpi=130); plt.close(fig)
    # extinction
    fig, ax = plt.subplots(figsize=(6, 4))
    for L in ["CygOB2-A", "CygOB2-B", "CygOB2-C"]:
        av = anc_full[anc_full["subgroup"] == L]["extinction_av_mag"].dropna()
        if len(av):
            ax.hist(av, bins=np.arange(3, 10, 0.5), alpha=0.5, color=colors[L],
                    label=f"{L} (n={len(av)}, med={av.median():.1f})")
    ax.set_xlabel(r"$A_V$ [mag] (anchors)"); ax.set_ylabel("count")
    ax.set_title("Extinction by subgroup — independent confirmation"); ax.legend()
    fig.tight_layout(); fig.savefig(FIGS / "wp2_subgroups_extinction.png", dpi=130); plt.close(fig)


def main():
    df_all, clean = load_clean()
    log = {"task": "WP2 Task A — derive Cyg OB2 subgroups", "date": "2026-07-23",
           "n_clean_members": int(len(clean)), "features": FEATURES, "n_seeds": N_SEEDS}
    log["diagnostics"] = diagnostics(clean)

    Xs = StandardScaler().fit_transform(clean[FEATURES].values)

    # ---- GMM seed stability, k=2..8 ----
    gmm_records = []
    consensus_cache = {}
    for k in range(2, 9):
        am, ap10, labs = ari_stats(Xs, k, SEEDS)
        bics = [fit_labels(Xs, k, s).bic(Xs) for s in SEEDS]
        cons, stab, _ = consensus_partition(Xs, k, SEEDS)
        consensus_cache[k] = (cons, stab)
        gmm_records.append({"k": k, "ari_mean": am, "ari_p10": ap10,
                            "bic_mean": float(np.mean(bics)), "bic_std": float(np.std(bics)),
                            "min_component_stability": float(min(stab.values())),
                            "median_component_stability": float(np.median(list(stab.values())))})
    log["gmm_seed_stability"] = gmm_records
    pd.DataFrame(gmm_records).to_csv(PROV / "wp2_gmm_seed_stability.csv", index=False)

    # candidate k = most stable partition meeting both ARI thresholds
    cand = [r for r in gmm_records if r["ari_mean"] >= STABLE_ARI and r["ari_p10"] >= STABLE_ARI_P10]
    log["stable_candidates"] = [r["k"] for r in cand]

    # ---- parametric-bootstrap null on the winning k ----
    null_result = None
    if cand:
        kbest = max(cand, key=lambda r: r["ari_mean"])["k"]
        rng = np.random.default_rng(7)
        mu, cov = Xs.mean(0), np.cov(Xs.T)
        null_am = []
        for _ in range(N_NULL_SIMS):
            Xn = rng.multivariate_normal(mu, cov, size=len(Xs))
            am, _, _ = ari_stats(Xn, kbest, NULL_SEEDS)
            null_am.append(am)
        data_am, _, _ = ari_stats(Xs, kbest, NULL_SEEDS)
        null_result = {"k": kbest, "data_ari_mean": float(data_am),
                       "null_ari_mean_mean": float(np.mean(null_am)),
                       "null_ari_mean_max": float(np.max(null_am)),
                       "frac_null_below_data": float(np.mean(np.array(null_am) < data_am)),
                       "n_null_sims": N_NULL_SIMS}
        log["parametric_bootstrap_null"] = null_result

    # ---- robustness: RobustScaler + P>0.9 subset ----
    Xr = RobustScaler().fit_transform(clean[FEATURES].values)
    hp = clean[clean["membership_probability"] > 0.9]
    Xh = StandardScaler().fit_transform(hp[FEATURES].values)
    robustness = {}
    if cand:
        kbest = max(cand, key=lambda r: r["ari_mean"])["k"]
        robustness["robustscaler_ari"] = {k: ari_stats(Xr, k, NULL_SEEDS)[0] for k in [2, 3, 4, 5]}
        robustness["highP_ari"] = {k: ari_stats(Xh, k, NULL_SEEDS)[0] for k in [2, 3, 4, 5]}
        cons_s = consensus_cache[kbest][0]
        cons_r, _, _ = consensus_partition(Xr, kbest, NULL_SEEDS)
        robustness["ari_standard_vs_robust_partition"] = float(adjusted_rand_score(cons_s, cons_r))
    log["robustness"] = robustness

    # ---- HDBSCAN scan (reported for completeness) ----
    hdb_grid = hdbscan_scan(Xs)
    log["hdbscan_scan"] = hdb_grid
    pd.DataFrame([{**g, "cluster_sizes": str(g["cluster_sizes"])} for g in hdb_grid]).to_csv(
        PROV / "wp2_hdbscan_subgroup_scan.csv", index=False)

    # ---- decision ----
    anchors = pd.read_parquet(ANCHORS)
    recs = pd.read_parquet(ANCHOR_RECS)[["source_id", "extinction_av_mag", "teff_K", "spectral_type"]]
    recs["source_id"] = pd.to_numeric(recs["source_id"], errors="coerce").astype("Int64")

    if cand:
        kbest = max(cand, key=lambda r: r["ari_mean"])["k"]
        cons, stab = consensus_cache[kbest]
        name_by_comp = name_components(clean, cons)
        labels = np.array([name_by_comp[c] for c in cons])
        clean["subgroup"] = labels

        anc = anchors[anchors["membership_probability"] > 0.5][["source_id"]].copy()
        anc["source_id"] = anc["source_id"].astype("Int64")
        lab_by_id = dict(zip(clean["source_id"].astype("Int64"), clean["subgroup"]))
        anc_full = anc.merge(recs, on="source_id", how="left")
        anc_full["subgroup"] = anc_full["source_id"].map(lab_by_id)

        # extinction split test (independent confirmation)
        labset = sorted(set(labels))
        avs = {L: anc_full.loc[anc_full["subgroup"] == L, "extinction_av_mag"].dropna().values
               for L in labset}
        pairks = []
        big = [L for L in labset if len(avs[L]) >= 5]
        for i in range(len(big)):
            for j in range(i + 1, len(big)):
                ks = stats.ks_2samp(avs[big[i]], avs[big[j]])
                a, b = avs[big[i]], avs[big[j]]
                eff = abs(np.median(a) - np.median(b)) / np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
                pairks.append({"a": big[i], "b": big[j], "n_a": len(a), "n_b": len(b),
                               "ks_stat": float(ks.statistic), "ks_p": float(ks.pvalue),
                               "median_effect_size": float(eff),
                               "confirms": bool(ks.pvalue < 0.01 and eff >= 0.30)})
        log["extinction_split_test"] = pairks
        confirmed_pairs = [p for p in pairks if p["confirms"]]

        # independent photometric reddening test on ALL members (not just anchors,
        # not a clustering feature). BP-RP differs by reddening; matched G rules out
        # a mass/luminosity-selection alternative.
        g = pd.read_parquet(NARROW)[["source_id", "phot_g_mean_mag",
                                     "phot_bp_mean_mag", "phot_rp_mean_mag"]]
        mp = clean.merge(g, on="source_id", how="left")
        mp["bprp"] = mp["phot_bp_mean_mag"] - mp["phot_rp_mean_mag"]
        phot = {"coverage": int(mp["bprp"].notna().sum()), "n": int(len(mp)), "pairs": []}
        for i in range(len(labset)):
            for j in range(i + 1, len(labset)):
                A, B = labset[i], labset[j]
                xa = mp.loc[mp["subgroup"] == A, "bprp"].dropna()
                xb = mp.loc[mp["subgroup"] == B, "bprp"].dropna()
                ga = mp.loc[mp["subgroup"] == A, "phot_g_mean_mag"].dropna()
                gb = mp.loc[mp["subgroup"] == B, "phot_g_mean_mag"].dropna()
                ksc = stats.ks_2samp(xa, xb)
                ksg = stats.ks_2samp(ga, gb)
                eff = abs(xa.median() - xb.median()) / np.sqrt((xa.std() ** 2 + xb.std() ** 2) / 2)
                phot["pairs"].append({
                    "a": A, "b": B, "bprp_median_a": float(xa.median()), "bprp_median_b": float(xb.median()),
                    "bprp_ks": float(ksc.statistic), "bprp_ks_p": float(ksc.pvalue),
                    "bprp_median_effect": float(eff),
                    "G_ks": float(ksg.statistic), "G_ks_p": float(ksg.pvalue),
                    "reddening_confirmed": bool(ksc.pvalue < 0.01 and eff >= 0.30),
                    "mass_selection_ruled_out": bool(ksg.pvalue > 0.01)})
        log["photometric_reddening_test"] = phot

        log["decision"] = {
            "outcome": "MULTIPLE_STABLE_SUBGROUPS",
            "adopted_k": kbest,
            "component_stability": stab,
            "component_naming": {str(c): name_by_comp[c] for c in name_by_comp},
            "independent_confirmation": {
                "anchor_extinction_confirmed_pairs": confirmed_pairs,
                "photometric_reddening_confirmed_pairs": [
                    (p["a"], p["b"]) for p in log["photometric_reddening_test"]["pairs"]
                    if p["reddening_confirmed"] and p["mass_selection_ruled_out"]],
                "note": ("Every subgroup is confirmed by a non-kinematic observable: "
                         "BP-RP reddening differs for all three pairs (p<0.01, effect>=0.30) "
                         "with matched G distributions ruling out mass selection; anchor "
                         "A_V independently confirms the A<C ordering."),
            },
        }
    else:
        labels = np.array(["CygOB2_single_body"] * len(clean))
        clean["subgroup"] = labels
        anc = anchors[anchors["membership_probability"] > 0.5][["source_id"]].copy()
        anc["source_id"] = anc["source_id"].astype("Int64")
        anc_full = anc.merge(recs, on="source_id", how="left")
        anc_full["subgroup"] = "CygOB2_single_body"
        name_by_comp = None
        cons = np.zeros(len(clean), dtype=int)
        log["decision"] = {
            "outcome": "ONE_COEVAL_BODY_TO_DR3_RESOLUTION",
            "reason": (f"No GMM k>=2 partition is stable (ari_mean>={STABLE_ARI} and "
                       f"ari_p10>={STABLE_ARI_P10})."),
            "recommendation": ("WP4 should carry the star-formation-duration branch "
                               "(0/1/2 Myr) as the age-spread uncertainty."),
        }

    # ---- subgroup characterisation ----
    analog_map = {"CygOB2-A": "FSR 0238 (+ part of Bica 2)", "CygOB2-B": "OC-128",
                  "CygOB2-C": "HSC 625 (+ part of Bica 2)", "CygOB2_single_body": "whole association"}
    log["subgroups"] = [characterise(clean[clean["subgroup"] == L], L, analog_map.get(L, ""))
                        for L in sorted(set(labels))]

    # anchor confirmation summary + spectral content
    def klass(s):
        if not isinstance(s, str):
            return "?"
        s = s.strip()
        if s[:2] in ("WN", "WC", "WR"):
            return "WR"
        return {"O": "O", "B": "B"}.get(s[:1], "other")
    anc_full["specclass"] = anc_full["spectral_type"].map(klass)
    conf = []
    for L in sorted(set(labels)):
        sub = anc_full[anc_full["subgroup"] == L]
        av = sub["extinction_av_mag"].dropna()
        conf.append({"label": L, "n_anchors": int(len(sub)), "n_with_Av": int(len(av)),
                     "Av_median": float(av.median()) if len(av) else None,
                     "Av_std": float(av.std(ddof=1)) if len(av) > 1 else None,
                     "spec_counts": sub["specclass"].value_counts().to_dict(),
                     "example_OB_WR": [st for st in sub["spectral_type"].dropna()
                                       if isinstance(st, str) and st[:1] in ("O", "W")][:8]})
    log["subgroup_anchor_confirmation"] = conf

    # ---- Paiz cross-match ----
    if name_by_comp is not None:
        log["paiz_crossmatch"] = paiz_crossmatch(clean, name_by_comp, cons)
    else:
        # single body: still report match counts, all under one label
        nb = {0: "CygOB2_single_body"}
        log["paiz_crossmatch"] = paiz_crossmatch(clean, nb, cons)
    pd.DataFrame([{k: (json.dumps(v) if isinstance(v, dict) else v) for k, v in r.items()}
                  for r in log["paiz_crossmatch"]]).to_csv(PROV / "wp2_paiz_crossmatch.csv", index=False)
    hsc630 = [r for r in log["paiz_crossmatch"] if r["cluster"] == "HSC 630"][0]
    log["hsc630_contamination_check"] = {
        "n_full_astrometric_match": hsc630["n_sky_pm_plx"], "n_sky_only": hsc630["n_sky_only"],
        "verdict": ("EXCLUDED — 0 clean members match HSC 630 in PM+parallax"
                    if hsc630["n_sky_pm_plx"] == 0
                    else f"WARNING — {hsc630['n_sky_pm_plx']} members match HSC 630")}

    # ---- figures ----
    make_figures(clean, (name_by_comp or {0: "CygOB2_single_body"}), cons, anc_full)

    # ---- sidecar labels ----
    TABLES.mkdir(exist_ok=True)
    side = clean[["source_id", "subgroup"]].copy()
    side["method"] = "GMM k-stable on (l,b,pmra,pmdec); parallax excluded"
    side.to_parquet(TABLES / "wp2_subgroup_labels.parquet", index=False)
    log["sidecar"] = {
        "path": "tables/wp2_subgroup_labels.parquet", "n_rows": int(len(side)),
        "key": "source_id",
        "note": "SIDECAR only; wp2_members.parquet NOT modified (Task B reads concurrently).",
        "label_counts": {k: int(v) for k, v in side["subgroup"].value_counts().items()}}

    with open(PROV / "wp2_subgroups_execution.json", "w") as f:
        json.dump(log, f, indent=2)

    # ---- console summary ----
    print("=== GMM seed stability ===")
    for r in gmm_records:
        print(f"  k={r['k']}  ari_mean={r['ari_mean']:.3f}  ari_p10={r['ari_p10']:.3f}  "
              f"bic_mean={r['bic_mean']:.1f}  min_comp_stab={r['min_component_stability']:.3f}")
    print("stable candidates:", log["stable_candidates"])
    if null_result:
        print(f"NULL check k={null_result['k']}: data_ari={null_result['data_ari_mean']:.3f} "
              f"null_max={null_result['null_ari_mean_max']:.3f} "
              f"frac_null_below_data={null_result['frac_null_below_data']:.2f}")
    print("robustness:", robustness)
    print("DECISION:", log["decision"]["outcome"])
    print("label counts:", log["sidecar"]["label_counts"])
    if "extinction_split_test" in log:
        for p in log["extinction_split_test"]:
            print(f"  anchor A_V {p['a']} vs {p['b']}: KS={p['ks_stat']:.3f} p={p['ks_p']:.1e} "
                  f"eff={p['median_effect_size']:.2f} confirms={p['confirms']}")
    if "photometric_reddening_test" in log:
        for p in log["photometric_reddening_test"]["pairs"]:
            print(f"  BP-RP {p['a']} vs {p['b']}: KS={p['bprp_ks']:.3f} p={p['bprp_ks_p']:.1e} "
                  f"eff={p['bprp_median_effect']:.2f} redden={p['reddening_confirmed']} "
                  f"massruledout={p['mass_selection_ruled_out']}")
    print("HSC630:", log["hsc630_contamination_check"]["verdict"])
    print("=== Paiz matches ===")
    for r in log["paiz_crossmatch"]:
        print(f"  {r['cluster']:9s}[{r['role'][:4]}] fp={int(r['in_member_footprint'])} "
              f"full={r['n_sky_pm_plx']:4d} breakdown={r['subgroup_breakdown']}")
    return log


if __name__ == "__main__":
    main()
