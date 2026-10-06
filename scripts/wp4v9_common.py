#!/usr/bin/env python3
"""repair_v9 WP4 age redesign (issue #21) -- shared method.

Hashed into provenance/wp4v9_age_prereg.json before any real-data posterior or
injection is computed; the drivers refuse to run against a changed copy unless
provenance/wp4v9_deviations.json records the change.

M1  each star's likelihood is integrated over its own WP3 A_V posterior
    (wp3_extinction_posterior_repair_v5.npz), represented by K = 32 equal-
    probability quantile nodes:
        L_i(t) = (1/K) sum_k L_phot( x_obs - A_k r | t ),  r = (k_BP - k_RP, k_G)
    with L_phot WP4's particle mixture (wp4_common.build_model_particles,
    unchanged) under photometric errors + MAG_FLOOR + SIGMA_INT only.
M2  window membership on the star's posterior-MEDIAN M_G0, modelled as a
    selection function inside the likelihood (the cut is a cut on observed g
    at fixed edges + k_G med(A_V)); the population density is NOT truncated.
    Primary edge M_G0 <= -1.0.
Old the WP4 model: av_err added in quadrature to G and colour INDEPENDENTLY,
    A_V fixed at the WP3 value (exactly wp4_common.star_loglike).
Corr the Gaussian correlated-error model of issue21_error_model_diagnostic.py
    (secondary): A_V at the posterior median, sigma_A the posterior 68 %
    half-width, propagated along r with full covariance.

Nothing here writes a file.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

import wp4_common as W4
import wp5_common as w
from wp3_extinction_law import band_coefficients
from wp4_common import DIST_MODULUS

PREREG_PATH = w.PROVENANCE / "wp4v9_age_prereg.json"
DEVIATIONS_PATH = w.PROVENANCE / "wp4v9_deviations.json"
THIS_FILE = Path(__file__).resolve()

SUBGROUPS = ("CygOB2-A", "CygOB2-B", "CygOB2-C")
FAMILIES = ("PARSEC", "MIST")
R_V_BRANCHES = (3.0, 3.1, 3.5)
F_BINS = (0.3, 0.4, 0.5)
DECISION = {"R_V": 3.1, "f_bin": 0.4}

K_NODES = 32
WINDOWS = {                          # name -> (lower edge exclusive, upper edge inclusive)
    "bright": (-np.inf, -1.0),       # PRIMARY
    "bright_m0.5": (-np.inf, -0.5),
    "bright_m1.5": (-np.inf, -1.5),
    "faint": (-1.0, W4.UMS_FAINT_EDGE),
    "full": (-np.inf, W4.UMS_FAINT_EDGE),
}
PRIMARY_WINDOW = "bright"

INPUTS = {
    "extinction_v5": "data/processed/wp3_extinction_repair_v5.parquet",
    "av_posterior_v5": "data/processed/wp3_extinction_posterior_repair_v5.npz",
    "subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    "isochrones_parsec": "data/processed/wp3_isochrones_parsec.parquet",
    "isochrones_mist": "data/processed/wp3_isochrones_mist.parquet",
    "age_posteriors_v5": "data/processed/wp4_age_posteriors_repair_v5.parquet",
    "anchor_hrd_v8": "data/processed/wp4_anchor_hrd_repair_v8.parquet",
    "spectroscopic_tests": "tables/issue20b_age_tests.csv",
    "spectroscopic_curves": "provenance/issue20b_run_execution.json",
    "issue21_diagnostic": "tables/issue21_error_model_diagnostic.csv",
}


class FrozenInputMoved(RuntimeError):
    pass


def prereg() -> dict:
    if not PREREG_PATH.exists():
        raise FileNotFoundError("run scripts/wp4v9_prereg.py first")
    return json.loads(PREREG_PATH.read_text())


def check_method_unchanged() -> dict:
    out = {}
    for rel, frozen_hash in prereg()["method_code"].items():
        now = w.sha256(w.ROOT / rel)
        if now == frozen_hash:
            out[rel] = {"sha256": now, "deviation": None}
            continue
        devs = json.loads(DEVIATIONS_PATH.read_text()) if DEVIATIONS_PATH.exists() else {}
        hit = [d for d in devs.get("method_code_changes", [])
               if d.get("file") == rel and d.get("new_sha256") == now]
        if not hit:
            raise FrozenInputMoved(f"{rel} changed since the wp4v9 preregistration; record a deviation")
        out[rel] = {"sha256": now, "deviation": hit[0]}
    return out


def frozen(name: str) -> Path:
    entry = prereg()["consumed_inputs"][name]
    path = w.ROOT / entry["path"]
    if w.sha256(path) != entry["sha256"]:
        raise FrozenInputMoved(f"{entry['path']} changed since the wp4v9 preregistration")
    return path


# ------------------------------------------------------------- star inputs
def load_stars(ext: pd.DataFrame, labels: pd.DataFrame, npz, rv: float,
               k: int = K_NODES) -> pd.DataFrame:
    """Labelled members with a finite A_V posterior at this R_V and finite
    G, BP, RP.  Stars whose posterior is missing (9 of 1,392: insufficient
    photometry, one WR, four broadband stars) are excluded from every model."""
    if not np.array_equal(npz["source_id"], ext["source_id"].to_numpy()):
        raise RuntimeError("A_V posterior row order differs from the extinction table")
    r = int(np.argmin(np.abs(npz["rv"] - rv)))
    prob = npz["probability"][:, r, :]
    grid = npz["av_grid"]
    kb = band_coefficients(rv)
    frame = ext.drop(columns=["subgroup"], errors="ignore").reset_index(drop=True)
    frame["_row"] = np.arange(len(frame))
    frame = frame.merge(labels, on="source_id", how="left")
    total = prob.sum(1)
    frame = frame[frame.subgroup.isin(SUBGROUPS)
                  & np.isfinite(total[frame._row]) & (total[frame._row] > 0)
                  & frame.G.notna() & frame.BP.notna() & frame.RP.notna()].reset_index(drop=True)
    pr = prob[frame._row.to_numpy()]
    cdf = np.cumsum(pr, 1)
    cdf /= cdf[:, -1:]
    frame["_nodes"] = [np.interp((np.arange(k) + 0.5) / k, c, grid) for c in cdf]
    med = np.array([np.interp(0.5, c, grid) for c in cdf])
    lo = np.array([np.interp(0.16, c, grid) for c in cdf])
    hi = np.array([np.interp(0.84, c, grid) for c in cdf])
    frame["g_obs"] = frame["G"] - DIST_MODULUS
    frame["c_obs"] = frame["BP"] - frame["RP"]
    frame["sig_g"] = np.sqrt(np.nan_to_num(frame["G_err"], nan=0.02) ** 2
                             + W4.MAG_FLOOR ** 2 + W4.SIGMA_INT ** 2)
    frame["sig_c"] = np.sqrt(np.nan_to_num(frame["BP_err"], nan=0.02) ** 2
                             + np.nan_to_num(frame["RP_err"], nan=0.02) ** 2
                             + W4.MAG_FLOOR ** 2 + W4.SIGMA_INT ** 2)
    frame["av_median"] = med
    frame["av_halfwidth68"] = 0.5 * (hi - lo)
    frame["MG0_median"] = frame["g_obs"] - kb["G"] * med
    frame["colour_median"] = frame["c_obs"] - (kb["BP"] - kb["RP"]) * med
    frame["P"] = frame["membership_probability"]
    frame["is_anchor"] = frame["av_method"].eq("intrinsic_color_spectroscopic")
    frame.attrs.update(kG=kb["G"], kc=kb["BP"] - kb["RP"], R_V=rv)
    return frame


WIDEN_FACTOR = 2.0          # double-counting sensitivity (owner request, 2026-10-06)


def widened(stars: pd.DataFrame, factor: float = WIDEN_FACTOR) -> pd.DataFrame:
    """Copy of the star table with every A_V posterior widened by `factor`
    about its median (nodes -> med + factor (node - med), clipped at 0).  Window
    membership is unchanged (it uses the median).  Sensitivity for the double
    use of the photometry: WP3's A_V posterior was fitted to the same G, BP,
    RP, so marginalising over it counts the optical colours twice."""
    out = stars.copy()
    out["_nodes"] = [np.clip(m + factor * (n - m), 0.0, None)
                     for n, m in zip(out["_nodes"], out["av_median"])]
    out.attrs.update(stars.attrs)
    return out


def in_window(stars: pd.DataFrame, window: str) -> pd.DataFrame:
    lo, hi = WINDOWS[window]
    m = stars["MG0_median"].to_numpy()
    return stars[(m > lo) & (m <= hi)]


# ------------------------------------------------------------- likelihoods
def _window_particles(particles, window, dmu=0.0):
    colour, mg, weight = particles
    mg = mg + dmu
    lo, hi = WINDOWS[window]
    sel = (mg > lo) & (mg <= hi)
    if sel.sum() < 3:
        return None
    wt = weight[sel] / weight[sel].sum()
    return colour[sel], mg[sel], np.log(wt)


SELECTION_MARGIN = 3.0     # particles up to edge + 3 mag enter the selection integral


def m1_loglike(stars: pd.DataFrame, particles, window: str, dmu: float = 0.0,
               outlier_eps: float = 0.0, outlier_density: float = 0.0) -> np.ndarray:
    """Per-star ln L under M1 with the window as a selection function.

    A star is in the window when g_obs - k_G * med(A_V) lies inside the window
    edges; for that star this is a cut on its OBSERVED g at fixed edges
    (g_lo, g_hi] = edges + k_G med(A_V).  The likelihood is the untruncated
    population density marginalised over the star's A_V nodes, divided by the
    probability that a model star carrying the same A_V distribution and
    photometric error passes the same cut:

      L_i = mean_k sum_j w_j N(g_i - k_G A_k - G_j) N(c_i - k_c A_k - C_j)
            / mean_k sum_j w_j [Phi((g_hi - k_G A_k - G_j)/s) - Phi((g_lo - ...)/s)]
    """
    from scipy.special import log_ndtr
    colour, mg, weight = particles
    mg = mg + dmu
    lo, hi = WINDOWS[window]
    sel = mg <= hi + SELECTION_MARGIN
    cm, gm = colour[sel], mg[sel]
    lw = np.log(weight[sel] / weight.sum())
    kG, kc = stars.attrs["kG"], stars.attrs["kc"]
    nodes = np.stack(stars["_nodes"].to_numpy())
    g, c = stars["g_obs"].to_numpy(), stars["c_obs"].to_numpy()
    sg, sc = stars["sig_g"].to_numpy(), stars["sig_c"].to_numpy()
    med = stars["av_median"].to_numpy()
    out = np.empty(len(stars))
    for i in range(len(stars)):
        gd = g[i] - kG * nodes[i]
        cd = c[i] - kc * nodes[i]
        q = ((gd[:, None] - gm[None, :]) / sg[i]) ** 2 + ((cd[:, None] - cm[None, :]) / sc[i]) ** 2
        num = logsumexp(lw[None, :] - 0.5 * q) - np.log(2 * np.pi * sg[i] * sc[i]) - np.log(len(gd))
        # selection probability: observed g within (g_lo, g_hi]
        g_hi = hi + kG * med[i]
        mean_g = gm[None, :] + kG * nodes[i][:, None]                 # (K, n_particles)
        upper = log_ndtr((g_hi - mean_g) / sg[i])
        if np.isfinite(lo):
            g_lo = lo + kG * med[i]
            lower = log_ndtr((g_lo - mean_g) / sg[i])
            log_p = upper + np.log1p(-np.exp(np.minimum(lower - upper, -1e-12)))
        else:
            log_p = upper
        den = logsumexp(lw[None, :] + log_p) - np.log(len(gd))
        out[i] = num - den
    if outlier_eps > 0:
        out = np.logaddexp(np.log1p(-outlier_eps) + out, np.log(outlier_eps * outlier_density))
    return out


def gauss_loglike(stars: pd.DataFrame, particles, window: str, correlated: bool,
                  dmu: float = 0.0, av_col: str = "av_median", sig_col: str = "av_halfwidth68") -> np.ndarray:
    """Gaussian A_V error at a fixed A_V: independent (WP4) or correlated along r."""
    wp = _window_particles(particles, window, dmu)
    if wp is None:
        return np.full(len(stars), -50.0)
    cm, gm, lw = wp
    kG, kc = stars.attrs["kG"], stars.attrs["kc"]
    a = stars[av_col].to_numpy()
    sa = np.nan_to_num(stars[sig_col].to_numpy(), nan=0.0)
    g = stars["g_obs"].to_numpy() - kG * a
    c = stars["c_obs"].to_numpy() - kc * a
    sg, sc = stars["sig_g"].to_numpy(), stars["sig_c"].to_numpy()
    out = np.empty(len(stars))
    for i in range(len(stars)):
        c11 = sg[i] ** 2 + (kG * sa[i]) ** 2
        c22 = sc[i] ** 2 + (kc * sa[i]) ** 2
        c12 = kG * kc * sa[i] ** 2 if correlated else 0.0
        det = c11 * c22 - c12 ** 2
        rg, rc = g[i] - gm, c[i] - cm
        q = (c22 * rg ** 2 - 2 * c12 * rg * rc + c11 * rc ** 2) / det
        out[i] = logsumexp(lw - 0.5 * q) - 0.5 * np.log(4 * np.pi ** 2 * det)
    return out


def particles_by_age(iso_family: pd.DataFrame, f_bin: float, drop_pms: bool = False):
    ages = np.sort(iso_family["age_Myr"].unique())
    out = []
    for a in ages:
        rows = iso_family[np.isclose(iso_family["age_Myr"], a)]
        if drop_pms:
            rows = rows[rows["phase"] > -1]          # MIST only; PARSEC labels are unreliable
        out.append(W4.build_model_particles(rows, f_bin))
    return ages, out


def posterior(ages: np.ndarray, total: np.ndarray, n_stars: int) -> dict:
    """WP4's posterior convention (wp4_common.posterior_from_loglike) plus
    95 % bounds, the WP4 measurability gate and grid_railed."""
    from wp4_fit_ages import grid_railed, exclusion_reason
    import issue20_common as I20
    p = W4.posterior_from_loglike(ages, total)
    q = I20.posterior_summary(ages, total)
    railed = grid_railed(p["map"], ages)
    reason = exclusion_reason(n_stars, p["map"], railed)
    return {"age_map": p["map"], "age_lo68": p["lo68"], "age_hi68": p["hi68"],
            "age_lo90": p["lo90"], "age_hi90": p["hi90"], "age_mean": p["mean"],
            "age_median": q["median"], "age_lo95": q["lo95"], "age_hi95": q["hi95"],
            "grid_railed": bool(railed), "measurable": reason == "", "exclusion_reason": reason}


# ----------------------------------------------------------------- injections
def synthetic(stars_all: pd.DataFrame, particles_t, window: str, n_target: int,
              rng: np.random.Generator) -> pd.DataFrame:
    """Synthetic subgroup at one true age.  Each synthetic star takes a random
    real star's photometric errors and A_V posterior; its TRUE A_V is a draw
    from that posterior (so the errors are correlated along r exactly as M1
    assumes), its intrinsic photometry a draw from the model particles.  The
    window is applied on the posterior-median M_G0, as for the real data, and
    n_target stars are kept (the real count in the window)."""
    colour, mg, weight = particles_t
    kG, kc = stars_all.attrs["kG"], stars_all.attrs["kc"]
    keep, total = [], 0
    while total < n_target:
        n = 4 * n_target + 200
        j = rng.choice(len(weight), n, p=weight / weight.sum())
        donor = stars_all.iloc[rng.integers(0, len(stars_all), n)].reset_index(drop=True)
        nodes = np.stack(donor["_nodes"].to_numpy())
        a_true = nodes[np.arange(n), rng.integers(0, nodes.shape[1], n)]
        a_true = a_true + rng.uniform(-0.5, 0.5, n) * np.median(np.diff(nodes, axis=1), axis=1)
        s = donor.copy()
        s["g_obs"] = mg[j] + kG * a_true + rng.normal(0, 1, n) * s["sig_g"].to_numpy()
        s["c_obs"] = colour[j] + kc * a_true + rng.normal(0, 1, n) * s["sig_c"].to_numpy()
        s["MG0_median"] = s["g_obs"] - kG * s["av_median"]
        s["P"] = 1.0
        s.attrs.update(stars_all.attrs)
        sel = in_window(s, window)
        keep.append(sel)
        total += len(sel)
    out = pd.concat(keep, ignore_index=True).iloc[:n_target].copy()
    out.attrs.update(stars_all.attrs)
    return out


# --------------------------------------------------------------- GA3 bound
SIG_LOGTE_GA3 = 0.03
SIG_MG0_GA3 = 0.40
CHI_GA3 = 2.0


def minimum_initial_mass(logte: float, mg0: float, iso_family: pd.DataFrame) -> float:
    """Smallest initial mass of ANY tabulated isochrone point (any native age,
    any phase) within chi <= 2 of the star (wp4_anchors_hrd metric).  A lower
    bound on the star's initial mass that assumes no age."""
    d2 = (((iso_family["logTe"] - logte) / SIG_LOGTE_GA3) ** 2
          + ((iso_family["G0"] - mg0) / SIG_MG0_GA3) ** 2)
    near = iso_family.loc[d2 <= CHI_GA3 ** 2, "Mini"]
    return float(near.min()) if len(near) else np.nan


# ------------------------------------------------------------- joint fallback
def bimodal(ages: np.ndarray, total: np.ndarray) -> bool:
    """Two local maxima of the posterior density (PCHIP in log age, as the
    summaries) whose separating minimum is < 0.5 x the smaller peak and whose
    smaller mode holds >= 10 % of the probability."""
    from scipy.interpolate import PchipInterpolator
    la = np.log10(ages)
    laf = np.linspace(la.min(), la.max(), 4001)
    f = PchipInterpolator(la, total - np.max(total))(laf)
    p = np.exp(f)
    p /= np.trapezoid(p, laf)
    peaks = [i for i in range(1, len(p) - 1) if p[i] >= p[i - 1] and p[i] > p[i + 1]]
    if len(peaks) < 2:
        return False
    for a, b in zip(peaks[:-1], peaks[1:]):
        trough = int(a + np.argmin(p[a:b + 1]))
        small = min(p[a], p[b])
        if p[trough] < 0.5 * small:
            left = np.trapezoid(p[:trough + 1], laf[:trough + 1])
            if min(left, 1 - left) >= 0.10:
                return True
    return False
