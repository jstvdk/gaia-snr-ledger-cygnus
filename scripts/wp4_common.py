#!/usr/bin/env python3
"""Shared machinery for WP4 (subgroup ages & per-star masses).

Design (every knob is labelled with its class from CUTS_AND_THRESHOLDS.md §1):

  * Age is read off the de-reddened CMD by a *forward-model likelihood*, not by
    eye and not by nearest-isochrone chi2.  For each trial age we synthesise the
    predicted CMD density of a single-age population (single stars + an
    unresolved-binary component) and evaluate the per-star marginal likelihood,
    weighting each star by its WP2 membership probability.

  * PARSEC and MIST are carried as parallel branches (Class E) and NEVER
    averaged.  R_V = 3.0/3.1/3.5 (Class E) enter through the de-reddened data.
    Binary fraction f_bin = 0.3/0.4/0.5 (Class D/E) is scanned.  Distance is
    fixed at the single-population posterior 1.6245 kpc and its +/-0.045 kpc is
    propagated as a coherent distance-modulus shift (Class A).

  * Two age indicators are fit independently on disjoint magnitude windows:
      - upper main sequence  (bright, M_G0 <= UMS_FAINT_EDGE)
      - pre-main-sequence turn-on / lower CMD (faint, M_G0 >= PMS_BRIGHT_EDGE)
    Agreement between them is the internal cross-check (plan WP4 step 2).

Nothing here fetches data; the isochrone grids are the frozen WP3 files and are
shared verbatim with WP4 (plan requirement).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

from wp3_extinction_law import band_coefficients, R_V_BRANCHES

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
TABLES = ROOT / "tables"

# ---------------------------------------------------------------------------
# Distance posterior (single population; WP2/WP3).  Class A.
# ---------------------------------------------------------------------------
D_KPC = 1.6245
D_KPC_ERR = 0.045
DIST_MODULUS = 5.0 * np.log10(D_KPC * 1e3) - 5.0                     # ~11.054
DIST_MODULUS_ERR = 5.0 / np.log(10) * (D_KPC_ERR / D_KPC)            # ~0.060 mag

SUBGROUPS = ["CygOB2-A", "CygOB2-B", "CygOB2-C"]

# ---------------------------------------------------------------------------
# Fit configuration constants (each labelled by class; logged to provenance).
# ---------------------------------------------------------------------------
IMF_SLOPE = 2.3          # Class B (Kroupa/Salpeter) - along-isochrone density weight only
F_BIN_BRANCHES = [0.3, 0.4, 0.5]   # Class D/E - unresolved-binary fraction scan
Q_MIN = 0.10             # Class B - minimum mass ratio for the binary component
SIGMA_INT = 0.03         # Class D - intrinsic model/isochrone scatter floor (mag)
MAG_FLOOR = 0.02         # Class D - photometric floor added in quadrature (mag)

# Indicator windows in de-reddened absolute G (Class C/D - set by data extent,
# see wp4 provenance; the faint edge of the *observable* sample sits near
# M_G0 ~ +3 at this distance+extinction, so the PMS window is intrinsically
# starved and an indicator is only reported where it holds >= N_MIN_INDICATOR
# stars).  The minimum is a Class-B reporting convention registered in
# CUTS_AND_THRESHOLDS.md; it is not tuned against the recovered ages.
UMS_FAINT_EDGE = 1.5     # upper-MS window: M_G0 <= 1.5
PMS_BRIGHT_EDGE = 2.5    # PMS/lower window: M_G0 >= 2.5
N_MIN_INDICATOR = 15     # Class B: minimum star count for a measurable age indicator

# CMD colour used for the fit: (BP-RP)0.  Fit plane is (colour, M_G0).
COLOUR_LO, COLOUR_HI = -1.6, 2.6   # global CMD colour box for model normalisation


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_members() -> pd.DataFrame:
    """WP3 extinction catalogue joined to the A/B/C subgroup labels (WP2)."""
    ext = pd.read_parquet(PROC / "wp3_extinction.parquet")
    lab = pd.read_parquet(TABLES / "wp2_subgroup_labels.parquet")
    if "subgroup" not in lab or "subgroup_label" in lab:
        raise RuntimeError("subgroup sidecar must expose only the canonical 'subgroup' column")
    ext = ext.drop(columns=["subgroup", "subgroup_label"], errors="ignore")
    m = ext.merge(
        lab[["source_id", "subgroup"]],
        on="source_id", how="left",
    )
    m["subgroup"] = m["subgroup"].fillna("unassigned")
    return m


def load_isochrones():
    par = pd.read_parquet(PROC / "wp3_isochrones_parsec.parquet")
    mist = pd.read_parquet(PROC / "wp3_isochrones_mist.parquet")
    return {"PARSEC": par, "MIST": mist}


def isochrone_ms_pms_flag(iso: pd.DataFrame, family: str) -> np.ndarray:
    """Boolean 'is pre-main-sequence' per isochrone row.

    PARSEC 'label': 0 = PMS, 1 = MS, >=2 = post-MS.
    MIST   'phase': -1 = PMS, 0 = MS, >=2 = post-MS.
    (Used only for figure annotation; the fit uses the full isochrone.)
    """
    if family == "PARSEC":
        return iso["label"].values == 0
    return iso["phase"].values == -1


# ---------------------------------------------------------------------------
# De-reddened data + error propagation for one R_V branch
# ---------------------------------------------------------------------------
def branch_photometry(m: pd.DataFrame, rv: float) -> pd.DataFrame:
    """Return per-star fit inputs for one R_V branch: colour, M_G0, and their
    1-sigma errors including the extinction-residual contribution.

    sigma_MG0^2   = G_err^2 + (kG * av_err)^2 + floor^2
    sigma_col^2   = BP_err^2 + RP_err^2 + ((kBP-kRP)*av_err)^2 + floor^2

    The coherent distance-modulus uncertainty is NOT folded in here (it is a
    shared shift, propagated by refitting at mu +/- sigma_mu, not per-star
    scatter).
    """
    tag = f"rv{rv:.1f}".replace(".", ".")  # column suffix uses e.g. 'rv3.1'
    suf = f"_rv{rv:.1f}"
    k = band_coefficients(rv)
    av_err = m[f"av_err_rv{rv:.1f}"].values

    col = m[f"BPRP0{suf}"].values
    mg = m[f"G0_abs{suf}"].values

    g_err = np.nan_to_num(m["G_err"].values, nan=0.02)
    bp_err = np.nan_to_num(m["BP_err"].values, nan=0.02)
    rp_err = np.nan_to_num(m["RP_err"].values, nan=0.02)
    av_err = np.nan_to_num(av_err, nan=0.0)

    sig_mg = np.sqrt(g_err ** 2 + (k["G"] * av_err) ** 2 + MAG_FLOOR ** 2)
    sig_col = np.sqrt(bp_err ** 2 + rp_err ** 2
                      + ((k["BP"] - k["RP"]) * av_err) ** 2 + MAG_FLOOR ** 2)

    out = pd.DataFrame({
        "source_id": m["source_id"].values,
        "subgroup": m["subgroup"].values,
        "P": m["membership_probability"].values,
        "colour": col,
        "MG0": mg,
        "sig_col": sig_col,
        "sig_mg": sig_mg,
        "av": m[f"av_rv{rv:.1f}"].values,
        "av_method": m["av_method"].values,
    })
    return out


# ---------------------------------------------------------------------------
# Forward synthetic-CMD model for a single age
# ---------------------------------------------------------------------------
def _combine_mag(m1, m2):
    """Flux-add two absolute magnitudes -> combined absolute magnitude."""
    return -2.5 * np.log10(10.0 ** (-0.4 * m1) + 10.0 ** (-0.4 * m2))


def build_model_particles(iso_age: pd.DataFrame, f_bin: float,
                          n_mass: int = 400, n_q: int = 24):
    """Synthesise the predicted CMD point-cloud for a single-age population.

    Returns (colour[j], MG[j], weight[j]) for a weighted set of 'particles':
      - a single-star sequence sampled along the isochrone in log-mass, weighted
        by the IMF (Class B slope IMF_SLOPE);
      - an unresolved-binary sequence: each primary paired with secondaries on a
        mass-ratio grid q in [Q_MIN, 1], q ~ uniform (Class B), fluxes added.

    Weights sum to 1; single component carries (1-f_bin), binary carries f_bin.
    """
    iso = iso_age.sort_values("Mini")
    mini = iso["Mini"].values
    g = iso["G0"].values
    bp = iso["BP0"].values
    rp = iso["RP0"].values
    good = np.isfinite(mini) & np.isfinite(g) & np.isfinite(bp) & np.isfinite(rp)
    mini, g, bp, rp = mini[good], g[good], bp[good], rp[good]
    # unique/monotone mass for interpolation
    mini, idx = np.unique(mini, return_index=True)
    g, bp, rp = g[idx], bp[idx], rp[idx]
    if mini.size < 5:
        return None

    # log-spaced primary-mass grid across the isochrone's mass range
    lm = np.linspace(np.log10(mini.min()), np.log10(mini.max()), n_mass)
    mgrid = 10.0 ** lm
    G = np.interp(mgrid, mini, g)
    BP = np.interp(mgrid, mini, bp)
    RP = np.interp(mgrid, mini, rp)

    # IMF weight per primary (density in mass): xi(m) dm ~ m^-slope * dm.
    # with a log grid, dm ~ m dln(m), so weight ~ m^(1-slope).
    w_primary = mgrid ** (1.0 - IMF_SLOPE)
    w_primary /= w_primary.sum()

    # --- single-star particles ---
    s_col = BP - RP
    s_mg = G
    s_w = (1.0 - f_bin) * w_primary

    # --- binary particles ---
    q = np.linspace(Q_MIN, 1.0, n_q)
    wq = np.ones_like(q) / q.size          # q ~ uniform
    m2 = np.outer(mgrid, q)                 # (n_mass, n_q)
    # secondary photometry by interpolation (below grid min -> clip to faint)
    G2 = np.interp(m2, mini, g)
    BP2 = np.interp(m2, mini, bp)
    RP2 = np.interp(m2, mini, rp)
    Gc = _combine_mag(G[:, None], G2)
    BPc = _combine_mag(BP[:, None], BP2)
    RPc = _combine_mag(RP[:, None], RP2)
    b_col = (BPc - RPc).ravel()
    b_mg = Gc.ravel()
    b_w = (f_bin * (w_primary[:, None] * wq[None, :])).ravel()

    colour = np.concatenate([s_col, b_col])
    mg = np.concatenate([s_mg, b_mg])
    weight = np.concatenate([s_w, b_w])
    weight /= weight.sum()
    return colour, mg, weight


# ---------------------------------------------------------------------------
# Per-star marginal likelihood and age posterior
# ---------------------------------------------------------------------------
def _window_mask_stars(fit: pd.DataFrame, window: str):
    if window == "ums":
        return fit["MG0"].values <= UMS_FAINT_EDGE
    elif window == "pms":
        return fit["MG0"].values >= PMS_BRIGHT_EDGE
    raise ValueError(window)


def _window_mask_model(colour, mg, window):
    if window == "ums":
        return mg <= UMS_FAINT_EDGE
    return mg >= PMS_BRIGHT_EDGE


def star_loglike(fit_sub: pd.DataFrame, particles, window: str,
                 dmu: float = 0.0) -> np.ndarray:
    """Per-star log-likelihood under one age's model, restricted+renormalised to
    the indicator window.  `dmu` shifts the *model* magnitudes to emulate a
    distance-modulus offset (distance branch).  Returns array over the stars in
    fit_sub (already window-masked)."""
    colour, mg, weight = particles
    mg = mg + dmu
    mm = _window_mask_model(colour, mg, window)
    if mm.sum() < 3:
        return np.full(len(fit_sub), -50.0)
    cm, gm, wm = colour[mm], mg[mm], weight[mm]
    wm = wm / wm.sum()

    ci = fit_sub["colour"].values[:, None]
    mi = fit_sub["MG0"].values[:, None]
    sc = np.sqrt(fit_sub["sig_col"].values[:, None] ** 2 + SIGMA_INT ** 2)
    sm = np.sqrt(fit_sub["sig_mg"].values[:, None] ** 2 + SIGMA_INT ** 2)

    # 2D independent-Gaussian mixture over model particles
    dc = (ci - cm[None, :]) / sc
    dm = (mi - gm[None, :]) / sm
    norm = 1.0 / (2.0 * np.pi * sc * sm)
    comp = wm[None, :] * norm * np.exp(-0.5 * (dc ** 2 + dm ** 2))
    dens = comp.sum(axis=1)
    dens = np.clip(dens, 1e-300, None)
    return np.log(dens)


def fit_age_grid(fit_sub: pd.DataFrame, iso_family: pd.DataFrame, family: str,
                 f_bin: float, window: str, dmu: float = 0.0):
    """Return (ages_Myr, loglike) over the isochrone age grid for one branch."""
    ages = np.sort(iso_family["age_Myr"].unique())
    smask = _window_mask_stars(fit_sub, window)
    stars = fit_sub[smask & np.isfinite(fit_sub["colour"]) & np.isfinite(fit_sub["MG0"])]
    P = stars["P"].values
    lls = []
    for a in ages:
        iso_age = iso_family[np.isclose(iso_family["age_Myr"], a)]
        parts = build_model_particles(iso_age, f_bin)
        if parts is None:
            lls.append(-1e9)
            continue
        ll = star_loglike(stars, parts, window, dmu=dmu)
        lls.append(float(np.sum(P * ll)))       # membership-weighted total loglike
    return ages, np.array(lls), int(len(stars))


def posterior_from_loglike(ages: np.ndarray, loglike: np.ndarray, n_fine: int = 2001):
    """Smooth the age-grid log-likelihood, form a posterior (uniform prior in
    log age since the grid is log-spaced), and return MAP + 68/90% credible
    intervals on a fine log-age grid."""
    la = np.log10(ages)
    good = np.isfinite(loglike) & (loglike > -1e8)
    if good.sum() < 4:
        return dict(map=np.nan, lo68=np.nan, hi68=np.nan, lo90=np.nan,
                    hi90=np.nan, mean=np.nan, ages_fine=None, post=None)
    la, ll = la[good], loglike[good]
    # cubic-ish interpolation of loglike in log-age
    laf = np.linspace(la.min(), la.max(), n_fine)
    llf = np.interp(laf, la, ll)              # linear; grid is fine enough
    try:
        from scipy.interpolate import PchipInterpolator
        llf = PchipInterpolator(la, ll)(laf)
    except Exception:
        pass
    llf -= llf.max()
    post = np.exp(llf)
    post /= np.trapezoid(post, laf)
    agef = 10.0 ** laf
    # MAP
    imap = int(np.argmax(post))
    amap = agef[imap]
    # credible intervals from the cumulative
    cdf = np.concatenate([[0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(laf))])
    cdf /= cdf[-1]
    def q(p):
        return 10.0 ** np.interp(p, cdf, laf)
    mean = 10.0 ** np.trapezoid(laf * post, laf)
    return dict(map=float(amap), lo68=float(q(0.16)), hi68=float(q(0.84)),
                lo90=float(q(0.05)), hi90=float(q(0.95)), mean=float(mean),
                ages_fine=agef, post=post)
