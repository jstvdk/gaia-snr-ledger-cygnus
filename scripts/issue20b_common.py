#!/usr/bin/env python3
"""Issue #20 Phase A' -- shared method: three independent age tests for A, B, C.

Hashed into provenance/issue20b_prereg.json before any test runs; the drivers
refuse to run against a changed copy unless a deviation records it.

  Test a  free-extinction CMD + near-IR likelihood.  Each star's A_V is a free
          parameter (flat prior), integrated analytically against its J-Ks
          colour, so the WP3 extinction map is never used.  The star's
          (G, BP-RP) is then compared with every IMF-weighted model particle of
          the trial age; the model particle supplies the intrinsic colours.
  Test b  Gaia DR3 ESP-HS hot-star HRD.  T_eff and A_G from Gaia's own XP
          spectra (zero points calibrated on the spectroscopic anchors), then
          the conditional HRD likelihood of issue #20 with an outlier term.
  Test c  the issue #20 A3 spectroscopic-HRD likelihood with an outlier term.

Nothing here writes a file.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from itertools import combinations

from scipy.interpolate import PchipInterpolator
from scipy.special import logsumexp

import issue20_common as I20
import wp5_common as w
from wp3_extinction_law import band_coefficients
from wp4_common import DIST_MODULUS

PREREG_PATH = w.PROVENANCE / "issue20b_prereg.json"
DEVIATIONS_PATH = w.PROVENANCE / "issue20b_deviations.json"
THIS_FILE = Path(__file__).resolve()

SUBGROUPS = I20.SUBGROUPS
FAMILIES = I20.FAMILIES
R_V_BRANCHES = I20.R_V_BRANCHES
F_BIN = 0.4
UMS_EDGE = 1.5                       # WP4 upper-MS window, M_G0 <= 1.5

INPUTS = {
    **{k: v for k, v in I20.INPUTS.items() if k in (
        "anchor_hrd", "spectroscopic_anchors", "extinction_v5", "subgroup_labels",
        "isochrones_parsec", "isochrones_mist", "age_posteriors_v5")},
    "esphs": "data/raw/gaia/issue20b_gaia_dr3_astrophysical_parameters.csv",
    "issue20_outcome": "provenance/issue20_phase_a_outcome.json",
    "issue20_per_star": "tables/issue20_hrd_per_star.csv",
}


class FrozenInputMoved(RuntimeError):
    pass


def prereg() -> dict:
    if not PREREG_PATH.exists():
        raise FileNotFoundError("run scripts/issue20b_prereg.py first")
    return json.loads(PREREG_PATH.read_text())


def check_method_unchanged() -> dict:
    rec = prereg()["method_code"]
    out = {}
    for rel, frozen_hash in rec.items():
        now = w.sha256(w.ROOT / rel)
        if now == frozen_hash:
            out[rel] = {"sha256": now, "deviation": None}
            continue
        devs = json.loads(DEVIATIONS_PATH.read_text()) if DEVIATIONS_PATH.exists() else {}
        hit = [d for d in devs.get("method_code_changes", [])
               if d.get("file") == rel and d.get("new_sha256") == now]
        if not hit:
            raise FrozenInputMoved(f"{rel} changed since the Phase A' preregistration "
                                   "and no deviation records it")
        out[rel] = {"sha256": now, "deviation": hit[0]}
    return out


def frozen(name: str) -> Path:
    entry = prereg()["consumed_inputs"][name]
    path = w.ROOT / entry["path"]
    if w.sha256(path) != entry["sha256"]:
        raise FrozenInputMoved(f"{entry['path']} changed since the Phase A' preregistration")
    return path


# ------------------------------------------------------------------ shared
EPS_OUTLIER = 0.05                   # outlier fraction, tests b and c
OUTLIER_LOGTE = (4.10, 4.75)         # uniform outlier density in log T_eff
MIN_STARS = 10                       # a test is valid for a subgroup only with >= 10 stars


def with_outliers(cond_ll: np.ndarray, eps: float = EPS_OUTLIER) -> np.ndarray:
    """ln[(1-eps) p_iso(logTe | M) + eps U(logTe)], per star and age."""
    lo, hi = OUTLIER_LOGTE
    return np.logaddexp(np.log1p(-eps) + cond_ll, np.log(eps) - np.log(hi - lo))


def subgroup_posterior(ages, ll, weights) -> dict:
    return I20.posterior_summary(ages, weights @ ll)


# ------------------------------------------------------------ test a (NIR)
JKS0_HOT = -0.20          # intrinsic J-Ks of O-B3 stars, used for calibration and selection only
SIG_JKS_INT = 0.03        # model intrinsic near-IR colour scatter
SIG_INT = 0.03            # as wp4_common.SIGMA_INT
MAG_FLOOR = 0.02          # as wp4_common.MAG_FLOOR
CMD_DMAG_STEP = 0.05
CMD_DLOGTE_STEP = 0.01
CMD_M_MIN = 1.0
NIR_QUAL = set("AB")


def nir_quality(frame: pd.DataFrame) -> np.ndarray:
    q = frame["ph_qual"].astype(str)
    ok = q.str.len().eq(3) & q.str[0].isin(NIR_QUAL) & q.str[2].isin(NIR_QUAL)
    return (ok & frame["J"].notna() & frame["Ks"].notna()
            & frame["J_err"].notna() & frame["Ks_err"].notna()).to_numpy()


def nir_scale(anchors: pd.DataFrame, ext: pd.DataFrame, rv: float) -> dict:
    """Calibrate E(J-Ks)/A_V on the spectroscopic anchors (O to B3, any class):
    s = median of [(J-Ks + 0.20)/(k_J - k_Ks)] / A_V,spec.  Validity: >= 30
    anchors, |s - 1| <= 0.3, robust scatter of (A_V,NIR/s - A_V,spec) <= 0.5 mag."""
    k = band_coefficients(rv)
    m = anchors.merge(ext[["source_id", "J", "Ks", "J_err", "Ks_err", "ph_qual",
                           f"av_rv{rv:.1f}"]], on="source_id")
    parsed = m.spectral_type.map(I20.parse_spectral_type)
    tn = np.array([p["type_number"] for p in parsed], float)
    keep = (nir_quality(m) & ~m.extreme_hot.to_numpy() & np.isfinite(tn) & (tn <= 13.0)
            & m[f"av_rv{rv:.1f}"].notna().to_numpy())
    m = m[keep]
    av_nir = (m.J - m.Ks - JKS0_HOT) / (k["J"] - k["Ks"])
    ratio = av_nir / m[f"av_rv{rv:.1f}"]
    s = float(np.median(ratio))
    resid = av_nir / s - m[f"av_rv{rv:.1f}"]
    sig = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
    return {"R_V": rv, "n": int(len(m)), "scale": s, "robust_sigma_mag": sig,
            "median_resid_mag": float(np.median(resid)),
            "valid": bool(len(m) >= 30 and abs(s - 1.0) <= 0.3 and sig <= 0.5)}


def _sample_bands(iso_age: pd.DataFrame, family: str, bands: list[str], m_min: float,
                  dlogte: float, dmag: float):
    """Fine IMF-weighted sampling along the isochrone (the issue20_common rule,
    including the MIST gap rule), carrying several bands."""
    mini = iso_age["Mini"].to_numpy(float)
    logte = iso_age["logTe"].to_numpy(float)
    vals = {b: iso_age[b].to_numpy(float) for b in bands}
    g = vals["G0"]
    phase = (iso_age["phase"].to_numpy(float) if family == "MIST" else np.zeros(len(mini)))
    out_m, out_v, out_w = [], {b: [] for b in bands}, []
    for i in range(len(mini) - 1):
        m_a, m_b = mini[i], mini[i + 1]
        if m_b <= m_min:
            continue
        if (family == "MIST" and phase[i] != phase[i + 1]
                and (abs(logte[i + 1] - logte[i]) > I20.MIST_GAP_DLOGTE
                     or abs(g[i + 1] - g[i]) > I20.MIST_GAP_DMAG)):
            continue
        lo = max(m_a, m_min)
        f_lo = (lo - m_a) / (m_b - m_a)
        n = int(np.clip(np.ceil(max(abs(logte[i + 1] - logte[i]) * (1 - f_lo) / dlogte,
                                    abs(g[i + 1] - g[i]) * (1 - f_lo) / dmag, 1.0)),
                        1, I20.MAX_SUBSTEPS))
        edges = lo + (m_b - lo) * np.arange(n + 1) / n
        mid = 0.5 * (edges[1:] + edges[:-1])
        frac = (mid - m_a) / (m_b - m_a)
        out_m.append(mid)
        for b in bands:
            out_v[b].append(vals[b][i] + frac * (vals[b][i + 1] - vals[b][i]))
        out_w.append(I20._imf_integral(edges[:-1], edges[1:]))
    return (np.concatenate(out_m), {b: np.concatenate(v) for b, v in out_v.items()},
            np.concatenate(out_w))


def cmd_particles(iso_age: pd.DataFrame, family: str, f_bin: float = F_BIN) -> dict:
    """Particles carrying G0, BP0-RP0 and J0-Ks0, singles + unresolved binaries
    (q ~ U(0.1, 1), 10 values, all bands flux-added)."""
    bands = ["G0", "BP0", "RP0", "J0", "Ks0"]
    m, v, wt = _sample_bands(iso_age, family, bands, CMD_M_MIN, CMD_DLOGTE_STEP, CMD_DMAG_STEP)
    mini = iso_age["Mini"].to_numpy(float)
    q = np.linspace(I20.Q_MIN, 1.0, I20.N_Q)
    comb = {}
    for b in bands:
        sec = np.interp(np.outer(m, q), mini, iso_age[b].to_numpy(float), left=99.0)
        comb[b] = np.concatenate([v[b], I20._flux_add(v[b][:, None], sec).ravel()])
    wts = np.concatenate([(1 - f_bin) * wt, np.repeat(f_bin * wt / I20.N_Q, I20.N_Q)])
    return {"G": comb["G0"], "c": comb["BP0"] - comb["RP0"], "x": comb["J0"] - comb["Ks0"],
            "w": wts, "log_w": np.log(wts)}


def nir_inputs(frame: pd.DataFrame, rv: float, scale: float) -> pd.DataFrame:
    """Observables of test a.  The star sample is fixed across ages: M_G0 from
    the calibrated near-IR excess with (J-Ks)_0 = -0.20 must be <= 1.5."""
    k = band_coefficients(rv)
    kx = (k["J"] - k["Ks"]) / scale
    f = frame[nir_quality(frame)].copy()
    f["g"] = f["G"] - DIST_MODULUS
    f["c"] = f["BP"] - f["RP"]
    f["x"] = f["J"] - f["Ks"]
    f["sig_g"] = np.sqrt(np.nan_to_num(f["G_err"], nan=0.02) ** 2 + MAG_FLOOR ** 2
                         + SIG_INT ** 2 + I20.SIG_DEPTH_MAG ** 2)
    f["sig_c"] = np.sqrt(np.nan_to_num(f["BP_err"], nan=0.02) ** 2
                         + np.nan_to_num(f["RP_err"], nan=0.02) ** 2 + MAG_FLOOR ** 2 + SIG_INT ** 2)
    f["sig_x"] = np.sqrt(f["J_err"] ** 2 + f["Ks_err"] ** 2 + SIG_JKS_INT ** 2)
    f["mg0_select"] = f["g"] - k["G"] * (f["x"] - JKS0_HOT) / kx
    f = f[(f["mg0_select"] <= UMS_EDGE) & np.isfinite(f["c"]) & np.isfinite(f["g"])]
    f.attrs.update(kG=k["G"], kc=k["BP"] - k["RP"], kx=kx)
    return f


def cmd_loglike(stars: pd.DataFrame, parts: dict, mode: str, rv: float) -> np.ndarray:
    """Per-star ln L at one age.  mode 'nir': A_V free, integrated against J-Ks
    (exact for Gaussian errors and a flat A_V prior).  mode 'wp3': A_V fixed at
    the repair_v5 WP3 value +- av_err (control).  Window-normalised over model
    particles with G0 <= 1.5, as wp4_common.star_loglike."""
    kG, kc, kx = stars.attrs["kG"], stars.attrs["kc"], stars.attrs["kx"]
    inside = parts["G"] <= UMS_EDGE
    gm, cm, xm, lw = parts["G"][inside], parts["c"][inside], parts["x"][inside], parts["log_w"][inside]
    log_norm = logsumexp(lw)
    out = np.empty(len(stars))
    g, c, x = stars["g"].to_numpy(), stars["c"].to_numpy(), stars["x"].to_numpy()
    sg, sc, sx = stars["sig_g"].to_numpy(), stars["sig_c"].to_numpy(), stars["sig_x"].to_numpy()
    if mode == "wp3":
        a_fix = stars[f"av_rv{rv:.1f}"].to_numpy()
        sa = np.nan_to_num(stars[f"av_err_rv{rv:.1f}"].to_numpy(), nan=0.0)
    for s in range(0, len(stars), 16):
        sl = slice(s, s + 16)
        if mode == "nir":
            a_hat = (x[sl, None] - xm[None, :]) / kx
            sa2 = (sx[sl] / kx) ** 2
        else:
            a_hat = np.broadcast_to(a_fix[sl, None], (len(g[sl]), len(gm)))
            sa2 = sa[sl] ** 2
        rg = g[sl, None] - kG * a_hat - gm[None, :]
        rc = c[sl, None] - kc * a_hat - cm[None, :]
        c11 = sg[sl] ** 2 + kG ** 2 * sa2
        c22 = sc[sl] ** 2 + kc ** 2 * sa2
        c12 = kG * kc * sa2
        det = c11 * c22 - c12 ** 2
        q = (c22[:, None] * rg ** 2 - 2 * c12[:, None] * rg * rc + c11[:, None] * rc ** 2) / det[:, None]
        out[sl] = (logsumexp(lw[None, :] - 0.5 * q, axis=1) - log_norm
                   - 0.5 * np.log((2 * np.pi) ** 2 * det))
    return out


# ------------------------------------------------------------ test b (ESP-HS)
ESPHS_TYPES = ("O", "B")


def esphs_calibration(anchors: pd.DataFrame, esphs: pd.DataFrame, ext: pd.DataFrame,
                      rv: float) -> dict:
    """Zero points of ESP-HS log T_eff and A_G against the anchors (O/B primary,
    not extreme hot, luminosity class III-V or unclassified, table T_eff, and
    the anchors' WP3 intrinsic-colour A_V).  Validity: >= 20 anchors, robust
    sigma(log T) <= 0.08 dex, robust sigma(A_G) <= 0.6 mag."""
    k = band_coefficients(rv)
    m = anchors.merge(esphs, on="source_id").merge(
        ext[["source_id", f"av_rv{rv:.1f}"]], on="source_id")
    parsed = m.spectral_type.map(I20.parse_spectral_type)
    m["letter"] = [p["letter"] for p in parsed]
    m["lum"] = [p["lum_class"] for p in parsed]
    m = m[m.teff_esphs.notna() & m.ag_esphs.notna() & ~m.extreme_hot & m.letter.notna()
          & ~m.lum.isin(["I", "II"]) & m.spectraltype_esphs.isin(ESPHS_TYPES)
          & m[f"av_rv{rv:.1f}"].notna()]
    d_t = np.log10(m.teff_esphs) - m.logTe_spec
    d_a = m.ag_esphs - k["G"] * m[f"av_rv{rv:.1f}"]

    def rob(v):
        return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else np.nan
    out = {"R_V": rv, "n": int(len(m)),
           "zp_logTe": float(np.median(d_t)) if len(m) else np.nan, "sig_logTe": rob(d_t),
           "zp_AG": float(np.median(d_a)) if len(m) else np.nan, "sig_AG": rob(d_a)}
    out["valid"] = bool(out["n"] >= 20 and out["sig_logTe"] <= 0.08 and out["sig_AG"] <= 0.6)
    return out


def esphs_inputs(members: pd.DataFrame, esphs: pd.DataFrame, cal: dict) -> pd.DataFrame:
    m = members.merge(esphs, on="source_id")
    m = m[m.teff_esphs.notna() & m.ag_esphs.notna() & m.spectraltype_esphs.isin(ESPHS_TYPES)].copy()
    m["logTe"] = np.log10(m.teff_esphs) - cal["zp_logTe"]
    unc_t = np.nan_to_num(m.teff_esphs_uncertainty / m.teff_esphs / np.log(10), nan=0.0)
    m["sig_logTe"] = np.sqrt(np.maximum(unc_t ** 2 + cal["sig_logTe"] ** 2, 0.03 ** 2))
    m["MG0"] = m.G - DIST_MODULUS - (m.ag_esphs - cal["zp_AG"])
    m["sig_MG0"] = np.sqrt(np.nan_to_num(m.ag_esphs_uncertainty, nan=0.0) ** 2 + cal["sig_AG"] ** 2
                           + I20.SIG_DEPTH_MAG ** 2 + np.nan_to_num(m.G_err, nan=0.02) ** 2)
    return m[m.MG0 <= UMS_EDGE]


# ------------------------------------------------------------ agreement rule
def overlap68(p, q) -> bool:
    return not (p["hi68"] < q["lo68"] or q["hi68"] < p["lo68"])


def adopt(results: dict) -> dict:
    """results: test -> posterior summary (with 'valid').  Adopt the unique
    largest set (size >= 2) of valid tests whose 68 % intervals pairwise
    overlap; otherwise NOT RESOLVED."""
    valid = [t for t, r in results.items() if r.get("valid")]
    best, sets = 0, []
    for size in (3, 2):
        for combo in combinations(valid, size):
            if all(overlap68(results[a], results[b]) for a, b in combinations(combo, 2)):
                sets.append(combo)
        if sets:
            best = size
            break
    if not sets:
        return {"status": "NOT RESOLVED", "reason": "fewer than two valid tests agree",
                "valid_tests": valid}
    if len(sets) > 1:
        return {"status": "NOT RESOLVED", "reason": f"ambiguous: {sets}", "valid_tests": valid}
    return {"status": "ADOPTED", "tests": list(sets[0]), "size": best, "valid_tests": valid}


def mixture_posterior(ages, loglikes: list[np.ndarray]) -> dict:
    """Equal-weight mixture of the agreeing tests' posteriors (on the log-age
    grid, each normalised first) -- the adopted age prior."""
    la = np.log10(ages)
    laf = np.linspace(la.min(), la.max(), 4001)
    dens = []
    for ll in loglikes:
        f = PchipInterpolator(la, ll - np.max(ll))(laf)
        p = np.exp(f)
        dens.append(p / np.trapezoid(p, laf))
    post = np.mean(dens, axis=0)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(laf))])
    cdf /= cdf[-1]

    def q(p):
        return float(10 ** np.interp(p, cdf, laf))
    return {"map": float(10 ** laf[np.argmax(post)]), "median": q(0.5), "lo68": q(0.16),
            "hi68": q(0.84), "lo95": q(0.025), "hi95": q(0.975)}
