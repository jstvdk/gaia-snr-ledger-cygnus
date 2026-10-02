#!/usr/bin/env python3
"""Issue #20 Phase A -- shared method: is CygOB2-C really young?

Everything the Phase A tests compute is defined here, and this file is hashed
into provenance/issue20_prereg.json before any test is run.  The analysis
scripts refuse to run if this file no longer matches that hash, unless a
deviation is recorded in provenance/issue20_deviations.json with its reason.

Contents
  * frozen(name)          -- read a preregistered input, verifying its SHA-256
  * spectral-type parsing -- spectral letter, subtype, luminosity class
  * T_eff assignment      -- the anchor table's T_eff, except that a type-derived
                             T_eff of a luminosity class I/II star is replaced by
                             a supergiant scale (Martins+2005 Table 6 for
                             O3-O9.5; Crowther+2006 Table 4 "mean" for O9.7-B3)
  * HRD model             -- an IMF-weighted particle cloud in (log T_eff, M_G0)
                             along the native isochrone, every evolutionary
                             phase included, with an unresolved-binary component
  * likelihoods           -- per star, per age: the conditional p(logTe | M_G0)
                             (primary) and the windowed joint p(logTe, M_G0)
                             (secondary), plus the two-age mixture

Nothing here writes a file.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.special import logsumexp

import wp5_common as w

PREREG_PATH = w.PROVENANCE / "issue20_prereg.json"
DEVIATIONS_PATH = w.PROVENANCE / "issue20_deviations.json"
THIS_FILE = Path(__file__).resolve()

SUBGROUPS = ("CygOB2-A", "CygOB2-B", "CygOB2-C")
FAMILIES = ("PARSEC", "MIST")
R_V_BRANCHES = (3.0, 3.1, 3.5)

# ------------------------------------------------------------------ inputs
# Logical name -> repository-relative path.  Hashed by the preregistration.
INPUTS = {
    "anchor_hrd": "data/processed/wp4_anchor_hrd_repair_v8.parquet",
    "spectroscopic_anchors": "data/processed/wp1_spectroscopic_anchors.parquet",
    "extinction_v5": "data/processed/wp3_extinction_repair_v5.parquet",
    "members": "data/processed/wp2_members.parquet",
    "subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    "isochrones_parsec": "data/processed/wp3_isochrones_parsec.parquet",
    "isochrones_mist": "data/processed/wp3_isochrones_mist.parquet",
    "age_posteriors_v5": "data/processed/wp4_age_posteriors_repair_v5.parquet",
    # A5 only: the attribution of the 3.98 -> 2.51 Myr jump needs the two
    # extinction runs either side of it.  They are read for that comparison
    # and nothing else; wp10_inputs forbids them as manuscript inputs.
    "extinction_pre_repair": "data/processed/wp3_extinction.parquet",
    "extinction_v1": "data/processed/wp3_extinction_repair_v1.parquet",
    "age_posteriors_pre_repair": "data/processed/wp4_age_posteriors.parquet",
    "age_posteriors_v1": "data/processed/wp4_age_posteriors_repair_v1.parquet",
    # A6 only: the stored repair_v8 closure inputs and result
    "wp5_normalization_v8": "data/processed/wp5_imf_normalization_repair_v8.parquet",
    "massive_census_v8": "tables/wp6_massive_census_repair_v8.csv",
    "closure_v8": "tables/wp6_closure_repair_v8.csv",
    "closing_slopes_v8": "tables/wp12_closing_slopes_repair_v8.csv",
    # the scope diagnostic (disclosed prior knowledge; I3 reproduces it)
    "diag_hrd_age_scan": "tables/issue20_hrd_age_scan.csv",
    "diag_c_age_by_version": "tables/issue20_c_age_by_version.csv",
    "diag_bright_members": "tables/issue20_bright_members.csv",
}


class FrozenInputMoved(RuntimeError):
    """An input on disk is not the one the preregistration hashed."""


def prereg() -> dict:
    if not PREREG_PATH.exists():
        raise FileNotFoundError(
            "provenance/issue20_prereg.json is missing -- run "
            "scripts/issue20_prereg.py first.  Phase A tests are meaningful "
            "only against a preregistration that predates them."
        )
    return json.loads(PREREG_PATH.read_text())


def check_method_unchanged() -> dict:
    """This file must match the hash frozen in the preregistration, or the
    change must be recorded as a deviation (with its reason) before use."""
    record = prereg()
    frozen_hash = record["method_code"]["scripts/issue20_common.py"]
    now = w.sha256(THIS_FILE)
    if now == frozen_hash:
        return {"method_hash": now, "deviation": None}
    if DEVIATIONS_PATH.exists():
        deviations = json.loads(DEVIATIONS_PATH.read_text())
        for entry in deviations.get("method_code_changes", []):
            if entry.get("new_sha256") == now:
                return {"method_hash": now, "deviation": entry}
    raise FrozenInputMoved(
        "scripts/issue20_common.py differs from the preregistered method and "
        "no deviation records this version.  Record the change and its reason "
        "in provenance/issue20_deviations.json first."
    )


def frozen(name: str) -> Path:
    """Resolve a preregistered input, verifying its SHA-256."""
    record = prereg()
    entry = record["consumed_inputs"].get(name)
    if entry is None:
        raise KeyError(f"{name!r} is not a preregistered issue-20 input")
    path = w.ROOT / entry["path"]
    digest = w.sha256(path)
    if digest != entry["sha256"]:
        raise FrozenInputMoved(
            f"{entry['path']} changed since the issue-20 preregistration "
            f"({entry['sha256'][:12]} -> {digest[:12]})"
        )
    return path


# ------------------------------------------------------- spectral types
_TYPE = re.compile(r"^\s*([OB])\s*(\d+(?:\.\d+)?)\s*(.*)$")
_LUM = re.compile(r"^(Ia\+|Iab|Ia|Ib|III|II|IV|V|I)")


def parse_spectral_type(text) -> dict:
    """Primary component's spectral letter, numeric subtype and luminosity
    class.  ``O5I+O3.5III`` -> O, 5.0, I, composite.  WR and unparseable types
    return letter None."""
    out = {"letter": None, "subtype": np.nan, "lum_class": "",
           "composite": False, "type_number": np.nan}
    if not isinstance(text, str) or not text.strip():
        return out
    out["composite"] = "+" in text
    primary = text.split("+")[0]
    match = _TYPE.match(primary)
    if not match:
        return out
    letter, subtype, rest = match.group(1), float(match.group(2)), match.group(3)
    lum = _LUM.match(rest.strip())
    lum_class = lum.group(1) if lum else ""
    if lum_class in ("Ia+", "Iab", "Ia", "Ib"):
        lum_class = "I"
    out.update(letter=letter, subtype=subtype, lum_class=lum_class,
               type_number=subtype if letter == "O" else 10.0 + subtype)
    return out


# ------------------------------------------------------------ T_eff scale
# Supergiant (luminosity class I) effective-temperature scale, in K, indexed by
# type number (O s -> s, B s -> 10 + s).
#   O3-O9.5 : Martins, Schaerer & Hillier 2005, A&A 436, 1049, Table 6
#             (observational T_eff scale, class I).
#   O9.7-B3 : Crowther, Lennon & Walborn 2006, A&A 446, 279, Table 4, column
#             "mean" (Galactic Ia/Iab, line-blanketed).
# Both read from the papers' arXiv versions (astro-ph/0503346,
# astro-ph/0509436) on 2026-10-02.  Applied to class II as well (brief A3).
SUPERGIANT_TEFF = {
    3.0: 42233, 4.0: 40422, 5.0: 38612, 5.5: 37706, 6.0: 36801, 6.5: 35895,
    7.0: 34990, 7.5: 34084, 8.0: 33179, 8.5: 32274, 9.0: 31368, 9.5: 30463,
    9.7: 28500, 10.0: 27500, 10.5: 26000, 10.7: 22500, 11.0: 21500,
    11.5: 20500, 12.0: 18500, 12.5: 16500, 13.0: 15500,
}
SIG_LOGTE_DWARF_GIANT = 0.03   # classes III, IV, V and unclassified
SIG_LOGTE_SUPERGIANT = 0.05    # classes I and II


def type_derived(teff: float) -> bool:
    """True when T_eff sits on the 0.01-dex log grid of the spectral-type
    scale that the anchor table used (Wright+15); measured values (Berlanas+20
    quantitative spectroscopy) do not."""
    x = 100.0 * np.log10(teff)
    return bool(abs(x - np.round(x)) < 1e-3)


def supergiant_teff(type_number: float) -> float:
    keys = np.array(sorted(SUPERGIANT_TEFF))
    if not (keys[0] <= type_number <= keys[-1]):
        return np.nan
    values = np.log10([SUPERGIANT_TEFF[k] for k in keys])
    return float(10.0 ** np.interp(type_number, keys, values))


def assign_teff(spectral_type: str, teff_table: float,
                sigma_mode: str = "primary") -> dict:
    """T_eff used by Phase A, its source, and its log-space sigma."""
    parsed = parse_spectral_type(spectral_type)
    supergiant = parsed["lum_class"] in ("I", "II")
    teff, source = float(teff_table), "table_measured"
    if type_derived(teff):
        source = "table_type_scale"
        if supergiant:
            replacement = supergiant_teff(parsed["type_number"])
            if np.isfinite(replacement):
                teff, source = replacement, "supergiant_scale"
            else:
                source = "table_type_scale_no_supergiant_scale"
    if sigma_mode == "primary":
        sigma = SIG_LOGTE_SUPERGIANT if supergiant else SIG_LOGTE_DWARF_GIANT
    elif sigma_mode == "all_0.05":
        sigma = 0.05
    else:
        raise ValueError(sigma_mode)
    return {**parsed, "teff_used": teff, "logTe_used": float(np.log10(teff)),
            "teff_source": source, "sigma_logTe": sigma, "supergiant": supergiant}


# --------------------------------------------------------------- M_G0 error
SIG_DEPTH_MAG = 5.0 / np.log(10.0) * (45.4 / 1624.5)   # WP2 depth 45.4 pc
SIG_CAL_MAG = 0.25    # intrinsic-colour / A_V calibration floor (primary)
SIG_CAL_MAG_WIDE = 0.40   # sensitivity: WP4's lumped SIG_MG0


def sigma_mg0(g_err, av_err, k_g, cal=SIG_CAL_MAG):
    g_err = np.nan_to_num(np.asarray(g_err, float), nan=0.02)
    av_err = np.nan_to_num(np.asarray(av_err, float), nan=0.0)
    return np.sqrt(g_err ** 2 + (k_g * av_err) ** 2 + SIG_DEPTH_MAG ** 2 + cal ** 2)


# ---------------------------------------------------------------- HRD model
IMF_SLOPE = 2.3        # along-isochrone weight, as wp4_common.IMF_SLOPE
M_MIN = 2.0            # lowest primary mass carried (anchors are M_G0 < 0)
M_BIRTH_MAX = 300.0    # birth normalisation upper limit (weights = fraction born)
Q_MIN = 0.10           # as WP4
N_Q = 10
DLOGTE_STEP = 0.005    # sub-sampling resolution along the isochrone
DMAG_STEP = 0.025
MAX_SUBSTEPS = 400
MIST_GAP_DLOGTE = 0.05     # MIST: a phase change with a jump this large in
MIST_GAP_DMAG = 0.5        # logTe or G between tabulated points is a gap


def _imf_integral(lo, hi, alpha=IMF_SLOPE):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    return (hi ** (1.0 - alpha) - lo ** (1.0 - alpha)) / (1.0 - alpha)


def _flux_add(m1, m2):
    return -2.5 * np.log10(10.0 ** (-0.4 * m1) + 10.0 ** (-0.4 * m2))


def isochrone_at(iso_family: pd.DataFrame, age: float) -> pd.DataFrame:
    ages = np.unique(iso_family["age_Myr"].to_numpy(float))
    native = float(ages[np.argmin(np.abs(ages - age))])
    frame = iso_family[np.isclose(iso_family["age_Myr"], native)]
    frame = frame.dropna(subset=["Mini", "logTe", "G0"]).sort_values("Mini")
    return frame.drop_duplicates("Mini", keep="first")


def hrd_particles(iso_age: pd.DataFrame, family: str, f_bin: float) -> dict:
    """IMF-weighted particle cloud in (logTe, M_G0) for one native age.

    Each pair of consecutive tabulated points (sorted by initial mass) is a
    segment of the isochrone; it is sub-sampled finely enough to resolve
    0.005 dex in logTe and 0.025 mag in G, and each sub-point carries the IMF
    integral over its initial-mass interval.  Every evolutionary phase is kept:
    post-main-sequence supergiants are the point of the test.  MIST tables omit
    some phases (e.g. between the end of the MS and the WR phase); a segment
    whose two ends differ in phase and jump by > 0.05 dex in logTe or > 0.5 mag
    in G is such a gap and carries no weight -- interpolating across it would
    invent stars.
    PARSEC tables are contiguous and are never cut.

    Unresolved binaries: a fraction f_bin of primaries has a secondary of mass
    q * m, q uniform on [0.1, 1] (10 values), whose G is read off the same
    isochrone at its initial mass and flux-added.  logTe is the primary's (the
    spectral type is the primary's).

    Weights are fractions of the stars born between 2 and 300 Msun, so two ages
    can be mixed with a birth fraction.
    """
    mini = iso_age["Mini"].to_numpy(float)
    logte = iso_age["logTe"].to_numpy(float)
    g = iso_age["G0"].to_numpy(float)
    phase = (iso_age["phase"].to_numpy(float) if family == "MIST"
             else np.zeros(len(mini)))
    seg_m, seg_t, seg_g, seg_w = [], [], [], []
    dropped = 0.0
    for i in range(len(mini) - 1):
        m_a, m_b = mini[i], mini[i + 1]
        if m_b <= M_MIN:
            continue
        gap = (family == "MIST" and phase[i] != phase[i + 1]
               and (abs(logte[i + 1] - logte[i]) > MIST_GAP_DLOGTE
                    or abs(g[i + 1] - g[i]) > MIST_GAP_DMAG))
        lo = max(m_a, M_MIN)
        if gap:
            dropped += float(_imf_integral(lo, m_b))
            continue
        frac_lo = (lo - m_a) / (m_b - m_a)
        t_a = logte[i] + frac_lo * (logte[i + 1] - logte[i])
        g_a = g[i] + frac_lo * (g[i + 1] - g[i])
        n = int(np.clip(np.ceil(max(abs(logte[i + 1] - t_a) / DLOGTE_STEP,
                                    abs(g[i + 1] - g_a) / DMAG_STEP, 1.0)),
                        1, MAX_SUBSTEPS))
        edges = lo + (m_b - lo) * np.arange(n + 1) / n
        mid = 0.5 * (edges[1:] + edges[:-1])
        frac = (mid - m_a) / (m_b - m_a)
        seg_m.append(mid)
        seg_t.append(logte[i] + frac * (logte[i + 1] - logte[i]))
        seg_g.append(g[i] + frac * (g[i + 1] - g[i]))
        seg_w.append(_imf_integral(edges[:-1], edges[1:]))
    m = np.concatenate(seg_m)
    t = np.concatenate(seg_t)
    gg = np.concatenate(seg_g)
    wt = np.concatenate(seg_w)
    birth = float(_imf_integral(M_MIN, M_BIRTH_MAX))
    represented = float(wt.sum())

    q = np.linspace(Q_MIN, 1.0, N_Q)
    g2 = np.interp(np.outer(m, q), mini, g, left=99.0)
    g_bin = _flux_add(gg[:, None], g2)
    logte_all = np.concatenate([t, np.repeat(t, N_Q)])
    g_all = np.concatenate([gg, g_bin.ravel()])
    w_all = np.concatenate([(1.0 - f_bin) * wt,
                            np.repeat(f_bin * wt / N_Q, N_Q)]) / birth
    return {"logTe": logte_all, "G": g_all, "w": w_all, "log_w": np.log(w_all),
            "alive_fraction": represented / birth,
            "dropped_gap_fraction": dropped / max(represented + dropped, 1e-300),
            "imf_check": abs(represented + dropped
                             - float(_imf_integral(M_MIN, max(mini.max(), M_MIN))))
                         / float(_imf_integral(M_MIN, max(mini.max(), M_MIN)))}


def star_terms(logte, mg, sig_t, sig_m, particles) -> tuple[np.ndarray, np.ndarray]:
    """Per star: ln p(logTe, M_G0 | age) (birth-normalised) and ln p(M_G0 | age).
    Stars are processed in chunks of 8 to bound memory; the result does not
    depend on the chunking."""
    logte, mg = np.asarray(logte, float), np.asarray(mg, float)
    sig_t, sig_m = np.asarray(sig_t, float), np.asarray(sig_m, float)
    lt = particles["logTe"][None, :]
    gm = particles["G"][None, :]
    lw = particles["log_w"][None, :]
    joint = np.empty(len(logte))
    marginal = np.empty(len(logte))
    for start in range(0, len(logte), 8):
        s = slice(start, start + 8)
        st, sm = sig_t[s, None], sig_m[s, None]
        log_nm = (-0.5 * ((mg[s, None] - gm) / sm) ** 2
                  - np.log(np.sqrt(2 * np.pi) * sm))
        log_nt = (-0.5 * ((logte[s, None] - lt) / st) ** 2
                  - np.log(np.sqrt(2 * np.pi) * st))
        joint[s] = logsumexp(lw + log_nm + log_nt, axis=1)
        marginal[s] = logsumexp(lw + log_nm, axis=1)
    return joint, marginal


def window_log_norm(particles, edge: float) -> float:
    inside = particles["G"] <= edge
    return float(np.log(particles["w"][inside].sum())) if inside.any() else -np.inf


# ---------------------------------------------------------------- posterior
def posterior_summary(ages: np.ndarray, loglike: np.ndarray, n_fine: int = 4001) -> dict:
    """wp4_common.posterior_from_loglike's construction (PCHIP in log age,
    uniform prior in log age over the native grid), with 95 % bounds added."""
    la = np.log10(np.asarray(ages, float))
    ll = np.asarray(loglike, float)
    good = np.isfinite(ll)
    la, ll = la[good], ll[good]
    laf = np.linspace(la.min(), la.max(), n_fine)
    llf = PchipInterpolator(la, ll)(laf)
    llf -= llf.max()
    post = np.exp(llf)
    post /= np.trapezoid(post, laf)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(laf))])
    cdf /= cdf[-1]

    def q(p):
        return float(10.0 ** np.interp(p, cdf, laf))

    age_map = float(10.0 ** laf[int(np.argmax(post))])
    step = float(np.median(np.diff(la)))
    railed = bool(np.log10(age_map) <= la[0] + step + 1e-3
                  or np.log10(age_map) >= la[-1] - step - 1e-3)
    return {"map": age_map, "median": q(0.5), "lo68": q(0.16), "hi68": q(0.84),
            "lo95": q(0.025), "hi95": q(0.975),
            "mean": float(10.0 ** np.trapezoid(laf * post, laf)),
            "railed": railed, "lnL_max": float(np.max(ll))}


# ----------------------------------------------------------------- mixture
PHI_GRID = np.round(np.arange(0.02, 0.981, 0.02), 2)
MIN_COMPONENT_STARS = 3.0


def single_best(cond: np.ndarray, weights: np.ndarray) -> tuple[int, float]:
    """cond: (n_star, n_age) conditional ln L.  Best single age index, lnL."""
    total = weights @ cond
    j = int(np.argmax(total))
    return j, float(total[j])


def mixture_best(joint: np.ndarray, marginal: np.ndarray, weights: np.ndarray,
                 ages: np.ndarray, statistic: str = "conditional",
                 log_norm: np.ndarray | None = None) -> dict:
    """Two-age mixture: phi born at ages[j1], 1-phi at ages[j2] (j1 < j2).

    conditional:  ln[(phi e^J1 + (1-phi) e^J2) / (phi e^M1 + (1-phi) e^M2)]
    joint      :  ln[(phi e^J1 + (1-phi) e^J2) / (phi N1 + (1-phi) N2)]
                  with N the windowed normalisation of each age.
    A component must hold >= 3 membership-weighted stars (responsibilities).
    """
    n_age = len(ages)
    best = {"lnL": -np.inf}
    for j1 in range(n_age):
        for j2 in range(j1 + 1, n_age):
            for phi in PHI_GRID:
                a, b = np.log(phi), np.log(1.0 - phi)
                num = np.logaddexp(a + joint[:, j1], b + joint[:, j2])
                if statistic == "conditional":
                    den = np.logaddexp(a + marginal[:, j1], b + marginal[:, j2])
                else:
                    den = np.full(len(num), np.logaddexp(a + log_norm[j1], b + log_norm[j2]))
                ll = float(weights @ (num - den))
                if ll <= best["lnL"]:
                    continue
                resp = np.exp(a + joint[:, j1] - num)
                n1 = float(weights @ resp)
                n2 = float(weights.sum() - n1)
                if min(n1, n2) < MIN_COMPONENT_STARS:
                    continue
                best = {"lnL": ll, "age_young": float(ages[j1]),
                        "age_old": float(ages[j2]), "phi_young": float(phi),
                        "n_young": n1, "n_old": n2}
    return best


def delta_bic(lnl_mix: float, lnl_single: float, n_eff: float) -> float:
    """BIC_single - BIC_mixture, k = 1 vs 3 parameters."""
    return 2.0 * (lnl_mix - lnl_single) - 2.0 * np.log(n_eff)
