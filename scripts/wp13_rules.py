#!/usr/bin/env python3
"""WP13 -- the pre-registered thresholds and scoring rules, as code.

Brief §5 (tasks/wp13_pooled_vs_resolved_ablation_brief.md), owner decisions of
2026-10-08 (provenance/decisions_2026_10_08.json), and the clarifications
A1-A6 of scripts/wp13_prereg.py.  The pre-registration hashes this file, and
the WP13 scorer imports every threshold and rule from here, so the record and
the scoring cannot disagree.

Nothing here reads data.  Importing it computes nothing.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import logsumexp

# ------------------------------------------------------------- design constants
EPS = 0.05                      # test-c outlier fraction (the repair_v9 headline)
POOLED_LABEL = "CygOB2-ALL"     # M0's single label; never written into the WP2 tables
BASELINE = {"family": "PARSEC", "R_V": 3.1, "alpha": 2.3, "sf_duration_Myr": 0.0}
SCORED_ALPHA = 2.3              # D3: T6 is scored on the 18 alpha = 2.3 branches
REPORTED_ALPHAS = (2.0, 2.3)    # D3: the 36 retained branches are reported
LEDGER_ITERATIONS = 2_000_000
RECENT_WINDOW_MYR = 0.1

# ------------------------------------------------------------- thresholds (brief §5, D5)
T1_LNBF_MIN = 2.3
T1_MIN_CELLS = 4                # of the 6 family x R_V cells, plus the baseline cell
T2_ABS = 3.0                    # |dN| >= 3 ...
T2_REL = 0.10                   # ... and |dN| / N_M1 >= 0.10 ...
T2_BAND_REL = 0.10              # ... and the paired 95 % interval excludes +-10 % of N_M1
T3_ABS = 0.10                   # |dP(last < 100 kyr)| >= 0.10, 95 % interval excluding 0
T4_D_MIN = 0.15                 # sup |F_M0 - F_M1| >= 0.15 and above the T7a null ...
T4_NULL_QUANTILE = 95.0         # ... 95th percentile
T4_FIRST_DEATH_SHIFT_MYR = 0.2  # ... or the first-death epoch moves by >= 0.2 Myr
T5_ABS = 2.0                    # |N_M1,sub - N_M0 k_sub / sum k| >= 2 for some subgroup
T6_FRACTION = 0.80              # sign and passing test hold on >= 80 % of the branches
T7_REALISATIONS = 200
T7_ITERATIONS = 20_000          # ledger iterations per model per realisation
T7A_MAX_FALSE_PASS = 0.10
T7B_COVERAGE = (0.58, 0.78)

# ------------------------------------------------------------- D6 grid (descriptive)
D6_AGE_GAP_MYR = (0.0, 0.5, 1.0, 1.5, 2.0)
D6_YOUNG_FRACTION = (0.0, 0.1, 0.2, 0.33, 0.5)
D6_ADEQUATE_N_REL = 0.10        # one-age estimate called adequate where |bias N| <= 10 % ...
D6_ADEQUATE_P_ABS = 0.05        # ... and |bias P(last < 100 kyr)| <= 0.05


# ------------------------------------------------------------- T1
def log_evidence(ages: np.ndarray, loglike: np.ndarray, n_fine: int = 2001) -> float:
    """ln Z under WP4's age prior: uniform in log10(age) over the native grid,
    with the ln L curve PCHIP-interpolated onto the same fine log-age grid that
    wp4_common.posterior_from_loglike uses.  Unlike that function, the curve is
    NOT shifted to its maximum, so Z is comparable between models."""
    la = np.log10(np.asarray(ages, float))
    ll = np.asarray(loglike, float)
    good = np.isfinite(ll) & (ll > -1e8)
    if good.sum() != ll.size:
        raise ValueError("T1 needs a finite ln L at every native age")
    laf = np.linspace(la.min(), la.max(), n_fine)
    llf = PchipInterpolator(la, ll)(laf)
    w = np.full(n_fine, laf[1] - laf[0])
    w[[0, -1]] *= 0.5                                   # trapezoid weights
    return float(logsumexp(llf, b=w) - np.log(la.max() - la.min()))


def t1_ln_bf(ages: np.ndarray, ll_a: np.ndarray, ll_c: np.ndarray) -> float:
    """A1: both models see the same data, the A and C test-c stars.
    M1 gives A and C their own age; M0 gives them one age."""
    return log_evidence(ages, ll_a) + log_evidence(ages, ll_c) - log_evidence(ages, ll_a + ll_c)


def t1_pass(ln_bf_by_cell: dict[tuple[str, float], float]) -> dict:
    if len(ln_bf_by_cell) != 6:
        raise ValueError("T1 needs all 6 family x R_V cells")
    base = ln_bf_by_cell[(BASELINE["family"], BASELINE["R_V"])]
    n_pass = sum(v >= T1_LNBF_MIN for v in ln_bf_by_cell.values())
    return {"baseline_ln_bf": base, "cells_passing": n_pass,
            "pass": bool(base >= T1_LNBF_MIN and n_pass >= T1_MIN_CELLS)}


# ------------------------------------------------------------- T2, T3 (A2)
def interval95(x: np.ndarray) -> tuple[float, float]:
    lo, hi = np.percentile(x, [2.5, 97.5])
    return float(lo), float(hi)


def t2(mu_m0: np.ndarray, mu_m1: np.ndarray) -> dict:
    """A2: on the per-iteration EXPECTED count mu = k * int_{M_to}^{120} M^-alpha
    (parameter uncertainty), paired by iteration index -- not on the Poisson
    realisation, whose 95 % interval could never exclude a +-10 % band at N ~ 7."""
    d = mu_m0 - mu_m1
    n1 = float(mu_m1.mean())
    dn = float(d.mean())
    lo, hi = interval95(d)
    band = T2_BAND_REL * n1
    return {"dN": dn, "dN_lo95": lo, "dN_hi95": hi, "N_M1": n1,
            "pass": bool(abs(dn) >= T2_ABS and abs(dn) / n1 >= T2_REL and (lo > band or hi < -band))}


def t3(recent_m0: np.ndarray, recent_m1: np.ndarray) -> dict:
    """A2: P(last < 100 kyr) is a predictive probability; its paired difference
    has the Monte Carlo interval dP +- 1.96 sd(I0 - I1) / sqrt(n)."""
    d = recent_m0.astype(float) - recent_m1.astype(float)
    dp = float(d.mean())
    se = float(d.std(ddof=1) / np.sqrt(d.size))
    lo, hi = dp - 1.96 * se, dp + 1.96 * se
    return {"dP": dp, "dP_lo95": lo, "dP_hi95": hi,
            "pass": bool(abs(dp) >= T3_ABS and (lo > 0 or hi < 0))}


# ------------------------------------------------------------- T4
def sup_distance(epochs_m0: np.ndarray, epochs_m1: np.ndarray) -> float:
    """D = sup_t |F_M0(t) - F_M1(t)| of the death-epoch distributions (all
    deaths of all iterations, lookback time in Myr)."""
    a, b = np.sort(epochs_m0), np.sort(epochs_m1)
    grid = np.concatenate([a, b])
    return float(np.max(np.abs(np.searchsorted(a, grid, side="right") / a.size
                               - np.searchsorted(b, grid, side="right") / b.size)))


def first_death_epoch(epochs: np.ndarray, iteration: np.ndarray, n_iter: int) -> float:
    """Median over iterations with >= 1 death of the earliest death's lookback time."""
    first = np.full(n_iter, -np.inf)
    np.maximum.at(first, iteration, epochs)
    first = first[np.isfinite(first)]
    return float(np.median(first)) if first.size else float("nan")


def t4(d: float, d_null95: float, first_m0: float, first_m1: float) -> dict:
    shift = abs(first_m0 - first_m1)
    return {"D": d, "D_null95": d_null95, "first_death_shift_Myr": shift,
            "pass": bool((d >= T4_D_MIN and d > d_null95) or shift >= T4_FIRST_DEATH_SHIFT_MYR)}


# ------------------------------------------------------------- T5
def t5(n_m1_sub: dict[str, float], n_m0: float, k_sub: dict[str, float]) -> dict:
    total_k = sum(k_sub.values())
    gap = {sg: n_m1_sub[sg] - n_m0 * k_sub[sg] / total_k for sg in n_m1_sub}
    return {"gap": gap, "pass": bool(max(abs(v) for v in gap.values()) >= T5_ABS)}


# ------------------------------------------------------------- T6
def t6(rows: list[dict], baseline: dict) -> dict:
    """rows: one per scored branch with keys family, dN, dP, t2, t3 (bools).
    The sign of dN and of dP, and whichever of T2/T3 passed on the baseline,
    hold on >= 80 % of the branches AND on both families."""
    def holds(r):
        ok = np.sign(r["dN"]) == np.sign(baseline["dN"]) and np.sign(r["dP"]) == np.sign(baseline["dP"])
        if baseline["t2"]:
            ok &= r["t2"]
        if baseline["t3"]:
            ok &= r["t3"]
        return bool(ok)
    flags = [holds(r) for r in rows]
    by_family = {}
    for fam in sorted({r["family"] for r in rows}):
        f = [holds(r) for r in rows if r["family"] == fam]
        by_family[fam] = float(np.mean(f))
    frac = float(np.mean(flags))
    return {"fraction": frac, "by_family": by_family,
            "pass": bool(frac >= T6_FRACTION and all(v >= T6_FRACTION for v in by_family.values()))}


# ------------------------------------------------------------- verdict (brief §5, D6)
def verdict(t: dict[str, bool]) -> str:
    """t: T1, T2, T3, T4, T6, T7a, T7b -> one of the four §6 outcomes."""
    if not (t["T7a"] and t["T7b"]):
        return "overfit"
    effect = t["T2"] or t["T3"] or t["T4"]
    if t["T1"] and effect and t["T6"]:
        return "material improvement"
    if t["T1"] and not effect:
        return "supported but immaterial"
    if not t["T1"] and not effect:
        return "equivalent"
    return "unclassified"   # e.g. T1 fails but an effect passes: reported as such, never re-mapped
