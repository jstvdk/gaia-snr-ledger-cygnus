#!/usr/bin/env python3
"""WP11 Part B — the 26Al / 60Fe forward prediction.

Pre-registered in provenance/wp11_isotope_prereg.json.  The three yield arms,
the interpolation rule, the estimator, the instrument values, predictions I1-I5
and the detectability rule are all fixed there and are APPLIED here, not
amended.

The population engine is IMPORTED from wp7_ledger, exactly as WP9 imports it, so
the forecast cannot drift from the ledger it is derived from.  Nothing upstream
is read except through the frozen artifacts.

What it computes, per branch and per yield arm:

  M_iso(now) = < sum over the branch's supernovae of
                 y_iso(m_progenitor) * exp(-t_explosion / tau_iso) >

evaluated at each supernova's OWN progenitor mass, then converted to a line
flux at the WP3/WP2 distance artifact.  The steady-state approximation the brief
sketches is reported alongside as a saturation ratio, never as the forecast --
Cyg OB2 is far too young for 60Fe to have saturated.

Outputs:
  tables/wp11_isotope_forecast.csv        per branch x arm
  tables/wp11_isotope_summary.csv         per arm, headline set rolled up
  provenance/wp11_isotope_forecast_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp11_isotope_forecast.py
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp3_common import D_KPC
from wp7_ledger import TurnoffRelation, draw_key, run_population
from wp7_ledger_prereg import SF_DURATIONS_MYR, SN_THRESHOLD_MSUN
from wp11_isotope_prereg import (
    HEADLINE_ALPHAS,
    INSTRUMENTS,
    ISOTOPES,
    LC06_AL26_TOTAL,
    LC06_AL26_WIND,
    LC06_FE60_TOTAL,
    LC06_MASS_MSUN,
    LC06L_AL26_SN,
    LC06L_FE60_SN,
    LC06L_MASS_MSUN,
    PRIMARY_EXPLODABILITY,
)

WP5_VERSION = "repair_v7"
ITERATIONS = 200_000
PRIMARY_ARM = "LC06_NL"

# LC18's Recommended scenario collapses everything above this mass.
LC18_COLLAPSE_ABOVE_MSUN = 25.0

# WP8 §5's comparison figure for the complex-wide 26Al inventory.  It is the
# denominator prediction I1 was pre-registered against, so it is reused VERBATIM
# and I1 is scored against it -- but see finding T4 below: this figure does not
# survive contact with the frozen WP1 flux, and the corrected, conversion-free
# comparison is computed alongside it.
COMPLEX_AL26_MSUN_WP8 = 1.0

# Physical constants.
MSUN_KG = 1.98847e30
ATOMIC_MASS_UNIT_KG = 1.66053906660e-27
KPC_CM = 3.0856775814913673e21
SECONDS_PER_MYR = 3.155815e13

# The most recent lookback bin used for the present-day rate, matching the
# 0.05 Myr binning of tables/wp7_rsn_curves.csv.
RATE_BIN_MYR = 0.05


def log_interp(mass: np.ndarray, grid: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Linear in log10(yield) against log10(mass), endpoint-clipped.

    The interpolation rule fixed in the pre-registration.  np.interp already
    holds the endpoint value outside the grid, which is the declared
    no-extrapolation behaviour.
    """
    return 10.0 ** np.interp(
        np.log10(mass), np.log10(grid), np.log10(values)
    )


def yields_for(arm: str, mass: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """SN-only (26Al, 60Fe) yields in Msun for each progenitor mass."""
    if arm == "LC06_NL":
        al = log_interp(mass, LC06_MASS_MSUN, LC06_AL26_TOTAL - LC06_AL26_WIND)
        fe = log_interp(mass, LC06_MASS_MSUN, LC06_FE60_TOTAL)
        return al, fe
    if arm == "LC06_Langer":
        al = log_interp(mass, LC06L_MASS_MSUN, LC06L_AL26_SN)
        fe = log_interp(mass, LC06L_MASS_MSUN, LC06L_FE60_SN)
        return al, fe
    if arm == "LC18_REC":
        # Zero above the collapse threshold; below it LC18 does explode, so the
        # LC06 curve stands in there.  The execution record reports how many
        # events fall below 25 Msun -- if none do, as the ledger implies, this
        # arm is identically zero and the fallback never fires.
        al, fe = yields_for("LC06_NL", mass)
        collapsed = mass > LC18_COLLAPSE_ABOVE_MSUN
        return np.where(collapsed, 0.0, al), np.where(collapsed, 0.0, fe)
    raise KeyError(f"unknown yield arm: {arm}")


def line_flux(mass_msun: float, isotope: str, line_keV: float) -> float:
    """Photon flux in one decay line, ph cm^-2 s^-1, at the WP3 distance."""
    spec = ISOTOPES[isotope]
    n_nuclei = mass_msun * MSUN_KG / (spec["mass_number"] * ATOMIC_MASS_UNIT_KG)
    tau_s = spec["mean_lifetime_Myr"] * SECONDS_PER_MYR
    p_gamma = spec["photons_per_decay"][f"{line_keV}"]
    area_cm2 = 4.0 * np.pi * (D_KPC * KPC_CM) ** 2
    return n_nuclei / tau_s * p_gamma / area_cm2


def mass_from_flux(flux_ph_cm2_s: float, isotope: str, line_keV: float) -> float:
    """Invert a measured LINE FLUX into an isotope MASS at the WP3 distance.

    Permitted and not the same thing as the prohibited operation: the WP8 ban is
    on inverting the 26Al flux into a SUPERNOVA COUNT, which would require a
    yield and an explodability assumption.  A flux-to-mass conversion needs only
    the decay constant and the distance.
    """
    spec = ISOTOPES[isotope]
    tau_s = spec["mean_lifetime_Myr"] * SECONDS_PER_MYR
    p_gamma = spec["photons_per_decay"][f"{line_keV}"]
    area_cm2 = 4.0 * np.pi * (D_KPC * KPC_CM) ** 2
    n_nuclei = flux_ph_cm2_s * area_cm2 / p_gamma * tau_s
    return n_nuclei * spec["mass_number"] * ATOMIC_MASS_UNIT_KG / MSUN_KG


def verdict(fluxes: np.ndarray, threshold: float) -> str:
    """The pre-registered detectability rule, applied mechanically."""
    if bool((fluxes >= threshold).all()):
        return "DETECTABLE"
    if bool((fluxes < threshold).all()):
        return "BELOW_REACH"
    return "MARGINAL"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=ITERATIONS)
    args = parser.parse_args()
    n_iter = int(args.iterations)

    draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{WP5_VERSION}.npz")
    relations = {family: TurnoffRelation(family) for family in w.FAMILIES}
    master = np.random.default_rng(w.SEED)

    tau26 = ISOTOPES["26Al"]["mean_lifetime_Myr"]
    tau60 = ISOTOPES["60Fe"]["mean_lifetime_Myr"]
    arms = ("LC06_NL", "LC06_Langer", "LC18_REC")

    rows = []
    events_below_collapse = 0
    total_events = 0
    convergence = []

    for family in w.FAMILIES:
        relation = relations[family]
        for rv in w.R_V_BRANCHES:
            for alpha in w.IMF_SLOPES:
                for delta in SF_DURATIONS_MYR:
                    # Per-iteration accumulators, summed over the three
                    # subgroups: this is one association, not three.
                    acc = {
                        arm: {
                            "al26": np.zeros(n_iter),
                            "fe60": np.zeros(n_iter),
                        }
                        for arm in arms
                    }
                    n_sn = np.zeros(n_iter)
                    recent = 0
                    mass_sum = 0.0
                    yield_sum = {arm: [0.0, 0.0] for arm in arms}
                    n_events = 0

                    for subgroup in w.SUBGROUPS:
                        key = draw_key(subgroup, family, rv, alpha)
                        k_all = draws[f"k__{key}"]
                        age_all = draws[f"truth_age_draws__{key}"]
                        rng = np.random.default_rng(
                            master.integers(0, 2 ** 63 - 1)
                        )
                        pick = rng.integers(0, k_all.size, n_iter)
                        res = run_population(
                            rng, k_all[pick], age_all[pick], alpha, delta,
                            relation,
                        )
                        mass = res["dead_masses"]
                        epoch = res["epochs"]
                        it = res["dead_iteration"]
                        explodes = mass >= SN_THRESHOLD_MSUN
                        if not explodes.any():
                            continue
                        m_sn = mass[explodes]
                        t_sn = epoch[explodes]
                        i_sn = it[explodes]

                        n_sn += np.bincount(i_sn, minlength=n_iter)
                        recent += int((t_sn < RATE_BIN_MYR).sum())
                        mass_sum += float(m_sn.sum())
                        n_events += int(m_sn.size)
                        events_below_collapse += int(
                            (m_sn <= LC18_COLLAPSE_ABOVE_MSUN).sum()
                        )

                        decay26 = np.exp(-t_sn / tau26)
                        decay60 = np.exp(-t_sn / tau60)
                        for arm in arms:
                            y26, y60 = yields_for(arm, m_sn)
                            acc[arm]["al26"] += np.bincount(
                                i_sn, weights=y26 * decay26, minlength=n_iter
                            )
                            acc[arm]["fe60"] += np.bincount(
                                i_sn, weights=y60 * decay60, minlength=n_iter
                            )
                            yield_sum[arm][0] += float(y26.sum())
                            yield_sum[arm][1] += float(y60.sum())

                    total_events += n_events
                    rate_now = recent / n_iter / RATE_BIN_MYR
                    mean_mass = mass_sum / n_events if n_events else float("nan")

                    for arm in arms:
                        m26 = float(acc[arm]["al26"].mean())
                        m60 = float(acc[arm]["fe60"].mean())
                        # Steady-state diagnostic: rate x lifetime x the
                        # population-mean yield on this arm.
                        ybar26 = yield_sum[arm][0] / n_events if n_events else 0.0
                        ybar60 = yield_sum[arm][1] / n_events if n_events else 0.0
                        ss26 = rate_now * tau26 * ybar26
                        ss60 = rate_now * tau60 * ybar60

                        f1809 = line_flux(m26, "26Al", 1808.65)
                        f1173 = line_flux(m60, "60Fe", 1173.2)
                        f1332 = line_flux(m60, "60Fe", 1332.5)

                        rows.append({
                            "family": family, "R_V": rv, "alpha": alpha,
                            "sf_duration_Myr": delta,
                            "explodability": PRIMARY_EXPLODABILITY,
                            "yield_arm": arm,
                            "in_headline_set": alpha in HEADLINE_ALPHAS,
                            "N_SN_mean": round(float(n_sn.mean()), 4),
                            "mean_progenitor_Msun": round(mean_mass, 2),
                            "rate_now_per_Myr": round(rate_now, 4),
                            "M_al26_Msun": m26,
                            "M_fe60_Msun": m60,
                            "M_al26_steady_state_Msun": ss26,
                            "M_fe60_steady_state_Msun": ss60,
                            "saturation_al26": (
                                m26 / ss26 if ss26 > 0 else float("nan")
                            ),
                            "saturation_fe60": (
                                m60 / ss60 if ss60 > 0 else float("nan")
                            ),
                            "F_1809_ph_cm2_s": f1809,
                            "F_1173_ph_cm2_s": f1173,
                            "F_1332_ph_cm2_s": f1332,
                            "F_fe60_combined_ph_cm2_s": f1173 + f1332,
                        })

                    # Split-half convergence on the primary arm's 60Fe mass.
                    half = acc[PRIMARY_ARM]["fe60"]
                    a, b = float(half[: n_iter // 2].mean()), float(
                        half[n_iter // 2:].mean()
                    )
                    convergence.append(
                        abs(a - b) / ((a + b) / 2.0) if (a + b) > 0 else 0.0
                    )
                print(f"  {family} R_V={rv} alpha={alpha} done", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(w.TABLES / "wp11_isotope_forecast.csv", index=False)

    head = table[table.in_headline_set]
    primary = head[head.yield_arm.eq(PRIMARY_ARM)]

    # ------------------------------------------------------------ per-arm roll-up
    spi_fe60 = INSTRUMENTS["SPI_Cygnus_60Fe_upper_limit"]["value_ph_cm2_s"]
    spi_al26 = INSTRUMENTS["SPI_Cygnus_26Al_complex_flux"]["value_ph_cm2_s"]
    cosi = INSTRUMENTS["COSI_narrow_line_3sigma_2yr"]["value_ph_cm2_s"]

    summary = []
    for arm in arms:
        block = head[head.yield_arm.eq(arm)]
        summary.append({
            "yield_arm": arm,
            "branches": len(block),
            "M_al26_min_Msun": block.M_al26_Msun.min(),
            "M_al26_median_Msun": block.M_al26_Msun.median(),
            "M_al26_max_Msun": block.M_al26_Msun.max(),
            "M_fe60_min_Msun": block.M_fe60_Msun.min(),
            "M_fe60_median_Msun": block.M_fe60_Msun.median(),
            "M_fe60_max_Msun": block.M_fe60_Msun.max(),
            "F_1809_min": block.F_1809_ph_cm2_s.min(),
            "F_1809_median": block.F_1809_ph_cm2_s.median(),
            "F_1809_max": block.F_1809_ph_cm2_s.max(),
            "F_fe60_line_min": block.F_1173_ph_cm2_s.min(),
            "F_fe60_line_median": block.F_1173_ph_cm2_s.median(),
            "F_fe60_line_max": block.F_1173_ph_cm2_s.max(),
            "F_fe60_combined_median": block.F_fe60_combined_ph_cm2_s.median(),
            "al26_fraction_of_complex_flux":
                block.F_1809_ph_cm2_s.median() / spi_al26,
            "fe60_over_spi_limit":
                block.F_fe60_combined_ph_cm2_s.max() / spi_fe60,
            "cosi_verdict_fe60": verdict(
                block.F_1173_ph_cm2_s.to_numpy(), cosi
            ),
            "cosi_verdict_al26_sn_only": verdict(
                block.F_1809_ph_cm2_s.to_numpy(), cosi
            ),
            "spi_verdict_fe60": (
                "BELOW_LIMIT"
                if bool((block.F_fe60_combined_ph_cm2_s < spi_fe60).all())
                else "EXCEEDS_LIMIT_ON_SOME_BRANCH"
            ),
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(w.TABLES / "wp11_isotope_summary.csv", index=False)

    # ------------------------------------------- T4: the WP8 denominator is wrong
    # The frozen WP1 marker gives a MEASURED Cygnus-complex 1809 keV flux.  At
    # the WP3 distance that flux corresponds to about 1e-2 Msun of 26Al, not the
    # ~1 Msun WP8 §5 compared against -- WP8 appears to have taken Martin et
    # al.'s GALACTIC stationary mass (1.7-2.0 Msun) for the complex.  Nothing
    # upstream is retuned: the defect is recorded, and the comparison is
    # restated in FLUX space, where no mass conversion of the measurement is
    # needed at all and the error cannot propagate.
    complex_al26_msun_measured = mass_from_flux(spi_al26, "26Al", 1808.65)
    al26_flux_fraction = primary.F_1809_ph_cm2_s / spi_al26
    t4 = {
        "id": "T4",
        "severity": "corrects a published WP8 number; changes no ledger value",
        "statement": (
            "WP8 §5 compares the ledger's SN-only 26Al mass against '~1 Msun "
            "inferred for the whole Cygnus complex' and reports the ordering as "
            "roughly 1000x.  The frozen WP1 marker's own measured complex flux "
            f"of {spi_al26:.1e} ph/cm2/s corresponds to only "
            f"{complex_al26_msun_measured:.2e} Msun of 26Al at "
            f"{D_KPC} kpc.  The 1 Msun figure is about two orders of magnitude "
            "too large and looks like Martin et al.'s GALACTIC stationary mass "
            "(1.7-2.0 Msun) used for the complex."
        ),
        "consequence": (
            "the true ordering is roughly one order of magnitude, not three: "
            "this ledger's supernovae alone supply "
            f"{100 * al26_flux_fraction.min():.0f}-"
            f"{100 * al26_flux_fraction.max():.0f}% of the measured "
            "complex-wide 1809 keV flux.  Supernovae are a sub-dominant but "
            "NOT negligible contributor, which is a stronger and more "
            "interesting statement than WP8 made -- and it is still a lower "
            "bound on the complex-wide signal, so it is consistent with the "
            "measurement and is not a tension with data."
        ),
        "measured_complex_al26_Msun": complex_al26_msun_measured,
        "wp8_assumed_complex_al26_Msun": COMPLEX_AL26_MSUN_WP8,
        "sn_only_flux_fraction_of_complex": {
            "min": float(al26_flux_fraction.min()),
            "median": float(al26_flux_fraction.median()),
            "max": float(al26_flux_fraction.max()),
        },
        "action_taken": (
            "NONE upstream.  WP11 is one-way validation: the defect is "
            "reported, WP8 is not edited, and the manuscript quotes the "
            "flux-space comparison, which needs no mass conversion of the "
            "measurement."
        ),
    }

    # ---------------------------------------------------- score the predictions
    # I1 -- SN-only 26Al at least 100x below the complex-wide inventory.
    # Scored against the denominator NAMED IN THE PRE-REGISTRATION, unamended.
    # T4 shows that denominator is wrong, which makes the PASS vacuous; both
    # facts are recorded rather than either being suppressed.
    al26_ratio_max = float(primary.M_al26_Msun.max() / COMPLEX_AL26_MSUN_WP8)
    i1 = bool(al26_ratio_max <= 0.01)
    i1_corrected = bool(float(al26_flux_fraction.max()) <= 0.01)

    # I2 -- 60Fe below the frozen SPI Cygnus upper limit on every branch.
    fe60_max_combined = float(primary.F_fe60_combined_ph_cm2_s.max())
    i2 = bool(fe60_max_combined < spi_fe60)

    # I3 -- alpha = 2.0 predicts at least 2x the alpha = 2.3 60Fe flux.
    med_20 = float(primary[primary.alpha.eq(2.0)].F_1173_ph_cm2_s.median())
    med_23 = float(primary[primary.alpha.eq(2.3)].F_1173_ph_cm2_s.median())
    alpha_ratio = med_20 / med_23 if med_23 > 0 else float("inf")
    i3 = bool(alpha_ratio >= 2.0)

    # I4 -- the between-arm spread exceeds the within-primary-arm branch spread.
    within = (
        float(primary.F_1173_ph_cm2_s.max() / primary.F_1173_ph_cm2_s.min())
        if primary.F_1173_ph_cm2_s.min() > 0 else float("inf")
    )
    non_null = [a for a in arms if a != "LC18_REC"]
    arm_medians = [
        float(head[head.yield_arm.eq(a)].F_1173_ph_cm2_s.median())
        for a in non_null
    ]
    between = max(arm_medians) / min(arm_medians) if min(arm_medians) > 0 else float("inf")
    i4 = bool(between > within)

    # I5 -- COSI verdict for 60Fe on the primary arm is MARGINAL.
    cosi_primary = verdict(primary.F_1173_ph_cm2_s.to_numpy(), cosi)
    i5 = bool(cosi_primary == "MARGINAL")

    predictions = [
        {"id": "I1", "outcome": "PASS" if i1 else "FAIL",
         "statement": "SN-only 26Al mass <= 1% of the ~1 Msun complex inventory",
         "measured": {"max_fraction_of_complex_Msun": al26_ratio_max,
                      "threshold": 0.01},
         "integrity_note": (
             "VACUOUS AS WRITTEN.  Scored against the denominator named in the "
             "pre-registration and not amended, but finding T4 shows that "
             "denominator is about 100x too large, so the PASS carries no "
             "information.  Re-scored against the measured complex-wide FLUX -- "
             "the same test with the error removed -- the outcome is "
             + ("PASS" if i1_corrected else "FAIL")
             + ".  The corrected comparison is what the manuscript quotes."
         ),
         "outcome_against_corrected_denominator":
             "PASS" if i1_corrected else "FAIL"},
        {"id": "I2", "outcome": "PASS" if i2 else "FAIL",
         "statement": "60Fe below the frozen SPI Cygnus upper limit everywhere",
         "measured": {"max_combined_flux_ph_cm2_s": fe60_max_combined,
                      "spi_upper_limit_ph_cm2_s": spi_fe60,
                      "margin_factor": spi_fe60 / fe60_max_combined
                      if fe60_max_combined > 0 else float("inf")}},
        {"id": "I3", "outcome": "PASS" if i3 else "FAIL",
         "statement": "alpha=2.0 median 60Fe flux >= 2x the alpha=2.3 median",
         "measured": {"median_alpha_2p0": med_20, "median_alpha_2p3": med_23,
                      "ratio": alpha_ratio, "threshold": 2.0}},
        {"id": "I4", "outcome": "PASS" if i4 else "FAIL",
         "statement": "between-arm yield spread > within-arm branch spread",
         "measured": {"between_arm_factor": between,
                      "within_arm_factor": within}},
        {"id": "I5", "outcome": "PASS" if i5 else "FAIL",
         "statement": "COSI verdict for 60Fe on the primary arm is MARGINAL",
         "measured": {"verdict": cosi_primary,
                      "cosi_sensitivity_ph_cm2_s": cosi}},
    ]

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp11_isotope_forecast.py",
        "status": "SUCCESS",
        "work_package": "WP11 Part B",
        "prereg": "provenance/wp11_isotope_prereg.json",
        "brief": "tasks/wp11_bowshock_isotope_brief.md §4",
        "post_hoc_disclosure": (
            "this comparison was NOT frozen at WP1.  It is a post-hoc addition, "
            "pre-registered before scoring but chosen after the ledger existed. "
            "It must be introduced as such in the manuscript."
        ),
        "iterations": n_iter,
        "wp5_version": WP5_VERSION,
        "distance_kpc": D_KPC,
        "distance_source": (
            "wp3_common.D_KPC <- provenance/wp2_membership_manifest.json "
            "one_component_mean_kpc"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "headline_branches": int(len(primary)),
        "yield_arms": list(arms),
        "primary_arm": PRIMARY_ARM,
        "lc18_null_arm_check": {
            "supernovae_sampled": total_events,
            "supernovae_at_or_below_25_Msun": events_below_collapse,
            "meaning": (
                "LC18's Recommended scenario collapses everything above "
                "25 Msun.  A count of zero here means the null arm is "
                "identically zero for this ledger, as the pre-registration "
                "asserted, rather than merely nearly so."
            ),
        },
        "convergence": {
            "statistic": "split-half relative difference in the primary arm's "
                         "60Fe mass, per branch",
            "max": float(np.max(convergence)),
            "median": float(np.median(convergence)),
        },
        "summary_by_arm": summary_df.to_dict(orient="records"),
        "predicted_isotope_ratio": {
            "fe60_combined_over_al26_sn_only": {
                "min": float((primary.F_fe60_combined_ph_cm2_s
                              / primary.F_1809_ph_cm2_s).min()),
                "median": float((primary.F_fe60_combined_ph_cm2_s
                                 / primary.F_1809_ph_cm2_s).median()),
                "max": float((primary.F_fe60_combined_ph_cm2_s
                              / primary.F_1809_ph_cm2_s).max()),
            },
            "galactic_measured_for_context":
                INSTRUMENTS["galactic_ratio_context"]["fe60_over_al26"],
            "why_ours_is_higher": (
                "our 26Al denominator is SUPERNOVA-ONLY while the Galactic "
                "ratio's is winds plus supernovae, and every progenitor here is "
                "above 30 Msun where LC06's 60Fe yield is largest.  The two "
                "numbers are not measuring the same quantity and the comparison "
                "is context, not a check."
            ),
        },
        "findings": [t4],
        "predictions": predictions,
        "prohibition_honoured": (
            "the measured complex-wide 26Al flux was never inverted into a "
            "supernova count.  The SN-only component is reported as a LOWER "
            "BOUND on the complex-wide signal."
        ),
        "outputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.TABLES / "wp11_isotope_forecast.csv",
                w.TABLES / "wp11_isotope_summary.csv",
            ]
        },
    }
    w.write_json(
        w.PROVENANCE / "wp11_isotope_forecast_execution.json", record
    )

    # ------------------------------------------------------------------ report
    print("\nWP11 Part B — 26Al / 60Fe forward prediction\n")
    print(f"  {len(primary)} headline branches x {len(arms)} yield arms, "
          f"{n_iter:,} iterations, d = {D_KPC} kpc\n")
    print(f"  {'arm':12s} {'M(26Al)/Msun':>22s} {'M(60Fe)/Msun':>22s}")
    for entry in summary:
        print(f"  {entry['yield_arm']:12s} "
              f"{entry['M_al26_min_Msun']:.2e}-{entry['M_al26_max_Msun']:.2e}   "
              f"{entry['M_fe60_min_Msun']:.2e}-{entry['M_fe60_max_Msun']:.2e}")
    print(f"\n  {'arm':12s} {'F(1809)':>10s} {'F(1173)':>10s}  "
          f"{'COSI 60Fe':>12s}  {'SPI 60Fe':>12s}")
    for entry in summary:
        print(f"  {entry['yield_arm']:12s} "
              f"{entry['F_1809_median']:.2e} {entry['F_fe60_line_median']:.2e}  "
              f"{entry['cosi_verdict_fe60']:>12s}  "
              f"{entry['spi_verdict_fe60']:>12s}")
    print(f"\n  COSI 3-sigma narrow line (2 yr): {cosi:.1e} ph/cm2/s")
    print(f"  SPI Cygnus 60Fe upper limit:     {spi_fe60:.1e} ph/cm2/s "
          f"(combined 1173+1332, 2 sigma)")
    print(f"  SPI Cygnus 26Al complex flux:    {spi_al26:.1e} ph/cm2/s "
          f"(winds + SNe; ours is a LOWER BOUND on it)")
    print(f"\n  saturation (primary arm, headline median): "
          f"26Al {primary.saturation_al26.median():.2f}, "
          f"60Fe {primary.saturation_fe60.median():.2f}")
    print(f"  supernovae at or below 25 Msun: {events_below_collapse} "
          f"of {total_events}")
    print(f"\n  T4  measured complex 26Al = "
          f"{complex_al26_msun_measured:.2e} Msun, not the 1 Msun WP8 §5 used;")
    print(f"      SN-only supplies {100 * al26_flux_fraction.min():.0f}-"
          f"{100 * al26_flux_fraction.max():.0f}% of the measured complex "
          f"1809 keV flux")
    print("\n  predictions:")
    for entry in predictions:
        print(f"    {entry['id']}  {entry['outcome']:4s}  {entry['statement']}")
        if "integrity_note" in entry:
            print(f"          ^ VACUOUS as written (see T4); corrected "
                  f"outcome {entry['outcome_against_corrected_denominator']}")
    print("\nwrote tables/wp11_isotope_forecast.csv")
    print("wrote tables/wp11_isotope_summary.csv")
    print("wrote provenance/wp11_isotope_forecast_execution.json")


if __name__ == "__main__":
    main()
