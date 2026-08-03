#!/usr/bin/env python3
"""Pre-registration: WP11 Part B, the 26Al / 60Fe forward prediction.

Written BEFORE any isotope mass or line flux is computed.  Per
tasks/wp11_bowshock_isotope_brief.md §4.2 step 1, fixing the YIELD BRANCH is
"the whole pre-registration for Part B": per-supernova 26Al and 60Fe yields for
>30 Msun progenitors are strongly model-dependent, so a declared literature
range is carried as a branch and is never adjusted after the fluxes are seen.

WHY THIS IS HARDER THAN IT LOOKS
--------------------------------
Every supernova in the WP7 ledger came from a progenitor above ~34 Msun on every
headline branch (WP9's C3 term), and for any black-hole threshold at or below
30 Msun the ledger returns exactly zero.  The ledger therefore lives entirely in
the mass range where the standard modern yield compilation refuses to explode
stars at all: Limongi & Chieffi (2018)'s Recommended scenario has every model
above 25 Msun collapse fully to the remnant, contributing only through its wind
-- so its explosive 26Al is zero and, because 60Fe reaches the ISM only through
the explosion, its 60Fe yield is IDENTICALLY ZERO above 25 Msun (Falla et al.
2025, §3).

That is not a reason to pick a friendlier table.  It is the single most
important fact about this forecast and it is carried as an explicit arm.  WP8
already excluded the zero-supernova branch on independent grounds -- PSR
J2032+4127 is a neutron star inside the association, so at least one star above
the local turnoff did explode -- which is why the null arm is reported rather
than adopted.  The 60Fe forecast is thereby a forward test of the very
proposition the paper turns on: whether very massive stars explode.

THE THREE ARMS
--------------
All three are published, none is invented here, and the two non-null arms come
from the SAME paper and differ only in the Wolf-Rayet mass-loss prescription --
so their spread is that paper's own stated systematic, not a spread this project
manufactured.

  LC06_NL      Limongi & Chieffi 2006 (ApJ 647, 483) Table 3.  Mass-resolved
               11-120 Msun, every mass explodes, Nugis & Lamers (2000) WR mass
               loss.  PRIMARY arm: the only published set giving explosive 26Al
               and 60Fe across the whole 30-120 Msun range the ledger populates.
  LC06_Langer  Same paper, Table 5: the Langer (1989) WNE+WCO mass-loss rate,
               40/60/80/120 Msun.  Stronger stripping leaves a small CO core, so
               the advanced-burning yields become nearly insensitive to initial
               mass.  LOW arm.
  LC18_REC     Limongi & Chieffi (2018) Recommended scenario as applied by Falla
               et al. (2025): M > 25 Msun collapses fully.  Explosive yields
               identically zero for every progenitor in this ledger.  NULL arm.

PARTIAL PRE-REGISTRATION, AND SAYING SO
---------------------------------------
The WP11 brief and the PROJECT_TRACE WP11 row already record a back-of-envelope
expectation: "a rough estimate puts the 60Fe flux near COSI's line sensitivity
with the alpha = 2.0 branches ~3x above the alpha = 2.3 ones".  That estimate was
made with a single yield value, no mass resolution, no decay accounting and no
arm structure, and none of the three arms above existed when it was written.  It
is nevertheless approximate prior knowledge of the answer and is disclosed here
rather than concealed.  Prediction I3's threshold is therefore set at 2x, BELOW
the remembered 3x, so that the prediction can fail against the estimate that
motivated the work.

POST-HOC DISCLOSURE (binding on the manuscript)
-----------------------------------------------
Unlike the WP8 markers, this comparison was NOT frozen at WP1.  It is a post-hoc
addition, pre-registered before scoring but chosen after the ledger existed.  The
manuscript must say so in the same sentence that introduces it (brief §5.3).

Output: provenance/wp11_isotope_prereg.json

Run:
  PYTHONPATH=scripts python3 scripts/wp11_isotope_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp8_crosschecks_prereg import AL26

# --------------------------------------------------------------- the isotopes
# 26Al mean lifetime is taken from the FROZEN WP8 constant so the two reports
# cannot disagree.  The exact half-life 0.717 Myr gives tau = 1.034 Myr, 1.5%
# lower; immaterial against a yield spread of order 50.
ISOTOPES = {
    "26Al": {
        "mass_number": 26,
        "half_life_Myr": 0.717,
        "mean_lifetime_Myr": AL26["mean_lifetime_Myr"],
        "mean_lifetime_source": (
            "frozen WP8 constant, provenance/wp8_crosschecks_prereg.json"
        ),
        "lines_keV": [1808.65],
        "photons_per_decay": {"1808.65": 0.9976},
        "produced_by_winds": True,
    },
    "60Fe": {
        "mass_number": 60,
        "half_life_Myr": 2.62,
        "half_life_source": (
            "Rugel et al. 2009, PRL 103, 072502; confirmed by Wallner et al. "
            "2015, PRL 114, 041101"
        ),
        "mean_lifetime_Myr": 2.62 / np.log(2.0),
        "mean_lifetime_source": "half-life / ln 2",
        "lines_keV": [1173.2, 1332.5],
        "photons_per_decay": {"1173.2": 0.999, "1332.5": 0.9998},
        "produced_by_winds": False,
        "why_it_matters": (
            "60Fe reaches the ISM only through the explosion.  Wolf-Rayet winds "
            "make none.  It is therefore the supernova-specific tracer and the "
            "ledger's cleanest future observable."
        ),
    },
}

# ------------------------------------------------------------- the yield arms
# Limongi & Chieffi 2006, ApJ 647, 483, Table 3 ("26Al and 60Fe yields"),
# transcribed from the published table.  Columns used:
#   26Al total, 26Al wind, 60Fe total.
# The SN-only 26Al yield is (total - wind) = the C-shell plus explosive Ne/C
# components; the wind component is excluded because it was released during the
# star's life, not by the explosion.  60Fe has no wind component at all.
LC06_MASS_MSUN = np.array(
    [11, 12, 13, 14, 15, 16, 17, 20, 25, 30, 35, 40, 60, 80, 120], dtype=float
)
LC06_AL26_TOTAL = np.array([
    1.60e-5, 2.11e-5, 2.45e-5, 1.04e-4, 1.32e-4, 6.80e-5, 6.87e-5, 5.43e-5,
    8.61e-5, 9.93e-5, 8.38e-5, 1.21e-4, 2.52e-4, 4.00e-4, 7.03e-4,
])
LC06_AL26_WIND = np.array([
    1.16e-11, 2.91e-11, 1.46e-10, 4.37e-10, 1.17e-09, 2.74e-09, 6.76e-09,
    4.32e-08, 3.93e-07, 2.39e-06, 1.14e-05, 2.06e-05, 6.94e-05, 1.32e-04,
    2.82e-04,
])
LC06_FE60_TOTAL = np.array([
    1.71e-6, 4.33e-6, 7.56e-5, 5.72e-6, 6.28e-6, 4.39e-6, 7.96e-6, 1.56e-5,
    3.69e-5, 1.49e-5, 4.03e-5, 6.23e-5, 2.27e-4, 7.55e-4, 9.93e-4,
])

# Limongi & Chieffi 2006 Table 5: models computed with the Langer (1989) mass
# loss rate in the WNE+WCO phases.  "26Al Cshell+expl." is already the SN-only
# quantity; "60Fe Total" is wholly explosive.
LC06L_MASS_MSUN = np.array([40, 60, 80, 120], dtype=float)
LC06L_AL26_SN = np.array([1.053e-4, 4.490e-5, 1.120e-4, 4.200e-5])
LC06L_FE60_SN = np.array([4.54e-5, 1.52e-5, 2.17e-5, 1.75e-5])

YIELD_ARMS = {
    "LC06_NL": {
        "role": "primary",
        "reference": "Limongi & Chieffi 2006, ApJ 647, 483, Table 3",
        "bibcode": "2006ApJ...647..483L",
        "mass_loss": (
            "Vink et al. 2000/2001 (BSG), de Jager et al. 1988 (RSG), "
            "Nugis & Lamers 2000 (WR)"
        ),
        "explodability": "every model 11-120 Msun explodes",
        "mass_grid_Msun": LC06_MASS_MSUN.tolist(),
        "al26_sn_only_Msun": (LC06_AL26_TOTAL - LC06_AL26_WIND).tolist(),
        "al26_total_Msun": LC06_AL26_TOTAL.tolist(),
        "al26_wind_Msun": LC06_AL26_WIND.tolist(),
        "fe60_sn_only_Msun": LC06_FE60_TOTAL.tolist(),
        "why_primary": (
            "the only published set that gives explosive 26Al AND 60Fe over the "
            "entire 30-120 Msun range this ledger populates"
        ),
    },
    "LC06_Langer": {
        "role": "low",
        "reference": "Limongi & Chieffi 2006, ApJ 647, 483, Table 5",
        "bibcode": "2006ApJ...647..483L",
        "mass_loss": "Langer 1989 rate in the WNE+WCO phases",
        "explodability": "every model explodes; only 40-120 Msun tabulated",
        "mass_grid_Msun": LC06L_MASS_MSUN.tolist(),
        "al26_sn_only_Msun": LC06L_AL26_SN.tolist(),
        "fe60_sn_only_Msun": LC06L_FE60_SN.tolist(),
        "why_carried": (
            "stronger WR stripping leaves a small CO core, so the "
            "advanced-burning yields become nearly independent of initial mass "
            "and fall by up to a factor 57 in 60Fe at 120 Msun.  Same authors, "
            "same code, same paper -- the spread is theirs, not ours."
        ),
    },
    "LC18_REC": {
        "role": "null",
        "reference": (
            "Limongi & Chieffi 2018, ApJS 237, 13, Recommended scenario, as "
            "applied by Falla et al. 2025, ApJ (arXiv:2508.07705) §2"
        ),
        "bibcode": "2018ApJS..237...13L",
        "mass_loss": "Vink/de Jager/Nugis & Lamers, as LC18",
        "explodability": "M > 25 Msun collapses fully; wind contribution only",
        "mass_grid_Msun": None,
        "al26_sn_only_Msun": 0.0,
        "fe60_sn_only_Msun": 0.0,
        "why_carried": (
            "it is the current standard compilation and it predicts EXACTLY "
            "ZERO supernova 26Al and 60Fe from this ledger, because every "
            "progenitor here is above its 25 Msun collapse threshold.  Reported, "
            "not adopted: WP8's pulsar already excludes a zero-supernova Cyg "
            "OB2.  Its role is to show what a COSI non-detection would mean."
        ),
    },
}

# Interpolation rule, fixed here so it cannot be tuned later.
INTERPOLATION = {
    "rule": "linear in log10(yield) against log10(initial mass)",
    "extrapolation": (
        "NONE.  Outside the tabulated mass range the endpoint value is held "
        "constant.  Progenitors in this ledger reach the 120 Msun IMF ceiling, "
        "which is also the top of the LC06 grid, so clipping bites only at the "
        "bottom of the LC06_Langer arm (below 40 Msun)."
    ),
    "why_log_log": (
        "both yields rise by more than an order of magnitude across the grid; "
        "linear interpolation in the raw values would be dominated by the top "
        "of the range"
    ),
}

# ------------------------------------------------------- the estimator itself
ESTIMATOR = {
    "primary": {
        "name": "decay-weighted accumulation",
        "formula": (
            "M_iso(now) = sum over the ledger's supernovae of "
            "y_iso(m_progenitor) * exp(-t_explosion / tau_iso), averaged over "
            "Monte Carlo iterations"
        ),
        "why_not_the_brief_formula": (
            "the brief (§4.2 step 2) writes 'rate x mean lifetime x yield', "
            "which is the STEADY-STATE limit.  Steady state requires the "
            "supernova rate to have been roughly constant for several mean "
            "lifetimes.  Cyg OB2's first supernova was about 1.4 Myr ago, "
            "against tau = 1.05 Myr for 26Al and 3.78 Myr for 60Fe, so 60Fe is "
            "nowhere near saturation and the steady-state formula would "
            "OVERSTATE it substantially.  The decay-weighted sum is the same "
            "physics without that assumption and reduces to the brief's formula "
            "in the limit the brief assumes."
        ),
        "mass_resolved": (
            "yes -- the yield is evaluated at each supernova's own progenitor "
            "mass, drawn by the WP7 engine, not at an ensemble average"
        ),
    },
    "secondary_reported": {
        "name": "steady-state approximation",
        "formula": "M_iso = R_SN(now) * tau_iso * <y_iso>",
        "purpose": (
            "reported alongside the primary as the SATURATION RATIO "
            "M_primary / M_steady_state, which quantifies how far from steady "
            "state this young association is.  It is a diagnostic, never the "
            "quoted forecast."
        ),
    },
    "flux": {
        "formula": "F = (M_iso / (A * u)) / tau_iso * p_gamma / (4 pi d^2)",
        "distance": (
            "the WP3/WP2 distance artifact, imported as wp3_common.D_KPC "
            "(1.6245 kpc, provenance/wp2_membership_manifest.json "
            "one_component_mean_kpc) -- the artifact, not a typed number"
        ),
        "point_source_assumption": (
            "Cyg OB2 subtends about 1 deg, well inside the ~3 deg resolution of "
            "SPI and COSI, so a point-source flux is the right comparison"
        ),
    },
}

# ------------------------------------------------- what it is compared against
INSTRUMENTS = {
    "SPI_Cygnus_60Fe_upper_limit": {
        "value_ph_cm2_s": 1.6e-5,
        "line": "1173 + 1332 keV, Cygnus region",
        "source": (
            "Martin et al. 2009, A&A 506, 703 -- FROZEN AT WP1, see "
            "wp1_sn_markers.md"
        ),
        "nature": "upper limit",
    },
    "SPI_Cygnus_26Al_complex_flux": {
        "value_ph_cm2_s": AL26["complex_flux_ph_cm2_s"],
        "error_ph_cm2_s": AL26["complex_flux_error"],
        "line": "1809 keV, component attributed to the Cygnus complex",
        "source": AL26["source"] + " -- FROZEN AT WP1",
        "nature": "measurement of winds PLUS supernovae, complex-wide",
    },
    "COSI_narrow_line_3sigma_2yr": {
        "value_ph_cm2_s": 3.0e-6,
        "lines_keV": [1173.0, 1333.0, 1809.0],
        "source": (
            "Tomsick et al. 2023, PoS(ICRC2023)745, Table 1 "
            "(arXiv:2308.12362) -- 3-sigma narrow-line point-source "
            "sensitivity in 2 years of survey observations"
        ),
        "launch": "planned 2027",
        "nature": "expected sensitivity, identical at all three lines",
    },
    "galactic_ratio_context": {
        "fe60_over_al26": 0.184,
        "error": 0.042,
        "source": "Wang et al. 2020, ApJ 889, 169",
        "use": "context only; never used to set or rescale any yield here",
    },
}

# ---------------------------------------------------------- the branch set
HEADLINE_ALPHAS = (2.0, 2.3)
EXCLUDED_ALPHA = 2.6
PRIMARY_EXPLODABILITY = "all_explode"

# ------------------------------------------------------- detectability rule
DETECTABILITY_RULE = {
    "DETECTABLE": "flux >= the instrument value on EVERY headline branch",
    "MARGINAL": "the headline branches straddle the instrument value",
    "BELOW_REACH": "flux < the instrument value on EVERY headline branch",
    "applied": (
        "mechanically, per isotope and per yield arm, to the numbers as "
        "computed.  A verdict of BELOW_REACH is a forecast and is reported as "
        "one; it is not a failure and is not to be rescued by re-picking a "
        "yield arm."
    ),
}


def main() -> None:
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp11_isotope_prereg.py",
        "status": "PREREGISTERED",
        "work_package": "WP11 Part B",
        "brief": "tasks/wp11_bowshock_isotope_brief.md §4",
        "objective": (
            "convert the WP8 26Al consistency band into a falsifiable forward "
            "prediction: the ledger's SN-only 26Al and 60Fe masses and line "
            "fluxes, per headline branch and per yield arm, against the frozen "
            "SPI Cygnus limits and COSI's expected narrow-line sensitivity"
        ),
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "partial_preregistration_disclosure": {
            "what_is_already_known": (
                "the WP11 brief and the PROJECT_TRACE WP11 row record a "
                "back-of-envelope expectation that the 60Fe flux lands near "
                "COSI's line sensitivity with the alpha = 2.0 branches about 3x "
                "above the alpha = 2.3 ones"
            ),
            "how_rough_it_was": (
                "a single yield value, no mass resolution, no decay accounting, "
                "no arm structure.  None of the three arms fixed above existed "
                "when it was written, and the dominant systematic identified "
                "here -- the WR mass-loss prescription, worth a factor up to 57 "
                "in 60Fe -- was not in it at all."
            ),
            "how_this_is_handled": (
                "prediction I3's threshold is set at 2x, BELOW the remembered "
                "3x, so that I3 can fail against the estimate that motivated "
                "the work.  Prediction I5's detectability verdict is left to a "
                "mechanical rule rather than being asserted."
            ),
            "what_would_have_been_illegitimate": (
                "choosing the yield arm, or the mass range over which it is "
                "interpolated, after seeing which choice puts the flux above "
                "COSI's sensitivity"
            ),
        },
        "isotopes": {
            name: {k: v for k, v in spec.items()} for name, spec in ISOTOPES.items()
        },
        "yield_branch": {
            "principle": (
                "a declared literature RANGE carried as a branch, never a "
                "single value, fixed before any flux is computed and not "
                "adjusted afterwards (brief §4.2 step 1)"
            ),
            "arms": YIELD_ARMS,
            "interpolation": INTERPOLATION,
            "sn_only_definition": (
                "26Al: total minus wind, i.e. the C-shell plus explosive Ne/C "
                "components, because the wind was released during the star's "
                "life and not by its explosion.  60Fe: the total, which is "
                "wholly explosive."
            ),
        },
        "estimator": ESTIMATOR,
        "instruments": INSTRUMENTS,
        "branch_set": {
            "headline_alphas": list(HEADLINE_ALPHAS),
            "excluded_alpha_reported_not_deleted": EXCLUDED_ALPHA,
            "explodability": PRIMARY_EXPLODABILITY,
            "expected_headline_branches": 36,
            "inherited_from": "provenance/wp7_alpha_headline_adoption_prereg.json",
        },
        "detectability_rule": DETECTABILITY_RULE,
        "predictions": [
            {
                "id": "I1",
                "statement": (
                    "on the primary arm the SN-only 26Al mass is at least 100x "
                    "below the ~1 Msun inferred for the whole Cygnus complex, "
                    "on every headline branch"
                ),
                "reasoning": (
                    "WP8 §5 already put the ordering near 1000x with a single "
                    "1e-4 Msun yield and the baseline rate.  Mass-resolved "
                    "yields at 30-120 Msun are several times larger and the "
                    "alpha = 2.0 branches raise the rate, so the margin should "
                    "shrink -- but nowhere near closing."
                ),
                "falsifies_if": (
                    "the SN-only 26Al mass exceeds 1% of the complex-wide "
                    "inventory on any headline branch"
                ),
                "if_falsified": (
                    "the WP8 consistency band would be under real strain and "
                    "the tension is REPORTED as a tension; nothing upstream is "
                    "retuned"
                ),
            },
            {
                "id": "I2",
                "statement": (
                    "the predicted 60Fe flux on the primary arm lies below the "
                    "frozen SPI Cygnus upper limit of 1.6e-5 ph/cm2/s on every "
                    "headline branch"
                ),
                "reasoning": (
                    "if it did not, existing INTEGRAL data would already "
                    "exclude the ledger, and this would be a result about the "
                    "ledger rather than a forecast"
                ),
                "falsifies_if": "any headline branch exceeds the SPI limit",
                "if_falsified": (
                    "recorded as a genuine tension between the ledger and an "
                    "existing measurement, and reported as such"
                ),
            },
            {
                "id": "I3",
                "statement": (
                    "on the primary arm the alpha = 2.0 headline branches "
                    "predict a median 60Fe flux at least 2x that of the "
                    "alpha = 2.3 branches"
                ),
                "reasoning": (
                    "alpha is the axis the WP9 verdict hinges on and the only "
                    "one Gaia DR4 will not settle.  If 60Fe separates the two "
                    "arms, COSI can discriminate what DR4 cannot -- which is "
                    "the whole reason Part B was adopted."
                ),
                "falsifies_if": "the ratio of medians is below 2",
                "if_falsified": (
                    "the discriminating-power claim is withdrawn and the "
                    "forecast is reported as a consistency check only"
                ),
            },
            {
                "id": "I4",
                "statement": (
                    "the spread in predicted 60Fe flux BETWEEN yield arms "
                    "exceeds the spread across headline branches WITHIN the "
                    "primary arm"
                ),
                "reasoning": (
                    "the brief's own risk register says the yield-model spread "
                    "may swamp the forecast, and that the spread IS the result. "
                    "This states it as a testable claim before the numbers "
                    "exist."
                ),
                "falsifies_if": (
                    "the within-arm branch spread is the larger of the two"
                ),
                "if_falsified": (
                    "better than expected: the forecast would then be limited "
                    "by the census rather than by nuclear astrophysics, and "
                    "the paper says so"
                ),
            },
            {
                "id": "I5",
                "statement": (
                    "the COSI detectability verdict for 60Fe on the primary arm "
                    "is MARGINAL under the mechanical rule -- the headline "
                    "branches straddle 3.0e-6 ph/cm2/s"
                ),
                "reasoning": (
                    "this is the prediction most exposed to the disclosed prior "
                    "estimate above, and is deliberately stated so that either "
                    "DETECTABLE or BELOW_REACH falsifies it"
                ),
                "falsifies_if": (
                    "the verdict comes out DETECTABLE or BELOW_REACH"
                ),
                "if_falsified": (
                    "the computed verdict stands and is reported; a stated "
                    "non-detectability is still a forecast (brief §4.2 step 3)"
                ),
            },
        ],
        "prohibitions": {
            "no_inversion": (
                "the WP8 prohibition stands in full: the measured complex-wide "
                "26Al flux is NEVER inverted into a supernova count.  The "
                "prediction runs forward only.  The SN-only component is stated "
                "as a LOWER BOUND on the complex-wide signal, because "
                "Wolf-Rayet winds dominate 26Al at this age and contribute no "
                "60Fe."
            ),
            "nothing_upstream_modified": (
                "WP11 consumes the frozen repair_v7 chain read-only.  Verified "
                "against the audit.py inventory before and after (gate G11d)."
            ),
            "failures_recorded_as_failures": (
                "project rule: a failed prediction is recorded as failed and is "
                "never reinterpreted."
            ),
            "no_yield_reselection": (
                "the three arms above are final.  No arm may be added, dropped "
                "or reweighted after a flux has been computed."
            ),
        },
        "manuscript_obligations": {
            "post_hoc_disclosure": (
                "MANDATORY.  The manuscript must state, in the same sentence "
                "that introduces this comparison, that it was not frozen at WP1 "
                "-- it is post-hoc, pre-registered before scoring but chosen "
                "after the ledger existed.  The credibility of the WP8 layer "
                "rests on its freeze and the two must not be blurred."
            ),
            "no_hand_typed_numbers": (
                "every quoted quantity enters as a macro via "
                "scripts/wp10_numbers.py, reading this work's artifacts through "
                "scripts/wp10_inputs.py"
            ),
            "placement": (
                "one Discussion paragraph adjacent to the DR4 outlook, plus one "
                "sentence in the 26Al cross-check paragraph"
            ),
        },
        "what_part_b_cannot_do": [
            "it cannot measure a yield; it carries published ones as a branch",
            "it cannot separate Cyg OB2 from the wider Cygnus complex in any "
            "existing 26Al measurement -- SPI/COMPTEL resolution is ~3 deg",
            "it cannot predict the WIND 26Al, which dominates the complex-wide "
            "signal at this age and is outside the ledger's scope",
            "it cannot decide between the three degenerate WP8 readings of the "
            "pulsar; it only shows that the null arm and the other two make "
            "observationally distinguishable predictions",
        ],
    }
    w.write_json(w.PROVENANCE / "wp11_isotope_prereg.json", record)

    print("WP11 Part B — 26Al / 60Fe forward prediction — pre-registered\n")
    print("  yield arms (fixed BEFORE any flux):")
    for name, arm in YIELD_ARMS.items():
        print(f"    {name:12s} [{arm['role']:7s}] {arm['reference']}")
    print(f"\n  estimator: {ESTIMATOR['primary']['name']}, mass-resolved")
    print(f"  distance:  wp3_common.D_KPC (artifact)")
    print(f"\n  compared against:")
    print(f"    SPI Cygnus 60Fe upper limit  1.6e-05 ph/cm2/s  (frozen at WP1)")
    print(f"    SPI Cygnus 26Al complex flux "
          f"{AL26['complex_flux_ph_cm2_s']:.1e} ph/cm2/s  (frozen at WP1)")
    print(f"    COSI 3-sigma narrow line     3.0e-06 ph/cm2/s  (2 yr survey)")
    print(f"\n  headline set: alpha in {list(HEADLINE_ALPHAS)}, "
          f"{PRIMARY_EXPLODABILITY}, 36 branches")
    print("\n  predictions:")
    for entry in record["predictions"]:
        print(f"    {entry['id']}  {entry['statement'][:66]}")
    print("\n  26Al is NEVER inverted into a supernova count (WP8 prohibition).")
    print("\nwrote provenance/wp11_isotope_prereg.json")


if __name__ == "__main__":
    main()
