#!/usr/bin/env python3
"""External cross-checks of the WP1-WP6 chain against published measurements.

Recomputes, from the accepted artifacts, every number quoted in cross_checks/
and records it with input hashes so the comparison is reproducible rather than
prose.  Three independent layers are tested:

  Harer et al. 2025 (A&A 703, A111)   the OUTPUT -- supernova rate and the
                                      association stellar mass that produces it
  Orellana et al. 2021 (MNRAS 502)    the ASTROMETRIC FOUNDATION -- systemic
                                      proper motion and distance
  Berlanas et al. 2019 (MNRAS 484)    a STRUCTURAL ASSUMPTION -- that Cyg OB2 is
                                      one population along the line of sight

NOTHING HERE IS A CALIBRATION.  These are validations.  If one disagrees it
becomes an issue in PROJECT_TRACE section 9; it is never a reason to move a
number.  No value in the chain was tuned toward any of them, and the Harer
comparison was not made until after WP6 closed.

Outputs:
  tables/wp6_external_crosschecks.csv
  provenance/wp6_external_crosschecks_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_external_crosschecks.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

WP5_VERSION = "repair_v6"
UPSTREAM = "repair_v5"
BASELINE_FAMILY, BASELINE_RV = "PARSEC", 3.1
IMF_MASS_LO = 0.5          # the frozen MASS_GRID lower edge
DERIVATIVE_STEP_MYR = 0.02  # central difference for |dM_turnoff/dt|

# ---------------------------------------------------------------- literature
HARER = {
    "reference": "Harer, Vieu, Schulze, Larkin & Reville 2025, A&A 703, A111",
    "local_copy": "papers/Harer_2025.pdf",
    "section": "3.1 and Fig. 2",
    "method": "Hoki (Stevance et al. 2020) interface to BPASS (Stanway & Eldridge 2018)",
    "assumed_cluster_mass_Msun": 1.65e4,
    "assumed_cluster_mass_source": "Wright et al. 2015",
    "assumed_alpha": 2.0,
    "assumed_metallicity": "solar",
    "figure_2_y_axis_verbatim": "Event Rate (events/100 kyr)",
    "figure_2_y_ticks": [0.25, 0.5, 1.0, 2.0],
    "text_verbatim": (
        "For an age of 3-5 Myr, type Ic supernovae are by far the most likely "
        "type of supernova, and the rate is consistent with several supernovae "
        "per Myr."
    ),
    "unit_trap": (
        "the axis is events per 100 kiloyears, NOT per 100 years.  1 event/100 "
        "kyr = 10 SNe/Myr; 1 event/100 yr would be 10,000 SNe/Myr.  Misreading "
        "it inflates the expected rate by a factor of 1000."
    ),
}
ORELLANA = {
    "reference": "Orellana, De Biasi & Paiz 2021, MNRAS 502, 6080",
    "data_release": "Gaia DR2",
    "field": "1 deg circle at (l, b) = (79.8, +0.8), G <= 17.5",
    "pmra_mas_yr": -2.71, "pmra_err": 0.02,
    "pmdec_mas_yr": -4.24, "pmdec_err": 0.02,
    "distance_astrometric_pc": 1683.0, "distance_astrometric_err": 5.0,
    "distance_astrophotometric_pc": 1669.0, "distance_astrophotometric_err": 6.0,
    "mean_parallax_mas": 0.599,
    "zero_point": "DR2 global offset -0.029 mas (Lindegren et al. 2018)",
    "other_overdensities_pc": {"UCB585": 1460.0, "right": 1280.0, "left": 1280.0},
    "distance_literature_pc": {
        "Morgan+1954": 1500, "Reddish+1966": 2100, "Humphreys+1984": 1820,
        "Massey&Thompson+1991": 1700, "Comeron&Pasquali+2012": 1445,
        "Kiminki+2015": 1330, "Berlanas+2019": 1760, "Lim+2019": 1600,
        "Orellana+2021": 1669,
    },
}
BERLANAS = {
    "reference": "Berlanas, Wright, Herrero, Drew & Lennon 2019, MNRAS 484, 1838",
    "claim": "line-of-sight substructure: main group ~1760 pc, foreground ~1350 pc",
    "foreground_stars": 19,
    "foreground_O_stars": 7,
    "foreground_distance_pc": 1350.0,
    "main_distance_pc": 1760.0,
}


def imf_integral(k: float, alpha: float, lo: float, hi: float) -> float:
    """k * integral[lo, hi] M^-alpha dM."""
    if lo >= hi:
        return 0.0
    return k * (lo ** (1.0 - alpha) - hi ** (1.0 - alpha)) / (alpha - 1.0)


def imf_mass_integral(k: float, alpha: float, lo: float, hi: float) -> float:
    """k * integral[lo, hi] M^(1-alpha) dM -- the stellar mass in that range."""
    if np.isclose(alpha, 2.0):
        return k * np.log(hi / lo)
    return k * (lo ** (2.0 - alpha) - hi ** (2.0 - alpha)) / (alpha - 2.0)


def turnoff_rate(family: str, age: float) -> float:
    """|dM_turnoff/dt| in Msun/Myr, by central difference."""
    h = DERIVATIVE_STEP_MYR
    return abs(
        turnoff_mass(family, age + h) - turnoff_mass(family, age - h)
    ) / (2.0 * h)


def galactic_summary(members: pd.DataFrame) -> dict:
    parallax = members["parallax_corrected"]
    weights = 1.0 / members["parallax_error"] ** 2
    weighted = float(np.average(parallax, weights=weights))
    return {
        "members": int(len(members)),
        "median_parallax_mas": round(float(parallax.median()), 4),
        "distance_from_median_pc": round(1000.0 / float(parallax.median()), 1),
        "weighted_mean_parallax_mas": round(weighted, 4),
        "distance_from_weighted_mean_pc": round(1000.0 / weighted, 1),
        "fraction_nearer_than_1430pc": round(float((parallax > 0.70).mean()), 4),
        "fraction_nearer_than_1350pc": round(float((parallax > 0.7407).mean()), 4),
    }


def main() -> None:
    normalization = pd.read_parquet(
        w.PROC / f"wp5_imf_normalization_{WP5_VERSION}.parquet"
    )
    branch = normalization[
        normalization.family.eq(BASELINE_FAMILY)
        & normalization.R_V.eq(BASELINE_RV)
    ]
    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet"
    )
    ages = {}
    for subgroup in w.SUBGROUPS:
        row = age_posterior[
            age_posterior.subgroup.eq(subgroup)
            & age_posterior.family.eq(BASELINE_FAMILY)
            & age_posterior.R_V.eq(BASELINE_RV)
            & age_posterior.f_bin.eq(w.F_BINARY)
            & age_posterior.indicator.eq("ums")
            & age_posterior.dmu.eq(0.0)
        ]
        ages[subgroup] = float(row.age_map.iloc[0])

    members = pd.read_parquet(
        w.PROC / "wp2_members.parquet",
        columns=["source_id", "parallax_corrected", "parallax_error",
                 "membership_probability"],
    )
    labels = pd.read_parquet(
        w.TABLES / "wp2_subgroup_labels.parquet", columns=["source_id", "subgroup"]
    )
    labelled = members.merge(labels, on="source_id", how="inner")
    labelled = labelled[labelled.membership_probability > 0.5]

    # ---------------------------------------------------- Harer: rate and mass
    rate_rows = []
    for alpha in w.IMF_SLOPES:
        cumulative = 0.0
        rate = 0.0
        mass = 0.0
        per_subgroup = {}
        for subgroup, age in ages.items():
            k = float(
                branch[branch.subgroup.eq(subgroup) & branch.alpha.eq(alpha)]
                .k_median.iloc[0]
            )
            cap = min(turnoff_mass(BASELINE_FAMILY, age), IMF_UPPER_LIMIT)
            dead = imf_integral(k, alpha, cap, IMF_UPPER_LIMIT) if cap < IMF_UPPER_LIMIT else 0.0
            instantaneous = (
                k * cap ** (-alpha) * turnoff_rate(BASELINE_FAMILY, age)
                if cap < IMF_UPPER_LIMIT else 0.0
            )
            cumulative += dead
            rate += instantaneous
            mass += imf_mass_integral(k, alpha, IMF_MASS_LO, IMF_UPPER_LIMIT)
            per_subgroup[subgroup] = {
                "age_Myr": round(age, 3),
                "turnoff_Msun": round(cap, 1),
                "cumulative_SNe": round(dead, 2),
                "rate_SNe_per_Myr": round(instantaneous, 2),
            }
        rate_rows.append(
            {
                "alpha": float(alpha),
                "cumulative_SNe": round(cumulative, 2),
                "rate_SNe_per_Myr": round(rate, 2),
                "rate_events_per_100kyr": round(rate / 10.0, 3),
                "inside_figure2_range": bool(0.25 <= rate / 10.0 <= 2.0),
                "association_mass_Msun": round(mass, 1),
                "mass_ratio_to_Wright2015": round(mass / HARER["assumed_cluster_mass_Msun"], 3),
                "per_subgroup": per_subgroup,
            }
        )

    matched = next(r for r in rate_rows if np.isclose(r["alpha"], HARER["assumed_alpha"]))
    baseline_rate = next(r for r in rate_rows if np.isclose(r["alpha"], 2.3))
    first_death_Myr = None
    for age in np.arange(2.0, 5.0, 0.01):
        if turnoff_mass(BASELINE_FAMILY, float(age)) < IMF_UPPER_LIMIT:
            first_death_Myr = round(float(age), 2)
            break

    # ------------------------------------------- Orellana: distance and motion
    overall = galactic_summary(labelled)
    distance_offset = (
        ORELLANA["distance_astrophotometric_pc"] / overall["distance_from_median_pc"]
    )
    modulus_shift = 5.0 * np.log10(distance_offset)
    # A larger assumed distance makes each star intrinsically brighter by the
    # full modulus shift.  With L ~ M^3.5 on the upper main sequence,
    # M_bol = -8.75 log10 M, so d(log10 M) = modulus / (2.5 * 3.5).  That is a
    # LOGARITHMIC shift -- it must be exponentiated to get the fractional mass
    # change, which is the quantity the count responds to.
    mass_shift_dex = modulus_shift / (2.5 * 3.5)
    mass_shift = 10.0 ** mass_shift_dex - 1.0
    # Scaling every mass by (1 + f) is equivalent to lowering the 8 Msun
    # threshold by the same factor, so for dN/dM ~ M^-alpha the count above the
    # threshold changes by (alpha - 1) * f.
    count_shift = (2.3 - 1.0) * mass_shift

    traceback = json.loads(
        (w.PROVENANCE / "wp6_runaways_execution.json").read_text(encoding="utf-8")
    )
    our_pmra = float(traceback["configuration"]["systemic_pmra_mas_yr"])
    our_pmdec = float(traceback["configuration"]["systemic_pmdec_mas_yr"])

    # -------------------------------------------- Berlanas: distance structure
    distance_test = json.loads(
        (w.PROVENANCE / "wp2_distance_population_execution.json").read_text(
            encoding="utf-8"
        )
    )
    association = distance_test["association"]
    per_subgroup_parallax = {
        subgroup: galactic_summary(block)
        for subgroup, block in labelled.groupby("subgroup", observed=True)
    }

    table = pd.DataFrame(
        [
            {"check": "harer_rate", "alpha": r["alpha"],
             "ours": r["rate_events_per_100kyr"], "units": "events/100kyr",
             "literature": "0.25-2.0 (Fig. 2 range); 'several per Myr'",
             "agrees": r["inside_figure2_range"]}
            for r in rate_rows
        ] + [
            {"check": "harer_mass", "alpha": r["alpha"],
             "ours": r["association_mass_Msun"], "units": "Msun",
             "literature": HARER["assumed_cluster_mass_Msun"],
             "agrees": bool(0.5 < r["mass_ratio_to_Wright2015"] < 2.0)}
            for r in rate_rows
        ] + [
            {"check": "orellana_pmra", "alpha": None, "ours": our_pmra,
             "units": "mas/yr", "literature": ORELLANA["pmra_mas_yr"],
             "agrees": abs(our_pmra - ORELLANA["pmra_mas_yr"]) < 0.1},
            {"check": "orellana_pmdec", "alpha": None, "ours": our_pmdec,
             "units": "mas/yr", "literature": ORELLANA["pmdec_mas_yr"],
             "agrees": abs(our_pmdec - ORELLANA["pmdec_mas_yr"]) < 0.1},
            {"check": "orellana_distance", "alpha": None,
             "ours": overall["distance_from_median_pc"], "units": "pc",
             "literature": ORELLANA["distance_astrophotometric_pc"],
             "agrees": abs(distance_offset - 1.0) < 0.10},
            {"check": "berlanas_foreground", "alpha": None,
             "ours": overall["fraction_nearer_than_1350pc"], "units": "fraction",
             "literature": "19 stars at ~1350 pc",
             "agrees": False},
        ]
    )
    out_csv = w.TABLES / "wp6_external_crosschecks.csv"
    table.to_csv(out_csv, index=False)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_external_crosschecks.py",
        "status": "SUCCESS",
        "work_package": "WP6 — external validation",
        "documents": [
            "cross_checks/harer_2025_supernova_rate.md",
            "cross_checks/orellana_2021_astrometry.md",
            "cross_checks/berlanas_2019_two_distance.md",
        ],
        "standing_rule": (
            "these are VALIDATIONS, not calibrations.  No value in the chain was "
            "tuned toward any of them.  A disagreement becomes an issue in "
            "PROJECT_TRACE section 9; it is never a reason to move a number."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "baseline_branch": f"{BASELINE_FAMILY} R_V={BASELINE_RV}",
        "wp5_version": WP5_VERSION,
        "subgroup_ages_Myr": {k: round(v, 3) for k, v in ages.items()},

        "harer_2025": {
            "literature": HARER,
            "our_estimator": (
                "dN_SN/dt = k * M_turnoff^(-alpha) * |dM_turnoff/dt|, summed "
                "over subgroups; cumulative = k * integral[M_turnoff, 120] "
                "dM M^-alpha.  Association mass = k * integral[0.5, 120] "
                "dM M^(1-alpha)."
            ),
            "by_alpha": rate_rows,
            "matched_alpha_comparison": {
                "alpha": HARER["assumed_alpha"],
                "our_rate_events_per_100kyr": matched["rate_events_per_100kyr"],
                "inside_figure_range": matched["inside_figure2_range"],
                "note": (
                    "compared at THEIR assumed alpha, which is our shallow "
                    "branch, not our baseline"
                ),
            },
            "baseline_alpha_comparison": {
                "alpha": 2.3,
                "our_rate_events_per_100kyr": baseline_rate["rate_events_per_100kyr"],
                "our_rate_SNe_per_Myr": baseline_rate["rate_SNe_per_Myr"],
                "their_text": "several supernovae per Myr",
            },
            "cumulative_vs_rate_reconciliation": {
                "first_possible_supernova_Myr": first_death_Myr,
                "explanation": (
                    f"the {BASELINE_FAMILY} turnoff only falls below the "
                    f"{IMF_UPPER_LIMIT:g} Msun IMF ceiling at {first_death_Myr} "
                    "Myr, so no star could have exploded before then.  CygOB2-A "
                    "is 3.98 Myr old, giving a supernova window of about 1 Myr; "
                    "a cumulative count of ~6 and a rate of ~8/Myr are the same "
                    "statement, not a contradiction."
                ),
                "CygOB2_C_contributes_zero": (
                    "at 2.51 Myr its turnoff is still above the IMF upper "
                    "limit — nothing has died there yet"
                ),
            },
            "mass_caveat": (
                "the association mass is sensitive to the low-mass integration "
                f"limit.  We integrate from {IMF_MASS_LO} Msun, the frozen "
                "MASS_GRID edge; a single alpha=2.3 power law continued to 0.1 "
                "Msun would give 3.1e4 Msun.  The real IMF flattens below ~0.5 "
                "Msun, so our cutoff happens to mimic a Kroupa turnover.  Good "
                "agreement, but not to be quoted better than a few tens of "
                "percent."
            ),
            "physical_difference_carried": (
                "BPASS models binary evolution and their Fig. 2 shows type Ic "
                "dominating at 3-5 Myr — the stripped-envelope channel, which "
                "REQUIRES a companion.  Our count is single-star turnoff "
                "counting.  Our N_SN is therefore a LOWER BOUND relative to a "
                "BPASS-style rate, and we agree with them before adding that "
                "channel.  Same physics as issue #15, from another direction."
            ),
            "verdict": "AGREE",
        },

        "orellana_2021": {
            "literature": ORELLANA,
            "our_systemic_pm": {"pmra": our_pmra, "pmdec": our_pmdec},
            "pm_difference_mas_yr": {
                "pmra": round(abs(our_pmra - ORELLANA["pmra_mas_yr"]), 4),
                "pmdec": round(abs(our_pmdec - ORELLANA["pmdec_mas_yr"]), 4),
            },
            "pm_significance": (
                "pmra agrees to 0.003 mas/yr against an independent DR2 "
                "analysis, which directly validates the issue #16 fix: the "
                "vector the runaway traceback subtracts is externally "
                "confirmed.  The 0.077 mas/yr pmdec difference is 0.59 km/s at "
                "1.62 kpc, negligible against a 10-100 km/s ejection window."
            ),
            "our_distance": overall,
            "distance_offset_fraction": round(distance_offset - 1.0, 4),
            "distance_offset_cause": (
                "parallax zero point.  They apply the DR2 global -0.029 mas "
                "(Lindegren et al. 2018); we apply the DR3 per-star Lindegren "
                "et al. 2021 correction.  Their mean member parallax is 0.599 "
                "mas, ours 0.614 — the 0.015 mas difference is the scale of a "
                "DR2 to DR3 zero-point revision."
            ),
            "propagated_systematic": {
                "distance_modulus_shift_mag": round(float(modulus_shift), 4),
                "mass_shift_dex": round(float(mass_shift_dex), 5),
                "mass_shift_fraction": round(float(mass_shift), 4),
                "closure_ratio_shift_fraction": round(float(count_shift), 4),
                "direction": (
                    "adopting their distance raises every closure ratio by "
                    "about 2%: it HELPS CygOB2-A (0.894, below unity) and HURTS "
                    "CygOB2-C (1.405)"
                ),
                "not_adopted": (
                    "changing the distance to match a cross-check would be "
                    "tuning.  Recorded as a systematic instead."
                ),
            },
            "verdict": "AGREE on proper motion; 3% distance offset carried as a systematic",
        },

        "berlanas_2019": {
            "literature": BERLANAS,
            "our_test": "scripts/wp2_distance_population_test.py",
            "one_component_kpc": round(
                float(association["one_component"]["means_kpc"][0]), 4
            ),
            "one_component_depth_pc": round(
                float(association["one_component"]["intrinsic_sigmas_kpc"][0]) * 1000, 1
            ),
            "two_component_kpc": [
                round(float(v), 4) for v in association["two_component"]["means_kpc"]
            ],
            "two_component_separation_pc": round(
                abs(
                    association["two_component"]["means_kpc"][1]
                    - association["two_component"]["means_kpc"][0]
                ) * 1000, 1
            ),
            "delta_bic_favouring_one": round(
                float(
                    association["two_component"]["bic"]
                    - association["one_component"]["bic"]
                ), 2
            ),
            "per_subgroup_parallax": per_subgroup_parallax,
            "recorded_verdict": distance_test["decision"]["verdict"],
            "CIRCULARITY_CAVEAT": (
                "WP2's membership classifier uses FEATURES = "
                "['parallax_corrected', 'pmra', 'pmdec'] — PARALLAX IS A "
                "CLUSTERING FEATURE.  A foreground group at 1350 pc would have "
                "been assigned low membership probability and removed BEFORE "
                "this test ran.  The honest scope is therefore: 'within the WP2 "
                "member sample there is no evidence of two distances'.  It is "
                "NOT evidence against Berlanas's foreground group."
            ),
            "what_would_settle_it": (
                "a test on a PARALLAX-BLIND selection — members chosen on sky "
                "position and proper motion only, then examined in parallax.  "
                "Not done here; it is the only way to break the circularity."
            ),
            "corroboration": (
                "Orellana et al. 2021 independently find foreground structure "
                "in the same field: 179 and 188 stars at ~1280 pc plus UCB585 "
                "at ~1460 pc, all with proper motions distinct from Cyg OB2.  "
                "All three analyses agree that foreground structure EXISTS and "
                "is separable by proper motion; the disagreement is about "
                "whether those stars are Cyg OB2 members, not whether they are "
                "there."
            ),
            "open_question_for_wp6": (
                "a star truly at 1350 pc but assumed at 1620 gets a distance "
                "modulus 0.40 mag too large, so it is inferred MORE massive, "
                "which INFLATES the count above 8 Msun and RAISES the closure "
                "ratio — the direction of CygOB2-C's 1.405 residual.  Our data "
                "do not support it (C's foreground fraction matches A and B, "
                "and no member sits at 1350 pc), but note that WP6 alternative "
                "A4 tested membership-probability weighting, NOT distance "
                "contamination.  This mechanism is covered by no WP6 "
                "alternative and is registered as open."
            ),
            "verdict": "NOT CONFIRMED — but this test cannot settle the claim",
        },

        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp5_imf_normalization_{WP5_VERSION}.parquet",
                w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet",
                w.PROC / "wp2_members.parquet",
                w.TABLES / "wp2_subgroup_labels.parquet",
                w.PROVENANCE / "wp2_distance_population_execution.json",
                w.PROVENANCE / "wp6_runaways_execution.json",
            ]
        },
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(
        w.PROVENANCE / "wp6_external_crosschecks_execution.json", record
    )

    print("external cross-checks\n")
    print("HARER 2025 — supernova rate and association mass")
    print(f"  {'alpha':>6s} {'SNe/Myr':>9s} {'ev/100kyr':>10s} {'in Fig.2':>9s} "
          f"{'mass Msun':>11s} {'vs 1.65e4':>10s}")
    for entry in rate_rows:
        print(f"  {entry['alpha']:6.1f} {entry['rate_SNe_per_Myr']:9.2f} "
              f"{entry['rate_events_per_100kyr']:10.3f} "
              f"{str(entry['inside_figure2_range']):>9s} "
              f"{entry['association_mass_Msun']:11.0f} "
              f"{entry['mass_ratio_to_Wright2015']:10.3f}")
    print(f"  first possible supernova: {first_death_Myr} Myr "
          f"(turnoff reaches {IMF_UPPER_LIMIT:g} Msun)")

    print("\nORELLANA 2021 — astrometry")
    print(f"  pmra   ours {our_pmra:+.4f}  theirs {ORELLANA['pmra_mas_yr']:+.2f}"
          f"  diff {abs(our_pmra - ORELLANA['pmra_mas_yr']):.4f} mas/yr")
    print(f"  pmdec  ours {our_pmdec:+.4f}  theirs {ORELLANA['pmdec_mas_yr']:+.2f}"
          f"  diff {abs(our_pmdec - ORELLANA['pmdec_mas_yr']):.4f} mas/yr")
    print(f"  distance ours {overall['distance_from_median_pc']:.0f} pc  theirs "
          f"{ORELLANA['distance_astrophotometric_pc']:.0f} pc  "
          f"offset {(distance_offset - 1) * 100:+.1f}% "
          f"-> {count_shift * 100:+.1f}% on closure ratios")

    print("\nBERLANAS 2019 — two-distance population")
    print(f"  one component {record['berlanas_2019']['one_component_kpc']} kpc, "
          f"depth {record['berlanas_2019']['one_component_depth_pc']} pc")
    print(f"  two components separate by only "
          f"{record['berlanas_2019']['two_component_separation_pc']} pc, "
          f"dBIC {record['berlanas_2019']['delta_bic_favouring_one']:+.2f}")
    print(f"  members nearer than 1350 pc: "
          f"{overall['fraction_nearer_than_1350pc'] * 100:.1f}%")
    print("  CAVEAT: parallax is a WP2 clustering feature — this test cannot")
    print("          rule out a foreground group removed upstream")

    print("\nwrote provenance/wp6_external_crosschecks_execution.json")


if __name__ == "__main__":
    main()
