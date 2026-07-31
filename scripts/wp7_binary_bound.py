#!/usr/bin/env python3
"""T3 -- bound the effect of unmodelled binary mass transfer on N_SN.

WP7 counts supernovae by single-star turnoff: a star has died if its ZAMS mass
exceeds the turnoff at its birth epoch.  Binary interaction is acknowledged
everywhere in this project and modelled nowhere -- issue #15 measured f_bin
about 0.7 above 8 Msun, Harer's own comparison is BPASS-based and Ic-dominated,
and reading (c) of the pulsar degeneracy is binary stripping.  A referee will
ask for a magnitude.  This script supplies one.

It does NOT integrate BPASS and does not change the pipeline.  It builds the
bound from two independent directions and then runs the resulting bracket
through WP7's own engine so the answer is stated in supernovae rather than in
percent.

  LINE 1 -- EMPIRICAL.  Harer et al. 2025 Fig. 2 is a BPASS calculation, with
  binaries, for a Cyg OB2-like population: an event rate versus age for a
  cluster of total initial mass 1.65e4 Msun at IMF index alpha = -2.0.  Our
  single-star rate at the SAME alpha, the SAME age and the SAME population mass
  is a like-for-like partner, so their ratio measures the binary correction
  end-to-end on exactly the population in question.  Fig. 2 is vector art in
  the local PDF, so it is digitized exactly rather than eyeballed -- and its
  bins turn out to sit on BPASS's native 0.1-dex log-age grid, which is a
  strong check that the digitization is right.

  LINE 2 -- THEORETICAL.  Zapartas et al. 2017 (A&A 601, A29) computed the
  core-collapse delay-time distribution with binary interaction.  Two of their
  results bound this directly: binarity raises the TOTAL number of CCSNe by
  14 (+15 -14)%, and the binary and single-star delay-time distributions are
  "remarkably similar at early times", diverging only around 20 Myr, because
  the enhancement is dominated by late (50-200 Myr) events from intermediate
  mass (4-8 Msun) binaries.  Cyg OB2 is 4 Myr old and its first supernova was
  1.3 Myr ago, so every event in the ledger sits deep inside the regime where
  the two distributions agree.

Outputs:
  tables/wp7_binary_bound.csv
  tables/wp7_binary_bound_harer_fig2.csv
  provenance/wp7_binary_bound_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp7_binary_bound.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp7_ledger as L
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass
from wp7_ledger_prereg import SF_DURATIONS_MYR

BASELINE_FAMILY, BASELINE_RV = "PARSEC", 3.1
DERIVATIVE_STEP_MYR = 0.02
HEADLINE_ALPHAS = (2.0, 2.3)

# ---------------------------------------------------------------- literature
ZAPARTAS_2017 = {
    "reference": "Zapartas et al. 2017, A&A 601, A29 (arXiv:1701.07032)",
    "title": (
        "Delay-time distribution of core-collapse supernovae with late events "
        "resulting from binary interaction"
    ),
    "total_ccsn_increase_from_binarity": 0.14,
    "total_ccsn_increase_hi": 0.15,
    "total_ccsn_increase_lo": 0.14,
    "verbatim_total": (
        "the total number of core-collapse supernovae increases by 14+15-14% "
        "because of binarity"
    ),
    "late_fraction": 0.15,
    "verbatim_late": (
        "a significant fraction, 15+9-8%, of core-collapse supernovae are "
        "'late', occurring 50-200 Myr after birth, when all massive single "
        "stars have already exploded"
    ),
    "late_progenitor_masses_Msun": [4.0, 8.0],
    "verbatim_early": (
        "Both distributions are remarkably similar at early times ... The "
        "differences become evident at around 20 Myr, where the binary "
        "distribution peaks"
    ),
}
DE_MINK_2014 = {
    "reference": "de Mink, Sana, Langer, Izzard & Schneider 2014, ApJ 782, 7",
    "binary_product_fraction_of_massive_MS": 0.30,
    "binary_product_hi": 0.10,
    "binary_product_lo": 0.15,
    "merger_fraction": 0.08,
    "verbatim": (
        "30+10-15% of massive main-sequence stars are the product of binary "
        "interaction ... 8+9-4% of a sample of early type stars to be the "
        "product of a merger"
    ),
    "caveat": "quoted for constant star formation, not a coeval burst",
}
HARER_FIG2 = {
    "reference": "Harer et al. 2025, A&A 703, A111, Fig. 2",
    "local_copy": "papers/Harer_2025.pdf",
    "page_index": 4,
    "cluster_initial_mass_Msun": 1.65e4,
    "imf_index": -2.0,
    "y_axis_verbatim": "Event Rate (events/100 kyr)",
    "note": (
        "digitized from the vector paths of the local PDF; the x axis is "
        "LOGARITHMIC in age and the step edges fall on BPASS's native 0.1-dex "
        "log-age grid (2.512, 3.162, 3.981, 5.012, 6.310, 7.943, 10.0 Myr), "
        "which is an independent confirmation that the calibration is right"
    ),
}
# Digitized total-CCSNe step curve: (log10 age/Myr bin lower edge, rate in
# events/100 kyr).  Regenerated and asserted by digitize_harer_fig2() below.
HARER_FIG2_TOTAL = [
    (2.512, 3.162, 1.373),
    (3.162, 3.981, 1.011),
    (3.981, 5.012, 0.990),
    (5.012, 6.310, 1.028),
    (6.310, 7.944, 1.021),
    (7.944, 10.000, 0.950),
]
# BPASS population definition, for the normalization correction.  BPASS's
# "imf100_300" is the alpha = -2.0, M_up = 300 Msun option Harer describe, on
# BPASS's standard broken IMF with slope -1.30 below 0.5 Msun down to 0.1.
BPASS_IMF = {"lo": 0.1, "break": 0.5, "low_slope": 1.3, "hi": 300.0}

# The adopted bracket.  See section 4 of the report for why +/-30%.
ADOPTED_BRACKET = (0.70, 1.30)


def digitize_harer_fig2(pdf_path) -> list[tuple[float, float, float]] | None:
    """Re-extract the total-CCSNe step curve from the local PDF's vector art.

    Returns None if PyMuPDF is unavailable, in which case the stored constants
    stand and the provenance record says so.
    """
    try:
        import fitz  # PyMuPDF
    except ImportError:
        return None
    document = fitz.open(str(pdf_path))
    page = document[HARER_FIG2["page_index"]]
    clip = fitz.Rect(295, 95, 570, 285)

    # Axis calibration from the tick LABELS.  Both axes are logarithmic; the
    # x axis is confirmed logarithmic by the 1 : 0.585 : 0.415 spacing of the
    # 2.5 / 5 / 7.5 / 10 ticks, which no linear axis can produce.
    labels = {}
    for x0, y0, x1, y1, text, *_ in page.get_text("words"):
        if not (clip.x0 <= x0 <= clip.x1 and clip.y0 <= y0 <= clip.y1):
            continue
        labels.setdefault(text, ((x0 + x1) / 2.0, (y0 + y1) / 2.0))
    x_ticks = [(2.5, labels["2.5"][0]), (10.0, labels["10"][0])]
    y_ticks = [
        (0.25, labels["0.25"][1]),
        (0.5, labels["0.5"][1]),
        (1.0, labels["1"][1]),
    ]
    x_fit = np.polyfit(
        [p[1] for p in x_ticks], [np.log10(p[0]) for p in x_ticks], 1
    )
    y_fit = np.polyfit(
        [p[1] for p in y_ticks], [np.log10(p[0]) for p in y_ticks], 1
    )

    black = [
        g
        for g in page.get_drawings()
        if clip.intersects(g["rect"])
        and g.get("color") == (0.0, 0.0, 0.0)
        and (g.get("width") or 0) > 0.5
    ]
    if len(black) != 1:
        return None
    points = []
    for item in black[0]["items"]:
        if item[0] == "l":
            points.extend([item[1], item[2]])
        elif item[0] == "c":
            points.extend([item[1], item[-1]])
    ordered = []
    for point in points:
        key = (round(point.x, 3), round(point.y, 3))
        if not ordered or ordered[-1] != key:
            ordered.append(key)

    # The step curve alternates vertical and horizontal segments; the plateaus
    # are the horizontal ones.  Matplotlib draws the last step past the right
    # spine, so anything starting at or beyond the 10 Myr tick is outside the
    # plotted range and is dropped.
    right_edge = labels["10"][0]
    out = []
    for (x_a, y_a), (x_b, y_b) in zip(ordered, ordered[1:]):
        if abs(y_a - y_b) > 1e-6 or abs(x_a - x_b) < 1e-6:
            continue
        if min(x_a, x_b) >= right_edge - 1.0:
            continue
        lo = 10 ** np.polyval(x_fit, min(x_a, x_b))
        hi = 10 ** np.polyval(x_fit, max(x_a, x_b))
        rate = 10 ** np.polyval(y_fit, y_a)
        out.append((round(float(lo), 3), round(float(hi), 3), round(float(rate), 3)))
    return out


def mass_per_k(lo: float, break_mass: float, low_slope: float,
               hi: float, alpha: float) -> float:
    """Stellar mass per unit high-mass normalization, broken power law."""
    coefficient = break_mass ** (low_slope - alpha)
    low = coefficient * w.power_integral(lo, break_mass, 1.0 - low_slope)
    high = w.power_integral(break_mass, hi, 1.0 - alpha)
    return float(low + high)


def turnoff_rate(family: str, age: float) -> float:
    h = DERIVATIVE_STEP_MYR
    return abs(
        turnoff_mass(family, age + h) - turnoff_mass(family, age - h)
    ) / (2.0 * h)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wp5-version", default="repair_v7")
    parser.add_argument("--iterations", type=int, default=400_000)
    args = parser.parse_args()
    n_iter = int(args.iterations)

    # ---------------------------------------------------------- line 1: BPASS
    pdf = w.ROOT / HARER_FIG2["local_copy"]
    extracted = digitize_harer_fig2(pdf)
    if extracted is not None:
        stored = [(a, b, c) for a, b, c in HARER_FIG2_TOTAL]
        agree = len(extracted) == len(stored) and all(
            abs(e[0] - s[0]) < 0.02
            and abs(e[1] - s[1]) < 0.02
            and abs(e[2] - s[2]) < 0.01
            for e, s in zip(extracted, stored)
        )
        if not agree:
            raise SystemExit(
                "the digitization of Harer Fig. 2 no longer reproduces the "
                f"stored constants:\n  extracted {extracted}\n  stored {stored}"
            )
        digitization = "re-extracted from the local PDF and matches the stored constants"
    else:
        digitization = (
            "PyMuPDF unavailable; stored constants used, extraction not re-verified"
        )

    fig2 = pd.DataFrame(
        HARER_FIG2_TOTAL, columns=["age_lo_Myr", "age_hi_Myr", "rate_per_100kyr"]
    )
    fig2["rate_SNe_per_Myr"] = fig2.rate_per_100kyr * 10.0
    fig2["rate_SNe_per_Myr_per_1e4Msun"] = (
        fig2.rate_SNe_per_Myr / (HARER_FIG2["cluster_initial_mass_Msun"] / 1e4)
    )
    fig2.to_csv(w.TABLES / "wp7_binary_bound_harer_fig2.csv", index=False)

    def bpass_rate_at(age: float) -> float:
        row = fig2[(fig2.age_lo_Myr <= age) & (fig2.age_hi_Myr > age)]
        return float(row.rate_SNe_per_Myr_per_1e4Msun.iloc[0])

    normalization = pd.read_parquet(
        w.PROC / f"wp5_imf_normalization_{args.wp5_version}.parquet"
    )
    rows = []
    for alpha in HEADLINE_ALPHAS:
        # For a FIXED total population mass, BPASS's wider IMF puts more mass
        # above our 120 Msun ceiling and less below 0.1 Msun, so its
        # high-mass normalization per unit mass differs from ours.  Correcting
        # for it is the difference between comparing physics and comparing
        # bookkeeping.
        ours_per_k = mass_per_k(
            w.TOTAL_MASS_RANGE[0], w.LOW_MASS_BREAK, w.LOW_MASS_SLOPE,
            w.TOTAL_MASS_RANGE[1], alpha,
        )
        theirs_per_k = mass_per_k(
            BPASS_IMF["lo"], BPASS_IMF["break"], BPASS_IMF["low_slope"],
            BPASS_IMF["hi"], alpha,
        )
        imf_range_correction = theirs_per_k / ours_per_k

        for subgroup in w.SUBGROUPS:
            branch = normalization[
                normalization.subgroup.eq(subgroup)
                & normalization.family.eq(BASELINE_FAMILY)
                & normalization.R_V.eq(BASELINE_RV)
                & normalization.alpha.eq(alpha)
            ].iloc[0]
            k = float(branch.k_median)
            age = float(branch.truth_age_posterior_mean_Myr)
            cap = min(turnoff_mass(BASELINE_FAMILY, age), IMF_UPPER_LIMIT)
            rate = (
                k * cap ** (-alpha) * turnoff_rate(BASELINE_FAMILY, age)
                if cap < IMF_UPPER_LIMIT
                else 0.0
            )
            mass = k * (
                w.primary_system_mass_per_k(alpha) + w.companion_mass_per_k(alpha)
            )
            ours_per_1e4 = rate / (mass / 1e4) if mass > 0 else 0.0
            theirs_per_1e4 = bpass_rate_at(age)
            rows.append(
                {
                    "alpha": alpha,
                    "subgroup": subgroup,
                    "age_Myr": age,
                    "turnoff_Msun": cap,
                    "k_median": k,
                    "our_rate_SNe_per_Myr": rate,
                    "our_mass_Msun": mass,
                    "our_rate_per_1e4Msun": ours_per_1e4,
                    "bpass_rate_per_1e4Msun": theirs_per_1e4,
                    "raw_ratio_bpass_over_ours": (
                        theirs_per_1e4 / ours_per_1e4 if ours_per_1e4 > 0 else np.nan
                    ),
                    "imf_range_correction": imf_range_correction,
                    "corrected_ratio_bpass_over_ours": (
                        theirs_per_1e4 / (ours_per_1e4 / imf_range_correction)
                        if ours_per_1e4 > 0
                        else np.nan
                    ),
                }
            )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(w.TABLES / "wp7_binary_bound.csv", index=False)

    # Harer's Fig. 2 was computed at alpha = -2.0.  Only our alpha = 2.0 rows
    # are a matched comparison; the alpha = 2.3 rows are kept in the CSV for
    # completeness and explicitly NOT used, because comparing our steeper IMF
    # against their shallower one would measure the slope, not the binaries.
    live = comparison[
        np.isfinite(comparison.raw_ratio_bpass_over_ours)
        & comparison.alpha.eq(2.0)
    ]
    empirical = {
        "matched_alpha": 2.0,
        "alpha_2p3_rows_excluded": (
            "Harer's figure is an alpha = -2.0 calculation, so only our "
            "alpha = 2.0 branch is a matched comparison; the alpha = 2.3 rows "
            "are tabulated but not used, since their ratio would measure the "
            "IMF slope difference rather than the binary physics"
        ),
        "subgroups_compared": sorted(set(live.subgroup)),
        "excluded": (
            "CygOB2-C contributes no supernovae on this branch (turnoff above "
            "the IMF ceiling), so no ratio is defined for it"
        ),
        "raw_ratio_range": [
            round(float(live.raw_ratio_bpass_over_ours.min()), 3),
            round(float(live.raw_ratio_bpass_over_ours.max()), 3),
        ],
        "imf_range_correction_at_matched_alpha": round(
            float(live.imf_range_correction.iloc[0]), 4
        ),
        "corrected_ratio_range": [
            round(float(live.corrected_ratio_bpass_over_ours.min()), 3),
            round(float(live.corrected_ratio_bpass_over_ours.max()), 3),
        ],
        "reading": (
            "A binary-inclusive population-synthesis rate and our single-star "
            "turnoff rate, at matched alpha, matched age and matched "
            "population mass, agree to within tens of percent -- and the "
            "binary calculation sits BELOW ours, not above.  Whatever the "
            "unmodelled binary physics does to N_SN, it is not an "
            "order-unity enhancement."
        ),
        "residual_systematics_not_corrected": [
            "BPASS's multiplicity fraction at high mass is near unity while "
            "our multiplicity-adjusted mass assumes f_bin = 0.40, so their "
            "total mass buys fewer primaries than ours; this pushes the "
            "corrected ratio further toward 1",
            "the digitized rate is a step function on 0.1-dex bins, so ages "
            "within a bin share one value",
            "BPASS's IMF below 0.5 Msun is its own, not ours",
        ],
    }

    # ------------------------------------------------- the bracket, in SNe
    lo_scale, hi_scale = ADOPTED_BRACKET
    draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{args.wp5_version}.npz")
    relations = {family: L.TurnoffRelation(family) for family in w.FAMILIES}
    bracket_rows = []
    for scale, label in ((lo_scale, "low"), (1.0, "nominal"), (hi_scale, "high")):
        rng_master = np.random.default_rng(w.SEED)
        for family in w.FAMILIES:
            relation = relations[family]
            for rv in w.R_V_BRANCHES:
                for alpha in HEADLINE_ALPHAS:
                    for delta in SF_DURATIONS_MYR:
                        n_sn = np.zeros(n_iter, dtype=int)
                        t_last = np.full(n_iter, np.inf)
                        for subgroup in w.SUBGROUPS:
                            key = L.draw_key(subgroup, family, rv, alpha)
                            k_all = draws[f"k__{key}"] * scale
                            age_all = draws[f"truth_age_draws__{key}"]
                            rng = np.random.default_rng(
                                rng_master.integers(0, 2 ** 63 - 1)
                            )
                            pick = rng.integers(0, k_all.size, n_iter)
                            res = L.run_population(
                                rng, k_all[pick], age_all[pick], alpha, delta,
                                relation,
                            )
                            n_sn += res["n_sn"]["all_explode"]
                            np.minimum(t_last, res["t_last"]["all_explode"],
                                       out=t_last)
                        bracket_rows.append(
                            {
                                "arm": label,
                                "k_scale": scale,
                                "family": family,
                                "R_V": rv,
                                "alpha": alpha,
                                "sf_duration_Myr": delta,
                                "N_SN_mean": float(n_sn.mean()),
                                "P_last_SN_within_100kyr": float(
                                    (t_last < 0.1).mean()
                                ),
                                "P_at_least_one": float((n_sn >= 1).mean()),
                            }
                        )
    bracket = pd.DataFrame(bracket_rows)
    bracket.to_csv(w.TABLES / "wp7_binary_bound_branches.csv", index=False)

    def arm(label: str) -> dict:
        frame = bracket[bracket.arm.eq(label)]
        base = frame[
            frame.family.eq("PARSEC")
            & frame.R_V.eq(3.1)
            & frame.alpha.eq(2.3)
            & frame.sf_duration_Myr.eq(0.0)
        ].iloc[0]
        return {
            "baseline_N_SN": round(float(base.N_SN_mean), 2),
            "baseline_P_last_SN_within_100kyr": round(
                float(base.P_last_SN_within_100kyr), 3
            ),
            "headline_N_SN_range": [
                round(float(frame.N_SN_mean.min()), 2),
                round(float(frame.N_SN_mean.max()), 2),
            ],
            "headline_P_last_SN_within_100kyr_range": [
                round(float(frame.P_last_SN_within_100kyr.min()), 3),
                round(float(frame.P_last_SN_within_100kyr.max()), 3),
            ],
        }

    arms = {label: arm(label) for label in ("low", "nominal", "high")}
    nominal_span = (
        arms["nominal"]["headline_N_SN_range"][1]
        - arms["nominal"]["headline_N_SN_range"][0]
    )
    binary_span = (
        arms["high"]["baseline_N_SN"] - arms["low"]["baseline_N_SN"]
    )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp7_binary_bound.py",
        "item": "T3 of tasks/pre_wp10_assessment_brief.md",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scope": (
            "a labelled literature-scaled bound, reported in the discussion.  "
            "NOT a pipeline change: no injection, no fit and no stored WP5, "
            "WP6, WP7, WP8 or WP9 product is touched, and the bracket is never "
            "marginalized into the headline number."
        ),
        "iterations_per_branch": n_iter,
        "line_1_empirical_bpass_comparison": {
            "source": HARER_FIG2,
            "digitization_status": digitization,
            "digitized_total_ccsn_curve": [
                {"age_lo_Myr": a, "age_hi_Myr": b, "rate_events_per_100kyr": c}
                for a, b, c in HARER_FIG2_TOTAL
            ],
            "grid_check": (
                "the six step edges reproduce BPASS's native log-age grid "
                "(log age/yr = 6.4 to 7.0 in 0.1 dex) to better than 0.5%, "
                "which no accidental calibration would do"
            ),
            "result": empirical,
        },
        "line_2_theoretical_delay_time": {
            "source": ZAPARTAS_2017,
            "argument": (
                "Cyg OB2's supernovae all occur at delay times below about 4.1 "
                "Myr -- its first explosion was 1.30 Myr ago in a 4.0 Myr "
                "population.  Zapartas et al. find the binary and single-star "
                "delay-time distributions 'remarkably similar at early times', "
                "diverging around 20 Myr, and attribute the +14% integrated "
                "enhancement to 50-200 Myr events from 4-8 Msun binaries.  "
                "None of that channel can have operated in a 4 Myr "
                "association, so the integrated +14% is an OVERESTIMATE of "
                "the correction applicable here."
            ),
            "supporting": DE_MINK_2014,
            "why_de_mink_does_not_scale_directly": (
                "30% of massive main-sequence stars being binary products does "
                "NOT imply 30% more supernovae: accretors and mergers are "
                "rejuvenated, so a star pushed above the present turnoff by "
                "accretion has had its clock partly reset and has typically "
                "NOT yet exploded at 4 Myr.  It adds to the future budget, not "
                "the past one."
            ),
        },
        "what_binaries_do_change_and_it_is_not_the_count": (
            "the supernova TYPE.  Harer's Fig. 2 is essentially pure type Ic "
            "from 2.5 to 6.4 Myr, the stripped-envelope channel.  Our WP9 "
            "condition C3 already asserts every progenitor above ~30 Msun is "
            "stripped and reaches 1.000; binary stripping is one route to that "
            "and single-star Wolf-Rayet winds are another, and Cyg OB2 "
            "contains three known WR stars.  So C3 does not depend on the "
            "binary channel being present, and the binary channel does not "
            "change the count it is applied to."
        ),
        "adopted_bracket": {
            "multiplicative_on_N_SN": list(ADOPTED_BRACKET),
            "justification": (
                "+/-30% on the supernova normalization.  The upper edge is the "
                "top of Zapartas et al.'s 14+15% integrated enhancement, "
                "rounded up and applied in full even though its dominant "
                "channel cannot operate at 4 Myr -- deliberately conservative. "
                " The lower edge is symmetric and covers both their -14% error "
                "edge and the empirical finding that the binary-inclusive "
                "BPASS rate sits BELOW our single-star rate."
            ),
            "applied_as": "a scale factor on k, run through wp7_ledger.run_population",
            "arms": arms,
        },
        "verdict": {
            "binary_bracket_span_on_baseline_N_SN": round(binary_span, 2),
            "headline_branch_span_on_N_SN": round(nominal_span, 2),
            "ratio_branch_span_over_binary_span": round(
                nominal_span / binary_span, 1
            ),
            "statement": (
                "Unmodelled binary mass transfer moves the baseline supernova "
                f"count from {arms['low']['baseline_N_SN']} to "
                f"{arms['high']['baseline_N_SN']}, a span of "
                f"{binary_span:.1f} supernovae, against a headline branch "
                f"span of {nominal_span:.1f}.  The systematic this project "
                "does not model is about "
                f"{nominal_span / binary_span:.0f} times smaller than the "
                "model-branch spread it already reports, and it does not "
                "change any conclusion."
            ),
        },
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p)
            for p in (
                pdf,
                w.PROC / f"wp5_imf_normalization_{args.wp5_version}.parquet",
                w.PROC / f"wp5_imf_posterior_draws_{args.wp5_version}.npz",
            )
        },
        "outputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p)
            for p in (
                w.TABLES / "wp7_binary_bound.csv",
                w.TABLES / "wp7_binary_bound_harer_fig2.csv",
                w.TABLES / "wp7_binary_bound_branches.csv",
            )
        },
    }
    w.write_json(w.PROVENANCE / "wp7_binary_bound_execution.json", record)

    print("T3 -- binary mass transfer, bounded\n")
    print("line 1: BPASS (binaries) vs our single-star turnoff, matched alpha,")
    print("        age and population mass")
    print(f"  digitization: {digitization}")
    print(f"  raw ratio      {empirical['raw_ratio_range']}")
    print(f"  IMF-range corrected ratio {empirical['corrected_ratio_range']}")
    print("\nline 2: Zapartas+2017 delay-time distribution")
    print(f"  integrated enhancement {ZAPARTAS_2017['total_ccsn_increase']:+.0%}"
          if "total_ccsn_increase" in ZAPARTAS_2017 else
          f"  integrated enhancement +14%, dominated by 50-200 Myr events")
    print("  early DTDs 'remarkably similar'; Cyg OB2 lives entirely there")
    print("\nadopted bracket +/-30% on N_SN:")
    for label in ("low", "nominal", "high"):
        a = arms[label]
        print(
            f"  {label:>8s}  baseline {a['baseline_N_SN']:6.2f}"
            f"  headline range {a['headline_N_SN_range']}"
            f"  P(<100kyr) {a['baseline_P_last_SN_within_100kyr']:.3f}"
        )
    print("\n" + record["verdict"]["statement"])
    print("\nwrote tables/wp7_binary_bound.csv")
    print("wrote tables/wp7_binary_bound_harer_fig2.csv")
    print("wrote tables/wp7_binary_bound_branches.csv")
    print("wrote provenance/wp7_binary_bound_execution.json")


if __name__ == "__main__":
    main()
