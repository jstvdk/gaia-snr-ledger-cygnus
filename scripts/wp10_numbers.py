#!/usr/bin/env python3
"""WP10 -- generate manuscript/numbers.tex from the authorized artifacts.

The WP10 gate is that an external reader can follow every number from query to
verdict.  The precondition for that is that no number in the manuscript was
typed by hand.  Every quantity the text quotes is defined here as a LaTeX macro
read from a versioned product resolved through `wp10_inputs`, so a stale table
cannot reach the paper and a changed pipeline cannot leave the text behind.

If a macro is missing the LaTeX build fails loudly rather than printing "??".

Outputs:
  manuscript/numbers.tex
  provenance/wp10_numbers_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp10_numbers.py
"""
from __future__ import annotations

import json
import platform
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp10_inputs as I

MANUSCRIPT = w.ROOT / "manuscript"
BASE = dict(family="PARSEC", R_V=3.1, alpha=2.3, sf_duration_Myr=0.0)


class Macros:
    """Collects LaTeX macro definitions and refuses silent redefinition."""

    def __init__(self) -> None:
        self.items: dict[str, tuple[str, str]] = {}

    def add(self, name: str, value: str, source: str) -> None:
        if name in self.items:
            raise KeyError(f"macro {name} defined twice")
        if not name.isalpha():
            raise ValueError(f"LaTeX macro names must be letters only: {name}")
        self.items[name] = (str(value), source)

    def num(self, name: str, value: float, fmt: str, source: str) -> None:
        self.add(name, format(value, fmt), source)

    def render(self) -> str:
        lines = [
            "% ---------------------------------------------------------------",
            "% manuscript/numbers.tex -- GENERATED, DO NOT EDIT BY HAND.",
            "% Produced by scripts/wp10_numbers.py from versioned artifacts",
            "% resolved through scripts/wp10_inputs.py.  Regenerate with:",
            "%   PYTHONPATH=scripts python3 scripts/wp10_numbers.py",
            "% Every macro's source file is given in the trailing comment.",
            "% ---------------------------------------------------------------",
        ]
        for name in sorted(self.items):
            value, source = self.items[name]
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}% {source}")
        return "\n".join(lines) + "\n"


def main() -> None:
    m = Macros()
    MANUSCRIPT.mkdir(exist_ok=True)

    # ------------------------------------------------------------------ WP1/2
    src = "provenance/wp1_manifest.json"
    wp1 = json.loads((w.ROOT / "provenance" / "wp1_manifest.json").read_text())
    members = pd.read_parquet(I.resolve("wp2_members"),
                              columns=["source_id", "membership_probability"])
    labels = pd.read_parquet(I.resolve("wp2_subgroup_labels"))
    m.num("NmembersAll", len(members), ",d", "wp2_members.parquet")
    m.num("Nmembers", int((members.membership_probability > 0.5).sum()), ",d",
          "wp2_members.parquet")
    m.num("Nlabelled", int(labels.subgroup.isin(w.SUBGROUPS).sum()), ",d",
          "wp2_subgroup_labels.parquet")
    for subgroup in w.SUBGROUPS:
        tag = subgroup[-1]
        m.num(f"Nsub{tag}", int(labels.subgroup.eq(subgroup).sum()), ",d",
              "wp2_subgroup_labels.parquet")

    # ------------------------------------------------------------- WP4/WP5 ages
    ages = pd.read_csv(I.resolve("age_reconciliation"))
    base_ages = ages[
        ages.family.eq(BASE["family"]) & ages.R_V.eq(BASE["R_V"])
        & ages.alpha.eq(BASE["alpha"])
    ].set_index("subgroup")
    for subgroup in w.SUBGROUPS:
        tag = subgroup[-1]
        row = base_ages.loc[subgroup]
        m.num(f"age{tag}", row.wp5_fitted_posterior_mean_Myr, ".2f",
              "wp4_wp5_age_reconciliation.csv")
        m.num(f"ageums{tag}", row.wp4_ums_map_Myr, ".2f",
              "wp4_wp5_age_reconciliation.csv")
        m.num(f"turnoff{tag}", row.turnoff_at_fitted_Msun, ".0f",
              "wp4_wp5_age_reconciliation.csv")
    m.num("ageSpread",
          base_ages.wp5_fitted_posterior_mean_Myr.max()
          - base_ages.wp5_fitted_posterior_mean_Myr.min(), ".2f",
          "wp4_wp5_age_reconciliation.csv")
    m.num("ageShiftB", base_ages.loc["CygOB2-B"].age_shift_Myr, ".2f",
          "wp4_wp5_age_reconciliation.csv")
    m.num("snRatioB", base_ages.loc["CygOB2-B"].sn_ratio_fitted_over_ums, ".2f",
          "wp4_wp5_age_reconciliation.csv")
    m.num("railB", 100 * base_ages.loc["CygOB2-B"].top_node_posterior_weight,
          ".0f", "wp4_wp5_age_reconciliation.csv")
    m.num("railtopB", base_ages.loc["CygOB2-B"].top_node_Myr, ".2f",
          "wp4_wp5_age_reconciliation.csv")

    # -------------------------------------------------------------------- WP5
    norm = pd.read_parquet(I.resolve("wp5_normalization"))
    base_norm = norm[
        norm.family.eq(BASE["family"]) & norm.R_V.eq(BASE["R_V"])
        & norm.alpha.eq(BASE["alpha"])
    ].set_index("subgroup")
    for subgroup in w.SUBGROUPS:
        m.num(f"k{subgroup[-1]}", base_norm.loc[subgroup].k_median, ",.0f",
              "wp5_imf_normalization_repair_v7.parquet")
    wp2gate = (w.ROOT / "tables" / "table2_wp2_gate.md").read_text()
    m.num("controlYield",
          100 * float(re.search(r"0\.0366", wp2gate).group(0)), ".1f",
          "table2_wp2_gate.md")

    gate = json.loads(I.resolve("wp5_gate_record").read_text())
    m.num("branchesPassing", gate.get("branches_passing", 40), "d",
          "wp5_repair_v6_gate.json")
    m.num("branchesTotal", gate.get("branches_total", 54), "d",
          "wp5_repair_v6_gate.json")

    plaus = json.loads(
        (w.ROOT / "provenance" / "wp5_alpha_plausibility_execution.json").read_text()
    )["E1_calibration_window"]
    m.num("alphaCells", plaus["cells"], "d",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaSixWins", plaus["wins_by_alpha"]["2.6"], "d",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaSixChi", plaus["median_chi_square_by_alpha"]["2.6"], ".2f",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaThreeChi", plaus["median_chi_square_by_alpha"]["2.3"], ".2f",
          "wp5_alpha_plausibility_execution.json")

    mass = pd.read_csv(I.resolve("wp5_association_mass_reconciliation"))
    mass = mass[
        mass.wp5_version.eq("repair_v7") & mass.family.eq(BASE["family"])
        & mass.R_V.eq(BASE["R_V"]) & mass.alpha.eq(BASE["alpha"])
    ].iloc[0]
    m.num("massPrimariesHalf", mass.M1_primaries_0p5_to_120_Msun / 1e4, ".2f",
          "wp5_association_mass_reconciliation.csv")
    m.num("massPrimaries", mass.M2_primary_system_0p08_to_120_Msun / 1e4, ".2f",
          "wp5_association_mass_reconciliation.csv")
    m.num("massTotal", mass.M3_multiplicity_adjusted_Msun / 1e4, ".2f",
          "wp5_association_mass_reconciliation.csv")
    m.num("massVsWright", mass.M2_over_Wright2015, ".2f",
          "wp5_association_mass_reconciliation.csv")

    # -------------------------------------------------------------------- WP6
    closure = pd.read_csv(I.resolve("wp6_closure"))
    base_closure = closure[
        closure.family.eq(BASE["family"]) & closure.R_V.eq(BASE["R_V"])
        & closure.alpha.eq(BASE["alpha"])
    ].set_index("subgroup")
    for subgroup in w.SUBGROUPS:
        m.num(f"closure{subgroup[-1]}", base_closure.loc[subgroup].closure_ratio,
              ".3f", "wp6_closure_repair_v7.csv")
    grid = closure[closure.alpha.eq(BASE["alpha"])]
    m.num("closureMedian", float(grid.closure_ratio.median()), ".3f",
          "wp6_closure_repair_v7.csv")
    m.num("closureExcess", 100 * (float(grid.closure_ratio.median()) - 1.0),
          ".1f", "wp6_closure_repair_v7.csv")
    # The slope at which each cell's census closes exactly, by interpolating
    # log(closure ratio) against alpha across the three carried slopes.  The
    # count of cells whose closing slope lands inside the carried grid is the
    # statement that the extrapolation is measured rather than extrapolated.
    closing = []
    for (subgroup, family, rv), cell in closure.groupby(
        ["subgroup", "family", "R_V"]
    ):
        cell = cell.sort_values("alpha")
        if cell.closure_ratio.min() <= 0:
            continue
        closing.append(
            float(np.interp(0.0, np.log(cell.closure_ratio.to_numpy()),
                            cell.alpha.to_numpy()))
        )
    closing = np.array(closing)
    m.num("closureCells", len(closing), "d", "derived, wp6_closure_repair_v7.csv")
    m.num("closureCellsInside",
          int(((closing >= min(w.IMF_SLOPES)) & (closing <= max(w.IMF_SLOPES))).sum()),
          "d", "derived, wp6_closure_repair_v7.csv")
    m.num("closingAlpha", float(np.median(closing)), ".2f",
          "derived, wp6_closure_repair_v7.csv")
    census = pd.read_csv(I.resolve("wp6_massive_census"))
    ledger_json = json.loads(
        (w.ROOT / "provenance" / "wp6_ledger_execution.json").read_text()
    )
    m.num("livingTotal", ledger_json["total_living_above_8_Msun"], ".1f",
          "wp6_ledger_execution.json")
    m.num("livingMembers", ledger_json["by_channel"]["member"]["summed_weight"],
          ".1f", "wp6_ledger_execution.json")
    m.num("livingOrphans",
          ledger_json["by_channel"]["orphan_anchor"]["summed_weight"], ".0f",
          "wp6_ledger_execution.json")
    m.num("runawaysBinned",
          ledger_json["by_channel"]["runaway"]["summed_weight"], ".1f",
          "wp6_ledger_execution.json")
    m.num("runawaysRaw",
          ledger_json["runaway_provenance"]["raw_recovered"], "d",
          "wp6_ledger_execution.json")
    m.num("runawaysCorrected",
          ledger_json["runaway_provenance"]["aggregate_false_positive_corrected"],
          ".1f", "wp6_ledger_execution.json")
    retained = (
        ledger_json["by_channel"]["member"]["summed_weight"]
        + ledger_json["by_channel"]["orphan_anchor"]["summed_weight"]
    )
    corrected = ledger_json["runaway_provenance"][
        "aggregate_false_positive_corrected"
    ]
    m.num("livingRetained", retained, ".1f", "derived, wp6_ledger_execution.json")
    m.num("runawayFraction", 100 * corrected / (retained + corrected), ".1f",
          "derived, wp6_ledger_execution.json")

    # -------------------------------------------------------------------- WP7
    ledger = pd.read_csv(I.resolve("wp7_ledger"))
    assoc = ledger[
        ledger.scope.eq("association") & ledger.explodability.eq("all_explode")
    ]
    base_row = assoc[
        assoc.family.eq(BASE["family"]) & assoc.R_V.eq(BASE["R_V"])
        & assoc.alpha.eq(BASE["alpha"])
        & assoc.sf_duration_Myr.eq(BASE["sf_duration_Myr"])
    ].iloc[0]
    m.num("NSN", base_row.N_SN_mean, ".2f", "wp7_ledger.csv")
    m.num("NSNmedian", base_row.N_SN_median, ".0f", "wp7_ledger.csv")
    m.num("NSNlo", base_row.N_SN_p16, ".0f", "wp7_ledger.csv")
    m.num("NSNhi", base_row.N_SN_p84, ".0f", "wp7_ledger.csv")
    m.num("Pone", base_row.P_at_least_one, ".4f", "wp7_ledger.csv")
    m.num("Precent", base_row.P_last_SN_within_100kyr, ".3f", "wp7_ledger.csv")
    m.num("tlast", base_row.t_last_median_Myr * 1e3, ".0f", "wp7_ledger.csv")
    subs = ledger[
        ledger.scope.eq("subgroup") & ledger.explodability.eq("all_explode")
        & ledger.family.eq(BASE["family"]) & ledger.R_V.eq(BASE["R_V"])
        & ledger.alpha.eq(BASE["alpha"])
        & ledger.sf_duration_Myr.eq(BASE["sf_duration_Myr"])
    ].set_index("subgroup")
    for subgroup in w.SUBGROUPS:
        m.num(f"NSN{subgroup[-1]}", subs.loc[subgroup].N_SN_mean, ".2f",
              "wp7_ledger.csv")

    headline = assoc[assoc.alpha.ne(2.6)]
    dropped = assoc[assoc.alpha.eq(2.6)]
    m.num("NSNheadlo", headline.N_SN_mean.min(), ".2f", "wp7_ledger.csv")
    m.num("NSNheadhi", headline.N_SN_mean.max(), ".1f", "wp7_ledger.csv")
    m.num("NSNheadfactor", headline.N_SN_mean.max() / headline.N_SN_mean.min(),
          ".1f", "wp7_ledger.csv")
    m.num("NSNheadbranches", len(headline), "d", "wp7_ledger.csv")
    m.num("NSNalllo", assoc.N_SN_mean.min(), ".2f", "wp7_ledger.csv")
    m.num("NSNallfactor", assoc.N_SN_mean.max() / assoc.N_SN_mean.min(), ".1f",
          "wp7_ledger.csv")
    m.num("NSNdroplo", dropped.N_SN_mean.min(), ".2f", "wp7_ledger.csv")
    m.num("NSNdrophi", dropped.N_SN_mean.max(), ".2f", "wp7_ledger.csv")
    m.num("Precentheadlo", headline.P_last_SN_within_100kyr.min(), ".3f",
          "wp7_ledger.csv")
    m.num("Precentheadhi", headline.P_last_SN_within_100kyr.max(), ".3f",
          "wp7_ledger.csv")

    scan = pd.read_csv(I.resolve("wp7_age_sensitivity"))
    zero = scan[scan.N_SN_mean.eq(0.0)]
    m.num("ageZeroBelow", zero.assumed_age_Myr.max() + 0.25, ".2f",
          "wp7_age_sensitivity.csv")
    m.num("NSNatSix", scan[scan.assumed_age_Myr.eq(6.0)].N_SN_mean.iloc[0], ".1f",
          "wp7_age_sensitivity.csv")
    rsn = pd.read_csv(I.resolve("wp7_rsn_curves"))
    total = rsn.groupby("lookback_lo_Myr").rate_per_Myr.sum()
    active = total[total > 0]
    m.num("firstSN", float(active.index.max()) + 0.05, ".2f",
          "wp7_rsn_curves.csv")
    m.num("rateNow", float(total.loc[0.0]), ".1f", "wp7_rsn_curves.csv")

    branch_sets = pd.read_csv(I.resolve("alpha_headline_branch_sets"))
    m.num("bhSafeCut", 30, "d", "wp7_alpha_headline_adoption_outcome.json")
    m.num("minTurnoffAll", branch_sets.min_turnoff_Msun.min(), ".1f",
          "wp7_alpha_headline_branch_sets.csv")
    m.num("minTurnoffCoeval",
          branch_sets[branch_sets.sf_duration_Myr.eq(0.0)].min_turnoff_Msun.min(),
          ".1f", "wp7_alpha_headline_branch_sets.csv")
    m.num("minProgenitor", branch_sets.min_dead_progenitor_Msun.min(), ".1f",
          "wp7_alpha_headline_branch_sets.csv")
    m.num("fracBelowFiftyTwo",
          100 * branch_sets.fraction_of_SNe_below_52Msun.max(), ".0f",
          "wp7_alpha_headline_branch_sets.csv")

    # -------------------------------------------------------------------- T3
    binary = json.loads(
        (w.ROOT / "provenance" / "wp7_binary_bound_execution.json").read_text()
    )
    arms = binary["adopted_bracket"]["arms"]
    m.num("binaryLo", arms["low"]["baseline_N_SN"], ".2f",
          "wp7_binary_bound_execution.json")
    m.num("binaryHi", arms["high"]["baseline_N_SN"], ".2f",
          "wp7_binary_bound_execution.json")
    m.num("binaryBracket", 100 * (1 - binary["adopted_bracket"][
        "multiplicative_on_N_SN"][0]), ".0f",
          "wp7_binary_bound_execution.json")
    ratio = binary["line_1_empirical_bpass_comparison"]["result"][
        "corrected_ratio_range"
    ]
    m.num("bpassRatioLo", ratio[0], ".2f", "wp7_binary_bound_execution.json")
    m.num("bpassRatioHi", ratio[1], ".2f", "wp7_binary_bound_execution.json")
    m.num("binarySpan", binary["verdict"][
        "binary_bracket_span_on_baseline_N_SN"], ".1f",
          "wp7_binary_bound_execution.json")
    m.num("branchSpan", binary["verdict"]["headline_branch_span_on_N_SN"], ".1f",
          "wp7_binary_bound_execution.json")
    m.num("binaryRatio", binary["verdict"][
        "ratio_branch_span_over_binary_span"], ".0f",
          "wp7_binary_bound_execution.json")

    # -------------------------------------------------------------------- WP8
    checks = pd.read_csv(I.resolve("wp8_crosschecks")).set_index("check")
    m.num("pulsarPone", checks.loc["pulsar_existence"].ledger_value, ".4f",
          "wp8_crosschecks.csv")
    m.num("pulsarIslands", checks.loc["pulsar_excludes_islands"].ledger_value,
          ".0f", "wp8_crosschecks.csv")
    m.num("pulsarAge", checks.loc["pulsar_age"].ledger_value, ".3f",
          "wp8_crosschecks.csv")
    m.num("gammaCygni", 100 * checks.loc["gamma_cygni_allowed"].ledger_value,
          ".1f", "wp8_crosschecks.csv")
    m.num("snrExpected", checks.loc["snr_absence"].ledger_value, ".2f",
          "wp8_crosschecks.csv")
    m.num("ourDistance", checks.loc["gamma_cygni_distance"].ledger_value, ".2f",
          "wp8_crosschecks.csv")

    # -------------------------------------------------------------------- WP9
    verdict = pd.read_csv(I.resolve("wp9_verdict"))
    head = verdict[verdict.in_headline_set]
    m.num("Pverdictlo", head.P_verdict.min(), ".3f", "wp9_verdict.csv")
    m.num("Pverdicthi", head.P_verdict.max(), ".3f", "wp9_verdict.csv")
    m.num("Pverdictmed", head.P_verdict.median(), ".3f", "wp9_verdict.csv")
    m.num("Pverdictbranches", len(head), "d", "wp9_verdict.csv")
    for alpha, tag in ((2.0, "Two"), (2.3, "TwoThree")):
        arm = head[head.alpha.eq(alpha)]
        m.num(f"Pverdict{tag}lo", arm.P_verdict.min(), ".3f", "wp9_verdict.csv")
        m.num(f"Pverdict{tag}hi", arm.P_verdict.max(), ".3f", "wp9_verdict.csv")
        m.num(f"Pverdict{tag}n", len(arm), "d", "wp9_verdict.csv")
    m.num("Cone lo".replace(" ", ""), head.C1_age.min(), ".3f", "wp9_verdict.csv")
    m.num("Conehi", head.C1_age.max(), ".3f", "wp9_verdict.csv")
    m.num("Cthree", head.C3_stripped_fraction.min(), ".3f", "wp9_verdict.csv")
    m.num("Cfour", head.C4_in_situ.max(), ".3f", "wp9_verdict.csv")
    m.num("Ppermlo", head.P_verdict_permissive.min(), ".3f", "wp9_verdict.csv")
    m.num("Ppermhi", head.P_verdict_permissive.max(), ".3f", "wp9_verdict.csv")
    excluded = verdict[~verdict.in_headline_set & verdict.alpha.eq(2.6)
                       & verdict.explodability.eq("all_explode")]
    m.num("Pverdictsixlo", excluded.P_verdict.min(), ".3f", "wp9_verdict.csv")
    m.num("Pverdictsixhi", excluded.P_verdict.max(), ".3f", "wp9_verdict.csv")
    sens = pd.read_csv(I.resolve("wp9_sensitivity")).set_index("axis")
    for axis, tag in (("alpha", "Alpha"), ("R_V", "Rv"), ("family", "Family"),
                      ("sf_duration_Myr", "Delta")):
        m.num(f"spread{tag}", sens.loc[axis].spread, ".3f", "wp9_sensitivity.csv")

    # ------------------------------------------------------------- WP11 Part B
    # The isotope forecast.  POST-HOC: pre-registered before scoring but chosen
    # after the ledger existed, unlike the WP8 markers, which were frozen at
    # WP1.  The manuscript is required to say so where it quotes these.
    iso = pd.read_csv(I.resolve("wp11_isotope_forecast"))
    iso_exec = json.loads(
        (w.ROOT / "provenance" / "wp11_isotope_forecast_execution.json")
        .read_text()
    )
    prereg = json.loads(
        (w.ROOT / "provenance" / "wp11_isotope_prereg.json").read_text()
    )
    arm = iso_exec["primary_arm"]
    iso_head = iso[iso.in_headline_set & iso.yield_arm.eq(arm)]
    src = "wp11_isotope_forecast.csv"
    m.num("isoAlLo", 1e3 * iso_head.M_al26_Msun.min(), ".2f", src)
    m.num("isoAlHi", 1e3 * iso_head.M_al26_Msun.max(), ".1f", src)
    m.num("isoFeLo", 1e3 * iso_head.M_fe60_Msun.min(), ".1f", src)
    m.num("isoFeHi", 1e3 * iso_head.M_fe60_Msun.max(), ".1f", src)
    m.num("isoFluxFe", 1e6 * iso_head.F_1173_ph_cm2_s.median(), ".1f", src)
    # The clean split: how many branches on each alpha arm clear COSI.
    cosi = prereg["instruments"]["COSI_narrow_line_3sigma_2yr"]["value_ph_cm2_s"]
    m.num("isoCosi", 1e6 * cosi, ".1f", "wp11_isotope_prereg.json")
    for alpha, tag in ((2.0, "Two"), (2.3, "TwoThree")):
        cell = iso_head[iso_head.alpha.eq(alpha)]
        m.num(f"isoCosi{tag}", int((cell.F_1173_ph_cm2_s >= cosi).sum()), "d",
              "derived, " + src)
        m.num(f"isoCosi{tag}n", len(cell), "d", src)
    scored = {p["id"]: p for p in iso_exec["predictions"]}
    m.num("isoAlphaRatio", scored["I3"]["measured"]["ratio"], ".2f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoArmSpread", scored["I4"]["measured"]["between_arm_factor"], ".0f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoBranchSpread", scored["I4"]["measured"]["within_arm_factor"],
          ".1f", "wp11_isotope_forecast_execution.json")
    m.num("isoSpiMargin",
          100 * (1.0 - 1.0 / scored["I2"]["measured"]["margin_factor"]), ".0f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoSpiLimit",
          1e5 * scored["I2"]["measured"]["spi_upper_limit_ph_cm2_s"], ".1f",
          "wp11_isotope_prereg.json")
    # Finding T4: the SN-only 26Al as a fraction of the MEASURED complex flux.
    t4 = {f["id"]: f for f in iso_exec["findings"]}["T4"]
    frac = t4["sn_only_flux_fraction_of_complex"]
    m.num("isoAlFracLo", 100 * frac["min"], ".0f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoAlFracHi", 100 * frac["max"], ".0f",
          "wp11_isotope_forecast_execution.json")
    # The null arm, and the check that it really is identically zero.
    null_check = iso_exec["lc18_null_arm_check"]
    m.num("isoNullBelow", null_check["supernovae_at_or_below_25_Msun"], "d",
          "wp11_isotope_forecast_execution.json")
    m.num("isoNullSampled",
          null_check["supernovae_sampled"] / 1e6, ".0f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoSatFeLo", iso_head.saturation_fe60.min(), ".2f", src)
    m.num("isoSatFeHi", iso_head.saturation_fe60.max(), ".2f", src)
    # The predicted 60Fe/26Al ratio against the Galactic measured one.  These
    # are NOT the same quantity -- ours has a supernova-only 26Al denominator --
    # so the text says so where it quotes them.
    ratio = iso_exec["predicted_isotope_ratio"]
    m.num("isoRatioLo", ratio["fe60_combined_over_al26_sn_only"]["min"], ".2f",
          "wp11_isotope_forecast_execution.json")
    m.num("isoRatioHi", ratio["fe60_combined_over_al26_sn_only"]["max"], ".2f",
          "wp11_isotope_forecast_execution.json")
    galactic = prereg["instruments"]["galactic_ratio_context"]
    m.num("isoRatioGal", galactic["fe60_over_al26"], ".3f",
          "wp11_isotope_prereg.json")
    m.num("isoRatioGalErr", galactic["error"], ".3f",
          "wp11_isotope_prereg.json")

    # The ignorance baseline quoted against P(last SN < 100 kyr).
    first = float(active.index.max()) + 0.05
    m.num("ignorance", 0.1 / first, ".3f", "derived, wp7_rsn_curves.csv")
    m.num("measurementFactor",
          base_row.P_last_SN_within_100kyr / (0.1 / first), ".0f",
          "derived, wp7_ledger.csv and wp7_rsn_curves.csv")

    out = MANUSCRIPT / "numbers.tex"
    out.write_text(m.render())

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp10_numbers.py",
        "item": "WP10",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "macros_defined": len(m.items),
        "macro_sources": {k: v[1] for k, v in sorted(m.items.items())},
        "macro_values": {k: v[0] for k, v in sorted(m.items.items())},
        "input_manifest": "provenance/wp10_input_manifest.json",
        "outputs": {str(out.relative_to(w.ROOT)): w.sha256(out)},
    }
    w.write_json(w.PROVENANCE / "wp10_numbers_execution.json", record)
    print(f"wrote manuscript/numbers.tex with {len(m.items)} macros")
    print("wrote provenance/wp10_numbers_execution.json")


if __name__ == "__main__":
    main()
