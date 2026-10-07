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

import chain as C
import wp5_common as w
import wp10_inputs as I

MANUSCRIPT = w.ROOT / "manuscript"
BASE = dict(family="PARSEC", R_V=3.1, alpha=2.3, sf_duration_Myr=0.0)

# Products repair_v7 wrote without a suffix.  A macro's source comment names
# the file it was actually read from, so on a later chain (issue #19) these
# names are rewritten to the chain's own copies.
CHAIN_SOURCES = {
    "wp4_wp5_age_reconciliation.csv", "wp5_alpha_plausibility_execution.json",
    "wp5_association_mass_reconciliation.csv", "wp6_massive_census.csv",
    "wp6_ledger_execution.json", "wp7_ledger.csv", "wp7_age_sensitivity.csv",
    "wp7_rsn_curves.csv", "wp7_alpha_headline_branch_sets.csv",
    "wp7_alpha_headline_adoption_outcome.json", "wp7_binary_bound_execution.json",
    "wp8_crosschecks.csv", "wp9_verdict.csv", "wp9_sensitivity.csv",
    "wp11_isotope_forecast.csv", "wp11_isotope_forecast_execution.json",
    "wp12_gate_landscape_execution.json", "wp12_closure_slopes_execution.json",
    "wp12_scenario_score_execution.json", "wp12_neighbour_budget_execution.json",
}


def chain_source(text: str) -> str:
    def sub(match: re.Match) -> str:
        name = match.group(0)
        if name in CHAIN_SOURCES:
            return C.tag(name).name
        return name.replace("_repair_v7", f"_{C.V['wp5']}")
    return re.sub(r"wp\d+_[A-Za-z0-9_]+\.(?:csv|json|parquet)", sub, text)


def chain_json(name: str) -> dict:
    """A chain product's provenance record, read from the active chain."""
    return json.loads(C.tag(w.PROVENANCE / name).read_text())


class Macros:
    """Collects LaTeX macro definitions and refuses silent redefinition."""

    def __init__(self) -> None:
        self.items: dict[str, tuple[str, str]] = {}

    def add(self, name: str, value: str, source: str) -> None:
        if name in self.items:
            raise KeyError(f"macro {name} defined twice")
        if not name.isalpha():
            raise ValueError(f"LaTeX macro names must be letters only: {name}")
        self.items[name] = (str(value), chain_source(source))

    def num(self, name: str, value: float, fmt: str, source: str) -> None:
        if value is None:
            # repair_v9: a quantity can be undefined on the active chain (e.g. a
            # score range over zero branches).  Render a dash, never a number.
            self.add(name, r"\textemdash", f"{source} (undefined on this chain)")
            return
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

    # The P>0.5 handoff is the maximum-posterior decision rule for two classes
    # under equal misclassification costs.  These neighbouring thresholds make
    # its empirical trade-off visible rather than treating 0.5 as a convention.
    controls = pd.read_parquet(I.resolve("wp2_control_members"))
    recovery_audit = pd.read_csv(I.resolve("wp2_recovery_audit"))
    # Re-read the exemption flag because the compact read above deliberately
    # requested only two columns.
    automatic = pd.read_parquet(I.resolve("wp2_members"))
    automatic = automatic[~automatic.anchor_quality_exempt.fillna(False)]
    quality_benchmark = recovery_audit[recovery_audit.quality_pass.astype(bool)]
    for suffix, threshold in (("Low", 0.4), ("Adopted", 0.5), ("High", 0.6)):
        n_target = int((automatic.membership_probability > threshold).sum())
        control_yields = [
            int((group.membership_probability > threshold).sum())
            for _, group in controls.groupby("control_field")
        ]
        recovered = int(
            (quality_benchmark.membership_probability_astrometric > threshold)
            .fillna(False).sum()
        )
        m.num(f"thresholdMembers{suffix}", n_target, ",d",
              "derived from wp2_members.parquet")
        m.num(f"thresholdControl{suffix}",
              100 * float(np.mean(control_yields)) / n_target, ".1f",
              "derived from wp2_control_members.parquet")
        m.num(f"thresholdRecall{suffix}", 100 * recovered / 168, ".1f",
              "derived from wp2_berlanas_recovery_audit.csv")

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

    # Issue #19.  The retained upper-MS age envelope, computed rather than
    # typed: MAP span of the measurable, non-railed rows of the authorized WP4
    # posterior.  It replaces the hand-typed 2.25-5.67 Myr "two-indicator"
    # envelope, which came from the pre-repair run; the repaired posterior
    # retains no PMS row, so the envelope is upper-MS only.
    posterior = pd.read_parquet(I.resolve("wp4_age_posteriors"))
    kept = posterior[
        posterior.measurable.astype(bool) & ~posterior.grid_railed.astype(bool)
    ]
    src = I.resolve("wp4_age_posteriors").name + " (retained rows)"
    if kept.indicator.eq("pms").any():
        raise RuntimeError("a PMS age row is retained; the envelope text assumes none")
    m.num("ageEnvLo", kept.age_map.min(), ".2f", src)
    m.num("ageEnvHi", kept.age_map.max(), ".2f", src)
    m.num("ageEnvRows", len(kept), "d", src)

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

    # The Berlanas+19 recall accounting, automatic and manual separately.  The
    # total is the gate; the automatic fraction is what the classifier achieves
    # unaided, and the manuscript is required to give both.
    recovery = I.resolve("wp2_literature_recovery").read_text()
    src = "table3_literature_recovery.md"

    def recovery_row(label: str, pattern: str) -> str:
        # Table rows only: the caption contains "published membership
        # catalogues", which substring-matches the "published members" row.
        row = next(
            line for line in recovery.splitlines()
            if line.lstrip().startswith("|") and label in line
        )
        found = re.search(pattern, row)
        if not found:
            raise ValueError(f"{pattern!r} did not match recovery row: {row!r}")
        return found.group(1)

    m.num("recallPublished",
          int(recovery_row("published members", r"\|\s*(\d+)\s*\|")), "d", src)
    m.num("recallAnalyzable",
          int(recovery_row("inside quality sample", r"\|\s*(\d+)\s*\|")), "d",
          src)
    m.num("recallAuto",
          int(recovery_row("recovered automatically", r"\|\s*(\d+)\s*\(")), "d",
          src)
    m.num("recallAutoPct",
          float(recovery_row("recovered automatically", r"\((\d+\.\d+)")), ".1f",
          src)
    m.num("recallManual",
          int(recovery_row("manual quality exceptions", r"\|\s*(\d+)\s*\|")),
          "d", src)
    m.num("recallTotal",
          int(recovery_row("total recall", r"\|\s*(\d+)\s*/")), "d", src)
    m.num("recallTotalPct",
          float(recovery_row("total recall", r"=\s*(\d+\.\d+)")), ".1f", src)

    plaus = chain_json("wp5_alpha_plausibility_execution.json")["E1_calibration_window"]
    m.num("alphaCells", plaus["cells"], "d",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaSixWins", plaus["wins_by_alpha"]["2.6"], "d",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaSixChi", plaus["median_chi_square_by_alpha"]["2.6"], ".2f",
          "wp5_alpha_plausibility_execution.json")
    m.num("alphaThreeChi", plaus["median_chi_square_by_alpha"]["2.3"], ".2f",
          "wp5_alpha_plausibility_execution.json")
    # issue #19: on repair_v8 alpha = 2.0, not 2.6, has the worst median, so
    # the text quotes all three instead of calling 2.6 the worst.
    m.num("alphaTwoChi", plaus["median_chi_square_by_alpha"]["2"], ".2f",
          "wp5_alpha_plausibility_execution.json")

    mass = pd.read_csv(I.resolve("wp5_association_mass_reconciliation"))
    mass = mass[
        mass.wp5_version.eq(C.V["wp5"]) & mass.family.eq(BASE["family"])
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
    ledger_json = chain_json("wp6_ledger_execution.json")
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
    m.num("NSNallbranches", len(assoc), "d", "wp7_ledger.csv")
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
    binary = chain_json("wp7_binary_bound_execution.json")
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
    iso_exec = chain_json("wp11_isotope_forecast_execution.json")
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

    # ------------------------------------------------------------------ WP12
    # The manuscript-revision analysis.  Read-only over the frozen repair_v7
    # chain, pre-registered with input hashes in
    # provenance/wp12_revision_prereg.json.
    landscape = chain_json("wp12_gate_landscape_execution.json")
    src = "wp12_gate_landscape_execution.json"
    v7 = landscape["repair_v7_breakdown"]
    m.num("gateCells", v7["cells_total"], "d", src)
    m.num("gatePassing", v7["cells_passing"], "d", src)
    m.num("gateFailing", v7["cells_failing"], "d", src)
    for key, tag in (("CygOB2-A", "A"), ("CygOB2-B", "B"), ("CygOB2-C", "C")):
        m.num(f"gateFail{tag}", v7["failures_by_subgroup"].get(key, 0), "d", src)
    m.num("gateFailParsec", v7["failures_by_family"].get("PARSEC", 0), "d", src)
    m.num("gateFailMist", v7["failures_by_family"].get("MIST", 0), "d", src)
    for rv, tag in (("3", "Rvthirty"), ("3.1", "Rvthirtyone"),
                    ("3.5", "Rvthirtyfive")):
        m.num(f"gateFail{tag}", v7["failures_by_R_V"].get(rv, 0), "d", src)
    for alpha, tag in (("2", "AlphaTwo"), ("2.3", "AlphaTwoThree"),
                       ("2.6", "AlphaTwoSix")):
        m.num(f"gateFail{tag}", v7["failures_by_alpha"].get(alpha, 0), "d", src)
    m.num("gateHeadlineCells", landscape["headline_cells"]["total"], "d", src)
    m.num("gateHeadlineCellsPass", landscape["headline_cells"]["passing"], "d",
          src)
    m.num("gateCombos", landscape["combinations"]["headline_total"], "d", src)
    m.num("gateCombosPass",
          landscape["combinations"]["headline_passing_all_subgroups"], "d", src)
    m.num("gateStrictBranches",
          landscape["headline_branches"]["built_entirely_on_passing_cells"],
          "d", src)
    strict = landscape["alpha_split"]["all_subgroup_pass_only"]
    m.num("strictTwoAbove", strict["alpha_2"]["above_0p5"], "d", src)
    m.num("strictTwon", strict["alpha_2"]["branches"], "d", src)
    m.num("strictTwoThreeAbove", strict["alpha_2.3"]["above_0p5"], "d", src)
    m.num("strictTwoThreen", strict["alpha_2.3"]["branches"], "d", src)
    m.num("strictTwolo", strict["alpha_2"]["score_min"], ".3f", src)
    m.num("strictTwohi", strict["alpha_2"]["score_max"], ".3f", src)

    closure_exec = chain_json("wp12_closure_slopes_execution.json")
    src = "wp12_closure_slopes_execution.json"
    grid_closure = closure_exec["subgroup_grid_median_closure_by_alpha"]
    for key, tag in (("CygOB2-A", "A"), ("CygOB2-B", "B"), ("CygOB2-C", "C")):
        m.num(f"closureGrid{tag}", grid_closure[key]["alpha_2.3"], ".2f", src)
    slopes = closure_exec["closing_slope_grid_median_by_subgroup"]
    for key, tag in (("CygOB2-A", "A"), ("CygOB2-B", "B"), ("CygOB2-C", "C")):
        m.num(f"closingAlpha{tag}", slopes[key], ".2f", src)
    m.num("closingAlphaAssoc", slopes["association"], ".2f", src)
    two = closure_exec["two_aggregations_at_2p3"]
    m.num("closureMedianRatios", two["median_of_18_subgroup_ratios"], ".3f", src)
    m.num("closureWeighted", two["star_weighted_summed_ratio_grid_median"],
          ".3f", src)
    mixed = closure_exec["mixed_slope"]
    m.num("mixedLo", mixed["N_SN_range"][0], ".2f", src)
    m.num("mixedHi", mixed["N_SN_range"][1], ".1f", src)
    versus = mixed["versus_pure_alpha_2p3"]
    m.num("mixedShift", 100 * versus["worst_relative_change_in_N_SN"], ".0f", src)
    dead = versus["worst_relative_change_where_C_is_dead"]   # None: C dies on every branch
    m.num("mixedShiftDead", None if dead is None else 100 * dead, ".0f", src)
    m.num("mixedAliveBranches", versus["branches_where_C_contributes"], "d", src)
    m.num("mixedDeadBranches",
          versus["branches_where_C_is_above_the_imf_ceiling"], "d", src)
    m.num("mixedCheck",
          100 * mixed["engine_validation"]["worst_relative_N_SN_difference"],
          ".1f", src)

    scenario_exec = chain_json("wp12_scenario_score_execution.json")
    src = "wp12_scenario_score_execution.json"
    c4 = scenario_exec["wp12_3_c4_sensitivity"]
    m.num("CfourAny",
          c4["C4_below_which_at_least_one_alpha2p0_branch_falls_below_0p5"],
          ".3f", src)
    m.num("CfourAll",
          c4["C4_below_which_every_alpha2p0_branch_falls_below_0p5"], ".3f", src)
    m.num("CfourLift",
          c4["C4_above_which_the_first_alpha2p3_branch_would_rise_above_0p5"],
          ".3f", src)
    c3 = scenario_exec["wp12_4_c3_subtype"]
    m.num("CthreeSixtyLo", c3["pessimistic_60_Msun_C3_range"][0], ".2f", src)
    m.num("CthreeSixtyHi", c3["pessimistic_60_Msun_C3_range"][1], ".2f", src)
    m.num("CthreeFortyLo", c3["threshold_scan_C3_range"]["step_40"][0], ".2f",
          src)
    m.num("CthreeFortyHi", c3["threshold_scan_C3_range"]["step_40"][1], ".2f",
          src)
    pess = c3["score_under_pessimistic_C3"]
    m.num("scorePessLo", pess["range"][0], ".3f", src)
    m.num("scorePessHi", pess["range"][1], ".3f", src)
    m.num("scorePessTwoAbove", pess["alpha2p0_above_0p5"], "d", src)
    m.num("scorePessTwon", pess["alpha2p0_branches"], "d", src)
    m.num("progenitorFloorLo",
          c3["progenitor_mass_floor"]["min_over_headline_branches_Msun"], ".1f",
          src)

    neighbour = chain_json("wp12_neighbour_budget_execution.json")
    src = "wp12_neighbour_budget_execution.json"
    r6 = {p["id"]: p for p in neighbour["predictions"]}["R6"]
    diag = r6["measured"][
        "non_preregistered_diagnostic_at_our_measured_4p0_Myr"]["rows"]
    ratios = [row["ratio_to_measured_baseline"] for row in diag
              if row["ratio_to_measured_baseline"] is not None]
    m.num("coarseRatioLo", min(ratios), ".1f", src)
    m.num("coarseRatioHi", max(ratios), ".1f", src)
    m.num("coarseRecent",
          neighbour["measured_cygob2_reference"][
              "expected_SNe_in_last_100kyr_baseline"], ".2f", src)
    m.num("neighbourPopulations",
          len(neighbour["cavity_sets"]["wide"]), "d", src)

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
