#!/usr/bin/env python3
"""Issue #20 Phase A -- pre-registration: is CygOB2-C really young?

Written and executed BEFORE any Phase A test (A2-A6) is run.  It fixes the
three hypotheses, every test statistic, the decision rule, the secondary
predictions, the integrity checks, and what each outcome does to N_death; it
records the SHA-256 of every input and of the method module
(scripts/issue20_common.py), so that a result cannot be computed from a moved
input or a quietly changed method.  It refuses to overwrite itself and refuses
to run if any Phase A output already exists.

Brief: tasks/issue20_subgroup_c_age_brief.md (sections 3-5) and
tasks/stage_after_issue19_c_age_and_wp13.md (step 1).

Run:
  PYTHONPATH=scripts python3 scripts/issue20_prereg.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy

import issue20_common as I
import wp5_common as w

OUT = I.PREREG_PATH

# Every product Phase A will write.  If any exists, the preregistration would
# postdate a result and is refused.
PHASE_A_OUTPUTS = [
    "tables/issue20_c_star_audit.csv",
    "tables/issue20_hrd_posteriors.csv",
    "tables/issue20_hrd_mixture.csv",
    "tables/issue20_hrd_per_star.csv",
    "tables/issue20_photometric_refits.csv",
    "tables/issue20_extinction_attribution.csv",
    "tables/issue20_closure_by_c_age.csv",
    "provenance/issue20_star_audit_execution.json",
    "provenance/issue20_hrd_likelihood_execution.json",
    "provenance/issue20_photometric_execution.json",
    "provenance/issue20_closure_execution.json",
    "provenance/issue20_phase_a_outcome.json",
    "reports/issue20_phase_a_report.md",
]

HYPOTHESES = {
    "H1": {
        "name": "C is young (~2.5 Myr)",
        "statement": (
            "C formed after A and B.  Its luminous supergiants are interlopers "
            "(A/B members in projection, or field) or binary/merger products."
        ),
    },
    "H2": {
        "name": "C is coeval with A/B (~3.5-4 Myr)",
        "statement": (
            "C's young photometric age is an artefact: evolved hot supergiants "
            "read as main-sequence stars by a colour-magnitude fit, possibly "
            "helped by the repair_v1 extinction change.  The O3 If star is the "
            "anomaly (rejuvenated merger product / blue straggler)."
        ),
    },
    "H3": {
        "name": "C has two components",
        "statement": (
            "A younger episode (the O3 stars) superposed on an older one (the B "
            "supergiants): an extended or multi-episode history, as Berlanas "
            "et al. (2019, 2020) report for the region."
        ),
    },
}

# ------------------------------------------------------------- A3 method
A3 = {
    "sample": (
        "Stars of wp4_anchor_hrd_repair_v8 labelled CygOB2-A, -B or -C, with a "
        "finite M_G0 on the branch's R_V, excluding extreme_hot (logTe > "
        "log10 52000) and excluding stars whose primary spectral type is not O "
        "or B (WR types).  B (5 anchors) is reported, never decided on."
    ),
    "teff": (
        "issue20_common.assign_teff.  The anchor table's T_eff, except that a "
        "TYPE-DERIVED T_eff (log10 T on the 0.01-dex grid of the Wright+15 "
        "spectral-type scale) of a luminosity class I or II star is replaced "
        "by the supergiant scale: Martins+2005 Table 6 (O3-O9.5, class I, "
        "observational) and Crowther+2006 Table 4 'mean' (O9.7-B3, Galactic "
        "Ia/Iab), linear in log T between tabulated subtypes.  Measured T_eff "
        "(Berlanas+2020 quantitative spectroscopy, off the 0.01-dex grid) are "
        "kept.  Primary component's type for composites."
    ),
    "sigma_logTe": {"primary": "0.03 (classes III-V, unclassified); 0.05 (I-II)",
                    "sensitivity": "0.05 for every star"},
    "M_G0": "MG0_obs_rv{R_V} of wp4_anchor_hrd_repair_v8 (= repair_v5 WP3 G0_abs).",
    "sigma_MG0": (
        "sqrt(G_err^2 + (k_G av_err)^2 + 0.0607^2 [WP2 depth 45.4 pc] + "
        "0.25^2 [intrinsic-colour / A_V calibration floor]); sensitivity with "
        "0.40 in place of 0.25 (WP4's lumped SIG_MG0).  G_err NaN -> 0.02 as WP4."
    ),
    "model": (
        "issue20_common.hrd_particles: every native isochrone age of the "
        "family (PARSEC 1.00-10.0 Myr, MIST 1.00-10.12 Myr, 21 each); all "
        "evolutionary phases; consecutive tabulated points joined linearly in "
        "initial mass and sub-sampled to 0.005 dex / 0.025 mag; IMF weight "
        "m^-2.3 per initial-mass interval, primaries >= 2 Msun; MIST gaps "
        "(phase change with a > 0.05 dex or > 0.5 mag jump) carry no weight; "
        "unresolved binaries with f_bin of the branch, q ~ U(0.1, 1) on 10 "
        "values, G flux-added, logTe of the primary.  No 120 Msun cap: the "
        "WP4 age fit uses the full isochrone and so does this."
    ),
    "statistic_primary": (
        "CONDITIONAL likelihood ln p(logTe_i | M_G0_i, t) = "
        "ln sum_j w_j N(logTe_i; T_j) N(M_i; G_j) - ln sum_j w_j N(M_i; G_j).  "
        "Reason (fixed before any likelihood was computed): the spectroscopic "
        "anchors are a luminosity-selected subset whose selection differs "
        "between subgroups (median M_G0: A -2.15, C -4.18), and a joint "
        "density rewards ages with more luminous stars, so it would bias the "
        "more luminous-selected subgroup (C) young.  Conditioning on M_G0 makes "
        "a magnitude-dependent selection ignorable."
    ),
    "statistic_secondary": (
        "JOINT likelihood within the window M_G0 <= 0.0 (just below the "
        "faintest anchor in A or C, -0.26), renormalised over the particles "
        "inside the window as wp4_common.star_loglike does.  Reported in full; "
        "if its verdict differs from the primary the verdict is flagged "
        "'statistic-dependent' but follows the primary."
    ),
    "subgroup_posterior": (
        "membership-weighted total lnL(t) = sum_i P_i ln L_i(t) on the native "
        "grid; posterior by issue20_common.posterior_summary (PCHIP in log "
        "age, uniform prior in log age over the grid, as "
        "wp4_common.posterior_from_loglike); equal-tailed 68 % and 95 % "
        "intervals.  A posterior whose MAP is within one native step of a grid "
        "edge is 'railed'."
    ),
    "mixture_H3": (
        "issue20_common.mixture_best for C (and A as specificity control): two "
        "ages t1 < t2 on the native grid, birth fraction phi in 0.02..0.98 "
        "(step 0.02); conditional of the mixture population "
        "[phi e^J1 + (1-phi) e^J2] / [phi e^M1 + (1-phi) e^M2]; each component "
        "must hold >= 3 membership-weighted stars by responsibility.  "
        "Delta BIC = 2 (lnL_mix - lnL_single) - 2 ln(sum P_i), k = 3 vs 1."
    ),
    "per_star_ages": (
        "each star's own ln L_i(t) through posterior_summary; 'defined' if "
        "max - min of ln L_i over the grid >= 2 and not railed.  Report only."
    ),
    "settings_grid": {
        "decision_cell": "R_V 3.1, f_bin 0.4, conditional, primary sigmas, dmu 0; PARSEC and MIST",
        "robustness": "R_V 3.0 and 3.5; f_bin 0.3 and 0.5 (both families)",
        "sensitivity_report_only": [
            "joint statistic (window M_G0 <= 0.0)",
            "sigma_logTe 0.05 for every star",
            "sigma_MG0 calibration floor 0.40",
            "distance modulus +/- 0.060 mag (coherent shift of the data)",
            "f_bin 0.7 (Sana+2012 close-binary fraction, issue #15)",
            "table T_eff throughout (no supergiant replacement)",
        ],
    },
}

DECISION_RULE = {
    "evaluated_on": (
        "the decision cell, separately for PARSEC and MIST.  A hypothesis is "
        "the verdict only if it holds in BOTH families."
    ),
    "H3": (
        "C's best two-age mixture beats C's best single age by Delta BIC > 10 "
        "AND age_old - age_young > 1.0 Myr AND both components hold >= 3 "
        "weighted stars."
    ),
    "H1": "C's 95 % upper bound < A's 95 % lower bound.",
    "H2": "C's 68 % interval overlaps A's 68 % interval.",
    "precedence": (
        "H3 is tested first: if it holds in both families the verdict is H3 "
        "whatever the single-age intervals say.  Otherwise H1, then H2 (H1 and "
        "H2 are mutually exclusive by construction).  Families disagreeing, "
        "or neither H1 nor H2 holding (68 % disjoint, 95 % overlapping), is "
        "INCONCLUSIVE."
    ),
    "measurability": (
        "if A's or C's posterior is railed in a family, that family supports "
        "neither H1 nor H2 (it can still support H3); a railed family makes "
        "the single-age verdict INCONCLUSIVE."
    ),
    "qualifiers_not_changing_the_label": {
        "robust / fragile": (
            "robust if the same per-family outcome obtains at R_V 3.0, 3.5 and "
            "f_bin 0.3, 0.5 in both families; otherwise fragile, with the "
            "dissenting cells listed."
        ),
        "statistic-dependent": "the secondary (joint) statistic gives a different verdict.",
        "H3 specificity": (
            "the same mixture test is run on A.  If A also passes the H3 "
            "criteria in both families, an H3 verdict is labelled 'mixture "
            "preference not specific to C' and Phase B for H3 needs a user "
            "decision (the test may be detecting model mis-specification, e.g. "
            "blue stragglers, rather than a second episode in C)."
        ),
    },
    "inconclusive_means": (
        "C is carried as a branch (young / coeval) through the ledger in a "
        "later, separately pre-registered chain, and no novelty claim rests on "
        "C's age."
    ),
}

CONSEQUENCES_FOR_N_DEATH = {
    "baseline_now": "repair_v8: N_death 8.36 (A 4.10, B 4.26, C 0.00), P(last < 100 kyr) 0.547",
    "H1": (
        "No change.  C stays at its repair_v5 age prior (2.52 Myr counts-based) "
        "and contributes 0 at baseline; the ledger stands; WP13's premise is "
        "strengthened and Phase C (focus study of C) opens."
    ),
    "H2": (
        "Direction UP.  In a separately pre-registered repair_v9, run after "
        "the WP13 decision (stage brief step 1), C's age prior moves toward "
        "A's (~3.5-4 Myr) from a hybrid CMD+HRD likelihood; C's turnoff drops "
        "below 120 Msun, C contributes deaths and the baseline rises above "
        "8.36; P(last < 100 kyr) is expected to rise toward the common-age "
        "values (WP13 brief section 2.1).  The pooled-vs-resolved difference "
        "is expected to shrink; WP13's 'equivalent' outcome becomes likely."
    ),
    "H3": (
        "Direction UP, by less than H2 for the same old age: the older "
        "component contributes deaths in proportion to its birth fraction, "
        "dated by its own age.  Requires either the existing formation-"
        "duration branch (delta up to 2 Myr) or a two-component C in repair_v9."
    ),
    "INCONCLUSIVE": (
        "Baseline unchanged; a coeval-C branch widens the headline range "
        "upward.  No novelty claim rests on C's age."
    ),
    "note": (
        "These are expected directions, not targets.  Phase A changes no chain "
        "product; repair_v8 remains the quoted chain whatever the verdict."
    ),
}

# --------------------------------------------- secondary predictions (scored)
SECONDARY_PREDICTIONS = [
    {
        "id": "S1",
        "from": "A2",
        "statement": (
            "Kinematic membership of C's luminosity class I/II spectroscopic "
            "stars.  H1 predicts they are closer in proper motion to A or B "
            "than to C; H2/H3 predict they are kinematic C members."
        ),
        "measured_quantity": (
            "for each C anchor of class I/II, the proper-motion Mahalanobis "
            "distance to the A, B and C centroids (membership-weighted mean "
            "and covariance of pmra, pmdec over the labelled members); the "
            "fraction whose nearest centroid is C"
        ),
        "scoring": "fraction >= 2/3 -> 'consistent with H2/H3'; <= 1/3 -> 'consistent with H1'; otherwise neutral",
        "caveat": (
            "partly circular: the labels come from a GMM in (l, b, pmra, "
            "pmdec), so labelled C stars are near C by construction.  A result "
            "'consistent with H1' would be informative; 'consistent with H2/H3' "
            "is weak evidence."
        ),
    },
    {
        "id": "S2",
        "from": "A4(b)",
        "statement": (
            "H1 predicts that removing C's class I/II stars keeps a young "
            "photometric age."
        ),
        "measured_quantity": (
            "WP4 upper-MS MAP of C (wp4_common machinery, repair_v5 "
            "extinction, f_bin 0.4, dmu 0) with C's class I/II anchors removed, "
            "at R_V 3.1, each family"
        ),
        "scoring": (
            "'young persists' if MAP_C < MAP_A(repair_v5, same family, R_V "
            "3.1, f_bin 0.4) - 0.5 Myr in BOTH families (consistent with H1); "
            "otherwise not"
        ),
    },
    {
        "id": "S3",
        "from": "A4(a)",
        "statement": "H2 predicts the hybrid CMD+HRD fit gives C >= 3.3 Myr in both families.",
        "measured_quantity": (
            "MAP of the hybrid likelihood for C at R_V 3.1, f_bin 0.4: "
            "photometric-only members by wp4_common.star_loglike in the upper-"
            "MS window (M_G0 <= 1.5); spectroscopic members by the joint "
            "(logTe, M_G0) density in the same window from "
            "issue20_common.hrd_particles, with the A3 T_eff and sigmas"
        ),
        "scoring": "'old' if MAP >= 3.3 Myr in both families (consistent with H2)",
    },
    {
        "id": "S4",
        "from": "A6",
        "statement": (
            "Expected direction: an older C lowers the turnoff, so fewer massive "
            "stars are predicted alive and C's closure excess grows."
        ),
        "measured_quantity": (
            "C's closure ratio at alpha 2.3 for C ages 2.51, 3.16, 3.55, 3.98 "
            "Myr on each family x R_V cell, and the 6-cell median closing slope "
            "at each age"
        ),
        "scoring": (
            "PASS if the alpha-2.3 ratio is non-decreasing in age on all six "
            "cells AND the median closing slope at 3.98 Myr is below that at "
            "2.51 Myr"
        ),
        "interpretation": (
            "recorded as evidence either way: if it passes, coeval C needs an "
            "even shallower top-end IMF in C (harder to explain); it does not "
            "by itself score H1/H2."
        ),
    },
]

REPORT_ONLY = {
    "A2": "star-by-star audit table, incl. RUWE, composite types, SIMBAD binarity object types if retrievable",
    "A3_per_star": "per-star ages; bimodality of C's defined per-star ages (H3 expectation) described, not scored",
    "A3_B": "B's 5-anchor posterior",
    "A5": (
        "attribution of the 3.98 -> 2.51 Myr jump (PARSEC, R_V 3.1): per-star "
        "Delta A_V and Delta M_G0 between pre-repair and repair_v1 for C's 43 "
        "anchors and 15 brightest members; counterfactual WP4 upper-MS fits of "
        "C on repair_v1 extinction with (i) the anchors' and (ii) the non-"
        "anchors' photometry taken from the pre-repair run.  No forecast is "
        "registered for A5."
    ),
    "A7": "literature ages (Wright+2015, Berlanas+2019/2020, O3 If studies) -- cross-check only, never a calibration",
}

A6_METHOD = {
    "what_varies": (
        "only the upper limit of C's closure integral: cap = min(turnoff(family, "
        "t), 120 Msun) with turnoff from wp6_mass_extension_decision.turnoff_mass, "
        "applied to every one of C's stored node responses."
    ),
    "held_fixed": (
        "C's stored repair_v7 node responses (WP5 nodes + WP6 extension; the "
        "repair_v8 chain reuses them) and their WP4 node weights, C's repair_v8 "
        "k_median per alpha, and C's repair_v8 observed census."
    ),
    "why_not_more": (
        "C's injection responses exist only at its own posterior nodes "
        "(2.0-3.2 Myr); varying C's age properly would need new injections at "
        "3.16-3.98 Myr and a refit of k at those ages.  That is a repair_v9 "
        "task, not a read-only Phase A step.  A6 therefore isolates the "
        "turnoff-window effect only; the age dependence of k and of the mass "
        "estimates of C's observed stars is not captured."
    ),
    "closing_slope": "wp6_closure_attribution.closing_alpha over alpha in {2.0, 2.3, 2.6}",
}

INTEGRITY_CHECKS = [
    {"id": "I1", "check": (
        "Re-running the WP4 upper-MS fit (wp4_common, unchanged) for C on "
        "repair_v5 extinction reproduces wp4_age_posteriors_repair_v5's age_map "
        "for both families and all three R_V at f_bin 0.4, dmu 0, within "
        "0.01 Myr.  Guards A4.")},
    {"id": "I2", "check": (
        "The same fit on the pre-repair and repair_v1 extinction reproduces "
        "C's stored PARSEC R_V 3.1 age_map in wp4_age_posteriors and "
        "wp4_age_posteriors_repair_v1 within 0.01 Myr.  Guards A5.")},
    {"id": "I3", "check": (
        "The scope diagnostic's nearest-point chi^2 scan, recomputed from the "
        "VERSIONED anchor table (table T_eff, R_V 3.1), reproduces "
        "tables/issue20_hrd_age_scan.csv chi2_all for A and C within 0.1.  "
        "Confirms the versioned table carries the same spectroscopy.")},
    {"id": "I4", "check": (
        "A6 with C's stored node prior and the stored turnoff caps reproduces "
        "C's closure_ratio rows of tables/wp6_closure_repair_v8.csv within "
        "1e-9 relative.")},
    {"id": "I5", "check": (
        "hrd_particles: kept + dropped IMF weight equals the analytic IMF "
        "integral over [2 Msun, max Mini] within 1e-9 at every age; dropped "
        "(gap) fraction reported.")},
    {"id": "I6", "check": "every input matches its preregistered SHA-256 at run time (frozen())."},
]

DISCLOSED_PRIOR_KNOWLEDGE = [
    "The scope diagnostic (scripts/issue20_c_age_diagnostic.py, 2026-10-01): a "
    "nearest-point spectroscopic-HRD scan (no IMF weighting, no binary model, "
    "table T_eff, reduced chi^2 ~3-4) prefers 3.55 Myr (PARSEC) / 4.01 Myr "
    "(MIST) for C over 2.51 (sum chi^2 144.0 vs 176.8 PARSEC; 134.8 vs 292.7 "
    "MIST); for C's 9 brightest, 16.2 at 3.55 vs 29.5 at 2.51 (PARSEC).  A's "
    "control prefers 3.55 / 3.57.",
    "C's upper-MS photometric MAP (PARSEC, R_V 3.1, f_bin 0.4): pre-repair "
    "3.98, repair_v1/v3/v5 2.51; MIST 3.57 -> 3.18.  repair_v5 grid: C "
    "2.03-3.17 Myr counts-based, younger than A (3.85-4.07) on every cell.",
    "The previous agent expected H2 or H3 to be plausible (stage brief step 1).",
    "While designing this test (2026-10-02, before any likelihood was "
    "computed) the present agent saw: the list of C's 43 anchors with "
    "spectral types, table T_eff and M_G0 (C's luminous stars are B0 Ia x2 at "
    "-7.51/-6.98, O5 I+O3.5 III -7.29, B1 Ib -7.13, O9.7 Iab -5.92, O3 If "
    "-5.86, B1 I -4.64); the A and B anchor lists likewise (A's brightest is "
    "O7 I+O6 I+O9 V at -8.53; B holds a WN7o+O7 V star, excluded as WR); that "
    "type-derived T_eff of B-type supergiants in the anchor table are DWARF "
    "values (B0 Ia 30.9 kK, B1 Ib 24.5 kK), which the supergiant scale "
    "lowers to 27.5 and 21.5 kK; that the PARSEC 3.98 Myr isochrone carries "
    "post-MS stars at logTe 4.40-4.45 and G -7.0 to -7.7 (initial mass "
    "43-47 Msun), i.e. near C's B0 Ia stars; that MIST 2.83-3.57 Myr tables "
    "omit post-MS phases between the end of the MS and the WR phase; that "
    "the anchors' tabulated av_err is ~0.008 mag (hence the 0.25 mag floor).",
    "The present agent's own expectation, recorded so it cannot be claimed "
    "as a prediction: the B supergiants pull C older and the O3 If star pulls "
    "it younger; A's anchors are mostly B dwarfs with weak age leverage, so "
    "A's posterior may be broad, which makes H2's overlap criterion easy to "
    "meet and H1's hard; H2 or INCONCLUSIVE judged more likely than H1, H3 "
    "possible.  These expectations are not scored.",
    "No Phase A likelihood, refit, closure scan or star audit had been run "
    "when this file was written.  The particle sampler was exercised on "
    "isochrones alone (weights, gaps, timing), never on stars.",
]


REUSED_CODE = [
    "scripts/wp4_common.py",                   # A4, A5: the WP4 upper-MS likelihood
    "scripts/wp4_anchors_hrd.py",              # I3: the diagnostic's nearest-point metric
    "scripts/wp3_extinction_law.py",           # k_G
    "scripts/wp6_mass_extension_decision.py",  # A6: turnoff_mass
    "scripts/wp6_closure_test.py",             # A6: closure estimator (replicated)
    "scripts/wp6_closure_attribution.py",      # A6: closing_alpha
    "scripts/wp5_joint_age_fit.py",            # A6: node prior and response paths
    "scripts/wp6_massive_injections.py",       # A6: extension response paths
]


def method_code_hashes() -> dict:
    return {"scripts/issue20_common.py": w.sha256(I.THIS_FILE)}


def reused_code_hashes() -> dict:
    return {rel: w.sha256(w.ROOT / rel) for rel in REUSED_CODE}


def main() -> None:
    if OUT.exists():
        raise SystemExit(
            f"{OUT.relative_to(w.ROOT)} already exists -- a preregistration is "
            "written once and never regenerated."
        )
    existing = [p for p in PHASE_A_OUTPUTS if (w.ROOT / p).exists()]
    if existing:
        raise SystemExit(
            "Phase A outputs already exist; a preregistration written after "
            f"them would be worthless: {existing}"
        )
    inputs = {}
    for name, rel in I.INPUTS.items():
        path = w.ROOT / rel
        if not path.exists():
            raise SystemExit(f"input missing at preregistration: {rel}")
        inputs[name] = {"path": rel, "sha256": w.sha256(path)}

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_prereg.py",
        "issue": "#20 -- is CygOB2-C really young?  Phase A",
        "briefs": ["tasks/issue20_subgroup_c_age_brief.md",
                   "tasks/stage_after_issue19_c_age_and_wp13.md"],
        "chain": "repair_v8 (unchanged by Phase A; Phase A is read-only)",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__,
                     "scipy": scipy.__version__},
        "hypotheses": HYPOTHESES,
        "a3_method": A3,
        "a3_constants": {
            "supergiant_teff_K": {str(k): v for k, v in I.SUPERGIANT_TEFF.items()},
            "sig_logTe_dwarf_giant": I.SIG_LOGTE_DWARF_GIANT,
            "sig_logTe_supergiant": I.SIG_LOGTE_SUPERGIANT,
            "sig_depth_mag": I.SIG_DEPTH_MAG, "sig_cal_mag": I.SIG_CAL_MAG,
            "sig_cal_mag_wide": I.SIG_CAL_MAG_WIDE, "imf_slope": I.IMF_SLOPE,
            "m_min": I.M_MIN, "m_birth_max": I.M_BIRTH_MAX, "q_min": I.Q_MIN,
            "n_q": I.N_Q, "dlogte_step": I.DLOGTE_STEP, "dmag_step": I.DMAG_STEP,
            "mist_gap_dlogte": I.MIST_GAP_DLOGTE, "mist_gap_dmag": I.MIST_GAP_DMAG,
            "phi_grid": [float(x) for x in I.PHI_GRID],
            "min_component_stars": I.MIN_COMPONENT_STARS,
            "joint_window_edge_mag": 0.0, "hybrid_window_edge_mag": 1.5,
        },
        "decision_rule": DECISION_RULE,
        "consequences_for_n_death": CONSEQUENCES_FOR_N_DEATH,
        "secondary_predictions": SECONDARY_PREDICTIONS,
        "report_only": REPORT_ONLY,
        "a6_method": A6_METHOD,
        "integrity_checks": INTEGRITY_CHECKS,
        "failed_predictions_rule": (
            "Failed predictions stay recorded as failed.  This file is not "
            "amended; corrections go in the outcome record with the reason.  "
            "A change to scripts/issue20_common.py after this point is recorded "
            "in provenance/issue20_deviations.json before it is used."
        ),
        "phase_b_rule": (
            "If H2 or H3 changes C's age prior, that is a separate, separately "
            "pre-registered chain version (repair_v9), run AFTER this stage's "
            "WP13 decision, not folded into repair_v8."
        ),
        "disclosed_prior_knowledge": DISCLOSED_PRIOR_KNOWLEDGE,
        "method_code": method_code_hashes(),
        "reused_code": reused_code_hashes(),
        "consumed_inputs": inputs,
        "planned_outputs": PHASE_A_OUTPUTS,
    }
    w.write_json(OUT, record)
    print(f"wrote {OUT.relative_to(w.ROOT)}")
    print(f"  inputs hashed: {len(inputs)}")
    print(f"  method code  : {record['method_code']}")


if __name__ == "__main__":
    main()
