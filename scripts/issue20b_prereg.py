#!/usr/bin/env python3
"""Issue #20 Phase A' -- pre-registration: three independent age tests for A, B and C.

Owner decision of 2026-10-02 (provenance/decisions_2026_10_02.json, item 4):
Phase A gave contradictory ages for C (photometric 2.51 Myr vs spectroscopic
3.98 / 4.01 Myr), so the age of every subgroup is re-measured by three tests
that do not share the disputed ingredient (the WP3 photometric extinction of
the stars without spectra), and the age adopted by a pre-registered agreement
rule is used throughout the study (repair_v9, before WP13).

Written and executed BEFORE any of tests a, b, c is run on real stars.
Refuses to overwrite itself or to run if any Phase A' output exists.

Run:
  PYTHONPATH=scripts python3 scripts/issue20b_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy

import issue20b_common as B
import wp5_common as w

OUT = B.PREREG_PATH

OUTPUTS = [
    "tables/issue20b_age_tests.csv",
    "tables/issue20b_calibrations.csv",
    "tables/issue20b_adopted_ages.csv",
    "provenance/issue20b_run_execution.json",
    "provenance/issue20b_outcome.json",
    "reports/issue20b_phase_a_prime_report.md",
]

TESTS = {
    "a": {
        "name": "free-extinction CMD + near-IR likelihood",
        "independent_of": "the WP3 extinction map and its anchor-based spatial prior (the "
                          "ingredient that moved C from 3.98 to 2.51 Myr at repair_v1)",
        "data": "labelled members with 2MASS J and Ks of quality A or B; Gaia G, BP, RP",
        "likelihood": (
            "issue20b_common.cmd_loglike, mode 'nir'.  For each model particle j of the trial "
            "age (IMF-weighted, all phases, unresolved binaries, f_bin 0.4, initial mass >= 1 "
            "Msun), the star's A_V is integrated analytically with a flat prior using its J-Ks "
            "against the particle's intrinsic J0-Ks0; the (G - mu, BP-RP) residual is then a 2-D "
            "Gaussian with the A_V uncertainty propagated (full covariance).  Window-normalised "
            "over particles with G0 <= 1.5 (as WP4).  Membership-weighted sum over stars."
        ),
        "sample": "fixed across ages: M_G0 <= 1.5 from the calibrated near-IR excess with (J-Ks)_0 = -0.20",
        "nir_law_calibration": (
            "issue20b_common.nir_scale per R_V: s = median over O-B3 anchors of "
            "[(J-Ks + 0.20)/(k_J - k_Ks)] / A_V,spec; the effective k_x = (k_J - k_Ks)/s.  "
            "The anchors' intrinsic-colour A_V is the project's reference extinction scale "
            "(WP3); it is not the disputed photometric A_V."
        ),
        "validity": "per R_V: >= 30 calibration anchors, |s - 1| <= 0.3, robust scatter <= 0.5 mag; "
                    "per subgroup: >= 10 stars; posterior not railed",
        "control_report_only": (
            "the same likelihood with A_V fixed at the repair_v5 WP3 value +- av_err (mode "
            "'wp3').  The difference test a - control is the effect of the extinction source "
            "alone.  The control does not gate anything."
        ),
    },
    "b": {
        "name": "Gaia DR3 ESP-HS hot-star HRD",
        "independent_of": "both the WP3 extinction and the literature spectral types of the "
                          "stars it measures (T_eff and A_G come from Gaia XP spectra)",
        "data": ("gaiadr3.astrophysical_parameters ESP-HS columns, queried 2026-10-02 for all "
                 "1,392 members (data/raw/gaia/issue20b_gaia_dr3_astrophysical_parameters.csv); "
                 "stars with teff_esphs, ag_esphs and spectraltype_esphs in {O, B}"),
        "calibration": (
            "issue20b_common.esphs_calibration per R_V: zero points of log T_eff (vs the anchors' "
            "table T_eff) and A_G (vs k_G x anchors' intrinsic-colour A_V) from anchors with "
            "class III-V or unclassified; robust scatters added in quadrature to every star"
        ),
        "likelihood": "conditional ln p(logTe | M_G0, t) of issue20_common with the outlier term "
                      "(eps 0.05, uniform in logTe over 4.10-4.75); window M_G0 <= 1.5",
        "validity": "per R_V: >= 20 calibration anchors, sigma(logTe) <= 0.08 dex, sigma(A_G) <= 0.6 mag; "
                    "per subgroup: >= 10 stars; posterior not railed",
        "coverage_seen_before_registration": "ESP-HS T_eff exists for 80 members: A 25 (14 anchors), "
                                             "B 22 (1), C 15 (5), unlabelled 18 (10); 30 anchors in total",
    },
    "c": {
        "name": "outlier-robust spectroscopic HRD",
        "likelihood": "the issue #20 A3 method (issue20_common, unchanged) with the outlier term "
                      "of test b applied to every star's conditional likelihood",
        "validity": "per subgroup >= 10 stars (B, with 4, is invalid); posterior not railed",
        "not_blind": ("test c re-uses Phase A's data and method; its result is expected to be "
                      "close to Phase A's (disclosed below).  The blind tests are a and b."),
    },
}

DECISION_RULE = {
    "cells": "every subgroup x family x R_V at f_bin 0.4; the headline reading is R_V 3.1",
    "agreement": "two tests agree when their 68 % intervals overlap",
    "adoption": (
        "issue20b_common.adopt: among the VALID tests, take the largest set (size 3, else 2) whose "
        "members pairwise agree.  If exactly one such set exists the subgroup's age is ADOPTED, "
        "and the adopted posterior is the equal-weight mixture of that set's posteriors "
        "(issue20b_common.mixture_posterior).  If none exists, or two different sets of the same "
        "size exist, the cell is NOT RESOLVED."
    ),
    "c_versus_a": (
        "from the adopted posteriors at R_V 3.1: 'C younger than A' iff C's 95 % upper bound < A's "
        "95 % lower bound in BOTH families; 'C coeval with A' iff the 68 % intervals overlap in "
        "BOTH families; otherwise 'C's age relative to A undetermined'.  Requires both A and C "
        "ADOPTED in both families."
    ),
    "use_throughout": (
        "ADOPTED posteriors (MAP, 68 %) replace the upper-MS age rows of the corresponding "
        "subgroup/family/R_V in a repair_v9 age table, which repair_v9 (separately pre-"
        "registered) propagates through WP5-WP12 before WP13 reads M0/M1.  A NOT RESOLVED cell "
        "keeps its repair_v5 photometric posterior as one branch and carries test c's posterior "
        "as a second, explicitly labelled branch."
    ),
    "robustness_qualifier": "the R_V 3.1 adoption is 'robust' if R_V 3.0 and 3.5 give the same "
                            "status and agreeing set in that family; otherwise 'R_V-fragile'",
}

SECONDARY = {
    "H3_with_guard": (
        "test c's two-age mixture for each subgroup with >= 10 stars (issue20_common.mixture_best "
        "on the outlier-robust likelihood): counts as a two-component detection only if Delta BIC "
        "> 10, age_old - age_young > 1 Myr, both components >= 3 stars, AND neither component age "
        "lies within one native grid step of the grid edges.  Reported for A and C, both "
        "families; it does not enter the adoption rule."
    ),
    "eps_sensitivity": "tests b and c repeated with eps 0.02 and 0.10 (report only)",
}

INTEGRITY = [
    {"id": "J1", "check": "synthetic populations drawn from the model at 2.51 and 3.98 Myr (both "
                          "families, A_V uniform 4-8 mag, realistic noise, 400 stars) are recovered by "
                          "test a's MAP within one native grid step"},
    {"id": "J2", "check": "test c with eps = 0 reproduces Phase A's decision-cell conditional MAPs "
                          "for A and C exactly (tables/issue20_hrd_posteriors.csv)"},
    {"id": "J3", "check": "every input matches its SHA-256 here; method code matches or a deviation "
                          "is recorded"},
]

CONSEQUENCES = {
    "C ADOPTED and coeval with A": (
        "repair_v9 moves C's age prior to ~A's; C's turnoff falls below 120 Msun, C contributes "
        "deaths, baseline N_death rises above 8.36; the surviving novelty claim (C below the "
        "first-death boundary) is withdrawn; WP13 is expected to come out 'equivalent' or "
        "'supported but immaterial'"
    ),
    "C ADOPTED and younger than A": (
        "C's young age stands, now independently confirmed; ledger direction unchanged; Phase C "
        "(focus study of C) opens; WP13's premise is strengthened"
    ),
    "C NOT RESOLVED": "C carried as two branches in repair_v9; no novelty claim on C's age",
    "A or B ages": "any adopted change to A or B also enters repair_v9; direction stated in the outcome",
}

DISCLOSED = [
    "Phase A (provenance/issue20_phase_a_outcome.json): spectroscopic-HRD MAP A 3.16 / 3.57, C 3.98 / "
    "4.01 Myr (PARSEC / MIST); photometric repair_v5 MAP A 3.98 / 4.01, C 2.51 / 3.18; "
    "the 3.98 -> 2.51 jump is carried by C's non-anchor members' repair_v1 extinction; exploratory "
    "(post-hoc, 2026-10-02): C's spectroscopic age without its two A-like supergiants 3.55 / 4.01 "
    "and without all class I/II stars 3.98 / 3.57.",
    "ESP-HS coverage counts (above) were inspected before registration; no ESP-HS T_eff or A_G was "
    "compared with any anchor or used in any age computation.",
    "GSP-Phot was inspected for coverage only: most hot members were fitted with cool-star "
    "libraries (PHOENIX 167, A 156, MARCS 64, OB 12), so GSP-Phot is not used.",
    "Test a's machinery was run on SYNTHETIC populations only before registration (recovered 2.51, "
    "3.16 and 3.98 Myr within one grid step, PARSEC and MIST).  It has not been run on real stars, "
    "and the near-IR scale s has not been computed.",
    "Isochrone check (no stars): on the MS the intrinsic J-Ks of hot stars is -0.24 to -0.17 for "
    "M_Ks0 < -0.5 and nearly age-independent; at the faint end of the window PMS stars are much "
    "redder, which test a handles by taking intrinsic colours from each age's own model particles.",
    "Expectation of the agent, not scored: test c ~ Phase A A3 (C at or above A's age); test a is the "
    "real arbiter of the photometric-vs-spectroscopic conflict; test b is thin (C 15 stars) and "
    "may be broad or invalid.",
]


def main() -> None:
    if OUT.exists():
        raise SystemExit("provenance/issue20b_prereg.json exists -- written once")
    existing = [p for p in OUTPUTS if (w.ROOT / p).exists()]
    if existing:
        raise SystemExit(f"Phase A' outputs already exist: {existing}")
    inputs = {n: {"path": rel, "sha256": w.sha256(w.ROOT / rel)} for n, rel in B.INPUTS.items()}
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20b_prereg.py",
        "issue": "#20 Phase A' -- independent ages for A, B and C",
        "decision_record": "provenance/decisions_2026_10_02.json",
        "python": sys.version.split()[0], "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "tests": TESTS,
        "constants": {"eps_outlier": B.EPS_OUTLIER, "outlier_logTe": list(B.OUTLIER_LOGTE),
                      "min_stars": B.MIN_STARS, "jks0_hot": B.JKS0_HOT, "sig_jks_int": B.SIG_JKS_INT,
                      "cmd_m_min": B.CMD_M_MIN, "cmd_dmag_step": B.CMD_DMAG_STEP,
                      "cmd_dlogte_step": B.CMD_DLOGTE_STEP, "f_bin": B.F_BIN, "ums_edge": B.UMS_EDGE,
                      "nir_quality": sorted(B.NIR_QUAL)},
        "decision_rule": DECISION_RULE,
        "secondary": SECONDARY,
        "integrity_checks": INTEGRITY,
        "consequences": CONSEQUENCES,
        "disclosed_prior_knowledge": DISCLOSED,
        "method_code": {rel: w.sha256(w.ROOT / rel) for rel in (
            "scripts/issue20b_common.py", "scripts/issue20_common.py",
            "scripts/issue20_hrd_likelihood.py")},
        "consumed_inputs": inputs,
        "planned_outputs": OUTPUTS,
        "failed_predictions_rule": "results stay as they come out; this file is never amended",
    }
    w.write_json(OUT, record)
    print(f"wrote {OUT.relative_to(w.ROOT)}; {len(inputs)} inputs hashed")


if __name__ == "__main__":
    main()
