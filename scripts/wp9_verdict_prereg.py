#!/usr/bin/env python3
"""Pre-registration: WP9, the verdict and result framing.

Written BEFORE the verdict is computed.  The verdict statistic, the
plausibility thresholds and the framing decision rule are fixed here.

WHAT WP9 IS
-----------
The project exists to test one hypothesis: Haerer et al. 2025's claim that a
supernova about 50 kyr ago in Cyg OB2 powers the Cygnus PeVatron.  WP0 extracted
that claim into FOUR conditions.  WP9 computes the probability that the ledger
supplies an event meeting them, chooses the paper's framing by a rule fixed in
advance, and survives a devil's-advocate pass.

A PARTIAL PRE-REGISTRATION, AND SAYING SO
-----------------------------------------
The AGE component is already published: WP7 reports P(last SN < 100 kyr) = 0.552
and a rate of 8.01/Myr, and WP8 reports P(SN within 7 kyr) = 0.055.  So the
thresholds below are being set with approximate knowledge of one of the three
computable terms.  This is stated rather than concealed.

Two things protect the exercise:

  * the thresholds are anchored to EXTERNAL conventional meanings -- 0.5 is
    "more likely than not", 0.1 is "strongly disfavoured" -- and not chosen
    relative to any measured value;
  * what is NOT yet known is the JOINT verdict, its behaviour across branches,
    and the in-situ and progenitor-type terms, which is what the framing rule
    actually turns on.

THE BRANCH SET (option B, chosen by the principal investigator)
--------------------------------------------------------------
Headline over alpha in {2.0, 2.3}.  alpha = 2.6 is REPORTED in the sensitivity
table as tested-and-excluded, not deleted: the WP5 calibration-window evidence
alone rejects it (best in 1 of 18 cells, worst median chi-square), so excluding
it does NOT spend the census closure, which stays free to serve as the
out-of-sample validation of the IMF.  See reports/wp5_alpha_plausibility.md.

Output: provenance/wp9_verdict_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp9_verdict_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

# ---- Haerer's four conditions, from wp0_requirements_table.md --------------
AGE_PREFERRED_KYR = 50.0
AGE_WINDOW_KYR = (10.0, 100.0)      # H25's stated floor to twice its preference
AGE_PERMISSIVE_KYR = (10.0, 500.0)  # H25's conclusion allows "few hundred kyr"
ENERGY_ERG = (3.0e51, 5.0e51)
SUPERBUBBLE_DIAMETER_PC = 150.0
STRIPPED_PROGENITOR_MSUN = 30.0     # above this, envelope loss -> type Ib/c

# ---- plausibility thresholds, anchored to conventional meanings ------------
PLAUSIBLE_THRESHOLD = 0.5     # "more likely than not"
IMPLAUSIBLE_THRESHOLD = 0.1   # "strongly disfavoured"

# ---- branch set -----------------------------------------------------------
HEADLINE_ALPHAS = (2.0, 2.3)
EXCLUDED_ALPHA = 2.6
PRIMARY_EXPLODABILITY = "all_explode"


def main() -> None:
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp9_verdict_prereg.py",
        "status": "PREREGISTERED",
        "work_package": "WP9",
        "objective": (
            "compute the exact statement the paper leads with, and choose the "
            "framing by a rule fixed before the number is known"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "partial_preregistration_disclosure": {
            "what_is_already_known": (
                "the AGE component is published: WP7 gives a rate of 8.01/Myr "
                "and P(last SN < 100 kyr) = 0.552; WP8 gives P(SN within "
                "7 kyr) = 0.055.  The thresholds below are therefore set with "
                "approximate knowledge of one of the three computable terms."
            ),
            "why_this_is_still_disciplined": (
                "the thresholds are anchored to EXTERNAL conventional meanings "
                "-- 0.5 is 'more likely than not', 0.1 is 'strongly "
                "disfavoured' -- rather than being chosen relative to a "
                "measured value.  What is NOT yet known is the joint verdict, "
                "its spread across branches, and the in-situ and "
                "progenitor-type terms, and the framing rule turns on those."
            ),
            "what_would_have_been_illegitimate": (
                "setting the threshold at, say, 0.45 after seeing 0.552, so "
                "that the verdict clears it.  The values 0.5 and 0.1 are fixed "
                "points that would have been chosen identically before WP7 ran."
            ),
        },
        "the_four_conditions": {
            "C1_age": {
                "requirement": (
                    f"H25 prefers about {AGE_PREFERRED_KYR:g} kyr, with a floor "
                    f"near 10-20 kyr because an older remnant should have "
                    f"faded, and allows an event within the last few hundred "
                    f"kyr in its conclusion"
                ),
                "primary_window_kyr": list(AGE_WINDOW_KYR),
                "permissive_window_kyr": list(AGE_PERMISSIVE_KYR),
                "statistic": (
                    "P(at least one supernova with look-back time inside the "
                    "window), computed from the WP7 explosion epochs -- NOT "
                    "P(the LAST supernova is in the window), because any event "
                    "in the window can power the PeVatron"
                ),
                "computable": True,
            },
            "C2_energy": {
                "requirement": f"{ENERGY_ERG[0]:.0e} to {ENERGY_ERG[1]:.0e} erg",
                "statistic": None,
                "computable": False,
                "how_it_is_handled": (
                    "this project does not compute explosion energies and will "
                    "NOT invent a probability for them.  WP0 requires that the "
                    "energy be REPORTED and COMPARED against the inferred "
                    "progenitor mass and type, and that the simulation energy "
                    "is not silently equated with a measured supernova energy.  "
                    "C2 therefore enters the verdict as a stated CONDITIONAL, "
                    "not as a multiplied factor."
                ),
            },
            "C3_progenitor_type": {
                "requirement": (
                    "H25 argues type-Ic explosions are the most likely outcome "
                    "at 3-5 Myr and that their energies may exceed the "
                    "canonical 1e51 erg"
                ),
                "statistic": (
                    f"the fraction of the ledger's supernovae whose progenitor "
                    f"exceeded {STRIPPED_PROGENITOR_MSUN:g} Msun, above which "
                    f"envelope loss makes a stripped-envelope (Ib/c) outcome "
                    f"the expected one"
                ),
                "computable": True,
                "expected_to_be_near_unity": (
                    "the smallest turnoff anywhere on the branch grid is about "
                    "52 Msun, so every star that has died in Cyg OB2 is already "
                    "far above the stripping threshold.  If this term is not "
                    "near 1 the progenitor bookkeeping is wrong."
                ),
            },
            "C4_location": {
                "requirement": (
                    f"a low-density superbubble roughly "
                    f"{SUPERBUBBLE_DIAMETER_PC:g} pc across associated with Cyg "
                    f"OB2, explicitly distinguished from the offset gamma Cygni "
                    f"remnant"
                ),
                "statistic": (
                    "the WP7 in-situ bound: at most 14.6% of the ledger's "
                    "supernovae occurred outside the association, from the "
                    "runaway fraction"
                ),
                "computable": "as a bound, not a distribution",
                "geometry_note": (
                    f"the three subgroups span about 1 deg, roughly 28 pc at "
                    f"1.62 kpc, comfortably inside a "
                    f"{SUPERBUBBLE_DIAMETER_PC:g} pc superbubble.  gamma Cygni "
                    f"sits 2.27 deg away, about 64 pc projected, which is why "
                    f"H25 distinguishes it."
                ),
            },
        },
        "verdict_statistic": {
            "definition": (
                "P_verdict = P(C1) x P(C3) x P(C4), evaluated per branch.  C2 "
                "is carried as a stated conditional and is NOT multiplied in."
            ),
            "why_a_product": (
                "the three terms are near-independent: the age of an event, the "
                "mass of its progenitor and whether it stayed in the "
                "association are set by different physics.  Any residual "
                "correlation is small against the branch spread, and treating "
                "them as independent is the conservative direction for C4, "
                "which is a bound rather than a probability."
            ),
            "reported_both_ways": (
                "marginalized over the headline branch set AND per branch, "
                "never only marginalized"
            ),
        },
        "framing_decision_rule": {
            "source": (
                "paper1_execution_plan.md WP9 step 3, written at planning time "
                "before any result existed"
            ),
            "verbatim": (
                "if the marginalized verdict is stable across branches and its "
                "68% interval does not straddle 'plausible/implausible' -- "
                "Letter framing (lead with the verdict).  Otherwise -- regular "
                "article framing (lead with the ledger + constraints + DR4 "
                "forecast)."
            ),
            "operationalization_fixed_here": {
                "SUPPORTED": (
                    f"P_verdict >= {PLAUSIBLE_THRESHOLD} on EVERY branch in the "
                    f"headline set"
                ),
                "DISFAVOURED": (
                    f"P_verdict <= {IMPLAUSIBLE_THRESHOLD} on EVERY branch in "
                    f"the headline set"
                ),
                "INCONCLUSIVE": (
                    "otherwise -- the branch range straddles a boundary"
                ),
                "framing": (
                    "Letter if SUPPORTED or DISFAVOURED, because either is a "
                    "clean verdict worth leading with.  Regular article if "
                    "INCONCLUSIVE."
                ),
            },
            "threshold_justification": (
                f"{PLAUSIBLE_THRESHOLD} is 'more likely than not' and "
                f"{IMPLAUSIBLE_THRESHOLD} is 'strongly disfavoured'.  Both are "
                f"conventional fixed points, not values tuned to the answer."
            ),
            "binding": (
                "the rule is applied mechanically to whatever the numbers are.  "
                "If the outcome is INCONCLUSIVE the paper is a regular article, "
                "and that is not a failure -- it is what the evidence supports."
            ),
        },
        "branch_set": {
            "headline_alphas": list(HEADLINE_ALPHAS),
            "excluded_alpha": EXCLUDED_ALPHA,
            "exclusion_justification": (
                "the WP5 calibration-window evidence alone rejects alpha = 2.6: "
                "best in 1 of 18 cells and worst median chi-square.  That "
                "evidence is internal to WP5 and free, so the exclusion does "
                "NOT spend the census closure, which remains available as the "
                "out-of-sample validation of the IMF.  alpha = 2.6 is REPORTED "
                "in the sensitivity table as tested-and-excluded, never deleted."
            ),
            "primary_explodability": PRIMARY_EXPLODABILITY,
            "explodability_justification": (
                "WP8: PSR J2032+4127 is a neutron star inside the association, "
                "and neutron stars require successful explosions.  The islands "
                "branch predicts exactly zero supernovae and is therefore "
                "observationally excluded.  It is REPORTED in the sensitivity "
                "table with that reason attached."
            ),
            "families_R_V_duration": "all carried, none excluded",
        },
        "sensitivity_table_required": {
            "axes": [
                "age treatment", "IMF slope", "explodability",
                "isochrone family",
            ],
            "requirement": (
                "identify the DOMINANT driver of the verdict, which need not be "
                "the dominant driver of N_SN"
            ),
            "note": (
                "N_SN's spread is 91.9% alpha, but the VERDICT depends on a "
                "probability that saturates: once the rate is high enough that "
                "P(C1) approaches 1, further increases in N_SN stop moving the "
                "verdict.  The dominant driver may therefore differ, and this "
                "is flagged in advance as something to look for rather than "
                "discovered afterwards."
            ),
        },
        "predictions": [
            {
                "id": "V1",
                "statement": (
                    "C3, the progenitor-type term, exceeds 0.95 on every "
                    "branch, because every star that has died in Cyg OB2 is "
                    "above 52 Msun and therefore far above the "
                    f"{STRIPPED_PROGENITOR_MSUN:g} Msun stripping threshold"
                ),
                "falsifies_if": "C3 < 0.95 on any headline branch",
                "if_falsified": (
                    "the progenitor mass bookkeeping disagrees with the turnoff "
                    "argument and must be reconciled before the verdict stands"
                ),
            },
            {
                "id": "V2",
                "statement": (
                    "the verdict is driven mainly by C1, the age term, and not "
                    "by C3 or C4, which are both near their ceilings"
                ),
                "falsifies_if": (
                    "C3 or C4 contributes more variation across branches than C1"
                ),
                "if_falsified": (
                    "the verdict is being set by a term this project measures "
                    "less well than the age, which weakens the result and must "
                    "be said plainly"
                ),
            },
            {
                "id": "V3",
                "statement": (
                    "the outcome is INCONCLUSIVE under the framing rule, so the "
                    "paper is a regular article.  The reasoning: P(C1) over the "
                    "10-100 kyr window is near 0.5 at the baseline rate, and "
                    "the alpha = 2.0 branches roughly triple the rate while "
                    "CygOB2-C contributes 0 or 7 supernovae depending on "
                    "isochrone family, so the branch range should straddle the "
                    "0.5 boundary"
                ),
                "falsifies_if": (
                    "the verdict is SUPPORTED or DISFAVOURED on every headline "
                    "branch"
                ),
                "if_falsified": (
                    "a Letter is justified and the framing rule says to write "
                    "one -- a better outcome than expected, and the rule is "
                    "followed either way"
                ),
            },
        ],
        "devils_advocate_required": (
            "the gate requires the strongest referee objection to be written "
            "and answered.  It is to be written as the objection a hostile "
            "expert would actually raise, not a soft one that is easy to "
            "rebut, and the answer must concede whatever cannot be defended."
        ),
        "dr4_forecast_required": (
            "quote what Gaia DR4 improves -- age precision from better "
            "astrometry and photometry, and 3D runaway traceback once radial "
            "velocities are available -- and the expected sharpening of the "
            "verdict, so the paper says what would settle what it cannot."
        ),
        "what_wp9_cannot_do": [
            "it cannot evaluate C2, the explosion energy, which this project "
            "does not compute; C2 is a stated conditional",
            "it cannot bound the alternative-source probability quantitatively "
            "-- WP8's neighbouring-association check is explicitly COARSE",
            "it cannot resolve whether gamma Cygni is associated with Cyg OB2",
            "it cannot distinguish high-mass explodability from binary "
            "stripping as the origin of the pulsar",
        ],
    }
    w.write_json(w.PROVENANCE / "wp9_verdict_prereg.json", record)

    print("WP9 — the verdict — pre-registered\n")
    print("  P_verdict = P(C1 age) x P(C3 type) x P(C4 in-situ)")
    print("  C2 energy is a stated CONDITIONAL, never multiplied in\n")
    print(f"  age window        {AGE_WINDOW_KYR[0]:g}-{AGE_WINDOW_KYR[1]:g} kyr"
          f"   (permissive to {AGE_PERMISSIVE_KYR[1]:g} kyr)")
    print(f"  headline alphas   {list(HEADLINE_ALPHAS)}   "
          f"(alpha={EXCLUDED_ALPHA} reported as excluded)")
    print(f"  explodability     {PRIMARY_EXPLODABILITY}   "
          f"(islands excluded by the pulsar)")
    print(f"\n  SUPPORTED    if P >= {PLAUSIBLE_THRESHOLD} on every branch")
    print(f"  DISFAVOURED  if P <= {IMPLAUSIBLE_THRESHOLD} on every branch")
    print("  INCONCLUSIVE otherwise -> regular article, applied mechanically")
    print("\n  predictions:")
    for entry in record["predictions"]:
        print(f"    {entry['id']}  {entry['statement'][:64]}")
    print("\nwrote provenance/wp9_verdict_prereg.json")


if __name__ == "__main__":
    main()
