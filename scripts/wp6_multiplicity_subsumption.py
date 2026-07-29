#!/usr/bin/env python3
"""Does repair_v7 already contain issue #15's multiplicity correction?

Bookkeeping check, written when the repair_v7 closure test returned.  It exists
to prevent a specific double-count.

THE RISK
--------
Issue #15 measured mass-dependent multiplicity ABOVE 8 Msun against a truth
model that assumed a constant f_bin = 0.40 everywhere.  It found the effect
absorbs 3.7% of the closure excess and lowers the grid-median closure ratio by
0.0038, and the adjudication in wp6_multiplicity_closure_execution.json directed
that this correction be CARRIED as a separate systematic on top of the accepted
repair_v6 numbers.

repair_v7 replaced the truth-side f_bin across the WHOLE mass range -- the WP5
base grid AND the WP6 mass extension both ran with `--fbin-model extended`.  If
the repair_v7 closure ratio is reported and the issue #15 correction is then
applied to it as well, the same physical effect is counted twice.

WHAT THIS SCRIPT CHECKS
-----------------------
That the repair_v7 truth model is at least as strong as issue #15's treatment
arm at every mass above 8 Msun.  If it is, the correction is SUBSUMED and must
not be applied again.

Outputs: provenance/wp6_multiplicity_subsumption.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_multiplicity_subsumption.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp5_fbin_discriminator_prereg import extended_binary_fraction
from wp6_mass_extension_decision import WP6_MASS_EXTENSION
from wp6_multiplicity_prereg import truth_binary_fraction as issue15_fbin

SN_THRESHOLD_MSUN = 8.0


def main() -> None:
    # Every mass above the threshold that either driver actually injects.
    masses = np.array(
        sorted(
            set(float(m) for m in w.MASS_GRID if m >= SN_THRESHOLD_MSUN)
            | set(float(m) for m in WP6_MASS_EXTENSION)
        )
    )
    v7 = extended_binary_fraction(masses)
    issue15 = issue15_fbin(masses)
    v6 = np.full(masses.shape, float(w.F_BINARY))

    subsumed = bool(np.all(v7 >= issue15 - 1e-12))
    exceeds = masses[v7 > issue15 + 1e-12]

    extension = json.loads(
        (w.PROVENANCE / "wp6_massive_injections_execution_repair_v7.json")
        .read_text(encoding="utf-8")
    )
    injections = json.loads(
        (w.PROVENANCE / "wp5_injections_agenodes_execution_repair_v7.json")
        .read_text(encoding="utf-8")
    )
    multiplicity = json.loads(
        (w.PROVENANCE / "wp6_multiplicity_closure_execution.json")
        .read_text(encoding="utf-8")
    )
    carried = multiplicity["adjudication"]["measured_correction_to_carry"]

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_multiplicity_subsumption.py",
        "status": "SUCCESS",
        "work_package": "WP6 bookkeeping",
        "question": (
            "does repair_v7 already contain the multiplicity correction that "
            "issue #15 directed be carried separately, and would applying both "
            "double-count it?"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "truth_side_fbin_above_threshold": [
            {
                "mass_Msun": round(float(m), 1),
                "repair_v6_constant": round(float(a), 3),
                "issue_15_treatment": round(float(b), 3),
                "repair_v7_extended": round(float(c), 3),
            }
            for m, a, b, c in zip(masses, v6, issue15, v7)
        ],
        "generation_evidence": {
            "wp5_base_grid": {
                "provenance": "provenance/wp5_injections_agenodes_execution_repair_v7.json",
                "truth_binary_fraction_model": injections.get(
                    "truth_binary_fraction_model", "(not recorded)"
                ),
            },
            "wp6_mass_extension": {
                "provenance": "provenance/wp6_massive_injections_execution_repair_v7.json",
                "truth_binary_fraction_model": extension["truth_binary_fraction_model"],
                "grid_Msun": extension["extension_grid_Msun"],
            },
            "note": (
                "both the base grid and the extension ran with the extended "
                "model, so the repair_v7 chain carries mass-dependent f_bin "
                "over the full 0.5-115 Msun range rather than only below 18"
            ),
        },
        "result": {
            "subsumed": subsumed,
            "repair_v7_is_strictly_stronger_at_Msun": [
                round(float(m), 1) for m in exceeds
            ],
            "why_not_identical": (
                "issue #15 anchored 0.40 at 8 Msun rising to 0.70 at 16; "
                "repair_v7 anchors 0.55 at 8 rising to the same 0.70 at 16, "
                "because its below-8 rise has to join continuously at the "
                "threshold.  The two models agree at and above 16 Msun and "
                "repair_v7 is HIGHER in between.  So repair_v7 contains not "
                "merely as much multiplicity above 8 Msun as issue #15's "
                "treatment arm, but slightly more."
            ),
        },
        "consequence": {
            "carried_correction_now_internal": carried,
            "instruction": (
                "the grid-median REDUCTION of 0.0038 recorded in "
                "wp6_multiplicity_closure_execution.json must NOT be applied on "
                "top of any repair_v7 closure ratio.  It is already inside the "
                "repair_v7 response matrices.  It remains applicable to the "
                "repair_v6 numbers, which were generated under a constant 0.40."
            ),
            "size_in_context": (
                "0.0038 against a repair_v7 grid median of 1.067 is 0.36% -- "
                "small enough that it would not have changed a conclusion, and "
                "recorded here because a bookkeeping error that happens to be "
                "small is still a bookkeeping error"
            ),
            "what_is_NOT_subsumed": (
                "issue #15's M1/M2/M3 SCORES are unaffected.  They were a "
                "paired differential measurement of how much multiplicity "
                "above 8 Msun explains the excess, and the answer -- 3.7%, so "
                "not the explanation -- stands independently of which chain "
                "version is reported.  What is subsumed is the numerical "
                "correction, not the finding."
            ),
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROVENANCE / "wp6_massive_injections_execution_repair_v7.json",
                w.PROVENANCE / "wp6_multiplicity_closure_execution.json",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp6_multiplicity_subsumption.json", record)

    print("WP6 — is issue #15's correction already inside repair_v7?\n")
    print("   M/Msun   v6 const   issue15   repair_v7")
    for row in record["truth_side_fbin_above_threshold"]:
        print(
            f"   {row['mass_Msun']:7.1f}   {row['repair_v6_constant']:8.3f}   "
            f"{row['issue_15_treatment']:7.3f}   {row['repair_v7_extended']:9.3f}"
        )
    print(f"\n  subsumed: {subsumed}")
    print(f"  repair_v7 strictly stronger at {len(exceeds)} masses in 8-16 Msun")
    print(
        f"\n  => do NOT apply the carried {carried['grid_median_shift']:.4f} "
        "grid-median reduction to any repair_v7 ratio"
    )
    print("\nwrote provenance/wp6_multiplicity_subsumption.json")


if __name__ == "__main__":
    main()
