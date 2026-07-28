#!/usr/bin/env python3
"""Issue #15: pre-register the mass-dependent multiplicity test.

Written BEFORE any multiplicity injection is generated.  WP6's closure test
found 45% more living massive stars than the WP5 normalization predicts, and
the excess could be read as a shallower-than-Salpeter IMF.  Before that reading
is allowed anywhere near the paper, the leading instrumental alternative has to
be measured: the injection truth model applies a CONSTANT binary fraction of
0.40 at every mass, and massive stars are far more multiple than that.

WHY THIS BIASES THE CENSUS UPWARD
---------------------------------
If real massive stars carry more unresolved companions than the truth model
assumes, they are brighter than the model predicts at fixed true mass, so their
inferred masses are biased high, so more of them cross 8 Msun.  That is an
apparent excess with no IMF change at all.  Gaia's angular resolution at
1.62 kpc is about 1000 AU, so essentially every spectroscopic binary and most
visual pairs blend into a single source.

MEASURED VALUES ADOPTED (not guessed)
-------------------------------------
Sana et al. 2012            close-binary fraction f_b = 0.70 for 15-60 Msun;
                            mass-ratio exponent kappa = -0.1 +- 0.6, i.e.
                            consistent with the uniform q the code already uses
Duchene & Kraus 2013        multiplicity >90% for O types, ~50% solar type
Caballero-Nieves et al.     Cyg OB2 ITSELF: 47% of 74 O/early-B stars have a
2020                        resolved companion at 0.08-10 arcsec; 48 of 74 are
                            multiple once spectroscopic binaries are included

THE EXPERIMENT
--------------
Only the TRUTH side changes.  The recovery side keeps its current assumption,
because that is the real-world situation: nature makes binaries at the true
rate and our estimator assumes 0.40.  Injecting truth at the measured rate and
recovering with the unchanged estimator reproduces exactly the bias under test.

Output: provenance/wp6_multiplicity_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_multiplicity_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

# Mass-dependent truth binary fraction.  Anchored at the current 0.40 for
# 2-8 Msun, where it is roughly right for B/A stars and where the WP5
# calibration lives, rising to the measured 0.70 for O stars and held flat
# above.  Log-linear in mass between the anchors: no new free parameter is
# tuned, both endpoints are measured values.
MULTIPLICITY_ANCHORS_MSUN = np.array([8.0, 16.0])
MULTIPLICITY_ANCHORS_FBIN = np.array([0.40, 0.70])
BASELINE_FBIN = 0.40


def truth_binary_fraction(mass: np.ndarray) -> np.ndarray:
    """f_bin(M): 0.40 at and below 8 Msun, 0.70 at and above 16 Msun."""
    mass = np.asarray(mass, dtype=float)
    return np.interp(
        np.log10(mass),
        np.log10(MULTIPLICITY_ANCHORS_MSUN),
        MULTIPLICITY_ANCHORS_FBIN,
        left=MULTIPLICITY_ANCHORS_FBIN[0],
        right=MULTIPLICITY_ANCHORS_FBIN[-1],
    )


def main() -> None:
    closure = pd.read_csv(w.TABLES / "wp6_closure.csv")
    baseline = closure[
        closure.family.eq("PARSEC") & closure.R_V.eq(3.1) & closure.alpha.eq(2.3)
    ]
    grid = np.array([2.0, 4.0, 8.0, 10.0, 12.0, 16.0, 20.0, 40.0, 80.0])

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_multiplicity_prereg.py",
        "status": "PREREGISTERED",
        "issue": "#15 — mass-dependent multiplicity in the injection truth model",
        "question": (
            "How much of WP6's 45% closure excess is produced by the truth "
            "model's constant binary fraction of 0.40, rather than by the IMF?"
        ),
        "why_it_matters": (
            "If multiplicity explains most of the excess, then WP6's "
            "shallower-than-Salpeter signal is an artefact and must not be "
            "reported as an IMF measurement.  If it explains little, the IMF "
            "reading survives its strongest instrumental challenge."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "literature_basis": [
            {
                "source": "Sana et al. 2012",
                "value": "close-binary fraction f_b = 0.70 for 15-60 Msun",
                "use": "the high-mass anchor",
            },
            {
                "source": "Sana et al. 2012",
                "value": "mass-ratio exponent kappa = -0.1 +- 0.6",
                "use": (
                    "consistent with the Uniform(0.1, 1) q distribution already "
                    "in the code, so q is NOT changed"
                ),
            },
            {
                "source": "Duchene & Kraus 2013",
                "value": "multiplicity >90% O type, ~50% solar type",
                "use": "confirms the direction and steepness of the mass trend",
            },
            {
                "source": "Caballero-Nieves et al. 2020",
                "value": (
                    "Cyg OB2 itself: 47% of 74 O/early-B stars have a resolved "
                    "companion at 0.08-10 arcsec; 48 of 74 multiple including "
                    "spectroscopic binaries"
                ),
                "use": (
                    "association-specific confirmation that 0.40 is too low at "
                    "the top; the 0.08-0.6 arcsec companions are unresolved by "
                    "Gaia and blend photometrically"
                ),
            },
        ],
        "model": {
            "changed": "truth-side binary fraction only",
            "unchanged": [
                "the mass-ratio distribution, Uniform(0.1, 1), which the "
                "measured kappa already supports",
                "the recovery/estimator side, which keeps f_bin = 0.40 — this "
                "is deliberate and is the point: nature makes binaries at the "
                "true rate while the estimator assumes 0.40, and that mismatch "
                "IS the bias under test",
                "every RNG draw count, so donors, extinction and QMC "
                "realizations are shared with the repair_v6 nodes",
            ],
            "f_bin_anchors": {
                "8_Msun": float(MULTIPLICITY_ANCHORS_FBIN[0]),
                "16_Msun": float(MULTIPLICITY_ANCHORS_FBIN[-1]),
                "interpolation": "linear in log mass, flat outside the anchors",
            },
            "f_bin_grid": {
                f"{m:g}": round(float(truth_binary_fraction(m)), 3) for m in grid
            },
            "no_new_free_parameter": (
                "both anchors are measured values; nothing is tuned to make the "
                "closure ratio move"
            ),
        },
        "scope": {
            "this_is_a_diagnostic_not_a_new_version": (
                "the test measures the SIZE of the multiplicity effect on the "
                "WP6 closure ratio while holding the accepted WP5 repair_v6 "
                "fixed.  A mass-dependent f_bin would also perturb the WP5 "
                "response inside 2-8 Msun through down-scatter from above, so "
                "ADOPTING it requires a full chain re-run (repair_v7).  That is "
                "deliberately not done here: first measure whether the effect "
                "is large enough to matter."
            ),
            "masses_injected": (
                "only M >= 8 Msun, the closure-test window, since the observed "
                "side and the 2-8 Msun calibration are held fixed"
            ),
        },
        "predictions": [
            {
                "id": "M1",
                "statement": (
                    "The multiplicity model raises the predicted observable "
                    "count above 8 Msun, lowering the closure ratio in every "
                    "subgroup."
                ),
                "falsifies_if": "any subgroup's closure ratio rises or is unchanged",
                "if_falsified": (
                    "the mechanism as reasoned is wrong; the excess is not "
                    "multiplicity and the IMF reading strengthens"
                ),
            },
            {
                "id": "M2",
                "statement": (
                    "The effect is mass-dependent, so it shrinks the excess "
                    "MORE in the subgroups with the higher turnoff: the "
                    "reduction should be largest for CygOB2-C and smallest for "
                    "CygOB2-A."
                ),
                "falsifies_if": (
                    "the ordering of the reduction does not follow the turnoff "
                    "ordering A < B < C"
                ),
                "if_falsified": (
                    "the excess and the multiplicity effect have different mass "
                    "dependences, so multiplicity cannot be the whole story "
                    "even if it lowers the ratio"
                ),
            },
            {
                "id": "M3",
                "statement": (
                    "Quantitative, and the decisive one: the multiplicity model "
                    "absorbs at least HALF of the baseline excess, i.e. the "
                    "grid-median closure ratio at alpha = 2.3 falls from 1.444 "
                    "to below 1.222."
                ),
                "falsifies_if": "the grid-median ratio stays at or above 1.222",
                "if_falsified": (
                    "multiplicity is a real but secondary effect.  The IMF "
                    "reading then survives its strongest instrumental "
                    "challenge and the closing slope should be reported, with "
                    "the residual multiplicity correction carried as a "
                    "systematic."
                ),
                "note": (
                    "this threshold is set BEFORE the run and is not a target: "
                    "either outcome is publishable and both are written down "
                    "here with their consequences"
                ),
            },
        ],
        "decision_rule": (
            "If M3 holds, WP6's shallower-than-Salpeter signal is NOT reported "
            "as an IMF measurement; the closure test is reported as consistent "
            "with Salpeter once multiplicity is modelled, and a repair_v7 full "
            "chain re-run is scheduled.  If M3 fails, the IMF reading is "
            "reported with the measured multiplicity correction applied and the "
            "remainder carried as a systematic.  In BOTH cases the "
            "disfavouring of alpha = 2.6 stands, because it is far larger than "
            "the multiplicity effect can plausibly be."
        ),
        "reference_state": {
            "closure_ratio_baseline_alpha2.3": [
                {
                    "subgroup": row.subgroup,
                    "closure_ratio": round(row.closure_ratio, 3),
                    "turnoff_Msun": round(row.turnoff_prior_mean_Msun, 1),
                }
                for row in baseline.itertuples()
            ],
            "grid_median_ratio_alpha2.3": round(
                float(closure[closure.alpha.eq(2.3)].closure_ratio.median()), 3
            ),
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.TABLES / "wp6_closure.csv",
                w.ROOT / "scripts" / "wp5_injections_repair.py",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp6_multiplicity_prereg.json", record)

    print("issue #15 — multiplicity test pre-registered\n")
    print("  truth-side binary fraction f_bin(M):")
    for mass, value in record["model"]["f_bin_grid"].items():
        print(f"    {mass:>4s} Msun  {value:.3f}")
    print(f"\n  reference: grid-median closure ratio at alpha=2.3 = "
          f"{record['reference_state']['grid_median_ratio_alpha2.3']:.3f}")
    print("  M3 threshold (pre-declared): must fall below 1.222 for "
          "multiplicity to be judged the dominant cause")
    print("\nwrote provenance/wp6_multiplicity_prereg.json")


if __name__ == "__main__":
    main()
