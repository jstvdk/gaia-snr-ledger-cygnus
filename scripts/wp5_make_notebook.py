#!/usr/bin/env python3
"""Create the explanatory WP5 notebook from the frozen diagnostic products."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "notebooks" / "wp5_imf_normalization_and_completeness.ipynb"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def main() -> None:
    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {
        "kernelspec": {
            "display_name": "Python (cygob2-gaia)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
        "wp": "WP5",
        "status": "BLOCKED_AT_VALIDATION_GATE",
    }
    notebook["cells"] = [
        markdown(
            """
# WP5 — IMF normalization and completeness

This notebook is the explanatory, read-only view of the executed WP5 pipeline.
It consumes only named WP5 artifacts.  It does not refit or overwrite upstream
WP2–WP4 products.

The intended downstream quantity is the birth-population normalization
\(k\) per subgroup.  The validation result is blocking: no absolute 95%
end-to-end completeness edge exists, and the response-aware IMF residuals are
not consistent with Poisson scatter.  The notebook therefore ends with a
blocked handoff rather than an accepted normalization.
"""
        ),
        markdown(
            """
## 1. Load the named artifacts and verify their shapes

Purpose: confirm that every mandatory family, extinction-law, slope, and
subgroup branch exists before interpreting any result.  Expected dimensions
are fixed by the registered grid: 2 families × 3 R_V values × 3 subgroups,
with 3 IMF slopes added for the normalization.

Downstream effect: a missing branch would invalidate the sensitivity envelope
and stop this notebook immediately.
"""
        ),
        code(
            """
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Image

ROOT = Path("..").resolve() if Path.cwd().name == "notebooks" else Path(".").resolve()
PROC = ROOT / "data" / "processed"

curves = pd.read_parquet(PROC / "wp5_completeness_curves.parquet")
response = pd.read_parquet(PROC / "wp5_injection_response.parquet")
normalization = pd.read_parquet(PROC / "wp5_imf_normalization.parquet")
mass_bins = pd.read_parquet(PROC / "wp5_mass_function_bins.parquet")
association = pd.read_parquet(PROC / "wp5_association_mass.parquet")

expected = {
    "completeness_rows": 2 * 3 * 3 * 31,
    "response_rows": 2 * 3 * 3 * 31 * 400,
    "normalization_rows": 2 * 3 * 3 * 3,
    "mass_bin_rows": 2 * 3 * 3 * 3 * 6,
    "association_rows": 2 * 3 * 3,
}
actual = {
    "completeness_rows": len(curves),
    "response_rows": len(response),
    "normalization_rows": len(normalization),
    "mass_bin_rows": len(mass_bins),
    "association_rows": len(association),
}
assert actual == expected, (actual, expected)
pd.DataFrame({"expected": expected, "actual": actual})
"""
        ),
        markdown(
            """
## 2. Inspect the end-to-end completeness curve

Purpose: distinguish Gaia-query visibility from the full chain.  A star can be
brighter than G=19 yet fail the astrometric-quality domain, the WP2
cluster-versus-field classifier, or the photometry needed for a WP4 mass.

Threshold: the execution plan requires the lower calibration edge to be at
least 95% complete.  The vertical dotted line marks the nominal 2 M_sun edge;
the horizontal dashed line is the absolute 95% requirement.

Limitation: this is a catalogue-level injection using real Gaia observational
states.  Gaia epoch images and AGIS cannot be rerun, so it is not an
image-level artificial-source test.
"""
        ),
        code(
            """
display(Image(filename=str(ROOT / "figures" / "wp5" / "wp5_completeness_curves.png")))

baseline_curve = curves[(curves.family == "PARSEC") & (curves.R_V == 3.1)]
summary = []
for subgroup, group in baseline_curve.groupby("subgroup"):
    at_two = group.loc[group.primary_mass.eq(2.0)].iloc[0]
    summary.append({
        "subgroup": subgroup,
        "full_recovery_at_2_Msun": at_two.recovery_fraction,
        "median_full_recovery_6_8_Msun":
            group.loc[group.primary_mass.between(6, 8), "recovery_fraction"].median(),
        "maximum_monotone_recovery_above_2":
            group.loc[group.primary_mass.ge(2), "recovery_isotonic"].max(),
    })
pd.DataFrame(summary)
"""
        ),
        markdown(
            """
## 3. Freeze the calibration-window verdict

Purpose: prevent a relative bright-plateau threshold from being silently
reported as absolute completeness.  `absolute_95_edge_Msun` must be finite to
pass the planned criterion.

Result: no subgroup/family/R_V branch has an absolute 95% edge.  The
`corrected_no_absolute95_edge` fit is retained only as a diagnostic because the
measured response can correct incomplete bins; it does not satisfy the original
window requirement.
"""
        ),
        code(
            """
edge_audit = normalization[
    ["subgroup", "family", "R_V", "absolute_95_edge_Msun",
     "relative_95_plateau_edge_Msun_diagnostic",
     "bright_plateau_completeness", "window_status"]
].drop_duplicates()
assert edge_audit.absolute_95_edge_Msun.isna().all()
edge_audit.head(18)
"""
        ),
        markdown(
            """
## 4. Inspect the response-aware Poisson IMF fit

Purpose: compare the membership-weighted observed counts with a forward model
that includes selection, binary brightening, extinction/photometric errors, and
migration through the exact WP4 nearest-isochrone mass estimator.

Gate: Pearson residuals must be compatible with Poisson scatter, with no
systematic bin behavior.  The shaded ±3 region is a visual diagnostic; the
registered numerical gate also requires chi-square p>=0.01 and no significant
monotonic residual trend.
"""
        ),
        code(
            """
display(Image(filename=str(ROOT / "figures" / "wp5" / "wp5_mass_function.png")))

baseline = normalization[
    (normalization.family == "PARSEC")
    & (normalization.R_V == 3.1)
    & (normalization.alpha == 2.3)
]
baseline[
    ["subgroup", "k_median", "k_lo68", "k_hi68",
     "poisson_chi_square_p", "max_abs_pearson_residual",
     "residual_gate_pass"]
]
"""
        ),
        markdown(
            """
## 5. Localize the residual failure

Purpose: show that the failure is not a vague fit-quality label.  In all three
subgroups the recovered 2.0–2.52 M_sun bin is substantially more populated
than the forward response predicts, while adjacent bins are deficient.

Downstream effect: changing only k cannot repair the shape.  The failure
persists through all mandatory IMF slopes, both isochrone families, and all
R_V branches.
"""
        ),
        code(
            """
baseline_bins = mass_bins[
    (mass_bins.family == "PARSEC")
    & (mass_bins.R_V == 3.1)
    & (mass_bins.alpha == 2.3)
]
first_bin = baseline_bins[baseline_bins.bin_index == 0][
    ["subgroup", "mass_lo", "mass_hi", "membership_weighted_count",
     "expected_count_at_k_median", "pearson_residual"]
]
display(first_bin)

gate_counts = normalization.residual_gate_pass.value_counts(dropna=False)
print("Residual-gate outcomes across all 54 fits:")
display(gate_counts)
assert not normalization.residual_gate_pass.any()
"""
        ),
        markdown(
            """
## 6. Association-mass sanity check

Purpose: integrate the diagnostic k values over a continuous Kroupa-like IMF
(slope 1.3 below 0.5 M_sun, branch slope above) and add unresolved companions
under f_bin=0.4, q~U[0.1,1].

Interpretation: the baseline mass lies within a factor two of the Wright+15
16,500 M_sun scale.  This is a useful order-of-magnitude check, but it cannot
override a failed mass-function shape: a biased distribution can still have a
plausible integral.
"""
        ),
        code(
            """
display(Image(filename=str(ROOT / "figures" / "wp5" / "wp5_association_mass.png")))

baseline_mass = association[
    (association.family == "PARSEC")
    & (association.R_V == 3.1)
    & (association.alpha == 2.3)
]
baseline_mass.T
"""
        ),
        markdown(
            """
## 7. Blocking handoff

WP5 is not accepted and WP6/WP7 must not consume the diagnostic normalization.
The evidence calls for a scoped revision of the WP4 2–5 M_sun mass inference:
replace nearest-point masses with per-star mass posteriors or fit the IMF
directly in CMD space, retaining binaries and extinction covariance.  The
current injection catalogue and response machinery should be preserved and
rerun against that revised estimator.

This conclusion is deliberately stronger than “large uncertainty”: the
registered residual gate failed in every branch.
"""
        ),
        code(
            """
fit_log = json.loads(
    (ROOT / "provenance" / "wp5_imf_fit_execution.json").read_text()
)
gate = pd.Series(fit_log["gate"], name="result")
display(gate)
assert not gate["absolute_95_edge_exists_all_branches"]
assert not gate["baseline_all_subgroups_residual_gate"]
print("FINAL WP5 VERDICT: BLOCKED; no downstream authority.")
"""
        ),
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(notebook, OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()

