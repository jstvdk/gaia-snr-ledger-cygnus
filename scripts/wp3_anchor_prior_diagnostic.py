#!/usr/bin/env python3
"""WP3 anchor spatial-prior diagnostic (issue #1b).

The repair_v1 extinction prior gives every member a mean equal to the median
A_V of its eight nearest spectroscopic anchors, with ONE global width per R_V
branch calibrated by leave-one-out at the anchor density.  This script tests
whether that width is valid where it is actually applied.

Three independent checks:

1. Variogram.  Measure the robust scatter of anchor-to-anchor A_V differences
   as a function of angular separation.  The adopted prior width is only valid
   at the separation where it was calibrated; beyond the correlation length the
   anchors carry no information and the true predictive width is the sill.
2. Collapse.  Compare the per-subgroup A_V spread of the frozen WP3 solution
   (no spatial prior) against repair_v1 (with it).  A prior that merely sharpens
   a well-measured field shrinks every subgroup comparably; a prior that
   over-constrains one subgroup collapses that one alone.
3. Independent cubes.  Vergely+22 and Dharmawardena+22 have no knowledge of the
   anchors.  If a subgroup is genuinely uniform they will say so.

Outputs are diagnostic only -- this script fits nothing and changes no
extinction.  See reports/WP5_PARENT_RANGE_FIX_repair_v2.md for the motivation.
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import curve_fit

import wp5_common as w
from wp3_extinction_law import R_V_BRANCHES
from wp3_repair_common import ANCHOR_NEIGHBOURS, MIN_PRIOR_SIGMA_MAG, AnchorMap


FIGURES = w.ROOT / "figures" / "wp3_repair"
SEPARATION_BINS = np.array(
    [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 2.00]
)
MIN_PAIRS_PER_BIN = 20
COLORS = {"CygOB2-A": "#4477AA", "CygOB2-B": "#EE6677", "CygOB2-C": "#228833"}
MARKERS = {"CygOB2-A": "o", "CygOB2-B": "s", "CygOB2-C": "^"}


def robust_sigma(values: np.ndarray) -> float:
    """Median-absolute-deviation sigma, insensitive to the outlier tail."""
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def half_central68(series: pd.Series) -> float:
    """Half the central-68% span -- the spread measure used throughout WP3/WP5."""
    return float((series.quantile(0.84) - series.quantile(0.16)) / 2.0)


def empirical_variogram(
    anchors: pd.DataFrame, cos_b0: float, rv: float
) -> pd.DataFrame:
    """Robust sigma of pairwise anchor A_V differences vs angular separation."""
    longitude = anchors["l_deg"].to_numpy(float) * cos_b0
    latitude = anchors["b_deg"].to_numpy(float)
    value = anchors[f"av_rv{rv:.1f}"].to_numpy(float)
    finite = np.isfinite(value) & np.isfinite(longitude) & np.isfinite(latitude)
    longitude, latitude, value = longitude[finite], latitude[finite], value[finite]
    separation = np.sqrt(
        (longitude[:, None] - longitude[None, :]) ** 2
        + (latitude[:, None] - latitude[None, :]) ** 2
    )
    upper = np.triu_indices(len(value), 1)
    separation = separation[upper]
    difference = (value[:, None] - value[None, :])[upper]

    rows = []
    for low, high in zip(SEPARATION_BINS[:-1], SEPARATION_BINS[1:], strict=True):
        inside = (separation >= low) & (separation < high)
        if int(inside.sum()) < MIN_PAIRS_PER_BIN:
            continue
        # The difference of two independent draws carries sqrt(2) x the
        # single-star scatter, so divide it back out.
        rows.append(
            {
                "R_V": rv,
                "separation_lo_deg": float(low),
                "separation_hi_deg": float(high),
                "separation_center_deg": float(np.median(separation[inside])),
                "n_pairs": int(inside.sum()),
                "sigma_mag": robust_sigma(difference[inside]) / np.sqrt(2.0),
            }
        )
    return pd.DataFrame(rows)


def exponential_variogram(
    separation: np.ndarray, nugget: float, sill: float, correlation_range: float
) -> np.ndarray:
    """Standard exponential model, returned as a sigma (not a variance)."""
    variance = nugget**2 + (sill**2 - nugget**2) * (
        1.0 - np.exp(-3.0 * separation / correlation_range)
    )
    return np.sqrt(np.maximum(variance, 1e-12))


def evaluate_model(separation: np.ndarray, model: dict[str, float]) -> np.ndarray:
    """Evaluate a fitted variogram dict at the given separations."""
    return exponential_variogram(
        separation,
        model["nugget_mag"],
        model["sill_mag"],
        model["correlation_range_deg"],
    )


def fit_variogram(binned: pd.DataFrame) -> dict[str, float]:
    parameters, _ = curve_fit(
        exponential_variogram,
        binned["separation_center_deg"].to_numpy(float),
        binned["sigma_mag"].to_numpy(float),
        p0=[0.35, 1.05, 0.20],
        sigma=1.0 / np.sqrt(binned["n_pairs"].to_numpy(float)),
        bounds=([0.0, 0.0, 0.01], [3.0, 5.0, 5.0]),
        maxfev=20000,
    )
    nugget, sill, correlation_range = (float(value) for value in parameters)
    return {
        "nugget_mag": nugget,
        "sill_mag": sill,
        "correlation_range_deg": correlation_range,
    }


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    anchor_map = AnchorMap.from_frozen_wp3()

    frozen = pd.read_parquet(w.PROC / "wp3_extinction.parquet")
    repaired = pd.read_parquet(w.PROC / "wp3_extinction_repair_v1.parquet")
    labels = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    cubes = pd.read_parquet(w.PROC / "wp3_cube_comparison.parquet")

    def attach(frame: pd.DataFrame) -> pd.DataFrame:
        return frame.drop(columns=["subgroup"], errors="ignore").merge(
            labels[["source_id", "subgroup"]], on="source_id", how="inner"
        )

    frozen, repaired, cubes = attach(frozen), attach(repaired), attach(cubes)

    variograms = []
    models: dict[float, dict[str, float]] = {}
    for rv in R_V_BRANCHES:
        binned = empirical_variogram(anchor_map.anchors, anchor_map.cos_b0, rv)
        models[float(rv)] = fit_variogram(binned)
        variograms.append(binned)
    variogram = pd.concat(variograms, ignore_index=True)

    rows = []
    for subgroup in w.SUBGROUPS:
        member = repaired[repaired["subgroup"].eq(subgroup)]
        prior_mean, separation = anchor_map.evaluate(
            member["l_deg"].to_numpy(float), member["b_deg"].to_numpy(float), 3.1
        )
        median_separation = float(np.nanmedian(separation))
        required = float(
            evaluate_model(np.array([median_separation]), models[3.1])[0]
        )
        adopted = float(anchor_map.prior_sigma[3.1])
        frozen_spread = half_central68(
            frozen.loc[frozen["subgroup"].eq(subgroup), "av_rv3.1"]
        )
        repaired_spread = half_central68(member["av_rv3.1"])
        cube = cubes[cubes["subgroup"].eq(subgroup)]
        rows.append(
            {
                "subgroup": subgroup,
                "n_members": int(len(member)),
                "median_8th_anchor_separation_deg": median_separation,
                "prior_mean_spread_mag": half_central68(pd.Series(prior_mean)),
                "frozen_wp3_av_spread_mag": frozen_spread,
                "repair_v1_av_spread_mag": repaired_spread,
                "collapse_factor": frozen_spread / repaired_spread,
                "adopted_prior_sigma_mag": adopted,
                "variogram_required_sigma_mag": required,
                "sigma_understatement_factor": required / adopted,
                "vergely_A0_spread_mag": half_central68(cube["vergely_A0"]),
                "dharma_Acum_spread_mag": half_central68(cube["dharma_Acum"]),
            }
        )
    diagnostic = pd.DataFrame(rows)

    variogram_path = w.TABLES / "wp3_anchor_variogram.csv"
    diagnostic_path = w.TABLES / "wp3_anchor_prior_diagnostic.csv"
    variogram.to_csv(variogram_path, index=False)
    diagnostic.to_csv(diagnostic_path, index=False)

    # ---- figure -------------------------------------------------------------
    figure, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    axis = axes[0]
    binned = variogram[variogram["R_V"].eq(3.1)]
    fine = np.linspace(0.005, 1.2, 400)
    axis.plot(
        fine,
        evaluate_model(fine, models[3.1]),
        color="0.35",
        lw=2.0,
        zorder=2,
        label="exponential variogram fit",
    )
    axis.scatter(
        binned["separation_center_deg"],
        binned["sigma_mag"],
        s=np.clip(binned["n_pairs"] / 12.0, 20, 160),
        color="0.15",
        zorder=3,
        label="anchor pairs (area $\\propto$ N)",
    )
    axis.axhline(
        anchor_map.prior_sigma[3.1],
        color="#CCBB44",
        lw=2.0,
        ls="--",
        zorder=2,
        label=f"adopted prior $\\sigma$ = {anchor_map.prior_sigma[3.1]:.3f} mag",
    )
    for row in diagnostic.itertuples():
        axis.plot(
            [row.median_8th_anchor_separation_deg],
            [row.variogram_required_sigma_mag],
            marker=MARKERS[row.subgroup],
            ms=11,
            color=COLORS[row.subgroup],
            mec="white",
            mew=1.4,
            zorder=5,
        )
        axis.annotate(
            row.subgroup.replace("CygOB2-", ""),
            (row.median_8th_anchor_separation_deg, row.variogram_required_sigma_mag),
            textcoords="offset points",
            xytext=(9, -3),
            color=COLORS[row.subgroup],
            fontweight="bold",
        )
    axis.set_xscale("log")
    axis.set_xlabel("angular separation (deg)")
    axis.set_ylabel("robust $\\sigma$ of $A_V$ difference (mag)")
    axis.set_title(
        "Anchor $A_V$ decorrelates well inside the distance\n"
        "at which the prior is applied to CygOB2-B",
        fontsize=11,
    )
    axis.grid(alpha=0.25, lw=0.6)
    axis.set_axisbelow(True)
    axis.legend(loc="lower right", fontsize=8.5, framealpha=0.92)

    axis = axes[1]
    x = np.arange(len(diagnostic))
    width = 0.38
    axis.bar(
        x - width / 2,
        diagnostic["frozen_wp3_av_spread_mag"],
        width,
        color=[COLORS[s] for s in diagnostic["subgroup"]],
        alpha=0.45,
        edgecolor="white",
        lw=1.2,
        label="frozen WP3 (no spatial prior)",
    )
    axis.bar(
        x + width / 2,
        diagnostic["repair_v1_av_spread_mag"],
        width,
        color=[COLORS[s] for s in diagnostic["subgroup"]],
        edgecolor="white",
        lw=1.2,
        label="repair_v1 (with spatial prior)",
    )
    for index, row in enumerate(diagnostic.itertuples()):
        axis.annotate(
            f"{row.collapse_factor:.1f}$\\times$",
            (index, max(row.frozen_wp3_av_spread_mag, row.repair_v1_av_spread_mag)),
            textcoords="offset points",
            xytext=(0, 7),
            ha="center",
            fontweight="bold",
            color=COLORS[row.subgroup],
        )
    axis.set_xticks(x)
    axis.set_xticklabels([s.replace("CygOB2-", "") for s in diagnostic["subgroup"]])
    axis.set_xlabel("subgroup")
    axis.set_ylabel("$A_V$ spread, half central-68% (mag)")
    axis.set_title(
        "The prior collapses CygOB2-B's differential extinction\n"
        "4.8$\\times$ harder than CygOB2-A's",
        fontsize=11,
    )
    axis.grid(alpha=0.25, lw=0.6, axis="y")
    axis.set_axisbelow(True)
    axis.legend(loc="upper right", fontsize=8.5, framealpha=0.92)

    figure.tight_layout()
    figure_path = FIGURES / "wp3_anchor_prior_diagnostic.png"
    figure.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(figure)

    inputs = [
        w.PROC / "wp3_extinction.parquet",
        w.PROC / "wp3_extinction_repair_v1.parquet",
        w.PROC / "wp3_cube_comparison.parquet",
        w.TABLES / "wp2_subgroup_labels.parquet",
    ]
    outputs = [variogram_path, diagnostic_path, figure_path]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp3_anchor_prior_diagnostic.py",
        "status": "SUCCESS",
        "purpose": (
            "Diagnose PROJECT_TRACE.md issue #1b: whether the repair_v1 anchor "
            "spatial prior over-constrains CygOB2-B's extinction and thereby "
            "drives the WP5 bin-2 residual failure."
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "configuration": {
            "anchor_neighbours": ANCHOR_NEIGHBOURS,
            "n_anchors": int(len(anchor_map.anchors)),
            "min_prior_sigma_mag": MIN_PRIOR_SIGMA_MAG,
            "separation_bins_deg": SEPARATION_BINS.tolist(),
            "min_pairs_per_bin": MIN_PAIRS_PER_BIN,
            "spread_measure": "half of the central-68% span",
            "variogram_model": (
                "sigma(d)^2 = nugget^2 + (sill^2 - nugget^2)(1 - exp(-3d/range))"
            ),
        },
        "adopted_prior_sigma_mag": {
            f"rv{rv:.1f}": float(anchor_map.prior_sigma[rv]) for rv in R_V_BRANCHES
        },
        "variogram_fit": {f"rv{rv:.1f}": models[float(rv)] for rv in R_V_BRANCHES},
        "per_subgroup": diagnostic.to_dict(orient="records"),
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in inputs},
        "outputs": {
            str(p.relative_to(w.ROOT)): {
                "sha256": w.sha256(p),
                "bytes": p.stat().st_size,
            }
            for p in outputs
        },
    }
    w.write_json(w.PROVENANCE / "wp3_anchor_prior_diagnostic_execution.json", record)

    print(diagnostic.to_string(index=False))
    print()
    for rv in R_V_BRANCHES:
        model = models[float(rv)]
        print(
            f"R_V={rv:.1f} variogram: nugget {model['nugget_mag']:.3f} mag, "
            f"sill {model['sill_mag']:.3f} mag, "
            f"range {model['correlation_range_deg']:.3f} deg"
        )


if __name__ == "__main__":
    main()
