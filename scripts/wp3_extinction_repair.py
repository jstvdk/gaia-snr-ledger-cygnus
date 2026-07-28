#!/usr/bin/env python3
"""Run the versioned WP3 extinction repair without overwriting frozen WP3."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn

from wp3_common import BANDS, DIST_MODULUS
from wp3_extinction_law import R_V_BRANCHES, band_coefficients
from wp3_repair_common import (
    ANCHOR_NEIGHBOURS,
    AV_GRID,
    PHOTOMETRIC_FLOOR_MAG,
    PROC,
    REPAIR_VERSION,
    ROOT,
    AnchorMap,
    fit_extinction_posterior,
    gaussian_grid,
    load_template_library,
    _grid_summary,
)


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    frozen = pd.read_parquet(PROC / "wp3_extinction.parquet")
    photometry = pd.read_parquet(PROC / "wp3_member_photometry.parquet")
    photometry = frozen[["source_id", "av_method"]].merge(
        photometry, on="source_id", how="left", validate="one_to_one"
    )
    templates, template_magnitudes, template_weights = load_template_library()
    anchor_map = AnchorMap.from_frozen_wp3()
    anchor_ids = set(anchor_map.anchors["source_id"].astype("int64"))
    anchor_photometry = photometry[photometry["source_id"].isin(anchor_ids)]
    frozen_index = frozen.set_index("source_id")
    branch_calibration: dict[float, dict[str, float]] = {}
    for rv in R_V_BRANCHES:
        residuals = []
        for _, anchor_row in anchor_photometry.iterrows():
            magnitudes = np.array(
                [anchor_row.get(band, np.nan) for band in BANDS], dtype=float
            )
            errors = np.array(
                [anchor_row.get(f"{band}_err", np.nan) for band in BANDS], dtype=float
            )
            _, summary = fit_extinction_posterior(
                magnitudes,
                errors,
                float(rv),
                7.0,
                100.0,
                template_magnitudes,
                template_weights,
                template_branch_sigma=0.0,
            )
            estimate = summary.get("av_q50", np.nan)
            truth = frozen_index.loc[
                int(anchor_row["source_id"]), f"av_rv{rv:.1f}"
            ]
            if np.isfinite(estimate) and np.isfinite(truth):
                residuals.append(float(estimate - truth))
        residual = np.asarray(residuals)
        centre = float(np.median(residual))
        rms_about_median = float(np.sqrt(np.mean((residual - centre) ** 2)))
        p16 = float(np.quantile(residual, 0.16))
        p84 = float(np.quantile(residual, 0.84))
        asymmetric_central68_extent = max(
            abs(p16 - centre), abs(p84 - centre)
        )
        adopted_branch_sigma = max(
            rms_about_median, asymmetric_central68_extent
        )
        branch_calibration[float(rv)] = {
            "n_anchors": int(len(residual)),
            "broadband_minus_spectroscopic_median_mag": centre,
            "rms_about_median_mag": rms_about_median,
            "p16_mag": p16,
            "p84_mag": p84,
            "asymmetric_central68_extent_mag": asymmetric_central68_extent,
            "adopted_template_branch_sigma_mag": adopted_branch_sigma,
            "adoption_rule": (
                "max(RMS about median, larger asymmetric central-68 endpoint); "
                "prevents a minority branch from being hidden by the MAD/RMS summary"
            ),
        }

    n_star = len(photometry)
    posterior_cube = np.full(
        (n_star, len(R_V_BRANCHES), len(AV_GRID)), np.nan, dtype=np.float32
    )
    rows: list[dict] = []

    for star_index, row in photometry.iterrows():
        source_id = int(row["source_id"])
        magnitudes = np.array([row.get(band, np.nan) for band in BANDS], dtype=float)
        errors = np.array([row.get(f"{band}_err", np.nan) for band in BANDS], dtype=float)
        record: dict[str, object] = {
            "source_id": source_id,
            "av_method_frozen": row["av_method"],
        }
        for rv_index, rv in enumerate(R_V_BRANCHES):
            centre, separation = anchor_map.evaluate(
                np.array([row.get("l_deg", np.nan)]),
                np.array([row.get("b_deg", np.nan)]),
                float(rv),
            )
            prior_mean = float(centre[0])
            prior_sigma = float(
                anchor_map.prior_sigma_at(separation, float(rv))[0]
            )
            method = str(row["av_method"])
            if method == "intrinsic_color_spectroscopic":
                frozen_row = frozen.loc[frozen["source_id"].eq(source_id)].iloc[0]
                av = float(frozen_row[f"av_rv{rv:.1f}"])
                av_error = float(frozen_row[f"av_err_rv{rv:.1f}"])
                posterior = gaussian_grid(av, av_error)
                summary = _grid_summary(posterior)
                summary.update(
                    {
                        "n_bands_fit": 0,
                        "conditional_av_sigma": av_error,
                        "chi2_min": np.nan,
                        "effective_error_floor_mag": np.nan,
                        "prior_mean": prior_mean,
                        "prior_sigma": prior_sigma,
                    }
                )
                estimator = "spectroscopic_anchor_frozen"
            elif method == "broadband_multiband":
                posterior, summary = fit_extinction_posterior(
                    magnitudes,
                    errors,
                    float(rv),
                    prior_mean,
                    prior_sigma,
                    template_magnitudes,
                    template_weights,
                    template_branch_sigma=branch_calibration[float(rv)][
                        "adopted_template_branch_sigma_mag"
                    ],
                )
                estimator = "broadband_floor_spatial_posterior"
            else:
                posterior = None
                summary = {"n_bands_fit": int(np.isfinite(magnitudes).sum())}
                estimator = "frozen_special_case"

            tag = f"rv{rv:.1f}"
            record[f"local_anchor_av_{tag}"] = prior_mean
            record[f"anchor_eighth_neighbour_deg_{tag}"] = float(separation[0])
            record[f"anchor_prior_sigma_{tag}"] = prior_sigma
            record[f"estimator_{tag}"] = estimator
            if posterior is not None:
                posterior_cube[star_index, rv_index] = posterior.astype(np.float32)
                for key, value in summary.items():
                    record[f"{key}_{tag}"] = value
            else:
                for key in [
                    "av_mean", "av_sd", "av_q16", "av_q50", "av_q84", "av_map",
                    "n_modes", "mode1_av", "mode2_av", "mode2_relative_height",
                    "chi2_min", "conditional_av_sigma",
                ]:
                    record[f"{key}_{tag}"] = np.nan
                record[f"n_bands_fit_{tag}"] = summary["n_bands_fit"]
        rows.append(record)

    summary_frame = pd.DataFrame(rows)
    summary_path = PROC / f"wp3_extinction_posterior_summary_{REPAIR_VERSION}.parquet"
    summary_frame.to_parquet(summary_path, index=False)
    posterior_path = PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz"
    np.savez_compressed(
        posterior_path,
        source_id=photometry["source_id"].to_numpy("int64"),
        rv=np.asarray(R_V_BRANCHES, dtype=float),
        av_grid=AV_GRID,
        probability=posterior_cube,
    )

    repaired = frozen.copy()
    repaired["av_method_pre_repair"] = repaired["av_method"]
    repaired["av_model_version"] = "frozen_special_case"
    broad = repaired["av_method"].eq("broadband_multiband")
    repaired.loc[broad, "av_model_version"] = "floor_spatial_full_posterior_repair_v1"
    summary_index = summary_frame.set_index("source_id")
    for rv in R_V_BRANCHES:
        tag = f"rv{rv:.1f}"
        replacement = repaired["source_id"].map(summary_index[f"av_q50_{tag}"])
        q16 = repaired["source_id"].map(summary_index[f"av_q16_{tag}"])
        q84 = repaired["source_id"].map(summary_index[f"av_q84_{tag}"])
        uncertainty = 0.5 * (q84 - q16)
        repaired.loc[broad, f"av_rv{rv:.1f}"] = replacement[broad]
        repaired.loc[broad, f"av_err_rv{rv:.1f}"] = uncertainty[broad]
        coefficient = band_coefficients(float(rv))
        for band in BANDS:
            repaired.loc[broad, f"{band}0_abs_rv{rv:.1f}"] = (
                repaired.loc[broad, band]
                - DIST_MODULUS
                - coefficient[band] * repaired.loc[broad, f"av_rv{rv:.1f}"]
            )
        repaired.loc[broad, f"A_G_rv{rv:.1f}"] = (
            coefficient["G"] * repaired.loc[broad, f"av_rv{rv:.1f}"]
        )
        repaired.loc[broad, f"BPRP0_rv{rv:.1f}"] = (
            repaired.loc[broad, f"BP0_abs_rv{rv:.1f}"]
            - repaired.loc[broad, f"RP0_abs_rv{rv:.1f}"]
        )
        repaired.loc[broad, f"GKs0_rv{rv:.1f}"] = (
            repaired.loc[broad, f"G0_abs_rv{rv:.1f}"]
            - repaired.loc[broad, f"Ks0_abs_rv{rv:.1f}"]
        )
    repaired_path = PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet"
    repaired.to_parquet(repaired_path, index=False)

    finite_broad = broad & repaired["av_rv3.1"].notna()
    multimodal = summary_frame.loc[
        summary_frame["av_method_frozen"].eq("broadband_multiband"),
        "n_modes_rv3.1",
    ]
    provenance = {
        "repair_version": REPAIR_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "trigger": {
            "brief": "tasks/wp3_extinction_repair_brief.md",
            "preserved_blocking_manifest": "provenance/wp5_manifest.json",
            "reason": "frozen WP3 optical-dominated extinction degeneracy caused WP4 mass migration",
        },
        "script": "scripts/wp3_extinction_repair.py",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "inputs": {
            str(path): sha256(ROOT / path)
            for path in [
                "data/processed/wp3_extinction.parquet",
                "data/processed/wp3_member_photometry.parquet",
                "data/processed/wp3_isochrones_parsec.parquet",
                "data/processed/wp3_isochrones_mist.parquet",
                "tasks/wp3_extinction_repair_brief.md",
            ]
        },
        "configuration": {
            "all_band_error_floor_mag": PHOTOMETRIC_FLOOR_MAG,
            "anchor_neighbours": ANCHOR_NEIGHBOURS,
            "anchor_prior_width_method": "leave-one-out robust sigma against median of 8 nearest other anchors",
            "anchor_variogram": {
                f"rv{rv:.1f}": anchor_map.variogram[float(rv)]
                for rv in R_V_BRANCHES
            },
            "template_branch_uncertainty_calibration": {
                f"rv{rv:.1f}": branch_calibration[float(rv)]
                for rv in R_V_BRANCHES
            },
            "av_grid_mag": {
                "minimum": float(AV_GRID.min()),
                "maximum": float(AV_GRID.max()),
                "step": float(AV_GRID[1] - AV_GRID[0]),
            },
            "template_prior": "equal PARSEC/MIST and age cells; d(log initial mass) quadrature within cell",
        },
        "counts": {
            "members": int(len(repaired)),
            "spectroscopic_anchors": int(repaired["av_method"].eq("intrinsic_color_spectroscopic").sum()),
            "broadband_rows": int(broad.sum()),
            "broadband_rows_with_repaired_rv31": int(finite_broad.sum()),
            "broadband_multimodal_rv31": int((multimodal.fillna(0) > 1).sum()),
            "templates": int(len(templates)),
        },
        "outputs": {
            str(path.relative_to(ROOT)): {
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in [summary_path, posterior_path, repaired_path]
        },
        "frozen_outputs_overwritten": False,
    }
    write_json(ROOT / "provenance" / "wp3_repair_execution.json", provenance)
    print(f"wrote {repaired_path.relative_to(ROOT)} ({len(repaired)} rows)")
    print(
        f"repaired broadband R_V=3.1: {int(finite_broad.sum())}/{int(broad.sum())}; "
        f"multimodal: {int((multimodal.fillna(0) > 1).sum())}"
    )
    print(
        "anchor prior sigmas: "
        + ", ".join(
            f"R_V={rv:.1f}: {anchor_map.prior_sigma[float(rv)]:.3f} mag"
            for rv in R_V_BRANCHES
        )
    )


if __name__ == "__main__":
    main()
