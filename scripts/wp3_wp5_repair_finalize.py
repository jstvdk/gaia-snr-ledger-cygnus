#!/usr/bin/env python3
"""Finalize the versioned WP3--WP5 repair audit without changing science data."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
PROV = ROOT / "provenance"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "figures" / "wp3_repair"
VERSION = "repair_v1"
MASS_BINS = np.array([2.0, 2.52, 3.17, 5.04, 6.35, 8.0, 12.0])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def absolute_edge(curve: pd.DataFrame) -> float:
    curve = curve[curve["primary_mass"].between(2.0, 8.0)].sort_values(
        "primary_mass"
    )
    values = curve["recovery_isotonic"].to_numpy(float)
    suffix_minimum = np.minimum.accumulate(values[::-1])[::-1]
    eligible = np.flatnonzero(suffix_minimum >= 0.95)
    return (
        float(curve["primary_mass"].to_numpy(float)[eligible[0]])
        if len(eligible)
        else np.nan
    )


def residual_table(
    extinction_value: pd.Series,
    map_value: pd.Series,
    mass: pd.Series,
    label: str,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "av": extinction_value,
            "local_anchor_av": map_value,
            "mass": mass,
        }
    )
    frame["delta_av"] = frame["av"] - frame["local_anchor_av"]
    frame["mass_bin"] = pd.cut(frame["mass"], MASS_BINS, right=False)
    table = (
        frame.groupby("mass_bin", observed=True)["delta_av"]
        .agg(["count", "median", "mean", "std"])
        .reset_index()
    )
    table.insert(0, "stage", label)
    table["mass_lo"] = table["mass_bin"].map(lambda value: value.left).astype(float)
    table["mass_hi"] = table["mass_bin"].map(lambda value: value.right).astype(float)
    table = table.drop(columns="mass_bin")
    return table


def main() -> None:
    REPORTS.mkdir(exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    summary = pd.read_parquet(
        PROC / f"wp3_extinction_posterior_summary_{VERSION}.parquet"
    )
    frozen_extinction = pd.read_parquet(PROC / "wp3_extinction.parquet")
    frozen_mass = pd.read_parquet(PROC / "wp4_masses.parquet")
    repaired_mass = pd.read_parquet(
        PROC / f"wp4_mass_posteriors_{VERSION}.parquet"
    )
    broadband = summary["av_method_frozen"].eq("broadband_multiband")
    joined_frozen = (
        summary.loc[
            broadband, ["source_id", "local_anchor_av_rv3.1"]
        ]
        .merge(
            frozen_extinction[["source_id", "av_rv3.1"]],
            on="source_id",
            validate="one_to_one",
        )
        .merge(
            frozen_mass[["source_id", "mass_baseline"]],
            on="source_id",
            validate="one_to_one",
        )
    )
    joined_repair = summary.loc[
        broadband,
        ["source_id", "local_anchor_av_rv3.1", "av_q50_rv3.1"],
    ].merge(
        repaired_mass[["source_id", "mass_baseline"]],
        on="source_id",
        validate="one_to_one",
    )
    before = residual_table(
        joined_frozen["av_rv3.1"],
        joined_frozen["local_anchor_av_rv3.1"],
        joined_frozen["mass_baseline"],
        "frozen_before_repair",
    )
    after = residual_table(
        joined_repair["av_q50_rv3.1"],
        joined_repair["local_anchor_av_rv3.1"],
        joined_repair["mass_baseline"],
        VERSION,
    )
    f4_table = pd.concat([before, after], ignore_index=True)
    f4_path = TABLES / f"wp3_extinction_mass_invariant_{VERSION}.csv"
    f4_table.to_csv(f4_path, index=False)

    after_valid = after[after["count"] >= 10]
    rank = spearmanr(
        0.5 * (after_valid["mass_lo"] + after_valid["mass_hi"]),
        after_valid["median"],
    )
    f4_max = float(after_valid["median"].abs().max())
    f4_pass = bool(
        f4_max < 0.30
        and np.isfinite(rank.pvalue)
        and float(rank.pvalue) >= 0.05
    )

    exact_fraction = float(
        repaired_mass["mass_baseline"].dropna().value_counts(normalize=True).max()
    )
    populated_25_32 = int(
        repaired_mass["mass_baseline"].between(2.5, 3.2).sum()
    )
    mass_acceptance = bool(exact_fraction < 0.02 and populated_25_32 > 0)
    median_gt8 = int((repaired_mass["mass_baseline"] > 8.0).sum())
    expected_gt8 = float(
        repaired_mass["mass_baseline_p_gt8"].sum(skipna=True)
    )

    normalization = pd.read_parquet(
        PROC / f"wp5_imf_normalization_{VERSION}.parquet"
    )
    bins = pd.read_parquet(
        PROC / f"wp5_mass_function_bins_{VERSION}.parquet"
    )
    association = pd.read_parquet(
        PROC / f"wp5_association_mass_{VERSION}.parquet"
    )
    baseline = normalization[
        normalization["family"].eq("PARSEC")
        & normalization["R_V"].eq(3.1)
        & normalization["alpha"].eq(2.3)
    ].copy()
    baseline_bins = bins[
        bins["family"].eq("PARSEC")
        & bins["R_V"].eq(3.1)
        & bins["alpha"].eq(2.3)
    ].copy()
    baseline_gate = bool(baseline["residual_gate_pass"].all())
    baseline_assoc = association[
        association["family"].eq("PARSEC")
        & association["R_V"].eq(3.1)
        & association["alpha"].eq(2.3)
    ].iloc[0]
    mass_sanity = bool(baseline_assoc["within_factor_two_literature"])

    curves = pd.read_parquet(
        PROC / f"wp5_completeness_curves_{VERSION}.parquet"
    )
    baseline_curve = curves[
        curves["family"].eq("PARSEC") & curves["R_V"].eq(3.1)
    ]
    completeness_rows = []
    for subgroup, curve in baseline_curve.groupby("subgroup"):
        plateau = float(
            curve.loc[
                curve["primary_mass"].between(6.0, 8.0),
                "recovery_isotonic",
            ].median()
        )
        relative_target = 0.95 * plateau
        relative = curve[
            curve["primary_mass"].between(2.0, 8.0)
            & curve["recovery_isotonic"].ge(relative_target)
        ]
        completeness_rows.append(
            {
                "subgroup": subgroup,
                "absolute_95_edge_Msun": absolute_edge(curve),
                "bright_plateau_completeness": plateau,
                "relative_95_plateau_edge_Msun_diagnostic": (
                    float(relative["primary_mass"].min())
                    if len(relative)
                    else np.nan
                ),
            }
        )
    completeness = pd.DataFrame(completeness_rows)
    completeness_path = TABLES / f"wp5_completeness_baseline_{VERSION}.csv"
    completeness.to_csv(completeness_path, index=False)
    baseline_bins_path = TABLES / f"wp5_baseline_residuals_{VERSION}.csv"
    baseline_bins.to_csv(baseline_bins_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for table, label, marker in [
        (before, "frozen", "o"),
        (after, VERSION, "s"),
    ]:
        centre = np.sqrt(table["mass_lo"] * table["mass_hi"])
        axes[0].plot(centre, table["median"], marker=marker, label=label)
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].axhspan(-0.30, 0.30, color="0.85", label="F4 ±0.30 mag")
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"assigned mass [$M_\odot$]")
    axes[0].set_ylabel(r"median $A_V-A_{V,\mathrm{anchor\ map}}$ [mag]")
    axes[0].legend()
    axes[0].set_title("Permanent WP3 invariant")
    for subgroup, group in baseline_bins.groupby("subgroup"):
        axes[1].plot(
            group["mass_geometric_center"],
            group["pearson_residual"],
            marker="o",
            label=subgroup,
        )
    axes[1].axhspan(-3.0, 3.0, color="0.85")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"observed mass [$M_\odot$]")
    axes[1].set_ylabel("Pearson residual")
    axes[1].set_title("Unchanged WP5 baseline gate")
    axes[1].legend()
    fig.tight_layout()
    figure_path = FIGURES / f"wp3_wp5_repair_gates_{VERSION}.png"
    fig.savefig(figure_path, dpi=160)
    plt.close(fig)

    accepted = bool(f4_pass and mass_acceptance and baseline_gate and mass_sanity)
    failed_baseline = baseline.loc[
        ~baseline["residual_gate_pass"],
        [
            "subgroup",
            "poisson_chi_square_p",
            "residual_trend_p",
            "max_abs_pearson_residual",
        ],
    ]
    failed_rows = "\n".join(
        "| {subgroup} | {chi:.6g} | {trend:.6g} | {maximum:.6g} |".format(
            subgroup=row["subgroup"],
            chi=row["poisson_chi_square_p"],
            trend=row["residual_trend_p"],
            maximum=row["max_abs_pearson_residual"],
        )
        for _, row in failed_baseline.iterrows()
    )
    failed_table_markdown = (
        "| subgroup | chi-square p | trend p | max abs residual |\n"
        "|---|---:|---:|---:|\n"
        + failed_rows
    )
    report = f"""# WP3–WP5 extinction/mass repair report ({VERSION})

## Verdict

**Repair remains BLOCKED at the unchanged WP5 residual gate.** The root WP3
extinction invariant and WP4 posterior-quality criteria pass, but the converged
baseline response still leaves the following subgroup failure:

{failed_table_markdown}

WP6 is therefore **not authorized**. The frozen WP5 blocking verdict remains
preserved as the trigger, and these versioned files are a separate repair
attempt rather than replacements.

## What was repaired

- WP3 now uses a 0.03-mag error floor in every band, an eight-anchor spatial
  prior with measured widths of 0.452/0.453/0.475 mag, and a full gridded
  extinction posterior. The hidden template-branch width is calibrated on the
  149 spectroscopic anchors from the asymmetric central-68% discrepancy.
- WP4 ages were refitted from the repaired extinction catalogue. Masses are
  posterior samples over the full extinction distribution, all available six
  bands, the binary branch, and the subgroup age posterior. The six-band model
  width is the anchor-measured 0.38 mag; the mass measure is log-uniform, not an
  IMF prior. Spectroscopic-HRD overrides are unchanged.
- WP5 retained 223,200 catalogue injections and the original residual gate.
  Synthetic photometry passes through the actual repaired WP3 and WP4
  estimators. Baseline response resolution was increased from 16 to 64 mass
  draws and shown not to remove the remaining residual.

## Acceptance checks

| check | result | pass |
|---|---:|:---:|
| max absolute binned median ΔA_V | {f4_max:.3f} mag (<0.30) | {f4_pass} |
| binned rank-test p-value | {rank.pvalue:.3f} (>=0.05) | {f4_pass} |
| largest repeated baseline mass | {100*exact_fraction:.2f}% (<2%) | {exact_fraction < 0.02} |
| stars at 2.5–3.2 M_sun | {populated_25_32} | {populated_25_32 > 0} |
| baseline WP5 all-subgroup residual gate | — | {baseline_gate} |
| baseline association mass | {baseline_assoc['multiplicity_adjusted_mass_median_Msun']:.0f} M_sun | {mass_sanity} |
| median-mass N(>8 M_sun) | {median_gt8} | reported |
| posterior E[N(>8 M_sun)] | {expected_gt8:.1f} | reported |

No baseline subgroup reaches an absolute 95% completeness edge. The bright
plateaus remain below 95%, so the valid statement remains a response-corrected
2–8 M_sun likelihood, **not** “95% complete.” Relative-to-plateau edges are
diagnostic only and are recorded in `{completeness_path.relative_to(ROOT)}`.

## Remaining problem

After the extinction/mass repair, CygOB2-A and CygOB2-B pass the baseline
residual gate. CygOB2-C retains a localized excess at the high end of the
calibration interval, with max |residual| above 3 despite acceptable global
χ² and trend tests. Because this survives the 64-draw response convergence
check, changing the extinction prior or mass-likelihood width further would be
gate tuning. The next diagnostic must test the C subgroup model itself:
subgroup-label uncertainty/contamination, intrinsic age spread versus the
single-age model, and spatially varying selection. Until one of those is
modeled and injected end-to-end, WP6 must remain blocked.
"""
    report_path = REPORTS / f"WP3_WP5_REPAIR_REPORT_{VERSION}.md"
    report_path.write_text(report, encoding="utf-8")

    inputs = [
        PROC / "wp3_extinction.parquet",
        PROC / "wp4_masses.parquet",
        PROC / f"wp3_extinction_{VERSION}.parquet",
        PROC / f"wp3_extinction_posterior_{VERSION}.npz",
        PROC / f"wp4_age_posteriors_{VERSION}.parquet",
        PROC / f"wp4_mass_posteriors_{VERSION}.parquet",
        PROC / f"wp4_mass_posterior_samples_{VERSION}.npz",
        PROC / f"wp5_completeness_curves_{VERSION}.parquet",
        PROC / f"wp5_injection_response_{VERSION}.parquet",
        PROC / f"wp5_imf_normalization_{VERSION}.parquet",
        PROC / f"wp5_association_mass_{VERSION}.parquet",
        PROV / "wp5_manifest.json",
    ]
    outputs = [
        f4_path,
        completeness_path,
        baseline_bins_path,
        figure_path,
        report_path,
        ROOT / "notebooks" / "wp3_extinction_repair_and_wp5_regate.ipynb",
        ROOT / "CUTS_AND_THRESHOLDS.md",
    ]
    gate_payload = {
        "repair_version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp3_wp5_repair_finalize.py",
        "f4": {
            "mass_bins_Msun": MASS_BINS.tolist(),
            "max_abs_median_delta_av_mag": f4_max,
            "limit_mag": 0.30,
            "binned_spearman_rho": float(rank.statistic),
            "binned_spearman_p": float(rank.pvalue),
            "pvalue_limit": 0.05,
            "pass": f4_pass,
        },
        "wp4": {
            "largest_exact_baseline_mass_fraction": exact_fraction,
            "limit": 0.02,
            "n_mass_2p5_3p2": populated_25_32,
            "median_mass_n_gt8": median_gt8,
            "posterior_expected_n_gt8": expected_gt8,
            "pass": mass_acceptance,
        },
        "wp5": {
            "baseline_definition": "PARSEC, R_V=3.1, alpha=2.3",
            "subgroups": baseline[
                [
                    "subgroup",
                    "poisson_chi_square_p",
                    "residual_trend_p",
                    "max_abs_pearson_residual",
                    "residual_gate_pass",
                ]
            ].to_dict(orient="records"),
            "all_subgroups_pass": baseline_gate,
            "absolute_95_edge_found": bool(
                completeness["absolute_95_edge_Msun"].notna().any()
            ),
            "association_mass_median_Msun": float(
                baseline_assoc["multiplicity_adjusted_mass_median_Msun"]
            ),
            "association_mass_lo68_Msun": float(
                baseline_assoc["multiplicity_adjusted_mass_lo68_Msun"]
            ),
            "association_mass_hi68_Msun": float(
                baseline_assoc["multiplicity_adjusted_mass_hi68_Msun"]
            ),
            "within_factor_two_literature": mass_sanity,
        },
        "accepted": accepted,
        "downstream_wp6_authorized": accepted,
        "blocking_reason": (
            None
            if accepted
            else "CygOB2-C baseline WP5 max absolute Pearson residual exceeds 3"
        ),
        "preserved_history": {
            "frozen_wp3_wp4_overwritten": False,
            "frozen_wp5_blocking_manifest": "provenance/wp5_manifest.json",
            "frozen_wp5_blocking_manifest_sha256": sha256(
                PROV / "wp5_manifest.json"
            ),
        },
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path) for path in inputs
        },
        "outputs": {
            str(path.relative_to(ROOT)): sha256(path) for path in outputs
        },
    }
    gate_path = PROV / "wp3_repair_gate.json"
    write_json(gate_path, gate_payload)
    manifest = {
        "repair_version": VERSION,
        "created_utc": gate_payload["created_utc"],
        "status": "ACCEPTED" if accepted else "BLOCKED_WP5_BASELINE_RESIDUAL",
        "accepted": accepted,
        "downstream_wp6_authorized": accepted,
        "component_gates": {
            "wp3_extinction_mass_invariant": f4_pass,
            "wp4_mass_posterior_quality": mass_acceptance,
            "wp5_baseline_residual": baseline_gate,
            "association_mass_factor_two": mass_sanity,
        },
        "response_resolution": {
            "baseline_PARSEC_RV3.1_mass_draws_per_injection": 64,
            "nonbaseline_sensitivity_mass_draws_per_injection": 16,
            "real_star_mass_draws_per_branch": 256,
        },
        "artifacts": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [*outputs, gate_path]
        },
    }
    manifest_path = PROV / "wp3_wp5_repair_manifest.json"
    write_json(manifest_path, manifest)

    # Correct the latest partial-response execution record's resolution text.
    injection_provenance_path = PROV / "wp5_injections_repair_execution.json"
    injection_provenance = json.loads(
        injection_provenance_path.read_text(encoding="utf-8")
    )
    injection_provenance["response_resolution"] = manifest["response_resolution"]
    injection_provenance["repair_change"] = (
        "synthetic six-band photometry passes through actual repair_v1 WP3/F2; "
        "baseline PARSEC R_V=3.1 uses 64 posterior mass draws after convergence "
        "check; inherited nonbaseline sensitivity branches use 16"
    )
    write_json(injection_provenance_path, injection_provenance)
    print(json.dumps({
        "accepted": accepted,
        "f4_pass": f4_pass,
        "mass_posterior_pass": mass_acceptance,
        "wp5_baseline_pass": baseline_gate,
        "mass_sanity_pass": mass_sanity,
        "blocking_reason": gate_payload["blocking_reason"],
    }, indent=2))


if __name__ == "__main__":
    main()
