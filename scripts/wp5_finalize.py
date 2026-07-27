#!/usr/bin/env python3
"""Validate and freeze the blocked WP5 diagnostic package."""
from __future__ import annotations

import json
from datetime import datetime, timezone

import nbformat
import numpy as np
import pandas as pd

import wp5_common as w


def file_record(path):
    return {"sha256": w.sha256(path), "bytes": path.stat().st_size}


def main() -> None:
    curves = pd.read_parquet(w.PROC / "wp5_completeness_curves.parquet")
    response = pd.read_parquet(w.PROC / "wp5_injection_response.parquet")
    normalization = pd.read_parquet(w.PROC / "wp5_imf_normalization.parquet")
    bins = pd.read_parquet(w.PROC / "wp5_mass_function_bins.parquet")
    association = pd.read_parquet(w.PROC / "wp5_association_mass.parquet")
    notebook_path = (
        w.NOTEBOOKS / "wp5_imf_normalization_and_completeness.ipynb"
    )
    notebook = nbformat.read(notebook_path, as_version=4)
    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    errors = [
        output
        for cell in code_cells
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]

    expected = {
        "completeness_rows": 2 * 3 * 3 * len(w.MASS_GRID),
        "response_rows": 2
        * 3
        * 3
        * len(w.MASS_GRID)
        * w.N_INJECT_PER_MASS,
        "normalization_rows": 2 * 3 * 3 * 3,
        "mass_bin_rows": 2 * 3 * 3 * 3 * w.N_IMF_BINS,
        "association_rows": 2 * 3 * 3,
    }
    actual = {
        "completeness_rows": len(curves),
        "response_rows": len(response),
        "normalization_rows": len(normalization),
        "mass_bin_rows": len(bins),
        "association_rows": len(association),
    }
    checks = {
        "row_counts_exact": actual == expected,
        "branch_grid_complete": (
            set(normalization["family"]) == set(w.FAMILIES)
            and set(normalization["R_V"]) == set(w.R_V_BRANCHES)
            and set(normalization["alpha"]) == set(w.IMF_SLOPES)
            and set(normalization["subgroup"]) == set(w.SUBGROUPS)
        ),
        "minimum_50_calibrators": bool(
            normalization["raw_calibration_sources"].ge(50).all()
        ),
        "no_absolute_95_edge_preserved": bool(
            normalization["absolute_95_edge_Msun"].isna().all()
        ),
        "all_residual_failures_preserved": bool(
            (~normalization["residual_gate_pass"]).all()
        ),
        "baseline_mass_factor_two": bool(
            association.loc[
                association["family"].eq("PARSEC")
                & association["R_V"].eq(3.1)
                & association["alpha"].eq(2.3),
                "within_factor_two_literature",
            ].iloc[0]
        ),
        "notebook_all_code_cells_executed": all(
            cell.get("execution_count") is not None for cell in code_cells
        ),
        "notebook_has_no_error_outputs": len(errors) == 0,
        "upstream_wp4_mass_hash_unchanged": (
            w.sha256(w.PROC / "wp4_masses.parquet")
            == json.loads(
                (w.PROVENANCE / "wp4_manifest.json").read_text(
                    encoding="utf-8"
                )
            )["outputs"]["data/processed/wp4_masses.parquet"]["sha256"]
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"WP5 validation failed unexpectedly: {checks}")

    notebook_log = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "notebook": str(notebook_path.relative_to(w.ROOT)),
        "sha256": w.sha256(notebook_path),
        "code_cells": len(code_cells),
        "executed_code_cells": sum(
            cell.get("execution_count") is not None for cell in code_cells
        ),
        "error_outputs": len(errors),
        "status": "EXECUTED_SUCCESSFULLY_REPORTS_BLOCKED_SCIENCE_GATE",
    }
    notebook_log_path = w.PROVENANCE / "wp5_notebook_execution.json"
    w.write_json(notebook_log_path, notebook_log)

    validation = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "VALIDATION_EXECUTED_EXPECTED_SCIENCE_GATE_FAILURE_PRESERVED",
        "expected_rows": expected,
        "actual_rows": actual,
        "checks": checks,
        "science_gate": {
            "accepted": False,
            "blocking_reasons": [
                "no absolute 95% completeness edge in any branch/subgroup",
                "0/54 response-aware Poisson IMF fits pass residual gate",
            ],
            "downstream_authority": False,
        },
    }
    validation_path = w.PROVENANCE / "wp5_validation.json"
    w.write_json(validation_path, validation)

    upstream = [
        w.PROC / "wp1_gaia_narrow.parquet",
        w.PROC / "wp1_2mass_join.parquet",
        w.PROC / "wp2_members.parquet",
        w.TABLES / "wp2_subgroup_labels.parquet",
        w.PROC / "wp3_extinction.parquet",
        w.PROC / "wp3_isochrones_parsec.parquet",
        w.PROC / "wp3_isochrones_mist.parquet",
        w.PROC / "wp4_age_posteriors.parquet",
        w.PROC / "wp4_masses.parquet",
    ]
    scripts = [
        w.ROOT / "scripts" / "wp5_common.py",
        w.ROOT / "scripts" / "wp5_injections.py",
        w.ROOT / "scripts" / "wp5_fit_imf.py",
        w.ROOT / "scripts" / "wp5_report.py",
        w.ROOT / "scripts" / "wp5_make_notebook.py",
        w.ROOT / "scripts" / "wp5_finalize.py",
    ]
    logs = [
        w.PROVENANCE / "wp5_injections_execution.json",
        w.PROVENANCE / "wp5_imf_fit_execution.json",
        w.PROVENANCE / "wp5_report_execution.json",
        notebook_log_path,
        validation_path,
    ]
    outputs = [
        w.PROC / "wp5_completeness_curves.parquet",
        w.PROC / "wp5_injection_response.parquet",
        w.PROC / "wp5_imf_normalization.parquet",
        w.PROC / "wp5_mass_function_bins.parquet",
        w.PROC / "wp5_association_mass.parquet",
        w.PROC / "wp5_imf_posterior_draws.npz",
        w.TABLES / "wp5_imf_norm.csv",
        w.TABLES / "wp5_imf_norm.md",
        w.TABLES / "wp5_association_mass.csv",
        w.FIGURES / "wp5_completeness_curves.png",
        w.FIGURES / "wp5_mass_function.png",
        w.FIGURES / "wp5_association_mass.png",
        notebook_path,
        w.ROOT / "wp5_imf_norm.md",
        w.ROOT / "wp5_completion_report.md",
        w.PROVENANCE / "wp5_provenance.md",
    ]
    manifest = {
        "work_package": "WP5 - IMF normalization and completeness",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "WP5_BLOCKED_RESIDUAL_AND_ABSOLUTE_95_COMPLETENESS_GATES",
        "accepted": False,
        "downstream_authority": False,
        "seed": w.SEED,
        "branches": {
            "subgroups": w.SUBGROUPS,
            "isochrone_family": w.FAMILIES,
            "R_V": w.R_V_BRANCHES,
            "imf_slope_alpha": w.IMF_SLOPES,
        },
        "upstream_inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in upstream
        },
        "scripts": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in scripts
        },
        "execution_logs": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in logs
        },
        "outputs": {
            str(path.relative_to(w.ROOT)): file_record(path)
            for path in outputs
        },
        "gate": validation["science_gate"],
        "checks": checks,
        "required_remediation": (
            "revise WP4 2--5 Msun mass inference to posterior/direct-CMD form, "
            "then rerun the preserved WP5 response and untuned residual gate"
        ),
    }
    manifest_path = w.PROVENANCE / "wp5_manifest.json"
    w.write_json(manifest_path, manifest)
    print(json.dumps({"manifest": str(manifest_path), "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
