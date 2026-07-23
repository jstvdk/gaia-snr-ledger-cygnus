#!/usr/bin/env python3
"""Freeze the read-only WP4 closure checks and WP6 proper-motion hand-off.

The only tabular output is a candidate list: P>0.5 stars more than five
standard deviations from the clean automatic-member mean in either proper
motion component.  This is a WP6 follow-up flag, not a membership cut.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
PROV = ROOT / "provenance"
PM_Z_THRESHOLD = 5.0
MASS_BRANCHES = [
    f"mass_{family}_rv{rv}"
    for family in ("PARSEC", "MIST")
    for rv in ("3.0", "3.1", "3.5")
]
SUBGROUP_PRODUCTS = [
    ROOT / "tables" / "wp2_subgroup_labels.parquet",
    PROC / "wp2_members.parquet",
    PROC / "wp2_anchor_assignments.parquet",
    PROC / "wp3_member_photometry.parquet",
    PROC / "wp3_extinction.parquet",
    PROC / "wp4_masses.parquet",
    PROC / "wp4_anchor_hrd.parquet",
    PROC / "wp4_age_posteriors.parquet",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    members = pd.read_parquet(PROC / "wp2_members.parquet")
    gate = members[members["membership_probability"] > 0.5].copy()
    clean = gate[
        ~gate["anchor_quality_exempt"].fillna(False).astype(bool)
    ].copy()

    pm_reference = {}
    outlier = np.zeros(len(gate), dtype=bool)
    for component in ("pmra", "pmdec"):
        mean = float(clean[component].mean())
        std = float(clean[component].std(ddof=1))
        gate[f"z_{component}"] = (gate[component] - mean) / std
        outlier |= gate[f"z_{component}"].abs().gt(PM_Z_THRESHOLD).fillna(False)
        pm_reference[component] = {
            "mean_mas_per_yr": mean,
            "std_mas_per_yr": std,
        }

    candidates = gate.loc[outlier, [
        "source_id", "pmra", "pmdec", "z_pmra", "z_pmdec",
        "membership_probability", "subgroup", "anchor_quality_exempt", "ruwe",
        "membership_basis",
    ]].copy()
    candidates["anchor_status"] = np.where(
        candidates["anchor_quality_exempt"],
        "spectroscopic_quality_exception",
        "automatic_member",
    )
    candidates["outlier_rule"] = (
        f"abs(z_pmra)>{PM_Z_THRESHOLD:g} or abs(z_pmdec)>{PM_Z_THRESHOLD:g}; "
        "z relative to 1,331 clean automatic P>0.5 members"
    )
    candidates = candidates.sort_values(
        ["z_pmdec", "z_pmra"], key=lambda values: values.abs(), ascending=False
    )
    candidate_path = PROV / "wp4_pm_outliers.csv"
    candidates.to_csv(candidate_path, index=False, float_format="%.6f")

    missing_astrometry = gate[gate["ra"].isna()].copy()
    expected_missing_basis = (
        "Berlanas2019_spectroscopic_member_manual_quality_exception"
    )
    if len(missing_astrometry) != 9:
        raise RuntimeError(f"expected nine astrometry-less gate members, found {len(missing_astrometry)}")
    if not missing_astrometry["membership_basis"].eq(expected_missing_basis).all():
        raise RuntimeError("an astrometry-less member lacks the documented Berlanas basis")
    if not missing_astrometry["membership_probability"].eq(1.0).all():
        raise RuntimeError("an astrometry-less manual exception does not have P=1")

    high_ruwe = gate[gate["ruwe"] > 1.4]
    if not high_ruwe["anchor_quality_exempt"].all():
        raise RuntimeError("a P>0.5 high-RUWE member is not an explicit anchor exemption")

    schema = {}
    for path in SUBGROUP_PRODUCTS:
        frame = pd.read_parquet(path)
        columns = [column for column in frame.columns if column.startswith("subgroup")]
        if columns != ["subgroup"]:
            raise RuntimeError(f"{path} has non-canonical subgroup columns: {columns}")
        placeholder = frame["subgroup"].astype(str).str.contains(
            "placeholder|distance_structure_unresolved|TODO|TBD|FIXME|dummy",
            case=False, regex=True, na=False,
        )
        if placeholder.any():
            raise RuntimeError(f"{path} contains placeholder subgroup values")
        schema[str(path.relative_to(ROOT))] = {
            "rows": int(len(frame)),
            "subgroup_columns": columns,
            "counts": {
                str(key): int(value)
                for key, value in frame["subgroup"].value_counts(dropna=False).items()
            },
        }

    masses = pd.read_parquet(PROC / "wp4_masses.parquet")
    if "mass_best" in masses or "mass_baseline" not in masses:
        raise RuntimeError("WP4 baseline mass schema is not repaired")
    null_sets = [
        set(masses.loc[masses[column].isna(), "source_id"].astype("int64"))
        for column in MASS_BRANCHES
    ]
    if not all(values == null_sets[0] for values in null_sets[1:]):
        raise RuntimeError("mass-branch null source_id sets differ")
    if len(null_sets[0]) != 55:
        raise RuntimeError(f"expected 55 shared massless stars, found {len(null_sets[0])}")

    ages = pd.read_parquet(PROC / "wp4_age_posteriors.parquet")
    required_age_columns = {"measurable", "grid_railed", "exclusion_reason"}
    if not required_age_columns.issubset(ages.columns):
        raise RuntimeError("age product lacks closure measurability columns")
    if ages.loc[ages["measurable"], "grid_railed"].any():
        raise RuntimeError("a measurable age branch is grid-railed")
    exclusions = {
        str(key): int(value)
        for key, value in ages.loc[
            ~ages["measurable"], "exclusion_reason"
        ].value_counts().items()
    }

    anchors = pd.read_parquet(PROC / "wp4_anchor_hrd.parquet")
    anchor_counts = {
        str(key): int(value)
        for key, value in anchors["subgroup"].value_counts().items()
    }
    if anchor_counts.get("CygOB2-B") != 5:
        raise RuntimeError("expected five subgroup-B HRD anchors")

    berlanas = pd.read_csv(PROV / "wp2_berlanas_recovery_audit.csv")
    recovered = int(berlanas["recovered_p_gt_0_5"].sum())
    if len(berlanas) != 229 or recovered != 189:
        raise RuntimeError("WP2 Berlanas gate changed during WP4 closure")

    log = {
        "script": "scripts/wp4_closure_audit.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "operation": "read-only closure audit plus WP6 candidate export; no member removed",
        "wp2_gate": {
            "soft_rows": int(len(members)),
            "p_gt_0_5_rows": int(len(gate)),
            "clean_automatic_rows": int(len(clean)),
            "anchor_quality_exempt_rows": int(gate["anchor_quality_exempt"].sum()),
            "berlanas_recovered": recovered,
            "berlanas_denominator": int(len(berlanas)),
        },
        "subgroup_schema": schema,
        "astrometry_less_gate_members": {
            "count": int(len(missing_astrometry)),
            "source_ids": [
                int(value) for value in missing_astrometry["source_id"].tolist()
            ],
            "mechanism": (
                "explicit Berlanas et al. 2019 spectroscopic membership manual "
                "quality exceptions; absent from the frozen narrow Gaia query, "
                "therefore P=1 by literature basis and no astrometric posterior"
            ),
        },
        "high_ruwe_gate_members": {
            "count": int(len(high_ruwe)),
            "all_anchor_quality_exempt": True,
            "maximum_ruwe": float(high_ruwe["ruwe"].max()),
        },
        "proper_motion_candidates": {
            "threshold_sigma": PM_Z_THRESHOLD,
            "reference": pm_reference,
            "count": int(len(candidates)),
            "source_ids": [int(value) for value in candidates["source_id"]],
            "output": str(candidate_path.relative_to(ROOT)),
            "sha256": sha256(candidate_path),
            "membership_action": "none; hand off to WP6 runaway analysis",
        },
        "mass_schema": {
            "branches": MASS_BRANCHES,
            "all_six_null_sets_identical": True,
            "shared_massless_source_ids": 55,
            "mass_baseline_definition": "mass_PARSEC_rv3.1; reporting only",
        },
        "age_measurability": {
            "rows": int(len(ages)),
            "measurable": int(ages["measurable"].sum()),
            "excluded": int((~ages["measurable"]).sum()),
            "grid_railed": int(ages["grid_railed"].sum()),
            "measurable_and_grid_railed": int(
                (ages["measurable"] & ages["grid_railed"]).sum()
            ),
            "exclusion_reason_counts": exclusions,
            "retained_map_envelope_myr": [
                float(ages.loc[ages["measurable"], "age_map"].min()),
                float(ages.loc[ages["measurable"], "age_map"].max()),
            ],
        },
        "anchor_counts": anchor_counts,
    }
    log_path = PROV / "wp4_closure_audit.json"
    log_path.write_text(json.dumps(log, indent=2) + "\n")
    print(f"wrote {candidate_path} ({len(candidates)} WP6 candidates)")
    print(f"wrote {log_path}")
    print("WP2 counts and Berlanas recall unchanged")


if __name__ == "__main__":
    main()
