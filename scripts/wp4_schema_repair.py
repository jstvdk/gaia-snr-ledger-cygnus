#!/usr/bin/env python3
"""Repair the canonical subgroup and WP4 mass-column schemas in place.

This is a migration, not a scientific re-analysis.  The authoritative subgroup
mapping is ``tables/wp2_subgroup_labels.parquet``.  Its one label column is
``subgroup``; every propagated WP2/WP3/WP4 copy is rebuilt from that sidecar.
Rows outside the 1,331 clean automatic-member sidecar are explicitly labelled
``unassigned``.  In particular, the 61 P>0.5 anchor-quality exceptions remain
unassigned by construction.

The migration also renames the reporting-only WP4 column ``mass_best`` to
``mass_baseline``.  No values are recalculated.

Run from the project root with the project environment:

    conda run -n cygob2-gaia --no-capture-output python scripts/wp4_schema_repair.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
TABLES = ROOT / "tables"
PROV = ROOT / "provenance"
SIDECAR = TABLES / "wp2_subgroup_labels.parquet"
SUBGROUP_VALUES = {"CygOB2-A", "CygOB2-B", "CygOB2-C", "unassigned"}
MASS_BRANCHES = [
    f"mass_{family}_rv{rv}"
    for family in ("PARSEC", "MIST")
    for rv in ("3.0", "3.1", "3.5")
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_parquet_atomic(frame: pd.DataFrame, path: Path) -> None:
    """Replace one Parquet product only after its complete successor exists."""
    temporary = path.with_suffix(path.suffix + ".schema-repair.tmp")
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, path)


def canonical_sidecar() -> pd.DataFrame:
    sidecar = pd.read_parquet(SIDECAR)
    if "subgroup_label" in sidecar.columns:
        if "subgroup" in sidecar.columns:
            raise RuntimeError("sidecar contains both subgroup schema names")
        sidecar = sidecar.rename(columns={"subgroup_label": "subgroup"})
    required = {"source_id", "subgroup"}
    if not required.issubset(sidecar.columns):
        raise RuntimeError(f"sidecar lacks required columns: {required - set(sidecar.columns)}")
    sidecar = sidecar[["source_id", "subgroup"]].copy()
    sidecar["source_id"] = pd.to_numeric(
        sidecar["source_id"], errors="raise"
    ).astype("int64")
    if sidecar["source_id"].duplicated().any():
        raise RuntimeError("authoritative subgroup sidecar has duplicate source_id values")
    if len(sidecar) != 1331:
        raise RuntimeError(f"expected 1,331 authoritative labels, found {len(sidecar)}")
    if not set(sidecar["subgroup"]).issubset(SUBGROUP_VALUES - {"unassigned"}):
        raise RuntimeError("authoritative sidecar contains an unexpected subgroup value")
    write_parquet_atomic(sidecar, SIDECAR)
    return sidecar


def propagate_subgroup(path: Path, sidecar: pd.DataFrame) -> dict:
    frame = pd.read_parquet(path)
    before_rows = len(frame)
    stale = [column for column in frame.columns if column in {"subgroup", "subgroup_label"}]
    frame = frame.drop(columns=stale, errors="ignore")
    frame["source_id"] = pd.to_numeric(frame["source_id"], errors="raise").astype("int64")
    frame = frame.merge(sidecar, on="source_id", how="left", validate="many_to_one")
    frame["subgroup"] = frame["subgroup"].fillna("unassigned")
    if len(frame) != before_rows:
        raise RuntimeError(f"row count changed while repairing {path}")
    if not set(frame["subgroup"]).issubset(SUBGROUP_VALUES):
        raise RuntimeError(f"unexpected subgroup value in {path}")
    write_parquet_atomic(frame, path)
    return {
        "rows": int(len(frame)),
        "label_counts": {
            str(key): int(value)
            for key, value in frame["subgroup"].value_counts(dropna=False).items()
        },
        "sha256": sha256(path),
    }


def repair_masses(sidecar: pd.DataFrame) -> dict:
    path = PROC / "wp4_masses.parquet"
    frame = pd.read_parquet(path)
    if "mass_best" in frame.columns:
        if "mass_baseline" in frame.columns:
            raise RuntimeError("mass product contains both baseline column names")
        frame = frame.rename(columns={"mass_best": "mass_baseline"})
    missing = [column for column in MASS_BRANCHES if column not in frame.columns]
    if missing:
        raise RuntimeError(f"WP4 mass branches missing: {missing}")
    null_sets = [
        set(frame.loc[frame[column].isna(), "source_id"].astype("int64"))
        for column in MASS_BRANCHES
    ]
    if not all(values == null_sets[0] for values in null_sets[1:]):
        raise RuntimeError("the six WP4 mass branches do not share one null source_id set")

    stale = [column for column in frame.columns if column in {"subgroup", "subgroup_label"}]
    frame = frame.drop(columns=stale, errors="ignore")
    frame["source_id"] = pd.to_numeric(frame["source_id"], errors="raise").astype("int64")
    frame = frame.merge(sidecar, on="source_id", how="left", validate="one_to_one")
    frame["subgroup"] = frame["subgroup"].fillna("unassigned")
    write_parquet_atomic(frame, path)
    return {
        "rows": int(len(frame)),
        "mass_branch_columns": MASS_BRANCHES,
        "shared_null_source_ids": int(len(null_sets[0])),
        "mass_baseline_definition": "mass_PARSEC_rv3.1",
        "label_counts": {
            str(key): int(value)
            for key, value in frame["subgroup"].value_counts(dropna=False).items()
        },
        "sha256": sha256(path),
    }


def main() -> None:
    before = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in [
            SIDECAR,
            PROC / "wp2_members.parquet",
            PROC / "wp2_anchor_assignments.parquet",
            PROC / "wp3_member_photometry.parquet",
            PROC / "wp3_extinction.parquet",
            PROC / "wp4_masses.parquet",
            PROC / "wp4_anchor_hrd.parquet",
        ]
    }
    sidecar = canonical_sidecar()
    outputs = {}
    for relative in [
        "data/processed/wp2_members.parquet",
        "data/processed/wp2_anchor_assignments.parquet",
        "data/processed/wp3_member_photometry.parquet",
        "data/processed/wp3_extinction.parquet",
        "data/processed/wp4_anchor_hrd.parquet",
    ]:
        path = ROOT / relative
        outputs[relative] = propagate_subgroup(path, sidecar)
    outputs["data/processed/wp4_masses.parquet"] = repair_masses(sidecar)

    members = pd.read_parquet(PROC / "wp2_members.parquet")
    gate = members["membership_probability"].gt(0.5)
    exempt = members["anchor_quality_exempt"].fillna(False).astype(bool)
    if int((gate & exempt & members["subgroup"].eq("unassigned")).sum()) != 61:
        raise RuntimeError("the 61 anchor-quality exemptions are not all explicitly unassigned")
    if int((gate & ~exempt & members["subgroup"].eq("unassigned")).sum()) != 0:
        raise RuntimeError("a clean automatic P>0.5 member lacks an authoritative subgroup")

    log = {
        "script": "scripts/wp4_schema_repair.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "operation": "schema-only migration; no membership, extinction, age, or mass re-analysis",
        "authoritative_subgroup_source": "tables/wp2_subgroup_labels.parquet",
        "canonical_subgroup_column": "subgroup",
        "sidecar_rows": int(len(sidecar)),
        "sidecar_label_counts": {
            str(key): int(value)
            for key, value in sidecar["subgroup"].value_counts().items()
        },
        "before_sha256": before,
        "outputs": outputs,
        "gate_member_anchor_exempt_unassigned": 61,
        "gate_member_clean_automatic_unassigned": 0,
    }
    log_path = PROV / "wp4_schema_repair_execution.json"
    log_path.write_text(json.dumps(log, indent=2) + "\n")
    print("schema repair complete")
    print("  authoritative labels:", len(sidecar))
    print("  P>0.5 anchor exemptions explicitly unassigned: 61")
    print("  shared massless source_ids across all six branches: 55")


if __name__ == "__main__":
    main()
