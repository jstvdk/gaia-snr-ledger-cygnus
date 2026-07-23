#!/usr/bin/env python
"""Freeze canonical anchor assignments and the row-level Berlanas gate audit."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MEMBERS = ROOT / "data" / "processed" / "wp2_members.parquet"
CANONICAL_ANCHORS = ROOT / "data" / "processed" / "wp1_spectroscopic_anchors.parquet"
EVIDENCE = ROOT / "data" / "processed" / "wp1_spectroscopic_anchor_records.parquet"
GAIA = ROOT / "data" / "processed" / "wp1_gaia_narrow.parquet"
ASSIGNMENTS = ROOT / "data" / "processed" / "wp2_anchor_assignments.parquet"
SUBGROUPS = ROOT / "tables" / "wp2_subgroup_labels.parquet"
AUDIT = ROOT / "provenance" / "wp2_berlanas_recovery_audit.csv"
EXECUTION = ROOT / "provenance" / "wp2_final_audits_execution.json"
MANIFEST = ROOT / "provenance" / "wp2_membership_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    members = pd.read_parquet(MEMBERS)
    canonical = pd.read_parquet(CANONICAL_ANCHORS)
    evidence = pd.read_parquet(EVIDENCE)
    gaia = pd.read_parquet(GAIA)
    subgroups = pd.read_parquet(SUBGROUPS)
    for frame in [members, canonical, evidence, gaia]:
        frame["source_id"] = pd.to_numeric(frame["source_id"], errors="coerce").astype("Int64")

    if canonical["source_id"].isna().any() or canonical["source_id"].duplicated().any():
        raise RuntimeError("canonical WP1 anchor source_id is not unique and complete")
    probabilities = members.set_index("source_id")["membership_probability"]
    bases = members.set_index("source_id")["membership_basis"]
    assignments = canonical[[
        "anchor_uid", "source_catalog", "object_name", "source_id", "spectral_type",
        "evidence_record_count", "object_aliases_json", "source_catalogs_json",
    ]].copy()
    assignments["membership_probability"] = assignments["source_id"].map(probabilities)
    assignments["membership_basis"] = assignments["source_id"].map(bases)
    if "subgroup" not in subgroups or "subgroup_label" in subgroups:
        raise RuntimeError("subgroup sidecar must expose only the canonical 'subgroup' column")
    subgroup_by_id = subgroups.set_index("source_id")["subgroup"]
    assignments["subgroup"] = assignments["source_id"].map(subgroup_by_id).fillna("unassigned")
    assignments["assignment_reason"] = np.where(
        assignments["membership_probability"].gt(0.5),
        "astrometric cluster-field mixture or explicit spectroscopic quality exception",
        "posterior P<=0.5 or outside the soft membership artifact",
    )
    assignments.to_parquet(ASSIGNMENTS, index=False)

    berlanas = evidence.loc[
        evidence["source_catalog"].eq("Berlanas et al. 2019"),
        ["catalog_record_id", "object_name", "source_id", "spectral_type", "ra_deg", "dec_deg"],
    ].copy()
    if len(berlanas) != 229 or berlanas["source_id"].isna().any() or berlanas["source_id"].duplicated().any():
        raise RuntimeError("Berlanas evidence denominator is not 229 unique non-null source_ids")
    gaia_fields = gaia[[
        "source_id", "l_deg", "b_deg", "parallax", "parallax_error", "pmra", "pmdec",
        "ruwe", "visibility_periods_used", "phot_bp_rp_excess_factor",
    ]].drop_duplicates("source_id")
    berlanas = berlanas.merge(gaia_fields, on="source_id", how="left", validate="one_to_one")
    member_fields = members[[
        "source_id", "quality_pass", "zero_point_boundary_flag", "membership_probability",
        "membership_probability_astrometric", "membership_probability_mc_se", "membership_basis",
    ]].drop_duplicates("source_id")
    berlanas = berlanas.merge(member_fields, on="source_id", how="left", validate="one_to_one")
    berlanas["recovered_p_gt_0_5"] = berlanas["membership_probability"].gt(0.5)

    reasons = []
    for row in berlanas.itertuples(index=False):
        if row.recovered_p_gt_0_5 and row.membership_basis == "full_covariance_astrometric_mixture":
            reason = "recovered automatically by 10000-draw full-covariance cluster-field posterior"
        elif row.recovered_p_gt_0_5:
            if pd.isna(row.parallax):
                reason = "manual spectroscopic quality exception: Gaia source absent from frozen narrow query"
            elif row.ruwe >= 1.4:
                reason = f"manual spectroscopic quality exception: RUWE={row.ruwe:.3f} >= 1.4"
            elif row.visibility_periods_used < 8:
                reason = f"manual spectroscopic quality exception: visibility_periods_used={row.visibility_periods_used} < 8"
            else:
                reason = "manual spectroscopic quality exception: failed one or more documented astrometric/photometric reliability flags"
        elif pd.isna(row.membership_probability):
            reason = "missed automatically: posterior below the P>0.05 soft-output floor"
        else:
            reason = f"missed automatically: cluster-field Monte Carlo posterior P={row.membership_probability:.4f} <= 0.5"
        reasons.append(reason)
    berlanas["gate_disposition_reason"] = reasons
    berlanas.to_csv(AUDIT, index=False)

    recovered = int(berlanas["recovered_p_gt_0_5"].sum())
    if recovered != 189:
        raise RuntimeError(f"Berlanas audit recovered {recovered}, expected frozen gate value 189")
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp2_finalize_audits.py",
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [MEMBERS, CANONICAL_ANCHORS, EVIDENCE, GAIA]
        },
        "canonical_anchor_assignments": {
            "rows": len(assignments),
            "unique_source_ids": int(assignments["source_id"].nunique()),
            "duplicate_source_ids": int(assignments["source_id"].duplicated().sum()),
        },
        "berlanas_gate_audit": {
            "denominator": len(berlanas),
            "recovered": recovered,
            "missed": int(len(berlanas) - recovered),
            "every_row_has_reason": bool(berlanas["gate_disposition_reason"].notna().all()),
        },
        "outputs": {
            str(ASSIGNMENTS.relative_to(ROOT)): sha256(ASSIGNMENTS),
            str(AUDIT.relative_to(ROOT)): sha256(AUDIT),
        },
    }
    EXECUTION.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest["outputs"][str(ASSIGNMENTS.relative_to(ROOT))] = {
        "sha256": sha256(ASSIGNMENTS), "bytes": ASSIGNMENTS.stat().st_size
    }
    manifest["outputs"][str(AUDIT.relative_to(ROOT))] = {
        "sha256": sha256(AUDIT), "bytes": AUDIT.stat().st_size
    }
    manifest["final_row_level_audits"] = {
        "execution": str(EXECUTION.relative_to(ROOT)),
        "sha256": sha256(EXECUTION),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
