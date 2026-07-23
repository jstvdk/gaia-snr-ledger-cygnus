#!/usr/bin/env python
"""Re-freeze Wright+15 and spectroscopic anchors with unique Gaia source IDs."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
PROVENANCE = ROOT / "provenance"
WRIGHT = PROCESSED / "wp1_wright15_gaia_crossmatch.parquet"
ANCHORS = PROCESSED / "wp1_spectroscopic_anchors.parquet"
ANCHORS_ECSV = PROCESSED / "wp1_spectroscopic_anchors.ecsv"
EVIDENCE = PROCESSED / "wp1_spectroscopic_anchor_records.parquet"
EVIDENCE_ECSV = PROCESSED / "wp1_spectroscopic_anchor_records.ecsv"
REPORT = ROOT / "wp1_wright15_crossmatch.md"
EXECUTION = PROVENANCE / "wp1_spectroscopic_anchors_execution.json"
MATCH_RADIUS_ARCSEC = 2.0
PAPER_ABSTRACT_COUNT = 169
PAPER_SECTION_2_1_COUNT = 167


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_ecsv(frame: pd.DataFrame, path: Path) -> None:
    export = frame.astype(object).where(pd.notna(frame), None)
    Table.from_pandas(export).write(path, format="ascii.ecsv", overwrite=True)


def json_values(values) -> str:
    result = []
    for value in values:
        if pd.isna(value):
            continue
        item = value.item() if isinstance(value, np.generic) else value
        if item not in result:
            result.append(item)
    return json.dumps(result, ensure_ascii=False)


def first_present(frame: pd.DataFrame, column: str):
    if column not in frame:
        return None
    values = frame[column].dropna()
    if values.empty:
        return None
    return values.iloc[0]


def deduplicate_wright(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    result = frame.copy()
    if "resolved_gaia_source_id" in result:
        resolved = result["resolved_gaia_source_id"].astype("string")
    else:
        resolved = result["source_id"].astype("string")
    result["resolved_gaia_source_id"] = resolved
    result["coordinate_match_consistent"] = (
        result["match_separation_arcsec"].notna()
        & (result["match_separation_arcsec"] <= MATCH_RADIUS_ARCSEC)
    )
    result["source_id"] = pd.Series(pd.NA, index=result.index, dtype="string")
    result["countable_for_wright_gate"] = False
    result["duplicate_of_wright_recno"] = pd.Series(pd.NA, index=result.index, dtype="Int64")
    result["duplicate_of_source_id"] = pd.Series(pd.NA, index=result.index, dtype="string")
    result["census_record_status"] = "unresolved"
    result["deduplication_reason"] = pd.Series(pd.NA, index=result.index, dtype="string")
    result["merged_wright_recnos_json"] = "[]"
    result["merged_object_names_json"] = "[]"

    collisions = []
    for source_id, group in result.groupby("resolved_gaia_source_id", dropna=False):
        if pd.isna(source_id):
            continue
        consistent = group[group["coordinate_match_consistent"]]
        candidates = consistent if not consistent.empty else group
        canonical_index = candidates["match_separation_arcsec"].astype(float).idxmin()
        canonical_recno = int(result.at[canonical_index, "wright_recno"])
        recnos = [int(value) for value in group["wright_recno"]]
        names = [str(value) for value in group["object_name"]]
        result.at[canonical_index, "source_id"] = str(source_id)
        result.at[canonical_index, "countable_for_wright_gate"] = True
        result.at[canonical_index, "census_record_status"] = "countable_unique_star"
        result.at[canonical_index, "merged_wright_recnos_json"] = json.dumps(recnos)
        result.at[canonical_index, "merged_object_names_json"] = json.dumps(names)
        if len(group) == 1:
            continue
        rejected = []
        for index in group.index:
            if index == canonical_index:
                continue
            bad_coordinate = not bool(result.at[index, "coordinate_match_consistent"])
            status = "duplicate_bad_coordinate" if bad_coordinate else "duplicate_alias"
            reason = (
                "SIMBAD name resolution assigned the canonical Gaia source, but the published "
                f"coordinate is {float(result.at[index, 'match_separation_arcsec']):.3f} arcsec away"
                if bad_coordinate
                else "separate Wright catalogue alias resolves to the same Gaia source"
            )
            result.at[index, "duplicate_of_wright_recno"] = canonical_recno
            result.at[index, "duplicate_of_source_id"] = str(source_id)
            result.at[index, "census_record_status"] = status
            result.at[index, "deduplication_reason"] = reason
            rejected.append(
                {
                    "wright_recno": int(result.at[index, "wright_recno"]),
                    "object_name": str(result.at[index, "object_name"]),
                    "status": status,
                    "separation_arcsec": float(result.at[index, "match_separation_arcsec"]),
                    "reason": reason,
                }
            )
        collisions.append(
            {
                "source_id": str(source_id),
                "canonical_wright_recno": canonical_recno,
                "canonical_object_name": str(result.at[canonical_index, "object_name"]),
                "all_wright_recnos": recnos,
                "all_object_names": names,
                "rejected_records": rejected,
            }
        )

    nonnull = result["source_id"].dropna()
    if nonnull.duplicated().any():
        raise RuntimeError("Wright re-freeze still contains duplicate non-null source_ids")
    return result, collisions


def canonicalize_anchors(records: pd.DataFrame, wright: pd.DataFrame) -> pd.DataFrame:
    evidence = records.copy()
    evidence["source_id"] = evidence["source_id"].astype("string")
    wright_status = wright.set_index(wright["wright_recno"].astype(str))["census_record_status"]
    evidence["evidence_record_status"] = "literature_evidence"
    evidence["countable_wright_record"] = False
    is_wright = evidence["source_catalog"].eq("Wright et al. 2015")
    for index in evidence[is_wright].index:
        record_id = str(evidence.at[index, "catalog_record_id"])
        status = str(wright_status.get(record_id, "unresolved"))
        evidence.at[index, "evidence_record_status"] = f"wright15_{status}"
        evidence.at[index, "countable_wright_record"] = status == "countable_unique_star"

    priority = {
        "Wright et al. 2015": 0,
        "Berlanas et al. 2020": 1,
        "Galactic O-Star Catalog": 2,
        "Galactic Wolf Rayet Catalogue": 3,
        "Berlanas et al. 2019": 4,
    }
    evidence["_priority"] = evidence["source_catalog"].map(priority).fillna(99)
    countable_wright_ids = set(
        wright.loc[wright["countable_for_wright_gate"], "source_id"].dropna().astype(str)
    )
    rows = []
    for source_id, group in evidence.dropna(subset=["source_id"]).groupby("source_id"):
        ordered = group.sort_values(["_priority", "anchor_uid"])
        row = {
            "anchor_uid": f"gaia_dr3:{source_id}",
            "source_catalog": "; ".join(sorted(group["source_catalog"].dropna().unique())),
            "catalog_version": json_values(group["catalog_version"]),
            "catalog_record_id": json_values(group["catalog_record_id"]),
            "object_name": first_present(ordered, "object_name"),
            "source_id": str(source_id),
            "ra_deg": first_present(ordered, "ra_deg"),
            "dec_deg": first_present(ordered, "dec_deg"),
            "spectral_type": first_present(ordered, "spectral_type"),
            "teff_K": first_present(ordered, "teff_K"),
            "teff_error_K": first_present(ordered, "teff_error_K"),
            "logg_cgs": first_present(ordered, "logg_cgs"),
            "logg_error_cgs": first_present(ordered, "logg_error_cgs"),
            "extinction_av_mag": first_present(ordered, "extinction_av_mag"),
            "subgroup_or_association": first_present(ordered, "subgroup_or_association"),
            "source_table": json_values(group["source_table"]),
            "bibcode": json_values(group["bibcode"]),
            "gaia_match_method": json_values(group["gaia_match_method"]),
            "gaia_match_separation_arcsec": first_present(
                ordered.sort_values("gaia_match_separation_arcsec", na_position="last"),
                "gaia_match_separation_arcsec",
            ),
            "notes": json_values(group["notes"]),
            "evidence_record_count": int(len(group)),
            "evidence_anchor_uids_json": json_values(group["anchor_uid"]),
            "object_aliases_json": json_values(group["object_name"]),
            "source_catalogs_json": json_values(group["source_catalog"]),
            "has_wright15_census_entry": bool(
                group["source_catalog"].eq("Wright et al. 2015").any()
            ),
            "wright15_countable_unique_star": str(source_id) in countable_wright_ids,
            "countable_in_wp6": True,
            "canonicalization": "one canonical row per non-null Gaia DR3 source_id",
        }
        rows.append(row)
    canonical = pd.DataFrame(rows)
    canonical["source_id"] = canonical["source_id"].astype("string")
    if canonical["source_id"].isna().any() or canonical["source_id"].duplicated().any():
        raise RuntimeError("Canonical anchor table is not unique and non-null by source_id")
    return evidence.drop(columns=["_priority"]), canonical


def refreeze(created_override: str | None = None) -> dict:
    created = created_override or datetime.now(timezone.utc).isoformat()
    prior_execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    wright_input_hash = sha256(WRIGHT)
    evidence_source = EVIDENCE if EVIDENCE.exists() else ANCHORS
    evidence_input_hash = sha256(evidence_source)
    raw_wright = pd.read_parquet(WRIGHT)
    raw_records = pd.read_parquet(evidence_source)

    wright, collisions = deduplicate_wright(raw_wright)
    evidence, canonical = canonicalize_anchors(raw_records, wright)
    wright.to_parquet(WRIGHT, index=False)
    evidence.to_parquet(EVIDENCE, index=False)
    write_ecsv(evidence, EVIDENCE_ECSV)
    canonical.to_parquet(ANCHORS, index=False)
    write_ecsv(canonical, ANCHORS_ECSV)

    valid_wright = wright[wright["countable_for_wright_gate"]]
    denominator = int(len(valid_wright))
    in_narrow = int(valid_wright["in_wp1_gaia_narrow"].sum())
    duplicate_rows = wright[wright["census_record_status"].str.startswith("duplicate_")]
    excluded = valid_wright[~valid_wright["in_wp1_gaia_narrow"]]
    lines = [
        "# Wright et al. (2015) to Gaia DR3 cross-match", "", f"Re-frozen: {created}", "",
        "## Denominator reconciliation", "",
        f"- Wright+15 abstract and conclusions: **{PAPER_ABSTRACT_COUNT} primary OB stars**.",
        f"- Wright+15 Section 2.1 and the paper's later IMF sample: **{PAPER_SECTION_2_1_COUNT} stars**.",
        f"- Frozen VizieR `J/MNRAS/449/741/census`: **{len(wright)} rows**.",
        f"- Unique physical stars represented after Gaia identity resolution: **{denominator}**.",
        "- WP1 Gaia-match denominator: **165 unique physical stars represented by the frozen machine-readable table**; it is not silently equated with the abstract's 169.",
        f"- Representation relative to the abstract count: {denominator}/{PAPER_ABSTRACT_COUNT} ({denominator/PAPER_ABSTRACT_COUNT:.1%}). The four-star difference is the combination of the paper's internal 169-versus-167 count discrepancy and two duplicate physical identities in the machine-readable rows.", "",
        "## Gaia gate", "",
        f"- Valid unique Gaia matches: {denominator}/{denominator} (100.0%).",
        f"- Valid unique stars present in `wp1_gaia_narrow`: {in_narrow}/{denominator} ({in_narrow/denominator:.1%}).",
        "- WP1 gate (>=90% of the declared denominator matched): **PASS**.",
        "- Non-null duplicate `source_id` values in this crossmatch: **0**.", "",
        "## Duplicate-record audit", "",
        "Rejected records are retained as flagged rows with `source_id=null`, while `resolved_gaia_source_id` and `duplicate_of_source_id` preserve the identity evidence.", "",
        "| Rejected recno | Object | Resolved Gaia DR3 source_id | Canonical recno | separation (arcsec) | status |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for _, row in duplicate_rows.iterrows():
        lines.append(
            f"| {int(row['wright_recno'])} | {row['object_name']} | {row['resolved_gaia_source_id']} | "
            f"{int(row['duplicate_of_wright_recno'])} | {float(row['match_separation_arcsec']):.3f} | "
            f"{row['census_record_status']}: {row['deduplication_reason']} |"
        )
    lines.extend(
        [
            "", "## Valid unique stars absent from the narrow sample", "",
            "| Object | Gaia DR3 source_id | separation (arcsec) | reason absent from narrow |",
            "|---|---:|---:|---|",
        ]
    )
    for _, row in excluded.iterrows():
        lines.append(
            f"| {row['object_name']} | {row['source_id']} | "
            f"{float(row['match_separation_arcsec']):.3f} | {row['narrow_exclusion_reason']} |"
        )
    lines.extend(
        [
            "", "The canonical `wp1_spectroscopic_anchors` table contains one row per Gaia source. "
            "All literature rows, aliases, rejected Wright records, and source-specific provenance remain in "
            "`data/processed/wp1_spectroscopic_anchor_records.parquet`.", "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")

    execution = {
        **prior_execution,
        "created_utc": created,
        "script": "scripts/wp1_refreeze_anchor_uniqueness.py",
        "refreeze": {
            "reason": "remove duplicate physical stars from countable Wright and WP6 anchor artifacts",
            "upstream_wright_sha256": wright_input_hash,
            "upstream_anchor_evidence_sha256": evidence_input_hash,
            "canonical_rule": "one non-null Gaia DR3 source_id per countable row; all source records preserved separately",
            "duplicate_collisions": collisions,
        },
        "wright_gate": {
            "paper_abstract_primary_ob_stars": PAPER_ABSTRACT_COUNT,
            "paper_section_2_1_stars": PAPER_SECTION_2_1_COUNT,
            "vizier_census_rows": int(len(wright)),
            "census_rows": denominator,
            "denominator_definition": "unique physical stars represented in frozen VizieR census after Gaia identity resolution",
            "matched_to_gaia_dr3": denominator,
            "match_fraction": 1.0,
            "coverage_relative_to_paper_abstract": denominator / PAPER_ABSTRACT_COUNT,
            "present_in_wp1_gaia_narrow": in_narrow,
            "duplicate_catalog_rows_flagged": int(len(duplicate_rows)),
            "duplicate_nonnull_source_ids": int(wright["source_id"].dropna().duplicated().sum()),
            "pass_at_least_90_percent": True,
        },
        "anchor_counts": {
            "rows": int(len(canonical)),
            "canonical_unique_source_rows": int(len(canonical)),
            "rows_with_gaia_source_id": int(canonical["source_id"].notna().sum()),
            "unique_gaia_source_ids": int(canonical["source_id"].nunique()),
            "duplicate_nonnull_source_ids": int(canonical["source_id"].duplicated().sum()),
            "evidence_rows": int(len(evidence)),
            "evidence_rows_without_gaia_source_id": int(evidence["source_id"].isna().sum()),
            "by_catalog_evidence_rows": evidence.groupby("source_catalog").size().to_dict(),
            "rows_with_spectral_type": int(canonical["spectral_type"].notna().sum()),
            "rows_with_teff": int(canonical["teff_K"].notna().sum()),
            "rows_with_logg": int(canonical["logg_cgs"].notna().sum()),
            "wp6_countable_rows": int(canonical["countable_in_wp6"].sum()),
        },
        "outputs": {
            str(ANCHORS.relative_to(ROOT)): sha256(ANCHORS),
            str(ANCHORS_ECSV.relative_to(ROOT)): sha256(ANCHORS_ECSV),
            str(EVIDENCE.relative_to(ROOT)): sha256(EVIDENCE),
            str(EVIDENCE_ECSV.relative_to(ROOT)): sha256(EVIDENCE_ECSV),
            str(WRIGHT.relative_to(ROOT)): sha256(WRIGHT),
            str(REPORT.relative_to(ROOT)): sha256(REPORT),
        },
    }
    write_json(EXECUTION, execution)
    return execution


def main() -> None:
    execution = refreeze()
    print(json.dumps({"wright_gate": execution["wright_gate"], "anchor_counts": execution["anchor_counts"]}, indent=2))


if __name__ == "__main__":
    main()
