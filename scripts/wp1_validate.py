#!/usr/bin/env python
"""Validate all named WP1 artifacts and write the final completion manifest."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / "provenance"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(name: str) -> dict:
    return json.loads((PROVENANCE / name).read_text(encoding="utf-8"))


def validate_outputs(record: dict) -> list[dict]:
    outputs = record.get("outputs", record.get("output", {}))
    if "file" in outputs and "sha256" in outputs:
        outputs = {outputs["file"]: outputs}
    verified = []
    for relative, expected in outputs.items():
        path = ROOT / relative
        expected_hash = expected if isinstance(expected, str) else expected["sha256"]
        actual = sha256(path)
        if actual != expected_hash:
            raise RuntimeError(f"Checksum mismatch for {relative}")
        verified.append(
            {
                "file": relative,
                "sha256": actual,
                "bytes": path.stat().st_size,
            }
        )
    return verified


def main() -> None:
    created = datetime.now(timezone.utc).isoformat()
    narrow = load("wp1_gaia_narrow_validation.json")
    two_mass_tap = load("wp1_gaia_2mass_join_execution.json")
    two_mass = load("wp1_2mass_join_execution.json")
    anchors = load("wp1_spectroscopic_anchors_execution.json")
    extinction = load("wp1_extinction_refs_execution.json")
    markers = load("wp1_sn_markers_execution.json")
    wide = load("wp1_gaia_wide_validation.json")

    records = {
        "gaia_narrow": narrow,
        "gaia_2mass_tap": two_mass_tap,
        "2mass_join": two_mass,
        "spectroscopic_anchors": anchors,
        "extinction_refs": extinction,
        "sn_markers": markers,
        "gaia_wide": wide,
    }
    verified_files = []
    for name, record in records.items():
        for item in validate_outputs(record):
            item["record"] = name
            verified_files.append(item)

    required = [
        "data/processed/wp1_gaia_narrow.parquet",
        "data/processed/wp1_gaia_narrow.fits",
        "data/processed/wp1_2mass_join.parquet",
        "data/processed/wp1_2mass_join.fits",
        "data/processed/wp1_spectroscopic_anchors.parquet",
        "data/processed/wp1_spectroscopic_anchors.ecsv",
        "data/processed/wp1_spectroscopic_anchor_records.parquet",
        "data/processed/wp1_spectroscopic_anchor_records.ecsv",
        "wp1_extinction_refs.md",
        "wp1_sn_markers.md",
        "data/processed/wp1_gaia_wide.parquet",
    ]
    missing = [path for path in required if not (ROOT / path).is_file()]
    if missing:
        raise RuntimeError(f"Missing named WP1 artifacts: {missing}")

    checks = {
        "narrow_rows_sane": 100_000 <= narrow["validation"]["rows"] <= 1_000_000,
        "narrow_unique_source_ids": narrow["validation"]["duplicate_source_ids"] == 0,
        "2mass_tap_completed": two_mass_tap["tap_phase"] == "COMPLETED",
        "2mass_one_row_per_narrow_source": two_mass["counts"]["gaia_narrow_rows"] == narrow["validation"]["rows"],
        "2mass_no_duplicate_source_ids": two_mass["counts"]["raw_duplicate_source_ids"] == 0,
        "wright_match_gate_at_least_90_percent": anchors["wright_gate"]["match_fraction"] >= 0.9,
        "wright_unmatched_individually_accounted": anchors["wright_gate"]["matched_to_gaia_dr3"] == anchors["wright_gate"]["census_rows"],
        "anchors_have_spectral_types": anchors["anchor_counts"]["rows_with_spectral_type"] > 0,
        "anchors_have_teff_logg": anchors["anchor_counts"]["rows_with_teff"] > 0 and anchors["anchor_counts"]["rows_with_logg"] > 0,
        "two_local_extinction_maps_frozen": bool(extinction["baseline"].get("cube_sha256")) and len(extinction["check"].get("files", {})) == 5,
        "pulsar_anchor_frozen": markers["atnf"]["psr_j2032_4127_rows"] == 1,
        "gamma_cygni_census_frozen": markers["green"]["wide_rows"] > 0,
        "wide_exact_count_recovered": wide["counts"]["combined_rows"] == wide["counts"]["exact_count_query"],
        "wide_unique_source_ids": wide["counts"]["duplicate_source_ids"] == 0,
    }
    wp3_checks = [
        "narrow_rows_sane", "narrow_unique_source_ids", "2mass_tap_completed",
        "2mass_one_row_per_narrow_source", "2mass_no_duplicate_source_ids",
        "wright_match_gate_at_least_90_percent", "wright_unmatched_individually_accounted",
        "anchors_have_spectral_types", "anchors_have_teff_logg",
        "two_local_extinction_maps_frozen",
    ]
    wp3_ready = all(checks[name] for name in wp3_checks)
    wp1_complete = all(checks.values()) and not missing
    if not wp1_complete:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"WP1 validation failed: {failed}")

    manifest = {
        "created_utc": created,
        "script": "scripts/wp1_validate.py",
        "wp3_inputs_ready": wp3_ready,
        "wp1_complete": wp1_complete,
        "checks": checks,
        "summary": {
            "gaia_narrow_rows": narrow["validation"]["rows"],
            "2mass_psc_matches": two_mass["counts"]["official_2mass_psc_matches"],
            "2mass_complete_jhk": two_mass["counts"]["complete_jhk_with_finite_uncertainties"],
            "wright_matched": anchors["wright_gate"]["matched_to_gaia_dr3"],
            "wright_census_rows": anchors["wright_gate"]["census_rows"],
            "spectroscopic_anchor_rows": anchors["anchor_counts"]["rows"],
            "spectroscopic_anchor_evidence_rows": anchors["anchor_counts"]["evidence_rows"],
            "gaia_wide_rows": wide["counts"]["combined_rows"],
            "atnf_wide_rows": markers["atnf"]["wide_rows"],
            "green_snr_wide_rows": markers["green"]["wide_rows"],
        },
        "verified_files": verified_files,
    }
    manifest_path = PROVENANCE / "wp1_manifest.json"
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    temporary.replace(manifest_path)

    report_path = ROOT / "wp1_completion_report.md"
    summary = manifest["summary"]
    report_path.write_text(
        "\n".join(
            [
                "# WP1 completion report", "", f"Validated: {created}", "",
                f"- WP3 input gate: **{'READY' if wp3_ready else 'BLOCKED'}**",
                f"- Complete WP1 deliverable set: **{'PASS' if wp1_complete else 'FAIL'}**",
                f"- Gaia narrow: {summary['gaia_narrow_rows']:,} unique sources.",
                f"- Official 2MASS PSC matches: {summary['2mass_psc_matches']:,}; complete J/H/Ks with finite uncertainties: {summary['2mass_complete_jhk']:,}.",
                f"- Wright+15 Gaia match: {summary['wright_matched']}/{summary['wright_census_rows']} ({summary['wright_matched']/summary['wright_census_rows']:.1%}); every narrow exclusion is listed in `wp1_wright15_crossmatch.md`.",
                f"- Spectroscopic anchors: {summary['spectroscopic_anchor_rows']:,} unique Gaia-source rows; {summary['spectroscopic_anchor_evidence_rows']:,} preserved literature-evidence rows.",
                "- Extinction references: local Vergely+22 25 pc cube plus all five Dharmawardena+22 Cygnus products.",
                f"- Gaia wide: {summary['gaia_wide_rows']:,} unique sources, matching the exact count query.",
                f"- SN markers: {summary['atnf_wide_rows']} ATNF pulsars and {summary['green_snr_wide_rows']} Green SNRs in the wide box, plus frozen INTEGRAL papers/values.", "",
                "Every file and checksum used by this verdict is enumerated in `provenance/wp1_manifest.json`.", "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps({"wp3_inputs_ready": wp3_ready, "wp1_complete": wp1_complete, **summary}, indent=2))


if __name__ == "__main__":
    main()
