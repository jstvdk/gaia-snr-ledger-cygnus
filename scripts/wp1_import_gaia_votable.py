#!/usr/bin/env python
"""Import a directly streamed completed Gaia TAP VOTable into frozen formats."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from astropy.table import Table

from gaia_download import fits_compatible_table


ROOT = Path(__file__).resolve().parents[1]


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    args = parser.parse_args()
    name = args.name
    if name not in {"wide_tile01", "wide_tile04"}:
        raise ValueError("This importer is restricted to the two directly streamed wide tiles")

    raw = ROOT / "data" / "raw" / "gaia"
    provenance = ROOT / "provenance"
    votable_path = raw / f"wp1_gaia_{name}.vot"
    fits_path = raw / f"wp1_gaia_{name}.fits"
    parquet_path = raw / f"wp1_gaia_{name}.parquet"
    query_path = ROOT / "queries" / f"gaia_{name}.adql"
    state_path = provenance / f"wp1_gaia_{name}_job.json"
    execution_path = provenance / f"wp1_gaia_{name}_execution.json"

    state = json.loads(state_path.read_text(encoding="utf-8"))
    query_hash = sha256(query_path)
    if state.get("query_sha256") != query_hash:
        raise RuntimeError("Saved TAP job and query hashes differ")
    table = Table.read(votable_path, format="votable")
    lowered = [name.lower() for name in table.colnames]
    if len(lowered) != len(set(lowered)):
        raise RuntimeError("Lower-casing VOTable columns would create duplicate names")
    table.rename_columns(table.colnames, lowered)
    table = fits_compatible_table(table)
    table.write(fits_path, format="fits", overwrite=True)
    table.write(parquet_path, format="parquet", overwrite=True)
    now = datetime.now(timezone.utc).isoformat()
    record = {
        "executed_utc": state["submitted_utc"],
        "imported_utc": now,
        "query_name": name,
        "tap_url": state["tap_url"],
        "tap_job_url": state["tap_job_url"],
        "tap_phase": "COMPLETED",
        "query_file": str(query_path.relative_to(ROOT)),
        "query_sha256": query_hash,
        "row_count": len(table),
        "column_count": len(table.colnames),
        "columns": list(table.colnames),
        "transport": {
            "method": "direct streamed GET of the completed TAP result after the pyvo in-memory transfer stalled",
            "file": str(votable_path.relative_to(ROOT)),
            "sha256": sha256(votable_path),
            "bytes": votable_path.stat().st_size,
        },
        "outputs": {
            str(fits_path.relative_to(ROOT)): sha256(fits_path),
            str(parquet_path.relative_to(ROOT)): sha256(parquet_path),
        },
    }
    write_json(execution_path, record)
    state.update(
        {
            "status": "COMPLETED",
            "last_phase": "COMPLETED",
            "completed_utc": now,
            "row_count": len(table),
            "column_count": len(table.colnames),
            "transport": "direct_votable_stream",
        }
    )
    write_json(state_path, state)
    print(json.dumps({"name": name, "rows": len(table), "outputs": record["outputs"]}, indent=2))


if __name__ == "__main__":
    main()
