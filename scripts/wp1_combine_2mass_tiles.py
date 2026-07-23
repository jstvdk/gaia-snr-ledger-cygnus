#!/usr/bin/env python
"""Combine six exact official Gaia/2MASS cross-match tiles."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "gaia"
PROVENANCE = ROOT / "provenance"
OUTPUT = RAW / "wp1_gaia_2mass_join.parquet"
EXECUTION = PROVENANCE / "wp1_gaia_2mass_join_execution.json"
NARROW = ROOT / "data" / "processed" / "wp1_gaia_narrow.parquet"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    tables = []
    inputs = []
    for index in range(1, 7):
        name = f"2mass_tile{index:02d}"
        path = RAW / f"wp1_gaia_{name}.parquet"
        execution_path = PROVENANCE / f"wp1_gaia_{name}_execution.json"
        execution = json.loads(execution_path.read_text(encoding="utf-8"))
        if execution.get("tap_phase") != "COMPLETED":
            raise RuntimeError(f"{name} is not COMPLETED")
        relative = str(path.relative_to(ROOT))
        actual_hash = sha256(path)
        if execution["outputs"].get(relative) != actual_hash:
            raise RuntimeError(f"{name} Parquet checksum differs from execution record")
        table = pq.read_table(path)
        if table.num_rows != execution["row_count"]:
            raise RuntimeError(f"{name} row count differs from execution record")
        tables.append(table)
        inputs.append(
            {
                "query_name": name,
                "file": relative,
                "rows": table.num_rows,
                "sha256": actual_hash,
                "tap_job_url": execution["tap_job_url"],
                "query_sha256": execution["query_sha256"],
            }
        )

    combined = pa.concat_tables(tables, promote_options="default")
    source_ids = combined["source_id"].to_pylist()
    duplicate_source_ids = len(source_ids) - len(set(source_ids))
    if duplicate_source_ids:
        raise RuntimeError(f"2MASS tiles contain {duplicate_source_ids} duplicate source_ids")
    narrow_source_ids = set(pq.read_table(NARROW, columns=["source_id"])["source_id"].to_pylist())
    outside_frozen_narrow = sorted(set(source_ids) - narrow_source_ids)
    pq.write_table(combined, OUTPUT, compression="zstd", use_dictionary=True)
    created = datetime.now(timezone.utc).isoformat()
    record = {
        "created_utc": created,
        "script": "scripts/wp1_combine_2mass_tiles.py",
        "tap_phase": "COMPLETED",
        "acquisition": "six non-overlapping Gaia-l tiles after the monolithic archive job aborted its transaction",
        "selection": {
            "l_deg": [77.0, 83.0],
            "b_deg": [-1.5, 4.0],
            "parallax_mas": [0.35, 1.10],
            "phot_g_mean_mag_lt": 19.0,
        },
        "input_tiles": inputs,
        "row_count": combined.num_rows,
        "duplicate_source_ids": duplicate_source_ids,
        "rows_outside_frozen_narrow": len(outside_frozen_narrow),
        "source_ids_outside_frozen_narrow": outside_frozen_narrow,
        "boundary_note": (
            "Gaia stored-l selection includes source 2070682352787315200, whose ICRS "
            "position transforms with Astropy to l=83.00000571998599 deg; the canonical "
            "processed join keys to frozen narrow source_ids and therefore excludes it."
        ),
        "outputs": {
            str(OUTPUT.relative_to(ROOT)): {
                "sha256": sha256(OUTPUT),
                "bytes": OUTPUT.stat().st_size,
            }
        },
        "failed_monolithic_job": {
            "tap_job_url": "https://gea.esac.esa.int/tap-server/tap/async/eba92310-85b5-11f1-982e-bc97e148b76b-O",
            "phase": "ERROR",
            "archive_error": "org.postgresql.util.PSQLException: ERROR: current transaction is aborted, commands ignored until end of transaction block",
        },
    }
    temporary = EXECUTION.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(EXECUTION)
    print(
        json.dumps(
            {
                "rows": combined.num_rows,
                "duplicate_source_ids": duplicate_source_ids,
                "rows_outside_frozen_narrow": len(outside_frozen_narrow),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
