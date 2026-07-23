#!/usr/bin/env python
"""Combine and validate the four exact WP1 Gaia wide-box tiles."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "gaia"
PROCESSED = ROOT / "data" / "processed"
PROVENANCE = ROOT / "provenance"
OUTPUT = PROCESSED / "wp1_gaia_wide.parquet"
EXPECTED_COUNT_FILE = RAW / "wp1_gaia_wide_exact_count.parquet"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    tiles = []
    inputs = []
    for index in range(1, 5):
        name = f"wide_tile{index:02d}"
        path = RAW / f"wp1_gaia_{name}.parquet"
        execution_path = PROVENANCE / f"wp1_gaia_{name}_execution.json"
        execution = json.loads(execution_path.read_text(encoding="utf-8"))
        if execution.get("tap_phase") != "COMPLETED":
            raise RuntimeError(f"{name} is not COMPLETED")
        if execution["outputs"].get(str(path.relative_to(ROOT))) != sha256(path):
            raise RuntimeError(f"{name} Parquet checksum differs from execution record")
        table = pq.read_table(path)
        if table.num_rows != execution["row_count"]:
            raise RuntimeError(f"{name} row count differs from execution record")
        table = table.append_column("tile_id", pa.array([f"tile{index:02d}"] * table.num_rows))
        tiles.append(table)
        inputs.append(
            {
                "file": str(path.relative_to(ROOT)),
                "rows": table.num_rows,
                "sha256": sha256(path),
                "tap_job_url": execution["tap_job_url"],
                "query_sha256": execution["query_sha256"],
            }
        )

    combined = pa.concat_tables(tiles, promote_options="default")
    expected = int(pq.read_table(EXPECTED_COUNT_FILE)["row_count"][0].as_py())
    if combined.num_rows != expected:
        raise RuntimeError(
            f"Wide tiles sum to {combined.num_rows}, exact count query returned {expected}"
        )
    unique_source_ids = pc.count_distinct(combined["source_id"]).as_py()
    duplicates = combined.num_rows - unique_source_ids
    if duplicates:
        raise RuntimeError(f"Wide tiles contain {duplicates} duplicate source_ids")
    for column, low, high in [
        ("parallax", 0.35, 1.10),
    ]:
        values = combined[column].to_numpy(zero_copy_only=False)
        if np.nanmin(values) < low or np.nanmax(values) > high:
            raise RuntimeError(f"{column} values violate [{low}, {high}]")
    gmag = combined["phot_g_mean_mag"].to_numpy(zero_copy_only=False)
    if np.nanmax(gmag) >= 19.0:
        raise RuntimeError("phot_g_mean_mag contains values >=19")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    pq.write_table(combined, OUTPUT, compression="zstd", use_dictionary=True)
    created = datetime.now(timezone.utc).isoformat()
    record = {
        "created_utc": created,
        "script": "scripts/wp1_finalize_wide.py",
        "selection": {
            "l_deg": [72.0, 88.0],
            "b_deg": [-5.0, 8.0],
            "parallax_mas": [0.35, 1.10],
            "phot_g_mean_mag_lt": 19.0,
        },
        "input_tiles": inputs,
        "counts": {
            "exact_count_query": expected,
            "combined_rows": combined.num_rows,
            "unique_source_ids": unique_source_ids,
            "duplicate_source_ids": duplicates,
        },
        "output": {
            "file": str(OUTPUT.relative_to(ROOT)),
            "sha256": sha256(OUTPUT),
            "bytes": OUTPUT.stat().st_size,
            "columns": combined.column_names,
        },
    }
    execution_path = PROVENANCE / "wp1_gaia_wide_validation.json"
    temporary = execution_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(execution_path)
    readme = PROVENANCE / "wp1_gaia_wide_README.md"
    readme.write_text(
        "\n".join(
            [
                "# wp1_gaia_wide", "",
                "Content: exact Gaia DR3 wide box for WP6 runaway/traceback work, assembled from four non-overlapping Galactic-longitude tiles.", "",
                "Selection: l=72-88 deg, b=-5-8 deg, parallax=0.35-1.10 mas, G<19. Raw parallax is retained; zero-point correction remains a downstream operation.", "",
                f"Frozen: {created}; rows: {combined.num_rows}; unique source_ids: {unique_source_ids}; SHA-256: `{sha256(OUTPUT)}`.", "",
                "Format: Parquet is the canonical analysis artifact. A FITS duplicate is omitted because this >3 million-row table is a WP6 support catalogue and Parquet preserves the nullable columns with substantially less storage.", "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(record["counts"], indent=2))


if __name__ == "__main__":
    main()
