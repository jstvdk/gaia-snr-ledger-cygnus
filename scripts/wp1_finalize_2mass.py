#!/usr/bin/env python
"""Validate the official Gaia DR3/2MASS cross-match and freeze a narrow join."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "gaia" / "wp1_gaia_2mass_join.parquet"
NARROW = ROOT / "data" / "processed" / "wp1_gaia_narrow.parquet"
OUTPUT = ROOT / "data" / "processed" / "wp1_2mass_join.parquet"
OUTPUT_FITS = ROOT / "data" / "processed" / "wp1_2mass_join.fits"
EXECUTION = ROOT / "provenance" / "wp1_2mass_join_execution.json"
README = ROOT / "provenance" / "wp1_2mass_join_README.md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    narrow = pd.read_parquet(NARROW, columns=["source_id"])
    raw = pd.read_parquet(RAW)
    narrow["source_id"] = narrow["source_id"].astype("int64")
    raw["source_id"] = raw["source_id"].astype("int64")
    duplicates = int(raw["source_id"].duplicated().sum())
    if duplicates:
        examples = raw.loc[raw["source_id"].duplicated(False), "source_id"].head().tolist()
        raise RuntimeError(f"2MASS query has {duplicates} duplicate source_ids; examples={examples}")
    narrow_ids = set(narrow["source_id"])
    outside_narrow = sorted(set(raw["source_id"]) - narrow_ids)
    joined = narrow.merge(raw, on="source_id", how="left", validate="one_to_one")
    joined["has_2mass_psc_match"] = joined["tmass_designation"].notna()
    joined["has_j"] = joined["j_m"].notna() & joined["j_msigcom"].notna()
    joined["has_h"] = joined["h_m"].notna() & joined["h_msigcom"].notna()
    joined["has_ks"] = joined["ks_m"].notna() & joined["ks_msigcom"].notna()
    joined["has_complete_jhk"] = joined[["has_j", "has_h", "has_ks"]].all(axis=1)
    quality = joined["ph_qual"].fillna("").astype(str)
    joined["ph_qual_aaa"] = quality.eq("AAA")

    joined.to_parquet(OUTPUT, index=False)
    fits_frame = joined.copy()
    for column in fits_frame.select_dtypes(include=["string", "object"]).columns:
        fits_frame[column] = fits_frame[column].astype("string").fillna("")
    Table.from_pandas(fits_frame).write(OUTPUT_FITS, format="fits", overwrite=True)

    rows = len(joined)
    psc_matches = int(joined["has_2mass_psc_match"].sum())
    complete = int(joined["has_complete_jhk"].sum())
    aaa = int(joined["ph_qual_aaa"].sum())
    created = datetime.now(timezone.utc).isoformat()
    record = {
        "created_utc": created,
        "script": "scripts/wp1_finalize_2mass.py",
        "upstream": {
            str(NARROW.relative_to(ROOT)): sha256(NARROW),
            str(RAW.relative_to(ROOT)): sha256(RAW),
            "provenance/wp1_gaia_2mass_join_execution.json": sha256(
                ROOT / "provenance" / "wp1_gaia_2mass_join_execution.json"
            ),
        },
        "counts": {
            "gaia_narrow_rows": rows,
            "official_2mass_psc_matches": psc_matches,
            "match_fraction": psc_matches / rows,
            "complete_jhk_with_finite_uncertainties": complete,
            "complete_jhk_fraction": complete / rows,
            "ph_qual_AAA": aaa,
            "ph_qual_AAA_fraction": aaa / rows,
            "raw_duplicate_source_ids": duplicates,
            "raw_rows": len(raw),
            "raw_rows_outside_frozen_narrow": len(outside_narrow),
            "raw_source_ids_outside_frozen_narrow": outside_narrow,
            "unmatched_rows_retained_with_null_photometry": rows - psc_matches,
        },
        "selection": {
            "l_deg": [77.0, 83.0],
            "b_deg": [-1.5, 4.0],
            "parallax_mas": [0.35, 1.10],
            "phot_g_mean_mag_lt": 19.0,
        },
        "notes": [
            "Uses gaiadr3.tmass_psc_xsc_best_neighbour and gaiadr3.tmass_psc_xsc_join.",
            "J/H/Ks values and uncertainties come from gaiadr1.tmass_original_valid, the archived 2MASS PSC table.",
            "All Gaia narrow rows are retained; missing/upper-limit photometry is represented by null values and explicit flags.",
            "One Gaia stored-l boundary row lies 0.00000572 deg outside the Astropy-defined frozen narrow catalogue and is explicitly excluded by source_id.",
            "XSC-only counterparts are not stellar point-source photometry and are not used by WP3.",
        ],
        "outputs": {
            str(OUTPUT.relative_to(ROOT)): {"sha256": sha256(OUTPUT), "bytes": OUTPUT.stat().st_size},
            str(OUTPUT_FITS.relative_to(ROOT)): {"sha256": sha256(OUTPUT_FITS), "bytes": OUTPUT_FITS.stat().st_size},
        },
    }
    temporary = EXECUTION.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(EXECUTION)

    README.write_text(
        "\n".join(
            [
                "# wp1_2mass_join", "",
                "Content: one row per frozen `wp1_gaia_narrow` source, with the official Gaia DR3 best-neighbour 2MASS PSC identifier, match diagnostics, J/H/Ks photometry and uncertainties, and explicit availability/quality flags.", "",
                "Upstream artifacts: `data/processed/wp1_gaia_narrow.parquet`, `queries/gaia_2mass_join.adql`, and `data/raw/gaia/wp1_gaia_2mass_join.parquet`.", "",
                f"Frozen: {created}", "",
                f"Rows: {rows}; PSC matches: {psc_matches} ({psc_matches/rows:.2%}); complete JHK: {complete} ({complete/rows:.2%}); AAA: {aaa} ({aaa/rows:.2%}).", "",
                "Use the Parquet artifact for analysis. The FITS copy is supplied for interoperability.", "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(record["counts"], indent=2))


if __name__ == "__main__":
    main()
