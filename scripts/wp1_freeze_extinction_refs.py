#!/usr/bin/env python
"""Freeze WP1 baseline/check 3D-extinction references and verify coverage."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import dustmaps
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from dustmaps.bayestar import BayestarWebQuery


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "extinction"
PROVENANCE = ROOT / "provenance"
OUTPUT_MD = ROOT / "wp1_extinction_refs.md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists():
        return
    temporary = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=300) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_content(1024 * 1024):
                if block:
                    handle.write(block)
    temporary.replace(path)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    created = datetime.now(timezone.utc).isoformat()

    longitudes = np.array([77.0, 80.0, 83.0])
    latitudes = np.array([-1.5, 1.0, 4.0])
    distances = np.array([1.0, 1.35, 1.60, 2.0])
    grid = np.array(
        [(lon, lat, distance) for lon in longitudes for lat in latitudes for distance in distances]
    )
    coords = SkyCoord(
        l=grid[:, 0] * u.deg,
        b=grid[:, 1] * u.deg,
        distance=grid[:, 2] * u.kpc,
        frame="galactic",
    )
    bayestar_error = None
    try:
        bayestar = BayestarWebQuery(version="bayestar2019")
        median, flags = bayestar(coords, mode="median", return_flags=True)
        percentiles = bayestar(coords, mode="percentile", pct=[16, 84])
        coverage = pd.DataFrame(
            {
                "l_deg": grid[:, 0], "b_deg": grid[:, 1],
                "distance_kpc": grid[:, 2], "bayestar19_median": median,
                "bayestar19_p16": percentiles[:, 0],
                "bayestar19_p84": percentiles[:, 1],
                "converged": flags["converged"],
                "reliable_dist": flags["reliable_dist"],
            }
        )
    except Exception as exc:
        bayestar_error = f"{type(exc).__name__}: {exc}"
        coverage = pd.DataFrame(
            {
                "l_deg": grid[:, 0], "b_deg": grid[:, 1],
                "distance_kpc": grid[:, 2], "bayestar19_median": np.nan,
                "bayestar19_p16": np.nan, "bayestar19_p84": np.nan,
                "converged": False, "reliable_dist": False,
            }
        )
    coverage_path = RAW / "wp1_bayestar19_coverage_grid.parquet"
    coverage.to_parquet(coverage_path, index=False)

    vergely_tables = Vizier(columns=["**"], row_limit=-1).get_catalogs(
        "J/A+A/664/A174/list"
    )
    if len(vergely_tables) != 1:
        raise RuntimeError("Expected exactly one Vergely+22 map-list table")
    vergely_list = vergely_tables[0]
    vergely_path = RAW / "vergely2022_extinction_cube_list.ecsv"
    vergely_list.write(vergely_path, format="ascii.ecsv", overwrite=True)

    vergely_cube_name = "explore_cube_density_values_025pc_v2.fits"
    vergely_cube = RAW / vergely_cube_name
    download(
        "https://cdsarc.cds.unistra.fr/ftp/cats/J/A+A/664/A174/fits/"
        + vergely_cube_name,
        vergely_cube,
    )

    dharma_base = "https://cdsarc.cds.unistra.fr/ftp/cats/J/A+A/658/A166/Cygnus/"
    dharma_files = [
        "gpy_dens_median.out", "ext_med_cube.pkl.npy", "l_bounds_pred.pkl.npy",
        "b_bounds_pred.pkl.npy", "d_bounds_pred.pkl.npy",
    ]
    dharma_paths = []
    for filename in dharma_files:
        path = RAW / f"dharmawardena2022_cygnus_{filename}"
        download(dharma_base + filename, path)
        dharma_paths.append(path)

    ebv_factor = 0.98 * 0.901
    lines = [
        "# WP1 extinction references", "",
        f"Frozen: {created}", "",
        "## Baseline: Vergely, Lallement & Cox (2022)", "",
        "- Product: Vergely, Lallement & Cox (2022), A&A 664, A174 (bibcode `2022A&A...664A.174V`; DOI `10.1051/0004-6361/202243319`).",
        "- Frozen archive: VizieR `J/A+A/664/A174`, catalogue DOI `10.26093/cds/vizier.36640174`, last catalogue modification reported as 2024-02-15. The five official cube choices are frozen in `data/raw/extinction/vergely2022_extinction_cube_list.ecsv`.",
        f"- Frozen WP3 cube: `data/raw/extinction/{vergely_cube_name}` ({vergely_cube.stat().st_size} bytes; SHA-256 `{sha256(vergely_cube)}`). It spans 6 x 6 x 0.8 kpc with 10 pc sampling and approximately 25 pc correlation length, so Cyg OB2 at about 1.35-1.6 kpc is inside the volume. The 5 kpc/20 pc cube is a fallback sensitivity check only.",
        "- Units: the cube is differential monochromatic extinction density at 550 nm in nanomagnitude per parsec. Integrate along the Sun-to-star ray to the sampled distance. Compare the resulting A0 with fitted A_V while retaining the paper's warning that A0 and A_V are close, not identical.", "",
        "## Check map: Dharmawardena et al. (2022) Cygnus X", "",
        "- Product: Dharmawardena et al. (2022), A&A 658, A166 (bibcode `2022A&A...658A.166D`; DOI `10.1051/0004-6361/202141298`), a non-negative Gaussian-process 3D dust reconstruction made specifically for Cygnus X.",
        "- Frozen archive: all five Cygnus products from VizieR `J/A+A/658/A166/Cygnus` (median density, median cumulative extinction, and l/b/d grid boundaries). Their hashes are in `provenance/wp1_extinction_refs_execution.json`.",
        "- Use the provided boundary arrays rather than assuming pixel centres. Compare its cumulative extinction to the integrated Vergely A0 and the WP3 fitted A_V; retain the different map methodology as a systematic branch.", "",
        "## Optional third diagnostic: Bayestar19", "",
        "- Product: Green et al. (2019), ApJ 887, 93 (bibcode `2019ApJ...887...93G`; DOI `10.3847/1538-4357/ab5362`; data DOI `10.7910/DVN/2EJ9TX`).",
        f"- Native-to-E(B-V) diagnostic conversion: `E(B-V) = 0.98 x 0.901 x E_B19 = {ebv_factor:.6f} x E_B19`.",
        f"- Access check at freeze time: {'successful' if bayestar_error is None else 'failed without blocking WP3 because two local maps are frozen: ' + bayestar_error}.",
        f"- The {len(coverage)}-point attempted coverage grid is retained in `data/raw/extinction/wp1_bayestar19_coverage_grid.parquet`; successful rows with reliable-distance flags: {int(coverage['reliable_dist'].sum())}.", "",
        "## Binding use in WP3 step 3b", "",
        "1. Query both maps at the same sampled member distance; never query only to infinity.",
        "2. Compare the Vergely and Dharmawardena maps with fitted extinction spatially and as residual distributions per subgroup; do not replace per-star Gaia+2MASS/spectral-type fits with either map.",
        "3. Treat angular-resolution mismatch explicitly: the maps cannot reproduce arcminute-scale cloud structure in the Cyg OB2 core.",
        "4. These maps are external pipelines but not photometrically independent of Gaia/2MASS. Agreement is a consistency check, not independent confirmation.",
        "5. Carry native quality flags and map choice as a labeled systematic branch; do not average the maps.", "",
    ]
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")

    record = {
        "created_utc": created,
        "script": "scripts/wp1_freeze_extinction_refs.py",
        "dustmaps_version": getattr(dustmaps, "__version__", "unknown"),
        "baseline": {
            "name": "Vergely et al. 2022 Gaia-2MASS map",
            "vizier_catalog": "J/A+A/664/A174",
            "cube_list_rows": len(vergely_list),
            "selected_resolution_pc": 25,
            "selected_sampling_pc": 10,
            "selected_extent_kpc": 3,
            "cube_file": str(vergely_cube.relative_to(ROOT)),
            "cube_sha256": sha256(vergely_cube),
        },
        "check": {
            "name": "Dharmawardena et al. 2022 Cygnus X map",
            "vizier_catalog": "J/A+A/658/A166/Cygnus",
            "files": {
                str(path.relative_to(ROOT)): sha256(path) for path in dharma_paths
            },
        },
        "optional_bayestar19": {
            "query_status": "success" if bayestar_error is None else "failed",
            "error": bayestar_error,
            "coverage_grid_rows": len(coverage),
            "converged": int(coverage["converged"].sum()),
            "reliable_dist": int(coverage["reliable_dist"].sum()),
        },
        "outputs": {
            str(OUTPUT_MD.relative_to(ROOT)): sha256(OUTPUT_MD),
            str(coverage_path.relative_to(ROOT)): sha256(coverage_path),
            str(vergely_path.relative_to(ROOT)): sha256(vergely_path),
            str(vergely_cube.relative_to(ROOT)): sha256(vergely_cube),
            **{str(path.relative_to(ROOT)): sha256(path) for path in dharma_paths},
        },
    }
    execution_path = PROVENANCE / "wp1_extinction_refs_execution.json"
    temporary = execution_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(execution_path)
    print(json.dumps(record["baseline"], indent=2))
    print(json.dumps(record["check"], indent=2))


if __name__ == "__main__":
    main()
