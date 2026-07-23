#!/usr/bin/env python
"""Freeze WP1 pulsar, SNR, and radioactive-nuclide marker inputs."""

from __future__ import annotations

import hashlib
import json
import re
import tarfile
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from bs4 import BeautifulSoup


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "markers"
PROCESSED = ROOT / "data" / "processed"
PROVENANCE = ROOT / "provenance"
OUTPUT_MD = ROOT / "wp1_sn_markers.md"
ATNF_VERSION = "2.8.1"
ATNF_URL = (
    "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/"
    f"psrcat_pkg.v{ATNF_VERSION}.tar.gz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> requests.Response:
    response = requests.get(url, timeout=180)
    response.raise_for_status()
    path.write_bytes(response.content)
    return response


def parameter_values(block: str) -> dict[str, list[str]]:
    values: dict[str, list[str]] = {}
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        fields = line.split()
        if len(fields) >= 2:
            values.setdefault(fields[0], []).append(fields[1])
    return values


def first(values: dict[str, list[str]], key: str) -> str | None:
    return values.get(key, [None])[0]


def number(values: dict[str, list[str]], key: str) -> float | None:
    value = first(values, key)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_atnf(database: str) -> pd.DataFrame:
    rows = []
    for block in database.split("@"):
        values = parameter_values(block)
        name = first(values, "PSRJ")
        ra = first(values, "RAJ")
        dec = first(values, "DECJ")
        if not (name and ra and dec):
            continue
        try:
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
        except ValueError:
            continue
        galactic = coord.galactic
        p0 = number(values, "P0")
        p1 = number(values, "P1")
        f0 = number(values, "F0")
        f1 = number(values, "F1")
        if p0 is None and f0 not in (None, 0):
            p0 = 1.0 / f0
        if p1 is None and f0 not in (None, 0) and f1 is not None:
            p1 = -f1 / f0**2
        age_yr = None
        if p0 is not None and p1 is not None and p1 > 0:
            age_yr = p0 / (2.0 * p1) / (365.25 * 86400.0)
        rows.append(
            {
                "psrj": name,
                "raj": ra,
                "decj": dec,
                "ra_deg": coord.ra.deg,
                "dec_deg": coord.dec.deg,
                "l_deg": galactic.l.deg,
                "b_deg": galactic.b.deg,
                "p0_s": p0,
                "p1": p1,
                "characteristic_age_yr": age_yr,
                "pmra_masyr": number(values, "PMRA"),
                "pmdec_masyr": number(values, "PMDEC"),
                "parallax_mas": number(values, "PX"),
                "distance_association_kpc": number(values, "DIST_A"),
                "distance_association_min_kpc": number(values, "DIST_AMN"),
                "distance_association_max_kpc": number(values, "DIST_AMX"),
                "distance_ymw17_kpc": number(values, "DIST_DM"),
                "distance_ne2001_kpc": number(values, "DIST_DM1"),
                "binary_model": first(values, "BINARY"),
                "binary_companion": first(values, "BINCOMP"),
                "association": first(values, "ASSOC"),
                "raw_block": block.strip(),
            }
        )
    return pd.DataFrame(rows)


def green_coordinates(name: str) -> tuple[float | None, float | None]:
    match = re.search(r"G?(\d+(?:\.\d+)?)([+-]\d+(?:\.\d+)?)", name)
    if not match:
        return None, None
    return float(match.group(1)), float(match.group(2))


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    created = datetime.now(timezone.utc).isoformat()

    tar_path = RAW / f"psrcat_pkg.v{ATNF_VERSION}.tar.gz"
    download(ATNF_URL, tar_path)
    with tarfile.open(fileobj=BytesIO(tar_path.read_bytes()), mode="r:gz") as archive:
        members = [member for member in archive.getmembers() if member.name.endswith("/psrcat.db")]
        if len(members) != 1:
            raise RuntimeError("Could not identify exactly one psrcat.db in ATNF archive")
        handle = archive.extractfile(members[0])
        if handle is None:
            raise RuntimeError("Could not read psrcat.db from ATNF archive")
        database = handle.read().decode("utf-8")
    db_path = RAW / f"psrcat.v{ATNF_VERSION}.db"
    db_path.write_text(database, encoding="utf-8")
    atnf = parse_atnf(database)
    atnf_wide = atnf[
        atnf["l_deg"].between(72, 88) & atnf["b_deg"].between(-5, 8)
    ].copy()
    atnf_path = PROCESSED / "wp1_atnf_pulsars_wide.parquet"
    atnf_wide.to_parquet(atnf_path, index=False)

    green_table = Vizier(columns=["**"], row_limit=-1).get_catalogs("VII/297")[0]
    green_raw = RAW / "green2024_vii_297_snrs.ecsv"
    green_table.write(green_raw, format="ascii.ecsv", overwrite=True)
    green = green_table.to_pandas()
    coordinates = [green_coordinates(str(name)) for name in green["SNR"]]
    green["l_deg"] = [item[0] for item in coordinates]
    green["b_deg"] = [item[1] for item in coordinates]
    green_wide = green[
        green["l_deg"].between(72, 88) & green["b_deg"].between(-5, 8)
    ].copy()
    green_path = PROCESSED / "wp1_green_snrs_wide.parquet"
    green_wide.to_parquet(green_path, index=False)

    chandra_url = "https://snrcat.cfa.harvard.edu/ChandraSNR/snrcat_gal.html"
    chandra_response = requests.get(chandra_url, timeout=120)
    chandra_response.raise_for_status()
    chandra_html = RAW / "chandra_snrcat_galactic.html"
    chandra_html.write_bytes(chandra_response.content)
    soup = BeautifulSoup(chandra_response.content, "html.parser")
    chandra_rows = []
    for row in soup.find_all("tr"):
        cells = [cell.get_text(" ", strip=True) for cell in row.find_all(["th", "td"])]
        if not cells:
            continue
        lon, lat = green_coordinates(cells[0])
        if lon is not None and 72 <= lon <= 88 and -5 <= lat <= 8:
            chandra_rows.append({"snr": cells[0], "l_deg": lon, "b_deg": lat, "cells": " | ".join(cells)})
    chandra = pd.DataFrame(chandra_rows)
    chandra_path = PROCESSED / "wp1_chandra_snrcat_wide.parquet"
    chandra.to_parquet(chandra_path, index=False)

    manitoba_url = "https://snrcat.physics.umanitoba.ca/"
    manitoba_html = RAW / "manitoba_snrcat_access_snapshot.html"
    manitoba_error = None
    manitoba_status = None
    try:
        manitoba_response = requests.get(manitoba_url, timeout=20)
        manitoba_status = manitoba_response.status_code
        manitoba_content = manitoba_response.content
        manitoba_text = BeautifulSoup(manitoba_content, "html.parser").get_text(" ", strip=True)
        manitoba_usable = manitoba_response.ok and "Error querying the database" not in manitoba_text
    except requests.RequestException as exc:
        manitoba_error = f"{type(exc).__name__}: {exc}"
        manitoba_usable = False
        manitoba_content = (
            "<!doctype html><meta charset='utf-8'><title>SNRcat access failure</title>"
            f"<p>{manitoba_error}</p>"
        ).encode("utf-8")
    manitoba_html.write_bytes(manitoba_content)

    martin09_pdf = RAW / "martin2009_integral_cygnus_arxiv1001.1521.pdf"
    martin10_pdf = RAW / "martin2010_population_synthesis_arxiv1001.1522.pdf"
    download("https://arxiv.org/pdf/1001.1521", martin09_pdf)
    download("https://arxiv.org/pdf/1001.1522", martin10_pdf)

    anchor = atnf_wide.loc[atnf_wide["psrj"].eq("J2032+4127")]
    if len(anchor) != 1:
        raise RuntimeError(f"Expected one PSR J2032+4127 entry, found {len(anchor)}")
    anchor = anchor.iloc[0]
    gamma = green_wide.loc[
        np.isclose(green_wide["l_deg"], 78.2)
        & np.isclose(green_wide["b_deg"], 2.1)
    ]
    if len(gamma) != 1:
        raise RuntimeError(f"Expected one Green gamma-Cygni row, found {len(gamma)}")
    gamma = gamma.iloc[0]
    gamma_minor = gamma.get("MinDiam")
    gamma_size = f"{gamma.get('MajDiam', '')} arcmin"
    if pd.notna(gamma_minor):
        gamma_size = f"{gamma.get('MajDiam', '')} x {gamma_minor} arcmin"
    header_count_match = re.search(r"# Number of pulsars is (\d+)", database)
    header_count = int(header_count_match.group(1)) if header_count_match else None

    lines = [
        "# WP1 supernova markers", "", f"Frozen: {created}", "",
        "## Pulsars: ATNF PSRCAT", "",
        f"- Frozen release: ATNF Pulsar Catalogue v{ATNF_VERSION}; the exact public archive and database are stored under `data/raw/markers/`.",
        f"- Wide box l=72-88 deg, b=-5-8 deg: {len(atnf_wide)} pulsars in `data/processed/wp1_atnf_pulsars_wide.parquet`.",
        "- PSR J2032+4127 / MT91 213 anchor:",
        f"  - P = {anchor['p0_s']:.12g} s; P-dot = {anchor['p1']:.6g}; characteristic age P/(2 P-dot) = {anchor['characteristic_age_yr']/1000:.1f} kyr.",
        f"  - Proper motion: mu_RA* = {anchor['pmra_masyr']} mas/yr, mu_Dec = {anchor['pmdec_masyr']} mas/yr; ATNF parallax = {anchor['parallax_mas']} mas.",
        f"  - Association distance field = {anchor['distance_association_kpc']} kpc; catalogue association string: `{anchor['association']}`.",
        "  - The characteristic age is not an explosion-age measurement; braking index, birth period, and binary timing systematics must be carried in WP8.", "",
        "## Supernova remnants", "",
        "- Primary census: Green's 2024 October catalogue, VizieR VII/297 (310 confirmed Galactic SNRs). "
        f"The wide box contains {len(green_wide)} entries, frozen in `data/processed/wp1_green_snrs_wide.parquet`.",
        f"- Gamma Cygni is Green {gamma['SNR']} ({gamma.get('Names', '')}), type `{gamma.get('type', '')}`, angular size {gamma_size}.",
        "- Physical anchor for Gamma Cygni: Leahy, Green & Ranasinghe (2013, MNRAS 436, 968; bibcode `2013MNRAS.436..968L`) infer d=1.7-2.6 kpc from H I absorption and a Sedov age of 6.8-10 kyr. Keep the older approximately 1.5 kpc class of estimates as an explicit literature branch; association with Cyg OB2 is unsettled.",
        f"- Chandra SNRcat snapshot: {len(chandra)} entries in the same wide box. The Ferrand/Safi-Harb Manitoba SNRcat endpoint was {'usable' if manitoba_usable else 'not usable (the server returned its database-query error)'} on the freeze date; its raw response is retained rather than silently replacing missing high-energy fields.", "",
        "## INTEGRAL radioactive-nuclide marker", "",
        "- Martin et al. (2009/2010), A&A 506, 703 (arXiv:1001.1521; DOI `10.1051/0004-6361/200912178`): total Cygnus-region 1809 keV flux = (6.0 +/- 1.0) x 10^-5 ph cm^-2 s^-1; component attributed to the Cygnus complex = (3.9 +/- 1.1) x 10^-5 ph cm^-2 s^-1; 60Fe upper limit = 1.6 x 10^-5 ph cm^-2 s^-1. The inferred 26Al morphology is centred near Cyg OB2 and extends about 9 deg or more.",
        "- Martin et al. (2010), A&A 511, A86 (arXiv:1001.1522; DOI `10.1051/0004-6361/200913385`) supplies the population-synthesis comparison. Both exact PDFs are frozen under `data/raw/markers/`.",
        "- WP8 use: 26Al is a combined winds-plus-SN consistency constraint with a mean lifetime near 1 Myr, not a direct count of recent SNe and not a Cyg OB2-only measurement.", "",
    ]
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")

    outputs = [
        tar_path, db_path, atnf_path, green_raw, green_path, chandra_html,
        chandra_path, manitoba_html, martin09_pdf, martin10_pdf, OUTPUT_MD,
    ]
    record = {
        "created_utc": created,
        "script": "scripts/wp1_freeze_sn_markers.py",
        "selection": {"l_deg": [72.0, 88.0], "b_deg": [-5.0, 8.0]},
        "atnf": {
            "version": ATNF_VERSION,
            "database_header_rows": header_count,
            "parsed_rows_with_coordinates": len(atnf),
            "wide_rows": len(atnf_wide),
            "psr_j2032_4127_rows": 1,
        },
        "green": {"version": "2024 October / VizieR VII/297", "all_rows": len(green), "wide_rows": len(green_wide)},
        "snrcat": {"chandra_wide_rows": len(chandra), "manitoba_endpoint_usable": manitoba_usable, "manitoba_http_status": manitoba_status, "manitoba_error": manitoba_error},
        "integral_26al": {"total_flux_1e-5_ph_cm-2_s-1": [6.0, 1.0], "complex_flux_1e-5_ph_cm-2_s-1": [3.9, 1.1], "fe60_upper_limit_1e-5_ph_cm-2_s-1": 1.6},
        "outputs": {str(path.relative_to(ROOT)): sha256(path) for path in outputs},
    }
    execution = PROVENANCE / "wp1_sn_markers_execution.json"
    temporary = execution.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(execution)
    print(json.dumps({"atnf_wide": len(atnf_wide), "green_wide": len(green_wide), "chandra_wide": len(chandra), "manitoba_usable": manitoba_usable}, indent=2))


if __name__ == "__main__":
    main()
