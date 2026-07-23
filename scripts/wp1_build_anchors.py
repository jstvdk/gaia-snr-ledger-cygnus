#!/usr/bin/env python
"""Freeze WP1 literature catalogues and build Gaia-linked spectroscopic anchors."""

from __future__ import annotations

import hashlib
import json
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pyvo
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from bs4 import BeautifulSoup
from urllib3.exceptions import InsecureRequestWarning


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "literature"
PROCESSED = ROOT / "data" / "processed"
PROVENANCE = ROOT / "provenance"
TAP_URL = "https://gea.esac.esa.int/tap-server/tap"
NARROW = ROOT / "data" / "processed" / "wp1_gaia_narrow.parquet"
MATCH_RADIUS_ARCSEC = 2.0


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


def scalar(value):
    if np.ma.is_masked(value) or value is None:
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def text_value(value) -> str | None:
    value = scalar(value)
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def float_value(value) -> float | None:
    value = scalar(value)
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def int_string(value) -> str | None:
    value = scalar(value)
    if value is None:
        return None
    try:
        return str(int(value))
    except (TypeError, ValueError):
        return None


def unique_headers(headers: list[str]) -> list[str]:
    result = []
    seen: dict[str, int] = {}
    for index, header in enumerate(headers):
        name = header.strip() or f"unnamed_{index}"
        count = seen.get(name, 0)
        seen[name] = count + 1
        result.append(name if count == 0 else f"{name}_{count + 1}")
    return result


def html_table(html: bytes) -> pd.DataFrame:
    table = BeautifulSoup(html, "html.parser").find("table")
    if table is None:
        raise RuntimeError("Expected an HTML table but none was present")
    rows = [
        [cell.get_text(" ", strip=True) for cell in row.find_all(["th", "td"])]
        for row in table.find_all("tr")
    ]
    headers = unique_headers(rows[0])
    body = [row for row in rows[1:] if len(row) == len(headers)]
    return pd.DataFrame(body, columns=headers)


def save_vizier_catalog(catalog: str, prefix: str) -> tuple[list[Table], list[dict]]:
    tables = list(Vizier(columns=["**"], row_limit=-1).get_catalogs(catalog))
    records = []
    for index, table in enumerate(tables, 1):
        table_name = str(table.meta.get("name", f"table{index}"))
        suffix = table_name.rsplit("/", 1)[-1]
        path = RAW / f"{prefix}_{suffix}.ecsv"
        table.write(path, format="ascii.ecsv", overwrite=True)
        records.append(
            {
                "catalog": catalog,
                "table": table_name,
                "rows": len(table),
                "file": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
        )
    return tables, records


def simbad_gaia_ids(names: list[str], raw_name: str) -> dict[str, str | None]:
    service = Simbad()
    service.add_votable_fields("ids")
    result = service.query_objects(names)
    if result is None:
        return {name: None for name in names}
    path = RAW / raw_name
    result.write(path, format="ascii.ecsv", overwrite=True)
    mapping = {name: None for name in names}
    for row in result:
        requested = text_value(row["user_specified_id"])
        match = re.search(r"(?:^|\|)Gaia DR3 (\d+)(?:\||$)", str(row["ids"]))
        if requested in mapping and match:
            mapping[requested] = match.group(1)
    return mapping


def tap_table(query: str) -> Table:
    return pyvo.dal.TAPService(TAP_URL).search(query).to_table()


def direct_coordinate_matches(
    coords: SkyCoord, labels: list[str], radius_arcsec: float = MATCH_RADIUS_ARCSEC
) -> dict[str, dict]:
    if not labels:
        return {}
    radius_deg = radius_arcsec / 3600.0
    conditions = " OR ".join(
        "1=CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{coord.ra.deg:.10f},{coord.dec.deg:.10f},{radius_deg:.10f}))"
        for coord in coords
    )
    query = (
        "SELECT source_id,ra,dec,l,b,phot_g_mean_mag,parallax,parallax_error,ruwe "
        "FROM gaiadr3.gaia_source WHERE " + conditions
    )
    candidates = tap_table(query)
    path = RAW / "wright2015_gaia_coordinate_candidates.ecsv"
    candidates.write(path, format="ascii.ecsv", overwrite=True)
    if len(candidates) == 0:
        return {}
    candidate_coords = SkyCoord(candidates["ra"], candidates["dec"])
    matches: dict[str, dict] = {}
    for label, coord in zip(labels, coords, strict=True):
        index, separation, _ = coord.match_to_catalog_sky(candidate_coords)
        index = int(np.asarray(index).item())
        separation_arcsec = float(np.asarray(separation.arcsec).item())
        if separation_arcsec <= radius_arcsec:
            matches[label] = {
                "source_id": int_string(candidates["source_id"][index]),
                "separation_arcsec": separation_arcsec,
            }
    return matches


def gaia_lookup(source_ids: list[str]) -> pd.DataFrame:
    unique_ids = sorted(set(source_ids))
    if not unique_ids:
        return pd.DataFrame()
    query = (
        "SELECT source_id,ra,dec,l,b,phot_g_mean_mag,parallax,parallax_error,ruwe "
        "FROM gaiadr3.gaia_source WHERE source_id IN (" + ",".join(unique_ids) + ")"
    )
    table = tap_table(query)
    path = RAW / "wright2015_gaia_source_lookup.ecsv"
    table.write(path, format="ascii.ecsv", overwrite=True)
    frame = table.to_pandas()
    frame["source_id"] = frame["source_id"].astype("int64").astype("string")
    return frame


def dr2_to_dr3(source_ids: list[str]) -> tuple[dict[str, str], Table]:
    unique_ids = sorted(set(source_ids))
    if not unique_ids:
        return {}, Table()
    query = (
        "SELECT dr2_source_id,dr3_source_id,angular_distance,magnitude_difference,"
        "proper_motion_propagation FROM gaiadr3.dr2_neighbourhood "
        "WHERE dr2_source_id IN (" + ",".join(unique_ids) + ")"
    )
    table = tap_table(query)
    frame = table.to_pandas()
    if frame.empty:
        return {}, table
    frame["dr2"] = frame["dr2_source_id"].astype("int64").astype("string")
    frame["dr3"] = frame["dr3_source_id"].astype("int64").astype("string")
    frame["abs_magdiff"] = frame["magnitude_difference"].abs().fillna(np.inf)
    frame = frame.sort_values(["dr2", "angular_distance", "abs_magdiff"])
    best = frame.drop_duplicates("dr2", keep="first")
    return dict(zip(best["dr2"], best["dr3"], strict=True)), table


def nearest_local_ids(
    ra_deg: list[float], dec_deg: list[float], local: pd.DataFrame
) -> tuple[list[str | None], list[float | None]]:
    local_coords = SkyCoord(local["ra"].to_numpy() * u.deg, local["dec"].to_numpy() * u.deg)
    coords = SkyCoord(np.asarray(ra_deg) * u.deg, np.asarray(dec_deg) * u.deg)
    indices, separations, _ = coords.match_to_catalog_sky(local_coords)
    ids: list[str | None] = []
    seps: list[float | None] = []
    for index, separation in zip(indices, separations.arcsec, strict=True):
        if separation <= MATCH_RADIUS_ARCSEC:
            ids.append(str(local.iloc[index]["source_id"]))
            seps.append(float(separation))
        else:
            ids.append(None)
            seps.append(float(separation))
    return ids, seps


def base_anchor(
    *, uid: str, catalog: str, version: str, record_id: str, name: str | None,
    source_id: str | None, ra: float | None, dec: float | None,
    spectral_type: str | None, teff: float | None = None,
    teff_error: float | None = None, logg: float | None = None,
    logg_error: float | None = None, av: float | None = None,
    group: str | None = None, source_table: str, bibcode: str,
    match_method: str | None = None, match_separation: float | None = None,
    notes: str | None = None,
) -> dict:
    return {
        "anchor_uid": uid,
        "source_catalog": catalog,
        "catalog_version": version,
        "catalog_record_id": record_id,
        "object_name": name,
        "source_id": source_id,
        "ra_deg": ra,
        "dec_deg": dec,
        "spectral_type": spectral_type,
        "teff_K": teff,
        "teff_error_K": teff_error,
        "logg_cgs": logg,
        "logg_error_cgs": logg_error,
        "extinction_av_mag": av,
        "subgroup_or_association": group,
        "source_table": source_table,
        "bibcode": bibcode,
        "gaia_match_method": match_method,
        "gaia_match_separation_arcsec": match_separation,
        "notes": notes,
    }


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    created = datetime.now(timezone.utc).isoformat()
    raw_records: list[dict] = []

    wright_tables, records = save_vizier_catalog("J/MNRAS/449/741", "wright2015")
    raw_records.extend(records)
    berlanas19_tables, records = save_vizier_catalog("J/MNRAS/484/1838", "berlanas2019")
    raw_records.extend(records)
    berlanas20_tables, records = save_vizier_catalog("J/A+A/642/A168", "berlanas2020")
    raw_records.extend(records)

    warnings.filterwarnings("ignore", category=InsecureRequestWarning)
    session = requests.Session()
    gosc_query_url = "https://gosc.cab.inta-csic.es/gosc-v3-query"
    query_page = session.get(gosc_query_url, timeout=60, verify=False)
    query_page.raise_for_status()
    query_soup = BeautifulSoup(query_page.content, "html.parser")
    attributes = [
        item.get("value") for item in query_soup.select('input[name="atributos[]"]')
    ]
    gosc_post = [("catalogo[]", "0")]
    gosc_post.extend(("atributos[]", attribute) for attribute in attributes)
    gosc_post.extend(
        [
            ("format", "html"), ("sel_crit1", "coor"),
            ("range1", "long"), ("range1min", "77"), ("range1max", "83"),
            ("range2", "lat"), ("range2min", "-1.5"), ("range2max", "4"),
            ("range3", "Bap"), ("range3min", "-99"), ("range3max", "99"),
            ("sort_var1", "long"), ("sort_var2", "lat"), ("path", "pub"),
        ]
    )
    gosc_response = session.post(
        "https://gosc.cab.inta-csic.es/gosc.php",
        data=gosc_post,
        timeout=180,
        verify=False,
    )
    gosc_response.raise_for_status()
    gosc_html = RAW / "gosc_v3_field_query.html"
    gosc_html.write_bytes(gosc_response.content)
    gosc = html_table(gosc_response.content)
    gosc = gosc[gosc["#"].str.fullmatch(r"\d+", na=False)].copy()
    gosc_csv = RAW / "gosc_v3_field_query.csv"
    gosc.to_csv(gosc_csv, index=False)
    raw_records.extend(
        [
            {
                "catalog": "GOSC v3 / GOSSS-DR2.2",
                "table": "main catalogue field query",
                "rows": len(gosc),
                "file": str(gosc_html.relative_to(ROOT)),
                "sha256": sha256(gosc_html),
                "tls_verification": False,
                "tls_note": "Site certificate chain is not trusted by the conda CA bundle; response frozen and hashed.",
            },
            {
                "catalog": "GOSC v3 / GOSSS-DR2.2",
                "table": "parsed field rows",
                "rows": len(gosc),
                "file": str(gosc_csv.relative_to(ROOT)),
                "sha256": sha256(gosc_csv),
            },
        ]
    )

    wr_url = "https://www.pacrowther.staff.shef.ac.uk/WRcat/index.php"
    wr_response = requests.get(wr_url, timeout=120, verify=False)
    wr_response.raise_for_status()
    wr_html = RAW / "galactic_wr_catalog_v1.33.html"
    wr_html.write_bytes(wr_response.content)
    wr_all = html_table(wr_response.content)
    wr_l = pd.to_numeric(wr_all["Galactic Longitude (deg)"], errors="coerce")
    wr_b = pd.to_numeric(wr_all["Galactic Latitude (deg)"], errors="coerce")
    wr = wr_all[wr_l.between(77, 83) & wr_b.between(-1.5, 4)].copy()
    wr_csv = RAW / "galactic_wr_catalog_v1.33_field.csv"
    wr.to_csv(wr_csv, index=False)
    raw_records.extend(
        [
            {
                "catalog": "Galactic Wolf Rayet Catalogue v1.33 (Gaia DR3), Aug 2025",
                "table": "full HTML snapshot",
                "rows": len(wr_all),
                "file": str(wr_html.relative_to(ROOT)),
                "sha256": sha256(wr_html),
                "tls_verification": False,
                "tls_note": "Site certificate chain is not trusted by the conda CA bundle; response frozen and hashed.",
            },
            {
                "catalog": "Galactic Wolf Rayet Catalogue v1.33 (Gaia DR3), Aug 2025",
                "table": "parsed field rows",
                "rows": len(wr),
                "file": str(wr_csv.relative_to(ROOT)),
                "sha256": sha256(wr_csv),
            },
        ]
    )

    local = pd.read_parquet(NARROW)
    local["source_id"] = local["source_id"].astype("int64").astype("string")
    local_ids = set(local["source_id"])

    wright = wright_tables[0]
    wright_names = [str(name).strip() for name in wright["SimbadName"]]
    wright_simbad = simbad_gaia_ids(wright_names, "wright2015_simbad_ids.ecsv")
    wright_coords = SkyCoord(
        wright["RAJ2000"], wright["DEJ2000"], unit=(u.hourangle, u.deg), frame="icrs"
    )
    missing_names = [name for name in wright_names if not wright_simbad.get(name)]
    missing_indices = [wright_names.index(name) for name in missing_names]
    coordinate_fallback = direct_coordinate_matches(
        wright_coords[missing_indices], missing_names
    )
    assigned_ids = [
        wright_simbad.get(name)
        or coordinate_fallback.get(name, {}).get("source_id")
        for name in wright_names
    ]
    unresolved = [name for name, source_id in zip(wright_names, assigned_ids) if source_id is None]
    if unresolved:
        raise RuntimeError(f"Unresolved Wright+15 Gaia matches: {unresolved}")
    outside_ids = [source_id for source_id in assigned_ids if source_id not in local_ids]
    outside = gaia_lookup(outside_ids)
    gaia = pd.concat(
        [
            local[["source_id", "ra", "dec", "l_deg", "b_deg", "phot_g_mean_mag", "parallax", "parallax_error", "ruwe"]]
            .rename(columns={"l_deg": "l", "b_deg": "b"}),
            outside,
        ],
        ignore_index=True,
    ).drop_duplicates("source_id", keep="first")
    gaia_by_id = gaia.set_index("source_id")

    wright_rows = []
    anchors: list[dict] = []
    for index, (row, name, source_id) in enumerate(
        zip(wright, wright_names, assigned_ids, strict=True)
    ):
        gaia_row = gaia_by_id.loc[source_id]
        gaia_coord = SkyCoord(float(gaia_row["ra"]) * u.deg, float(gaia_row["dec"]) * u.deg)
        separation = float(wright_coords[index].separation(gaia_coord).arcsec)
        in_narrow = source_id in local_ids
        reasons = []
        parallax = float_value(gaia_row["parallax"])
        gmag = float_value(gaia_row["phot_g_mean_mag"])
        if not in_narrow:
            if not (77 <= float(gaia_row["l"]) <= 83):
                reasons.append("outside l=[77,83]")
            if not (-1.5 <= float(gaia_row["b"]) <= 4):
                reasons.append("outside b=[-1.5,4]")
            if parallax is None:
                reasons.append("missing Gaia parallax")
            elif parallax < 0.35:
                reasons.append("parallax < 0.35 mas")
            elif parallax > 1.10:
                reasons.append("parallax > 1.10 mas")
            if gmag is None:
                reasons.append("missing Gaia G")
            elif gmag >= 19:
                reasons.append("G >= 19")
        method = "SIMBAD Gaia DR3 cross-identifier" if wright_simbad.get(name) else "Gaia coordinate match <=2 arcsec"
        wright_rows.append(
            {
                "wright_recno": int(index + 1),
                "object_name": name,
                "source_id": source_id,
                "match_method": method,
                "match_separation_arcsec": separation,
                "in_wp1_gaia_narrow": in_narrow,
                "narrow_exclusion_reason": "; ".join(reasons) if reasons else None,
                "gaia_ra": float(gaia_row["ra"]),
                "gaia_dec": float(gaia_row["dec"]),
                "gaia_l": float(gaia_row["l"]),
                "gaia_b": float(gaia_row["b"]),
                "gaia_g": gmag,
                "gaia_parallax": parallax,
                "gaia_parallax_error": float_value(gaia_row["parallax_error"]),
                "gaia_ruwe": float_value(gaia_row["ruwe"]),
            }
        )
        logt = float_value(row["logT"])
        anchors.append(
            base_anchor(
                uid=f"wright2015:{index + 1}", catalog="Wright et al. 2015",
                version="VizieR J/MNRAS/449/741", record_id=str(index + 1),
                name=name, source_id=source_id, ra=wright_coords[index].ra.deg,
                dec=wright_coords[index].dec.deg, spectral_type=text_value(row["SpType"]),
                teff=10**logt if logt is not None else None, av=float_value(row["AV"]),
                source_table="J/MNRAS/449/741/census", bibcode="2015MNRAS.449..741W",
                match_method=method, match_separation=separation,
                notes="Primary WP1 Wright-census gate sample",
            )
        )

    wright_frame = pd.DataFrame(wright_rows)
    wright_path = PROCESSED / "wp1_wright15_gaia_crossmatch.parquet"
    wright_frame.to_parquet(wright_path, index=False)

    dr2_ids: list[str] = []
    for table in [*berlanas19_tables, berlanas20_tables[3]]:
        if "GaiaDR2" in table.colnames:
            dr2_ids.extend(filter(None, (int_string(value) for value in table["GaiaDR2"])))
    dr2_mapping, dr2_table = dr2_to_dr3(dr2_ids)
    dr2_path = RAW / "berlanas_dr2_to_dr3_neighbourhood.ecsv"
    dr2_table.write(dr2_path, format="ascii.ecsv", overwrite=True)
    raw_records.append(
        {
            "catalog": "Gaia DR3 dr2_neighbourhood",
            "table": "Berlanas DR2 identifiers",
            "rows": len(dr2_table),
            "file": str(dr2_path.relative_to(ROOT)),
            "sha256": sha256(dr2_path),
        }
    )

    for table_index, table in enumerate(berlanas19_tables[:2], 1):
        coords = SkyCoord(table["_RA"], table["_DE"], unit="deg")
        fallback_ids, fallback_sep = nearest_local_ids(
            coords.ra.deg.tolist(), coords.dec.deg.tolist(), local
        )
        for index, row in enumerate(table):
            dr2 = int_string(row["GaiaDR2"])
            source_id = dr2_mapping.get(dr2) or fallback_ids[index]
            method = "Gaia DR2 neighbourhood" if dr2_mapping.get(dr2) else (
                "local Gaia coordinate match <=2 arcsec" if source_id else None
            )
            anchors.append(
                base_anchor(
                    uid=f"berlanas2019:a{table_index}:{index + 1}",
                    catalog="Berlanas et al. 2019",
                    version="VizieR J/MNRAS/484/1838", record_id=str(index + 1),
                    name=text_value(row["Name"]), source_id=source_id,
                    ra=float(coords[index].ra.deg), dec=float(coords[index].dec.deg),
                    spectral_type=None, group=text_value(row["Group"]) if "Group" in table.colnames else None,
                    source_table=f"J/MNRAS/484/1838/tablea{table_index}",
                    bibcode="2019MNRAS.484.1838B", match_method=method,
                    match_separation=fallback_sep[index] if method and method.startswith("local") else None,
                    notes=f"Published Gaia DR2 source_id={dr2}" if dr2 else None,
                )
            )

    b20_frames = [table.to_pandas() for table in berlanas20_tables[:4]]
    b20 = b20_frames[0]
    for table_number, frame in enumerate(b20_frames[1:], 2):
        b20 = b20.merge(
            frame, on="ID", how="left", suffixes=("", f"_tablea{table_number}")
        )
    b20_coords = SkyCoord(b20["RAJ2000"], b20["DEJ2000"], unit=(u.hourangle, u.deg))
    b20_fallback, b20_sep = nearest_local_ids(
        b20_coords.ra.deg.tolist(), b20_coords.dec.deg.tolist(), local
    )
    for index, row in b20.iterrows():
        dr2 = int_string(row.get("GaiaDR2"))
        source_id = dr2_mapping.get(dr2) or b20_fallback[index]
        method = "Gaia DR2 neighbourhood" if dr2_mapping.get(dr2) else (
            "local Gaia coordinate match <=2 arcsec" if source_id else None
        )
        anchors.append(
            base_anchor(
                uid=f"berlanas2020:{int(row['ID'])}", catalog="Berlanas et al. 2020",
                version="VizieR J/A+A/642/A168", record_id=str(int(row["ID"])),
                name=text_value(row.get("Name")), source_id=source_id,
                ra=float(b20_coords[index].ra.deg), dec=float(b20_coords[index].dec.deg),
                spectral_type=text_value(row.get("SpType")) or text_value(row.get("SpTypeP")),
                teff=float_value(row.get("Teff")), teff_error=float_value(row.get("e_Teff")),
                logg=float_value(row.get("logg")), logg_error=float_value(row.get("e_logg")),
                av=float_value(row.get("Av")), group=text_value(row.get("group")),
                source_table="J/A+A/642/A168 tables A1-A4 merged on ID",
                bibcode="2020A&A...642A.168B", match_method=method,
                match_separation=b20_sep[index] if method and method.startswith("local") else None,
                notes=f"Published Gaia DR2 source_id={dr2}" if dr2 else None,
            )
        )

    gosc_names = [str(value).strip() for value in gosc["Simbad"]]
    gosc_resolvable = [name for name in gosc_names if name]
    gosc_simbad = simbad_gaia_ids(gosc_resolvable, "gosc_field_simbad_ids.ecsv")
    gosc_ra = pd.to_numeric(gosc["RA ( o )"], errors="coerce").tolist()
    gosc_dec = pd.to_numeric(gosc["dec ( o )"], errors="coerce").tolist()
    gosc_fallback, gosc_sep = nearest_local_ids(gosc_ra, gosc_dec, local)
    for index, (_, row) in enumerate(gosc.iterrows()):
        name = str(row["Simbad"]).strip() or str(row["unnamed_0"]).strip()
        source_id = gosc_simbad.get(name) or gosc_fallback[index]
        method = "SIMBAD Gaia DR3 cross-identifier" if gosc_simbad.get(name) else (
            "local Gaia coordinate match <=2 arcsec" if source_id else None
        )
        spectral_type = "".join(
            str(row.get(column, "")).strip() for column in ["STv3", "LCv3", "qualv3"]
        ) or None
        anchors.append(
            base_anchor(
                uid=f"gosc:{row['GOS ID']}", catalog="Galactic O-Star Catalog",
                version="GOSC v3 / GOSSS-DR2.2 current web catalogue",
                record_id=str(row["GOS ID"]), name=name, source_id=source_id,
                ra=float(gosc_ra[index]), dec=float(gosc_dec[index]),
                spectral_type=spectral_type, group=text_value(row.get("associations")),
                source_table="GOSC main catalogue web query", bibcode="2013MSAIS..33...55M",
                match_method=method,
                match_separation=gosc_sep[index] if method and method.startswith("local") else None,
            )
        )

    for _, row in wr.iterrows():
        alias = str(row["Alias1"])
        match = re.search(r"DR3\s+(\d+)", alias)
        source_id = match.group(1) if match else None
        coord = SkyCoord(
            str(row["Right Ascension J2000"]), str(row["Declination J2000"]),
            unit=(u.hourangle, u.deg), frame="icrs",
        )
        current_association = text_value(row.get("Association"))
        wright_member = str(row["WR#"]) in {"144", "145", "146"}
        notes = (
            "Included in Wright+15 Cyg OB2 census; current WR catalogue association field="
            f"{current_association or 'blank'}"
            if wright_member else None
        )
        anchors.append(
            base_anchor(
                uid=f"wrcat:{row['ID']}", catalog="Galactic Wolf Rayet Catalogue",
                version="v1.33 (Gaia DR3), Aug 2025", record_id=str(row["WR#"]),
                name=f"WR {row['WR#']}", source_id=source_id, ra=float(coord.ra.deg),
                dec=float(coord.dec.deg), spectral_type=text_value(row["Spectral Type"]),
                group=current_association, source_table="WRcat v1.33 web catalogue",
                bibcode="2015MNRAS.447.2322R", match_method="catalogue-provided Gaia DR3 alias" if source_id else None,
                notes=notes,
            )
        )

    anchor_frame = pd.DataFrame(anchors)
    anchor_frame["source_id"] = anchor_frame["source_id"].astype("string")
    anchor_path = PROCESSED / "wp1_spectroscopic_anchors.parquet"
    anchor_frame.to_parquet(anchor_path, index=False)
    anchor_ecsv = PROCESSED / "wp1_spectroscopic_anchors.ecsv"
    ecsv_frame = anchor_frame.astype(object).where(pd.notna(anchor_frame), None)
    Table.from_pandas(ecsv_frame).write(anchor_ecsv, format="ascii.ecsv", overwrite=True)

    wright_matched = int(wright_frame["source_id"].notna().sum())
    wright_in_narrow = int(wright_frame["in_wp1_gaia_narrow"].sum())
    wright_gate = wright_matched / len(wright_frame)
    report_path = ROOT / "wp1_wright15_crossmatch.md"
    excluded = wright_frame[~wright_frame["in_wp1_gaia_narrow"]]
    lines = [
        "# Wright et al. (2015) to Gaia DR3 cross-match", "",
        f"Frozen: {created}", "",
        f"- Wright census rows: {len(wright_frame)}",
        f"- Rows assigned a Gaia DR3 source_id: {wright_matched} ({wright_gate:.1%})",
        f"- Rows present in `wp1_gaia_narrow`: {wright_in_narrow} ({wright_in_narrow/len(wright_frame):.1%})",
        f"- WP1 gate (>=90% matched to Gaia): **{'PASS' if wright_gate >= 0.9 else 'FAIL'}**", "",
        "The Gaia-match gate and narrow-selection inclusion are deliberately separate. "
        "The following matched stars are absent from the narrow sample:", "",
        "| Object | Gaia DR3 source_id | separation (arcsec) | reason absent from narrow |",
        "|---|---:|---:|---|",
    ]
    for _, row in excluded.iterrows():
        lines.append(
            f"| {row['object_name']} | {row['source_id']} | "
            f"{row['match_separation_arcsec']:.3f} | {row['narrow_exclusion_reason']} |"
        )
    lines.extend(
        [
            "", "The machine-readable audit is `data/processed/wp1_wright15_gaia_crossmatch.parquet`. "
            "All catalog rows, including bright/problematic Gaia sources, remain in the anchor table.", "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")

    execution = {
        "created_utc": created,
        "script": "scripts/wp1_build_anchors.py",
        "selection": {"l_deg": [77.0, 83.0], "b_deg": [-1.5, 4.0]},
        "catalog_inputs": raw_records,
        "wright_gate": {
            "census_rows": len(wright_frame),
            "matched_to_gaia_dr3": wright_matched,
            "match_fraction": wright_gate,
            "present_in_wp1_gaia_narrow": wright_in_narrow,
            "pass_at_least_90_percent": bool(wright_gate >= 0.9),
        },
        "anchor_counts": {
            "rows": len(anchor_frame),
            "rows_with_gaia_source_id": int(anchor_frame["source_id"].notna().sum()),
            "unique_gaia_source_ids": int(anchor_frame["source_id"].dropna().nunique()),
            "by_catalog": anchor_frame.groupby("source_catalog").size().to_dict(),
            "rows_with_spectral_type": int(anchor_frame["spectral_type"].notna().sum()),
            "rows_with_teff": int(anchor_frame["teff_K"].notna().sum()),
            "rows_with_logg": int(anchor_frame["logg_cgs"].notna().sum()),
        },
        "outputs": {
            str(anchor_path.relative_to(ROOT)): sha256(anchor_path),
            str(anchor_ecsv.relative_to(ROOT)): sha256(anchor_ecsv),
            str(wright_path.relative_to(ROOT)): sha256(wright_path),
            str(report_path.relative_to(ROOT)): sha256(report_path),
        },
    }
    write_json(PROVENANCE / "wp1_spectroscopic_anchors_execution.json", execution)
    print(json.dumps(execution["wright_gate"], indent=2))
    print(json.dumps(execution["anchor_counts"], indent=2))

    # Enforce the downstream counting invariant after every full catalogue rebuild.
    from wp1_refreeze_anchor_uniqueness import refreeze

    refreeze(created_override=created)


if __name__ == "__main__":
    main()
