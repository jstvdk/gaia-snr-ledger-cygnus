# WP1 provenance log

## Project and environment

- Project root: `/Users/vdk/science/gaia_snr_history_cygnus`
- WP: WP1 — Data acquisition
- Data release: Gaia DR3
- TAP endpoint: `https://gea.esac.esa.int/tap-server/tap`
- Query interface: ADQL over IVOA TAP, submitted with `pyvo`
- Conda environment: `cygob2-gaia`
- Environment verification date: 2026-07-20
- Verified packages: `pyvo 1.9.1`, `astropy 8.0.1`, `pandas 3.0.3`
- Raw query files: `queries/gaia_narrow.adql`, `queries/gaia_wide.adql`
- Submission script: `scripts/gaia_download.py`

## Acquisition design

- Narrow candidate region: Galactic `l=[77,83] deg`, `b=[-1.5,4] deg`.
- Wide candidate region: Galactic `l=[72,88] deg`, `b=[-5,8] deg`.
- Because the Gaia archive query uses ICRS geometry, both ADQL files use safe ICRS cones around RA `308.3 deg`, Dec `41.2 deg`; exact Galactic-coordinate filtering will be applied locally with Astropy.
- Parallax preselection: `0.35 <= parallax <= 1.10 mas`.
- Photometric preselection: `phot_g_mean_mag < 19`.
- Gaia table: `gaiadr3.gaia_source`.
- The raw parallax is retained. The Lindegren et al. zero-point correction will be applied downstream, with raw parallax, zero point, and corrected parallax stored separately.
- No `SELECT *` queries are permitted.

## 2026-07-20 execution record

- The user confirmed that the pilot query in `notebooks/gaia_pilot_query.ipynb` and `scripts/gaia_pilot_query.py` worked.
- The default Python interpreter lacked `pyvo`; the `cygob2-gaia` environment was verified instead.
- The environment verification succeeded with Python 3.11.15, `pyvo 1.9.1`, `astropy 8.0.1`, and `pandas 3.0.3`.
- Full-query execution is performed by `scripts/gaia_download.py`.
- Each completed TAP job writes FITS and Parquet outputs under `data/raw/gaia/` and a JSON execution record with row count, columns, job URL, and SHA-256 checksums under `provenance/`.

## Planned next acquisitions

1. Execute and validate the narrow Gaia query.
2. Execute the wide Gaia query for runaway/traceback support.
3. Apply the exact local Galactic-coordinate cuts and record pre-/post-filter row counts.
4. Add the Gaia–2MASS precomputed cross-match.
5. Build the spectroscopic-anchor table.
6. Freeze pulsar, SNR, 26Al, and extinction-map reference inputs.


## Incomplete execution record — 2026-07-20/21

- Query attempted: `queries/gaia_narrow.adql`.
- Submission method: `conda run -n cygob2-gaia --no-capture-output python scripts/gaia_download.py narrow --timeout 7200`.
- Archive status observed repeatedly: `EXECUTING`.
- The local process was interrupted before the TAP job reached `COMPLETED` and before `fetch_result()` returned.
- Result: no Gaia FITS/Parquet output was written; no `wp1_gaia_narrow_execution.json` completion record exists.
- Interpretation: WP1 narrow acquisition remains pending. Do not treat the attempted job as a downloaded dataset.


## Downloader failure and remediation — 2026-07-21

- The narrow query was submitted successfully and reported `EXECUTING` repeatedly.
- After approximately 27 polling cycles, the Gaia TAP connection closed during a status request with `RemoteDisconnected: Remote end closed connection without response`.
- This was a transient network/polling failure, not an ADQL syntax error. The old downloader had no retry logic and had not persisted the TAP job URL before completion.
- No Gaia data files were downloaded from that attempt.
- `scripts/gaia_download.py` was replaced with a resumable implementation that saves `provenance/wp1_gaia_<name>_job.json` immediately, preserves the server-side job, retries transient polling/result-download failures with exponential backoff, and supports `--new-job` for an explicitly requested fresh submission.
- Validation performed: Python byte-compilation and `--help` completed successfully in conda environment `cygob2-gaia`.
- No new Gaia query was submitted during remediation.

## Gaia retrieval sizing and tiling (2026-07-21T08:31:21.287322+00:00)

The initial synchronous `COUNT(*)` over the full 5 deg cone timed out at the Gaia TAP service. A cheaper, reproducible enclosing-rectangle count was then run over five contiguous RA strips with the same parallax and G-magnitude filters. These counts are upper bounds because the final tile queries also retain the original 5 deg cone predicate.

| Tile | RA interval (deg) | Dec interval (deg) | Enclosing-strip count |
|---|---:|---:|---:|
| 01 | [301.6, 304.3) | [36.2, 46.2] | 362,902 |
| 02 | [304.3, 307.0) | [36.2, 46.2] | 241,580 |
| 03 | [307.0, 309.7) | [36.2, 46.2] | 160,951 |
| 04 | [309.7, 312.4) | [36.2, 46.2] | 227,953 |
| 05 | [312.4, 315.0) | [36.2, 46.2] | 276,944 |
| **Total upper bound** | 301.6–315.0 | 36.2–46.2 | **1,270,330** |

Generated query files:
- `queries/gaia_narrow_tile01.adql`
- `queries/gaia_narrow_tile02.adql`
- `queries/gaia_narrow_tile03.adql`
- `queries/gaia_narrow_tile04.adql`
- `queries/gaia_narrow_tile05.adql`

Sizing command used:
```bash
conda run -n cygob2-gaia --no-capture-output python - <<'PY'
import pyvo
tap = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
for i, (ra0, ra1) in enumerate([(301.6,304.3),(304.3,307.0),(307.0,309.7),(309.7,312.4),(312.4,315.0)], 1):
    q = f"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE ra >= {ra0} AND ra < {ra1} AND dec BETWEEN 36.2 AND 46.2 AND parallax BETWEEN 0.35 AND 1.10 AND phot_g_mean_mag < 19"
    print(i, tap.search(q).to_table()[0][0])
PY
```

## Tiled downloader handoff (2026-07-21T08:34:20.639382+00:00)

- `scripts/gaia_download.py` now accepts any validated query name whose file exists as `queries/gaia_<name>.adql`; this supports `narrow_tile01` through `narrow_tile05` while preserving resumable TAP job state and retry behavior.
- A follow-up exact `COUNT(*)` attempt on tile 01 with the cone predicate also timed out at the archive. No data retrieval was started in this sizing step.
- The five tile intervals are contiguous and half-open in RA, so a source cannot be duplicated at a tile boundary. The cone predicate is retained in every tile.

- 2026-07-21T08:36:25.451968+00:00: User-requested local stop of the original `narrow` downloader (PIDs 62016/62018/62057). No output files existed; saved TAP job state was preserved for possible resume.

## WP1 closure execution — 2026-07-22

This section supersedes the earlier "Planned next acquisitions" list. The aggregate verdict is recorded in `provenance/wp1_manifest.json` and `wp1_completion_report.md`; every file used by that verdict was re-hashed at validation time.

### Gaia narrow

- Canonical artifacts: `data/processed/wp1_gaia_narrow.parquet` and `.fits`.
- Exact frozen selection: Galactic `l=[77,83] deg`, `b=[-1.5,4] deg`, `0.35 <= parallax <= 1.10 mas`, `G < 19`.
- Validated rows: **245,843**, all unique by Gaia DR3 `source_id`.
- The five raw TAP tiles and exact local Galactic-coordinate validation are recorded in `provenance/wp1_gaia_narrow_validation.json`.

### Official Gaia DR3–2MASS PSC join

- Query family: `queries/gaia_2mass_join.adql` and the six non-overlapping `queries/gaia_2mass_tile01.adql` through `tile06.adql` files.
- Tables: `gaiadr3.tmass_psc_xsc_best_neighbour`, `gaiadr3.tmass_psc_xsc_join`, and `gaiadr1.tmass_original_valid`.
- The monolithic job `eba92310-85b5-11f1-982e-bc97e148b76b-O` ended in server phase `ERROR` with PostgreSQL `current transaction is aborted`; the same query was therefore partitioned into six one-degree Galactic-longitude tiles.
- Tile row counts: 58,873; 36,516; 30,400; 35,895; 30,252; 36,719. Union: **228,655** rows, zero duplicate `source_id`s.
- One query-side boundary row (`source_id=2070682352787315200`) is not in the frozen narrow catalogue: Gaia's stored `l` passes the upper bound, while the ICRS position transforms with Astropy to `l=83.00000571998599 deg` (0.0206 arcsec outside). The canonical product keys to the frozen narrow IDs and explicitly excludes this row.
- Canonical `wp1_2mass_join` has one row for every narrow Gaia source: **245,843** rows, **228,654** official PSC matches (93.01%), **203,646** complete J/H/Ks measurements with finite uncertainties (82.84%), and **170,757** `ph_qual=AAA` rows (69.46%). All 17,189 unmatched sources are retained with null photometry and explicit flags.
- The archive returned variable-length object strings for 2MASS identifiers. `scripts/gaia_download.py` now normalizes only object-string columns before FITS/Parquet export; the normalized export was tested in both formats before the completed jobs were resumed.

### Spectroscopic anchors and Wright+15 gate

- Builder: `scripts/wp1_build_anchors.py`; outputs: `data/processed/wp1_spectroscopic_anchors.parquet`, `.ecsv`, and `data/processed/wp1_wright15_gaia_crossmatch.parquet`.
- Frozen long-form rows: **540** (Wright+15 167; Berlanas+19 229; Berlanas+20 78; current GOSC field entries 56; current Galactic WR field entries 10).
- Rows with Gaia DR3 IDs: 536; unique Gaia DR3 IDs: 252; rows with spectral type: 311; with Teff: 227; with log g: 63.
- Wright+15 gate: **167/167 = 100%** matched to Gaia DR3, exceeding the required 90%. Of these, 160 are in the frozen narrow catalogue. The seven matched narrow exclusions and their individual parallax/missing-parallax reasons are listed in `wp1_wright15_crossmatch.md`; bright or astrometrically problematic anchors remain in the anchor table.

### Extinction references

- Baseline: Vergely et al. (2022) 25 pc correlation-length differential-extinction cube, frozen locally as `data/raw/extinction/explore_cube_density_values_025pc_v2.fits` (81 x 601 x 601; 6 x 6 x 0.8 kpc coverage).
- Independent Cygnus check: all five Dharmawardena et al. (2022) products, covering `l=[73,87] deg`, `b=[-4,6] deg`, and distance 1.28–2.20 kpc.
- A 36-point Bayestar web-API coverage request was attempted and frozen, but the service returned HTTP 500. It is not used as a completion dependency because two local 3D-map products are present and checksummed.
- Interpretation and units are documented in `wp1_extinction_refs.md`; these maps are consistency checks, not substitutes for WP3 per-star fits.

### SN and feedback markers

- ATNF catalogue v2.8.1 frozen: 4,149 coordinate-parsed pulsars, 80 in the wide box. PSR J2032+4127 is present with its period, period derivative, characteristic age, proper motion, parallax, and association distance fields.
- Green 2024 catalogue: 310 total SNRs, 9 in the wide box, including G078.2+02.1 (gamma Cygni). The Chandra SNRcat snapshot contributes six wide-box entries.
- The Manitoba SNRcat endpoint timed out; the access failure is preserved rather than silently substituted.
- Martin et al. Cygnus 26Al primary PDFs are frozen. The extracted 1809-keV values and errors are recorded in `wp1_sn_markers.md`.

### Gaia wide

- Independent exact count query: `queries/gaia_wide_exact_count.adql` returned **3,133,326**.
- Four non-overlapping exact Galactic-longitude tiles returned 1,001,200; 691,996; 563,495; and 876,635 rows.
- The two largest completed archive results were streamed directly to local VOTables after long-lived `pyvo` in-memory transfers stalled. `scripts/wp1_import_gaia_votable.py` lower-cased the archive's FITS-style column names, exported FITS/Parquet, and recorded transport-file hashes and the original TAP job URLs.
- Canonical `data/processed/wp1_gaia_wide.parquet`: **3,133,326 rows, 3,133,326 distinct int64 source IDs, zero duplicates**, exactly equal to the count query. The validator uses Arrow-native distinct counting so Gaia IDs are never coerced through imprecise float64 values.

### Commands and final gate

```bash
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_build_anchors.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_freeze_extinction_refs.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_freeze_sn_markers.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_combine_2mass_tiles.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_finalize_2mass.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_finalize_wide.py
conda run -n cygob2-gaia --no-capture-output python scripts/wp1_validate.py
```

Final aggregate result: `wp3_inputs_ready=true`; `wp1_complete=true`. The final validator re-opened the execution records, re-hashed all named outputs, and passed every row-count, uniqueness, Wright-match, anchor-content, extinction-reference, marker, and wide-count check.

## Wright+15 identity correction and anchor re-freeze — 2026-07-22

- Trigger: a post-freeze identity audit found that the 167 machine-readable Wright+15 census rows contain two pairs resolving to the same Gaia DR3 sources. The earlier `167/167` row-level statement was therefore not a valid physical-star denominator and is superseded by this section.
- Collision 1: Wright recno 46 and 57 are both named `Schulte 41` and both resolved through SIMBAD to Gaia DR3 `2067782936320586240`. Recno 46 is **951.685 arcsec** from that Gaia position; it is retained as `duplicate_bad_coordinate` with `source_id=null`, `resolved_gaia_source_id` populated, and recno 57 identified as the canonical record.
- Collision 2: Wright recno 69 (`2MASS J20323949+4052475`) and recno 81 (`[CPR2002] A31`) resolve to Gaia DR3 `2067768195993175936`. Recno 81 is retained as a flagged `duplicate_alias` record with recno 69 canonical.
- Denominator reconciliation: Wright+15's abstract and conclusions state 169 primary OB stars, while Section 2.1 states that the spectroscopic compilation produced 167 stars and the paper later describes a 167-star IMF sample. The frozen VizieR table contains 167 rows. Gaia identity resolution reduces those rows to **165 unique physical stars**. The operational WP1 denominator is therefore explicitly 165 unique stars represented by the ingested machine-readable table; coverage relative to the abstract statement is 165/169 (97.6%), and the paper/table count discrepancy is not hidden.
- Corrected crossmatch: 167 evidence rows, 165 countable non-null `source_id`s, 165 unique non-null `source_id`s, two flagged duplicate records, 165/165 valid Gaia matches, and 158/165 present in the frozen narrow catalogue.
- Corrected canonical anchor: **252 rows and 252 unique non-null Gaia DR3 source IDs**, with zero duplicates. This is the only anchor artifact downstream census code may count.
- Preserved evidence: all **540** literature records are retained in `data/processed/wp1_spectroscopic_anchor_records.parquet` and `.ecsv`, including four records without Gaia IDs and both flagged Wright duplicates. Cross-catalogue aliases are encoded in JSON provenance columns of the canonical row.
- Re-freeze implementation: `scripts/wp1_refreeze_anchor_uniqueness.py`. `scripts/wp1_build_anchors.py` invokes the same invariant after every future full rebuild.
- Updated records: `provenance/wp1_spectroscopic_anchors_execution.json`, `wp1_wright15_crossmatch.md`, `provenance/wp1_manifest.json`, and `wp1_completion_report.md`. The manifest update was produced by rerunning `scripts/wp1_validate.py`, which re-hashed the corrected products.
