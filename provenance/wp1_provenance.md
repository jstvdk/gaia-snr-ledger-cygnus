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
