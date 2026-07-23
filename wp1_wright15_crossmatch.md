# Wright et al. (2015) to Gaia DR3 cross-match

Re-frozen: 2026-07-22T13:54:33.752887+00:00

## Denominator reconciliation

- Wright+15 abstract and conclusions: **169 primary OB stars**.
- Wright+15 Section 2.1 and the paper's later IMF sample: **167 stars**.
- Frozen VizieR `J/MNRAS/449/741/census`: **167 rows**.
- Unique physical stars represented after Gaia identity resolution: **165**.
- WP1 Gaia-match denominator: **165 unique physical stars represented by the frozen machine-readable table**; it is not silently equated with the abstract's 169.
- Representation relative to the abstract count: 165/169 (97.6%). The four-star difference is the combination of the paper's internal 169-versus-167 count discrepancy and two duplicate physical identities in the machine-readable rows.

## Gaia gate

- Valid unique Gaia matches: 165/165 (100.0%).
- Valid unique stars present in `wp1_gaia_narrow`: 158/165 (95.8%).
- WP1 gate (>=90% of the declared denominator matched): **PASS**.
- Non-null duplicate `source_id` values in this crossmatch: **0**.

## Duplicate-record audit

Rejected records are retained as flagged rows with `source_id=null`, while `resolved_gaia_source_id` and `duplicate_of_source_id` preserve the identity evidence.

| Rejected recno | Object | Resolved Gaia DR3 source_id | Canonical recno | separation (arcsec) | status |
|---:|---|---:|---:|---:|---|
| 46 | Schulte 41 | 2067782936320586240 | 57 | 951.685 | duplicate_bad_coordinate: SIMBAD name resolution assigned the canonical Gaia source, but the published coordinate is 951.685 arcsec away |
| 81 | [CPR2002] A31 | 2067768195993175936 | 69 | 0.382 | duplicate_alias: separate Wright catalogue alias resolves to the same Gaia source |

## Valid unique stars absent from the narrow sample

| Object | Gaia DR3 source_id | separation (arcsec) | reason absent from narrow |
|---|---:|---:|---|
| WR 146 | 2067867358199590400 | 0.117 | missing Gaia parallax |
| [MT91] 516 | 2067778542573330304 | 0.315 | parallax < 0.35 mas |
| [CPR2002] A20 | 2067755379810392320 | 0.361 | parallax > 1.10 mas |
| [MT91] 140 | 2067828467273010176 | 0.754 | parallax > 1.10 mas |
| [MT91] 605 | 2067878357612771200 | 0.636 | parallax < 0.35 mas |
| [MT91] 129 | 2067848219824399616 | 0.297 | parallax > 1.10 mas |
| [MT91] 459 | 2067785070923660416 | 0.095 | parallax < 0.35 mas |

The canonical `wp1_spectroscopic_anchors` table contains one row per Gaia source. All literature rows, aliases, rejected Wright records, and source-specific provenance remain in `data/processed/wp1_spectroscopic_anchor_records.parquet`.
