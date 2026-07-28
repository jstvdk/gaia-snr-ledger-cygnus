# WP5 parent-range fix (repair_v2)

*Executed 2026-07-27. Applies the correction specified in
`reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md`. All `repair_v1`
artifacts are untouched; this is a WP5-only re-run consuming the `repair_v1`
WP3/WP4 upstream unchanged.*

## Verdict

**The diagnosed bug is fixed and the fix behaved exactly as predicted. WP5
remains BLOCKED on a different, independent failure. WP6 is still not
authorized.**

## What was changed

| file | change |
|---|---|
| `scripts/wp5_common.py` | `MASS_GRID` extended from 0.5–8.0 to **0.5–18.0 M☉** (45 points); new `PARENT_MASS_HI = 18.0`. `CALIBRATION_HI` unchanged at 8.0 — the *observed* window is identical. |
| `scripts/wp5_fit_imf.py` | forward-response clip changed from `CALIBRATION_HI` to `MASS_GRID.max()`; provenance now records `parent_mass_range_Msun` separately from `nominal_window_Msun`. |
| `scripts/wp5_injections_repair.py` | new `--output-version`; `scripts/wp5_fit_imf.py` new `--wp5-version`. Lets a WP5-only re-run write versioned products without overwriting an earlier repair version. |
| `CUTS_AND_THRESHOLDS.md` | new §7.2 (window edge vs parent range) and a §7.1 amendment recording the unreachable 95% edge. |

Grid spacing: 0.25 M☉ through 0.5–8.0, 0.5 through 8.5–12.0, 1.0 through
13–18. Non-uniform trapezoid integration verified accurate to **0.17%** against
a 20,001-point reference at every IMF slope.

**Ceiling justification (measured, not chosen).** The >8 M☉ contribution to the
top-bin rate reaches 95.8% of its converged value at a 12 M☉ ceiling, 99.9% at
16, and 100.0% at 18, identically for α ∈ {2.0, 2.3, 2.6}.

## Pre-flight checks

| check | result |
|---|---|
| Isochrone MS coverage at fitted ages, both families | Mini range extends to 47–96 M☉ after the `label`/`phase` ≤ 1 cut — 18 M☉ is nowhere near extrapolation |
| Donor pool covers the bright end | pool spans G = 8.23–19.0 with 284 donors brighter than G = 12.5; brightest synthetic in the extended grid is G ≈ 11.5 — no railing |
| `repair_v1` products preserved | confirmed unmodified |

## Execution

324,000 injections (45 masses × 400 × 18 branches; was 223,200), all 18
branches, 64 mass-posterior draws each.
Provenance: `provenance/wp5_injections_repair_execution_repair_v2.json`,
`provenance/wp5_imf_fit_execution_repair_v2.json`.

## Result: the fix did what it claimed

Top-bin (6.35–8.0 M☉) Pearson residual, baseline PARSEC R_V = 3.1, α = 2.3:

| subgroup | repair_v1 | repair_v2 | change |
|---|---:|---:|---:|
| CygOB2-A | +2.63 | **+1.47** | −1.16 |
| CygOB2-B | +1.08 | **+0.25** | −0.82 |
| CygOB2-C | +3.26 | **+1.89** | −1.37 |

The top bin improved in all three subgroups, in the predicted direction and by a
comparable magnitude. **CygOB2-C — the subgroup that blocked repair_v1 — now
passes on every criterion** (χ² p 0.040 → 0.285, max |residual| 3.26 → 1.89).
The diagnosis in the correction report is confirmed.

## Result: the gate still fails, on a different bin

**Correction to the estimate in the diagnosis report.** That report projected
baseline PASS and 43/54 branches. The delivered numbers are **baseline FAIL and
33/54** (up from 31/54). The estimate was optimistic because it added only the
top-bin term while holding every other bin's response fixed; the real re-run
re-estimated the response in all six bins and lowered k by 2–5%, which raised
the mid-window residuals. The top-bin prediction itself was accurate.

Baseline gate, repair_v2:

| subgroup | χ² p | trend p | max abs residual | worst bin | pass |
|---|---:|---:|---:|---|:--:|
| CygOB2-A | 0.384 | 0.872 | 1.47 | 5 (6.35–8.0) | yes |
| CygOB2-B | 0.036 | 0.208 | **3.02** | **2 (3.17–4.0)** | **no** |
| CygOB2-C | 0.285 | 0.787 | 1.89 | 5 (6.35–8.0) | yes |

CygOB2-B's bin-2 residual is **not new** — it was already the worst bin in B at
+2.75 under repair_v1. It is a pre-existing, independent anomaly in the middle
of the window that the truncation fix never addressed; the ~5% reduction in k
pushed it from just under the 3.0 threshold to just over.

Across all 54 branches the 21 failures now split into two clean populations:

| worst bin | count | interpretation |
|---|---:|---|
| bin 2 (3.17–4.0 M☉) | 9 | the new blocker, overwhelmingly CygOB2-B |
| bin 5 (6.35–8.0 M☉) | 8 | residual top-bin excess, confined to R_V = 3.5 and α = 2.6 corners |
| bins 0/1/4 | 4 | scattered |

## Characterization of the CygOB2-B bin-2 anomaly

It is **not** a lower-edge or completeness-ramp artifact. Raising the
calibration lower edge does not fix it, it only relocates the failure:

| lower edge | branches passing | baseline max abs residual A / B / C |
|---:|---:|---|
| 2.0 M☉ | 33/54 | 1.47 / **3.02** / 1.89 |
| 2.5 M☉ | 30/54 | 1.87 / **4.44** / 1.85 |
| 3.0 M☉ | 28/54 | 2.52 / **4.71** / 1.78 |
| 3.5 M☉ | 33/54 | 2.94 / 1.63 / 1.96* |

(*fails on χ²/trend rather than the residual bound.) This rules out the remedy
prescribed in plan WP5 step 2 and `CUTS_AND_THRESHOLDS.md` §7.1 — raising the
lower edge makes CygOB2-A and CygOB2-B worse before it makes them better.

**It is strongly and monotonically R_V-dependent**, in both isochrone families.
CygOB2-B bin-2 residual at α = 2.3:

| R_V | PARSEC | MIST |
|---:|---:|---:|
| 3.0 | +1.72 | +0.55 |
| 3.1 | +3.02 | +1.53 |
| 3.5 | +4.20 | +4.26 |

That is an **extinction-law signature**, not an IMF or membership signature. An
IMF or contamination problem would not order itself by R_V; a systematic error
in the adopted extinction law would.

**The likely mechanism.** CygOB2-B has the narrowest extinction distribution of
the three subgroups and the highest median:

| subgroup | N | median A_V | half of central-68% spread |
|---|---:|---:|---:|
| CygOB2-A | 476 | 5.59 | 0.917 |
| CygOB2-B | 426 | **6.56** | **0.196** |
| CygOB2-C | 429 | 6.51 | 0.575 |

B's members are extinguished almost uniformly (σ ≈ 0.20 mag, 4.7× narrower than
A). A coherent A_V distribution means any systematic error in the extinction law
displaces B's entire mass function **rigidly**, piling stars into one bin,
instead of being smeared across bins as it is in A. That is exactly the observed
signature: one sharp bin excess, scaling with R_V, in the subgroup with the
tightest A_V.

## Recommended next diagnostic

The question is whether B's very narrow A_V distribution is **physical** or is
an artifact of the repair_v1 WP3 extinction prior. repair_v1 introduced an
eight-anchor spatial prior with measured widths 0.452/0.453/0.475 mag. If
CygOB2-B's members sit close to a single anchor, that prior could be collapsing
their A_V posterior toward one value, manufacturing the coherence that makes B
R_V-sensitive.

Concrete checks, in order of cost:

1. Map CygOB2-B members against the eight anchor positions; compute the
   effective prior weight per star. If B is dominated by one anchor, the narrow
   A_V is a prior artifact.
2. Compare B's A_V distribution against the Vergely+22 and Dharmawardena+22
   cubes at B's footprint and distance — the WP3 Step 7 machinery already does
   this and can be run per subgroup.
3. Refit B's extinction with the spatial prior widened or removed, and re-run
   the WP5 baseline for B only.

Only if all three come back clean does the subgroup-model class of hypothesis
(label contamination, intrinsic age spread) become the live one — and then for
CygOB2-B, not CygOB2-C.

## Artifacts

| kind | path |
|---|---|
| injection response | `data/processed/wp5_injection_response_repair_v2.parquet` (324,000 rows) |
| completeness curves | `data/processed/wp5_completeness_curves_repair_v2.parquet` |
| normalization | `data/processed/wp5_imf_normalization_repair_v2.parquet` (54 fits) |
| mass-function bins | `data/processed/wp5_mass_function_bins_repair_v2.parquet` |
| association mass | `data/processed/wp5_association_mass_repair_v2.parquet` |
| posterior draws | `data/processed/wp5_imf_posterior_draws_repair_v2.npz` |
| branch table | `tables/wp5_imf_norm_repair_v2.csv` |
| baseline residuals | `tables/wp5_baseline_residuals_repair_v2.csv` |
| baseline completeness | `tables/wp5_completeness_baseline_repair_v2.csv` |
| provenance | `provenance/wp5_injections_repair_execution_repair_v2.json`, `provenance/wp5_imf_fit_execution_repair_v2.json` |

Baseline association mass 30,453 M☉ (multiplicity-adjusted), within a factor 2
of the 16,500 M☉ literature scale — unchanged in character from repair_v1's
31,293 M☉.

The headline deliverables `tables/wp5_imf_norm.csv` / `.md` and
`wp5_imf_norm.md` were **not** regenerated. `scripts/wp5_report.py` is hardwired
to the unversioned frozen products, and promoting a still-blocked fit into the
headline tables would be wrong. They remain stale pre-repair artifacts; see
open issue #2 in `PROJECT_TRACE.md`.
