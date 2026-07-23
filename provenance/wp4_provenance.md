# WP4 provenance — subgroup ages & per-star masses

**Generated:** 2026-07-23 (UTC) · DR3 · single-population distance 1.6245 ± 0.045 kpc (μ = 11.054).
**Manifest:** `provenance/wp4_manifest.json` (SHA-256 per input / script / execution-log / output).
**Status:** WP4_COMPLETE_GATE_PASSED.

This file is the reproducibility appendix for WP4. It records the pipeline, every
threshold with its class (per `CUTS_AND_THRESHOLDS.md §1`), and the sensitivity
register (`§9`, `§10`).

---

## 1. Pipeline (scripts, in execution order)

| # | Script | Role | Execution log |
|---|---|---|---|
| 1 | `scripts/wp4_common.py` | constants, loaders, extinction-error propagation, forward synthetic-CMD model, per-star weighted likelihood, posterior fitter | (imported) |
| 2 | `scripts/wp4_fit_ages.py` | full branch-grid age fit → `wp4_age_posteriors.parquet`, `wp4_posterior_curves.npz` | `wp4_fit_ages_execution.json` |
| 3 | `scripts/wp4_anchors_hrd.py` | spectroscopic anchors on HRD, masses, gate | `wp4_anchors_hrd_execution.json` |
| 4 | `scripts/wp4_clump.py` | high-A_V clump classification + age-robustness test | `wp4_clump_execution.json` |
| 5 | `scripts/wp4_masses.py` | per-star masses (6 branches) + spectroscopic override | `wp4_masses_execution.json` |
| 6 | `scripts/wp4_figures.py` | CMD / HRD / posterior / money-plot figures | (deterministic) |
| 7 | `scripts/wp4_closure_audit.py` | schema/count gates + WP6 PM-candidate export | `wp4_closure_audit.json` |
| 8 | `scripts/wp4_report.py` | tables, `wp4_masses.cat`, manifest, regenerated report blocks | (deterministic) |

Reproduce end-to-end from `scripts/`:
```
conda run -n cygob2-gaia --no-capture-output env PYTHONPATH=scripts \
  python scripts/wp4_fit_ages.py
# Then run anchors_hrd, clump, masses, figures, closure_audit and report with
# the same Conda/PYTHONPATH prefix, in that order.
```
All steps are deterministic (no RNG; the synthetic population is a deterministic
IMF-weighted mass grid × mass-ratio grid). Total runtime ≈ 35 s.

## 2. Frozen inputs (SHA-256, see manifest for full digests)

| Input | SHA-256 (first 16) |
|---|---|
| `data/processed/wp3_extinction.parquet` | `f83a612215e9e0e0` |
| `data/processed/wp3_isochrones_parsec.parquet` | `9e0444045acf3e93` |
| `data/processed/wp3_isochrones_mist.parquet` | `3b1437d791e61a5d` |
| `tables/wp2_subgroup_labels.parquet` | `44a5ee7b9190703d` |
| `data/processed/wp2_anchor_assignments.parquet` | `22d7fc622dc0fc99` |
| `data/processed/wp1_spectroscopic_anchors.parquet` | `9bb30e91d98568b3` |

The two isochrone grids are **byte-identical to the WP3 files** (digests match
`provenance/wp3_manifest.json` `shared_with_WP4`) — the plan's WP3/WP4 isochrone
consistency requirement is satisfied by construction: WP4 reads, never refetches.

**source_id join note:** `wp1_spectroscopic_anchors.parquet` stores `source_id`
as a string; it is coerced to int64 before joining to the int64 member keys.
190 anchors match into the 1,392 members; 156 carry a spectroscopic Teff.

## 3. Threshold register (class-labelled)

| Cut / knob | Value | Class | Justification | Range tested | Effect on age |
|---|---|---|---|---|---|
| Isochrone family | PARSEC, MIST | **E** | genuine model ambiguity; carried in parallel, never averaged | both | dominant systematic, MIST younger by 0.4–0.9 Myr |
| R_V | 3.0 / 3.1 / 3.5 | **E** | dust-law ambiguity (WP3) | all three | ≤ 1 grid step (~0.4 Myr) |
| f_bin (binary fraction) | 0.3 / 0.4 / 0.5 | **D/E** | O-star multiplicity ~70%; scanned | all three | negligible on upper-MS |
| Distance modulus | 11.054 ± 0.060 | **A** | derived from WP2 posterior; propagated at μ ± σ_μ | ±0.060 mag | ±0.4–0.9 Myr (coherent across A/B/C) |
| IMF slope (density weight) | 2.3 | **B** | Kroupa/Salpeter; sets along-isochrone density only | — | not a driver (morphology sets age) |
| q distribution | U[0.1, 1] | **B** | standard binary mass-ratio prior | — | negligible |
| σ_int (model floor) | 0.03 mag | **D** | isochrone/calibration imperfection | — | negligible |
| Mag floor | 0.02 mag | **D** | photometric systematic floor | — | negligible |
| ums window | M_G0 ≤ +1.5 | **C/D** | bright, well-populated regime | — | defines indicator |
| pms window | M_G0 ≥ +2.5 | **C/D** | faint contracting regime; data-limited | — | measurable only for A, B |
| N_min per indicator | 15 stars | **B** | population-level measurability convention; registered in CUTS §5.4 | 10/20/30 | full retained MAP span is 2.25–5.67 at 10/15 and 2.25–4.50 at 20/30 |
| anchor χ tolerance | 2.5 | **B** | ~2.5σ HRD match incl. binary brightening | — | gate pass/fail |
| logTe_max_trust | log(52000) | **B** | above this = WR/stripped, off normal isochrones | — | 2 anchors excluded |
| clump box | see below | **C** | reproduces the WP3 flag | — | age shift 0.00 Myr when removed |

**clump box:** (BP−RP)₀ ∈ [−0.4, 0.5], M_G0 ∈ [−1, 2.5], ≥ 0.5 mag above the
single-star MS.

`mass_baseline = mass_PARSEC_rv3.1` (PARSEC, R_V = 3.1) is a reporting
convenience only and has no downstream authority; WP5 must carry all six mass
branches separately.

The HRD age lookup is algorithmically symmetric between families: for each
subgroup, PARSEC and MIST each use their own independently fitted upper-MS MAP.
All three baseline MIST fits select the same 3.571-Myr native grid point, so the
constant `age_used_MIST` is a fit result, not a fixed-age input. Consequently,
χ_PARSEC versus χ_MIST compares each family's best-age morphology, not two
families evaluated at one common numerical age.

## 4. Sensitivity register (downstream-relevant, per `§9`)

- **Age vs isochrone family:** MIST 0.4–0.9 Myr younger than PARSEC at fixed data. **Dominant systematic.** Carried as a branch into WP5/WP7.
- **Age vs distance:** ±0.4–0.9 Myr for Δμ = ±0.060 mag (A most sensitive in MIST, ±0.93). Coherent across subgroups → does not manufacture a relative age spread.
- **Age vs R_V:** ≤ one grid step (~0.4 Myr).
- **Age vs f_bin:** negligible on the upper MS (< 0.05 Myr); binaries are captured by the explicit binary component, not by clipping.
- **Age vs clump removal:** 0.00 Myr in every subgroup/family — the high-A_V clump does not drive the age.
- **Upper-MS vs PMS turn-on:** B agrees at baseline, but its retained R_V=3.5
  MIST PMS branch reaches 5.67 Myr with N=19; A's upper-MS/PMS offset is
  documented; C's PMS is starved and excluded.
- **Age-row exclusions:** 104/132 retained; 28 excluded. Nineteen fail only
  N≥15, and nine fail N≥15 and are grid-railed. No measurable row is railed.

## 5. Gate

The gate is satisfied with explicit limitations — see `wp4_ages.md §10`.
A/B/C remain consistent with a coeval upper MS (retained MAPs 3.16–4.50 Myr),
while the honest two-indicator span is 2.25–5.67 Myr. The association-wide HRD
gate passes at 131/150 (87%), but subgroup B has only five anchors and is marked
limited. The indicator gate is satisfied by **documented disagreement**, not a
blanket agreement claim. The star-formation-duration branch (0/1/2 Myr) remains
mandatory for WP7.

<!-- BEGIN GENERATED:REPORT_REGENERATION -->
Final deterministic report regeneration: `2026-07-23T14:05:25.562369+00:00` UTC, after all WP4 Parquet, NPZ, figure, table, audit and manifest products.

Diff against the pre-closure reports: headline values changed from an unqualified 3.5–4.5 Myr statement to a 3.16–4.50 Myr upper-MS envelope plus the honest 2.25–5.67 Myr two-indicator envelope; the indicator gate changed from PASS to documented disagreement; subgroup-B anchor coverage (N=5), the MIST equal-age fit result, 28 age-row exclusions, the 55 common massless rows, nine astrometry-less exceptions and seven WP6 PM candidates were added. Membership counts, Berlanas recall, anchor-HRD aggregate 131/150, baseline mass summaries and the coeval upper-MS verdict did not move.
<!-- END GENERATED:REPORT_REGENERATION -->

## 6. Known limitations (carried forward honestly)

- The upper-MS **statistical** credible intervals (±0.05 Myr) are internal
  precision only; the true uncertainty is the **branch envelope** (~±1 Myr),
  systematic-dominated. Not tuned narrower.
- The **PMS turn-on is data-starved** at this distance+extinction (faint edge
  near M_G0 ≈ +3); measurable only for A and B. A DR4 fainter sample would
  strengthen the second age indicator.
- The forward model uses **independent (colour, M_G0) Gaussian errors**; the
  small colour↔magnitude covariance from shared A_V error is not modelled
  (sub-dominant to the branch envelope).
- Photometric masses for the **hottest (M_G0 < −4) stars** are colour-degenerate;
  those are superseded by the spectroscopic-HRD masses where an anchor exists.
- All six mass columns share the same **55 null source_ids**: 46 optical-only
  Gaia members have no 2MASS JHK and received no frozen WP3 A_V/de-reddened CMD
  coordinate; nine manual Berlanas quality exceptions are absent from the
  frozen narrow Gaia query and have no absolute G. Their missing masses are
  input-driven and not branch-specific.
- Anchor HRD coverage is lopsided: A/B/C/unassigned = **60/5/43/48**. B's five
  anchors are all consistent, but N=5 is not an independently strong subgroup
  validation; 48 informative anchors are unassigned because the subgroup fit
  deliberately excluded all quality exceptions.
- Nine P>0.5 members have no Gaia astrometry or photometry because they are
  explicit Berlanas+19 spectroscopic manual quality exceptions absent from the
  frozen narrow query. Their `membership_probability=1` is literature-based;
  `membership_probability_astrometric` remains null. Counts and Berlanas recall
  are unchanged.
- All 52 P>0.5 members with RUWE>1.4 are explicit
  `anchor_quality_exempt=True` rows (maximum RUWE 24.03), consistent with the
  documented binary-bias exception.
- Seven >5σ proper-motion candidates are frozen in
  `provenance/wp4_pm_outliers.csv` for WP6. None is cut from membership.
