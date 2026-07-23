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
| 7 | `scripts/wp4_report.py` | tables, `wp4_masses.cat`, manifest | (deterministic) |

Reproduce end-to-end from `scripts/`:
```
python3 wp4_fit_ages.py && python3 wp4_anchors_hrd.py && python3 wp4_clump.py \
  && python3 wp4_masses.py && python3 wp4_figures.py && python3 wp4_report.py
```
All steps are deterministic (no RNG; the synthetic population is a deterministic
IMF-weighted mass grid × mass-ratio grid). Total runtime ≈ 35 s.

## 2. Frozen inputs (SHA-256, see manifest for full digests)

| Input | SHA-256 (first 16) |
|---|---|
| `data/processed/wp3_extinction.parquet` | `278dce68647f8319` |
| `data/processed/wp3_isochrones_parsec.parquet` | `9e0444045acf3e93` |
| `data/processed/wp3_isochrones_mist.parquet` | `3b1437d791e61a5d` |
| `tables/wp2_subgroup_labels.parquet` | `4cccab14878c1f87` |
| `data/processed/wp2_anchor_assignments.parquet` | `26290a7e543b0a0a` |
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
| N_min per indicator | 15 stars | **D** | measurability floor | — | flags C-pms as starved |
| anchor χ tolerance | 2.5 | **B** | ~2.5σ HRD match incl. binary brightening | — | gate pass/fail |
| logTe_max_trust | log(52000) | **B** | above this = WR/stripped, off normal isochrones | — | 2 anchors excluded |
| clump box | see below | **C** | reproduces the WP3 flag | — | age shift 0.00 Myr when removed |

**clump box:** (BP−RP)₀ ∈ [−0.4, 0.5], M_G0 ∈ [−1, 2.5], ≥ 0.5 mag above the
single-star MS.

`mass_baseline = mass_PARSEC_rv3.1` (PARSEC, R_V = 3.1) is a reporting
convenience only and has no downstream authority; WP5 must carry all six mass
branches separately.

## 4. Sensitivity register (downstream-relevant, per `§9`)

- **Age vs isochrone family:** MIST 0.4–0.9 Myr younger than PARSEC at fixed data. **Dominant systematic.** Carried as a branch into WP5/WP7.
- **Age vs distance:** ±0.4–0.9 Myr for Δμ = ±0.060 mag (A most sensitive in MIST, ±0.93). Coherent across subgroups → does not manufacture a relative age spread.
- **Age vs R_V:** ≤ one grid step (~0.4 Myr).
- **Age vs f_bin:** negligible on the upper MS (< 0.05 Myr); binaries are captured by the explicit binary component, not by clipping.
- **Age vs clump removal:** 0.00 Myr in every subgroup/family — the high-A_V clump does not drive the age.
- **Upper-MS vs PMS turn-on:** agree for B; A's upper-MS is ~1.5 Myr older than its PMS (extended-SF signature, documented); C's PMS is starved (n=3, not measurable).

## 5. Gate

Three criteria, all passed — see `wp4_ages.md §10`. Headline: A/B/C **consistent
with coeval at ~3.5–4.5 Myr**; the coeval result is stated explicitly rather than
forced into a spread. Anchor HRD 87% consistent (median χ < 1); all outliers are
over-luminous binaries. The **star-formation-duration branch (0/1/2 Myr) is carried
into WP7** regardless of the coeval finding.

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
