# Project trace — WP0 to WP5

*Master index for "The Supernova History of Cygnus OB2 from Gaia DR3" (Paper 1).
Created 2026-07-27. Every artifact produced so far is linked from here, with its
gate verdict and the numbers that downstream work depends on. This is a
navigation and audit document — it duplicates no content, it points at it.*

**Governing documents:** [paper1_execution_plan.md](paper1_execution_plan.md)
(what to do) · [method_explained.md](method_explained.md) (why) ·
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) (how every number was chosen) ·
[GLOSSARY.md](GLOSSARY.md) (nomenclature) · [AUDIT.txt](AUDIT.txt) (generated
file inventory with checksums, produced by [audit.py](audit.py))

---

## 1. Status board

| WP | Objective | Gate verdict | Downstream authorized |
|---|---|---|---|
| WP0 | Requirements + dedup | **PASS** — no direct duplicate; one Gaia DR3 partial overlap to cite | yes |
| WP1 | Data acquisition | **PASS** | yes |
| WP2 | Membership + substructure | **PASS** — 8/8 criteria, after one rejected run | yes |
| WP3 | Per-star extinction | **PASS**, then repaired (repair_v1) | yes |
| WP4 | Ages + masses | **SATISFIED WITH CONSTRAINTS**, then repaired (repair_v1) | yes, with constraints |
| WP5 | IMF normalization + completeness | **PASS — accepted 2026-07-28 (`repair_v6`).** Gate G3 clears under the **strict** per-branch reading with `A_or_C_regressions: []`; all three subgroups pass the baseline residual gate under both trend statistics; branch grid **40/54**; association mass 29,122 M☉ | **yes** |
| WP6 | Census closure + runaways | **AUTHORIZED, not started.** `downstream_wp6_authorized: true`. Its one known design defect is already fixed: issue #3 closed, closure estimator binding in [CUTS §16](CUTS_AND_THRESHOLDS.md), bright-mass response measured on the accepted version | — |
| WP7–WP10 | Ledger, cross-checks, verdict, manuscript | not started | — |

**Current status — updated 2026-07-28 (evening). The CygOB2-B problem is
solved.** It had two independent causes, each found, separately pre-registered,
and confirmed:

1. **The injection truth model used a single age** (the upper-MS MAP), and at
   B's fitted age the PMS/Henyey isochrone fold sat on the 3.17–4.0 M☉ bin
   boundary, manufacturing a localized bump. Fixed in `repair_v4` by truth-side
   joint age–k marginalization.
2. **B's extinction was biased high by the anchor spatial prior.** The prior's
   width grew with distance from the calibrating anchors, correctly declaring
   them uninformative, while its *mean* still took those same anchors at full
   weight — the two moments contradicted each other. B's nearest anchors sit
   0.377° away against 0.089° for A and 0.139° for C. Fixed in `repair_v5` by
   replacing the neighbour median with a simple-kriging mean using the
   *already-fitted* variogram; no new parameter, and repair_v3 behaviour is
   bit-preserved under the old mode.

Fix 2 was pre-registered before it was run
([wp3_kriging_prior_prereg.json](provenance/wp3_kriging_prior_prereg.json)) and
**all four predictions were confirmed**
([wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json)):
B's A_V fell 0.359 mag; B's age rose 2.818 → **3.548 Myr** with a 68% interval
[3.389, 3.991] that reaches the age at which WP5's mass function independently
passed; A and C moved by **0.000 Myr**; and B's residual tilt flattened, the
replacement trend statistic going |T| 2.40 → **1.31**.

B's baseline residuals across the three versions:

| version | bin 0 | bin 1 | bin 2 | bin 3 | bin 4 | bin 5 | T | gate |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| repair_v3 | +1.45 | +2.41 | **+3.62** | −0.36 | −1.30 | −0.92 | −3.22 | fail |
| repair_v4 | +1.20 | +2.08 | **+2.51** | −0.24 | −0.79 | −0.54 | −2.40 | fail |
| **repair_v5** | +0.87 | +2.27 | **+0.36** | −0.24 | −0.03 | +0.17 | **−1.31** | **pass** |

**`repair_v5` is the first version whose baseline passes for all three
subgroups under *both* the incumbent and the better-calibrated replacement
trend statistic.** Branch grid **26 → 38 of 54** (37 under the replacement).

**WP5 is nevertheless still not accepted.** Under the strict per-branch reading
of gate G3 adopted on 2026-07-27, three A/C non-baseline cells regress against
repair_v4: two sit at trend p = 0.040 and 0.048, inside the indeterminate band
around the 0.05 threshold that **issue #11** identified and that criterion R3
has not yet resolved; the third, CygOB2-C MIST R_V = 3.5 α = 2.0, is a genuine
max-residual failure at 3.25. `downstream_wp6_authorized` stays `false`. What
blocks WP5 is now a threshold-methodology question, not an astrophysical one.

---

## 2. WP0 — Requirements extraction and dedup

| | |
|---|---|
| **Verdict** | No direct duplicate. One relevant Gaia DR3 partial overlap identified and must be cited. |
| **Deliverables** | [wp0_requirements_table.md](wp0_requirements_table.md) · [wp0_prior_art_summary.md](wp0_prior_art_summary.md) · [wp0_dedup_report.md](wp0_dedup_report.md) |

The requirements table defines the statement WP9 must eventually evaluate. Dedup
covered ADS 2023–2026 citers of Martin+10 / Berlanas+20 / Härer+25 plus ICRC and
TeVPA 2025 materials.

---

## 3. WP1 — Data acquisition

| | |
|---|---|
| **Verdict** | **PASS** · WP3 input gate READY |
| **Report** | [wp1_completion_report.md](wp1_completion_report.md) |
| **Provenance** | [wp1_provenance.md](provenance/wp1_provenance.md) · [wp1_manifest.json](provenance/wp1_manifest.json) |

**Headline numbers**

| quantity | value |
|---|---:|
| Gaia DR3 narrow-box unique sources | 245,843 |
| Gaia DR3 wide-box unique sources | 3,133,326 |
| Official 2MASS PSC matches | 228,654 |
| Complete J/H/Ks with finite errors | 203,646 |
| Wright+15 → Gaia match rate | **165 / 165 (100%)** |
| Spectroscopic anchors (unique Gaia rows) | 252 (540 literature-evidence rows) |
| ATNF pulsars / Green SNRs in wide box | 80 / 9 |

**Scripts** — [gaia_download.py](scripts/gaia_download.py) ·
[wp1_import_gaia_votable.py](scripts/wp1_import_gaia_votable.py) ·
[wp1_combine_2mass_tiles.py](scripts/wp1_combine_2mass_tiles.py) ·
[wp1_finalize_2mass.py](scripts/wp1_finalize_2mass.py) ·
[wp1_finalize_wide.py](scripts/wp1_finalize_wide.py) ·
[wp1_build_anchors.py](scripts/wp1_build_anchors.py) ·
[wp1_refreeze_anchor_uniqueness.py](scripts/wp1_refreeze_anchor_uniqueness.py) ·
[wp1_freeze_extinction_refs.py](scripts/wp1_freeze_extinction_refs.py) ·
[wp1_freeze_sn_markers.py](scripts/wp1_freeze_sn_markers.py) ·
[wp1_validate.py](scripts/wp1_validate.py)

**Data** — `data/processed/wp1_gaia_narrow.parquet` ·
`wp1_gaia_wide.parquet` · `wp1_2mass_join.parquet` ·
`wp1_spectroscopic_anchors.parquet` · `wp1_wright15_gaia_crossmatch.parquet` ·
`wp1_atnf_pulsars_wide.parquet` · `wp1_green_snrs_wide.parquet` ·
`wp1_chandra_snrcat_wide.parquet`

**Supporting docs** — [wp1_wright15_crossmatch.md](wp1_wright15_crossmatch.md) ·
[wp1_sn_markers.md](wp1_sn_markers.md) ·
[wp1_extinction_refs.md](wp1_extinction_refs.md) ·
[table1_wp1_inventory.md](tables/table1_wp1_inventory.md)

**Notebooks** — [gaia_pilot_query.ipynb](notebooks/gaia_pilot_query.ipynb) ·
[wp1_step1_build_narrow_catalogue.ipynb](notebooks/wp1_step1_build_narrow_catalogue.ipynb)

**Per-step execution records** — `provenance/wp1_*_execution.json` (one per tile
and per step), `provenance/wp1_*_job.json` (archive job IDs),
`provenance/wp1_gaia_narrow_validation.json`,
`provenance/wp1_gaia_wide_validation.json`

---

## 4. WP2 — Membership and substructure

| | |
|---|---|
| **Verdict** | **PASS** — all 8 gate criteria |
| **Deliverable** | [wp2_subgroups.md](wp2_subgroups.md) |
| **Gate table** | [table2_wp2_gate.md](tables/table2_wp2_gate.md) |
| **Provenance** | [wp2_provenance.md](provenance/wp2_provenance.md) · [wp2_membership_manifest.json](provenance/wp2_membership_manifest.json) |

**Gate results**

| criterion | threshold | achieved |
|---|---|---:|
| Berlanas+19 recall at P > 0.5 | ≥ 0.80 | 0.825 |
| Control / target yield ratio | ≤ 0.10 | 0.0366 |
| Member count | 10²–10⁴ | 1,392 |
| Central-90% l span | < 4.8° | 1.07 |
| Central-90% b span | < 4.4° | 0.89 |
| Convex-hull area | < 16.5 deg² | 2.08 |
| Largest seed / analysis sample | ≤ 0.10 | 0.0080 |
| Published-structure comparison | documented | yes |

**Sample structure** — `data/processed/wp2_members.parquet` holds 2,112 rows
(everything above P > 0.05). Of these 1,392 have P > 0.5; 1,331 carry an A/B/C
subgroup label; 61 are anchor-quality-exempt spectroscopic members; 781 are
`unassigned` (median P = 0.25) and are carried for provenance only.

| subgroup | N labelled | median P | mean P |
|---|---:|---:|---:|
| CygOB2-A | 476 | 0.97 | 0.89 |
| CygOB2-B | 426 | 0.91 | 0.84 |
| CygOB2-C | 429 | 0.91 | 0.85 |
| *unassigned* | 781 | 0.25 | 0.30 |

**Method notes** — clustering in (l, b, μα*, μδ); **parallax deliberately
excluded** (one population at 1.62 kpc, intrinsic depth 45 pc, exhausted at DR3
precision). Gaussian mixtures k = 2…8 over 50 deterministic seeds; acceptance by
**seed stability (ARI), not BIC**. Only k = 3 is stable. HDBSCAN scanned for
completeness, not used for labels. The Berlanas two-distance hypothesis is
**not confirmed** — see [wp2_subgroups.md](wp2_subgroups.md) §6.

**Scripts** — [wp2_membership_pipeline.py](scripts/wp2_membership_pipeline.py) ·
[wp2_finalize_membership.py](scripts/wp2_finalize_membership.py) ·
[wp2_derive_subgroups.py](scripts/wp2_derive_subgroups.py) ·
[wp2_distance_population_test.py](scripts/wp2_distance_population_test.py) ·
[wp2_finalize_audits.py](scripts/wp2_finalize_audits.py)

**Data** — `wp2_members.parquet` · `wp2_control_members.parquet` ·
`wp2_anchor_assignments.parquet` ·
`wp2_distance_component_probabilities.parquet` ·
[wp2_subgroup_labels.parquet](tables/wp2_subgroup_labels.parquet)

**Figures** — [wp2_membership_sky_pm.png](figures/wp2/wp2_membership_sky_pm.png) ·
[wp2_subgroups_sky.png](figures/wp2/wp2_subgroups_sky.png) ·
[wp2_subgroups_vpd.png](figures/wp2/wp2_subgroups_vpd.png) ·
[wp2_subgroups_extinction.png](figures/wp2/wp2_subgroups_extinction.png) ·
[wp2_distance_population_test.png](figures/wp2/wp2_distance_population_test.png) ·
[kdistance_min15.png](figures/wp2/kdistance_min15.png) ·
[target_kdistance_min15.png](figures/wp2/target_kdistance_min15.png)

**Paper figures** — [fig1_membership_literature.pdf](figures/paper/fig1_membership_literature.pdf) ·
[fig2_control_fields.pdf](figures/paper/fig2_control_fields.pdf), built by
[make_paper_figures_wp1_wp2.py](scripts/make_paper_figures_wp1_wp2.py)

**Notebook** — [wp2_membership_and_substructure.ipynb](notebooks/wp2_membership_and_substructure.ipynb)

**Rejected / superseded** — the 2026-07-22 baseline run is preserved, not
deleted: `data/processed/wp2_members_failed_20260722.parquet`,
[wp2_membership_manifest_failed_20260722.json](provenance/wp2_membership_manifest_failed_20260722.json),
[wp2_membership_attempt_deconvolved_failure.json](provenance/wp2_membership_attempt_deconvolved_failure.json),
[wp2_distance_population_attempt_failure.json](provenance/wp2_distance_population_attempt_failure.json).
Hyperparameter scans: `wp2_dbscan_scan.csv`, `wp2_hdbscan_scan.csv`,
`wp2_hdbscan_leaf_scan.csv`, `wp2_hdbscan_subgroup_scan.csv`,
`wp2_gmm_seed_stability.csv`, `wp2_mixture_prior_control_scan.csv`,
`wp2_literature_footprint_density_scan.csv`, all under `provenance/`.

**Literature recovery** — [table3_literature_recovery.md](tables/table3_literature_recovery.md) ·
[wp2_berlanas_recovery_audit.csv](provenance/wp2_berlanas_recovery_audit.csv) ·
[wp2_paiz_crossmatch.csv](provenance/wp2_paiz_crossmatch.csv)

---

## 5. WP3 — Per-star extinction and de-reddened CMD

| | |
|---|---|
| **Verdict** | **PASS** · WP4 input gate READY |
| **Report** | [wp3_completion_report.md](wp3_completion_report.md) |
| **Provenance** | [wp3_provenance.md](provenance/wp3_provenance.md) · [wp3_manifest.json](provenance/wp3_manifest.json) |

**Headline numbers** — distance posterior 1.6245 ± 0.045 kpc (μ = 11.05).
A_V median 5.85, span ~2–10.6 mag, core 4–8 populated. Spatial coherence
ρ = 0.66 at 18σ (Vergely cube ρ = 0.29). Two marginally negative extinctions,
none beyond 2σ.

**Method** — anchors get A_V from fixed-Teff intrinsic colours (avoids the
Teff/A_V degeneracy for hot stars); all other members get a χ²-posterior over
PARSEC+MIST templates with μ fixed. Isochrones are **shared with WP4** by
construction. R_V branches 3.0 / 3.1 / 3.5 carried in parallel from the start.

**Scripts** — [wp3_fetch_parsec.py](scripts/wp3_fetch_parsec.py) ·
[wp3_fetch_mist.py](scripts/wp3_fetch_mist.py) ·
[wp3_extinction_law.py](scripts/wp3_extinction_law.py) ·
[wp3_assemble_photometry.py](scripts/wp3_assemble_photometry.py) ·
[wp3_anchor_extinction.py](scripts/wp3_anchor_extinction.py) ·
[wp3_broadband_extinction.py](scripts/wp3_broadband_extinction.py) ·
[wp3_build_catalogue.py](scripts/wp3_build_catalogue.py) ·
[wp3_cube_comparison.py](scripts/wp3_cube_comparison.py) ·
[wp3_figures_and_coherence.py](scripts/wp3_figures_and_coherence.py) ·
[wp3_common.py](scripts/wp3_common.py)

**Data** — `wp3_extinction.parquet` · `wp3_member_photometry.parquet` ·
`wp3_anchor_extinction.parquet` · `wp3_broadband_extinction.parquet` ·
`wp3_cube_comparison.parquet` · `wp3_isochrones_parsec.parquet` ·
`wp3_isochrones_mist.parquet` (see
[wp3_isochrones_README.md](provenance/wp3_isochrones_README.md))

**Tables** — [wp3_extinction_rv30.cat](tables/wp3_extinction_rv30.cat) ·
[wp3_extinction_rv31.cat](tables/wp3_extinction_rv31.cat) ·
[wp3_extinction_rv35.cat](tables/wp3_extinction_rv35.cat)

**Figures** — [wp3_dereddened_cmd.png](figures/wp3/wp3_dereddened_cmd.png) ·
[wp3_extinction_map.png](figures/wp3/wp3_extinction_map.png) ·
[wp3_av_distribution.png](figures/wp3/wp3_av_distribution.png) ·
[wp3_wp5_truncation.png](figures/wp3/wp3_wp5_truncation.png)

**Note for WP5** — [wp3_provenance.md](provenance/wp3_provenance.md) Step 8
anticipates the WP5 truncation bias: the G < 19 limit removes 2 M☉ stars above
A_V ≈ 6.5, 3 M☉ above ≈ 8.7, 5 M☉ above ≈ 10.1 (4 Myr PARSEC). Completeness must
therefore be a function of both mass and position.

**Extinction method as adopted (2026-07-28).** The working extinction is
`repair_v5`: the anchor spatial prior's **mean** is a simple-kriging estimate
over the 8 nearest anchors using the repair_v3 variogram, and its **width** is
the unconditional variogram sigma. See
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §15 and
[WP3_KRIGED_PRIOR_repair_v5.md](reports/WP3_KRIGED_PRIOR_repair_v5.md).
Median A_V by subgroup: A 5.501, B **6.031** (was 6.391), C 6.435.
Carried systematic: broadband photometry sits ~0.5 mag below spectroscopically
calibrated anchors at matched position, so the absolute extinction and mass
scales rest on the anchor calibration (issue #1d, obligation O2).

---

## 6. WP4 — Subgroup ages and per-star masses

| | |
|---|---|
| **Verdict** | **SATISFIED** — one association-wide PASS with a subgroup-B limitation and one documented indicator disagreement · WP5/WP7 input gate READY WITH CONSTRAINTS |
| **Report** | [wp4_completion_report.md](wp4_completion_report.md) · [wp4_ages.md](wp4_ages.md) |
| **Provenance** | [wp4_provenance.md](provenance/wp4_provenance.md) · [wp4_manifest.json](provenance/wp4_manifest.json) |

**Headline result** — the three kinematic subgroups remain consistent with a
coeval upper main sequence. Retained upper-MS branch MAPs span
**3.16–4.50 Myr**; the honest envelope across both retained indicators is
**2.25–5.67 Myr**. ≈3.5–4.5 Myr is a central upper-MS summary, *not* the full
envelope. The A/B/C age ordering is unresolved, so the star-formation-duration
branch (0/1/2 Myr) remains mandatory for WP7.

| audit quantity | value |
|---|---|
| Age branches measurable | 104 / 132 (28 excluded, 9 grid-railed) |
| Mass branches retained | 6 / 6 · `mass_baseline` = PARSEC R_V=3.1, reporting only |
| HRD anchors A/B/C/unassigned | 60 / 5 / 43 / 48 |
| WP6 proper-motion hand-off | 7 candidates, none removed from membership |

**Scripts** — [wp4_fit_ages.py](scripts/wp4_fit_ages.py) ·
[wp4_masses.py](scripts/wp4_masses.py) ·
[wp4_anchors_hrd.py](scripts/wp4_anchors_hrd.py) ·
[wp4_clump.py](scripts/wp4_clump.py) ·
[wp4_figures.py](scripts/wp4_figures.py) ·
[wp4_closure_audit.py](scripts/wp4_closure_audit.py) ·
[wp4_report.py](scripts/wp4_report.py) ·
[wp4_schema_repair.py](scripts/wp4_schema_repair.py) ·
[wp4_common.py](scripts/wp4_common.py)

**Data** — `wp4_age_posteriors.parquet` · `wp4_masses.parquet` ·
`wp4_anchor_hrd.parquet` · `wp4_clump.parquet` · `wp4_posterior_curves.npz`

**Tables** — [wp4_ages_table.md](tables/wp4_ages_table.md) ·
[wp4_ages_envelope.md](tables/wp4_ages_envelope.md) ·
[wp4_ages_summary.csv](tables/wp4_ages_summary.csv) ·
[wp4_masses.cat](tables/wp4_masses.cat)

**Figures** — [wp4_age_posteriors.png](figures/wp4/wp4_age_posteriors.png) ·
[wp4_age_summary.png](figures/wp4/wp4_age_summary.png) ·
[wp4_cmd_subgroups.png](figures/wp4/wp4_cmd_subgroups.png) ·
[wp4_hrd_anchors.png](figures/wp4/wp4_hrd_anchors.png)

**Audits** — [wp4_closure_audit.json](provenance/wp4_closure_audit.json) ·
[wp4_pm_outliers.csv](provenance/wp4_pm_outliers.csv) ·
[task_B_wp4_closure_1.md](tasks/task_B_wp4_closure_1.md)

**Ages as adopted (repair_v5, 2026-07-28).** Upper-MS MAP ages on the baseline
branch: **A 3.981, B 3.548, C 2.512 Myr**. CygOB2-B moved +0.730 Myr from
repair_v3's 2.818 when its extinction was corrected; A and C did not move at
all. The total spread is unchanged at 1.469 Myr, but the association now reads
as **two older subgroups (A, B) plus one younger (C)** rather than one older
plus two younger — a revised star-formation history that WP7 inherits
(obligation O1). B's PMS indicator remains unmeasurable (n = 4, grid-railed),
so issue #9's indicator disagreement is documented rather than resolved.

---

## 7. WP5 — IMF normalization and completeness · **BLOCKED**

| | |
|---|---|
| **Verdict** | **BLOCKED** — mass-function residual gate fails |
| **Report** | [wp5_completion_report.md](wp5_completion_report.md) · [wp5_imf_norm.md](wp5_imf_norm.md) (copy at [tables/wp5_imf_norm.md](tables/wp5_imf_norm.md)) |
| **Provenance** | [wp5_provenance.md](provenance/wp5_provenance.md) · [wp5_manifest.json](provenance/wp5_manifest.json) · [wp5_validation.json](provenance/wp5_validation.json) |

**Gate definition** ([wp5_fit_imf.py:297-299](scripts/wp5_fit_imf.py#L297-L299)) —
a branch passes only if `chi_p ≥ 0.01` **and** `trend_p ≥ 0.05` **and**
`max_abs_pearson_residual ≤ 3.0`, for all three subgroups.

**Current state after repair_v5 (adopted)** — baseline (PARSEC, R_V = 3.1, α = 2.3):

| subgroup | χ² p | trend p (incumbent) | trend p (replacement) | max abs residual | pass |
|---|---:|---:|---:|---:|:--:|
| CygOB2-A | 0.448 | 0.468 | 0.590 | 1.34 | yes |
| CygOB2-B | 0.293 | 0.111 | 0.186 | **2.27** (was 3.62) | **yes** |
| CygOB2-C | 0.350 | 0.957 | 0.858 | 1.56 | yes |

This is the **first version to pass the baseline under both trend statistics**.

Branch grid across versions: **0/54** (original) → **31/54** (repair_v1) →
**33/54** (repair_v2) → **26/54** (repair_v3) → **29/54** (repair_v4) → **38/54** (repair_v5). The
repair_v3 dip is real and explained in §8 — the corrected extinction prior
unmasked the CygOB2-B residual rather than causing it, and repair_v4 then
removed it at the source.

Association mass at baseline: **29,185 [28,116, 30,275] M☉** (repair_v5), within a
factor two of the 16,500 M☉ literature scale.

**WP5 is still not accepted.** Under the incumbent trend statistic the block is
issue #11 (Monte-Carlo instability) via the strict reading of gate G3; under
the pre-declared replacement statistic B's baseline fails outright on a
surviving mass-dependent tilt (issue #12, trend p = 0.017). The table above
therefore reads "pass" only against a trend test now known to be underpowered
against exactly this residual shape. Gate record:
[wp5_repair_v4_gate.json](provenance/wp5_repair_v4_gate.json); verdict report:
[wp5_completion_report_repair_v4.md](wp5_completion_report_repair_v4.md).

324,000 catalogue injections. Association mass 30,696 M☉ (within a factor 2 of
the 16,500 M☉ literature scale). No subgroup reaches an absolute 95%
completeness edge; bright-plateau completeness is ~0.79–0.80.

**Blocker 1 — top-bin parent truncation: FIXED in repair_v2.** The forward
response integrated the parent mass distribution only to 8 M☉ while the mass
estimator has 22–25% scatter, so the ~300 living stars above 8 M☉ that scatter
down into the top bin had no term in the model. Diagnosis:
[WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md).
Execution and outcome:
[WP5_PARENT_RANGE_FIX_repair_v2.md](reports/WP5_PARENT_RANGE_FIX_repair_v2.md).
Top-bin residuals fell in all three subgroups (A −1.16, B −0.82, C −1.37) and
CygOB2-C now passes on every criterion.

**Blocker 2 — CygOB2-B mid-window excess: OPEN.** Pre-existing (it was already
B's worst bin at +2.75 under repair_v1) and independent of the truncation; the
~5% drop in k from the fix pushed it over the bound. It is **not** a lower-edge
artifact — raising the calibration lower edge to 2.5 / 3.0 / 3.5 M☉ only
relocates the failure (33 → 30 → 28 → 33 branches passing). It **is** strongly
and monotonically R_V-dependent in both isochrone families (B bin-2 residual at
α = 2.3: PARSEC 1.72 / 3.02 / 4.20 and MIST 0.55 / 1.53 / 4.26 at
R_V = 3.0 / 3.1 / 3.5), which is an extinction-law signature rather than an IMF
or membership one. Likely mechanism: CygOB2-B has the narrowest A_V distribution
of the three subgroups (half-central-68% spread 0.196 mag vs 0.917 for A) and
the highest median (6.56), so an extinction-law error displaces B's mass
function rigidly instead of smearing it. Next diagnostic is the WP3 anchor
spatial prior — see the report's "Recommended next diagnostic".

**Scripts** — [wp5_common.py](scripts/wp5_common.py) ·
[wp5_injections.py](scripts/wp5_injections.py) ·
[wp5_fit_imf.py](scripts/wp5_fit_imf.py) ·
[wp5_report.py](scripts/wp5_report.py) ·
[wp5_finalize.py](scripts/wp5_finalize.py) ·
[wp5_make_notebook.py](scripts/wp5_make_notebook.py)

**Data** — `wp5_injection_response.parquet` ·
`wp5_completeness_curves.parquet` · `wp5_imf_normalization.parquet` ·
`wp5_mass_function_bins.parquet` · `wp5_imf_posterior_draws.npz` ·
`wp5_association_mass.parquet`

**Tables** — [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) **(stale — pre-repair,
0/54 passing; must be regenerated)** ·
[wp5_association_mass.csv](tables/wp5_association_mass.csv)

**Figures** — [wp5_completeness_curves.png](figures/wp5/wp5_completeness_curves.png) ·
[wp5_mass_function.png](figures/wp5/wp5_mass_function.png) ·
[wp5_association_mass.png](figures/wp5/wp5_association_mass.png)

**Notebook** — [wp5_imf_normalization_and_completeness.ipynb](notebooks/wp5_imf_normalization_and_completeness.ipynb)

---

## 8. Repair history

Repairs are versioned side-by-side with the frozen originals, never overwriting
them. `repair_v1` covers WP3, WP4 and WP5 as one coupled attempt.

| | |
|---|---|
| **Brief** | [wp3_extinction_repair_brief.md](tasks/wp3_extinction_repair_brief.md) |
| **Open handoff brief** | [wp5_cygob2b_massfunction_brief.md](tasks/wp5_cygob2b_massfunction_brief.md) — the current blocker (#1c), written for a fresh agent |
| **Report** | [WP3_WP5_REPAIR_REPORT_repair_v1.md](reports/WP3_WP5_REPAIR_REPORT_repair_v1.md) *(frozen; diagnosis superseded)* |
| **Correction** | [WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md) *(2026-07-27)* |
| **Gate record** | [wp3_repair_gate.json](provenance/wp3_repair_gate.json) — `accepted: false`, `downstream_wp6_authorized: false` |
| **Manifest** | [wp3_wp5_repair_manifest.json](provenance/wp3_wp5_repair_manifest.json) |
| **Notebook** | [wp3_extinction_repair_and_wp5_regate.ipynb](notebooks/wp3_extinction_repair_and_wp5_regate.ipynb) |

**What repair_v1 changed**

- **WP3** — 0.03-mag error floor in every band; eight-anchor spatial prior with
  measured widths 0.452 / 0.453 / 0.475 mag; full gridded extinction posterior;
  hidden template-branch width calibrated on the 149 spectroscopic anchors.
- **WP4** — ages refitted from the repaired extinction catalogue; masses are
  posterior samples over the full extinction distribution, all six bands, the
  binary branch and the subgroup age posterior; six-band model width 0.38 mag;
  log-uniform mass measure, **not** an IMF prior; spectroscopic-HRD overrides
  unchanged.
- **WP5** — injections re-run through the actual repaired WP3/WP4 estimators;
  baseline response resolution raised from 16 to 64 mass draws.

**Effect** — the branch grid went from **0 / 54** passing to **31 / 54**.
The baseline still fails on CygOB2-C.

### repair_v2 — WP5 parent-range fix (2026-07-27)

| | |
|---|---|
| **Diagnosis** | [WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md) |
| **Report** | [WP5_PARENT_RANGE_FIX_repair_v2.md](reports/WP5_PARENT_RANGE_FIX_repair_v2.md) |
| **Provenance** | [wp5_injections_repair_execution_repair_v2.json](provenance/wp5_injections_repair_execution_repair_v2.json) · [wp5_imf_fit_execution_repair_v2.json](provenance/wp5_imf_fit_execution_repair_v2.json) |
| **Scope** | WP5 only. Consumes the `repair_v1` WP3/WP4 upstream unchanged; every `repair_v1` artifact is preserved byte-identical. |

**What changed** — `MASS_GRID` extended 0.5–8.0 → **0.5–18.0 M☉** (45 points),
new `PARENT_MASS_HI = 18.0`; the forward-response clip in `wp5_fit_imf.py` moved
from `CALIBRATION_HI` to `MASS_GRID.max()`. `CALIBRATION_HI` stays 8.0 — the
*observed* window is unchanged. New `--output-version` / `--wp5-version` flags
keep WP5-only re-runs from overwriting earlier versions.
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) gained §7.2 (window edge vs
parent range) and a §7.1 amendment.

**Effect** — 324,000 injections (was 223,200); branch grid **31 → 33 / 54**;
top-bin residuals fell in all three subgroups; CygOB2-C cleared. Baseline still
blocked, now on CygOB2-B bin 2.

**Tables** — [wp5_imf_norm_repair_v2.csv](tables/wp5_imf_norm_repair_v2.csv) ·
[wp5_baseline_residuals_repair_v2.csv](tables/wp5_baseline_residuals_repair_v2.csv) ·
[wp5_completeness_baseline_repair_v2.csv](tables/wp5_completeness_baseline_repair_v2.csv)

### repair_v3 — WP3 variogram anchor prior (2026-07-27) · **negative result**

| | |
|---|---|
| **Report** | [WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md](reports/WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md) |
| **Chain runner** | [run_repair_v3_chain.sh](scripts/run_repair_v3_chain.sh) — WP3 → WP4 ages → WP4 masses → WP5 injections → WP5 fit |
| **Diagnostic** | [wp3_anchor_prior_diagnostic.py](scripts/wp3_anchor_prior_diagnostic.py) → [wp3_anchor_prior_diagnostic_execution.json](provenance/wp3_anchor_prior_diagnostic_execution.json) |
| **Lower-edge scan** | [wp5_lower_edge_scan.py](scripts/wp5_lower_edge_scan.py) → [wp5_lower_edge_scan_execution.json](provenance/wp5_lower_edge_scan_execution.json) |
| **Scope** | Full WP3→WP5 chain. repair_v1 and repair_v2 artifacts untouched. |

**Diagnosis (confirmed).** The repair_v1 prior gives each member the median A_V
of its 8 nearest spectroscopic anchors with one global width per R_V, calibrated
by leave-one-out at the anchor density (0.071°). CygOB2-B's members sit 0.377°
from their anchors, where the fitted variogram (nugget 0.000, sill 1.228 mag,
range 0.853°) says the required width is 1.052 mag — **understated 2.32×**. The
prior collapsed B's differential extinction **9.5×** (1.860 → 0.196) versus 2.0×
for A, and the independent Vergely/Dharmawardena cubes rank B's patchiness at or
above A's.

**Fix.** New `ANCHOR_PRIOR_MODE` in
[wp3_repair_common.py](scripts/wp3_repair_common.py): per-star prior width read
off the fitted variogram at each star's own 8th-anchor separation. Both
`WP_REPAIR_VERSION` and `WP3_ANCHOR_PRIOR_MODE` are environment variables
defaulting to repair_v1 behaviour, so the chain re-runs at any version with no
script edits and repair_v1 stays exactly reproducible. Extinction outcome was as
intended: B's spread restored 3.1× (0.196 → 0.603) while A and C moved −7%/+4%.

**Outcome — the fix did not work.** Branch grid **33 → 26 / 54**; CygOB2-B
baseline residual **3.02 → 3.62**, χ² p 0.036 → 0.0002. B's residuals sharpened
into a monotone tilt (+1.45, +2.41, +3.62, −0.36, −1.30, −0.92): its mass
function is steeper than the model across the entire window. The over-tight
repair_v1 prior had been partly *masking* this.

**Lower-edge scan also negative.** Raising B's calibration lower edge makes B
monotonically worse (3.62 → 5.11 → 5.23 at edges 2.0 / 2.5 / 2.75), which is the
opposite of a completeness-ramp artifact.

**Version recommendation (SUPERSEDED 2026-07-28 by repair_v5, which keeps the
variogram width but fixes the prior mean — see below).** Keep repair_v3 as the working extinction on the
merits of the variogram argument, *not* on gate score; picking repair_v2 because
it scores better would be selecting the extinction model by its residuals.
Caveat carried: with a ~1.05 mag prior width, B's A_V is now set almost entirely
by broadband photometry, reintroducing the Teff/A_V degeneracy the anchors were
meant to break. This is an explicit science decision, not a default.

**Repair scripts** — [wp3_extinction_repair.py](scripts/wp3_extinction_repair.py) ·
[wp3_repair_common.py](scripts/wp3_repair_common.py) ·
[wp3_repair_make_notebook.py](scripts/wp3_repair_make_notebook.py) ·
[wp3_wp5_repair_finalize.py](scripts/wp3_wp5_repair_finalize.py) ·
[wp4_mass_posteriors_repair.py](scripts/wp4_mass_posteriors_repair.py) ·
[wp4_repair_common.py](scripts/wp4_repair_common.py) ·
[wp5_injections_repair.py](scripts/wp5_injections_repair.py)

**Repair artifacts** — `data/processed/*_repair_v1.parquet` / `.npz` ·
[wp3_extinction_mass_invariant_repair_v1.csv](tables/wp3_extinction_mass_invariant_repair_v1.csv) ·
[wp5_baseline_residuals_repair_v1.csv](tables/wp5_baseline_residuals_repair_v1.csv) ·
[wp5_completeness_baseline_repair_v1.csv](tables/wp5_completeness_baseline_repair_v1.csv) ·
[wp3_wp5_repair_gates_repair_v1.png](figures/wp3_repair/wp3_wp5_repair_gates_repair_v1.png) ·
execution records `provenance/wp*_repair*_execution.json`

### Steps 1–2 of the #1c gated plan (2026-07-27, evening) · **cause identified**

| | |
|---|---|
| **Plan** | [wp5_cygob2b_age_caustic_fix_brief.md](tasks/wp5_cygob2b_age_caustic_fix_brief.md) (status header updated with results) |
| **Report** | [WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md) |
| **Gate records** | [wp5_stable_label_refit_execution.json](provenance/wp5_stable_label_refit_execution.json) (G1) · [wp5_age_scan_execution.json](provenance/wp5_age_scan_execution.json) (G2) · [wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json) (step-3 forecast) |
| **Scope** | Diagnostics only; every stored repair_v3 artifact preserved byte-identical. Snapshots `data/processed/wp5_age_scan_B_*.parquet` are new files, hashed in the scan JSON. |

**Step 1 (G1 PASSED, issue #6 closed).** Per-star label confusion over the 50
frozen GMM seeds ([wp2_label_stability.py](scripts/wp2_label_stability.py)):
1327/1331 members ≥90% seed-stable, cross-assignment ≤ 0.16%, no mass trend in
B (ρ = 0.011, p = 0.84). Stable-label refit of the full 54-branch grid moves
B's baseline bin-2 residual by **0.009** (limit 0.5); grid stays 26/54.
Contamination by mislabeling excluded.

**Step 2 (G2 PASSED).** B-only baseline injections with the truth age forced
to eight native isochrone ages, recovery side untouched, official `fit_one`
throughout ([wp5_age_conditional_scan.py](scripts/wp5_age_conditional_scan.py)):
bin-2 residual strictly monotone in truth age — **5.49 / 4.47 / 3.77 / 2.85 /
2.29 / 1.81 / 1.22 / 0.74 at 2.24 / 2.51 / 2.82 / 3.16 / 3.55 / 3.98 / 4.47 /
5.01 Myr** — with the entire tilt flattening in step; the control node (2.82 =
production truth) reproduces the stored failure (MC floor ~0.15). Nodes
**3.162 Myr (inside B's 68% interval, 2/9 posterior weight) and 3.548 Myr
(1/9) pass all three gate statistics** → the anomaly is the single-age
injection truth model sitting on the PMS/Henyey fold, displaced by an age
error, exactly finding F2/F3. B's PMS MAP itself is unmeasurable (n = 2,
grid-railed at 10 Myr); the older probes stand in for it.

**Plain truth-age marginalization is insufficient**: both posterior mixtures
fail (bin-2 3.40 / 2.89 with 1 Myr SF spread) because B's WP4 posterior is
bottom-heavy (σ_lo ≈ 0.02 vs σ_hi ≈ 0.43 Myr; 6/9 nodes on the MAP). The
brief-§4-authorized **joint age–k fit (WP4 posterior as prior, Jeffreys k
integrated analytically) passes under both prior variants** — max|res| 2.48
(χ²p 0.014) and 2.06 (χ²p 0.063), with B's joint age posterior concentrating
at 3.5–4.0 Myr, toward the PMS-indicator direction of issue #9. **repair_v4
must therefore implement the joint age–k fit uniformly for all three
subgroups**, with the anti-tuning rule intact (no gate statistic enters the
node weights).

**Steps 3a–3b (2026-07-27) — machinery implemented and validated.**
[wp5_joint_age_fit.py](scripts/wp5_joint_age_fit.py) implements the joint fit;
`inject_curve` gained an optional `truth_age_override` whose default preserves
the frozen repair_v1–v3 behaviour exactly. The node rule adds **no new
parameter** — it is the recovery side's own `age_posterior_nodes`
discretization snapped to native isochrone ages, needing **36 node responses**
for the whole 54-branch grid. Validation
([wp5_joint_fit_baseline_check_execution.json](provenance/wp5_joint_fit_baseline_check_execution.json),
zero new injections): **single-node equivalence is exact** — where a subgroup
spans one native age (A and C on the baseline) `fit_joint` reproduces the
official `fit_one` bit-for-bit, max k-draw difference **0.0**, so the new
machinery provably cannot move such a subgroup; and the **baseline branch now
passes for all three subgroups** — B's bin-2 residual 3.62 → **2.51**
(χ²p 0.0002 → 0.023) with k moving only −3.6%, A 1.36 and C 1.53 unchanged to
within Monte-Carlo noise (≤0.012 residual, ≤0.16% in k, from RNG stream
position). B's fitted truth-age posterior is 2.818: 0.02 · 3.162: 0.36 ·
3.548: 0.62.

**Step 3c is written and ready but NOT yet run** —
[run_repair_v4_chain.sh](scripts/run_repair_v4_chain.sh) →
[wp5_injections_agenodes.py](scripts/wp5_injections_agenodes.py) (33 node
injections; 3 of 36 reused from the G2 scan; measured 49 s per node) then
[wp5_fit_imf_joint.py](scripts/wp5_fit_imf_joint.py) (54-branch fit,
repair_v4 products, gate G3 with an explicit A/C per-branch regression check).
repair_v4 is **WP5-only**: WP3/WP4 repair_v3 products are consumed unchanged,
as repair_v2 consumed repair_v1.

### repair_v4 — WP5 truth-side joint age–k fit (2026-07-27) · **baseline passes**

| | |
|---|---|
| **Chain runner** | [run_repair_v4_chain.sh](scripts/run_repair_v4_chain.sh) — node injections → 54-branch joint fit |
| **Scripts** | [wp5_joint_age_fit.py](scripts/wp5_joint_age_fit.py) · [wp5_injections_agenodes.py](scripts/wp5_injections_agenodes.py) · [wp5_fit_imf_joint.py](scripts/wp5_fit_imf_joint.py) |
| **Provenance** | [wp5_injections_agenodes_execution_repair_v4.json](provenance/wp5_injections_agenodes_execution_repair_v4.json) · [wp5_imf_fit_execution_repair_v4.json](provenance/wp5_imf_fit_execution_repair_v4.json) · [wp5_report_execution_repair_v4.json](provenance/wp5_report_execution_repair_v4.json) |
| **Report** | [wp5_imf_norm_repair_v4.md](wp5_imf_norm_repair_v4.md) (data-driven; `wp5_report.py --wp5-version`) |
| **Scope** | **WP5 only.** WP3 extinction and WP4 ages/masses are the repair_v3 products, consumed unchanged. Every repair_v1–v3 artifact preserved byte-identical. |

**The fix.** The injection truth age is no longer a single upper-MS MAP. The
forward response is the WP4-posterior mixture over truth-age nodes, with node
weights updated by the Poisson likelihood of the observed counts (Jeffreys k
integrated analytically). Node rule adds **no new parameter** — it is the
recovery side's own `age_posterior_nodes` discretization snapped to native
isochrone ages. Identical machinery for all three subgroups; no gate threshold
moved; no per-subgroup age chosen. With a single node the estimator reduces to
the unmodified `fit_one` **bit-for-bit**.

**Result — the baseline branch passes for all three subgroups:**

| subgroup | max abs residual v3 → v4 | χ² p v3 → v4 | k v3 → v4 | pass |
|---|---:|---:|---:|:--:|
| CygOB2-A | 1.37 → 1.45 | 0.350 → 0.380 | 1724 → 1708 | yes |
| **CygOB2-B** | **3.62 → 2.51** | **0.0002 → 0.023** | 1925 → 1855 | **yes (was no)** |
| CygOB2-C | 1.52 → 1.52 | 0.341 → 0.492 | 1884 → 1873 | yes |

Branch grid **26 → 29 / 54**. Association mass 30,155 [28,987, 31,312] M☉,
within a factor two of the 16,500 M☉ literature scale. B's fitted truth-age
posterior mean is 3.40 Myr on the baseline (3.16–3.96 across branches), above
its 2.82 Myr upper-MS MAP and toward the PMS-indicator direction of issue #9.

**Open verdict question — the "A and C must not regress" clause of gate G3.**
Two CygOB2-C cells (MIST, R_V = 3.1, α = 2.0 and 2.3) flip pass → fail, so the
per-branch reading of that clause evaluates **false** while the baseline
reading passes. The flips are **not** caused by the fix, and this is measured,
not asserted ([wp5_trend_stability_check_execution.json](provenance/wp5_trend_stability_check_execution.json)):

- Of the 7 pass/fail flips between v3 and v4, **3 are the model change** (all
  fail → pass, driven by large max-residual improvements: B baseline α = 2.3
  3.63 → 2.51, α = 2.6 3.77 → 2.66, C MIST R_V = 3.5 α = 2.0 3.03 → 2.27).
- The other **4 are single-truth-age-node cells where the estimator is provably
  identical to repair_v3** — 2 flip each way (A MIST R_V = 3.1 α = 2.0 and
  A MIST R_V = 3.5 α = 2.6 gain; C MIST R_V = 3.1 α = 2.0 and α = 2.3 lose).
  In all four the χ² and max residual barely move; only the rank-based trend
  statistic crosses 0.05.
- **Paired test:** fitting the *identical* estimator to the two independent
  Monte-Carlo realizations of the *same* model (the repair_v3 injection
  response and the repair_v4 node response, generated at the same truth age)
  gives **opposite gate verdicts** in all four cells, and `trend_p` is
  completely insensitive to the fit RNG seed. The flips are injection
  Monte-Carlo noise.
- **Mechanism:** with 6 bins the two-sided Spearman p-value lives on a coarse
  lattice (ρ = 0.771 → p = 0.072; 0.829 → 0.042; 0.943 → 0.005), so a small
  realization change in a *well-fitting* residual vector can cross p = 0.05
  while χ² and max residual are unchanged or better. The instability is
  unbiased — it moved two cells each way.

This is a genuine weakness of the trend diagnostic, not of repair_v4. Resolving
it is a science decision recorded in §10; **no threshold was moved and
`downstream_wp6_authorized` has not been flipped.**

### repair_v5 — kriged anchor prior (2026-07-28) · **ADOPTED WORKING VERSION**

| | |
|---|---|
| **Decision record** | [wp3_kriging_adoption.json](provenance/wp3_kriging_adoption.json) |
| **Report** | [WP3_KRIGED_PRIOR_repair_v5.md](reports/WP3_KRIGED_PRIOR_repair_v5.md) |
| **Method spec** | [CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §15 |
| **Pre-registration / outcome** | [wp3_kriging_prior_prereg.json](provenance/wp3_kriging_prior_prereg.json) · [wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json) |
| **Downstream impact** | [wp3_kriging_downstream_impact.json](provenance/wp3_kriging_downstream_impact.json) |
| **Gate record** | [wp5_repair_v5_gate.json](provenance/wp5_repair_v5_gate.json) · report [wp5_completion_report_repair_v5.md](wp5_completion_report_repair_v5.md) |
| **Chain runner** | [run_repair_v5_chain.sh](scripts/run_repair_v5_chain.sh) — full WP3 → WP5 |
| **Scope** | Full chain. repair_v1–v4 artifacts preserved byte-identical; `WP3_ANCHOR_PRIOR_MODE=variogram\|global` still reachable. |

**The change.** The anchor spatial prior's **mean** is now a simple-kriging
estimate using the variogram *already fitted in repair_v3*, replacing the plain
median of the 8 nearest anchors. The **width is unchanged**. No new parameter.
The prior previously widened its uncertainty toward the sill for distant-anchor
stars — declaring those anchors uninformative — while still centring on them at
full weight; the two moments contradicted each other.

**Why the prior and not the photometry was at fault** — four independent lines:
the misfit is a *pure* A_V offset (per-band residual ÷ extinction coefficient
constant to <0.1 mag); optical-only and near-IR-only A_V agree to 0.04–0.07 mag
with B *not* the outlier; the prior-vs-photometry disagreement grows with anchor
distance **only in B** (ρ = +0.42, p = 1×10⁻¹⁸, versus −0.13 for A and −0.08
n.s. for C); and both 3D dust cubes rank B *least* extinguished while the
uncorrected prior ranked it *most*.

**Pre-registered, all four predictions confirmed:** B's A_V −0.359 mag; B's age
2.818 → **3.548 Myr** (68% [3.389, 3.991], reaching the age at which WP5's mass
function independently passed); A and C ages moved **0.000 Myr**; B's tilt
statistic |T| 2.40 → **1.31**.

**Result.** B's bump-bin residual **3.62 → 2.51 → 0.36** across v3/v4/v5 and the
high-mass deficit −1.30/−0.92 → −0.03/+0.17. Branch grid **26 → 29 → 38 of 54**.
`repair_v5` is the **first version whose baseline passes for all three subgroups
under both trend statistics**. Together with repair_v4 this closes finding F1 of
the fix brief: the *bump* was the isochrone fold under a single-age truth model,
the *mass-scale displacement* was this extinction error.

**Still not accepted.** Three A/C non-baseline cells regress (A PARSEC R_V=3.1
α=2.0 at trend p 0.040; A MIST R_V=3.5 α=2.6 at 0.048; C MIST R_V=3.5 α=2.0 at
max|r| 3.25). Two sit inside issue #11's indeterminate band.
`downstream_wp6_authorized` stays **false**.

**Downstream effects.** Association mass 30,155 → 29,185 M☉ (**−3.2%**, inside
the ±18% branch spread, so N_SN moves ~3%). k: A +2.2%, **B −12.2%**, C +1.0%.
The consequential change is structural, not numerical: the age spread is
unchanged (1.469 Myr) but B moves from near C to near A, so the association
reads as **two older subgroups plus one younger** rather than one older plus
two younger — a different star-formation history, which WP7 inherits.
Anchor selection bias was tested and is **absent** (anchors read +0.500 mag
*higher* than co-located members' prior-free photometry, p = 4.3×10⁻¹⁶).

**Obligations accepted with adoption:** O1 report the revised SF history as a
result; O2 carry the anchor absolute-scale systematic (photometry sits ~0.5 mag
below anchors, so the mass scale and N_SN inherit the anchor calibration);
O3 carry B's calibration asymmetry (4 anchors vs 59 and 42; Gaia XP spectra the
natural remedy); ~~O4 resolve issue #11's R3 before WP5 acceptance~~ —
**O4 discharged 2026-07-28**, see
[WP5_TREND_STATISTIC_R3_RESOLUTION.md](reports/WP5_TREND_STATISTIC_R3_RESOLUTION.md).

**Correction to the "three regressing cells" line above (2026-07-28).** Issue
#11's resolution re-classified them. The two CygOB2-A cells are **indeterminate
in both `repair_v4` and `repair_v5`** — their verdict is a coin flip on the
injection realization (π 0.46 → 0.30 and 0.79 → 0.78) — and are discounted under
[CUTS §14.7](CUTS_AND_THRESHOLDS.md). **One determinate regression survives:**
CygOB2-C MIST R_V = 3.5 α = 2.0, determinate pass in v4 (π = 1.000), and it is
issue #13's node-snapping artifact rather than a physical failure.

### repair_v6 — truth-age isochrone interpolation (2026-07-28) · **ACCEPTED WORKING VERSION**

| | |
|---|---|
| **Pre-registration** | [wp5_node_interpolation_prereg.json](provenance/wp5_node_interpolation_prereg.json) — written before any injection |
| **Mechanism evidence** | [wp5_node_rule_continuity_execution.json](provenance/wp5_node_rule_continuity_execution.json) · [wp5_node_rule_continuity.py](scripts/wp5_node_rule_continuity.py) |
| **Report** | [WP5_NODE_INTERPOLATION_repair_v6.md](reports/WP5_NODE_INTERPOLATION_repair_v6.md) |
| **Outcome** | [wp5_node_interpolation_outcome.json](provenance/wp5_node_interpolation_outcome.json) — all four predictions confirmed |
| **Gate record** | [wp5_repair_v6_gate.json](provenance/wp5_repair_v6_gate.json) · report [wp5_completion_report_repair_v6.md](wp5_completion_report_repair_v6.md) |
| **Chain runner** | [run_repair_v6_chain.sh](scripts/run_repair_v6_chain.sh) — WP5 only |
| **Backward compatibility** | [wp5_v6_backward_compatibility.json](provenance/wp5_v6_backward_compatibility.json) — all four checks pass |
| **Scope** | **WP5-only.** WP3 extinction and WP4 ages and masses consumed from `repair_v5` unchanged, exactly as `repair_v4` consumed `repair_v3`. repair_v1–v5 artifacts byte-identical. |

**The change.** The nine WP4 age-posterior quantiles are no longer snapped to
native isochrone ages. Each keeps weight 1/9 where the posterior puts it, and
`wp5_common.load_isochrone_between_ages` builds the truth isochrone there using
the **same bracketing and linear age blend the recovery side already used**
(`wp4_repair_common._interpolate_age_sequence`). The two sides had been
inconsistent — the recovery side interpolated while the truth side snapped. No
new parameter: the node count is still `N_AGE_NODES = 9`.

**Why.** Measured, not inferred. Translating a branch's WP4 posterior by
**0.0005 Myr** moved the snapped node distribution by up to **0.0546 Myr** in
1-Wasserstein distance — a **109× amplification** — and the node count varied
between 1 and 4 within a single branch under a ±0.1 Myr translation. The
unsnapped rule is **Lipschitz-1**: it moves by exactly the shift that caused it,
the tightest continuity a rule can have. The interpolated loader was verified to
reproduce the native isochrone to **0.00e+00 mag** at every native age, so at a
node that happens to land on the grid nothing changes at all.

**Cost.** 37 → **162** node injections. No reuse is possible because every node
now sits at a new age.

**Pre-registered predictions — all four confirmed.** P1 continuity (established
before the run, Lipschitz-1); **P2** the surviving regression clears, max|r|
**3.25 → 2.94** — explicitly *not* the adoption criterion, since it is the cell
that motivated the work; **P3** the baseline passes all three subgroups under
both statistics with B's |T| = **1.43**; **P4** the change is local — worst
truth-age posterior mean shift **0.118 Myr** (limit 0.15, passing but not by a
wide margin), association mass **−0.22%**, grid +3. **Adoption rested on P3 and
P4, which bound the size of the change without prescribing its direction.**

**Result — WP5 accepted.** Gate G3 passes with `A_or_C_regressions: []` under
the **strict** per-branch reading; §14.7's refinement discounted nothing and was
not needed (`readings_agree: true`). `accepted: true`,
`blocking_reason: null`, **`downstream_wp6_authorized: true`**. Branch grid
38 → **40 of 54**. The baseline barely moved — a nine-node mixture reproduced
what a two-node mixture found, so repair_v5's physics was not an artifact of the
coarse approximation. k moved −1.1% / +0.6% / −0.3% (A/B/C), association mass
29,185 → **29,122 M☉**. Truth-age nodes per baseline branch 2/4/1 → **9/9/9**.

**Unanticipated benefit.** The indeterminate fraction of the branch grid — 41–44%
in every earlier version, which had triggered §14.7(5) as an underpowered-grid
finding — falls to **17%** (9 of 54), below the 25% limit. Mechanism, stated
rather than treated as luck: a nine-node mixture averages the forward response
over nine independent injection realizations instead of two, cutting the
response's Dirichlet uncertainty by roughly √(9/2). The finding stands for
repair_v3–v5 and **does not apply to the accepted version**.

**Backward compatibility proven** ([wp5_v6_backward_compatibility.json](provenance/wp5_v6_backward_compatibility.json)):
the default node rule reproduces the exact repair_v4/v5 node sets on disk; the
interpolating loader equals the native table to **0.00e+00 mag** at all 35
native ages; the vectorized Jeffreys k is bitwise identical; the interpolation
flag is off for every pre-v6 version.

---

## 9. Open issues register

| # | Issue | Severity | Owner WP | Status |
|---|---|---|---|---|
| 1 | Forward response truncated parent mass distribution at 8 M☉; top bin under-predicted in all three subgroups | ~~blocking~~ | WP5 | **CLOSED 2026-07-27** by repair_v2; top-bin residuals fell in all three subgroups |
| 1b | WP3 anchor spatial prior over-constrains CygOB2-B: width calibrated at 0.071° applied at 0.377°, understated 2.32×, collapsing B's A_V spread 9.5× | high | WP3 | **CLOSED 2026-07-27** by repair_v3 variogram prior. Real bug, correctly fixed — but *not* the cause of #1c |
| 1c | **CLOSED 2026-07-28.** Both causes found and fixed: the single-age injection truth model on the PMS/Henyey fold (`repair_v4`) and the anchor-prior extinction-scale error (`repair_v5`). B's baseline worst residual went **3.62 → 2.51 → 0.36** in the bump bin, the tilt statistic \|T\| **3.22 → 2.40 → 1.31**, and the grid **26 → 29 → 38** of 54. `repair_v5` is the first version whose baseline passes for all three subgroups under both trend statistics. WP5 acceptance is now held by issue #11, not by #1c. Earlier status: **PARTIALLY SOLVED — the bump is fixed, a tilt remains.** ~~CAUSE FOUND AND FIXED~~ (that claim, made 2026-07-27, was **corrected on 2026-07-28** and was too strong). One real cause was identified and removed: the single-age injection truth model sitting on the PMS/Henyey isochrone fold. The fix (truth-side joint age–k marginalization, `repair_v4`) takes CygOB2-B's baseline residual **3.62 → 2.51**, χ² p 0.0002 → 0.023, grid 26 → 29/54. **But a properly powered trend test shows B's residuals still drift systematically with mass** (+1.20, +2.08, +2.51, −0.24, −0.79, −0.54; replacement trend p = **0.017**, T = −2.40, versus A −0.39 and C +0.19 which are flat). The incumbent rank test scored this 0.111 only because the sequence is not *monotone*. This matches finding F1 of the fix brief, which decomposed B's anomaly into **a localized bump plus a global mass-scale displacement** — the age fix removed the bump; the mass-scale displacement was never addressed. Next hypothesis is therefore well-posed and already written down: issue #1d (B's extinction is set by broadband photometry alone, so a scale error displaces B's masses coherently) and step 4b of the gated plan. Original statement: CygOB2-B's observed mass function is steeper than the model across the whole 2–8 M☉ window (bin-2 residual 3.62; monotone tilt +1.45/+2.41/+3.62/−0.36/−1.30/−0.92). Robust to every instrumental explanation tested. **2026-07-27 shape diagnostic** ([wp5_bump_shape_diagnostic_execution.json](provenance/wp5_bump_shape_diagnostic_execution.json)): no single IMF slope explains B (bin-2 excess survives all α at R_V ≥ 3.1); the bump sits on the PMS/Henyey isochrone fold (2.77–3.22 M☉ at B's MAP age); injections generate truth at a single UMS MAP age and are structurally blind to an age error; contamination signatures absent in B's bin-2 stars | **blocking** | WP5←WP4 | **cause identified 2026-07-27** — steps 1–2 of the gated plan executed ([WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md)): **G1 passed** (labels exonerated, shift 0.009; contamination excluded), **G2 passed** — B's residual is strictly monotone in the injection truth age (bin-2: 5.49 → 0.74 over 2.24 → 5.01 Myr) and ages inside B's WP4 posterior (3.162, 3.548 Myr) pass the full gate ([wp5_age_scan_execution.json](provenance/wp5_age_scan_execution.json)). Plain truth-age marginalization fails (posterior bottom-heavy); the §4-authorized **joint age–k fit passes both prior variants** ([wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json)). **Steps 3a–3c executed** — `repair_v4` chain complete ([wp5_repair_v4_gate.json](provenance/wp5_repair_v4_gate.json), [wp5_completion_report_repair_v4.md](wp5_completion_report_repair_v4.md)); `accepted: false`, `downstream_wp6_authorized: false` under the adopted strict G3 reading. The blocker is now **issue #11**, not the CygOB2-B mass function |
| 1d | **FIXED 2026-07-28 in `repair_v5`** by the kriged prior mean (`ANCHOR_PRIOR_MODE=kriging`): the prior now shrinks toward the field mean as anchor information decays, exactly as the fitted variogram says it should. Kriging weights sum to 1.000 for A and 0.992 for C but only **0.772 for B**, so the correction is intrinsically B-specific and no per-subgroup choice was made. B's A_V fell 0.359 mag and its age rose to 3.548 Myr while A and C did not move. The underlying observational limitation stands as a carried caveat — B still has only 4 spectroscopic anchors against A's 59 and C's 42, so its extinction rests on weaker calibration than its siblings', and Gaia XP spectra remain the natural way to remove that asymmetry. Original statement: with the variogram prior, CygOB2-B's A_V is set almost entirely by broadband photometry, reintroducing the Teff/A_V degeneracy the anchors were meant to break. No DR3 spectroscopic anchors lie near B's footprint | **high** (was medium) | WP3 | **CONFIRMED WITH A MEASUREMENT 2026-07-28** ([wp3_nir_extinction_crosscheck_execution.json](provenance/wp3_nir_extinction_crosscheck_execution.json), [wp3_nir_extinction_crosscheck.py](scripts/wp3_nir_extinction_crosscheck.py)). Anchor counts inside the member sample are **A 59, C 42, B 4** — the anchor famine is real and B-specific. And B's WP3 A_V sits **+0.385 mag above** what its near-IR colour excess implies, relative to A and C (Mann-Whitney p = 2.4×10⁻¹², n=334 vs 684 in the 2–8 M☉ window). The offset keeps its sign and significance at *every* candidate age for B (+0.385 at 2.82 Myr → +0.265 at 4.47 Myr), so it is not an artifact of B's disputed age. **Mechanism corrected 2026-07-28** ([wp3_band_tension_diagnostic_execution.json](provenance/wp3_band_tension_diagnostic_execution.json)): my first reading — that the optical bands were outvoting the near-IR ones — is **wrong and is retracted**. Inside the fit, the A_V implied by G/BP/RP alone and by J/H/Ks alone agree to 0.04–0.07 mag in every subgroup, and B is *not* the outlier (B +0.044 vs A +0.074). The misfit is instead a **pure A_V offset**: per-band residuals are proportional to the extinction coefficient (residual/k constant to <0.1 mag across all six bands), meaning the model has the colours right and only the extinction wrong. Its source is the **anchor spatial prior**, whose mean sits above what the photometry alone prefers by **+0.48 mag (A), +0.71 mag (B), −0.10 mag (C)**, pulling the adopted A_V up by +0.415 / +0.485 / −0.080 respectively. B's eighth-nearest anchor is **0.374°** away against 0.089° for A and 0.139° for C, so B's prior is the least locally informed. **Caveat, stated plainly:** A is pulled almost as hard as B and still passes the WP5 gate, so the prior pull is necessary-but-not-sufficient and does not by itself establish causation for B |
| 2 | [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) / `.md` are the pre-repair run (0/54) and contradict the current fit. `wp5_report.py` is hardwired to the unversioned frozen products and cannot regenerate from a repair version | high | WP5 | **CLOSED 2026-07-27** — `wp5_report.py` gained `--wp5-version`, and the versioned path is fully **data-driven**: every pass/fail statement is read from the run's own gate record instead of the hardcoded blocked-gate prose. Versioned products exist for repair_v3 and repair_v4 ([wp5_imf_norm_repair_v4.md](wp5_imf_norm_repair_v4.md), [tables/wp5_imf_norm_repair_v4.csv](tables/wp5_imf_norm_repair_v4.csv), matching figures). The unversioned files are left untouched as the frozen historical record of the 0/54 run |
| 3 | **CLOSED 2026-07-28 — measured, corrected, and the correction turned out to differ from the flagged estimate in both size and sign.** The defect was real: recovery is flat in mass out to 18 M☉ (plateau **0.787** grid median), so the loss is the WP2 quality filter and does not vanish at the bright end. But the naive repair — dividing the observed count by that plateau — is **forbidden**, because the response also scatters mass estimates across the 8 M☉ threshold. The correct forward estimator is `k · ∫ dM M^(−α) R(observed above threshold \| M)`, the same construction WP5 already uses inside 2–8 M☉ extended to one open-ended bin. Its **effective completeness is 0.871** (grid median), range 0.808–**1.041**; in **6 of 54 cells it exceeds 1**, so net up-scatter beats recovery loss and the scalar divisor has the *wrong sign*. Spurious deficit if uncorrected: **12.9%**, not the ~17–20% originally flagged. A per-subgroup systematic also emerged and is now binding: **CygOB2-B sits at 0.969 against A 0.840 and C 0.852**, because B's wider mass posteriors (obligation O3) scatter more stars up — a single association-wide completeness would under-count B's explosions. Specification is binding in [CUTS_AND_THRESHOLDS.md §16](CUTS_AND_THRESHOLDS.md); must be re-run on whichever WP5 version is accepted. Original statement: bright-mass completeness plateaus at ~0.83; WP6 step 2(a) assumes ~1.0 and will manufacture a false ~17% deficit → biases N_SN high | ~~high~~ | WP6 | **CLOSED** — [wp6_bright_completeness_execution.json](provenance/wp6_bright_completeness_execution.json) · [wp6_bright_completeness.py](scripts/wp6_bright_completeness.py) |
| 4 | No subgroup reaches an absolute 95% completeness edge; `CUTS_AND_THRESHOLDS.md` §7.1 procedure is unachievable on this field | medium | WP5 | **CLOSED 2026-07-27** — formal supersession recorded in §7.1; `corrected_no_absolute95_edge` fallback in force |
| 4b | Plan WP5 step 2 / §7.1 prescribe raising the lower edge when completeness fails. Measured twice — global edges on repair_v2 (30/28/33 of 54) and per-subgroup completeness-driven edges on repair_v3 (25/28/29 of 54). Neither works; raising B's edge makes B monotonically worse. The prescribed remedy does not apply on this field | medium | WP5 | measured and documented in [wp5_lower_edge_scan_execution.json](provenance/wp5_lower_edge_scan_execution.json); do not retry blind |
| 5 | Branches failing on the top bin now confined to the R_V = 3.5 and α = 2.6 corners; plan §1.4 forbids silently dropping branches | medium | WP5 | **CLOSED 2026-07-27** — retention policy adopted as [CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §13, written *before* the repair_v4 grid was evaluated so it could not be tailored to the outcome: no branch is ever dropped; the baseline must pass; non-baseline failures are carried as a systematic statement and reported by their concentration; downstream consumes the ensemble, never the best-scoring branch |
| 6 | Subgroup-label (A vs B vs C) uncertainty never quantified — `membership_probability` measures cluster-vs-field only | medium | WP2 | **CLOSED 2026-07-27** — quantified per star over the 50 frozen GMM seeds ([wp2_label_stability.py](scripts/wp2_label_stability.py) → [wp2_label_stability_execution.json](provenance/wp2_label_stability_execution.json), per-star CSV alongside): 1327/1331 clean members are ≥90% seed-stable, seed-averaged cross-assignment ≤ 0.16%, and B shows **no mass-dependent instability** (window Spearman ρ = 0.011, p = 0.84; 341/342 window stars stable). The stable-label WP5 refit ([wp5_stable_label_refit_execution.json](provenance/wp5_stable_label_refit_execution.json)) moves B's baseline bin-2 residual by **0.009** (G1 limit 0.5) and leaves the grid at 26/54 — **labels exonerated; gate G1 of the #1c plan passed**, contamination-by-mislabeling excluded as the bump cause |
| 7 | A/B/C age ordering unresolved → star-formation-duration branch (0/1/2 Myr) stays mandatory | medium | WP4→WP7 | carried by design |
| 8 | Berlanas two-distance population not confirmed; single population at 1.62 kpc, depth 45 pc | low | WP2 | documented result |
| 9 | Upper-MS vs PMS indicator disagreement in subgroup B (N = 19 tail) | low | WP4 | documented; **2026-07-27**: the #1c joint age–k forecast independently pulls B's 2–4 M☉ truth age to 3.5–4.0 Myr — the same direction as the PMS indicator ([wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json)) |
| 10 | Wide-box definition inconsistency noted in `CUTS_AND_THRESHOLDS.md` §4.2 | low | WP1→WP6 | resolve before runaway traceback |
| 11 | **CLOSED 2026-07-28 — the criterion was measuring the wrong thing.** Criterion R3 required both realizations of an identical model to return the same verdict. That is **unsatisfiable by any statistic**: a cell with pass-probability π disagrees across independent realizations with probability 2π(1−π), so only a degenerate test could comply. R3 is **withdrawn** and replaced by R3a/R3b/R3c ([CUTS §14.6–14.7](CUTS_AND_THRESHOLDS.md), pre-declared before measurement); R1, R2, R4 stand. **R3a PASS** — the incumbent has no achievable p between **0.0416 and 0.0724** (gap 0.031 straddling the threshold), the replacement's gap is 0.00005. **R3b PASS, and it independently convicts the incumbent** — using each cell's own injection uncertainty (Dirichlet response posterior, 400 replicates, no new parameter) the implied flip count is 1.43 [0,3], and the incumbent flipped **4 of 4**, more often than its own noise can explain; the replacement flipped 1 against 0.74 expected. **R3c: implemented, and it triggers §14.7(5)** — **41–44% of the 54-cell grid is indeterminate** (0.05 < π < 0.95) in *every* version, so the branch grid is underpowered at N_INJECT_PER_MASS = 400 and per-cell verdicts are not independent evidence. Determinate failures nevertheless fall **15 → 11 → 5** across v3/v4/v5, so repair_v5's improvement is not an artifact. **Consequence for G3:** the two A cells at trend p 0.040 and 0.048 are indeterminate in *both* v4 and v5 and are discounted; **one regression survives — CygOB2-C MIST R_V=3.5 α=2.0**, determinate pass in v4 (π = 1.000) and indeterminate in v5 (π = 0.900). That cell is issue #13, exactly as forecast, and `repair_v6` cleared it. **Postscript:** the underpowered-grid finding this raised (§14.7(5)) does **not** apply to the accepted version — `repair_v6`'s indeterminate fraction is **17%**, against 41–44% for v3–v5, because a nine-node truth-age mixture averages the response over nine injection realizations instead of two. Prior status: **NOW THE SOLE WP5 BLOCKER**. Original statement: **The WP5 residual-trend gate statistic is Monte-Carlo unstable.** With 6 bins the two-sided Spearman p-value lives on a coarse lattice (ρ = 0.771 → p = 0.072; 0.829 → 0.042; 0.943 → 0.005), so an independent injection realization of an *identical* model can move a well-fitting cell across p = 0.05 while χ² and max residual are unchanged or better. Measured directly: 4 of the 7 repair_v3→v4 gate flips are single-truth-age-node cells whose estimator is provably identical, and a paired refit on the two realizations gives opposite verdicts (2 flips each way, so the instability is unbiased). `trend_p` is insensitive to the fit RNG seed — the driver is injection noise | ~~blocking~~ | WP5 | **RESOLVED** — [WP5_TREND_STATISTIC_R3_RESOLUTION.md](reports/WP5_TREND_STATISTIC_R3_RESOLUTION.md) · [wp5_verdict_stability_execution.json](provenance/wp5_verdict_stability_execution.json) · [wp5_verdict_stability.py](scripts/wp5_verdict_stability.py) · prior evidence: [wp5_trend_stability_check_execution.json](provenance/wp5_trend_stability_check_execution.json) · [wp5_trend_stability_check.py](scripts/wp5_trend_stability_check.py). Fix by **replacing the diagnostic, never by moving the threshold** ([CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §6.4): done and evaluated — see the status column. Raising `N_INJECT_PER_MASS` reduces the noise but not the lattice coarseness |
| 13 | **CLOSED 2026-07-28 by `repair_v6`.** The fix was pre-registered and all four predictions confirmed: the blocking cell CygOB2-C MIST R_V=3.5 α=2.0 cleared (max\|r\| **3.25 → 2.94**), the baseline held under both statistics (B's \|T\| = 1.43), and the change stayed local (worst truth-age posterior mean shift 0.118 Myr against a 0.15 limit, association mass −0.22%, grid 37 → 40). Gate G3 then passed with `A_or_C_regressions: []` and **WP5 was accepted**. Report: [WP5_NODE_INTERPOLATION_repair_v6.md](reports/WP5_NODE_INTERPOLATION_repair_v6.md). Prior status: **The truth-age node set changes discontinuously under age-grid snapping.** The repair_v4 node rule snaps the nine WP4 posterior quantiles to native isochrone ages, whose grid is coarse (0.05 dex ≈ 12%), so an arbitrarily small change in the age posterior can add or delete a whole node carrying tens of percent of the weight. Measured: CygOB2-C MIST R_V=3.5 α=2.0 kept an identical fitted age (2.002 → 2.002 Myr) and moved only −0.059 mag in A_V, but its 68% upper edge fell 2.084 → 2.062 Myr, which deleted the 2.248 Myr node holding **41.6%** of the weight — and bin 0's residual jumped +2.27 → **+3.25**. This is the same class of defect as issue #11: a coarse discretization turning a continuous input into a discontinuous output | ~~blocking~~ | WP5 | **DIAGNOSED, FIXED, PRE-REGISTERED, RUNNING.** The discontinuity is now measured rather than inferred ([wp5_node_rule_continuity_execution.json](provenance/wp5_node_rule_continuity_execution.json), [wp5_node_rule_continuity.py](scripts/wp5_node_rule_continuity.py)): translating the WP4 posterior by **0.0005 Myr** moves the snapped node distribution by up to **0.0546 Myr** in 1-Wasserstein distance — a **109× amplification** — and the node count varies between 1 and 4 within a single branch under a ±0.1 Myr translation. Dropping the snap makes the rule **Lipschitz-1**: it moves in W1 by exactly the shift that caused it, the tightest continuity a rule can have. Fix implemented as `repair_v6`: `wp5_common.load_isochrone_between_ages` builds the truth isochrone by the same bracketing and linear age blend the recovery side already uses (verified to reproduce the native table to **0.00e+00 mag** at every native age), and `truth_age_nodes(snap=False)` keeps the nine quantiles where the posterior puts them. No new parameter; repair_v1–v5 are bit-preserved behind default arguments. Predictions **pre-registered before the run** in [wp5_node_interpolation_prereg.json](provenance/wp5_node_interpolation_prereg.json) — P2 the regressing cell clears (explicitly **not** the adoption criterion), P3 the baseline still passes under both statistics, P4 the change is local (posterior mean age < 0.15 Myr, mass ±5%, grid ±6 cells). 162 node injections running (37 → 162; no reuse possible, every node is at a new age). Confirmed by issue #11's independent analysis: this cell is the one surviving G3 regression, determinate pass in v4 (π = 1.000) and indeterminate in v5 (π = 0.900). Prior status: open — [wp5_node_snapping_discontinuity.json](provenance/wp5_node_snapping_discontinuity.json). Candidate fix: interpolate the truth isochrone between native ages instead of snapping. Note the two sides are currently **inconsistent** — `wp4_repair_common._interpolate_age_sequence` already interpolates on the *recovery* side while the truth side snaps. **Re-specifying R3 will not clear this cell**, so #11 and #13 are two separate prerequisites for a fair G3 re-evaluation |
| 12 | **CLOSED 2026-07-28 — cause confirmed and fixed in `repair_v5`.** The extinction-scale error was real and was the cause: correcting the anchor prior's mean (kriging) flattened B's tilt, \|T\| 2.40 → **1.31** (trend p 0.017 → 0.186), with bin-2 residual 2.51 → **0.36** and the high-mass deficit −0.79/−0.54 → −0.03/+0.17. Pre-registered and confirmed ([wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json)). Original statement: **an independently measured extinction-scale error in B.** The near-IR cross-check (issue #1d) finds B's optical A_V high by +0.385 mag relative to A and C, which rescales B's whole mass ladder by **12.6%**; by contrast the disputed age shift moves assigned mass by only 0.4%, because it mainly relocates the PMS/Henyey fold rather than rescaling masses. That is precisely finding F1's decomposition of B's anomaly into *a localized bump plus a global mass-scale displacement*: `repair_v4` removed the bump (age/fold), and the displacement is the extinction offset. The near-IR reddening-law **shape** is the same for all three subgroups (E(J−H)/E(H−Ks) = 1.21–1.28), so B does not need its own R_V — this disfavours step 4b and favours a scale error. **Not yet demonstrated:** that correcting the extinction actually flattens the tilt; the injections pass through the same biased estimator so the bias partly cancels, and the sign of the net effect must be established by running the chain, not by argument. Original statement: **CygOB2-B's mass function still has a mass-dependent tilt after the age fix.** Under the better-powered replacement trend statistic, B's `repair_v4` baseline residuals (+1.20, +2.08, +2.51, −0.24, −0.79, −0.54) give trend p = **0.017** (T = −2.40), against A −0.39 and C +0.19 which are flat — so B's baseline **fails** the gate under that statistic even at `repair_v4`. The incumbent rank test scored it 0.111 only because the sequence rises then falls instead of increasing monotonically. This is the *global mass-scale displacement* half of finding F1 in [wp5_cygob2b_age_caustic_fix_brief.md](tasks/wp5_cygob2b_age_caustic_fix_brief.md): the age fix removed the localized bump, but the displacement was never addressed | **blocking** | WP5←WP3 | open — leading hypothesis is issue #1d: with no spectroscopic anchors near B, its A_V is set by broadband photometry alone, so an extinction-scale error moves B's whole mass scale coherently. Step 4b of the gated plan (continuous R_V / extinction-scale scan for B, to be justified by **independent extinction evidence, never gate score**) is the pre-written next experiment |

---

## 10. Next actions

*Updated 2026-07-28 after WP5 acceptance. Every item that stood here this
morning is closed: issues #3, #11 and #13 are resolved, obligation O4 is
discharged, and gate G3 passed under the strict reading.*

**WP5 is accepted and WP6 is authorized.** The work below is WP6's, not WP5's.

1. **Start WP6 census closure.** The estimator is already binding in
   [CUTS §16](CUTS_AND_THRESHOLDS.md) and its inputs are measured on the
   accepted version: expected observed count above threshold is
   `k · ∫ dM M^(−α) R(observed above threshold \| M)`, **per subgroup**, never a
   scalar completeness. Effective completeness above 8 M☉ is **0.872** grid
   median, and CygOB2-B's **0.962** against A's 0.834 and C's 0.846 is a real
   systematic — a single association-wide value would under-count B's
   explosions.
2. **Resolve issue #10 before the runaway traceback** — the wide-box definition
   inconsistency in [CUTS §4.2](CUTS_AND_THRESHOLDS.md) sets a velocity ceiling
   that the 10–100 km/s runaway window may exceed.
3. **Report the 14 branches that still fail** under the
   [§13 retention policy](CUTS_AND_THRESHOLDS.md). The grid is 40/54; the
   failures are concentrated in the R_V = 3.5 and α = 2.0/2.6 corners and must
   be shown, not dropped.
4. **Discharge the remaining adoption obligations** from
   [wp3_kriging_adoption.json](provenance/wp3_kriging_adoption.json): **O1**
   report the revised star-formation history as a result — the association now
   reads as two older subgroups (A 4.00, B 4.07 Myr) plus one younger
   (C 2.52 Myr), which WP7's supernova timeline inherits; **O2** carry the
   anchor absolute-scale systematic (photometry sits ~0.5 mag below
   spectroscopic anchors, so the absolute mass scale and N_SN rest on the
   anchor calibration); **O3** carry CygOB2-B's calibration asymmetry (4 anchors
   against 59 and 42 — Gaia DR3 XP spectra are the natural remedy and are **not
   yet ingested**, WP1 carries no astrophysical-parameter columns). **O4 is
   discharged.**
5. **Carried methodological finding, for the paper's systematics section.**
   At `N_INJECT_PER_MASS = 400` a single-node truth model leaves **41–44%** of
   the 54-cell branch grid with an indeterminate verdict — a pass/fail that is a
   coin flip on the injection realization. `repair_v6`'s nine-node mixture cuts
   this to **17%** by averaging the response over nine realizations. Any future
   branch-grid comparison must report the indeterminate set rather than treating
   per-cell verdicts as independent evidence.
6. **Open decision, unforced.** [CUTS §14.7](CUTS_AND_THRESHOLDS.md) refines the
   strict G3 reading chosen on 2026-07-27 so that a cell indeterminate in *both*
   versions counts as neither regression nor improvement. It was **not needed**
   for `repair_v6` — both readings agree and nothing was discounted — so it
   stands as a pre-declared rule for future comparisons rather than a decision
   taken under pressure.

## 11. Conventions used throughout

- **Nothing is overwritten.** Failed runs and superseded artifacts are kept with
  `_failed_<date>` or `_repair_v<n>` suffixes and stay referenced from the
  provenance logs.
- **Every gate verdict is a JSON record**, not just prose — `provenance/*_gate.json`,
  `*_validation.json`, `*_manifest.json` carry sha256 digests of every input and
  output.
- **Model branches are carried, never averaged** — isochrones {PARSEC, MIST} ×
  R_V {3.0, 3.1, 3.5} × IMF slope {2.0, 2.3, 2.6}, per plan §1.4.
- **`P` is overloaded** — membership probability (WP2, per-star, [0,1]),
  p-values (gate statistics, low = bad), response probabilities (WP5 injections),
  and pulsar spin period (WP1/WP8). See [GLOSSARY.md](GLOSSARY.md).
