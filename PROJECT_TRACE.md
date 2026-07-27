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
| WP5 | IMF normalization + completeness | **BLOCKED** — mass-function residual gate | **no** |
| WP6 | Census closure + runaways | not started | — |
| WP7–WP10 | Ledger, cross-checks, verdict, manuscript | not started | — |

**Current blocker.** WP5's baseline residual gate fails on one bin. The cause is
diagnosed and scoped but **not yet fixed** — see §7 and
[WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md).
WP6 must not start until the re-run passes.

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

**Current state after repair_v1** — baseline (PARSEC, R_V = 3.1, α = 2.3):

| subgroup | χ² p | trend p | max abs residual | pass |
|---|---:|---:|---:|:--:|
| CygOB2-A | 0.090 | 0.872 | 2.63 | yes |
| CygOB2-B | 0.037 | 0.156 | 2.75 | yes |
| CygOB2-C | 0.040 | 0.623 | **3.26** | **no** |

Full branch grid: **31 / 54 passing**. Association mass 31,293 M☉ (within a
factor 2 of the 16,500 M☉ literature scale). 223,200 catalogue injections.
No subgroup reaches an absolute 95% completeness edge; bright-plateau
completeness is 0.833 / 0.854 / 0.831.

**⚠ The blocker is diagnosed but not fixed.** The repair report's verdict is
right and its diagnosis is wrong. See
[WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md)
for the full argument. In short: the forward response integrates the parent mass
distribution only to 8 M☉ while the mass estimator has 22–25% scatter, so the
~300 stars above 8 M☉ that scatter down into the top bin have no term in the
model. Fixing it is estimated to take the baseline to PASS and the grid to
43 / 54.

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

---

## 9. Open issues register

| # | Issue | Severity | Owner WP | Status |
|---|---|---|---|---|
| 1 | Forward response truncates parent mass distribution at 8 M☉; top bin under-predicted in all three subgroups | **blocking** | WP5 | diagnosed, fix not applied |
| 2 | [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) is the pre-repair run (0/54) and contradicts the repair parquet (31/54) | high | WP5 | stale, regenerate |
| 3 | Bright-mass completeness plateaus at ~0.83 (WP2 quality filter, not magnitude limit); WP6 step 2(a) assumes ~1.0 and will manufacture a false ~17% deficit → biases N_SN high | high | WP6 | flagged, not yet handled |
| 4 | No subgroup reaches an absolute 95% completeness edge; `CUTS_AND_THRESHOLDS.md` §7.1 procedure is unachievable on this field | medium | WP5 | fallback in use, needs formal supersession |
| 5 | Branches failing after the fix cluster at R_V = 3.5 and α = 2.6; plan §1.4 forbids silently dropping branches | medium | WP5 | policy decision needed |
| 6 | Subgroup-label (A vs B vs C) uncertainty never quantified — `membership_probability` measures cluster-vs-field only | medium | WP2 | open, non-blocking |
| 7 | A/B/C age ordering unresolved → star-formation-duration branch (0/1/2 Myr) stays mandatory | medium | WP4→WP7 | carried by design |
| 8 | Berlanas two-distance population not confirmed; single population at 1.62 kpc, depth 45 pc | low | WP2 | documented result |
| 9 | Upper-MS vs PMS indicator disagreement in subgroup B (N = 19 tail) | low | WP4 | documented |
| 10 | Wide-box definition inconsistency noted in `CUTS_AND_THRESHOLDS.md` §4.2 | low | WP1→WP6 | resolve before runaway traceback |

---

## 10. Next actions

1. Apply the WP5 parent-range fix (§7 and the correction report), re-run
   injections, re-fit and re-gate all 54 branches.
2. Regenerate [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) and refresh
   [AUDIT.txt](AUDIT.txt) via [audit.py](audit.py).
3. Record the branch-retention policy for issue #5 and the §7.1 supersession for
   issue #4 in [CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md).
4. Only then write a WP5 completion report that supersedes
   [wp5_completion_report.md](wp5_completion_report.md) and flips
   `downstream_wp6_authorized` in the gate record.
5. Before WP6 starts, resolve issue #3 — the closure test must divide by the
   injection response, not assume high-mass completeness ≈ 1.

---

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
