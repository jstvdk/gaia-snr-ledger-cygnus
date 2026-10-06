# Stage — WP4 age redesign (issue #21), the repair_v9 chain, then WP13

*Written 2026-10-06 for an agent starting cold. Self-contained: read §0–§2, then work
§3 → §6 in order. Owner decision: option (b) of the Phase A′ report §4, recorded in
[provenance/decisions_2026_10_06.json](../provenance/decisions_2026_10_06.json).
Nothing in this brief has been executed.*

---

## 0. What this project is, and why this stage exists

**The project.** A Gaia DR3 census of Cygnus OB2 (d = 1.62 kpc), split into three kinematic subgroups A/B/C. Each is normalised through the IMF and turned into a branch-resolved ledger of how many massive stars have died and when (`N_death`; all-explode unless stated). Every step has a pre-registered gate, and failed gates stay recorded as failed.

**The surviving novelty claim** was that resolving the subgroups changes the death history *because C is young* (photometric age 2.52 Myr, so 0 deaths).

**This stage exists because the ages are wrong in a known way.**
- **Issue #20** (Phase A, Phase A′): every age measurement that does not depend on a known defect puts C at about 3.2–4 Myr. These are the spectroscopic HRD and the bright upper main sequence.
- **Issue #21:** the WP4 photometric age fit itself is defective. Its extinction-error model is wrong, and the faint part of its fitting window is mis-modelled.

Every WP5–WP12 product inherits those ages, including `N_death = 8.36`. This stage:
1. redesigns the WP4 age fit, pre-registered and checked against spectroscopy;
2. rebuilds the chain on it as `repair_v9`;
3. only then runs WP13, the pooled-versus-resolved ablation that decides what the paper claims.

## 1. Read first, and the rules

**Read in this order:**
- `CLAUDE.md`
- `PROJECT_TRACE.md` §1, plus §9 issues 19, 20 and 21, and the Berlanas two-distance / foreground row
- `reports/issue20_phase_a_report.md`
- `reports/issue20b_phase_a_prime_report.md` (**the key document for this stage**)
- `scripts/issue21_error_model_diagnostic.py` with `tables/issue21_error_model_diagnostic.csv`
- `wp4_ages.md` §1, for the method only: its generated tables are pre-repair and stale (issue #19)
- `tasks/wp13_pooled_vs_resolved_ablation_brief.md`
- `reports/monte_carlo_methodology_explained.md`, for how ages enter the ledger

**Never read `AUDIT.txt`** (16 MB, ~4.3 M tokens; regenerate it with `audit.py`).

**Rules:**
1. **No bare numbers in the manuscript.** Never type a number into `manuscript/main.tex`. Add it to `scripts/wp10_numbers.py` and use the macro; `wp10_validate.py` V7/V2 fail the build otherwise.
2. **Resolve inputs properly.** Manuscript inputs go through `scripts/wp10_inputs.py`. WP12-style analyses go through `wp12_common.frozen()`, which is hash-checked.
3. **Nothing is overwritten or retuned.** A cross-check that disagrees becomes an issue in PROJECT_TRACE §9, never a reason to move a number.
4. **Predictions and thresholds are written, hashed and committed before the result is read.** A pre-registration is written once and never regenerated. If one must be voided, rename it (`…_void_<reason>.json`) and reference it from its replacement.
5. **Chain selection.** `repair_v8` stays the quoted chain until `repair_v9` passes its adoption rule (§4.4). It is selected in `scripts/chain.py`, overridable with `CYGOB2_CHAIN`.
6. **Commit only when the owner asks.** Before running any stage, show the owner the pre-registration for sign-off (§3.1, §4.1, §5.2).

**Environment.** Use the conda env `cygob2-gaia` (`python` is first on PATH). Run scripts as `PYTHONPATH=scripts python scripts/<name>.py` from the repo root.

## 2. State at the start of this stage

**Quoted chain: `repair_v8`.** Baseline = PARSEC, R_V 3.1, α 2.3, coeval, all-explode.

| quantity | value |
|---|---|
| N_death (A / B / C) | **8.36** (4.10 / 4.26 / 0.00) |
| P(≥1) / P(last < 100 kyr) | 0.9996 / 0.547 |
| WP4 photometric upper-MS ages (repair_v5) | A 3.98 · B 3.55 · C 2.51 Myr (PARSEC, R_V 3.1) |
| WP5 counts-based ages (used by the ledger) | A 4.01 · B 4.09 · C 2.52 |
| headline rule (decision of 2026-10-02) | α = 2.3 only (18 branches); α 2.0 and 2.6 are sensitivity branches |
| closing slopes | A 2.29 · B 2.25 · C 2.05 · association 2.19 |

**What issue #21 established** (post-hoc diagnostic, reproduces every stored WP4 MAP exactly). PARSEC / MIST, R_V 3.1, f_bin 0.4:

| | WP4 model, full window | A_V error correlated, full window | bright window M_G0 ≤ −1, either error model | spectroscopic HRD, Phase A′ test c |
|---|---|---|---|---|
| A | 3.98 / 4.01 | 5.62 / 6.37 | 2.82 / 4.01 (broad) | 3.16 / 3.57 |
| B | 3.55 / 4.01 | 5.62 / 8.03 | 3.98 / 2.83 (correlated: 2.82 / 4.01) | invalid (4 stars) |
| C | **2.51** / 3.18 | 5.01 / 4.50 | **3.98** / 3.18 | 3.55 / 4.01 |

**Reading of this table:**
- `wp4_common.branch_photometry` (line ~134) adds `(k_G·σ_AV)²` to σ²(M_G0) and `((k_BP−k_RP)·σ_AV)²` to σ²(colour) as *independent* terms. A real A_V error moves a star along the reddening vector, so the two are fully correlated.
- Propagated correctly, the full window gives 4.5–8 Myr. That is impossible: A contains O5–O7 supergiants, which need a turnoff of at least ~40 M☉, i.e. ≲ 4.5 Myr. So **the faint part of the window (−1 < M_G0 ≤ 1.5) is mis-modelled whatever the error model.**
- The bright window is nearly insensitive to the error model.

**Resources available:**
- **The full per-star A_V posterior:** `data/processed/wp3_extinction_posterior_repair_v5.npz`, with `source_id` (1392), `rv` (3), `av_grid` (301) and `probability` (1392 × 3 × 301). Its row order must equal the repair_v5 extinction table; `wp4_mass_posteriors_repair.py` already asserts this.
- **Bright-window counts** (M_G0 ≤ −1, labelled members):

| R_V | A | B | C |
|---|---:|---:|---:|
| 3.0 | 94 | 114 | 157 |
| 3.1 | 99 | 120 | 162 |
| 3.5 | 114 | 147 | 193 |

  The full upper-MS window holds about 400 stars per subgroup.
- **Spectroscopic-HRD posteriors** (outlier-robust, Phase A′ test c): `tables/issue20b_age_tests.csv`. Spectroscopic stars per subgroup: A 59, C 43, B 4.

---

## 3. Stage 1 — WP4 age redesign (read-only with respect to the chain)

### 3.1 Pre-register first

Write `scripts/wp4v9_prereg.py` → `provenance/wp4v9_age_prereg.json`. Mirror `scripts/issue20b_prereg.py`: it must refuse to overwrite, hash every input, and refuse to run if any output already exists. It must contain:
- the method (§3.2), with every constant;
- the injection validation (§3.3);
- the acceptance rule and fallback (§3.4);
- the predictions (§3.5);
- a **disclosed prior knowledge** section, listing:
  - the issue #21 table above;
  - Phase A and Phase A′ results;
  - that the previous agent expects C at about 3.2–4 Myr;
  - that A's bright window was broad in the diagnostic, so acceptance for A-PARSEC may fail (§3.4).

**Show the pre-registration to the owner, get sign-off, commit it, then run.**

### 3.2 Method (what changes in WP4, and what does not)

**M1 — Extinction propagation by marginalisation (primary).**
- **Before:** each star contributed a Gaussian error in (colour, M_G0).
- **Now:** integrate each star's likelihood over its own A_V posterior:
$$
\mathcal{L}_i(t) = \sum_{j} p_i(A_{V,j})\;\mathcal{L}^{\rm phot}_i\!\left(\mathbf{x}_i^{\rm obs} - A_{V,j}\,\mathbf{r}\;\middle|\;t\right),
\qquad \mathbf{r} = \bigl(k_{BP}-k_{RP},\; k_G\bigr),
$$
  where $\mathbf{x}^{\rm obs}$ = (observed colour, observed $G - \mu$), and $\mathcal{L}^{\rm phot}$ is WP4's particle-mixture likelihood using photometric errors plus the 0.03 mag floor only.
- Take the band coefficients from the same `band_coefficients(rv)` WP3 used.
- Anchors keep their spectroscopic A_V; their posterior is narrow, which is handled automatically.
- **Secondary (sensitivity only):** the Gaussian correlated-error model of `issue21_error_model_diagnostic.py`.

**M2 — Indicator window.**
- **Primary:** the bright upper main sequence, **M_G0 ≤ −1.0**. This is where the turnoff information is, and where the issue #21 diagnostic showed the result does not depend on the error model.
- **Sensitivity edges:** −0.5 and −1.5.
- **Selection consistency (pre-register the choice).** M_G0 depends on the marginalised A_V, so define window membership on the star's posterior-median M_G0, and apply the same cut to the model particles. Alternatively, model the window as a selection function inside the likelihood. State which; do not switch after seeing results.

**M3 — Unchanged from WP4:**
- the forward model: IMF-weighted particles, the explicit unresolved-binary component, f_bin ∈ {0.3, 0.4, 0.5};
- the native isochrone age grid (1–10 Myr, 0.05 dex), PARSEC and MIST, R_V ∈ {3.0, 3.1, 3.5};
- membership-probability weights;
- distance fixed at 1.6245 kpc, with ±σ_μ refits;
- the measurability gate (≥ 15 stars, no MAP within one grid step of the grid edges).

**M4 — Faint-end study (diagnostic, never adopted from this step).** Fit the faint window (−1 < M_G0 ≤ 1.5) alone with M1, then test four candidate causes one at a time:
- (i) PMS–MS transition modelling: compare against the issue #1c Henyey-fold finding;
- (ii) the binary fraction: f_bin scan beyond 0.5;
- (iii) contamination: restrict to membership probability P > 0.9;
- (iv) A_V outliers: a heavy-tailed or outlier mixture term.

Report which cause, if any, reconciles the faint window with the bright window. This result feeds future work; it does not gate this stage.

### 3.3 Injection validation (do before reading real-data posteriors)

For each subgroup, build synthetic populations at 2.5, 3.2, 4.0 and 5.0 Myr:
- matching the real star counts in each window and the real per-star A_V posterior widths;
- drawing true A_V from each real star's posterior, so the errors are correlated along r.

Fit them with M1 and with the old WP4 independent-error model.

**Required:** M1 recovers the true age with |bias| < 0.2 Myr and 68 % coverage of 0.55–0.80, in both families. Record how biased the old model is; that is the quantitative statement of issue #21. Budget about 1–2 h of compute. Mirror `issue20b` J1, which already showed that recovery machinery works.

### 3.4 Acceptance rule and fallback (freeze in §3.1)

All checks are evaluated per subgroup × family at the decision cell (R_V 3.1, f_bin 0.4). Robustness across the other cells is a qualifier: robust or fragile.

- **GA1 (validation):** §3.3 passes.
- **GA2 (agreement with spectroscopy):** for A and C, the M1 bright-window 68 % interval overlaps the Phase A′ test-c spectroscopic 68 % interval.
- **GA3 (physical plausibility):** less than 5 % of the age posterior lies at ages whose turnoff is below the initial mass of the subgroup's most massive confirmed main-sequence or giant spectroscopic member. Use `wp6_mass_extension_decision.turnoff_mass`.
- **B:** there is no spectroscopic check (4 stars). B is accepted on GA1 and GA3 only, and is flagged "photometry-only" wherever quoted.
- **Fallback, per subgroup × family where GA2 fails.**
  - **Pre-registered:** the adopted age is the **joint** likelihood, bright-window photometry (M1) × spectroscopic HRD (the test-c model); the two data sets are disjoint where anchors are removed from the photometric term.
  - If the joint posterior is itself bimodal, or fails GA3, carry the two ages as **separate branches** through repair_v9. In that case no claim may rest on the difference between subgroup ages.
- **Known risk, disclosed before running:** in the diagnostic, A-PARSEC's bright window was broad (68 % 1.1–2.8 Myr) against its spectroscopic 2.8–3.5. GA2 may fail for A in PARSEC.

### 3.5 Predictions (directional; write before running)

- **W1:** C's adopted age is ≥ 3.2 Myr in both families.
- **W2:** the |A − C| age difference is < 1 Myr in both families.
- **W3:** the old independent-error model is biased young in the faint window on injections.
- **W4:** the faint-window study (M4) does not reconcile the faint window under any single tested cause.

These are expectations to be scored, not targets.

### 3.6 Outputs

- `data/processed/wp4_age_posteriors_repair_v9.parquet`: the same schema as repair_v5, plus `window`, `error_model`, `adopted`, `fallback` columns.
- `tables/wp4v9_injection_validation.csv`
- `tables/wp4v9_faint_end_study.csv`
- `provenance/wp4v9_age_outcome.json`, scoring GA1–GA3 and W1–W4
- `reports/wp4v9_age_redesign.md`
- An update of issue #21 in §9.

---

## 4. Stage 2 — the repair_v9 chain (ages change; nothing else does)

### 4.1 Pre-register

Write `scripts/repair_v9_prereg.py` → `provenance/repair_v9_prereg.json` (mirror `scripts/issue19_prereg.py`), before any chain step. Contents:
- the input hashes, including the adopted WP4 posteriors from §3;
- the integrity checks (§4.3);
- the chain predictions:
  - **V1:** baseline N_death rises if C's age rises above 3.0 Myr. The estimate is +1.1 at 3.2 Myr, +2.8 at 3.6 and +4.6 at 4.0, holding A and B fixed. A and B will also move; state the expected direction per subgroup from §3's result.
  - **V2:** the closure ratios move in the direction implied by each turnoff change.
- the adoption rule (§4.4).

Get owner sign-off before running.

### 4.2 Steps, in order (model the runner on `scripts/run_repair_v8_chain.sh`)

`export CYGOB2_CHAIN=repair_v9 WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging`.

0. **`scripts/chain.py`:** add `"repair_v9"` with `wp3_extinction: repair_v5` (unchanged), `wp4_ages: repair_v9`, `wp4_masses: repair_v9`, `wp5: repair_v9`, `responses: repair_v9`. Leave `ADOPTED` unchanged until §4.4.
1. **Anchor masses at the v9 ages:**
   `wp4_anchors_hrd.py --age-version repair_v9 --extinction-version repair_v5 --output-version repair_v9`.
   Keep the present-day `Mass` as v8 did. F3(a), switching to `Mini`, stays a separate later change (§6.4).
2. **Photometric mass posteriors:**
   `wp4_mass_posteriors_repair.py --anchor-version repair_v9 --output-version repair_v9`, reading the v9 age posterior. Check that the script reads `wp4_age_posteriors_{version}` from an argument rather than a constant, and add the argument if needed. Also fix the stale `mass_method` label (`photometric_posterior_repair_v1`).
3. **WP5 injections at the new truth-age nodes, then the joint age–k fit.**
   - The node responses depend on the ages, so **they must be regenerated; they cannot be reused** as in v8. Use `wp5_injections_agenodes.py` / `wp5_joint_age_fit.py`; read `run_repair_v6_chain.sh` for how nodes were built. About 1–2 h.
   - Then `wp5_fit_imf_joint.py --mass-version repair_v9 --age-version repair_v9 --response-version repair_v9 --wp5-version repair_v9 --compare-version repair_v8`.
   - Re-score the WP5 gates; G3 stays recorded as failed on v8.
4. **WP6:**
   - mass-extension injections at the new node ages: `wp6_mass_extension_decision.py` plan, then `wp6_multiplicity_injections.py` / the injection driver;
   - the closure test (`wp6_closure_test.py`, 4.0 M☉ floor);
   - the massive census and orphan anchors;
   - `wp12_closure_slopes.py`.
5. **WP7:**
   - `wp7_ledger.py` at **2,000,000** iterations; the iteration count must match the published one;
   - age scan, BH scan, convergence scan, `wp7_binary_bound.py`;
   - the headline set under decision 2 (α = 2.3).
6. **WP8–WP12:**
   - WP9 verdict, plus the WP12 gate map, combination gate, branch-gate table, mixed-slope ledger and scenario score;
   - write a **new** `provenance/wp12_revision_prereg_repair_v9.json` hash record; never edit the v8 one;
   - WP11 runs at **500,000** iterations (its script default is 200,000).
7. **Manuscript plumbing:**
   - `wp10_inputs.py`: bump versions; forbid the superseded v8 WP4/WP5 age products from being quoted;
   - `wp10_numbers.py` → `numbers.tex`, including the age envelope, which is now the bright-window envelope;
   - `wp12_tables.py`, `wp12_figures.py`, `wp10_validate.py`.

### 4.3 Integrity checks (pre-registered)

- **I1, replay:** feeding the repair_v5 WP4 posteriors through the v9 code path reproduces the repair_v8 WP5 normalisation and the WP7 baseline to within Monte Carlo noise. This proves the code changes are inert apart from the ages.
- **I2:** the anchor masses are per-branch, and their ages equal the v9 posteriors.
- **I3:** every v9 input hash matches the pre-registration.
- **I4:** no v9 script reads an unversioned or v8 age product: grep, and check against `wp10_inputs.FORBIDDEN`.

### 4.4 Adoption rule

repair_v9 is adopted (`ADOPTED = "repair_v9"` in `chain.py`) if I1–I4 pass and §3's acceptance or fallback is documented. Predictions V1–V2 are scored but do not gate adoption.

Then:
- write `reports/repair_v9_completion_report.md` with a before/after table for every headline number;
- close issue #21, and update issue #20 with C's adopted age;
- update the `CLAUDE.md` headline numbers.

### 4.5 Sensitivity items run on repair_v9 (each reported separately; never the headline)

- **S1 — IMF upper limit m_max ∈ {100, 120, 150} M☉,** baseline branch.
  - `IMF_UPPER_LIMIT = 120` is a convention, not fitted (`wp6_mass_extension_decision.py:54`).
  - `wp7_ledger.TurnoffRelation` only inverts up to 1.02 × 120. Extending it to 150 needs the turnoff grid to bracket 150 (PARSEC: about 2.8 Myr), and must respect issue #14: above about 120 M☉ the tables' maximum is a table ceiling at young ages.
  - Analytic expectation on v8 inputs: 100 → 7.0, 120 → 8.4, 150 → 9.7 deaths.
- **S2 — Subgroup-assignment uncertainty.**
  - The labels are hard, but a refit of the WP2 Gaussian mixture gives a median maximum responsibility of 0.90. Only 51 % of stars exceed 0.9; B is clean, and A and C mix by about 15 % each way.
  - Recompute k and N_death with responsibilities as weights, or resample labels per Monte Carlo iteration.
  - The responsibilities are not stored. Pre-register the refit (k = 3, full covariance, StandardScaler on l, b, μα\*, μδ, 50 seeds, consensus mapping to the stored labels).
- **S3 — Parallax-blind membership check** (open register item from the Berlanas row): foreground groups at about 1.3 kpc are separable by proper motion. This is optional in this stage, and recorded if not done.

---

## 5. Stage 3 — WP13 (pooled vs resolved), on repair_v9

1. **Re-base the WP13 brief on repair_v9.** Use `scripts/wp13_brief_rebase.py` as the pattern for the v8 re-basing. Add a dated note, and strike or annotate the old numbers rather than overwriting them. Recompute:
   - §2.1: M1, plus the common-age rows (forced 4.0, k-weighted mean age, and the age that reproduces the resolved count) from `wp7_age_sensitivity_repair_v9.csv`;
   - §2.2: the age ranges.
2. **Freeze before M0 is run.** Transcribe §5's pass criteria into `scripts/wp13_prereg.py` with input hashes through the repair_v9 chain, and commit it.
   - **T6 is scored on the 18 headline branches** (α = 2.3; decision 2).
   - **M0, the pooled model, must use the same redesigned age method** (M1 marginalisation, bright window), applied to the pooled sample. Otherwise the ablation compares two age methods instead of pooled versus resolved.
3. **Run M0:** the pooled WP4 fit, pooled injections, the pooled joint fit and the ledger. Score T1–T6, then run the ledger-level injections (T7, 200 realisations per design).
4. **Write up:** `reports/wp13_ablation.md` and the PROJECT_TRACE row. Apply the brief's §6 framing for whichever outcome occurs.
5. **Disclosed expectation (record it so it can't be claimed as a prediction):** if repair_v9 puts C at about A's age, the "equivalent" outcome is likely. The paper then becomes a careful Gaia census and death ledger, plus the methodological finding: photometric ages of hot massive stars need correct extinction-error propagation and spectroscopic anchoring.

## 6. Stage 4 — after WP13

1. **Manuscript.** Rewrite the abstract, results and conclusions around the WP13 outcome. Also:
   - remove every claim built on the PMS indicator or on a young C;
   - replace "supernova history" with conditional wording;
   - every number goes through `numbers.tex`.
2. **Figure 4.** Adopt the CMD-only redesign (`scripts/draft_fig04_cmd_redesign.py --cmd-only`) into `wp12_figures.py`, which needs owner approval. It must draw the 61 unlabelled anchors.
3. **Talk materials.**
   - The deck (https://claude.ai/artifact/3EWLoXBDL97KkSKGEs5xHM, private) and `reports/deck_slide_guide.md`.
   - `reports/group_meeting_brief_2026-09-30.md`.
   - `slides/cygob2_supernova_history_talk.pptx`: still repair_v7 values on slides 13, 16 and 17. Check for the lock file `slides/~$…pptx` first, and do not write to a deck the owner has open. `slides/make_talk.py` still carries the note "This is a genuine discovery about the system"; remove it.
4. **Separate, pre-registered changes left for later:**
   - F3(a): anchor masses from `Mass` to `Mini` (affects 40–60 M☉ anchors by 10–25 %);
   - `wp6_external_crosschecks.py` still reads the repair_v6 normalisation.

## 7. Done criteria

| stage | done when |
|---|---|
| 1 | the pre-registration is committed before any posterior is read; GA1–GA3 and W1–W4 are scored in `wp4v9_age_outcome.json`; the report is written; issue #21 is updated |
| 2 | I1–I4 pass; repair_v9 is adopted or explicitly not adopted, with reasons; the completion report has before/after tables; S1 and S2 are reported; the issues are updated; CLAUDE.md is updated |
| 3 | the WP13 pre-registration is committed before M0; T1–T7 are scored; the report and register row are written; the manuscript framing is chosen by the frozen rule |
| 4 | the manuscript builds (`wp10_validate.py` clean); figures, deck, brief and pptx are consistent with the adopted chain |

## 8. Pitfalls learned in earlier stages

- **Iteration counts must match the published run** (ledger 2,000,000; WP11 500,000). A wrong count once pinned the wrong file into a hash record.
- **Checks can be wrong, not just data.** Two integrity checks have failed on bugs in the checks themselves. Fix the check, record why, and never loosen the science.
- **`chain.tag(path)`** inserts `_<chain>` before the extension for products repair_v7 wrote unsuffixed. New scripts that read WP6–WP12 products must go through it.
- **Generated report tables can be stale.** `wp4_report.py` was hard-wired to the unversioned posteriors (issue #19). Any report generator touched here must take a version argument.
- **Window selection and A_V marginalisation interact** (§3.2 M2). Decide the selection rule once, in the pre-registration.
- **Tiny posterior intervals (±0.05 Myr) usually signal a mis-specified model,** not precision (Phase A′ §2).
- **Effort:** Stage 1 about 2–4 days; Stage 2 about 1 day of compute plus checks; Stage 3 about 1–3 weeks to a submittable manuscript (WP13 brief §8).
