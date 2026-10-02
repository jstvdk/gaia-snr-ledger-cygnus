# Stage — after issue #19: settle C's age (issue #20), re-base and run WP13

*Written 2026-10-02 for an agent starting cold. Self-contained: read §0–§2, then
work §3 in order. Nothing here has been executed except what §2 says is done.*

## 0. What this project is, in five lines

Gaia DR3 census of Cygnus OB2 (d = 1.62 kpc), split into three kinematic
subgroups A/B/C, normalised through the IMF, converted into a branch-resolved
ledger of how many massive stars have died and when. Work packages WP0–WP12 are
complete and every step has a pre-registered gate; failed gates stay recorded as
failed. The 12 Aug 2026 novelty audit found most emphasised results already
established. The surviving claim — **resolving the subgroups materially changes
the death history because the young subgroup C sits below the first-death
boundary** — is a hypothesis under test, not a result. This stage decides it.

## 1. Read first, and the rules

Read in this order: `CLAUDE.md` · `PROJECT_TRACE.md` §1 (status board, with the
2026-10-01 repair_v8 notice) and §9 issues 19–20 · `reports/issue19_completion_report.md`
· `tasks/issue20_subgroup_c_age_brief.md` · `tasks/wp13_pooled_vs_resolved_ablation_brief.md`
· `reports/novelty_prior_art_and_incremental_value_audit_2026-08-12.md` §5–§9, §12.

**Never read `AUDIT.txt`** (16 MB, ~4.3 M tokens; regenerate with `audit.py`).

Rules not enforced by code:
1. Never type a number into `manuscript/main.tex`; add it to `scripts/wp10_numbers.py`
   and use the macro. `wp10_validate.py` V7/V2 fail the build otherwise.
2. Resolve manuscript inputs through `scripts/wp10_inputs.py`; WP12-style analyses
   through `wp12_common.frozen()` (hash-checked).
3. Nothing is overwritten, nothing is retuned. A cross-check that disagrees becomes
   an issue in PROJECT_TRACE §9, never a reason to move a number.
4. Predictions and thresholds are written, hashed and committed **before** the
   result is read. A pre-registration is written once and never regenerated.
5. The chain is **`repair_v8`** (selected in `scripts/chain.py`; override with
   `CYGOB2_CHAIN`). repair_v7 products stay on disk but are forbidden to quote.

Environment: the conda env `cygob2-gaia` is already first on PATH
(`python` = `/Users/vdk/miniforge3/envs/cygob2-gaia/bin/python`). Run scripts as
`PYTHONPATH=scripts python scripts/<name>.py` from the repo root. Commit only when
asked.

## 2. State at the start of this stage (repair_v8, adopted 2026-10-01)

Baseline branch = PARSEC, R_V 3.1, α 2.3, coeval, all-explode.

| quantity | value |
|---|---|
| N_death, association (A / B / C) | **8.36** (4.10 / 4.26 / **0.00**) |
| headline range, 36 branches (α ∈ {2.0, 2.3}) | 5.57–28.7; all 54: 1.94–28.74 |
| P(≥1 event) / P(last < 100 kyr) | 0.9996 / 0.547 |
| counts-based ages A / B / C (Myr) | 4.01 / 4.09 / **2.52** |
| closure ratio A / B / C, baseline | 1.006 / 1.068 / **1.36** |
| closing slopes A / B / C / association | 2.29 / 2.25 / 2.05 / 2.19 |
| WP5 | 40/54 cells pass; **gate G3 fails** its no-regression clause |
| WP9 verdict; α split | INCONCLUSIVE; 18/18 at α 2.0 vs 0/18 at α 2.3 |

Issue #19 did **not** touch C (its N_death, age and closure are identical on v7
and v8). Everything is committed in `cb16f7e` except one file; check
`git status`. Open items from #19 that this stage depends on: §3 step 3.

## 3. Steps (in order; each leaves a versioned artifact plus a provenance JSON)

### Step 1 — Issue #20 Phase A: is C young? (read-only, ~2–3 days; do this first)

**Why first.** If C is coeval with A, C contributes deaths, the baseline rises
above 8.36, the pooled-vs-resolved difference shrinks, and the WP13 "equivalent"
outcome becomes likely. WP13 may read M0/M1 only after this verdict, because M1's
C age depends on it.

**State of the evidence — do not treat it as settled.** One diagnostic
(`scripts/issue20_c_age_diagnostic.py`; tables `issue20_c_age_by_version.csv`,
`issue20_bright_members.csv`, `issue20_hrd_age_scan.csv`): a nearest-point
spectroscopic-HRD scan of C's 43 anchors prefers ~3.55 (PARSEC) / 4.01 (MIST) Myr
over 2.51 (Σχ² 144.0 vs 176.8 PARSEC). Its own caveats: no IMF weighting, no binary
model, reduced χ² ~3–4 (Δχ² overstates), luminous-biased sample, uncertain
supergiant T_eff. Against it: C's photometric age is younger than A's on every
grid cell (2.03–3.17 vs A 3.85–4.07, WP13 brief §2.2, repair_v7 grid, C unchanged
on v8); C holds an O3 If star; the young age appeared with the extinction repair,
which was validated against anchors. The live hypotheses are H1 young, H2 coeval,
**H3 two components**.

Tasks (full text in the issue-20 brief §4):
- **A1 Pre-register** `scripts/issue20_prereg.py` → `provenance/issue20_prereg.json`,
  with input SHA-256s, hypotheses, every test statistic, the decision rule
  (proposal in brief §5: H1 if C's 95 % upper bound < A's 95 % lower bound; H2 if
  the 68 % intervals overlap; H3 if a two-age mixture beats the best single age by
  ΔBIC > 10 and the ages differ by > 1 Myr, in **both** families; else
  inconclusive) and **what each outcome does to N_death** (state the expected
  direction, not a target). Add a *disclosed prior knowledge* section listing the
  diagnostic above and that the previous agent expected H2/H3 to be plausible.
  Mirror `scripts/issue19_prereg.py` (refuses to overwrite; hashes inputs;
  refuses if outputs already exist).
- **A2** Star-by-star audit of C's 43 spectroscopic stars and 15 brightest members
  → `tables/issue20_c_star_audit.csv` (spectral type/class, RUWE and binarity
  flags, membership probability, proper-motion Mahalanobis distance to the A/B/C
  centroids, position, spectroscopic vs photometric A_V, literature IDs).
- **A3** Proper spectroscopic-HRD age likelihood per subgroup, both families, three
  R_V: IMF-weighted, with an unresolved-binary component, σ_logTe 0.03–0.05 and a
  supergiant T_eff scale for class I/II. Then single-age vs two-age mixture for C.
- **A4** Photometric sensitivity refits of C (hybrid CMD+HRD likelihood; upper-MS
  fit excluding class I–II). Sensitivity only; never adopted from this step.
- **A5** Attribute the 3.98 → 2.51 jump: compare A_V and M_G,0 of C's brightest
  stars between pre-repair and repair_v1 extinction.
- **A6** Closure cross-check of C at ages 2.51, 3.16, 3.55, 3.98 Myr. Expected
  direction: older C → lower turnoff → fewer predicted massive stars → C's excess
  (slope 2.05) *grows*. Record as evidence either way.
- **A7** Literature check (Wright+2015, Berlanas+2019/2020, O3 If studies):
  cross-check only, never a calibration.

**Input hygiene.** The diagnostic reads the unversioned `data/processed/wp4_anchor_hrd.parquet`,
which `wp10_inputs` now **forbids**. Phase A must use the versioned
`data/processed/wp4_anchor_hrd_repair_v8.parquet` (legacy columns `logTe_spec`,
`teff_spec`, `spectral_type`, `MG0_obs`, `age_used_*`, `mass_*` equal the R_V 3.1
branch; per-branch columns carry `_rv3.0/_rv3.1/_rv3.5`). Do not use the
unversioned file. Note: the closure test (A6) reads isochrone-node responses built
at repair_v7 ages; check `scripts/wp6_closure_test.py` for how to vary C's age
without regenerating injections, and say so in the pre-registration if it can't.

**Acceptance.** `provenance/issue20_prereg.json` committed before A3–A6 run; outcome
JSON scores each hypothesis by the frozen rule; `reports/issue20_phase_a_report.md`;
issue #20 updated in PROJECT_TRACE §9 with the verdict; failed predictions stated
as failed.

**Phase B/C** (act on the verdict) are in the issue-20 brief §4. If H2 or H3 changes
C's age prior, that is a separate, separately pre-registered chain version
(`repair_v9`), run **after** this stage's WP13 decision, not folded into v8.

### Step 2 — Re-base the WP13 brief on repair_v8 (≈½ day; can run alongside Step 1)

The WP13 brief's §2 "disclosed prior knowledge" was computed on repair_v7. Add a
dated *re-basing note* (do not silently overwrite the old numbers; strike or
annotate them) and regenerate from repair_v8 products:
- §2.1 table: M1 = 8.43 / 0.552 → **8.36 / 0.547**; recompute the common-age rows
  (age forced to 4.0; k-weighted mean age; the age that reproduces the resolved
  count) from `tables/wp7_age_sensitivity_repair_v8.csv`.
- §2.2 C/A/B age ranges: recompute from `data/processed/wp5_imf_normalization_repair_v8.parquet`.
- Every other "8.43", "0.552", "repair_v7" in the brief and in
  `issue20_subgroup_c_age_brief.md` ("baseline rises above 8.43").
- Check the §5 thresholds still make sense against the new numbers (T2 needs
  |ΔN| ≥ 3 and ≥ 10 % of N_M1; nothing else is numeric in the baseline).
- Add the audit item: `slides/make_talk.py` still carries the note "This is a
  genuine discovery about the system".

### Step 3 — Decide the two issue-#19 judgement calls (user decisions; prepare options)

These must be settled **before** WP13 thresholds are frozen, because T6 requires
the result to hold on ≥ 80 % of the 36 headline branches, which are defined by
dropping α = 2.6.
1. **WP5 gate G3** fails on repair_v8 only through its "no regression in A or C"
   clause (A-PARSEC α = 2.0 at R_V 3.0/3.1, residual-trend p = 0.0048). Does a
   clause written for model changes apply to a data-defect fix?
2. **α headline set.** On repair_v8 the window χ² median is α 2.0 = 10.4, 2.3 =
   7.00, 2.6 = 10.06, so the set keeps a slope that fits worse than the one it
   drops. Registered test (2.3 vs 2.6) still passes. Options for the user: keep the
   36-branch set as registered and say so; add an α = 2.0 sensitivity; or revise
   the headline — which is a new pre-registered decision, not a retune.
3. **D1-P4** in `wp7_alpha_headline_adopt.py` fails literally (threshold hard-codes
   8.43); the chain reading governs and both are recorded. Confirm.
Prepare a one-page options table; do not decide for the user.

### Step 4 — WP13: pre-register, run M0, compare to M1 (after Steps 1–3)

Follow the WP13 brief §8: transcribe §5 into `scripts/wp13_prereg.py` with input
hashes through the repair_v8 chain, commit, **then** run M0 (pooled WP4 fit, pooled
injections, pooled joint fit, ledger), score T1–T6, run ledger-level injections
(T7, 200 realisations per design), write `reports/wp13_ablation.md` and the
PROJECT_TRACE row. M1's numbers are the stored repair_v8 ones, reproduced not
recomputed. Verdict rule and per-outcome manuscript framing are in the brief §5–§6.
Expected effort ≈ 3 weeks to a submittable manuscript. The brief's own expectation
(recorded so it can't be claimed as a prediction afterwards): the count test T2 is
marginal; the time-history tests T3/T4 are where a difference would live.

### Step 5 — Housekeeping (any time)

- `slides/cygob2_supernova_history_talk.pptx` slides **13** (closing α 2.25, 6.7 %),
  **16** (8.43, 0.9997, 0.552) and **17** (8.43, 5.63–28.7) still show repair_v7
  values. Check for the PowerPoint lock file `slides/~$cygob2_supernova_history_talk.pptx`
  first; do not write to a deck the user has open.
- fig04: user choice — draw the 61 unlabelled anchors, or adopt
  `scripts/draft_fig04_cmd_redesign.py`.
- F3(a): switch anchor masses from present-day `Mass` to initial `Mini` as its own
  pre-registered, separately scored change (affects 40–60 M☉ anchors by 10–25 %).
- `scripts/wp6_external_crosschecks.py` still reads the repair_v6 normalization
  (predates #19; not rerun to avoid mixing changes).
- Group-meeting brief and slides tables above the issue-19 bullet are still repair_v7.

## 4. Operational pitfalls learned in the previous stage

- **Iteration counts must match the published run** when re-running a script with
  a non-default setting (WP11 forecast: 500,000, script default is 200,000). A wrong
  count silently pinned the wrong file into a hash record; it was renamed (not
  deleted) to `provenance/wp12_revision_prereg_repair_v8_void_wp11_at_200k.json`.
  If you must void a record, rename it, reference it from its replacement.
- **Checks can be wrong, not just data.** Two integrity checks failed on bugs in the
  checks themselves (CSV round-trip at 1.8e-15; NaN comparison on an all-NaN
  column). Fix the check, record why, never loosen the science.
- **`scripts/chain.py`**: `C.tag(path)` inserts `_<chain>` before the extension for
  products repair_v7 wrote unsuffixed; `C.V[...]` gives the chain's upstream
  versions. New scripts that read WP6–WP12 products must go through it, or they
  will read repair_v7 files that `wp10_inputs` now forbids.
- **zsh**: a bare `echo =====` is a command-substitution error and aborts the rest
  of the command line. Use `printf '%s\n' '-----'`.
- **Auto-mode classifier** denied a Python one-liner that truncated a tracked
  script. Make such edits with the Edit/Write tools (git keeps the old version).
- Long jobs: run with `run_in_background` and wait on completion; do not poll.
- Convention for new analyses (mirror `scripts/issue19_prereg.py`,
  `issue19_integrity.py`, `issue19_score.py`): prereg script (hashes, predictions,
  integrity checks, adoption rule, disclosed prior knowledge) → run → integrity
  JSON → scorer that writes an outcome JSON and a before/after table.
- Never claim a number in prose that is not in a versioned artifact; quote
  `reports/issue19_completion_report.md` for the repair_v7 → v8 deltas.

## 5. Definition of done for this stage

1. `provenance/issue20_prereg.json` committed before any Phase A test ran; verdict
   recorded; PROJECT_TRACE §9 issue 20 updated.
2. WP13 brief re-based on repair_v8 with the old numbers still visible.
3. User decisions on §3 step 3 recorded in PROJECT_TRACE.
4. `scripts/wp13_prereg.py` hashed and committed before M0 was read; WP13 verdict
   and `reports/wp13_ablation.md` written; the manuscript framing chosen by the
   verdict, per the brief's §6 table.
