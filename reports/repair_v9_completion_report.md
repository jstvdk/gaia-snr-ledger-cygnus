# repair_v9: completion report

Date: 2026-10-07 (the chain ran 18:57–20:37 CEST on the owner's AlmaLinux machine). Status: **ADOPTED.** `scripts/chain.py` now has `ADOPTED = "repair_v9"`.

Pre-registration: `provenance/repair_v9_prereg.json`, committed in `d947f31` before any chain step. Integrity record: `provenance/repair_v9_integrity.json`. Outcome: `provenance/repair_v9_outcome.json`. Before/after table: `tables/repair_v9_before_after.csv`.

## 0. Handoff: start here (for the next agent)

This report replaces [tasks/HANDOFF_repair_v9.md](../tasks/HANDOFF_repair_v9.md), whose work is all done. Read this section, then `CLAUDE.md`. "The owner" is the project's human lead; use they/them.

### 0.1 Check before anything else

1. **Code:** `git log --oneline -1` should show the handoff commit, which comes after `146698a`.
2. **Data:** `data/` is gitignored, and repair_v9 added 1,300 files (3.6 GB), so `data/` is now about 13 GB. They are listed with SHA-256 in `provenance/repair_v9_data_manifest.sha256`.
   - Verify from the repository root:
     ```sh
     (cd data && sha256sum -c ../provenance/repair_v9_data_manifest.sha256 --quiet) && echo "v9 data present"
     ```
     On macOS use `shasum -a 256 -c`.
   - If the files are missing, **do not regenerate them.** Ask the owner to copy `data/` from the AlmaLinux machine (`/home/vvoitsek/science/gaia-snr-ledger-cygnus/data`), all of it or just the 1,300 manifest files.
   - Then run `PYTHONPATH=scripts python scripts/wp10_inputs.py`; it must report `audit: PASS`.
3. **Environment:** follow `INSTALL.md`. On Linux the env also needs the pip packages from the pinned file, including `gaiadr3-zeropoint`, or `wp5_common` fails to import.
4. **Prompt log:** the owner wants every prompt they give recorded in `AI_PROMPTS.md` with the date and time. A `UserPromptSubmit` hook in `.claude/settings.json` does this automatically. If a prompt is missing from the log, add it by hand, marked "(backfilled)".

### 0.2 Where things stand (2026-10-08)

- **Quoted chain:** `repair_v9`, adopted 2026-10-07. The ages are spectroscopic; C is coeval with A.
  - Baseline N_death **6.66** (A 2.09, B 2.85, C 1.72).
  - P(last < 100 kyr) 0.697.
  - 36-branch range 6.29–34.5.
- **Manuscript:** `numbers.tex`, the tables and figures are regenerated, and `wp10_validate` passes. The prose in `main.tex` is **not** updated.
- **Done:** S1 (IMF ceiling) and S2 (subgroup labels); see [repair_v9_sensitivity_s1_s2.md](repair_v9_sensitivity_s1_s2.md). WP13 is re-based, with its pre-registration **drafted, not signed**.
- **Not done, recorded:** S3 (the parallax-blind membership check).

### 0.3 Next steps, in order. Each needs the owner where marked.

1. **Issue #22 (owner decision).** The WP5 low-mass counts disagree with the spectroscopic ages: the gate passes in 28 of 54 cells, C in 6 of 18 (§2, §5; PROJECT_TRACE §9 #22). Two options:
   - (a) accept it as a stated caveat;
   - (b) investigate it first. Candidate causes from #22: the faint-window mis-modelling of #21; low-mass membership contamination (Stage 1 M4 (iii)); a real age spread in C.
   
   Any investigation needs its own pre-registration.
2. **WP13 sign-off (owner).** [tasks/wp13_prereg_proposal.md](../tasks/wp13_prereg_proposal.md) lists decisions D1–D6. The recommended age method for M0 is a pooled test-c spectroscopic fit. After sign-off:
   - transcribe it into `scripts/wp13_prereg.py`, hash it and commit it **before** computing anything on the M0 side, including T1, which is computable in seconds and deliberately not computed;
   - run M0, using `scripts/repair_v9_injections.py` with a pooled-label age table;
   - score T1–T7 and write `reports/wp13_ablation.md`.
   
   The disclosed expectation is the "equivalent" outcome.
3. **Stage 4, the manuscript** (stage brief §6), after WP13 decides the framing:
   - rewrite the abstract, results and conclusions; remove every young-C, PMS and below-the-first-death-boundary claim;
   - macro names such as `\ageumsA` and `\ageumsC` now hold spectroscopic ages; rename them in `wp10_numbers.py`;
   - two macros render as a dash on repair_v9 (`\strictTwolo` / `\strictTwohi` and `\mixedShiftDead`). The sentences that use them must change;
   - add S1's ceiling systematic (−32 % / +30 %) and S2's ±0.5 per-subgroup label uncertainty to the sensitivity section;
   - Figure 4: adopt the CMD-only redesign (owner approval);
   - compile with Tectonic once `manuscript/aa.cls` and `aa.bst` are fetched;
   - regenerate `AUDIT.txt` with `audit.py`. **Never open `AUDIT.txt`.**
4. **Talk materials (stage brief §6.3):**
   - the deck and `reports/deck_slide_guide.md`;
   - `reports/group_meeting_brief_2026-09-30.md`;
   - `slides/cygob2_supernova_history_talk.pptx`, which still has repair_v7 values. Check for an open lock file first.
   - `slides/make_talk.py` still says "This is a genuine discovery about the system"; remove that.
5. **Separate changes, each pre-registered later:**
   - F3(a), anchor masses from `Mass` to `Mini`;
   - `wp6_external_crosschecks.py` still reads the repair_v6 normalization;
   - optionally S3.

### 0.4 Rules learned in this stage (in addition to CLAUDE.md)

- **Order of work:** pre-register, commit, then run.
  - Commit only when the owner asks. The overnight delegation (`provenance/decisions_2026_10_07_overnight.json`) covered 2026-10-07 only.
  - A defect in a check is fixed and recorded in `provenance/repair_v9_deviations.json`; a threshold is never loosened.
- **Long runs:** use `tmux` on Linux (there is no `screen` or `caffeinate` there), with `OMP_NUM_THREADS=1`.
  - `scripts/repair_v9_injections.py` parallelises injections safely: every node starts a fresh seed. With 28 workers, 162 nodes take about 20 min.
- **Chain registry:** `scripts/chain.py` registers `repair_v9`, `repair_v9_replay` and `repair_v9_scan00`–`10`. Use `CYGOB2_CHAIN` to address a non-adopted chain.
- **Generated figures can carry hard-coded prose.** Check each claim against the active chain after any rerun. Three such claims were stale this time (§7).

## Morning summary for the owner

1. **repair_v9 is adopted.** All the integrity checks pass, I1–I5.
   - Two parts of the code-replay check reproduce repair_v8 exactly (I1a, I1e).
   - Three agree with it to the last floating-point bits across platforms (I1b, I1c, I1d).
   - One check (I4) first failed on a **bug in the check itself**. It was fixed and recorded as deviation D1; no threshold changed (§3).
2. **New baseline: N_death = 6.66** (A 2.09 · B 2.85 · C 1.72), down from 8.36 (A 4.10 · B 4.26 · C 0.00). This is still the all-explode count, PARSEC, R_V 3.1, α 2.3, coeval.
   - **C is not below the first-death boundary.** It has 1.72 deaths and P(≥ 1) = 0.79. The paper's surviving novelty claim is therefore false on the adopted chain.
3. **P(last death < 100 kyr) rises from 0.547 to 0.697**, because C now contributes recent deaths.
   - The 36-branch headline range is 6.29–34.5 (was 5.57–28.7).
   - All 54 branches span 2.02–34.5.
4. **Most important caveat: the WP5 mass-function gate is much worse.**
   - It passes in 28 of 54 cells (was 40), 18 of 36 headline cells (was 25), and 1 of 12 headline combinations (was 5).
   - Almost all the loss is **C, which drops from 15 of 18 passing cells to 6**.
   - The age scan shows why: C's low-mass counts pass the gate only at ages ≤ 3.16 Myr (PARSEC) or ≤ 3.57 Myr (MIST). They fail at C's spectroscopic age (3.55 / 4.01 Myr).
   - The WP5 counts also pull every subgroup off its prior: C's WP5 age falls to 3.38 Myr (prior 3.55), while A and B rail at the **top** of their age nodes (A 15 of 18 cells, B 18 of 18).
   - So the spectroscopic O-star ages and the Gaia low-mass counts disagree, in C most of all. This is the same tension Stage 1 found and it is **unresolved**. It does not gate adoption (the pre-registered rule is I1–I5), but every number below inherits it.
5. **Predictions:** V1a, V1b, V1c, V2, V3 and V5 came true. **V4 failed** (§4).
6. **The age scan** (§5): at equal ages the three subgroups have almost the same death curve, within 10 %. On the adopted ages any difference between subgroups therefore comes from age alone.
7. **Tonight's follow-on work, all done** (S1 and S2 pre-registered in `3c571f6` while this chain ran, before any v9 result; report in [repair_v9_sensitivity_s1_s2.md](repair_v9_sensitivity_s1_s2.md)):
   - **S1, IMF ceiling:** 100 / 120 / 150 M☉ gives 4.55 / 6.66 / 8.66 deaths (−32 % / +30 %). The ceiling is now one of the largest systematics, because the turnoffs sit close to it.
   - **S2, subgroup labels:** the association count moves +1.2 %, but A and C trade about 0.4–0.5 deaths each.
   - All checks pass and all predictions are true for both.
   - **WP13:** the brief is re-based. A single common age reproduces the resolved model (6.64 against 6.66 deaths; P(last < 100 kyr) 0.697 against 0.697), so the outcome to expect is "equivalent". The pre-registration is drafted with six decisions for you: [tasks/wp13_prereg_proposal.md](../tasks/wp13_prereg_proposal.md). **M0 was not run.**
8. **Needs you:**
   - the WP5 gate degradation for C (item 4) — accept it as a caveat, or investigate;
   - signing off the WP13 pre-registration draft;
   - the Stage 4 prose in `main.tex`, which still argues for a young C;
   - pushing the overnight commits (nothing was pushed).

## 1. What changed and what did not

**Changed: the WP4 age table only.** It is now `wp4_age_posteriors_repair_v9_headline.parquet`:
- **A and C** carry their Phase A′ test-c spectroscopic-HRD posteriors. These are the frozen Stage 1 fallback.
- **B** borrows the A + C combined posterior: the product of the two likelihood curves. It is flagged "not measured for B".

| R_V 3.1 | PARSEC MAP (68 %) | MIST MAP (68 %) | source |
|---|---|---|---|
| A | 3.16 (2.80–3.50) | 3.57 (3.31–3.70) | test c, 59 stars |
| B | 3.55 (3.16–3.70) | 3.57 (3.52–3.88) | A + C product, **not measured for B** |
| C | 3.55 (3.32–3.97) | 4.01 (3.76–4.27) | test c, 43 stars |

**Regenerated downstream at the new ages:**
- anchor masses, read per R_V branch at the new MAPs;
- photometric mass posteriors;
- 162 WP5 truth-age node injections and 162 WP6 mass-extension injections, with the repair_v7 truth model (interpolated nodes, mass-dependent f_bin);
- the WP5 joint age–k fit;
- WP6 to WP12, with WP7 at 2,000,000 iterations and WP11 at 500,000.

**Unchanged:** WP1–WP3 (membership, labels, the repair_v5 extinction), the isochrones, the IMF grid, the explodability rules, and the turnoff and lifetime relation.

**Code changes:** all are inert for earlier versions, and I1 tests this.
- `chain.py`: new chain entries.
- `wp5_joint_age_fit.py`: the interpolated-version set.
- `wp4_mass_posteriors_repair.py`: `--age-version`.
- `wp4_anchors_hrd.py`: `--preregistration`.
- `wp12_prereg.py`: the chain record re-points the age table.
- `repair_v9_injections.py`: a new parallel driver around the unchanged `inject_curve`.
- The "stale `mass_method` label" from the proposal had already been fixed in repair_v8. Nothing changed for it.

## 2. Before / after (baseline PARSEC, R_V 3.1, α 2.3, coeval, all-explode)

| quantity | repair_v8 | **repair_v9** |
|---|---:|---:|
| N_death, association | 8.36 | **6.66** |
| N_death A / B / C | 4.10 / 4.26 / 0.00 | **2.09 / 2.85 / 1.72** |
| P(≥ 1 death) A / B / C | 0.98 / 0.98 / 0.00 | 0.81 / 0.93 / **0.79** |
| P(≥ 1), association | 0.9996 | 0.9974 |
| P(last < 100 kyr), association | 0.547 | **0.697** |
| median t_last (kyr) | 88 | 57 |
| headline range, 36 branches | 5.57–28.7 | **6.29–34.5** |
| α = 2.3 branches only (18) | 5.57–11.16 | 6.30–13.41 |
| all 54 branches | 1.94–28.74 (median 8.72) | 2.02–34.53 (median 9.57) |
| α = 2.6 branches | 1.94–4.24 | 2.02–5.05 |
| WP5 truth age A / B / C (Myr, posterior mean) | 4.01 / 4.09 / 2.52 | 3.48 / 3.70 / 3.38 |
| WP5 k_median A / B / C | 1,692 / 1,647 / 1,894 | 1,769 / 1,693 / 1,757 |
| turnoff at the baseline age A / B / C (M☉) | 58 / 56 / 279 | 77 / 68 / 82 |
| minimum progenitor mass (M☉) | 33.9 | 32.8 |
| closure ratio, baseline cell A / B / C | 1.006 / 1.068 / 1.360 | 0.924 / 1.036 / **1.498** |
| closure grid median (median of 18 ratios) | 1.073 | 1.049 |
| closing slope A / B / C / association | 2.29 / 2.25 / 2.05 / 2.19 | 2.33 / 2.27 / 2.00 / 2.19 |
| living ledger above 8 M☉ | 388.6 | 386.5 |
| WP5 residual gate, cells passing (of 54) | 40 | **28** (A 13, B 9, **C 6**) |
| WP5 headline cells passing (of 36) | 25 | 18 |
| headline combinations passing in all subgroups (of 12) | 5 | 1 |
| headline branches built only on passing cells | 15 | 3 |
| WP5 gate G3 (no A/C regression) | fails | fails |
| WP9 verdict | INCONCLUSIVE (P 0.32–0.73) | INCONCLUSIVE (P 0.38–0.81) |
| verdict score above 0.5: α = 2.0 / α = 2.3 branches (of 18 each) | 18 / 0 | 18 / **7**; the clean α split no longer holds |

The full macro-level comparison has 150 of 216 manuscript macros changed. It is in `git diff d0bd29d~1 d0bd29d -- manuscript/numbers.tex`.

## 3. Integrity checks (adoption gate)

| check | verdict | what it showed |
|---|---|---|
| I1a anchors replay | **PASS_EXACT** | 41 columns identical |
| I1b masses replay | PASS_WITHIN_TOLERANCE | sample cube bit-identical; largest summary difference 3.6×10⁻¹⁵ |
| I1c injection driver | PASS_WITHIN_TOLERANCE | 27 regenerated repair_v7 nodes; ≥ 99.9999 % of draws equal to 10⁻⁶; curves identical |
| I1d WP5 replay | PASS_WITHIN_TOLERANCE | k within 1.9×10⁻¹⁶ in all 54 cells; 0 gate flips |
| I1e WP7 replay | **PASS_EXACT** | the baseline reproduces 8.359927 exactly |
| I2 anchors at the v9 ages | PASS | 0 age mismatches; 96 % of anchors have masses that depend on the branch |
| I3 input and code hashes | PASS | 24 inputs, 41 scripts; one code change, covered by deviation D1 |
| I4 no stale age product | PASS (on rerun) | see below |
| I5 scan vs chain | PASS | weighted-scan versus ledger differences ≤ 0.03 in all 6 cells |

Each "within tolerance" verdict reflects last-bit floating-point differences between macOS arm64 (where repair_v8 was made) and Linux x86_64. Every value reproduces.

**Deviation D1 (`provenance/repair_v9_deviations.json`):**
- I4's first run failed. It substring-searched each execution record for the paths on `wp10_inputs.FORBIDDEN`.
- Three records contain fixed prose naming legacy tables, for example "All 54 branches remain in tables/wp7_ledger.csv". The same sentences appear verbatim in the repair_v8 records.
- The scripts read the `_repair_v9` files through `chain.tag()`.
- The fix: a forbidden path now counts only where a record lists it as a file, meaning a JSON key or a whole string value.
- The fixed check was tested and still catches a forbidden path listed as an input.
- The first-run result is kept in the integrity record as `I4_first_run_defective_check`. No threshold changed.

## 4. Predictions (scored as written; they do not gate)

| | statement | result |
|---|---|---|
| V1a | baseline below repair_v8 | **true** (6.66 < 8.36) |
| V1b | baseline in [4, 9] | **true** |
| V1c | A falls, B falls, C rises | **true** (4.10 → 2.09, 4.26 → 2.85, 0 → 1.72) |
| V2 | closure ratio moves A down, B down, C up (PARSEC) | **true**: 1.006 → 0.924, 1.068 → 1.036, 1.360 → 1.498; MIST moves the same way (reported only) |
| V3 | C's baseline above 0.5 | **true** (1.72) |
| V4 | scan N_death(t) non-decreasing, and below 0.1 for t ≤ 2.53 Myr on the coeval branches | **false** |
| V5 | subgroups agree within ×1.25 at t ≥ 3.15 Myr | **true** |

**Why V4 failed.** Both reasons are physics, not code.
- **24 small drops**, all on MIST between 3.18 and 3.57 Myr. MIST's tabulated turnoff is flat at 69.0 M☉ there: the running-minimum correction documented in `turnoff_sequence`, applied to a 69.0 → 74.9 inversion. With the turnoff fixed, the refitted k decreases slightly and N_death dips by 3–4 %.
- **27 cells above 0.1 at 2.52 Myr**, all MIST. MIST's turnoff at 2.52 Myr is 108 M☉, below the 120 M☉ ceiling, so a few deaths are possible: 0.15–1.27.
- The pre-registration had already noted this risk when it set the 0.1 operating point. The prediction stays failed.

V2's turnoff-only term (2–5 %) pointed the same way as the net change. The refitted k and the census added to it rather than reversing it.

## 5. Age scan (sensitivity, not the headline)

`tables/repair_v9_age_scan.csv` · `data/processed/repair_v9_age_scan_draws.npz` (full N_death distributions, for convolving any (t_A, t_B, t_C)) · `figures/repair_v9/repair_v9_age_scan.png`

N_death mean, PARSEC, R_V 3.1, α 2.3, coeval, with k refitted at each age:

| age (Myr) | ≤ 2.82 | 3.16 | 3.55 | 3.98 | 4.47 | 5.01 | 5.62 | 6.31 |
|---|---|---|---|---|---|---|---|---|
| A | 0 | 0.93 | 2.36 | 4.05 | 5.84 | 7.88 | 9.76 | 11.77 |
| B | 0 | 0.91 | 2.30 | 3.91 | 5.55 | 7.80 | 9.98 | 11.92 |
| C | 0 | 0.93 | 2.34 | 4.09 | 5.86 | 8.25 | 10.52 | 12.71 |
| C passes the WP5 gate? | yes | yes | **no** | no | no | no | no | no |

**Reading:**
- At equal ages the subgroups differ by at most 7 %, because their k values are similar.
- The death history is therefore controlled entirely by the ages. With A, B and C coeval to within 0.4 Myr, resolving the population barely changes the count.
- k falls steadily with assumed age (A: 2,009 at 2.0 Myr, 1,476 at 6.3 Myr). The mechanism was not diagnosed here.
- The full table covers all 54 branches and both families.

## 6. What this does to the project's claims

- **Withdrawn on adoption** (now in CLAUDE.md's withdrawn table):
  - every repair_v8 number above;
  - "C 0.00 deaths";
  - "C sits below the first-death boundary".
- **The 12 August novelty audit's surviving claim does not survive.** Resolving the population changes the inferred history only through age differences, and the adopted ages make the subgroups coeval within their errors. WP13 decides how the paper is framed; its re-based brief now expects the "equivalent" outcome. That expectation was disclosed in brief §5 before anything ran.
- **What survives** is the Gaia census: membership, the subgroups, the distance, the extinction, and the IMF normalization k from the counts. This is the owner's decision 1 of 2026-10-07. The ages behind the ledger are now spectroscopic (Gaia G luminosity with literature spectral-type temperatures), not Gaia photometry alone.
- **New open question (#22 in PROJECT_TRACE §9):** the WP5 low-mass counts are in tension with the spectroscopic ages, in C most of all.

## 7. Manuscript state

- `numbers.tex` regenerated (216 macros), and so are `tables_generated.tex` and figures 1–11.
- The `wp10_inputs` audit passes; `wp10_validate` passes V1–V7.
- Plumbing fixes made on adoption:
  - WP3 extinction is now resolved through the chain's own `wp3_extinction` entry;
  - every repair_v8 product is forbidden on repair_v9;
  - two WP12 quantities that are undefined on repair_v9 render as a dash. The strict α = 2.0 subset has 0 branches, and C is never dead.
- Figure fixes:
  - fig08's hard-coded "Cyg OB2-C contributes nothing on this branch" is now data-driven; it was false on repair_v9;
  - the fig04 and fig09 envelope labels now say "adopted WP4 age envelope".
- **Not done, deliberately:**
  - the `main.tex` prose is Stage 4. It still argues for a young C, and some macro names (`\ageumsA`, `\ageumsC`) now carry spectroscopic ages;
  - no compile (the A&A class files are not on this machine);
  - `AUDIT.txt` not regenerated.
