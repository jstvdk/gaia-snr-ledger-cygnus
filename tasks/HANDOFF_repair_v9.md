# Handoff: repair_v9 (WP4 age redesign → repair_v9 chain → WP13)

> **SUPERSEDED 2026-10-08.** Everything below was executed on 2026-10-07: repair_v9 is adopted. The current handoff is §0 of [reports/repair_v9_completion_report.md](../reports/repair_v9_completion_report.md). Kept for the record.

Written 2026-10-07 for an agent continuing this work on another machine. Read it top to bottom before touching anything.

## 0. Read first

1. `CLAUDE.md`: project rules. **Never read `AUDIT.txt`** (16 MB). Never type a number into `manuscript/main.tex`. Nothing is overwritten or retuned. Resolve inputs through `scripts/wp10_inputs.py`.
2. The brief: [stage_repair_v9_wp4_age_redesign_and_wp13.md](stage_repair_v9_wp4_age_redesign_and_wp13.md). Rule 6 there: **commit only when the owner asks; show the owner each pre-registration for sign-off before running a stage.**
3. Stage 1 outcome: [../reports/wp4v9_age_redesign.md](../reports/wp4v9_age_redesign.md).
4. The next pre-registration, awaiting sign-off: [repair_v9_prereg_proposal.md](repair_v9_prereg_proposal.md).
5. `PROJECT_TRACE.md`: §1 notices (top) and §9 issues #20 and #21.

The owner is called "the owner"; use they/them. Reports are written for someone who skims: lead with the answer.

## 1. Check this before anything else: the data tree

`data/` (about 8.4 GB) is **gitignored and exists only on the owner's original machine**. It has no archive. Everything below that runs code needs it, including:
- `data/processed/wp4_age_posteriors_repair_v9.parquet`;
- `data/processed/wp4v9_injection_parts/`;
- the WP3 posterior `wp3_extinction_posterior_repair_v5.npz`.

If it is missing on this machine, you can still read and edit documents and write scripts. **Do not regenerate inputs to work around it.** Tell the owner and ask them to copy `data/`.

Also check `git log --oneline -5`. The last commit on the original machine was `2e4088b` (Stage 1 pre-registration). Everything from Stage 1's run onward is **uncommitted**, so the files in §6 exist only on that machine's working tree until the owner commits and pushes.

## 2. Where the project stands

**Question under test:** is subgroup CygOB2-C young (so it has had no deaths yet)? That was the project's one surviving novelty claim.

**Chain history:**
- The quoted chain is `repair_v8` (N_death 8.36: A 4.10, B 4.26, C 0.00).
- Issue #21: WP4's photometric age fit treated the A_V error as independent in G and in colour. It should be correlated along the reddening vector.

**Phase A′ (2026-10-06)** ran three independent age tests. Only C-MIST was resolved, at 4.45 Myr. It exposed #21.

**Stage 1** was a pre-registered redesign of the WP4 age fit ("M1"): marginalise over each star's WP3 A_V posterior, use a bright window M_G0 ≤ −1, and treat the window as a selection function. The pre-registration is `provenance/wp4v9_age_prereg.json`, committed in `2e4088b`. `scripts/wp4v9_common.py` is hash-frozen: **do not edit it.**

**Stage 1 result (2026-10-07):**
- GA1 (injection bias < 0.2 Myr at every true age, coverage 0.55–0.80) **fails in all six subgroup × family cells.**
  - Mean bias is −0.2 to −0.5 Myr at true ages of 3–5 Myr, from a young-biased tail; coverage is fine (0.65–0.71).
  - The redesign is still about three times less biased than the old WP4 fit in the bright window.
- On real data the M1 age is unstable: the MAP jumps between about 1 Myr and 5–10 Myr under R_V, anchor removal, the ×2 A_V widening and an outlier term.
  - In A, one unresolved triple (O7 I + O6 I + O9 V, Gaia DR3 2067830941174418048) supplies about 11 of the 14 log-likelihood units favouring 1.26 Myr.
  - These post-hoc checks (`scripts/wp4v9_posthoc_diagnostics.py`) are not pre-registered and gate nothing.
- GA2 fails for A and C; GA3 passes but has no power against too-young ages.
- Predictions: W1 true, W2 true, W3 false, W4 PASS. W1 and W2 are scored on the spectroscopic fallback, so they add no independent evidence.
- **The Gaia optical CMD does not measure these subgroups' ages.** This is a result and is documented as such.

**Frozen fallback** (`if_GA1_or_GA3_fails`): repair_v9 carries the Phase A′ test-c spectroscopic-HRD posteriors (eps 0.05):

| | PARSEC | MIST |
|---|---|---|
| A, 59 stars | 3.16 Myr (68 %: 2.80–3.50) | 3.57 (3.31–3.70) |
| C, 43 stars | 3.55 (3.32–3.97) | 4.01 (3.76–4.27) |
| B (flagged M1 fit, **unusable**) | 5.01 | 4.50 |

The B row varies from 1.4 to 5.7 Myr across R_V, so B cannot use it. **C is coeval with A, not young, so the "C below the first-death boundary" claim does not survive.**

## 3. Owner decisions (recorded in `provenance/decisions_2026_10_07.json`)

1. **Accept A and C spectroscopic ages, with full documentation.**
   - The owner worried this diminishes the Gaia work.
   - The answer already given: membership, subgroups, distance, extinction, the census that sets the IMF normalisation k, and the kinematics all stay Gaia products. Only the claim of a Gaia-photometry-only *age* is lost.
   - Keep the documentation honest on this point: the spectroscopic ages use Gaia G on the luminosity axis but literature spectral-type temperatures on the other.
2. **Every subgroup is scanned over all plausible ages:** the 11 native isochrone ages from 2.0 to 6.3 Myr per family, each subgroup independently, with k refitted at each age.
   - PARSEC 2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01, 5.62, 6.31.
   - MIST 2.00, 2.25, 2.52, 2.83, 3.18, 3.57, 4.01, 4.50, 5.05, 5.67, 6.37.
3. **Headline weighting:** A and C use their own spectroscopic posteriors. B borrows the A + C combined spectroscopic posterior (product of the two test-c likelihood curves, same family and R_V), flagged "not measured for B".

**Not yet decided or signed off by the owner:**
- the repair_v9 pre-registration proposal (§4);
- committing any of the work in §6.

## 4. Your next steps, in order

1. **Get the owner's sign-off on [repair_v9_prereg_proposal.md](repair_v9_prereg_proposal.md).** Do not run any chain step before that. The proposal covers:
   - the headline ages (B = A + C combined);
   - the repair_v8-style chain with a new `--age-version` argument for `wp4_mass_posteriors_repair.py`;
   - the age scan (198 WP5 node injections plus 11 mini-chains);
   - integrity checks I1–I5;
   - predictions V1–V5, which state that **baseline N_death is expected to FALL, not rise**.
2. Derive the expected sign of each closure-ratio change from `wp6_closure_test` (V2 in the proposal) before writing the pre-registration.
3. Write `scripts/repair_v9_headline_ages.py`.
   - It builds B's A + C product curve and writes `data/processed/wp4_age_posteriors_repair_v9_headline.parquet`.
   - It must not edit the Stage 1 table.
   - Then write `scripts/repair_v9_prereg.py`, modelled on `scripts/issue19_prereg.py`. It must hash inputs, refuse to overwrite, and refuse to run if outputs exist.
4. **Ask the owner to commit** the pre-registration (and the Stage 1 files, §6) before running. Commit only when asked.
5. Run the chain, modelled on `scripts/run_repair_v8_chain.sh`, as `scripts/run_repair_v9_chain.sh`:
   - `export CYGOB2_CHAIN=repair_v9 WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging`;
   - add `"repair_v9"` to `scripts/chain.py` and leave `ADOPTED` as `repair_v8` until the adoption rule passes (proposal §6);
   - node responses must be regenerated, not reused; read `run_repair_v6_chain.sh` for how;
   - WP7 at 2,000,000 iterations, WP11 at 500,000.
6. Score I1–I5 and V1–V5, write `reports/repair_v9_completion_report.md` with a before/after table, close #21, update #20, and update `CLAUDE.md` (headline numbers and the withdrawn table).
7. Only then, WP13 (brief §5) with its own pre-registration and sign-off.

**Rough estimates, before running:**
- repair_v8's k at fixed best-fit ages (PARSEC, R_V 3.1, α 2.3) gives about A 0.9, B 2.2, C 2.6, total about 5.7.
- Averaging over age posteriors will raise this, because N_death is convex in age.
- Deaths per subgroup go from 0 at 2.5 Myr or younger to about 1 at 3.2, 4 at 4.0 and 8 at 5.0 Myr.
- Runtime is about 9–10 h, unattended.

## 5. How to run long jobs

- Detached: `screen -dmS <name> caffeinate -i <script>`. It survives a closed session.
- Make steps resumable in blocks: write a part file per block through `.tmp` then rename, and skip existing parts on restart. `scripts/wp4v9_injections.py` is the working example, and `data/processed/wp4v9_logs/run_stage1.sh` is a runner example.
- Use about 13 workers on a 14-core machine.
- A session that ends kills any non-detached job. This happened once already and lost an unsaved run.
- The owner once asked for the laptop to sleep after a run (`pmset sleepnow`). That is not a standing request; ask first.

## 6. Uncommitted work from this session

Nothing here is committed.

- scripts: `wp4v9_injections.py`, `wp4v9_fit.py`, `wp4v9_score.py` (edited after the pre-registration, see below), `wp4v9_posthoc_diagnostics.py`, `decisions_2026_10_07.py`
- tables: `tables/wp4v9_*.csv`
- provenance: `wp4v9_injection_execution.json`, `wp4v9_fit_execution.json`, `wp4v9_age_outcome.json`, `decisions_2026_10_07.json`
- reports and tasks: `reports/wp4v9_age_redesign.md`, `tasks/repair_v9_prereg_proposal.md`, `tasks/RESUME_repair_v9_stage1.md`, this file
- `PROJECT_TRACE.md` (§1 notice and §9 issue #21 status)

Unrelated uncommitted files (draft figures, `reports/deck_slide_guide.md`, `reports/monte_carlo_methodology_explained.md`, and similar) are the owner's; leave them.

**One deviation to know about, already disclosed in the report §2:**
- Before any real-data fit, `wp4v9_score.py` was changed so the fallback spectroscopic age is written as the adopted `ums` row. The draft stored it only as a "candidate", which left keys with no adopted row.
- The pre-registration's output schema and `if_GA1_or_GA3_fails` require the adopted row. No criterion or threshold changed.

## 7. Rules for this stage

- Pre-registered gates and predictions are scored as written. Failed predictions stay failed (here W3).
- Do not loosen GA1 or any threshold to rescue the redesign.
- Do not quote any repair_v9 number as a result until adoption.
- Do not write that C is young or free of deaths.
- Superseded products stay on disk.
