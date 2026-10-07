# Overnight plan, 2026-10-07 → 08 (unattended)

The owner left the machine running overnight and asked that nothing wait for their input. The decisions are in [provenance/decisions_2026_10_07_overnight.json](../provenance/decisions_2026_10_07_overnight.json). Any session that resumes this work reads this file first and ticks items off here.

## State at hand-off (19:05 CEST)

- Pre-registration `provenance/repair_v9_prereg.json` was committed in `d947f31` before any chain step.
- The chain runs in tmux session `repair_v9` (`bash scripts/run_repair_v9_chain.sh`). Its log is `data/processed/repair_v9_logs/chain.log`, and the line `EXIT=<code>` marks the end.
- Python: `/home/vvoitsek/software/installations/miniconda3/envs/cygob2-gaia/bin/python`. Run with `PYTHONPATH=scripts`.

## Steps

1. [x] **Watch the chain.** Done 20:37; one crash-free run. Phases: A (anchors and masses), B (injections), C (fits and ledgers), D (WP6–WP12, checks, scoring).
   - On a crash caused by a code bug:
     - fix the bug;
     - add an entry to `provenance/repair_v9_deviations.json` (`path`, `old_sha256`, `new_sha256`, `reason`, `date`);
     - relaunch the runner. It is resumable, and finished steps are skipped.
   - Never loosen a threshold in `repair_v9_integrity.TOL` or the scoring constants.
2. [x] **When the chain ends:** I4 failed on a check defect, fixed as deviation D1 and rerun; all of I1–I5 pass. read `provenance/repair_v9_integrity.json` and `provenance/repair_v9_outcome.json`.
   - If a check failed because of a bug in the check or in the plumbing, fix it, record the deviation, rerun that check and rescore.
3. [x] **If I1–I5 all pass (decision O2):** done (`d0bd29d` and the docs commit).
   - set `ADOPTED = "repair_v9"` in `scripts/chain.py`;
   - write `reports/repair_v9_completion_report.md` with a before/after table;
   - in `PROJECT_TRACE.md`: add a §1 notice, close #21, update #20;
   - in `CLAUDE.md`: update the headline numbers and add the repair_v8 values to the withdrawn table;
   - manuscript plumbing: `wp10_inputs` versions, `wp10_numbers`, `wp12_tables`, `wp12_figures`, `wp10_validate`. Leave `main.tex` prose alone.
   - Commit (decision O1: no push).
4. [—] **If not adopted:** not applicable. report it, keep repair_v8, skip S1 and S2, and still do the WP13 re-base draft. Commit.
5. [ ] **S1, IMF ceiling 100/120/150 M☉:** pre-register (script, JSON, commit), run, report, commit.
6. [ ] **S2, subgroup-assignment uncertainty:** pre-register (WP2 mixture refit: k = 3, full covariance, StandardScaler on l, b, μα\*, μδ, 50 seeds, consensus mapping to the stored labels), run, report, commit.
7. [ ] **WP13:**
   - re-base `tasks/wp13_pooled_vs_resolved_ablation_brief.md` on repair_v9 (pattern: `scripts/wp13_brief_rebase.py`);
   - draft `tasks/wp13_prereg_proposal.md` for the owner's sign-off.
   - **Do not run M0.** Commit.
8. [ ] **Morning summary** for the owner, at the top of the completion report: what ran, what passed or failed, what needs them.

S3 (the parallax-blind membership check) is not done; that is recorded in decision O3.
