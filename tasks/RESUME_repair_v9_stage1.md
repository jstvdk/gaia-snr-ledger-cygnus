# repair_v9 Stage 1 (WP4 age redesign): done 2026-10-07 14:15

Full brief: [stage_repair_v9_wp4_age_redesign_and_wp13.md](stage_repair_v9_wp4_age_redesign_and_wp13.md).

Outcome: [reports/wp4v9_age_redesign.md](../reports/wp4v9_age_redesign.md).

## State

- **Stage 1 has been run and scored.**
  - GA1 fails in all 6 cells, so no age is adopted from the redesign.
  - Under the frozen fallback, A and C carry the spectroscopic test-c ages, and B carries a flagged M1 fit.
- **Owner decisions made 2026-10-07** (`provenance/decisions_2026_10_07.json`):
  - A/C spectroscopic ages accepted, with documentation;
  - every subgroup scanned over the native ages 2.0–6.3 Myr;
  - headline: A and C spectroscopic, B borrows A+C.
- **Next:** owner sign-off on [repair_v9_prereg_proposal.md](repair_v9_prereg_proposal.md), then `scripts/repair_v9_prereg.py`, commit, run.

## Uncommitted (commit only when the owner asks)

- `scripts/wp4v9_injections.py`, `wp4v9_fit.py`, `wp4v9_score.py`, `wp4v9_posthoc_diagnostics.py`
- `tables/wp4v9_*.csv`
- `provenance/wp4v9_{injection,fit}_execution.json`, `provenance/wp4v9_age_outcome.json`
- `reports/wp4v9_age_redesign.md`
- the PROJECT_TRACE §1 and §9 #21 edits
- `scripts/decisions_2026_10_07.py`, `provenance/decisions_2026_10_07.json`, `tasks/repair_v9_prereg_proposal.md`
- this note

The data products (`data/processed/wp4_age_posteriors_repair_v9.parquet`, the injection parts and the log) are gitignored.
