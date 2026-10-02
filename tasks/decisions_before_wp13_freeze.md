# Decisions needed before the WP13 thresholds are frozen

*Prepared 2026-10-02 for the project owner (stage brief
[stage_after_issue19_c_age_and_wp13.md](stage_after_issue19_c_age_and_wp13.md),
step 3). Each row lists the options and what each one does. Nothing here has
been decided. Record the choices in PROJECT_TRACE §9 and they become inputs to
`scripts/wp13_prereg.py`.*

Why these must come first: T6 requires the WP13 result to hold on ≥ 80 % of
"the 36 headline branches". That set is defined by dropping α = 2.6, so
decisions 1–2 fix the denominator WP13 is scored against. Decision 4 fixes how
the WP13 result can be framed.

## 1. WP5 gate G3 on repair_v8: does the no-regression clause apply to a data-defect fix?

**Facts** ([wp5_imf_fit_execution_repair_v8.json](../provenance/wp5_imf_fit_execution_repair_v8.json)):
- G3 fails *only* through `no_A_or_C_branch_regression`;
- the baseline passes in all three subgroups, masses are within a factor of two, and the grid is 40/54;
- the two regressions are both A-PARSEC α = 2.0, at R_V 3.0 and 3.1, where the max |residual| rises 2.13 → 2.36 and 1.80 → 1.98;
- in the other direction, two α = 2.6 cells (A, R_V 3.5) went fail → pass;
- [CUTS §14.7](../CUTS_AND_THRESHOLDS.md) says every G3 evaluation reports the strict reading *and* the refined reading, in which a cell indeterminate in both versions is not a regression. The repair_v8 record reports only the strict reading. `wp5_verdict_stability.py` supports only repair_v4–v6, so the refined reading has not been computed.

| option | what it means | consequence |
|---|---|---|
| 1a. Clause applies; G3 stays FAILED | WP5 is reported as failing G3 on repair_v8, and repair_v8 remains adopted (its adoption never depended on G3) | Manuscript and PROJECT_TRACE say "G3 fails its no-regression clause on repair_v8". No number moves |
| 1b. Clause does not apply to defect fixes | Needs a written, dated amendment to the G3 scope: "regressions are judged against the defect-free predecessor, which does not exist for a defect fix" | G3 recorded as PASS-with-amendment. Being post hoc, the amendment must be labelled as such |
| 1c. First compute the refined §14.7 reading (prerequisite to 1a/1b) | Extend `wp5_verdict_stability.py` to repair_v7/v8 with versioned outputs; about ½ day plus compute | If both A cells are indeterminate in both versions, the refined reading passes and the question may become moot. If either is determinate, 1a or 1b is still needed |

## 2. The α headline set (α ∈ {2.0, 2.3}, dropping 2.6)

**Facts** ([wp5_alpha_plausibility_execution_repair_v8.json](../provenance/wp5_alpha_plausibility_execution_repair_v8.json)):

| statistic | α = 2.0 | α = 2.3 | α = 2.6 |
|---|---:|---:|---:|
| calibration-window median χ² (E1) | 10.4 | 7.00 | 10.06 |
| closure median σ from unity (E2, the unspent out-of-sample line) | 6.55 | 1.03 | 5.63 |

- The registered test (A1: 2.3 vs 2.6) still passes.
- On both E1 and E2, α = 2.0 now fits association-wide no better than 2.6, which the set drops.
- C alone prefers 2.0 on both statistics (E1 5.34, E2 1.07σ).
- That record's own "reading" text still says 2.6 "is the worst on the median". The sentence was not regenerated and is now wrong.

| option | consequence |
|---|---|
| 2a. Keep the 36-branch set as registered, and say so | No rerun. The text states plainly that the set keeps a slope fitting as badly as the one it drops. WP13's T6 is scored on the same 36 branches |
| 2b. Keep 36 as the headline and add α = 2.0 as an explicit sensitivity | The headline range is unchanged (5.57–28.7). The α = 2.0 half is shown separately: it spans 15.2–28.7 against 5.6–11.2 at α = 2.3, and it carries the WP9 split (18/18 vs 0/18). Source: `tables/wp7_ledger_repair_v8.csv` |
| 2c. Revise the headline set (e.g. α = 2.3 only, or all 54) | A new pre-registered decision, not a retune. It changes the headline range and T6's denominator. It must be frozen before any WP13 number is read |

## 3. α-headline adoption D1-P4 (threshold hard-codes 8.43)

**Facts** ([wp7_alpha_headline_adoption_outcome_repair_v8.json](../provenance/wp7_alpha_headline_adoption_outcome_repair_v8.json)):
- the literal check (`baseline = 8.43`) FAILS on repair_v8 (8.36);
- the chain reading ("the restriction does not move the baseline") PASSES and is recorded as governing;
- both readings are already stored.

| option | consequence |
|---|---|
| 3a. Confirm the chain reading governs | Nothing to rerun. Record the confirmation in PROJECT_TRACE |
| 3b. Treat the literal failure as binding | D1 is re-opened on repair_v8: a new pre-registration of the α restriction with a chain-relative threshold |

## 4. Issue #20 Phase B (new, from today's Phase A)

**Facts** ([issue20_phase_a_report.md](../reports/issue20_phase_a_report.md)):
- the frozen rule gives **H3, flagged "not specific to C"**, because the same mixture test passes on A;
- the mixture's second component sits at the 10 Myr grid edge in both subgroups;
- C's single-age spectroscopic age is 3.98 / 4.01 Myr, at least as old as A's;
- C's young photometric age is carried by its photometric-only members' repair_v1 extinction.

| option | consequence for WP13 |
|---|---|
| 4a. Accept H3; two-component C in repair_v9 after WP13 | WP13 runs on repair_v8, and its premise ("C below the first-death boundary") is flagged as contested |
| 4b. Treat as INCONCLUSIVE (single-age); carry C as a young/coeval branch | WP13 runs on repair_v8. No novelty claim rests on C's age, so a "material improvement" outcome cannot be attributed to C's youth |
| 4c. Pre-register Phase A′ first (outlier-robust HRD likelihood with a grid-edge guard, plus a near-IR extinction check of C's photometric upper-MS stars) | About 1–2 days, read-only. It decides "photometry vs spectroscopy" before WP13's framing is fixed |
