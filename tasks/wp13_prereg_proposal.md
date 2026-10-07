# WP13 pre-registration proposal on repair_v9 (draft, for the owner's sign-off)

Drafted on the night of 2026-10-07 under decision O3 ([provenance/decisions_2026_10_07_overnight.json](../provenance/decisions_2026_10_07_overnight.json)). **M0 has not been run.** Nothing in this file is binding until the owner signs it off and it is transcribed into `scripts/wp13_prereg.py`, hashed and committed. That must happen before any M0 number exists.

The basis is the WP13 brief ([wp13_pooled_vs_resolved_ablation_brief.md](wp13_pooled_vs_resolved_ablation_brief.md)), re-based on repair_v9 on 2026-10-07, together with stage brief §5 and the repair_v9 completion report.

## 1. What has changed since the brief was written

- **The premise is gone.** On the adopted chain, C is coeval with A and contributes 1.72 deaths. The brief's falsifiable sentence ("…because Cyg OB2-C lies below the first-death boundary…") is kept and tested as written, but its "because" clause is already false.
- **The outcome is all but known.** This is disclosed prior knowledge from brief §2 as re-based. A single k-weighted common age of 3.516 Myr reproduces the resolved model:
  - N_death 6.644 against 6.657;
  - P(last < 100 kyr) 0.697 against 0.697.
  
  T2 and T3 are expected to fail, and the expected outcome is "equivalent", or "supported but immaterial" if T1 passes.
- **Most decisions below concern the age method.** On repair_v9 the WP4 ages are spectroscopic-HRD (test c), not a Gaia CMD fit, and B's age is not measured.

## 2. Decisions the owner needs to take

**D1. How M0's age is measured (stage brief §5.2: "M0 must use the same age method").**
- **(a) Recommended:** M0's age is the test-c spectroscopic-HRD fit on the pooled spectroscopic sample: A 59 + C 43 + B 4 = 106 stars, eps 0.05, same family and R_V.
  - This is the method M1 uses for A and C.
  - B's 4 stars, unusable on their own, enter the pool.
- (b) The brief's original step 1, the WP4 photometric upper-MS fit on the pooled CMD. **Not recommended:** Stage 1 showed this method does not measure these ages, and it would compare two different age methods.

**D2. T1 (do the data support separate ages?), as computed under D1(a).**
- ln BF(M1 : M0) = ln Z_A + ln Z_C − ln Z_pooled, using the same grid prior.
- B is excluded on the M1 side because it has no measurement of its own. M1 therefore has **two** measured ages, not three, and the BF penalises one extra parameter, not two.
- **Disclosure:** T1 can be computed in seconds from the stored test-c likelihood curves. It has deliberately **not** been computed, because that would be reading an M0-side number.
- The threshold is ln BF ≥ 2.3, as in the brief, on the baseline cell and on ≥ 4 of the 6 family × R_V cells.

**D3. Which branch set T6 uses.**
- The stage brief says T6 is scored on the **18 α = 2.3 headline branches** (decision 2 of 2026-10-02).
- But WP7's adopted headline set (`wp7_alpha_headline_adopt`) is still the **36 retained branches** (α 2.0 and 2.3), and so is CLAUDE.md.
- **Recommendation:** score T6 on the 18 α = 2.3 branches, as the stage brief says, and report the 36 as well. The owner should confirm which set "headline" means for the paper, since the two documents disagree.

**D4. Whether WP13 waits for issue #22.**
- #22 is the tension between the WP5 low-mass counts and the spectroscopic ages: the gate passes in 28 of 54 cells, and in C in 6 of 18.
- M0's k comes from the same counts, so M0 inherits the same tension.
- **Recommendation:** do not wait. The ablation compares M0 with M1 on identical counts, so #22 affects both sides equally. State it as a shared caveat.

**D5. Thresholds T2–T7.** Keep the brief's §5 values unchanged. Re-checked against the repair_v9 baseline:
- T2 (|ΔN| ≥ 3 and ≥ 10 %) is now 45 % of the baseline.
- T3 (|ΔP| ≥ 0.10) is unchanged.
- Changing thresholds after this disclosure would be tuning, so none is proposed.

**D6. The verdict rule and §6's manuscript framing.** Unchanged. The expected "equivalent" row means the following, and the owner should confirm it is acceptable before the run:
- the title becomes "When is a fixed-age supernova count adequate? Cygnus OB2 as a benchmark";
- this requires T7 extended over an (age-spread × young-fraction) grid.

## 3. What the pre-registration will contain once signed

- The input hashes of the repair_v9 chain: ages, masses, WP5 draws, the ledger, and the stored test-c curves.
- M0 under the chosen D1:
  - the pooled age fit;
  - pooled WP5 node injections, using `scripts/repair_v9_injections.py` with a pooled-label age table;
  - the pooled joint fit;
  - the ledger on all 36 headline branches at 2,000,000 iterations, paired with M1 by iteration index.
- T1–T7 with the brief's thresholds. T7 is ledger-level, with 200 realisations per design, plus the (age-spread × young-fraction) grid if D6 is accepted.
- §2 of the brief (as re-based) as disclosed prior knowledge, including the expectation of "equivalent".
- The anti-tuning rules of brief §9.

## 4. Effort once signed

| step | time |
|---|---|
| pre-registration script and commit | about 1 h |
| pooled spectroscopic fit (D1a) | minutes |
| pooled injections (6 family × R_V × 9 nodes = 54) | about 10 min on this machine |
| pooled fit and ledger | about 15 min |
| T1–T6 scoring, the figure | about 1 day of work |
| T7 ledger-level injections, plus the grid if D6 is accepted | 1–2 days |
| report and register row | half a day |
