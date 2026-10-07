# repair_v9 sensitivity items S1 and S2

Date: 2026-10-07; both ran 20:45–20:50 CEST. Pre-registration: `provenance/repair_v9_s1_s2_prereg.json`, committed in `3c571f6` while the repair_v9 chain was still running, before any v9 WP5 or WP7 result existed. Sign-off was delegated by the owner (decision O3). Outcome: `provenance/repair_v9_s1_s2_outcome.json`. **Every check passes and every prediction is true.** Both are sensitivity results, never the headline.

## S1: the IMF upper limit, m_max ∈ {100, 120, 150} M☉

`tables/repair_v9_s1_imf_ceiling.csv` covers the 18 α = 2.3 branches. Each branch has 2,000,000 iterations and uses paired seeds across the three ceilings.

Baseline N_death (R_V 3.1, α 2.3, coeval, all-explode):

| | PARSEC 100 | **120** | 150 | MIST 100 | 120 | 150 |
|---|---:|---:|---:|---:|---:|---:|
| A | 1.39 | 2.09 | 2.77 | 2.45 | 3.15 | 3.81 |
| B | 2.16 | 2.85 | 3.50 | 3.13 | 3.80 | 4.42 |
| C | 1.00 | 1.72 | 2.39 | 2.64 | 3.34 | 3.99 |
| association | **4.55** | **6.66** | **8.66** | 8.22 | 10.28 | 12.22 |

**Checks and predictions:**
- **S1-I1, pass.** At 120 M☉ the run reproduces the repair_v9 ledger: 6.655 against 6.657. The ceiling patch does nothing else.
- **S1-P1, true.** The count rises with m_max in every branch and subgroup.
- **S1-P2, true.**
  - The ratios are ×0.683 at 100 M☉ and ×1.301 at 150 M☉.
  - The analytic expectation at the MAP ages was ×0.622 and ×1.356; the tolerance was ±30 %.
  - As predicted, the spread of the age posteriors pulls the measured ratios towards 1.

**Reading:**
- On the repair_v9 ages, the IMF ceiling is worth −32 % / +30 % on the death count. That is larger than the brief's −17 % / +15 % expectation, which used the repair_v8 inputs.
- The reason is that the turnoffs (68–82 M☉ at the baseline ages) now sit close to the ceiling.
- **The 120 M☉ convention is therefore one of the larger systematics in the ledger.** It is comparable to the PARSEC-versus-MIST difference (6.66 against 10.28). The paper's sensitivity section should carry it.

## S2: subgroup-assignment uncertainty

`tables/repair_v9_s2_responsibilities.csv` · `tables/repair_v9_s2_label_uncertainty.csv` · `data/processed/wp5_imf_normalization_repair_v9_s2soft.parquet`

**Method:**
- The frozen WP2 mixture was replayed: 1,331 clean labelled members, the 50 frozen seeds, components named by the frozen rule.
- Responsibilities were averaged over the seeds.
  - Median maximum responsibility: **0.903**.
  - Stars with maximum responsibility above 0.9: **50.6 %**. This reproduces the brief's 0.90 / 51 %.
  - For every star, the most likely component is its hard label.
- Each star then enters every subgroup with weight membership × responsibility. Membership weight moves from A to C (A 422 → 405, C 366 → 386; B 360 → 358).
- WP5 was refitted on the repair_v9 ages and responses, and the WP7 engine run on the hard and the soft draws with identical seeds.

Baseline (PARSEC, R_V 3.1, α 2.3, coeval):

| | k hard → soft | change | N_death hard → soft |
|---|---|---:|---|
| A | 1,769 → 1,690 | −4.5 % | 2.09 → 1.69 |
| B | 1,693 → 1,682 | −0.7 % | 2.85 → 2.86 |
| C | 1,757 → 1,877 | +6.8 % | 1.72 → 2.20 |
| association | | | **6.66 → 6.74 (+1.2 %)** |

**Checks and predictions:**
- **S2-I1, pass.** Total membership weight is conserved exactly (1,148.39).
- **S2-P1, true.** The association count is within 10 % (+1.2 %).
- **S2-P2, true.** Every subgroup's k is within 15 % (largest change: C, 6.8 %).
- **S2-P3, true.** B's k changes least.

**Reading:**
- Label uncertainty hardly moves the association's death history.
- It does move deaths between A and C, by about 0.4–0.5 each.
- Any per-subgroup death count should carry roughly ±0.5 for label uncertainty, which is comparable to the A–C difference itself.
- **The approximation stated in the pre-registration:** a relabelled star keeps the mass posterior computed at its hard subgroup's age. That age differs by at most 0.4 Myr.
