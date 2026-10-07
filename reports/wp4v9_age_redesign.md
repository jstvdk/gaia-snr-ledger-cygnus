# repair_v9 Stage 1: WP4 age redesign: outcome

Date: 2026-10-07.

The pre-registration was committed in `2e4088b` (`provenance/wp4v9_age_prereg.json`) before any real-data fit was run.

The brief is [stage_repair_v9_wp4_age_redesign_and_wp13.md](../tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md). The issue is #21 in `PROJECT_TRACE.md` §9.

## 1. Bottom line

1. **GA1 fails in all six subgroup × family cells.**
   - This is the injection-validation gate: the bias limit is 0.2 Myr at every true age.
   - Coverage is fine everywhere (0.65–0.71). The failures are a young-biased tail at true ages of 3–5 Myr, worth −0.2 to −0.5 Myr.
   - **By the frozen rule, no subgroup's age is adopted from the redesign.**
2. **On real data, the photometric bright-window age is unstable, not just biased.** The fitted age jumps between about 1 Myr and 5–10 Myr under:
   - the R_V branch;
   - removing the spectroscopic members;
   - the ×2 widening of the A_V posterior;
   - an outlier term;
   - dropping the composite spectroscopic stars.

   One star dominates A's fit (§4). The ×2 widening flags 13 of 18 subgroup × family × R_V cells as "double-counting sensitive".

   **Conclusion: the Gaia optical CMD of the upper main sequence does not measure these subgroups' ages, given A_V posteriors 0.5–1 mag wide.** This is a negative result about the method, not a new age.
3. **What repair_v9 carries under the frozen rule (`if_GA1_or_GA3_fails`):**
   - for A and C, the Phase A′ test-c spectroscopic-HRD posteriors;
   - for B, the M1 bright-window posterior, flagged "failed GA1/GA3".

   | | PARSEC | MIST |
   |---|---|---|
   | A (59 spectroscopic stars) | 3.16 Myr (68 %: 2.80–3.50) | 3.57 Myr (3.31–3.70) |
   | C (43 spectroscopic stars) | 3.55 Myr (3.32–3.97) | 4.01 Myr (3.76–4.27) |
   | B, photometric, **flagged** | 5.01 Myr (4.02–5.03) | 4.50 Myr (4.06–4.50) |

   R_V 3.1, f_bin 0.4. A and C come from `tables/issue20b_age_tests.csv` (test c, ε 0.05). B is the M1 bright-window fit.
4. **C is not younger than A.** It is the same age within 0.4–0.45 Myr, slightly older in both families. This is the same answer as Phase A and Phase A′ test c, because it is the same measurement. W1 and W2 come out true, but they add no independent evidence (§6).
5. **B's carried age is not usable as it stands.** At f_bin 0.4 it reads 1.41 / 4.50 / 5.67 Myr (MIST) and 5.01 / 5.01 / 3.55 Myr (PARSEC) across R_V 3.0 / 3.1 / 3.5. B has only 4 spectroscopic stars, so it has no fallback. **This needs an owner decision** before repair_v9 (§8).

## 2. What ran

| step | script | output | runtime |
|---|---|---|---|
| injections | `scripts/wp4v9_injections.py` | `tables/wp4v9_injection_validation.csv`, `provenance/wp4v9_injection_execution.json` | 2 h 12 min (6 resumable blocks) |
| real-data fits | `scripts/wp4v9_fit.py` | `tables/wp4v9_real_fits.csv`, `tables/wp4v9_faint_end_study.csv`, `provenance/wp4v9_fit_execution.json` | 22 min |
| scoring | `scripts/wp4v9_score.py` | `data/processed/wp4_age_posteriors_repair_v9.parquet`, `provenance/wp4v9_age_outcome.json` | 3 s |
| post-hoc diagnostics (**not pre-registered**) | `scripts/wp4v9_posthoc_diagnostics.py` | `tables/wp4v9_posthoc_star_contributions.csv`, `tables/wp4v9_posthoc_refits.csv` | 3 min |

**Method code:** `scripts/wp4v9_common.py` matches its frozen hash, and every script checks it. There is no method deviation.

**One change to the scoring bookkeeping**, made on 2026-10-07 before any real-data fit existed:
- The draft scorer stored the fallback spectroscopic age only as a "candidate" row. That left keys with no adopted row.
- The pre-registration's `output_schema` requires exactly one adopted `ums` row per key, and `if_GA1_or_GA3_fails` says repair_v9 carries the test-c posterior.
- The scorer now writes it that way. The failed M1 row is kept as `ums_not_adopted`.
- No criterion or threshold changed.
- The scorer was dry-run on mock inputs first.

**The run was interrupted once:**
- On 2026-10-06 a session end killed the run with nothing written.
- On 2026-10-07 the owner paused it at 10:28 to charge the laptop, and it resumed at 11:41.
- Each block is seeded, so a resumed block is identical to an uninterrupted one.

## 3. GA1: injection validation

M1, bright window (M_G0 ≤ −1), R_V 3.1, f_bin 0.4. 100 realisations per true age.

| subgroup | family | bias at 2.5 / 3.2 / 4.0 / 5.0 Myr | coverage 68 % | GA1 |
|---|---|---|---|---|
| A | PARSEC | +0.16 / −0.22 / −0.39 / −0.48 | 0.70 | FAIL |
| B | PARSEC | +0.20 / +0.08 / −0.29 / −0.23 | 0.71 | FAIL |
| C | PARSEC | +0.16 / −0.06 / −0.34 / −0.28 | 0.71 | FAIL |
| A | MIST | −0.01 / −0.28 / −0.13 / −0.34 | 0.67 | FAIL |
| B | MIST | +0.14 / −0.09 / −0.03 / −0.22 | 0.65 | FAIL |
| C | MIST | −0.05 / −0.28 / −0.11 / −0.14 | 0.68 | FAIL |

**The bias comes from a tail, not a shift of the whole distribution:**
- In 45–74 % of realisations the MAP lands on the true native age.
- 10–20 % land at 1.3–2.5 Myr for true ages of 3–5 Myr.
- Few land too old, so the mean of the medians is pulled young.

**Compared with the old WP4 model:**
- The redesign is a large improvement in the bright window. Mean bias is −0.13 (MIST) and −0.14 (PARSEC) Myr, against −0.38 and −0.57 for the old independent-error model.
- In the faint and full windows both models are unbiased on synthetic data (|mean| ≤ 0.06 Myr).

## 4. Real-data fits: why they are unstable

### M1 bright window, f_bin 0.4: MAP (Myr)

| | R_V 3.0 | 3.1 | 3.5 | 3.1, ×2 A_V width | 3.1, non-anchor stars only |
|---|---|---|---|---|---|
| A PARSEC | 1.26 | 1.26 | 3.98 | 1.26 | 10.0 (railed) |
| A MIST | 1.12 | 1.12 | 4.50 | 1.12 | 10.1 (railed) |
| B PARSEC | 5.01 | 5.01 | 3.55 | 1.26 | 5.01 |
| B MIST | 1.41 | 4.50 | 5.67 | 1.41 | 4.50 |
| C PARSEC | 5.01 | 5.01 | 1.00 | 1.00 | 1.00 |
| C MIST | 2.83 | 3.57 | 1.00 | 1.00 | 1.00 |

A's 1.1–1.3 Myr contradicts its own O-star spectroscopy (3.2–3.6 Myr). GA3 cannot catch this: it only rejects ages that are too *old* for the most massive dwarf or giant (§5).

### Post-hoc diagnostics: per-star ln L contributions (`tables/wp4v9_posthoc_star_contributions.csv`)

**A:** a single star dominates.
- The unresolved multiple system O7 I + O6 I + O9 V (Gaia DR3 2067830941174418048, M_G0 = −8.3) contributes Δln L = +11.3 (PARSEC) / +11.8 (MIST) towards 1.26 Myr. A's total is +13.7 / +11.9.
- No single-star or binary model reaches that luminosity at about 3.5 Myr.
- The spectroscopic members' A_V posteriors are narrow (68 % half-width about 0.03 mag), so such a star cannot hide in its extinction. It sets the age.

**C:** the pull runs the other way.
- C's spectroscopic stars favour about 3.6 Myr (Δln L −6 / −18).
- The photometric-only members just inside the window edge favour 1.0 Myr (+10 / +12). These are stars with M_G0 between −1.1 and −1.6, A_V 5.7–7.7 and wide posteriors.
- The fit is a tug of war between the two sets of stars. The R_V branch and the posterior width decide who wins.

**B:** the photometric stars favour the older solution. One composite supergiant (O7 I + O9 I) pulls young.

### Post-hoc refits (`tables/wp4v9_posthoc_refits.csv`)

MAP in Myr, at R_V 3.0 / 3.1 / 3.5.

| | with outlier term ε 0.05 | without composite spectroscopic stars |
|---|---|---|
| A PARSEC | 1.26 / 3.98 / 10.0 | 2.82 / 7.08 / 8.91 |
| A MIST | 2.00 / 2.52 / 10.1 | 4.50 / 6.37 / 10.1 |
| B PARSEC | 1.58 / 1.41 / 1.00 | 5.01 / 5.01 / 3.55 |
| B MIST | 1.59 / 1.41 / 1.00 | 4.50 / 4.50 / 6.37 |
| C PARSEC | 1.12 / 1.12 / 7.94 | 5.01 / 5.01 / 1.00 |
| C MIST | 1.00 / 1.12 / 10.1 | 2.83 / 3.57 / 1.00 |

Neither repair stabilises the fit. These refits explain the instability; they are not candidates for adoption.

**Why the injections could not see this.** The synthetic subgroups are drawn from the model itself. They contain no unresolved triples and no mis-specified A_V posteriors, and their extinction is consistent with their colours by construction. The likelihood is informative on synthetic data and fragile on real data. That is exactly the double-counting and model-misspecification risk the owner asked to be tested.

## 5. GA2, GA3 and the adoption rule

**GA2** (A and C: does the M1 68 % interval overlap the test-c 68 % interval?) fails in all four cells at the decision cell:
- A PARSEC 1.13–1.46 vs 2.80–3.50;
- A MIST 1.12–1.41 vs 3.31–3.70;
- C PARSEC 4.02–4.97 vs 3.32–3.97;
- C MIST 3.18–3.57 vs 3.76–4.27.

GA2 is robust over R_V × f_bin except for A-PARSEC, which overlaps at R_V 3.5 (f_bin 0.4 and 0.5).

**GA3 passes everywhere,** with posterior fraction below the turnoff bound equal to 0.000. The bound m_lo at R_V 3.1 is set by these stars:

| | PARSEC | MIST |
|---|---|---|
| A, from an O7 III star | 22.8 M☉ | 23.6 M☉ |
| B, from an O9.7 III star | 15.0 M☉ | 15.3 M☉ |
| C, from an O5 III star | 33.3 M☉ | 31.5 M☉ |

These bounds only exclude ages well above 5 Myr, so GA3 had no power against the young failure seen here.

**Adoption:** GA1 fails, so every cell is "not adopted from the redesign" (`if_GA1_or_GA3_fails`). The joint and branch fallbacks are never reached, because they apply only after GA1 and GA3 pass.

**The table `wp4_age_posteriors_repair_v9.parquet`** has 680 rows:

| indicator | rows | content |
|---|---|---|
| `ums` (adopted) | 66 | exactly one per (subgroup, family, R_V, f_bin, dmu). A and C carry test c; B carries the flagged M1 fit. |
| `ums_not_adopted` | 44 | the failed M1 rows for A and C |
| `ums_sensitivity` | 570 | every other fit |

Test c has no distance or f_bin variants. So the dmu ±0.06 and f_bin 0.3 / 0.5 keys of A and C repeat the same spectroscopic curve, marked `no_dmu_refit` in `fallback`.

## 6. Predictions W1–W4 (frozen in the pre-registration and scored as written)

| | statement | result | comment |
|---|---|---|---|
| W1 | C's adopted age ≥ 3.2 Myr in both families | **true** (3.55 / 4.01) | scored on the test-c fallback, so it repeats Phase A′ rather than confirming it independently |
| W2 | \|A − C\| < 1 Myr in both families | **true** (0.39 / 0.44) | same caveat |
| W3 | the old model is biased young in the faint window on injections | **false** (−0.02 / −0.02 Myr) | failed as predicted; the old model's young bias is a bright-window effect (−0.38 / −0.57) |
| W4 | no single cause reconciles the faint window | **PASS** | each cause overlaps the bright window in at most 1 of 4–5 cells |

**A side observation from M4 (iii):** keeping only stars with P > 0.9 moves C's faint-window age to 2.5 Myr in both families (PARSEC 2.28–2.52, MIST 2.45–2.60). The same cut moves A to 5.6 / 7.7 Myr. The faint window in C is therefore sensitive to the membership cut. This is reported only, and it is not evidence that C is young: the faint window is the region issue #21 showed to be mis-modelled.

## 7. What this changes

- **Issue #21's fix does not yield usable photometric ages.** Correctly propagating the A_V uncertainty, the redesign's purpose, shows the Gaia optical photometry carries too little age information to resist a handful of over-luminous or over-reddened stars.
- **The only age measurements that survive are spectroscopic-HRD ages.** These are A, with 59 stars, and C, with 43. B has none it can use.
- **C is not young on any surviving measurement.** The paper's remaining novelty claim (that the C subgroup sits below the first-death boundary, CLAUDE.md "Where it actually stands") **does not survive**: C is coeval with A to within the errors.
  - The brief's V1 estimate puts the baseline N_death at +2.8 (at 3.6 Myr) to +4.6 (at 4.0 Myr) above repair_v8's 8.36, holding A and B fixed.
  - A and B will move too, so repair_v9 computes the actual numbers.

## 8. Decisions needed before repair_v9 (brief §4.1)

1. **Accept the frozen outcome for A and C** (the test-c spectroscopic ages), or reject it. Rejecting would mean that no age is adopted and repair_v9 waits.
2. **B's age.** The carried flagged M1 posterior varies from 1.4 to 5.7 Myr across R_V branches, which is not usable. Options:
   - **(a) borrow** the association's spectroscopic age for B: the A + C test-c curves combined, flagged "not measured for B". This is coherent with the Wright+15 picture and simple. **(Recommended.)**
   - **(b) branches:** carry B at fixed ages 3.2 / 4.0 / 5.0 Myr through repair_v9.
   - **(c) keep the frozen flagged M1 row** as written. Not recommended: it makes B's death count depend on R_V by construction.

   Any choice other than (c) is a decision record before the repair_v9 pre-registration, not a change to this stage.
3. **Then:** pre-register the repair_v9 chain (brief §4.1) with V1's direction stated: C up and A about unchanged, while B depends on decision 2. Get sign-off.

## 9. Files

- Pre-registration: `provenance/wp4v9_age_prereg.json` (`2e4088b`)
- Outcome: `provenance/wp4v9_age_outcome.json`
- Adopted table: `data/processed/wp4_age_posteriors_repair_v9.parquet`
- Injections: `tables/wp4v9_injection_validation.csv`; parts in `data/processed/wp4v9_injection_parts/`
- Real fits: `tables/wp4v9_real_fits.csv`, `tables/wp4v9_faint_end_study.csv`; lnL curves in `provenance/wp4v9_fit_execution.json`
- Post-hoc: `scripts/wp4v9_posthoc_diagnostics.py` → `tables/wp4v9_posthoc_*.csv`
- Log: `data/processed/wp4v9_logs/stage1.log`
