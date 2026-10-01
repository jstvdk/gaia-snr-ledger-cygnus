# WP13 — the pooled-versus-resolved ablation (M0 vs M1): the experiment that decides the paper

**Date:** 2026-09-14
**Status:** brief. Nothing in it has been executed against a pre-registration.
Thresholds in §5 are *proposed*; they become binding only when transcribed into
`scripts/wp13_prereg.py` and hashed, **before** any M0 number is read.
**Governing document:** [reports/novelty_prior_art_and_incremental_value_audit_2026-08-12.md](../reports/novelty_prior_art_and_incremental_value_audit_2026-08-12.md) §5–§9 and §12.
**Extension to the plan:** WP13 is not in `paper1_execution_plan.md`, like WP11 and WP12.

---

## 0. Why the project is stuck, in three sentences

The chain is finished, compiled and validated, and the 12 August audit found
that most of what the manuscript emphasises is either already established or a
direct consequence of known physics. The one claim that could carry a paper —
that resolving Cygnus OB2 into subgroups materially changes the inferred death
history because Cyg OB2-C sits below the first-death boundary — is a hypothesis
that has never been tested against the best *one*-population model of the same
stars. Nothing in the repository has changed since the audit was written, so the
project cannot move until that test is run; once it has run, the paper is
publishable under **every** outcome, because each branch of the decision tree in
§6 has a framing already written for it.

## 1. The falsifiable sentence

> Relative to the best single-age, single-normalization model of the same
> *Gaia*-selected stars (**M0**), the subgroup-resolved model (**M1**) changes
> the inferred number of stellar deaths by at least T_N and the probability of an
> event within the last 100 kyr by at least T_P, because Cyg OB2-C lies below
> the first-death boundary while Cyg OB2-A and -B lie above it.

If this sentence is false, the paper says so and becomes a census paper with a
supernova ledger as an application (§6, outcomes 2–4). If it is true, it is the
result.

## 2. Disclosed prior knowledge

Everything below was computed while writing this brief, from tables that already
exist. It is disclosed so that the thresholds in §5 are set *knowing* it, as WP11
did for I3. None of it is the M0 result, because none of it fits a pooled age.

### 2.1 What the existing common-age scan already implies

`tables/wp7_age_sensitivity.csv` forces one age on all three subgroups while
keeping each subgroup's own `k`. Interpolating it:

| model | age (Myr) | N_death | P(last < 100 kyr) |
|---|---:|---:|---:|
| resolved baseline (M1) | A 4.00 / B 4.09 / C 2.52 | **8.43** | **0.552** |
| common age forced to 4.0 | 4.00 | 12.73 | 0.727 |
| common age = k-weighted mean of the three counts-based ages | 3.494 | 6.46 | 0.700 |
| common age that reproduces the resolved count | 3.659 | 8.43 | 0.701 |

**Reading.** A single age of about 3.66 Myr reproduces the resolved *count*.
No single age between 3.25 and 6 Myr reproduces the resolved *recency
probability*: the scan sits at 0.70–0.77 throughout while M1 gives 0.55. The
honest expectation is therefore that the count test (T2) will be marginal or
fail and the time-history tests (T3, T4) are where the difference lives.
Thresholds are set anyway; the expectation is recorded here so it cannot be
claimed as a prediction afterwards.

### 2.2 The structural finding is robust on the authoritative grid

`data/processed/wp5_imf_normalization_repair_v7.parquet`, counts-based
posterior-mean ages over all 18 family × R_V × α cells:

| subgroup | range (Myr) |
|---|---|
| Cyg OB2-A | 3.85–4.07 |
| Cyg OB2-B | 3.30–4.30 |
| Cyg OB2-C | **2.03–3.17** |

C is younger than A on every cell. The "two older, one younger" structure is not
an artefact of one branch. What *is* branch-dependent is whether C's turnoff has
crossed the 120 M☉ ceiling (it has on the 9 PARSEC-and-MIST-3.5 cells, giving C
exactly zero; it has not on MIST 3.0/3.1, giving C 2.5–7.3).

> ⚠️ `tables/wp4_ages_table.md` and `tables/wp4_ages_envelope.md` are
> **pre-repair** products: they list C's upper-main-sequence age as 3.98 Myr,
> which is the age WP4 found before the `repair_v1`–`repair_v5` extinction fixes.
> The repair_v5 parquet gives 2.51 (PARSEC 3.0/3.1), 2.00 (PARSEC 3.5), 3.18
> (MIST 3.0/3.1), 2.00 (MIST 3.5). These two markdown tables are not in
> `wp10_inputs.py`'s forbidden list and should be added to it, or regenerated.

### 2.3 M_prior: Menchiari et al. (2024) Appendix B, reproduced with our relations

Their recipe, verbatim from the appendix: Kroupa (2001) IMF sampled over
0.08–150 M☉, total mass 16,500 (+3,800/−2,800) M☉ from Wright et al. (2015),
1,000 mock populations, "counting how many massive stars have a main sequence
phase lasting less than the cluster age". Their lifetime prescription is **not
stated** in the appendix (it is inherited from Menchiari 2023). Results:
7 ± 2.5 at 3 Myr, 26 ± 5 at 5 Myr.

Expected deaths from the same IMF, mass and mass range, using *this project's*
`TurnoffRelation` (the relation every WP7 number uses), computed analytically
(the 1,000-mock scatter is Poisson around these means):

| mass (M☉) | age (Myr) | PARSEC, ≤150 | PARSEC, ≤120 | MIST, ≤150 | MIST, ≤120 |
|---:|---:|---:|---:|---:|---:|
| 16,500 | 3.0 | 1.50 | 0.14 | 6.23 | 4.92 |
| 16,500 | 5.0 | 18.75 | 17.57 | 18.44 | 17.26 |
| 20,300 | 3.0 | 1.85 | 0.17 | 7.67 | 6.05 |
| 20,300 | 5.0 | 23.06 | 21.62 | 22.68 | 21.24 |

The common age at which our relations give their 7 is 3.65 Myr (PARSEC) or
3.16 Myr (MIST); for their 26 it is 5.9 Myr on either family.

**Reading.** Their 3 Myr number is reproduced by MIST within their own error bar
and *not* by PARSEC, and their 5 Myr number is reproduced by neither. The
turnoff at 3 Myr is 118 M☉ (PARSEC) against 73 M☉ (MIST), so **the lifetime
prescription alone is worth the same factor of four at 3 Myr as the 3-to-5 Myr
age change they emphasise.** Two consequences for the paper: (i) the numerical
proximity of our 8.43 to their 7 ± 2.5 is coincidental, not a validation, and
must not be presented as one; (ii) "a fixed-age estimate is sensitive to the
lifetime prescription at the same level as to the age" is a legitimate, cheap
finding for the sensitivity section. Our association mass is 1.47× theirs
like-for-like, which would scale their number to ~10 at 3 Myr under MIST.

## 3. The model hierarchy

| model | what it is | status |
|---|---|---|
| **M_prior** | Menchiari's recipe with our turnoff relation | pilot done (§2.3); finalise once the lifetime source is obtained from Menchiari (2023), else report as "reproduced to within the lifetime-prescription systematic" |
| **M0** | the same 1,331 labelled members treated as one population: one age, one `k`, per family × R_V × α cell | **to run** |
| **M0-lite** | M0's pooled age with `k_ALL = k_A + k_B + k_C` draw-wise | fallback only; may not be quoted as "the best one-population model" |
| **M1** | the accepted `repair_v7` ledger, unchanged, hash-verified | done |

**M0 definition, step by step.** The pooled label is constructed inside the WP13
script and never written into `wp2_subgroup_labels.parquet`.

1. *Age.* Run the WP4 upper-main-sequence likelihood of `wp4_common` on the
   pooled de-reddened CMD, same f_bin ∈ {0.3, 0.4, 0.5}, R_V ∈ {3.0, 3.1, 3.5},
   family ∈ {PARSEC, MIST}, same membership-probability weighting, same
   N_MIN_INDICATOR, same distance systematic. Record MAP, 68/90 % intervals and
   the full grid posterior; the grid posterior is what T1 needs.
2. *Completeness response.* Run `wp5_injections_agenodes` for the pooled
   footprint at the nine pooled posterior nodes × 6 family/R_V cells. Each node
   run took ~55 s on `repair_v7`, so this is ~1 h of compute.
3. *Normalization.* Run the `wp5_joint_age_fit` machinery on the pooled 2–8 M☉
   counts with the pooled WP4 posterior as prior, producing paired
   `(k_ALL, truth_age_ALL)` draws per cell, exactly as for A, B, C.
4. *Ledger.* Run `wp7_ledger.run_population` on those draws for all 18
   family × R_V × α cells at δ = 0 (and δ ∈ {1, 2} Myr if the M1 comparison is
   to be made on the full 36-branch headline set), 2 × 10⁶ iterations, same seed
   recipe. Pair M0 and M1 by iteration index so the difference has a posterior.

**Held fixed between M0 and M1**, and asserted in the execution record: input
catalogue and membership probabilities; the 61 unlabelled members excluded from
both; extinction and de-reddening; isochrone tables; IMF family, slopes and the
120 M☉ ceiling; calibration window 2–8 M☉; f_bin convention; formation-duration
definition; explodability prescription (all-explode for the primary comparison,
hard cutoff as a secondary row); turnoff relation; RNG recipe.

## 4. Quantities compared

Per cell, paired by iteration, for M0 and M1:

- `N_death` (mean, median, 68 %, 95 %) and the paired difference ΔN with its
  95 % interval;
- `P(last < 100 kyr)` and Δ;
- median `t_last` and the first-death epoch;
- the cumulative death-time distribution F(t), and the sup-distance
  D = sup_t |F_M0 − F_M1|;
- attribution: M1's per-subgroup counts against M0's pro-rata attribution
  `N_M0 × k_sub / Σk` (a pooled model's only possible attribution);
- the WP4 grid marginal likelihood of the pooled fit against the product of the
  three subgroup fits (T1).

## 5. Proposed pass criteria — freeze before M0 is run

| id | question | criterion (baseline PARSEC / 3.1 / 2.3 / δ = 0 unless stated) |
|---|---|---|
| **T1** | do the data support three ages over one? | ln BF(M1 : M0) ≥ 2.3 from the WP4 grid marginal likelihoods, on the baseline cell and on ≥ 4 of the 6 family × R_V cells. The evidence penalises M1's two extra parameters automatically; no ELPD machinery is needed because the WP4 posterior is grid-based. |
| **T2** | does resolution change the count? | \|ΔN\| ≥ 3 **and** \|ΔN\|/N_M1 ≥ 0.10, with the paired 95 % interval of ΔN excluding the ±10 % equivalence band |
| **T3** | does it change the recency inference? | \|ΔP(last < 100 kyr)\| ≥ 0.10, paired 95 % interval excluding zero |
| **T4** | does it change the history, not only the integral? | D ≥ 0.15 **and** D exceeds the 95th percentile of the coeval-injection null (T7a); *or* the first-death epoch shifts by ≥ 0.2 Myr |
| **T5** | is attribution informative? | for at least one subgroup \|N_M1,sub − N_M0 × k_sub/Σk\| ≥ 2 |
| **T6** | is it robust to branches? | the sign of ΔN and of ΔP, and whichever of T2/T3 passed on the baseline, hold on ≥ 80 % of the 36 headline branches and on both families |
| **T7a** | does M1 manufacture a false advantage? | coeval injections (one true age, the three real label sets): T2 or T3 "passes" in < 10 % of realizations |
| **T7b** | does M1 reduce bias? | multi-age injections (true ages A/B/C): M1's bias in N and in P is smaller than M0's, and M1's 68 % interval covers the truth within 58–78 % |

**Verdict rule.** *Material improvement* ⇔ T1 ∧ (T2 ∨ T3 ∨ T4) ∧ T6 ∧ T7a ∧ T7b.
*Supported but immaterial* ⇔ T1 ∧ T7 ∧ ¬(T2 ∨ T3 ∨ T4).
*Equivalent* ⇔ ¬T1 ∧ ¬(T2 ∨ T3 ∨ T4).
*Overfit* ⇔ ¬T7a ∨ ¬T7b, whatever else passes.

**Injection design for T7.** The full-chain injection (mock CMD → WP4 → WP5 →
WP7) is expensive. The required version is ledger-level: draw true ages, apply
the *measured* WP5 posterior widths to produce fitted ages for M0 and M1, run
`run_population` for each, and score bias and coverage; 200 realizations per
design. The full-chain version is optional and deferred to Paper 2 unless T7 is
close to its boundary.

## 6. What each outcome does to the manuscript

Each row is written now so that the narrative is *chosen by the result*, not
after it.

| outcome | title | abstract result sentence | §1.3 "What this paper adds" |
|---|---|---|---|
| **material improvement** | unchanged | "A single-age treatment of the same stars changes N_death by [ΔN] and the probability of an event within 100 kyr by [ΔP], because Cyg OB2-C lies below the first-death boundary; the resolved model is preferred by ln BF = […] and is unbiased in injection." | rewritten around aggregation bias near the first-death boundary; the portable rule of audit §6.3 becomes the closing paragraph |
| **supported but immaterial** | "A *Gaia* DR3 subgroup-resolved census of Cygnus OB2 and its massive-star death ledger" | "The subgroups have distinct ages (ln BF = […]) but the death history is equivalent to a single-age fit within [band]; the ledger is an application." | the census, closure and heterogeneous closing slopes lead; the ledger is Section 6, conditional throughout |
| **equivalent** | "When is a fixed-age supernova count adequate? Cygnus OB2 as a benchmark" | "The resolved and pooled histories agree to [band]; injection maps the domain in which a one-age estimate is unbiased." | methods/benchmark; requires T7 extended over an (age-spread × young-fraction) grid |
| **overfit** | as "supported but immaterial" | the negative finding stated | simplified paper; the subgroup ledger is not the headline |

In all four the cocoon score, the pulsar cross-check and the COSI forecast stay
as applications, exactly as WP12 placed them.

## 7. Manuscript changes required whatever the outcome

From audit §8, checked against the current `main.tex`:

- [ ] The Menchiari paragraph (§1.2) is fair as written; add the §2.3 result
      and state explicitly that the proximity of 8.43 to 7 ± 2.5 is not a
      validation.
- [ ] Add the lifetime-prescription lever to the sensitivity section and to
      Fig. 9 (turnoff at 3 Myr: 118 vs 73 M☉).
- [ ] `slides/make_talk.py` still carries the speaker note "This is a genuine
      discovery about the system" on the mass-floor slide; the audit demotes
      this to an expected consequence of the ages. Fix the note and the slide
      kicker.
- [x] Add the two pre-repair WP4 markdown tables to the forbidden list in
      `wp10_inputs.py` (§2.2). *Done 2026-10-01 under issue #19, with the
      unversioned WP4 posterior and anchor file.*
- [x] Verify which version the quoted 2.25–5.67 Myr envelope comes from; the
      repair_v5 parquet must reproduce it or the macro must be regenerated.
      *Done 2026-10-01 (issue #19): it was the pre-repair run; repair_v5 gives
      2.00–4.01 Myr, now the macros `\ageEnvLo`/`\ageEnvHi`.*
- [ ] Refresh the prior-art sweep immediately before submission (already listed
      in `manuscript/README.md`).

## 8. Work sequence and effort

| step | what | effort |
|---|---|---|
| 1 | transcribe §5 into `scripts/wp13_prereg.py` with input hashes; run it; commit | ½ day |
| 2 | finalise M_prior (obtain Menchiari 2023's lifetime source; otherwise report as in §2.3) | ½ day |
| 3 | M0: pooled WP4 fit, pooled injections, pooled joint fit, ledger | 2–3 days, ~2 h compute |
| 4 | paired comparison, T1–T6 scoring, one figure (M0 vs M1 death-time curves and the ΔN, ΔP posteriors) | 1–2 days |
| 5 | T7 ledger-level injections | 1–2 days |
| 6 | execution record, `reports/wp13_ablation.md`, PROJECT_TRACE §1 row | ½ day |
| 7 | manuscript revision per §6 row selected by the verdict (WP14) | 3–5 days |

About three weeks of focused work to a submittable manuscript, in every branch
of §6.

## 9. Anti-tuning rules

Identical to WP12: no WP13 script alters any upstream artifact, threshold, gate
or branch set; every input is consumed through a hash-verified `frozen()`
resolver; M1's numbers are the stored ones and are reproduced, not recomputed;
if M0 is unwelcome it is reported. The §2 prior knowledge is part of the
pre-registration record. Failed criteria stay failed.
