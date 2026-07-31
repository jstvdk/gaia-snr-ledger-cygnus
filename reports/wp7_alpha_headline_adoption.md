# D1 — the headline branch set, adopted; and a structural claim, corrected

*Item **D1** of [pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
Pre-registered in
[wp7_alpha_headline_adoption_prereg.json](../provenance/wp7_alpha_headline_adoption_prereg.json)
**before** the adoption script ran; scored in
[wp7_alpha_headline_adoption_outcome.json](../provenance/wp7_alpha_headline_adoption_outcome.json).
Per-branch table: [wp7_alpha_headline_branch_sets.csv](../tables/wp7_alpha_headline_branch_sets.csv).
Scripts: [wp7_alpha_headline_prereg.py](../scripts/wp7_alpha_headline_prereg.py)
· [wp7_alpha_headline_adopt.py](../scripts/wp7_alpha_headline_adopt.py).*

> **Three of four predictions passed. D1-P1 FAILED and is recorded as failed.
> Its failure was not caused by the decision under test — it exposed a
> pre-existing over-generalization in WP7's prose, which is corrected in §4.**

---

## 1. The decision, adopted

The manuscript's **headline** branch set for N_SN, the explosion timeline and
the last-supernova posterior is **α ∈ {2.0, 2.3}** — 36 branches. α = 2.6 is
reported in the sensitivity table and is never deleted.

| set | branches | N_SN | factor | P(last SN < 100 kyr) |
|---|---:|---|---:|---|
| as published | 54 | 1.93 – 28.74 | 14.9 | 0.183 – 0.889 |
| **headline (adopted)** | **36** | **5.63 – 28.74** | **5.1** | **0.411 – 0.889** |
| sensitivity only (α = 2.6) | 18 | 1.93 – 4.24 | 2.2 | 0.183 – 0.320 |

**The criterion is E1** — the WP5 calibration-window Poisson χ² over 2–8 M☉, a
gate statistic fixed long before WP7 existed. α = 2.6 is the best fit in 1 of 18
cells and has the worst median χ² (10.80, against 6.86 for α = 2.3 and 10.30 for
α = 2.0). **E2, the >8 M☉ census closure, was deliberately not spent**: it
independently agrees (α = 2.6 wins 0 of 18 cells) but is this project's only
out-of-sample validation of the IMF extrapolation, and it is the answer to
WP9's devil's-advocate objection 1. Spending it to select α would make that
answer circular. Restricting further to α = 2.3 alone would require E2 and is
**not** adopted — the standing decision of PROJECT_TRACE §10 item 2 stands.

## 2. A correction the adoption forced: the retained set's median is 13.3, not 9

[wp5_alpha_plausibility.md](wp5_alpha_plausibility.md) §8 and item D1 of the
assessment brief both state the retained set as *"5.63–28.74, median ≈ 9"*.
**That median is wrong.** 8.79 is the ensemble median of the **full 54-branch**
set; the retained set's is **13.29**, which the plausibility run's own
`candidate_branch_sets` block records correctly. The prose read the wrong row.

The consequence matters for how the result is framed: **dropping α = 2.6 does
not merely trim a low tail, it raises the ensemble median by about 50%.** A
decision that narrows a reported range while raising its centre is exactly the
kind of move a referee should interrogate, so the defence has to be explicit:

- the criterion is E1, which existed before the branch grid was evaluated and
  was never chosen with the resulting range in view;
- the excluded branch is reported with its own numbers, not deleted;
- **the baseline branch does not move.** N_SN = 8.43 on PARSEC, R_V = 3.1,
  α = 2.3, coeval, all-explode, before and after (prediction D1-P4, PASS).

**Therefore the manuscript leads with the baseline branch value and the carried
range beside it — never with an ensemble median.** Branches are carried, never
averaged (plan §1.4); the median of a non-probability-weighted branch set is not
a posterior summary and should not be quoted as one. This is now a reporting
obligation attached to the adoption.

## 3. Predictions, scored

| | statement | outcome |
|---|---|---|
| **D1-P1** | zero SNe for any BH threshold ≤ 40 M☉ on every retained branch | **FAIL** — see §4 |
| **D1-P2** | the restriction does not weaken the WP8 pulsar window | **PASS** — retained minimum 0.818 against 0.510 for the full set |
| **D1-P3** | nothing is deleted; `wp7_ledger.csv` unchanged, all 18 α = 2.6 branches present | **PASS** — sha256 identical before and after |
| **D1-P4** | the baseline branch does not move | **PASS** — 8.43 |

Two of these needed per-branch quantities that WP7 stored **only for the
baseline branch** (the black-hole scan and R_SN(t) are baseline-only products).
Rather than weaken the tests to what happened to be on disk, they were
recomputed over all 54 branches with WP7's own engine and the same frozen WP5
draws at 400,000 iterations. That reduced run reproduces the published
2,000,000-iteration N_SN means to **0.22% worst case**, which is the validation
that licenses using it.

## 4. D1-P1 failed, and it caught something real

**The claim as published:** WP7 §2a and the PROJECT_TRACE status line say *"for
any black-hole threshold ≤ 40 M☉ the ledger returns exactly zero on every
branch"*. Measured over all 54 branches, that is **false at the 40 M☉ cut**.

**Why.** Under the star-formation-duration branch a star may be born up to δ/2
*before* the subgroup's fitted age, so the lowest turnoff any star was compared
against is not the turnoff at the fitted age. Over the grid it falls to
**33.9 M☉** at δ = 2 Myr, against **45.4 M☉** on the coeval branches. Twenty-one
of 54 cells then place a small non-zero count below 40 M☉ — worst case
**1.8%** of that branch's N_SN.

**It is not caused by the restriction.** The affected cells are split across
retained and dropped sets in the same proportion (14 retained, 7 dropped, i.e.
exactly the 2:1 branch ratio). The pre-registration's failure clause requires
the structural claim to be **re-derived** before adoption, not the restriction
to be abandoned; that is what §4.1 does.

### 4.1 The claim, re-derived

> **For any black-hole threshold at or below 30 M☉ the ledger returns exactly
> zero supernovae on every one of the 54 branches and in every iteration.** At a
> 40 M☉ threshold it is exactly zero on every coeval branch and at most 1.8% of
> N_SN on the 1–2 Myr formation-window branches.

| quantity | value |
|---|---:|
| largest BH cut that is zero on every branch | **30 M☉** |
| lowest turnoff, all branches | 33.9 M☉ |
| lowest turnoff, coeval branches | 45.4 M☉ |
| lowest progenitor that actually exploded, all branches | 33.9 M☉ |
| lowest progenitor that actually exploded, baseline branch | **52.1 M☉** |

**Prediction L2 is unaffected and stands as PASS.** L2 was written against the
*islands* cut of **25 M☉**, not 40, and the count below 25 M☉ is exactly zero on
every branch and in every one of the iterations run. What was over-general is
only the prose that generalized L2 from 25 to 40 M☉.

**The paper's actual argument is unchanged.** The Sukhbold+2016 / Ertl+2016
islands of implosion sit between roughly 15 and 25 M☉. The corrected floor of
30 M☉ lies above that entire structure, so *"the island pattern is irrelevant at
this age, and the whole budget is conditional on whether very massive stars
explode"* holds exactly as before — now stated at a mass the grid actually
supports rather than one it supports only at δ = 0.

### 4.2 The "> 52 M☉" statements must also be qualified

WP7 §2a and WP9 §2 say every supernova came from a star **above 52 M☉**. That is
exactly true on the **baseline** branch (minimum dead progenitor 52.1 M☉, 0%
below 52) and on coeval branches generally (≤ 6.5% below 52 M☉), but on the
δ = 2 Myr branches up to **19%** of supernovae come from progenitors between
34 and 52 M☉.

**WP9's C3 = 1.000 stands as computed**: C3 requires progenitors above ~30 M☉ to
be envelope-stripped, and the minimum dead progenitor anywhere on the grid is
33.9 M☉ — still above 30, on every branch. What must be corrected is the
*characterisation* "far above the ~30 M☉ stripping threshold": on the widest
formation-window branches the margin is a few solar masses, not a factor. The
sentence the paper can defend is: *every supernova in the ledger, on every
branch, came from a progenitor above 34 M☉ — above the envelope-stripping
threshold on every branch, and above 52 M☉ on the coeval ones.*

## 5. What the manuscript quotes

> Cyg OB2 has produced **N_SN = 8.4** supernovae on the baseline branch (median
> 8, 68% [5, 11]), with a carried range of **5.6–28.7** across the 36 headline
> branches (α ∈ {2.0, 2.3}); the α = 2.6 branches give 1.9–4.2 and are
> disfavoured by the calibration-window mass function, reported in the
> sensitivity table.

with, in the same section, the four obligations accepted at adoption: name E1 as
the criterion, state that E2 was not spent and why, keep α = 2.6's numbers
visible, and never lead with an ensemble median.
