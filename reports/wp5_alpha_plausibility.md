# How plausible is each IMF slope branch?

*Pre-registered in [wp5_alpha_plausibility_prereg.json](../provenance/wp5_alpha_plausibility_prereg.json)
before the statistics were computed. Executed 2026-07-29. Machine-readable
record: [wp5_alpha_plausibility_execution.json](../provenance/wp5_alpha_plausibility_execution.json).*

> **Nothing was removed, regenerated or reweighted.** All 54 branches remain in
> `tables/wp7_ledger.csv` and every published number stands unaltered. This is a
> plausibility measurement for a decision that has **not been taken**. If it is
> ever acted on, that adoption requires its own pre-registration.

---

## 1. Why only α matters

A variance decomposition of log N_SN across the 54 all-explode branches:

| axis | share |
|---|---:|
| **IMF slope α** | **91.9%** |
| R_V | 4.0% |
| isochrone family | 2.0% |
| star-formation duration | 0.2% |

The WP7 spread is one question about α, not four questions about four axes.
Narrowing R_V or the family would buy essentially nothing.

## 2. Two independent lines of evidence

**E1 — the 2–8 M☉ calibration window.** `k` is refitted at every α, so the
normalization is free and cannot absorb a slope error; the Poisson χ² measures
whether the assumed *shape* matches the observed low-mass mass function.
**Cost of using it: none.** It is internal to WP5 and is already a published
gate statistic.

**E2 — the >8 M☉ census closure.** `k` is fitted from 2–8 M☉ counts alone and
the >8 M☉ census never enters the WP5 likelihood, so this is a genuine
out-of-sample test of the extrapolation. **Cost of using it: real** — see §5.

The two probe **different mass ranges**, which is what makes their agreement
informative.

## 3. What they say

| | α = 2.0 | α = 2.3 | α = 2.6 |
|---|---:|---:|---:|
| **E1** median χ² (5 dof) | 10.30 | **6.86** | 10.80 |
| **E1** cells won (of 18) | 7 | **10** | **1** |
| **E2** median \|ratio − 1\| in 68% half-widths | 6.55 | **1.65** | 5.38 |
| **E2** cells won (of 18) | 6 | **12** | **0** |

**α = 2.6 is disfavoured by both lines**, winning 1 of 18 cells on E1 and **none
at all** on E2. **α = 2.0 versus 2.3 is genuinely competitive** on E1, and less
so on E2.

## 4. The striking result: the subgroups disagree, and both lines agree about *how*

| subgroup | E1 prefers | E2 prefers |
|---|---:|---:|
| CygOB2-A | **2.3** | **2.3** |
| CygOB2-B | **2.3** | **2.3** |
| CygOB2-C | **2.0** | **2.0** |

**Two independent mass ranges reach the same per-subgroup verdict, cell by
cell.** A and B want Salpeter; C wants something shallower — and the 2–8 M☉ fit
and the >8 M☉ extrapolation say so *separately*.

This was pre-registered as A2 and it **passed**. It is the fourth independent
appearance of the "CygOB2-C is different" signal:

1. WP6's closure ratio — C at 1.360 against A's 0.865
2. WP6's M2 multiplicity prediction — **failed** because C behaved unlike B
3. WP6's closing α — C at 2.056 against A's 2.336
4. **this diagnostic** — C prefers a shallower slope on both E1 and E2

Four appearances across three different analyses is no longer a curiosity. It is
a property of CygOB2-C that the paper has to state.

**E2 discriminates about 4× more sharply than E1** (prediction A3, passed), which
is expected: a slope error compounds over the lever arm from the calibration
window to the census.

## 5. Candidate branch sets — costed, not adopted

| candidate set | branches | N_SN | factor | what it costs |
|---|---:|---|---:|---|
| **as published** | 54 | 1.93 – 28.74 | 14.9 | nothing — current reported range |
| **drop α = 2.6** | 36 | 5.63 – 28.74 | **5.1** | **nothing** — E1 alone supports it, and E1 is free |
| α = 2.3 only | 18 | 5.63 – 11.16 | 2.0 | **spends E2** — see below |
| α ∈ {2.0, 2.3}, drop R_V = 3.5 | 24 | 8.43 – 28.74 | 3.4 | needs extinction evidence this diagnostic does not provide |

### The cheapest defensible move

**Dropping α = 2.6 costs nothing.** E1 alone supports it — 1 win in 18 cells, and
the worst median χ² — and E1 is a WP5-internal gate statistic that consumes no
validation. E2 independently agrees without needing to be spent. This cuts the
reported spread from a factor of **14.9 to 5.1**, and removes the "maybe only 2
supernovae" tail.

### The move with a price

**Restricting to α = 2.3 would cut the spread to a factor of 2.0** — but E1 alone
does not support it, since α = 2.0 wins 7 of 18 cells. It therefore requires E2,
and **E2 is currently this analysis's only out-of-sample validation of the IMF.**
Spend it to select α and it can no longer be reported as an independent
confirmation. The project gets one or the other, not both.
[PROJECT_TRACE §10c item 2](../PROJECT_TRACE.md) already records the standing
decision that WP7 must not refit α, for exactly this reason.

## 6. What no choice of branch set can fix

**CygOB2-C contributes 0 or about 7 supernovae depending on the isochrone
family**, because its fitted age straddles the age at which the turnoff crosses
the 120 M☉ IMF ceiling. That is independent of α and survives every candidate
set above.

And a standing caveat: E1 constrains the slope at **2–8 M☉**, while N_SN depends
on the slope above **52 M☉**. The link between them is the **single power-law
assumption** from 2 to 120 M☉. This diagnostic assumes it and does not test it.

## 7. Predictions, scored

| | statement | outcome |
|---|---|---|
| **A1** | E1 and E2 agree in direction, both disfavouring 2.6 | **PASS** |
| **A2** | the subgroups disagree, C preferring a shallower slope than A | **PASS** |
| **A3** | E2 discriminates more sharply than E1 | **PASS** — 4.9σ vs 1.25 in own units |

## 8. The decision, left open

The principal investigator is being asked one question: **report N_SN over all
three α branches, or over a subset?**

My reading — offered, not applied — is that **dropping α = 2.6 is free and
worth doing**, and that **restricting further to α = 2.3 is not worth the price
of the project's only out-of-sample check**. That leaves

> **N_SN = 5.6 – 28.7, median ≈ 9**, with CygOB2-C's contribution ranging 0 to 7
> depending on isochrone family.

Nothing in the pipeline has been changed to reflect this.
