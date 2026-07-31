# B2 — the two CygOB2-B ages, reconciled

*Written 2026-07-30 to discharge item **B2** of
[pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
Diagnostic only: nothing was refitted and no stored number moved. Evidence:
[wp4_wp5_age_reconciliation_execution.json](../provenance/wp4_wp5_age_reconciliation_execution.json)
· [wp4_wp5_age_reconciliation.csv](../tables/wp4_wp5_age_reconciliation.csv)
· [wp4_wp5_age_reconciliation.py](../scripts/wp4_wp5_age_reconciliation.py).*

---

## 1. The complaint

| where | A | B | C |
|---|---:|---:|---:|
| WP4 ages as adopted ([PROJECT_TRACE §6](../PROJECT_TRACE.md)) | 3.981 | **3.548** | 2.512 |
| obligation O1 ([PROJECT_TRACE §10 item 6](../PROJECT_TRACE.md)) | 4.00 | **4.07** | 2.52 |

Same document, two sets, no stated distinction. B differs by 0.5 Myr — and at
Cyg OB2's age half a megayear is not a rounding difference.

## 2. Both are real, and neither is stale

**3.548 Myr is the WP4 upper-main-sequence MAP** — a photometric age, from
isochrone fitting to B's de-reddened CMD on the baseline branch (PARSEC,
R_V = 3.1, f_bin = 0.40, no distance-modulus offset), on the repair_v5
extinction.

**4.07 Myr is the WP5 fitted truth-age posterior mean** — a *counts-based* age.
WP5's joint age–k fit takes the WP4 posterior as a prior over nine unsnapped
quantile nodes and reweights those nodes by the Poisson likelihood of the
observed 2–8 M☉ counts, with `k` integrated out under a Jeffreys prior
(`repair_v4` machinery, `repair_v6` interpolation). It is the age at which the
injection truth model reproduces the observed mass function.

The identification is exact, not approximate. Fitted posterior means on the
baseline branch:

| version | A | B | C |
|---|---:|---:|---:|
| repair_v6 | 3.997 | **4.066** | 2.516 |
| repair_v7 (accepted chain) | 4.000 | **4.086** | 2.517 |

O1's "A 4.00, B 4.07, C 2.52" is the repair_v6 row, rounded. Nothing is stale.

## 3. Which one is load-bearing

**The counts-based age.** [wp7_ledger.py](../scripts/wp7_ledger.py) and
[wp9_verdict.py](../scripts/wp9_verdict.py) both draw paired `(k, truth_age)`
samples from `truth_age_draws__*` in the WP5 posterior archive — deliberately
paired, because `k` and the age are correlated through the WP5 fit. The
upper-MS MAP enters the ledger only as the *prior* of that posterior. So every
downstream number — N_SN, R_SN(t), P(last SN < 100 kyr), the WP9 verdict — rests
on ~4.09 Myr for B, never on 3.548.

## 4. What the difference is worth, in supernovae

N_SN is essentially the IMF integral above the turnoff, so an age difference is
a turnoff difference is a supernova-count difference:

| subgroup | UMS MAP | fitted | shift | M_turnoff(UMS) | M_turnoff(fitted) | SN count ratio |
|---|---:|---:|---:|---:|---:|---:|
| CygOB2-A | 3.981 | 4.000 | +0.019 | 58.3 M☉ | 57.9 M☉ | **1.02** |
| **CygOB2-B** | **3.548** | **4.086** | **+0.538** | 73.7 M☉ | 55.8 M☉ | **1.92** |
| CygOB2-C | 2.512 | 2.517 | +0.005 | 281.7 M☉ | 278.7 M☉ | n/a (turnoff above the 120 M☉ ceiling — nothing has died) |

**Quoting B at 3.548 Myr would halve its supernova contribution.** B supplies
4.26 of the baseline ledger's 8.43 SNe, so the wrong age would take the
association total to roughly 6.4 — a 24% error in the paper's headline number,
produced entirely by a labelling ambiguity. This is why B2 was correctly
classified as blocking.

Across the whole 18-cell baseline-family grid the same pattern holds and is
B-specific: median age shift A +0.036, **B +0.448**, C +0.009 Myr.

## 5. A finding that fell out of the check: B's fitted age rails

B's counts-based posterior does not sit comfortably inside its prior support —
it piles up against the top of it. On the baseline branch **53.2% of the
posterior weight is on the topmost of nine nodes and 83.7% on the top two**,
against 0.166 and 0.124 for A and C. Grid medians for the top-node weight:
**B 0.435, A 0.133, C 0.134.**

A posterior that concentrates on the edge of its prior grid is reporting a
**bound, not a measurement**. The direction of that bound is known and one-sided:
older → lower turnoff → *more* supernovae. If B's true age were the top node
itself (4.26 Myr), its supernova count would be **+15%** on the value carried.

Three consequences:

1. **B's contribution to N_SN is a lower bound at the quoted precision**, and
   the paper must say so. It is not symmetric with the branch spread.
2. It is the same direction as issue #9 (B's PMS indicator wants B older than
   its upper-MS indicator) and the same direction as the repair_v5 extinction
   correction. Three independent lines now point the same way, which is why the
   *structural* result — two older subgroups plus one younger — is more robust
   than B's central value.
3. It is a real limitation of the joint-fit design, stated here rather than
   discovered by a referee: the truth-age node set is WP4's nine posterior
   quantiles, so the counts can never pull the age outside WP4's own support,
   however strongly they want to. Widening the node support is a Paper 2 / DR4
   item, not a Paper 1 re-run — B's photometric age is exactly the quantity
   Gaia DR4 astrometry and XP spectra improve (obligation O3).

## 6. Resolution — what the paper quotes

**One age per subgroup, with its definition attached, and it is the
counts-based one:**

> CygOB2-A 4.00, CygOB2-B 4.09, CygOB2-C 2.52 Myr *(WP5 joint age–k truth-age
> posterior means, baseline branch, repair_v7; the underlying WP4 upper-MS CMD
> MAPs are 3.98 / 3.55 / 2.51 Myr and enter as the prior)*.

The upper-MS MAPs appear once, in the methods section, labelled as the prior.
They are never quoted as "the age" of a subgroup. Where the star-formation
history is stated as a result (obligation O1), the counts-based set is the one
used, with B's railing carried as a one-sided systematic.

The wider honest envelope from WP4 — 2.25–5.67 Myr across both retained
indicators — is unchanged by any of this and remains the age uncertainty the
paper reports; §4's table is the *central-value* reconciliation, not a claim of
0.02 Myr precision.
