# WP12 — manuscript revision analysis

*Executed 2026-08-03 against
[tasks/manuscript_reframing_revision_brief.md](../tasks/manuscript_reframing_revision_brief.md).
Pre-registration:
[provenance/wp12_revision_prereg.json](../provenance/wp12_revision_prereg.json),
written and hashed **before** any WP12 result existed.*

---

## 1. What WP12 is

WP12 adds no observation and refits nothing upstream. Every input is a frozen
`repair_v7`-chain product, consumed read-only through `wp12_common.frozen()`,
which verifies each file against the SHA-256 recorded in the preregistration and
raises if it moved. What WP12 adds is presentation-grade structure the accepted
chain always contained but never exposed, plus the two sensitivities the brief
required and one attempted analysis that failed its own gate.

**Gate summary: R1–R5 and R7 PASS, R6 FAIL.** R6 was the declared gate on the
neighbouring-association budget, and its failure has been honoured: the budget
is not reported as a number and the manuscript scopes the cavity-source claim
out instead.

| id | statement | outcome |
|---|---|---|
| **R1** | fewer than half of the 36 headline branches rest entirely on passing WP5 cells | **PASS** — 15/36 (*not blind; disclosed as prior knowledge*) |
| **R2** | the α split survives restriction to all-subgroup-pass branches | **PASS** — 3/3 at α=2.0 above 0.5, 0/12 at α=2.3 |
| **R3** | the three subgroups do not share a closing slope (spread > 0.15) | **PASS** — spread 0.322 at the baseline family and R_V |
| **R4** | the mixed-slope ledger stays inside the carried headline N_SN range | **PASS** — 5.63–15.19 inside 5.63–28.74 |
| **R5** | C3 is not robust to its mapping (< 0.95 at a 60 M☉ threshold) | **PASS** — C3 falls to 0.68–0.96 |
| **R6** | the coarse neighbour estimator reproduces our own ledger to within 10× | **FAIL** — returns 0 at Knödlseder's 2.5 Myr |
| **R7** | Cyg OB2 does not dominate the recent budget (f < 0.9) | PASS, but computed from an estimator R6 rejected — **not reportable** |

---

## 2. WP12.1 — the residual-gate landscape

`scripts/wp12_gate_landscape.py` ·
[execution record](../provenance/wp12_gate_landscape_execution.json)

The manuscript reported "40 of 54 mass-function fits pass" beside a headline
branch set of 36, which invites the reading that the headline set is validated.
It is not: the headline set is selected on α and on nothing else.

**Three statements that must not be merged:**

| statement | value |
|---|---|
| the baseline branch passes in all three subgroups | yes |
| subgroup fit cells passing the residual gate | **40 / 54** |
| headline cells passing | **27 / 36** |
| headline family–R_V–α combinations passing in all three subgroups | **5 / 12** |
| headline ledger branches built entirely on passing cells | **15 / 36** |

**The failure map, on `repair_v7`** — the chain every downstream number was
computed on:

| axis | distribution |
|---|---|
| subgroup | A 4, B 7, C 3 |
| family | PARSEC 4, MIST 10 |
| R_V | 4 at 3.0, 5 at 3.1, 5 at 3.5 |
| α | 7 at 2.0, 2 at 2.3, 5 at 2.6 |

The manuscript's previous claim that failures concentrate in the `R_V = 3.5` and
`α = 2.6` corners is **false** and has been replaced by this distribution. The
dominant axes are isochrone family and, on the headline arm, the shallow slope.

**A version discrepancy worth recording.** The *accepted* WP5 gate record is
`repair_v6`; the chain the ledger consumed is `repair_v7`. Both give 40/54, but
they disagree about *which* 14 cells fail (v6: A 5, B 5, C 4; PARSEC 5, MIST 9).
`repair_v7` is authoritative for the paper because it produced every downstream
number, and the execution record carries both breakdowns.

**The α split survives the strict restriction.** Among the 15 all-subgroup-pass
headline branches, 3/3 at α = 2.0 stay above 0.5 (0.685–0.721) and 0/12 at
α = 2.3 rise above it. The qualitative conclusion is robust; the branch
representation is badly imbalanced, and the "18/18 versus 0/18" phrasing uses
branches that include failed cells.

---

## 3. WP12.2 — subgroup closure, closing slopes, mixed slopes

`scripts/wp12_closure_slopes.py` ·
[execution record](../provenance/wp12_closure_slopes_execution.json)

**Two aggregations, previously conflated.** At α = 2.3 the median of the 18
subgroup closure ratios is **1.067**; the star-weighted ratio of summed observed
to summed predicted counts is **1.154**. A median of ratios is not the ratio of
sums, and the subgroups carry very different numbers of observed massive stars.
The manuscript now quotes both and names which is which.

**Grid-median closure by subgroup at α = 2.3:** A **0.95**, B **1.07**,
C **1.39** — reproducing the brief's numbers exactly.

**Closing slopes (grid medians):** A **2.34**, B **2.25**, C **2.06**,
association aggregate **2.20**. All 18 subgroup cells close *inside* the carried
[2.0, 2.6] range, so no extrapolation is involved. A and B want Salpeter; C
wants distinctly shallower. The heading "The census closes at Salpeter" has been
replaced by "Association-wide closure is near Salpeter, but the subgroups
differ".

**The mixed-slope ledger.** Each subgroup was assigned the carried slope its own
closure prefers — uniformly A 2.3, B 2.3, C 2.0 on every family and R_V — and
the ledger re-run on the frozen posterior draws at those slopes. The engine is
WP7's own `run_population`, unmodified, and it reproduces the stored
single-slope WP7 numbers to **0.2 %** in N_SN and 0.002 absolute in
P(last SN < 100 kyr), well inside the declared 3 % and 0.02 tolerances.

Result: **5.63–15.19 supernovae**, inside the carried headline range (R4 PASS).
But the effect splits sharply on whether Cyg OB2-C's turnoff has crossed the
120 M☉ IMF ceiling:

- on the **9 branches where C is above the ceiling** (PARSEC at every R_V, MIST
  at R_V = 3.5) the change is at most **1 %**;
- on the **9 branches where C contributes** (MIST at R_V = 3.0 and 3.1, where C
  supplies 2.5–7.3 SNe) the change reaches **44 %**.

So the subgroup heterogeneity is not presentational: where it acts it is worth
as much as an R_V step.

### Defects found and fixed inside WP12

Both were caught by disagreement with an independently computed quantity, not by
inspection. Neither affects any upstream number, and both are recorded in the
execution JSON.

| id | defect | how it was caught |
|---|---|---|
| **D1** | `closing_slope()` assumed the closure ratio *decreases* with α. It increases — k is refitted to the same 2–8 M☉ counts at every slope and over-compensates. The bug returned 2.0 for all three subgroups and would have made **R3 read FAIL**. | disagreed with `wp10_numbers.py`'s independent `closingAlpha = 2.25` |
| **D2** | the star-weighted association aggregate was computed as summed *predicted* over summed *observed*. WP6 defines `closure_ratio` the other way round (`wp6_closure_test.py:237`). | the aggregate ran *opposite* to all three subgroup curves, which is impossible for a weighted combination of them |

---

## 4. WP12.3 / 12.4 — the cocoon quantity, relabelled and made conditional

`scripts/wp12_scenario_score.py` ·
[execution record](../provenance/wp12_scenario_score_execution.json)

**Renamed.** `P_verdict` → **conditional scenario-availability score `S`**. The
arithmetic `S = C1·C3·C4` is unchanged and the stored WP9 numbers are
reproduced, not replaced. What changed is the label, because only one of the
three factors is a sampled probability:

| term | what it actually is |
|---|---|
| `C1` | a sampled Monte-Carlo probability |
| `C3` | a **deterministic model mapping** — an indicator that the progenitor exceeds 30 M☉ |
| `C4` | an **upper bound**, held fixed at 0.854 |

Independence of the three is asserted, not demonstrated, and the energy
condition is deliberately not multiplied in.

**Engine validation:** worst |ΔC1| = 0.0018 against the stored 2×10⁶-iteration
WP9 table (tolerance 0.01) and ΔC3 = 0 exactly.

### C4 sensitivity (WP12.3)

The manuscript's sentence "C4 would have to fall below 0.60 to move **any**
α = 2.0 branch under 0.5" conflated two thresholds. Corrected:

- **C4 < 0.720** moves *at least one* α = 2.0 branch below 0.5;
- **C4 < 0.581** moves *every* α = 2.0 branch below 0.5.

The quoted 0.60 was close to the *all-branches* value, not the *any-branch* one.
In the other direction, **C4 would have to exceed 0.901** before even the
highest α = 2.3 branch rose above 0.5 — so the α split is not an artefact of the
adopted C4, but neither is it immune to it. `tables/wp12_c4_scan.csv` gives the
full scan and Fig. 10 plots it.

### C3 as a mapping (WP12.4)

C3 = 1.000 exactly because every sampled progenitor exceeds the 30 M☉ threshold
— a fact about the ledger's mass floor (33.9 M☉ at its lowest), not evidence
that every explosion would have been observed as type Ib/c. Alternative
mappings, all bounded and all reported as a bracket:

| mapping | C3 range over headline branches |
|---|---|
| step at 30 M☉ (baseline, unchanged) | 1.00 |
| step at 40 M☉ | 0.98–1.00 |
| step at 60 M☉ (explicitly pessimistic) | **0.68–0.96** |
| logistic, M₀ ∈ {30, 45, 60}, ΔM ∈ {5, 10} | see `tables/wp12_c3_subtype.csv` |

Under the pessimistic 60 M☉ mapping S spans 0.304–0.691 and 18/18 α = 2.0
branches still clear 0.5 — the qualitative split survives, but **C3 = 1 does not
survive as a statement about nature** (R5 PASS). The manuscript now cites
primary sources for the dependence of envelope stripping on wind mass loss
([Vink et al. 2001](https://ui.adsabs.harvard.edu/abs/2001A%26A...369..574V/abstract);
[Björklund et al. 2021](https://doi.org/10.1051/0004-6361/202038384);
[Smith 2014](https://doi.org/10.1146/annurev-astro-081913-040025)) and on
binarity (Sana+2012, de Mink+2014).

---

## 5. WP12.5 — the neighbouring-association budget: **attempted, gate failed**

`scripts/wp12_neighbour_budget.py` ·
[execution record](../provenance/wp12_neighbour_budget_execution.json)

The original plan's WP8.5 required a coarse literature-based supernova budget
for the neighbouring Cygnus populations. WP8 declined to claim one. The brief
required either completing it or scoping the claim out. **We attempted it and
the pre-registered gate rejected it, so the claim is scoped out.**

**Source (primary, single tabulation):** Martin et al. 2010, A&A 511, A86,
Table 1, reproduced there from Knödlseder et al. 2002, A&A 390, 945. Each row
gives an *observed* massive-star count in a stated ZAMS mass interval, with a
distance and an age, for 13 Cygnus populations.

**Method:** normalize a power law to the observed count truncated at the
turnoff; integrate above the turnoff for cumulative deaths; difference over
100 kyr for the recent rate. The turnoff relation is *this project's own*, so
the coarse estimator and the measured ledger differ in inputs, not in stellar
physics.

**The gate (R6) and its failure.** The estimator was applied to Martin+2010's
own Cyg OB2 row (120 stars in [20, 120] M☉ at 2.5 Myr) and required to land
within a factor of 10 of our measured baseline `N_SN = 8.43`.

| | result |
|---|---|
| at Knödlseder's adopted **2.5 Myr** | PARSEC **0.00** (turnoff 282 M☉ > the 120 M☉ IMF ceiling — nothing has died); MIST 0.99–2.01 |
| ratio to measured | **0.00–0.24** → outside [0.1, 10] → **R6 FAIL** |
| *diagnostic, not part of R6*: at **our measured 4.0 Myr** | 17.8–32.8, ratio **2.1–3.9** — within a factor of 10, but biased high by 2–4× |

Both failures have the same cause: **the estimator is dominated by the assumed
age**, and heterogeneous 2002-era literature ages for the neighbours are exactly
what no homogeneous census yet supplies. A 1.5 Myr error in Cyg OB2's assumed
age moves its coarse budget from zero to tens.

**Consequence, per the preregistration.** The budget is not reported as a
number. `f_OB2` was computed (0.33 on the wide cavity set at baseline) and is
**not quoted in the manuscript**, because it derives from an estimator its own
gate rejected. The paper states instead that it estimates only the availability
of an event *from Cygnus OB2* and cannot quantify whether Cygnus OB2 was the
source of an event elsewhere in the cavity. The attempt, its inputs and its
failed gate are retained in the released records so the decision can be audited.

---

## 6. Outputs

**Scripts** — [wp12_prereg.py](../scripts/wp12_prereg.py) ·
[wp12_common.py](../scripts/wp12_common.py) ·
[wp12_gate_landscape.py](../scripts/wp12_gate_landscape.py) ·
[wp12_closure_slopes.py](../scripts/wp12_closure_slopes.py) ·
[wp12_scenario_score.py](../scripts/wp12_scenario_score.py) ·
[wp12_neighbour_budget.py](../scripts/wp12_neighbour_budget.py) ·
[wp12_tables.py](../scripts/wp12_tables.py) ·
[wp12_figures.py](../scripts/wp12_figures.py)

**Tables** — `wp12_wp5_gate_map.csv` · `wp12_combination_gate.csv` ·
`wp12_branch_gate_table.csv` · `wp12_appendix_a_full.csv` ·
`wp12_closure_by_alpha.csv` · `wp12_closing_slopes.csv` ·
`wp12_mixed_slope_ledger.csv` · `wp12_c4_scan.csv` · `wp12_c3_subtype.csv` ·
`wp12_scenario_score.csv` · `wp12_neighbour_budget.csv` ·
`wp12_cavity_share.csv`

**Figures** — `figures/paper/fig01..fig11`, replacing the six-figure set.

**Manuscript** — `manuscript/tables_generated.tex` (Tables 1–3 and Appendix A,
generated) and a rewritten `manuscript/main.tex`.

**Nothing upstream was modified.** The frozen-input hash check in
`wp12_common.frozen()` enforces this mechanically rather than by discipline.
