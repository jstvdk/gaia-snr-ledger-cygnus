# WP6 — Massive-star census closure and runaway correction

*2026-07-28. Consumes the accepted WP5 `repair_v6`.
Evidence: [wp6_closure_test_execution.json](provenance/wp6_closure_test_execution.json) ·
[wp6_closure_attribution_execution.json](provenance/wp6_closure_attribution_execution.json) ·
[wp6_runaways_execution.json](provenance/wp6_runaways_execution.json) ·
[wp6_runaway_crossmatch_execution.json](provenance/wp6_runaway_crossmatch_execution.json) ·
[wp6_ledger_execution.json](provenance/wp6_ledger_execution.json).*

**Headline: the census very nearly closes at the Salpeter branch.** Observed
living massive stars exceed the WP5 prediction by **11%** at α = 2.3 (grid
median), and the census closes at **α ≈ 2.23** — within the carried branch grid
and close to Salpeter's 2.35. Because α is a Class E branch carried in parallel,
the correct reading is that **WP6 discriminates between the branches**; it does
not show that WP5 was wrong.

> **These numbers replace an earlier version of this document.** It reported a
> 45% excess closing at α ≈ 2.07 and read that as evidence for a
> shallower-than-Salpeter IMF. **That reading was mostly an estimator bug**
> (issue #17, §2): the forward integral truncated at 8 M☉ while the observed
> side counted every member's `P(M > 8)` regardless of true mass, so stars
> scattered up across the threshold were counted but never predicted. Fixing the
> bound absorbed about three quarters of the excess. The withdrawn figures —
> **1.087 / 1.448 / 1.706, grid median 1.444, closing α = 2.070** — should not
> be quoted.
>
> **Issue #15 is now resolved** (§6): unmodelled massive-star multiplicity above
> 8 M☉ absorbs only **3.7%** of the remaining excess, so the residual IMF reading
> survives its strongest pre-registered instrumental challenge. It is reported
> with that correction applied and the remainder carried as a systematic. The
> disfavouring of α = 2.6 stands. **Caveat:** multiplicity *below* 8 M☉ is
> untested and would require a full `repair_v7` chain re-run.

---

## 1. What was fixed before the test could run

**Issue #3 — the estimator.** WP6's closure step assumed bright-mass
completeness ≈ 1. It is not. The binding estimator is now
[CUTS §16](CUTS_AND_THRESHOLDS.md): a forward comparison,
`k · ∫[M_floor, M_turnoff] dM M^(−α) R(observed | M)`, applied **per subgroup**.
Dividing an observed count by a scalar completeness is forbidden — the response
also scatters masses across the threshold, and in 6 of 54 cells the scalar
correction has the wrong sign.

**Issue #14 — the turnoff was a step function of age.** Reading it as
`load_isochrone_between_ages(...)['Mini'].max()` returns the *older* bracket's
turnoff across each native interval: 48.0 M☉ at both 4.00 and 4.20 Myr. Since
N_SN is essentially the IMF integral above the turnoff, this would have made the
paper's headline number a step function of age. Fixed by log-log interpolation;
it also removed a spurious 22% PARSEC/MIST disagreement at 4 Myr (48.0 vs 58.6 →
57.9 vs 58.7).

**Step 0a — the response ceiling.** WP5's response stopped at 18 M☉ while the
closure window runs to the turnoff (60–120 M☉). The two sides were not
comparable. 162 extension injections were run rather than assuming flatness
(decision: [wp6_mass_extension_decision.json](provenance/wp6_mass_extension_decision.json)).

**That decision was vindicated by the measurement.** Recovery is **not** flat
above 18 M☉ — it falls from ~0.80 at 8–20 M☉ to **~0.57 near 48 M☉** before
recovering to ~0.75 at 115 M☉. Assuming flatness would have *understated* the
predicted observable count and made the excess larger still.

**V1 passed**: all five accepted `repair_v6` artifacts are byte-identical after
the extension ([wp6_verify_wp5_untouched.json](provenance/wp6_verify_wp5_untouched.json)).
`MASS_GRID` was deliberately left frozen.

## 2. Issue #17 — the integral truncated at the threshold it smears across

*Pre-registered before the grid run:
[wp6_closure_floor_prereg.json](provenance/wp6_closure_floor_prereg.json) ·
scored: [wp6_closure_floor_score_execution.json](provenance/wp6_closure_floor_score_execution.json).*

The predicted side integrated from 8 M☉. The observed side counts every member's
`P(M > 8)` **whatever that member's true mass is**. Those are not the same
quantity: a star with true mass 7 M☉ whose estimated mass lands above 8 is
counted as observed and was never predicted.

The omitted up-scatter is not small:

| true mass | R(estimated > 8 M☉ \| M) |
|---:|---:|
| 7.0 | 0.231 |
| 6.0 | 0.089 |
| 5.0 | 0.029 |
| 4.0 | 0.004 |

**The specification was already right and the code contradicted it.** CUTS §16.2
writes the integral with *no* lower limit, and its stated reason for forbidding a
scalar completeness is that "the response also scatters mass estimates across the
threshold" — it even reports six cells where net up-scatter exceeds the recovery
loss. This is **issue #3's error committed a second time, at the other end of the
same integral**.

Why it survived review: the *upper* limit absorbed all the scrutiny — it is
physically meaningful, it drove the mass extension, and it produced issue #14 —
while the lower limit read as a restatement of the 8 M☉ supernova threshold
rather than as an integration bound. And the bias has the **same sign** as the
excess WP6 was already reporting, so nothing downstream contradicted it.

**The floor is set by convergence, not by its effect** — and the convergence scan
is reported in full precisely because that claim is checkable:

| floor | grid median ratio |
|---:|---:|
| 8.0 | 1.444 |
| 7.0 | 1.304 |
| 6.0 | 1.203 |
| 5.0 | 1.140 |
| **4.0 (adopted)** | **1.105** |
| 3.0 | 1.089 |
| 2.0 | 1.079 |

**Regression check:** run at floor 8.0, the refactored estimator reproduces the
withdrawn numbers exactly, so the bound is the only behaviour change. No new
injections were needed — the frozen `MASS_GRID` already spans 0.5–18 M☉.

### Scored against the pre-registration

| | prediction | verdict |
|---|---|---|
| **F1** | the ratio falls in all 54 cells | **PASS** — 54/54 |
| **F2** | the reduction is larger at steeper α | **PASS** |
| **F3** | subgroup spread under 5 points | **FAIL** — 6.0 points |
| **F4** | grid median falls below 1.25 | **PASS** — 1.105 |

**F3 is recorded as failed, not reinterpreted.** What the data show is that its
*substantive* claim holds: across all 54 cells the mean reduction is **A −18.3%,
C −18.5%** despite turnoffs of 59.7 and 120.0 M☉, so the effect is flatly
uncorrelated with the turnoff. The spread comes from **B alone**, whose effective
completeness [CUTS §16.4](CUTS_AND_THRESHOLDS.md) already documents as differing
from A's and C's by ~0.12 — which is why the response is applied per subgroup in
the first place. The 5-point threshold was set too tight for a quantity whose
per-subgroup spread was already on the record.

Convergence is slower at grid level than the single probe node suggested (1.4%
from 4.0 → 3.0, a further 1.0% to 2.0). **The floor was not moved in response** —
it was pre-declared — and the residual ~2.4% is carried as a systematic.

## 3. Counting convention

`N_obs = Σ p_membership · P(M > 8)` per branch, not a threshold on point-estimate
masses. Measured Eddington factor: **0.936–1.083** across the grid — *smaller
than anticipated, and it under-counts in A and B rather than over-counting
everywhere*, which is the opposite of what a naive argument predicts.

**62 of 252 countable anchors are absent from the member sample** — 30 O stars,
13 B, 6 Wolf-Rayets, 21 supergiants. These are bright stars that failed the WP2
quality filter, which is exactly the loss channel the response models. They are
therefore an **independent check on the response, not a correction added to the
census** — folding them in would double-correct. After a 1° footprint cut (19
belong to other associations) and marking 13 with no spectral type as *unknown*
rather than sub-threshold: **27 orphans above 8 M☉**.

## 4. The closure test

Baseline PARSEC, R_V = 3.1, α = 2.3:

| subgroup | turnoff | predicted true | predicted observed | observed | **ratio** | 68% |
|---|---:|---:|---:|---:|---:|---|
| CygOB2-A | 59.7 | 82.3 | 66.5 | 59.4 | **0.894** | [0.84, 0.95] |
| CygOB2-B | 71.1 | 79.3 | 72.0 | 79.6 | **1.106** | [1.03, 1.19] |
| CygOB2-C | 120.0 | 94.6 | 76.5 | 107.4 | **1.405** | [1.31, 1.51] |

Across all 54 cells: 0.582–2.273, median **1.105**; **6 of 54** consistent with
unity at 68% (was 2 before issue #17).

Note the two prediction columns have deliberately different lower bounds:
`predicted true` counts stars **physically** above 8 M☉ and starts there;
`predicted observed` counts stars **estimated** above 8 M☉ and starts at the
4 M☉ convergence floor. A is now *below* unity — the corrected estimator does
not simply shrink the excess everywhere, it changes its sign in one subgroup.

## 5. Attribution

**Two of the three channels the plan lists cannot produce this sign.**
Extinction-hidden stars and escaped runaways both *remove* stars from view and
can only push the ratio down. Only a genuine IMF deviation — or a
mis-specification the plan did not anticipate, because it assumed a deficit —
can produce an excess.

**The largest single contributor turned out to be the third option.** Issue #17
was a mis-specification of exactly that unanticipated kind, and it accounted for
roughly three quarters of what the earlier version of this document attributed
to the IMF. The residual is much smaller and, in CygOB2-A, has changed sign.

**The slope dominates everything else:**

| α | grid median ratio |
|---|---:|
| 2.0 | **0.722** |
| 2.3 | 1.105 |
| 2.6 | **1.677** |

Slope that would close the census:

| subgroup | closing α | range |
|---|---:|---|
| CygOB2-A | 2.314 | [2.248, 2.379] |
| CygOB2-B | 2.230 | [2.198, 2.239] |
| CygOB2-C | 2.031 | [2.000, 2.073] |
| **grid median** | **2.230** | 17/18 cells inside the carried branch grid |

The consistency across subgroups argues this is structure, not noise. **All
three closing slopes now sit inside the carried grid** and within 0.3 of
Salpeter, where before issue #17 they straddled its lower edge (12/18 cells, and
C extrapolating to 1.92).

**This is an out-of-sample test.** `k` is fitted from 2–8 M☉ counts alone; the
> 8 M☉ census never enters the WP5 likelihood.

Alternatives tested rather than assumed away:

| | verdict | evidence |
|---|---|---|
| A1 mass scale | not supported | injected stars recovered without bias at the top (27.0 → 25.9; 100.0 → 100.2) |
| A2 response shape | works *against* the excess | recovery falls to ~0.57 at 48 M☉; flatness would enlarge the excess |
| A3 turnoff | a real lever, reported | ratio rises with turnoff across subgroups; issue #14 removed the step-function error |
| A4 contamination | too small | membership-weighted vs hard-cut counts differ by a few percent, not 11% |
| A5 integration floor | **the dominant contributor** | issue #17 (§2): grid median 1.444 → 1.105, closing α 2.070 → 2.230 |
| A6 multiplicity | under test | issue #15 (§6), pre-registered |

## 6. Multiplicity — the leading instrumental alternative

*Pre-registered before running:
[wp6_multiplicity_prereg.json](provenance/wp6_multiplicity_prereg.json).*

The truth model injects a constant f_bin = 0.40 at every mass. Massive stars are
not like that:

| source | measurement |
|---|---|
| Sana+2012 | close-binary fraction **0.70** for 15–60 M☉; mass-ratio exponent κ = −0.1 ± 0.6 |
| Duchêne & Kraus 2013 | multiplicity **>90%** for O types vs ~50% solar |
| Caballero-Nieves+2020 | **Cyg OB2 itself**: 47% of 74 O/early-B stars have a resolved companion at 0.08–10″; 48/74 multiple with spectroscopic binaries included |

At 1.62 kpc Gaia resolves ~1000 AU, so essentially every spectroscopic binary
and most visual pairs blend into one source. Adopted model: **0.40 at ≤8 M☉,
rising log-linearly to 0.70 at ≥16 M☉**, both endpoints measured, nothing tuned.
κ already supports the Uniform(0.1, 1) mass ratio in the code, so **q is not
changed**.

Only the **truth** side changes. The estimator keeps assuming 0.40 — that
mismatch *is* the bias under test.

**Paired design.** The published M ≥ 8 response is assembled from two files
generated on different mass grids, so a grid-restricted re-run draws a different
RNG stream even at the same seed; comparing the treatment straight against it
would mix multiplicity with a change of Monte Carlo realization. Each node is
therefore injected **twice on the identical grid** — control at the frozen 0.40,
treatment at f_bin(M) — and the per-star threshold consumes the same single
`rng.random(n_injected)` draw in both arms, so donor, extinction, photometric and
QMC realizations are bit-identical. The arms differ in exactly one thing: which
stars got a companion.

**V4 passed** before any of it ran: with `truth_binary_fraction=None` the
regenerated response is byte-identical to the published one, so the code change
cannot have perturbed anything already accepted.

### Result

*Scored: [wp6_multiplicity_closure_execution.json](provenance/wp6_multiplicity_closure_execution.json)
· injections: [wp6_multiplicity_injections_execution.json](provenance/wp6_multiplicity_injections_execution.json)
(162 paired nodes, 3.3 h, realized f_bin: control 0.396 ± 0.005, treatment 0.611 ± 0.011).*

| subgroup | control | treatment | drop | excess absorbed |
|---|---:|---:|---:|---:|
| CygOB2-A | 0.901 | 0.900 | 0.001 | — (below unity) |
| CygOB2-B | 1.106 | 1.101 | 0.005 | 4% |
| CygOB2-C | 1.416 | 1.413 | 0.004 | 1% |
| **grid median** | **1.103** | **1.099** | **0.004** | **3.7%** |

**Multiplicity above 8 M☉ is not the explanation.** It moves the grid median by
0.004 and absorbs **3.7%** of the excess, against the ≥50% M3 required.

**The paired design validates itself.** The control arm reproduces the published
grid median to **0.002** (1.103 vs 1.105) despite a different RNG realization and
a spliced sub-8 segment — so realization noise is about half the size of the
effect being measured, which is what makes a 0.004 shift interpretable at all.

| | prediction | verdict |
|---|---|---|
| **M1** | ratio falls in every subgroup | **PASS** — all three fell; 48 of 54 cells |
| **M2** | reduction follows turnoff ordering A < B < C | **FAIL** — A 0.0029, B 0.0035, **C 0.0011** |
| **M3** | absorbs ≥ half the excess (median < 1.222) | **PASS literally, FAIL on the governing reading** |

**M2 failed, and not merely for lack of power.** CygOB2-C's reduction is smaller
than B's by ~2.7σ — the *opposite* of the predicted turnoff ordering. The
pre-registered consequence applies as written: *the excess and the multiplicity
effect have different mass dependences, so multiplicity cannot be the whole story
even if it lowers the ratio.*

**Why M3's literal pass does not count.** Its 1.222 threshold was derived
arithmetically from a grid median of 1.444 that **issue #17 has since withdrawn**.
The corrected control arm already sits at 1.103 — below 1.222 *before the
treatment arm changes anything*. Applying the decision rule to that would
conclude "the shallow-IMF signal is a multiplicity artefact" while the
measurement shows multiplicity absorbed 3.7% of it. **The pre-registration was
not amended** — both readings are recorded, and the relative form, which is what
M3's own sentence states, governs.

Consistency check: the effect grows with α (0.0011 / 0.0022 / 0.0042 at
α = 2.0 / 2.3 / 2.6), the same low-mass-weighting mechanism that drove F2 in §2.

### Decision rule applied

The **"if M3 fails"** branch: the IMF reading is reported with the measured
multiplicity correction applied and the remainder carried as a systematic.
**No repair_v7 is triggered by issue #15.** The disfavouring of α = 2.6 stands —
as it does under either branch.

### Scope limit — stated because it is load-bearing

This tests multiplicity **above 8 M☉ only**. A star already above 8 M☉ cannot be
scattered *into* the census by being made brighter, so this is precisely the
range where the mechanism has **least** room to act. The 4–8 M☉ up-scatter
channel — where unresolved companions would matter most, and which issue #17
showed carries real weight — is held at f_bin = 0.40 in **both** arms.

The correct reading is **"multiplicity above 8 M☉ does not explain the excess"**,
not "multiplicity does not explain the excess".

### The sub-8 M☉ channel was then measured — and it is not negligible

*Pre-registered: [wp5_fbin_discriminator_prereg.json](provenance/wp5_fbin_discriminator_prereg.json)
· scored: [wp5_fbin_discriminator_execution.json](provenance/wp5_fbin_discriminator_execution.json).*

A three-node paired discriminator extended f_bin below 8 M☉ (0.40 at 2 → 0.55 at
8 → 0.70 at 16):

| | measured | threshold |
|---|---:|---:|
| **D1** — recovery over the 2–8 M☉ calibration window | **+0.16% ± 0.09%** | 2% |
| **D2** — R(estimated > 8 M☉ \| M) over 4–8 M☉ | **+9.87% ± 1.58%** | 2% |

**`repair_v7` is justified.** The control arm reproduced the accepted `repair_v6`
node byte for byte, so the shift is the model change and nothing else.

The 4–8 M☉ segment carries **23.5%** of the predicted observable count above
8 M☉, so a +9.87% shift there raises the total predicted count by **2.3%** and
would move the grid median **1.099 → ~1.074** — in the direction of better
closure.

**Why the contrast between D1 and D2 matters more than either number.** Inside
the calibration window most stars are recovered either way, so the recovery
fraction barely moves. At the 8 M☉ boundary a small brightness shift converts
directly into a crossing probability. **A truth-model error shows up ~25× more
strongly at a threshold than in a window-averaged recovery fraction** — the same
structural lesson as issues #3 and #17, and the reason this one was nearly
missed.

## 7. Runaways

Standard 2D proper-motion traceback with full covariance sampling (256 draws per
star). **3,337 OB candidates outside the footprint; 119 recovered raw.**

The number that matters is the chance-alignment rate, and it is **measured, not
modelled** — control fields traced back by identical code. A first attempt gave
a *negative* corrected count because control stars sat at systematically smaller
separations than real candidates, and chance recovery falls steeply with
separation (20.7% at 1–1.5°, 0.4% beyond 8°). Corrected by measuring the
false-positive rate **per separation bin** and applying it at each candidate's
own separation: **1.9% effective chance rate → 54.9 genuine runaways**
(binned-and-clipped in the ledger: 58.2; both reported, and WP7 should carry the
difference as a systematic).

### Issue #16 — the traceback used the wrong proper motions

**A previous version of this section, reporting 260 raw and 109 corrected, is
withdrawn.** It traced back **absolute** Gaia proper motions. Cyg OB2's systemic
motion is **(−2.707, −4.317) mas/yr — larger than a typical ejection
signature**, so the traceback was measuring the association's bulk drift rather
than the ejection. The fix subtracts the systemic motion, rotated into galactic
coordinates at each star's own position, with control fields treated identically
so the false-positive rate stays comparable.

**No internal check could have found this.** It was caught by the external gate
in §10, on a star whose answer is known independently.

**Issue #10 quantified.** At the WP2 distance of 1.62 kpc the ±8° box caps
recovery at **110.6 km/s over 2 Myr and 44.2 km/s over 5 Myr**. (CUTS §4.2
quoted 96/38 using a round 1.4 kpc; these supersede it.) The plan's stated
10–100 km/s range is **not achievable** across the full age baseline. Adopted
wording: *the runaway search is complete to v ≲ 44 km/s for ejections up to
5 Myr ago; faster or older ejections lie outside the footprint and make N_SN a
lower bound.*

The traceback is **2D** — no radial velocities are used — so line-of-sight
ejections are missed and the count is a lower bound for that reason too.

## 8. The ledger handed to WP7

[tables/wp6_massive_census.cat](tables/wp6_massive_census.cat)

| channel | entries | summed weight |
|---|---:|---:|
| member | 632 | 295.40 |
| orphan_anchor | 27 | 27.00 |
| runaway | 119 | 58.18 |
| **total living above 8 M☉** | | **380.58** |

By subgroup: A 66.41, B 87.55, C 119.44, unassigned (runaways) 58.18.

Runaway subgroup is **unassigned** by design: a 2D traceback constrains the
origin footprint, not which subgroup inside it, and assigning them per subgroup
would invent precision.

The binned figure (58.18) and the aggregate (54.9) now agree to 6%. Under the
withdrawn absolute-PM traceback they disagreed by 37% (149 vs 109) — the
corrected traceback is markedly more internally consistent, which is weak
independent support for the fix.

## 9. What WP7 must do with this

1. **Carry α = 2.6 as disfavoured** by the massive-star census, and quote N_SN
   per branch with that noted. Do **not** refit α to make the census close —
   that would convert the only out-of-sample check the analysis has into a
   fitted parameter.
2. **Do not quote the withdrawn figures.** Any downstream text still carrying a
   45% excess, a grid median of 1.444, or a closing α of 2.070 predates issue
   #17 and must be corrected.
3. **Report the residual IMF reading with the multiplicity correction applied**
   (§6): multiplicity above 8 M☉ accounts for 3.7% of the excess, measured, and
   the remainder is a systematic. Multiplicity *below* 8 M☉ is untested and is
   the one instrumental channel that could still move this.
4. **Subtract runaways** from the missing-equals-dead bookkeeping; they are
   living stars that left, so they reduce N_SN.
5. **Treat the runaway count as a lower bound** on both the velocity ceiling and
   the 2D-projection grounds.
6. **Carry the orphan-anchor systematic**: their masses come from spectral type,
   not the WP4 posterior, and 13 have no spectral type at all.
7. **CygOB2-C's residual is the open science question.** After both fixes, C
   still sits at **1.405** while A is at 0.894 — they disagree in *direction*, so
   no single mechanism explains both. Three untested candidates, in priority
   order:
   - **subgroup-label uncertainty** (issue #6);
   - **distance contamination** — a star truly at 1350 pc but assumed at 1620
     gets a distance modulus 0.40 mag too large, is inferred more massive, and
     inflates the count above 8 M☉. This is C's direction, and it is covered by
     **no** WP6 alternative: A4 tested membership weighting, not distance. See
     [berlanas_2019_two_distance.md](cross_checks/berlanas_2019_two_distance.md);
   - the **+2.4% Orellana distance systematic**, which pushes C further up.
8. **Multiplicity below 8 M☉ is untested** and is the one remaining instrumental
   channel that could still move the result. It requires a `repair_v7` chain
   re-run because it perturbs the accepted WP5 calibration.

## 10. Gate — the external cross-match

*[wp6_runaway_crossmatch_execution.json](provenance/wp6_runaway_crossmatch_execution.json)
· [tables/wp6_runaway_crossmatch.csv](tables/wp6_runaway_crossmatch.csv)*

| criterion | status |
|---|---|
| closure ratio explained within quantified contributions | **met** — the dominant contributor was issue #17 (~3/4 of the excess); six alternatives tested and bounded, multiplicity measured at 3.7% |
| runaway search reproduces literature candidates | **met** — 1/1 testable candidate recovered |

**WP6 is closed.** Both gate criteria are met and no method question remains
open. What remains is a *science* question: CygOB2-C's residual ratio of 1.405,
discussed in §9.

Nine published candidates were resolved from SIMBAD and frozen into the script
with the identifier queried, so the check reruns offline.

**Eight of the nine are not testable.** A15, A20, A26, A37, A46, MT91-516,
Cyg OB2 #10 and #4B all sit **0.15–0.70° from the centroid — inside the
association**. They are stars with anomalous proper motions that have not left.
The escaped-star search excludes them by construction, not by failure.

**The gate rests on BD+43 3654** (Comerón & Pasquali 2007), the one candidate
genuinely outside the footprint at 2.78°:

| | measured here | literature |
|---|---:|---:|
| recovery probability | **1.000** | — |
| implied velocity | **38.8 km/s** | ~40 |
| flight time | **1.36 Myr** | 1.6 Myr (age) |
| position match | 0.04″ | — |

Both quantities agree, from an estimator that knows nothing about the published
values. **This is the check that caught issue #16**: under the absolute-PM
traceback the same star scored exactly 0.000.

**What this gate does and does not establish.** It rests on a single star, so it
is a strong test of the method's sign and scale and a weak test of its
completeness. That limit is carried in the provenance record.

**Incidental finding.** Only 2 of the 8 in-footprint literature candidates are
WP2 members. Berlanas+2020 flag most of the rest as RUWE > 1.4, which is exactly
what the WP2 quality filter removes — the same loss channel the orphan-anchor
count describes, independently confirmed. They are **not** added to the census:
the response already models that loss, so folding them in would double-correct.
Note also that unresolved multiplicity inflates RUWE, so this population
overlaps the one §6 is testing; the two are not independent evidence.

## 11. External cross-checks

Full documents in [cross_checks/](cross_checks/); record:
[wp6_external_crosschecks_execution.json](provenance/wp6_external_crosschecks_execution.json).

| source | quantity | ours | theirs |
|---|---|---|---|
| Härer+2025 | SN rate | 0.78 ev/100 kyr (α=2.3) | "several per Myr", Fig. 2 range 0.25–2 |
| Wright+2015 | association mass | 1.74 × 10⁴ M☉ | 1.65 × 10⁴ M☉ |
| Orellana+2021 | systemic μ_α* | −2.7067 | −2.71 ± 0.02 |
| Orellana+2021 | distance | 1616 pc | 1669 ± 6 pc (**+2.4%** systematic) |

The proper-motion agreement independently validates the issue #16 fix — the
vector §7's traceback subtracts is confirmed by an external DR2 analysis to
0.003 mas/yr.

**Caution on the SN rate.** Härer's Fig. 2 axis is *events per 100 **kyr***, not
per 100 years. Misreading it inflates the expected rate by 1000×.

**Our N_SN is a lower bound relative to BPASS**: their rate is dominated by type
Ic at 3–5 Myr, the stripped-envelope channel, which requires a binary companion.
Our count is single-star turnoff counting. This is issue #15's physics reached
from another direction.

## 12. Reproduction

```bash
export PYTHONPATH=scripts WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging
python3 scripts/wp6_mass_extension_decision.py
python3 scripts/wp6_massive_injections.py
python3 scripts/wp6_verify_wp5_untouched.py
python3 scripts/wp6_massive_census.py
python3 scripts/wp6_closure_floor_prereg.py   # BEFORE the corrected run
python3 scripts/wp6_closure_test.py           # issue #17 floor, 4.0 Msun
python3 scripts/wp6_closure_floor_score.py    # scores F1-F4
python3 scripts/wp6_closure_attribution.py
python3 scripts/wp6_runaways.py
python3 scripts/wp6_runaway_crossmatch.py
python3 scripts/wp6_ledger.py

# issue #15 — the multiplicity diagnostic (about 3 h)
python3 scripts/wp6_multiplicity_prereg.py       # BEFORE the injections
python3 scripts/wp6_multiplicity_injections.py   # V4 preflight, then both arms
python3 scripts/wp6_multiplicity_closure.py      # scores M1/M2/M3

# external validation
python3 scripts/wp6_external_crosschecks.py
```
