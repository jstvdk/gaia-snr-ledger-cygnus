# Cuts and Thresholds — How to Choose Every Number in This Pipeline

*Companion to `GLOSSARY.md` (what the terms mean) and `paper1_execution_plan.md` (what the steps are). This document is about the numbers: where each one comes from, whether it is derived, conventional, measured or free, and what you owe a referee for each.*

---

## 1. The central idea: every number belongs to a class

The most common way analyses like this go wrong is not picking a bad number — it is **not knowing what kind of number it is**. A value that should have been measured from the data gets guessed; a value that should have been scanned gets fixed; a genuine model ambiguity gets silently resolved by picking a favourite.

Five classes. The class determines what you do:

| Class | Definition | What you owe | Failure if misclassified |
|---|---|---|---|
| **A — Derived** | Follows from a calculation. Change the physics, the number changes. | Show the derivation in the paper. | Number drifts out of sync with the science it encodes. |
| **B — Conventional** | Standard in the literature; defensible by citation. | Cite it, and demonstrate the result survives moving it. | Referees wave it through — so *errors hide here*. |
| **C — Measured** | The data determines it. Guessing it is simply wrong. | Show the measurement (curve, knee plot, injection test). | You invent a value the data contradicts. |
| **D — Free** | No principled value exists. | Scan it; require the result to be stable across the scan. | A result that exists only at one setting is not a result. |
| **E — Branch** | Genuine unresolved model ambiguity. | Carry all values in parallel to the end; never average. | Model systematic gets buried inside a quoted error bar. |

Your plan already handles Class E well — the mandatory branches in §1.4 are exactly right. The gap is that **A, B, C and D are all currently being treated the same way**: chosen once, written into a notebook, and not revisited. The WP2 failure is precisely this — `eps` is a Class C/D quantity that was treated as Class B.

---

## 2. Master table — every number currently in the pipeline

| # | Cut | Current value | Class | Status | Action |
|---|---|---|---|---|---|
| 1 | Narrow sky box | *l* ∈ [77,83], *b* ∈ [−1.5,4] | A | derived, generous | keep; document derivation |
| 2 | Wide sky box | *l* ∈ [72,88], *b* ∈ [−5,8] | A | **under-derived** | see §4.2 — sets a velocity ceiling |
| 3 | Parallax window | 0.35–1.10 mas | A | derived but **truncating** | see §4.1 — must be modelled downstream |
| 4 | Magnitude limit | G < 19 | A/C | derived; completeness unmeasured | see §4.3 |
| 5 | RUWE | < 1.4 | B | convention, **binary-biased** | see §5.1 — sensitivity test mandatory |
| 6 | Visibility periods | ≥ 8 | B | convention, low risk | test ≥ 8 vs ≥ 9, expect null |
| 7 | BP/RP excess | non-null only | C | **inadequate** | see §5.2 — implement C* cut |
| 8 | Zero-point recipe | Lindegren+21 | B | correct | audit boundary flags |
| 9 | DBSCAN `eps` | 0.42 | **C+D** | **wrong — percolating** | see §6.1 — k-distance knee, then scan |
| 10 | DBSCAN `min_samples` | 15 | B/D | plausible, untested | see §6.2 |
| 11 | Min cluster size | 20 | D | arbitrary | scan with `min_samples` |
| 12 | MC draws | N = 100 | A | **too few** | see §6.3 — σ<sub>P</sub> = 0.05 |
| 13 | Acceptance locus | 3 robust scales | D | arbitrary **and circular** | replace with field model, §6.4 |
| 14 | Soft membership | P > 0.05 | B | fine as a soft floor | keep; never hard-cut for IMF |
| 15 | Gate membership | P > 0.5 | B | fine | keep |
| 16 | BIC preference | ΔBIC < 0 | B | **invalid at this N** | see §6.5 |
| 17 | R<sub>V</sub> | 3.1 base; 3.0/3.5 | B+E | correct | keep |
| 18 | Isochrone families | PARSEC, MIST | E | correct | keep |
| 19 | Isochrone age grid | 1–10 Myr log | A | fine | keep |
| 20 | Binary fraction | 0.3–0.5 | E | ~~correct~~ **range too narrow above ~10 M<sub>☉</sub>** | see §18 — measured O-star value is 0.70, outside the carried range |
| 21 | Calibration window | 2–8 M<sub>☉</sub> | **C** | **must be measured** | see §7.1 — likely needs raising |
| 22 | Completeness floor | ≥ 95% | B | fine | keep |
| 23 | Min stars/subgroup | ≥ 50 | B | fine | keep |
| 24 | IMF slope | 2.0 / 2.3 / 2.6 | E | correct | keep |
| 25 | SN mass threshold | 8 M<sub>☉</sub> | B+E | convention | pair with explodability branch |
| 26 | Explodability | all-explode / islands | E | correct | keep |
| 27 | SF duration | 0 / 1 / 2 Myr | E | correct | keep |
| 28 | Runaway velocity | 10–100 km/s | A | **inconsistent with box** | see §4.2 |
| 29 | Recent-SN window | 100 kyr | A | external (Härer+25) | quote verbatim, never tune |
| 30 | WP1 gate | ≥ 90% Wright+15 | B | fine | keep |
| 31 | WP2 gate | ≥ 80% Berlanas | B | **recall-only** | see §8 |
| 32 | Mass sanity | factor ~2 | B | fine | keep |
| 33 | WP4 age-indicator minimum | ≥ 15 stars | B | registered at closure | exclude smaller samples; report 10/20/30-star sensitivity |

---

## 3. Four cross-cutting principles

### 3.1 Irreversible cuts must be generous; reversible cuts must be flags

A cut applied **in the ADQL query** (parallax window, G < 19, sky box) is irreversible — undoing it means re-querying the archive and re-downloading. A cut applied **in analysis** is free to undo.

> **Rule: make query-time cuts wider than you think you need, and make analysis-time cuts flags rather than deletions.**

Your plan states this correctly ("These are analysis flags, not deletions from WP1") and your WP1 window is appropriately generous. Keep this discipline — it is why you can fix WP2 without touching the archive.

### 3.2 Any cut correlated with magnitude is a cut correlated with mass

This is the one that silently destroys IMF work. Brightness maps to mass. So **every** cut that preferentially removes faint stars also preferentially removes low-mass stars, and therefore tilts the mass function you are trying to measure.

Affected cuts: G < 19 (directly), parallax S/N (faint stars have worse parallaxes), RUWE (see §5.1), BP/RP excess (crowding hits faint stars harder), and the parallax window itself (large errors scatter faint stars out).

> **Rule: for WP5 you do not need the cut, you need the *selection function* — the probability of surviving the full chain as a function of mass and position.** That is what injection testing measures. Never substitute a cut value for a measured recovery curve.

### 3.3 Truncation must be carried into every distribution fit

Your parallax window 0.35–1.10 mas is a **hard truncation**. The parent distribution extends beyond it; you observe only the interior.

This matters enormously for the Berlanas two-population test. Fitting a Gaussian mixture to truncated data without modelling the truncation misfits the tails, and the fit will often compensate by adding a spurious component. Combined with the BIC-at-large-N problem (§6.5), that is very likely part of what produced your ΔBIC = −25,265.

> **Rule: any likelihood fitted to parallaxes must include the window as part of the model, or be restricted to a region where truncation is negligible.**

### 3.4 Distinguish "the cut changes the result" from "the cut changes the sample"

A sensitivity test that changes the sample size will change error bars — that is expected and uninteresting. What you are testing is whether the **central value moves by more than the statistical uncertainty**. Report both, and say which you mean.

---

## 4. Class A — the derived numbers, with derivations

### 4.1 Parallax window: 0.35–1.10 mas

Cyg OB2 spans 1.35–1.75 kpc → ϖ = **0.571–0.741 mas**. The window 0.35–1.10 mas corresponds to **909–2857 pc**, giving generous margin on both sides. Correct in spirit.

Two things to record:

- **Margin justification.** The window must be wide enough that stars scattered by measurement error still land inside. At G = 19, σ<sub>ϖ</sub> ≈ 0.4 mas — comparable to the entire window half-width. So for faint stars the window is *not* generous; it is comparable to their error. State this.
- **Truncation.** See §3.3. The window enters any parallax likelihood.

### 4.2 Wide box: the number that is quietly inconsistent

Your plan asks WP6 to recover runaways at **10–100 km/s**. The wide box is *l* ∈ [72,88], *b* ∈ [−5,8] — a half-width of ~8° about the association.

At 1.4 kpc:

| Half-width | Physical | Max velocity recoverable over 2 Myr | over 5 Myr |
|---|---|---|---|
| 4° | 98 pc | 48 km/s | 19 km/s |
| 6° | 147 pc | 72 km/s | 29 km/s |
| **8°** | **196 pc** | **96 km/s** | **38 km/s** |
| 10° | 244 pc | 120 km/s | 48 km/s |

So the current box recovers a 100 km/s runaway only if it was ejected within the last ~2 Myr, and over a 5 Myr baseline it caps out near **38 km/s**. The stated 10–100 km/s range is not achievable with this footprint.

**This is not necessarily wrong — but it must be stated as a bound, not left implicit.** Either widen the box (costs query volume and contamination) or state explicitly: *"our runaway search is complete to v ≲ 38 km/s for ejections up to 5 Myr ago; faster or older ejections are outside the search footprint and make our N<sub>SN</sub> a lower bound."* The second is honest and cheap. The plan already flags this qualitatively in WP6's caveats — this quantifies it.

### 4.3 Magnitude limit G < 19, and what it costs you

Distance modulus at 1.4 kpc: **μ = 10.73 mag**.

For a 2 M<sub>☉</sub> star (M<sub>G</sub> ≈ +1.4, the bottom of your calibration window), with A<sub>G</sub> ≈ 0.86 A<sub>V</sub>:

| A<sub>V</sub> | A<sub>G</sub> | Apparent G | Inside G < 19? |
|---|---|---|---|
| 4 | 3.4 | 15.6 | yes |
| 6 | 5.2 | 17.3 | yes |
| **8** | **6.9** | **19.0** | **borderline** |
| 10 | 8.6 | 20.7 | **no** |
| 15 | 12.9 | 25.0 | no |
| 20 | 17.2 | 29.3 | no |

Cygnus extinction runs A<sub>V</sub> = 4–20. So **the 2 M<sub>☉</sub> lower edge of your calibration window is only reachable in the least obscured sightlines**, and is lost entirely above A<sub>V</sub> ≈ 8.

Three consequences:

1. The calibration window lower edge almost certainly needs raising — your plan anticipates this ("e.g. 3–8 at Cygnus distance/extinction; document the choice"). Now you know roughly why and by how much.
2. **Completeness is a function of position, not just mass**, because extinction is patchy. A single completeness curve for the whole association is not enough — it must be per subgroup at minimum, and ideally account for the extinction distribution within each.
3. Raising the window's lower edge shrinks the number of calibration stars, which weakens the IMF normalisation. There is a real trade-off here between completeness and counting statistics, and it should be shown, not hidden.

*(These magnitudes are order-of-magnitude guides using an approximate A<sub>G</sub>/A<sub>V</sub> and a rough M<sub>G</sub>. Replace with actual isochrone + extinction-law lookups in WP5 — but the conclusion that A<sub>V</sub> ≳ 8 kills 2 M<sub>☉</sub> at G < 19 is robust.)*

---

## 5. Class B — conventions, and where they bite

### 5.1 RUWE < 1.4 — the most dangerous "safe" number

1.4 is an empirical threshold from Gaia's own technical notes, not a physical boundary. It is universally accepted, which is exactly the problem: a referee will not question it, so if it biases you, nobody catches it.

**The bias:** RUWE is elevated by unresolved binaries. O-star multiplicity is ~70%. So `RUWE < 1.4` preferentially deletes the massive binaries that are central to your census — and the deletion rate is *mass-dependent*, which is the §3.2 failure mode.

**What to do:**

- Exempt the spectroscopic anchors entirely (your plan already says this — good).
- Run the census at RUWE < 1.4 and RUWE < 2.0 (and no cut, with excess-noise flagging instead) and compare the high-mass counts. If the closure ratio in WP6 moves, that is a real systematic and belongs in the error budget.
- Report the fraction removed by RUWE **as a function of magnitude**. If it rises steeply at bright magnitudes, you are cutting massive binaries.

### 5.2 BP/RP excess factor — currently not actually a cut

Your filter checks only that the value is *present*. That is not a quality cut.

The correct treatment (Riello et al. 2021) is the **colour-corrected** excess factor C\*, which removes the intrinsic colour dependence, followed by a magnitude-dependent sigma clip: keep |C\*| < N·σ<sub>C\*</sub>(G), with N typically 3–5.

This matters more here than in most fields: Cygnus is crowded and nebulous, which is exactly what inflates BP/RP excess. Without this cut, contaminated colours propagate into WP3 extinctions and WP4 masses. Treat N as Class D and scan it.

### 5.3 The ones that are fine

`visibility_periods_used ≥ 8`, R<sub>V</sub> = 3.1 baseline, the 8 M<sub>☉</sub> SN threshold, ≥95% completeness floor, ≥50 stars per subgroup, the factor-2 mass sanity check. Cite, test cheaply where free, move on. Do not spend effort here.

### 5.4 WP4 age-indicator measurability: N ≥ 15

An age branch is reported as measurable only when its CMD indicator window
contains at least **15 stars**. This is a Class-B reporting convention: below
that count, an apparent MAP can be set by one or a few stars rather than a
population sequence. It is a measurability flag, not a membership cut; no star
is removed from any catalogue.

The closure sensitivity check repeats the envelope selection at floors of
10/20/30 stars. The 15-star baseline retains 104/132 branch rows; the excluded
rows are enumerated in `tables/wp4_ages_summary.csv`. Boundary-pinned MAPs are
excluded independently, even if they meet the count floor.

---

## 6. Class C and D — the numbers that must be measured or scanned

### 6.1 DBSCAN `eps` — measure it, then scan it

**Do not guess this.** The standard method:

1. For every point, compute the distance to its *k*-th nearest neighbour, where *k* = `min_samples`.
2. Sort those distances and plot them.
3. The curve has a **knee**: below it, points are in dense regions; above it, they are isolated. The knee is the natural ε.

Then scan ε across a range spanning the knee and watch the right diagnostic. **Cluster count is the wrong diagnostic.** Track instead:

- **largest-cluster fraction** — must be small; if it approaches the field size, you are percolating
- **member spatial extent** — must be smaller than the selection box
- **control-field yield at the same ε** — your false-positive rate

Your scan of 0.34/0.42/0.50 failed to reveal the problem because all three values sit in the percolated regime. A useful scan must bracket the transition — start near 0.03 and go up geometrically.

**Consider HDBSCAN instead.** It builds a hierarchy across all density scales rather than fixing one ε, which suits a region whose density varies by orders of magnitude. Paíz+25 used it here successfully.

**Also reconsider the feature space.** Clustering jointly on (*l*, *b*, ϖ, μ<sub>α*</sub>, μ<sub>δ</sub>) after RobustScaler lets sky position trade against kinematics — two stars far apart on the sky can be "close" if their proper motions match. Many groups cluster in proper motion + parallax only, then use spatial coherence as an *independent* confirmation. That is methodologically cleaner: it keeps one axis free for validation instead of consuming all of them in the fit.

### 6.2 `min_samples` — heuristic floor, then scan

Common guidance: `min_samples` ≥ D + 1, often 2D, for D dimensions. With D = 5 that suggests ~10–15, so your 15 is defensible as a starting point. But it interacts strongly with ε, so scan the two **jointly on a grid**, not one at a time. Report the stability region.

### 6.3 Monte Carlo draws — this one is pure arithmetic

Membership probability from N draws is a binomial proportion; its standard error at p = 0.5 is √(0.25/N):

| N | σ<sub>P</sub> |
|---|---|
| 100 | 0.050 |
| 500 | 0.022 |
| 1000 | 0.016 |
| 2500 | 0.010 |
| 10000 | 0.005 |

At N = 100 your probabilities carry ±0.05 noise — larger than the P > 0.05 threshold itself. Fine for development, not for publication.

> **Rule: choose N so σ<sub>P</sub> is small compared to the smallest probability difference that matters.** If P > 0.05 is a meaningful boundary, you need σ<sub>P</sub> ≲ 0.01, so **N ≥ 2500**. N = 10,000 costs little and removes the question.

### 6.4 The acceptance locus — replace, don't tune

"Within 3 robust scales of the candidate centre" has two problems, and only one is about the number 3.

The number is arbitrary (Class D). But the deeper issue is **circularity**: the centre and scale are computed from the candidate set, then that same set is tested against them. Nearly everything passes by construction — hence 85% of your members at P > 0.95.

**No value of "3" fixes this.** You need a model with two components:

- a **cluster** model (compact in parallax and proper motion), and
- a **field** model (broad, from control fields or from the data outside the cluster locus),

with membership as the posterior odds P(cluster)/[P(cluster) + P(field)]. Then there is no arbitrary radius at all — the probability falls out of the likelihood ratio. This converts a Class D free parameter into a Class C measured quantity, which is always the right direction to move.

### 6.5 BIC — why your ΔBIC is not evidence

Conventional reading (Kass & Raftery): ΔBIC of 2–6 is positive evidence, 6–10 strong, >10 very strong. **These thresholds were never meant for N = 160,000.**

BIC's complexity penalty grows as *k* ln N, but the log-likelihood gain from adding a component grows as *N*. At large N the likelihood term always wins. So BIC will prefer more components for **any** distribution that is not exactly Gaussian — merely being skewed is sufficient. A field parallax distribution inside a truncated window along a spiral arm is emphatically not Gaussian.

Your ΔBIC = −25,265 is therefore consistent with "the sample is large and not Gaussian," which is guaranteed *a priori*. It is not evidence for two physical populations.

**What to do instead:**

1. Run the test on a **clean member sample**, not 160k field stars.
2. Use **extreme deconvolution** so per-star parallax errors are removed rather than being free to manufacture components.
3. Include the **truncation window** in the likelihood (§3.3).
4. Work in the space where the claim lives — Berlanas' claim is about **distance**, and a Gaussian in parallax is not a Gaussian in distance.
5. Validate on **control fields**: run the identical test where you *know* there is no association. If control fields also "prefer" two components, the test is measuring the field, not the object.
6. Cross-check any preferred split against **independent** information — sky position, proper motion, extinction, spectroscopy. A split that appears only in parallax is a statistical artefact until something else confirms it.

Step 5 is the one that would have caught this immediately, and it costs almost nothing.

---

## 7. Class C in WP5 — the calibration window

### 7.1 The lower edge is measured, never chosen

The nominal 2–8 M<sub>☉</sub> window has a principled *upper* edge (Class A: 8 M<sub>☉</sub> is where stars start dying, so above it the surviving population no longer traces the birth population) but its **lower edge is Class C — it is whatever the completeness curve says it is.**

Procedure:

1. Injection test: insert synthetic stars of known mass, run the *entire* selection chain (query cuts, quality cuts, clustering, membership), measure recovery fraction vs mass.
2. Find the mass where recovery ≥ 95%.
3. That is your lower edge. If it is 3.5 M<sub>☉</sub>, the window is 3.5–8, and you say so.

From §4.3, expect this to land above 2 M<sub>☉</sub>, and to differ between subgroups because extinction differs.

**The injection must go through the full chain, not just the magnitude limit.** A star can be bright enough to detect yet still lost by the RUWE cut, the BP/RP cut, or by DBSCAN failing to assign it. Only end-to-end injection captures that.

**Amendment, 2026-07-27 — the 95% edge was not reachable on this field.** No subgroup reaches an absolute 95% recovery edge anywhere; the bright plateau tops out at 0.83–0.85 because ~15% of *bright* injected stars are lost by the WP2 quality filter, which is magnitude-independent and therefore never clears with increasing mass. Steps 2–3 above are consequently unexecutable as written. The predeclared fallback is in force: the fit uses the full nominal 2–8 M<sub>☉</sub> window with the injection curve inside the Poisson intensity, labelled `corrected_no_absolute95_edge`, and is **explicitly not** described as 95% complete. Relative-to-plateau edges are diagnostic only. This supersedes steps 2–3 for Cyg OB2; it does not license skipping the measurement on another field.

### 7.2 The window edge and the parent integration range are two different numbers

*Added 2026-07-27, after this conflation blocked the WP5 gate. See `reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md`.*

§7.1 gives the upper edge a Class A justification — 8 M<sub>☉</sub> is where stars start dying, so above it the surviving population no longer traces the birth population. That is correct, and it is a statement about **which stars you count**. It says nothing about **which stars the forward model is allowed to contain**, and the two must not be tied together.

The mass estimator has a finite width (~23% in log mass, measured from the WP4 injection response). Any star whose *true* mass lies within roughly 3σ of a window edge can be *measured* onto the other side of it. So the parent population the response integrates has to be wider than the observed window **on both sides**:

| | value | class | justification |
|---|---|---|---|
| Observed counting window | 2–8 M<sub>☉</sub> | A (upper), C (lower) | §7.1 |
| Injected parent range | 0.5–18 M<sub>☉</sub> | A | must cover every true mass that can be measured into the window |

The downward extension to 0.5 M<sub>☉</sub> was present from the start. The upward extension was not: the parent was truncated at the same 8.0 M<sub>☉</sub> value as the window, so the ~300 living members above 8 M<sub>☉</sub> had no term in the model, and the top observed bin was under-predicted in every subgroup. The ceiling of 18 M<sub>☉</sub> is measured, not chosen — it is where the >8 M<sub>☉</sub> contribution to the top-bin rate converges to 100.0% for every IMF slope branch (99.9% at 16, 95.8% at 12).

> **Rule: whenever a fit bins a quantity that carries measurement error, the model's parent range must exceed the binned range by several times the error width at every edge. Truncating the parent at the bin boundary is the same error as ignoring truncation (§3.3), applied in the opposite direction.**

---

## 8. Gates are numbers too — and yours has a hole

Every gate threshold is itself a Class B choice and deserves the same scrutiny.

**The WP2 gate is currently one-sided.** "Recover ≥ 80% of Berlanas+19 members at P > 0.5" measures **recall** only. A pipeline that labels every star a member scores 100% recall and sails through. That is not hypothetical — it is what your current run does, and the gate would have passed it.

A gate must be a **conjunction of recall and precision**:

| Criterion | Type | Suggested threshold |
|---|---|---|
| Berlanas+19 recovery at P > 0.5 | recall | ≥ 80% (as written) |
| Control-field yield / association yield | precision | ≤ 5–10% |
| Total member count | sanity | 10²–10⁴, not 10⁵ |
| Member spatial extent vs selection box | sanity | strictly smaller |
| Largest cluster / analysis sample | sanity | ≪ 1 |

The last three are nearly free to compute and would each have caught the current failure independently. **Cheap sanity checks that can only fail loudly are worth more than expensive ones that fail quietly.**

---

## 9. How to run a sensitivity test that means something

1. **Fix everything else.** One knob at a time, or a proper grid — never drift multiple values together.
2. **Choose the range around the default**, spanning both the literature range and the point where you expect breakdown. A scan that never breaks tells you nothing about where the edge is.
3. **Record the downstream quantity, not the intermediate.** Nobody cares how the member count varies with ε; they care how **N<sub>SN</sub>** varies with ε. Push each variation all the way to the end.
4. **Compare the shift to the statistical error.** A change smaller than the error bar is not a systematic worth reporting; a change larger than it is, and belongs in the budget.
5. **Report the driver.** WP9 asks for the dominant sensitivity — that answer comes from this exercise, so structure it to produce a ranked list.
6. **Log every run**, including the ones you rejected. The provenance appendix should show what you tried, not just what you kept.

---

## 10. What to record for every number

For the reproducibility appendix, each cut gets one row:

| Field | Example |
|---|---|
| Cut name | DBSCAN `eps` |
| Value used | 0.075 |
| Class | C + D |
| Justification | k-distance knee at 0.07; stable over 0.06–0.09 |
| Range tested | 0.03–0.50 |
| Effect on N<sub>SN</sub> | +0.3 / −0.2 across stability region |
| Where set | `notebooks/wp2_membership_and_substructure.ipynb` §3 |
| Date / version | 2026-07-22, DR3 |

If you cannot fill in the Class and Justification fields for a number, that number is not yet ready to be in the pipeline.

---

## 11. Priority order

Given where the project stands:

1. **§6.1 `eps`** — nothing downstream is meaningful until membership works.
2. **§6.4 field model** — the membership probabilities are currently tautological.
3. **§8 gate revision** — do this *before* rerunning, or the rerun cannot fail informatively.
4. **§6.5 BIC redo** — the current two-population result would be an embarrassing error to publish, precisely because it agrees with the literature.
5. **§5.2 BP/RP excess** — cheap, and it protects WP3.
6. **§6.3 MC draws** — trivial fix, N = 10,000.
7. **§7.1 calibration window** — needs the injection machinery; comes with WP5.
8. **§4.2 wide-box statement** — one honest sentence in the manuscript; no code required.

---

## 12. Versioned WP3 extinction / WP4 mass repair (`repair_v1`, 2026-07-24)

These settings apply only to the versioned repair products. The frozen WP3,
WP4, and blocking WP5 artifacts remain unchanged.

| Quantity | Adopted value | Class | Justification and effect |
|---|---:|---|---|
| Per-band calibration/model floor in the WP3 likelihood | 0.03 mag in G/BP/RP/J/H/Ks | D | Added in quadrature to prevent millimag Gaia errors from suppressing 2MASS. The 0.02–0.03-mag range was specified before the repair; the conservative endpoint is used. |
| Spatial extinction neighbours | 8 spectroscopic anchors | B/C | Fixed by the repair brief and evaluated on the same angular scale as the smoking-gun diagnostic. Median eighth-neighbour separation is 0.071°. |
| Spatial-prior width | 0.452, 0.453, 0.475 mag at R_V=3.0, 3.1, 3.5 | C | Measured as the robust leave-one-out scatter of each anchor around the median of its eight nearest other anchors. This preserves real differential extinction rather than replacing it with a smooth map. |
| Template-branch uncertainty | 1.71, 1.67, 1.70 mag at R_V=3.0, 3.1, 3.5 | C | Calibrated on all 149 spectroscopic anchors as the larger asymmetric central-68% endpoint of broadband-minus-spectroscopic A_V. This exposes the hidden PMS/ZAMS branch ambiguity instead of reporting a 0.05-mag local-curvature error. |
| Template quadrature resolution | at most 64 log-mass representatives per family/age cell | D | Computational quadrature for applying the identical estimator to 223,200 injections. PARSEC/MIST and age cells receive equal total weight; rows within a cell carry d(log M) quadrature weight. |
| Real-star mass posterior draws | 256 per family/R_V branch | D | Stratified propagation of the full A_V grid, all available six-band errors, subgroup age posterior, and unresolved-binary branch. |
| Mass-posterior measure | log-uniform in initial mass | B | Scale-invariant inference measure; no IMF slope is imposed before WP5 tests α=2.0/2.3/2.6. |
| Six-band mass-model width | 0.38 mag per band | C/D | Median robust G/BP/RP/J/H/Ks residual scale from 107 non-extreme spectroscopic-HRD anchors (individual band scales 0.35–0.39 mag). The binary tail is marginalized explicitly. |
| Injection mass-posterior draws | 64 for baseline PARSEC R_V=3.1; 16 for nonbaseline sensitivity branches | D | The baseline was increased from 16 to 64 as a response-convergence check. The remaining CygOB2-C gate failure persisted, so it is not a Monte Carlo-resolution artefact. |
| Permanent WP3 F4 bins | [2,2.52), [2.52,3.17), [3.17,5.04), [5.04,6.35), [6.35,8), [8,12) M_sun | B | Frozen before final evaluation and aligned with the original smoking-gun/mass-function boundaries. |
| Permanent WP3 F4 amplitude | max absolute bin-median ΔA_V < 0.30 mag | B | Predeclared in the repair brief. `repair_v1` gives 0.083 mag. |
| Permanent WP3 F4 rank test | Spearman test on the six bin medians, p >= 0.05 | B | Tests a monotonic mass trend without letting N≈10^3 make negligible individual-star correlations formally significant. `repair_v1` gives p=0.329. |

The WP5 Poisson-residual gate is **unchanged**: chi-square p ≥ 0.01,
residual-trend p ≥ 0.05, and max absolute Pearson residual ≤ 3 in every
subgroup. The converged baseline repair still fails only the last condition for
CygOB2-C (3.257), so WP6 remains unauthorized.

---

## 13. Branch-retention policy (open issue #5, adopted 2026-07-27)

*Written before the `repair_v4` grid was evaluated, so it cannot be tailored to
which branches happened to fail. Plan §1.4 requires model branches to be
carried and reported, never averaged; this section states what "carried" means
when some of them fail the residual gate.*

The branch grid is 2 isochrone families × 3 R_V × 3 IMF slopes = 18 model
branches, each fit for 3 subgroups (54 subgroup-branch fits).

**The rules.**

1. **No branch is ever dropped.** Every one of the 54 fits is reported with its
   gate statistics, whatever the outcome. A failing branch is a measurement of
   a systematic, not a defective run to be discarded.
2. **The baseline branch is the headline and must pass.** PARSEC, R_V = 3.1,
   α = 2.3 is the predeclared baseline (`CUTS_AND_THRESHOLDS.md` §2). WP5
   acceptance requires the baseline to pass for all three subgroups. A version
   whose baseline fails is blocked regardless of how many other branches pass.
3. **Non-baseline failures are carried as a systematic statement**, in the form
   "N/54 pass; the failures are concentrated at <corner>". They do not block
   acceptance by themselves, and they must never be repaired by moving a
   threshold (§6.4).
4. **A concentrated failure pattern is a result and must be reported as one.**
   If failures cluster in a corner of the grid (historically R_V = 3.5 and
   α = 2.6), that clustering is evidence about the extinction law or the slope
   prior and belongs in the systematics discussion, not in a footnote.
5. **Downstream consumes the ensemble, not the winner.** WP6/WP7 propagate
   every branch that passes its own gate and report each derived quantity with
   its across-branch spread. Selecting the single best-scoring branch for a
   headline number is the same error as tuning a threshold.
6. **A failing branch may not be silently re-run to a better outcome.** Any
   re-run is a new version with its own provenance record; the earlier result
   stays on disk and referenced.

**Consequence for acceptance.** The WP5 gate reads: baseline passes for all
three subgroups **and** the association mass is within a factor two of the
literature scale **and** no previously passing subgroup-branch regressed. The
all-54 pass rate is reported, not required.

---

## 14. Replacement of the WP5 residual-trend statistic (pre-declared 2026-07-28)

*Written and committed **before** the replacement was evaluated on any pipeline
version, so it cannot be reverse-engineered to a desired verdict. This is a
§6.4 action — **replace the diagnostic, do not move the threshold**. The
threshold stays 0.05 and the statistic keeps the same role in the same
three-way conjunction.*

### 14.1 Why the incumbent is being replaced

The incumbent trend test is a Spearman rank correlation between
log₁₀(bin centre) and the Pearson residual over the 6 observed-mass bins. It is
**not stable against Monte-Carlo resampling of an unchanged model** (open issue
#11). Two properties combine:

- It responds to the residuals' *rank order*, not their magnitude, so noise in
  residuals that are individually consistent with zero can reorder them and
  swing the statistic.
- With n = 6 the two-sided p-value is confined to a lattice — near the
  threshold the only achievable values are 0.0048, 0.0188, 0.0416, 0.0724.
  **No p-value exists between 0.042 and 0.072**, so a cell in that region has
  no stable verdict; it sits on a cliff.

Measured consequence: of the 7 gate flips between `repair_v3` and `repair_v4`,
4 occurred in cells whose estimator is *provably identical* between the two
versions, and a paired refit on two independent Monte-Carlo realizations of the
same model returns opposite verdicts (2 flips in each direction).
Evidence: `provenance/wp5_trend_stability_check_execution.json`.

Exact null over all 720 orderings: false-positive rate 0.0583 at nominal 0.05.

### 14.2 The replacement statistic

For each subgroup × family × R_V × α cell, with x_i = log₁₀ of the geometric
centre of bin i and res_i the Pearson residual:

    b    = Σ (x_i − x̄) res_i / Σ (x_i − x̄)²        (least-squares slope)
    T    = b · sqrt( Σ (x_i − x̄)² )                  (test statistic)

Pearson residuals are already standardized, so T needs no estimate of the
residual variance from 6 points. T measures **how large** the mass-dependent
drift is, not merely how the bins rank, and is continuous in the residuals.

### 14.3 Null distribution — parametric bootstrap

The p-value is **not** taken from an asymptotic formula. For each cell, with
λ_i the fitted expected counts and r_i the per-k rates:

1. simulate ñ_i ~ Poisson(λ_i), M = 20,000 times;
2. refit the normalization on each simulated dataset by the same Jeffreys rule
   used in production, k̃ = median of Gamma(Σñ_i + ½, rate Σr_i);
3. recompute res̃_i = (ñ_i − k̃ r_i)/sqrt(k̃ r_i) and hence T̃;
4. p = fraction of draws with |T̃| ≥ |T_observed|.

This is exact for the actual bin counts: it carries Poisson non-normality in
low-count bins and the constraint imposed by fitting k from the same counts,
neither of which an asymptotic null handles. Seeded deterministically and
recorded in provenance.

**Gate:** `trend_p ≥ 0.05`, unchanged, inside the unchanged conjunction
`chi2_p ≥ 0.01` **and** `trend_p ≥ 0.05` **and** `max|Pearson residual| ≤ 3.0`.

### 14.4 Acceptance criteria for the replacement itself

The replacement is adopted **only if all four hold**. If any fails it is
rejected and the incumbent stands.

| # | Criterion | Requirement |
|---|---|---|
| R1 | **Calibration** | False-positive rate on simulated null data drawn from the real fitted λ of every cell lies in [0.04, 0.06] at nominal 0.05. Reported alongside the incumbent's rate on the same data. |
| R2 | **Power** | Against residual drifts of injected known slope, detection rate ≥ the incumbent's **at matched false-positive rate**. A test that merely passes more cells is a weakened gate, not a better one, and is rejected under this criterion. |
| R3 | **Stability** | On the four identical-model cell pairs of issue #11, both Monte-Carlo realizations must return the **same** verdict. This is the specific defect being repaired. |
| R4 | **Equal footing** | Every version — frozen, repair_v1, v2, v3, v4 — is re-scored and its grid change reported. A replacement that improves only the newest version is tuning in disguise and must be reported as such. |

### 14.5 Binding consequence

Gate G3 for `repair_v4` is then re-evaluated under the replaced statistic with
the **same strict per-branch reading** already adopted. The outcome is binding
in both directions: if `repair_v4` fails under a correctly calibrated test,
that is the result and WP5 stays blocked.

### 14.6 Criterion R3 re-specified (pre-declared 2026-07-28, before measurement)

R3 in §14.4 required that on the four identical-model cell pairs of issue #11,
**both** Monte-Carlo realizations return the same verdict. That criterion is
**wrong in principle, not merely strict**, and it is withdrawn.

The reason is elementary and applies to every statistic, not just this one. For
a cell whose model gives it probability π of passing, two independent injection
realizations disagree with probability **2π(1 − π)**. A test that never
disagrees would need π ∈ {0, 1} for every cell — that is, its verdict would have
to be a deterministic function of the model with no sampling variation at all.
No test of finite data has that property. R3 as written was therefore
unsatisfiable by *any* candidate, which is precisely why the replacement failed
it on the one cell whose true p sits at ≈ 0.05. The failure was in the
criterion.

What is genuinely wrong with the incumbent is narrower, and R3 is re-specified
to test exactly that and nothing more. **These three replace R3; R1, R2 and R4
are unchanged and their recorded outcomes stand.**

| # | Criterion | Requirement |
|---|---|---|
| **R3a** | **No forbidden region at the threshold** | Simulating null data from every cell's own fitted λ, the widest gap between achievable p-values inside [0.01, 0.20] must be **≤ 0.005**. The incumbent's lattice has no achievable value between **0.0416 and 0.0724** — a gap of 0.031 straddling the 0.05 threshold — so a cell landing there has no verdict the data can support. This is the defect being repaired. |
| **R3b** | **Flip rate consistent with calibration, not zero** | On the four identical-model pairs, the observed number of verdict flips must lie inside the 95% interval of the Poisson-binomial implied by each cell's own pass-probability π. Both too many flips (unstable) and too few (degenerate) fail. |
| **R3c** | **Indeterminacy declared, not hidden** | Every cell carries a computed pass-probability π under its own injection uncertainty. Cells with **0.05 < π < 0.95** are labelled **indeterminate** and reported as such — never counted as a clean pass or a clean failure. |

**How π is computed.** The noise source measured in issue #11 is the injection
realization, and the fit already carries a model of it: the Dirichlet posterior
over the response matrix's category counts. For each cell, M = 400 response
replicates are drawn from that same posterior, the normalization is refitted on
each by the production Jeffreys rule against the **unchanged** observed counts,
and the full three-way gate is re-evaluated. π is the fraction that pass. No new
parameter and no new noise model: this is the fit's own response uncertainty,
already propagated into every published k.

### 14.7 Binding consequence of R3c for gate G3

Under the strict per-branch reading of G3 adopted 2026-07-27, any CygOB2-A or
CygOB2-C cell moving pass → fail blocks acceptance. R3c refines *what counts as
a move*, and the refinement is declared here **before** the version it will be
applied to (`repair_v6`) exists:

1. A cell that is **indeterminate in both versions** contributes to neither a
   regression nor an improvement. Its verdict is a coin flip on the injection
   realization, so reading it either way is reading noise.
2. The rule is **symmetric**: an indeterminate cell flipping fail → pass is
   likewise not counted as an improvement, and the grid count is reported with
   indeterminate cells broken out separately.
3. The rule applies to **every version**, retroactively, not only the newest.
4. A regression that is **determinate in both versions** blocks, exactly as
   before. The strict reading is untouched for cells whose verdict the data
   actually supports.
5. If more than **25%** of the 54-cell grid is indeterminate, that is itself a
   finding — the branch grid is underpowered — and must be reported as such
   rather than absorbed.

**This is a change to a reading the project owner chose explicitly on
2026-07-27.** Every G3 evaluation from here reports the outcome under **both**
readings, and WP6 authorization is not flipped on the refined reading alone
without that difference being stated in the completion report.

---

## 15. Adopted WP3 anchor prior: kriged mean (`repair_v5`, 2026-07-28)

*Adopted after the pre-registered test in §14's spirit — predictions written
before the run, all four confirmed. Decision record:
`provenance/wp3_kriging_adoption.json`; evidence:
`reports/WP3_KRIGED_PRIOR_repair_v5.md`.*

### 15.1 What changed

| Quantity | Adopted value | Class | Justification |
|---|---|---|---|
| Anchor-prior **mean** | simple-kriging estimate over the 8 nearest anchors, shrinking to the anchor field median as separation grows | **A (derived)** | The previous plain neighbour median gave distant anchors weights summing to one while `prior_sigma_at` simultaneously widened the same star's uncertainty toward the sill — the prior's first and second moments contradicted each other. Kriging removes the contradiction using the variogram already fitted in repair_v3. |
| Variogram parameters | nugget ≈ 0, sill 1.228 mag, range 0.853° | C | **Unchanged** — fitted in repair_v3 and reused as-is. |
| Anchor-prior **width** | unchanged (unconditional variogram sigma, *not* the smaller kriging variance) | C | Deliberately conservative: the change corrects where the prior sits, not how confident it is. |
| Kriging jitter | 1e-6 added to the 8×8 covariance diagonal | D | Numerical only; the fitted nugget is ~0 so the matrix is near-singular when anchors nearly coincide. |

**No new free parameter is introduced.** The size of the correction follows
from each star's anchor geometry: kriging weights sum to 1.000 for CygOB2-A,
0.992 for CygOB2-C and 0.772 for CygOB2-B, so the correction is intrinsically
B-specific without any per-subgroup choice having been made.

### 15.2 Why this is not tuning

Under §6.4 the test is whether the change was motivated by the diagnostic or by
the gate. Four things establish the former:

1. The defect is an internal inconsistency identifiable **without reference to
   any gate outcome** — the prior's width and mean disagreed.
2. The predictions were **pre-declared** (`wp3_kriging_prior_prereg.json`) with
   explicit falsification conditions, including one — P4 — that could have
   failed and was committed in advance to mean "the extinction error is real
   but is not the cause, and CygOB2-B being genuinely distinct returns as the
   leading hypothesis."
3. Two **independent** 3D dust maps rank CygOB2-B as the least extinguished
   subgroup, contradicting the uncorrected prior's ranking of it as the most.
4. The correction's magnitude is fixed by the fitted variogram, not chosen.

### 15.3 Gate status

The WP5 residual gate is **unchanged**: χ² p ≥ 0.01 **and** trend p ≥ 0.05
**and** max |Pearson residual| ≤ 3.0, in every subgroup. `repair_v5` passes the
baseline for all three subgroups under **both** the incumbent rank statistic
and the §14 replacement statistic — the first version to do so — with the
branch grid at 38/54. It nevertheless **does not** pass gate G3 under the
strict per-branch reading, because three A/C non-baseline cells regress; two of
those sit at trend p = 0.040 and 0.048, inside the indeterminate band of issue
#11. **WP5 is not accepted and WP6 is not authorized.**

### 15.4 Systematics this makes explicit

Broadband photometry sits about **0.5 mag below** spectroscopically calibrated
anchors at matched sky position (+0.500 mag, Wilcoxon p = 4.3×10⁻¹⁶, measured
against the members' own prior-free photometric A_V). The project's absolute
extinction scale — and therefore its absolute mass scale and N_SN — rests on
the anchor calibration rather than on the photometry. This predates `repair_v5`
and is unchanged by it, but it belongs in the systematics budget.

---

## 16. WP6 census-closure estimator (binding, 2026-07-28)

*Resolves open issue #3, which had flagged that WP6 step 2(a) assumes bright-mass
completeness ≈ 1.0. Evidence:
`provenance/wp6_bright_completeness_execution.json`,
`scripts/wp6_bright_completeness.py`.*

### 16.1 The defect

WP6's closure test compares the observed massive-star count against the number
the WP5 normalization predicts, and reads the shortfall as the population that
has already exploded. As designed it compares against the **true** number. The
catalogue does not contain every massive star: the injection experiment recovers
a fraction that is **flat in mass out to 18 M☉** (plateau 0.787 median across
the grid), so the loss is the WP2 quality filter, not a magnitude limit, and it
does not disappear at the bright end. Uncorrected, that fraction is read as a
real deficit and inflates N_SN by the same factor.

### 16.2 The estimator WP6 must use

With `k` the WP5 normalization, α the IMF slope, and `R(observed above
threshold | M)` the injection response:

    expected observed count above threshold
        = k · ∫[M_floor, M_turnoff] dM M^(−α) · R(observed above threshold | M)

This is the construction WP5 already uses inside 2–8 M☉, extended to one
open-ended bin. The deficit is `observed − expected_observed`, and **only that
difference** is attributable to stars that have left the main sequence.

> **The lower bound is an integration bound, not the supernova threshold**
> (issue #17, added 2026-07-28 after the original implementation set them
> equal). The observed side counts every member's `P(M > 8)` *whatever that
> member's true mass is*, so the predicted side must integrate over every true
> mass the response can push above 8 M☉. Setting `M_floor = 8` omits that
> up-scatter entirely and makes the two sides different quantities — it is the
> same error this section forbids one paragraph below, committed at the other
> end of the integral. **`M_floor` is driven to convergence, never chosen for
> its effect**: the adopted 4.0 M☉ leaves the integral changing by ~1.4% per
> further 1 M☉ of extension, and the residual is carried as a systematic.
> Measured up-scatter: R(estimated > 8 | true M) = 0.23 at 7 M☉, 0.09 at 6,
> 0.03 at 5, 0.004 at 4.
>
> `predicted_true_living` keeps its lower bound at 8 M☉, because it counts stars
> *physically* above 8. The two bounds differ on purpose, and their ratio is an
> effective completeness that may exceed 1.

**Dividing the observed count by a scalar completeness is forbidden**, and not
merely as a matter of style. The response also scatters mass estimates across
the threshold, and the two treatments disagree in size and sometimes in sign:

| quantity | value |
|---|---|
| scalar plateau recovery (the naive divisor) | **0.785** |
| forward effective completeness, grid median | **0.872** |
| forward effective completeness, grid range | 0.810 – **1.040** |
| cells where the forward value **exceeds 1** | **6 of 54** |

In those 6 cells net up-scatter across 8 M☉ more than compensates the recovery
loss, so the scalar correction has the **wrong sign**. Dividing by 0.787
everywhere would over-correct by ~10% and turn a modest real deficit into a
manufactured one in the opposite direction.

### 16.3 Size of the bias avoided

A closure test assuming completeness 1.0 reports a spurious deficit of **12.8%**
(grid median) and inflates N_SN by that factor. Note this is materially smaller
than the ~20% that issue #3 originally estimated from the raw plateau — the
up-scatter partly cancels it, and only the forward calculation reveals that.

### 16.4 Per-subgroup systematic that WP6 must carry

Effective completeness is **not** common across subgroups:

| subgroup | min | median | max |
|---|---:|---:|---:|
| CygOB2-A | 0.810 | 0.834 | 0.872 |
| **CygOB2-B** | 0.912 | **0.962** | 1.040 |
| CygOB2-C | 0.822 | 0.846 | 0.921 |

CygOB2-B sits ~0.12 higher because its mass posteriors are wider — the same
weaker extinction calibration recorded as obligation O3 in §15 — so more of its
stars scatter up across the threshold. A single association-wide completeness
would bias B's closure relative to A's and C's, in the direction of
under-counting B's explosions. **WP6 must apply the response per subgroup.**

### 16.5 Re-measurement obligation

The numbers above are measured on **`repair_v6`, the accepted version**
(2026-07-28); they moved by <0.01 from the `repair_v5` measurement, so the
specification is insensitive to which of the two is used. Should any later WP5
version supersede `repair_v6`, `scripts/wp6_bright_completeness.py --version
<accepted>` must be re-run and its output re-read into WP6 before the closure
test is executed.

---

## 17. Runaway traceback — peculiar motion is binding (issue #16, 2026-07-28)

### 17.1 The rule

**A traceback must run on peculiar proper motions, with the association's
systemic motion subtracted.** Absolute Gaia proper motions are dominated by the
association's own bulk motion plus Galactic rotation. For Cyg OB2 the systemic
motion is **(−2.707, −4.317) mas/yr**, which is *larger than a typical ejection
signature* — so an absolute-PM traceback measures the common drift and not the
ejection.

The subtraction is rotated into galactic coordinates **at each star's own
position**, not applied as a constant offset in (l, b): over an ±8° box the
rotation between equatorial and galactic axes varies enough to matter.

**Control fields receive identical treatment.** If controls were traced back on
absolute motions while candidates used peculiar ones, the measured
false-positive rate would not describe the candidate sample.

### 17.2 The false-positive rate is measured per separation bin

Chance recovery falls steeply with angular separation from the centroid — 20.7%
at 1–1.5°, 0.4% beyond 8°. A pooled control rate is therefore meaningless: the
control and candidate samples do not share a separation distribution. **Applying
a pooled rate produced a negative corrected count.** The rate must be measured in
bins of separation (edges 1, 1.5, 2, 3, 4, 5, 6, 8, 12°) and applied to each
candidate at its own separation.

Purity per bin is `(n_selected − expected_chance) / n_selected`, with the
expected chance count summed over **all** candidates in the bin, not the selected
ones. Bins are clipped at zero, so the binned total can exceed the unclipped
aggregate; **both must be reported** and WP7 carries the difference as a
systematic.

### 17.3 The external gate is not optional

The runaway result is gated on recovering published candidates
(`scripts/wp6_runaway_crossmatch.py`). This is what caught issue #16: every
internal diagnostic passed while BD+43 3654 scored exactly 0.000.

Two rules follow.

1. **Positions are frozen into the script** with the identifier queried and the
   query date, so the check reruns offline and reproducibly.
2. **In-footprint literature candidates cannot pass or fail the gate.** Most
   published "runaway candidates" are stars still inside the association with
   anomalous proper motions; an escaped-star search excludes them by
   construction. Counting them as gate failures would be wrong, and counting
   them as passes would be worse.

---

## 18. Massive-star multiplicity in the injection truth model (issue #15, pre-declared 2026-07-28)

Full pre-registration, with the M1/M2/M3 predictions and the decision rule:
`provenance/wp6_multiplicity_prereg.json`, written **before** any multiplicity
injection was generated.

### 18.0 A correction to master-table row 20

Row 20 of §2 carries the binary fraction as a Class E branch spanning **0.3–0.5**
and marks it "correct — keep". Issue #15 shows that verdict is wrong at the top
end, in two separate ways.

1. **The range does not span the measured value.** Sana+2012 give f_b ≈ 0.70 for
   15–60 M☉ and Caballero-Nieves+2020 find 48 of 74 Cyg OB2 O/early-B stars
   multiple. **0.70 lies outside 0.3–0.5 entirely**, so no branch in the carried
   grid represents the massive-star regime. The range was specified from
   field/low-mass literature and then applied to the whole mass axis.
2. **WP5 does not carry the branch at all.** WP4 genuinely fits ages at
   f_bin ∈ {0.3, 0.4, 0.5}, but the WP5 injections select `f_bin == F_BINARY`
   and inject a **single constant 0.40 at every mass**. On the response side the
   Class E branch is not carried; it is fixed.

Neither point affects WP5's accepted result, whose calibration window stops at
8 M☉ where 0.40 is defensible. Both bear directly on WP6, whose closure window
starts exactly where the specification stops being right.

**Consequence for row 20**: the correct reading is *Class E, and the carried
range is under-specified above ~10 M☉*. §18 below specifies the replacement.

### 18.1 The truth-side model

| mass | truth f_bin |
|---|---:|
| ≤ 8 M☉ | 0.400 |
| 10 M☉ | 0.497 |
| 12 M☉ | 0.575 |
| ≥ 16 M☉ | 0.700 |

Linear in log mass between the anchors, flat outside. **Both anchors are
measured values** — 0.40 is the frozen WP5 value where the 2–8 M☉ calibration
lives, 0.70 is Sana+2012 for 15–60 M☉ — so no free parameter is tuned to make
the closure ratio move.

`q ~ Uniform(0.1, 1)` is **unchanged**: Sana+2012's mass-ratio exponent
κ = −0.1 ± 0.6 already supports it.

### 18.2 Only the truth side changes

The recovery side keeps assuming f_bin = 0.40. This is deliberate and is the
whole point: nature makes binaries at the true rate while the estimator assumes
0.40, and **that mismatch is the bias under test**. Changing both sides would
measure nothing.

### 18.3 The comparison must be paired

The published M ≥ 8 response is assembled from two files generated on different
mass grids, so a grid-restricted re-run draws a different RNG stream even at the
same seed. **Comparing a treatment run directly against the published ratio
would mix the multiplicity change with a change of Monte Carlo realization.**

Each node is therefore injected twice on the identical grid, control and
treatment, from a fresh `default_rng(SEED)`. The per-star binary threshold
consumes the same single `rng.random(n_injected)` draw as the constant version,
so donor, extinction, photometric and QMC realizations are bit-identical between
the arms.

Both readings are reported and **never averaged**: the pre-registered absolute
threshold binds; the paired reading is reported alongside so the size of the
realization noise is visible rather than hidden.

### 18.4 Bit-preservation is verified, not asserted (V4)

`truth_binary_fraction=None` must leave `inject_curve` byte-identical. This is
checked by regenerating a published WP6 extension node with the current code and
comparing sha256 against the stored artifact. **If V4 fails, nothing else runs.**

### 18.5 Scope — this is a diagnostic, not a new version

The test measures the size of the multiplicity effect on the WP6 closure window
while holding accepted WP5 `repair_v6` fixed. **Adopting** a mass-dependent
f_bin would also perturb the WP5 response inside 2–8 M☉ through down-scatter
from above, and therefore requires a full `repair_v7` chain re-run. Masses below
8 M☉ are deliberately not re-injected here.

### 18.6 Outcome (measured 2026-07-28)

Issue #15 is **resolved**. 162 paired nodes, 3.3 h of injection, realized
f_bin control 0.396 ± 0.005 and treatment 0.611 ± 0.011.

| | grid median at α = 2.3 |
|---|---:|
| published | 1.105 |
| control arm | 1.103 |
| treatment arm | **1.099** |
| **multiplicity effect** | **0.004 — 3.7% of the excess** |

**M1 PASS** (all subgroups fell; 48 of 54 cells). **M2 FAIL** — CygOB2-C's
reduction is *smaller* than B's by ~2.7σ, the opposite of the predicted turnoff
ordering, so the pre-registered consequence applies: the excess and the
multiplicity effect have different mass dependences. **M3 PASS literally, FAIL on
the governing relative reading** — its 1.222 threshold was derived from a 1.444
baseline that §16.2's issue #17 correction has withdrawn, and the control arm
already sits below it before the treatment changes anything.

**Binding consequence, replacing the previous one:** WP6's residual excess is
reported **with the measured 3.7% multiplicity correction applied** and the
remainder carried as a systematic. **No `repair_v7` is triggered by issue #15.**
The disfavouring of α = 2.6 stands.

**The scope limit is binding on how this may be quoted.** The test covers
multiplicity **above 8 M☉ only** — the range where the mechanism has least room
to act, since a star already above the threshold cannot be scattered into the
census by being brightened. The 4–8 M☉ up-scatter channel is held at
f_bin = 0.40 in both arms. Permitted claim: *"multiplicity above 8 M☉ does not
explain the excess."* **Not permitted:** *"multiplicity does not explain the
excess."*

### 18.7 The sub-8 M☉ channel is not negligible (measured 2026-07-29)

A pre-registered go/no-go discriminator
(`provenance/wp5_fbin_discriminator_prereg.json`) extended f_bin below 8 M☉
(0.40 at 2 M☉ → 0.55 at 8 → 0.70 at 16) on three paired nodes of the reporting
branch:

| | measured | pre-declared threshold |
|---|---:|---:|
| **D1** — recovery over the 2–8 M☉ calibration window | **+0.16% ± 0.09%** | 2% |
| **D2** — R(estimated > 8 M☉ \| M) over 4–8 M☉ | **+9.87% ± 1.58%** | 2% |

**`repair_v7` is JUSTIFIED.** The control arm reproduced the accepted `repair_v6`
node byte for byte, so the treatment shift is the model change and nothing else.

**The lesson generalizes beyond f_bin.** The calibration window is insensitive
because most stars there are recovered either way; the 8 M☉ boundary is highly
sensitive because a small brightness shift converts directly into a crossing
probability. **Any truth-model error will therefore show up ~25× more strongly
at a threshold than in a recovery fraction.** A diagnostic that checks only
recovery inside the fitting window will systematically under-report truth-model
error, which is precisely how this one was nearly missed.

This is the same structural lesson as issues #3 and #17: **thresholds are where
the response does its damage**, and they must be probed directly rather than
through window-averaged quantities.
