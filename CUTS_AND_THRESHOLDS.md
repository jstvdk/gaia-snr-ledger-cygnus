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
| 20 | Binary fraction | 0.3–0.5 | E | correct | keep |
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
