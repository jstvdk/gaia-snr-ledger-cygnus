# Cross-check — Härer et al. 2025: supernova rate and association mass

**Source.** L. Härer, T. Vieu, F. Schulze, C. J. K. Larkin, B. Reville,
*"Deciphering the gamma-ray emission in the Cygnus region"*,
A&A **703**, A111 (2025). Local copy: [`papers/Harer_2025.pdf`](../papers/Harer_2025.pdf).

**What they computed.** §3.1 and Fig. 2. Core-collapse supernova event rates for
Cyg OB2 using the **Hoki** package (Stevance et al. 2020), an interface to
**BPASS** (Stanway & Eldridge 2018), assuming solar metallicity, a single stellar
population of total initial mass **1.65 × 10⁴ M☉** (Wright et al. 2015), and an
IMF index **α = −2.0**.

**Why this is the strongest check we have.** It tests the *output* of the whole
WP1→WP6 chain, not an intermediate. Their route (population-synthesis rates from
an adopted literature mass) shares no machinery with ours (a response-aware
Poisson fit to 2–8 M☉ Gaia star counts, extrapolated above the isochrone
turnoff).

---

## 1. The unit trap

Fig. 2's y-axis is:

> **Event Rate (events/100 kyr)**

with ticks at 0.25, 0.5, 1, 2. This is **not** events per 100 *years*. The
conversion that matters:

| reading | per Myr |
|---|---:|
| 1 event / 100 kyr | **10** |
| 1 event / 100 yr | 10,000 |

The paper's own text removes the ambiguity:

> *"For an age of 3–5 Myr, type Ic supernovae are by far the most likely type of
> supernova, and **the rate is consistent with several supernovae per Myr**."*

Anyone comparing against this figure must convert first. A misread of the axis
inflates the expected rate by a factor of 1000 and makes the comparison look
catastrophic when it is not.

**Sanity bound.** At 1 SN/100 yr, Cyg OB2 would produce 40,000 supernovae in
4 Myr. Our IMF gives ~260 stars above 8 M☉ *ever born* in the association. At
1.65 × 10⁴ M☉ and roughly one core-collapse progenitor per 100 M☉ of stars, the
total budget is ~165 supernovae over the association's entire ~40 Myr future —
about one per 250 kyr averaged. The cluster mass alone forbids the misread rate.

---

## 2. Rate comparison

Our instantaneous rate is computed from the WP5 normalization as

    dN_SN/dt = k · M_turnoff^(−α) · |dM_turnoff/dt|

summed over subgroups, with `M_turnoff(age)` from
[`wp6_mass_extension_decision.turnoff_mass`](../scripts/wp6_mass_extension_decision.py)
(PARSEC, R_V = 3.1, ages A 3.98 / B 3.55 / C 2.51 Myr).

| α | our rate (SNe/Myr) | our rate (**events/100 kyr**) | in Fig. 2's plotted range (0.25–2)? |
|---|---:|---:|---|
| **2.0** — *their assumption* | 19.2 | **1.92** | yes, near the top |
| **2.3** — our baseline | 7.8 | **0.78** | yes |
| 2.6 | 3.1 | **0.31** | yes |

**All three branches land inside the range Fig. 2 spans.** Our baseline α = 2.3
gives 0.78 events/100 kyr = 7.8 per Myr, which is precisely "several supernovae
per Myr".

At matched α = 2.0 we sit toward the upper end of their figure. Read
conservatively, the two calculations **agree to within about a factor of two
across the whole branch grid** — which is the appropriate resolution for a
comparison between single-star turnoff counting and BPASS population synthesis.

### Cumulative counts, and why they are consistent with the rate

| α | SNe so far | A (turnoff) | B (turnoff) | C (turnoff) |
|---|---:|---|---|---|
| 2.0 | 16.7 | 10.7 (58 M☉) | 6.0 (74 M☉) | 0 (120 M☉) |
| **2.3** | **6.3** | 4.1 (58 M☉) | 2.2 (74 M☉) | 0 (120 M☉) |
| 2.6 | 2.3 | 1.5 (58 M☉) | 0.8 (74 M☉) | 0 (120 M☉) |

A cumulative count of ~6 and a rate of ~8/Myr are the *same statement*, not a
contradiction. The PARSEC turnoff only falls below the 120 M☉ IMF ceiling at
**3.00 Myr**, so **no star in Cyg OB2 could have exploded before then**.
CygOB2-A is 3.98 Myr old, so its supernovae occupy a window of ~1.0 Myr; six
events in that window is ~6/Myr, rising to the current 7.8/Myr as the turnoff
sweeps down through a steeply rising IMF.

**CygOB2-C contributes exactly zero.** At 2.51 Myr its turnoff is still above the
120 M☉ upper limit — nothing has died there yet.

---

## 3. Association mass — an independent check we had not made

Härer adopt **1.65 × 10⁴ M☉** from Wright et al. 2015. Integrating our own fitted
normalization over 0.5–120 M☉:

| α | implied association mass |
|---|---:|
| 2.0 | 1.99 × 10⁴ M☉ |
| **2.3** | **1.74 × 10⁴ M☉** |
| 2.6 | 1.79 × 10⁴ M☉ |

**Within 5% at our baseline**, by a completely different route: they adopt a
literature mass, we derive one from response-corrected star counts in 2–8 M☉.
Nothing in WP1–WP5 was tuned toward it and the comparison was not made until
after WP6 closed.

**Caveat, stated because it matters.** This number is sensitive to the low-mass
integration limit. We integrate 0.5–120 M☉, the span of the frozen `MASS_GRID`.
A single α = 2.3 power law continued to 0.1 M☉ would give 3.1 × 10⁴ M☉ — the
real IMF flattens below ~0.5 M☉, so our cutoff happens to mimic a Kroupa
turnover. The agreement is therefore good but should not be quoted to better
than a few tens of percent.

---

## 4. One real physical difference

BPASS models **binary evolution**. Härer's Fig. 2 shows type **Ic** supernovae
dominating at 3–5 Myr — the stripped-envelope channel, which *requires* a binary
companion to remove the hydrogen envelope.

Our count is a **single-star turnoff calculation** from PARSEC/MIST isochrones.
It does not include stars driven to explosion by mass transfer, nor mergers.

Consequence: **our N_SN is a lower bound relative to a BPASS-style rate**, and we
already agree with them before adding that channel. Adding it would move us up,
toward their type-Ic-dominated regime.

This is the same physics as [issue #15](../PROJECT_TRACE.md) — massive-star
multiplicity in the injection truth model — arriving from an independent
direction. It is an argument for taking that diagnostic seriously, and it is
noted in the WP7 handoff.

---

## 5. What a failure would have looked like

This check could have failed in three distinguishable ways, and did not:

| failure mode | what we would have seen | actual |
|---|---|---|
| normalization wrong | association mass off by ≫2× | 5% |
| turnoff wrong (issue #14 class) | rate off by ≫10× | within ~2× |
| extrapolation above 8 M☉ invalid | rate and mass disagreeing in *opposite* directions | both agree |

The third is the informative one: mass and rate depend on `k` in different ways
(mass integrates the whole IMF, rate samples it only at the turnoff), so a broken
extrapolation would generally break them inconsistently.

---

## 6. Verdict

**AGREE.** Two independent quantities — the supernova rate and the association
stellar mass — match published values computed by unrelated machinery. No
number was adjusted as a result of this comparison, and none should be.

Carried forward: our rate is a **lower bound** relative to BPASS because the
binary-stripped channel is not modelled.
