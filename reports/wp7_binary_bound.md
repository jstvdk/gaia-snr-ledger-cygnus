# T3 — binary mass transfer, converted from "acknowledged" to "bounded"

*Item **T3** of [pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
Machine-readable record:
[wp7_binary_bound_execution.json](../provenance/wp7_binary_bound_execution.json).
Tables: [wp7_binary_bound.csv](../tables/wp7_binary_bound.csv) ·
[wp7_binary_bound_harer_fig2.csv](../tables/wp7_binary_bound_harer_fig2.csv) ·
[wp7_binary_bound_branches.csv](../tables/wp7_binary_bound_branches.csv).
Script: [wp7_binary_bound.py](../scripts/wp7_binary_bound.py).*

> **Scope, stated first.** This is a labelled, literature-scaled bound for the
> discussion section. It is **not** a pipeline change: no injection was run, no
> fit was refitted, no stored WP5–WP9 product was touched, and the bracket is
> never marginalized into the headline number.

---

## 1. The vulnerability

WP7 counts supernovae by single-star turnoff: a star has died if its ZAMS mass
exceeds the turnoff at its birth epoch. Binary interaction is acknowledged
throughout this project and modelled nowhere. Three independent places point at
the same gap:

- **issue #15** measured a massive-star binary fraction of ≈ 0.7 above 8 M☉
  (Sana+2012; Caballero-Nieves+2020 find 48 of 74 Cyg OB2 O/early-B stars
  multiple), against the truth model's constant 0.40;
- **Härer's own comparison** is BPASS-based and type-Ic-dominated — a channel
  that requires a companion;
- **reading (c)** of the WP8 pulsar degeneracy is binary stripping.

A referee will ask for a magnitude. Here is one, from two independent
directions.

## 2. Line 1 — the empirical bound: BPASS versus us, matched

Härer et al. 2025 Fig. 2 is a BPASS calculation **with binaries**, for a
Cyg OB2-like population: a core-collapse event rate versus age for a cluster of
total initial mass 1.65 × 10⁴ M☉ at IMF index α = −2.0. Our single-star rate at
the **same α, the same age and the same population mass** is its like-for-like
partner, so the ratio measures the binary correction end-to-end, on exactly the
population in question.

Fig. 2 is vector art in the local PDF, so it was **digitized exactly** rather
than read by eye. Two features confirm the calibration: the x axis is
logarithmic (the 2.5 / 5 / 7.5 / 10 tick spacing is 1 : 0.585 : 0.415, which no
linear axis produces), and the resulting step edges land on **BPASS's native
0.1-dex log-age grid** — 2.512, 3.162, 3.981, 5.012, 6.310, 7.944, 10.0 Myr — to
better than 0.5%. The y calibration reproduces the 0.25 / 0.5 / 1 ticks exactly.

| age bin (Myr) | total CCSN rate (events/100 kyr) |
|---|---:|
| 2.512 – 3.162 | 1.373 |
| 3.162 – 3.981 | 1.011 |
| **3.981 – 5.012** | **0.990** |
| 5.012 – 6.310 | 1.028 |
| 6.310 – 7.944 | 1.021 |
| 7.944 – 10.00 | 0.950 |

The comparison at α = 2.0 (their assumption), per 10⁴ M☉ of population:

| | age | turnoff | our rate | BPASS rate | ratio |
|---|---:|---:|---:|---:|---:|
| CygOB2-A | 4.01 Myr | 57.5 M☉ | 9.28 SNe/Myr | 6.00 | 0.647 |
| CygOB2-B | 4.12 Myr | 55.1 M☉ | 9.45 SNe/Myr | 6.00 | 0.635 |

One bookkeeping correction is required before these are comparable as physics.
BPASS's population is defined over **0.1–300 M☉** and ours over 0.08–120 M☉, so
for a fixed total mass their normalization buys **1.130×** fewer stars near the
turnoff. Applying it:

> **corrected ratio BPASS(binary) / ours(single-star) = 0.72 – 0.73.**

**The binary-inclusive calculation sits ~28% BELOW our single-star count, not
above it.** Whatever the unmodelled binary physics does to N_SN at this age, it
is not an order-unity enhancement — and its sign is not even established as
positive.

Three residuals are *not* corrected and all push the ratio back toward 1: BPASS's
high-mass multiplicity is near unity while our multiplicity-adjusted mass
assumes f_bin = 0.40, so their total mass buys still fewer primaries; the
digitized rate is a step function, so ages inside a bin share one value; and
BPASS's sub-0.5 M☉ IMF is its own. The honest reading is *agreement at the
few-tens-of-percent level*, not a measurement of 0.72.

*(Only α = 2.0 is used. Comparing our α = 2.3 branch against their α = −2.0
figure would measure the IMF slope, not the binaries; those rows are tabulated
and explicitly unused.)*

## 3. Line 2 — the theoretical bound: delay times

Zapartas et al. 2017 (A&A 601, A29) computed the core-collapse delay-time
distribution *with* binary interaction. Two of their results bound this case
directly:

1. **"the total number of core-collapse supernovae increases by 14⁺¹⁵₋₁₄%
   because of binarity"** — integrated over all delay times;
2. **"a significant fraction, 15⁺⁹₋₈%, of core-collapse supernovae are 'late',
   occurring 50–200 Myr after birth"**, from binaries with one or both stars
   initially of *intermediate* mass (4–8 M☉) — and the binary and single-star
   distributions are **"remarkably similar at early times"**, diverging only
   around 20 Myr.

**Cyg OB2 is 4.0 Myr old and its first supernova was 1.30 Myr ago.** Every event
in the ledger sits at a delay time below ~4.1 Myr — deep inside the regime where
the two distributions agree, and five times younger than the divergence. The
channel that produces the +14% *cannot have operated here at all*. Applying the
integrated +14% to this population is therefore a deliberate **overestimate**.

**Why de Mink et al. 2014's 30% does not scale directly.** That 30⁺¹⁰₋₁₅% of
massive main-sequence stars are binary products (8⁺⁹₋₄% mergers) is a statement
about the *living* population, not the dead one. Accretors and mergers are
**rejuvenated**: a star pushed above the present turnoff by accretion has had its
clock partly reset and has typically *not* yet exploded at 4 Myr. It adds to the
future budget, not the past one — which is the quantity this paper reports.

## 4. The adopted bracket: ±30% on N_SN

| edge | source |
|---|---|
| **+30%** | the top of Zapartas et al.'s 14 + 15% integrated enhancement, rounded up and applied **in full** even though its dominant channel cannot operate at 4 Myr |
| **−30%** | symmetric; covers their −14% error edge and the empirical finding that the BPASS rate sits ~28% below ours |

Applied as a scale factor on `k` and run through WP7's own engine
([wp7_ledger.run_population](../scripts/wp7_ledger.py), 400,000 iterations per
branch, all 36 headline branches):

| arm | baseline N_SN | headline range | P(last SN < 100 kyr), baseline |
|---|---:|---|---:|
| −30% | **5.90** | 3.94 – 20.13 | 0.430 |
| nominal | **8.43** | 5.62 – 28.74 | 0.552 |
| +30% | **10.96** | 7.31 – 37.36 | 0.648 |

## 5. The answer to the referee

> **Unmodelled binary mass transfer moves the baseline supernova count from 5.9
> to 11.0 — a span of 5.1 supernovae, against a headline model-branch span of
> 23.1. The systematic this project does not model is about five times smaller
> than the model-branch spread it already reports, and it changes no
> conclusion.**

On the verdict it is likewise sub-dominant: the bracket moves
P(last SN < 100 kyr) from 0.43 to 0.65 about a nominal 0.55, a swing of 0.22 —
comparable to the α axis (0.264 on P_verdict) but applied to a quantity where
α also dominates, and it does not move the baseline across the 0.5 framing
boundary in the direction that matters. The +30% arm *helps* the Härer
scenario, so the reported verdict is conservative with respect to this
systematic.

## 6. What binaries **do** change — and it is not the count

**The supernova type.** Härer's Fig. 2 is essentially pure type Ic from 2.5 to
6.4 Myr: the stripped-envelope channel. WP9's condition C3 asserts that every
progenitor above ~30 M☉ is stripped and reaches **1.000**.

That is not a hidden dependence on binaries. Binary stripping is one route to an
envelope-free progenitor; **single-star Wolf–Rayet winds are another**, and
Cyg OB2 contains three known WR stars. At the 34–58 M☉ progenitor masses this
ledger actually deals with, wind stripping is the standard expectation. So C3
does not require the binary channel to be present, and the binary channel does
not change the count C3 is applied to.

## 7. A new dependency

The digitization needs **PyMuPDF** (`pip install pymupdf`, 1.28.0 here) to read
the vector paths out of the local PDF. It is imported defensively: if it is
absent the script falls back to the stored constants and says so in the
provenance record rather than failing. Nothing else in the chain depends on it.

## 8. What this bound is not

- It is **not** a BPASS integration. Doing one properly would mean adopting
  BPASS's own IMF, metallicity grid and remnant prescription, which would
  replace this project's measurement rather than bound it. That is Paper 2 work.
- It is **not** marginalized in. The headline N_SN stays 8.43 with its
  36-branch range; the ±30% is quoted as a stated systematic beside it.
- It does **not** cover the 4–8 M☉ up-scatter channel in the WP5 calibration
  window, which issue #15 flagged as untestable without a chain re-run and which
  `repair_v7` subsequently measured at +0.54% on `k` — three orders of magnitude
  below this bracket, and already inside the accepted chain.
