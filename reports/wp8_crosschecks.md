# WP8 — external cross-checks of the supernova ledger

*Pre-registered in [wp8_crosschecks_prereg.json](../provenance/wp8_crosschecks_prereg.json)
before any comparison was scored. Markers frozen at WP1 on 2026-07-22, before
WP5, WP6 and WP7 existed. Machine-readable record:
[wp8_crosschecks_execution.json](../provenance/wp8_crosschecks_execution.json).*

**All five predictions passed. Gate G8a — pulsar consistency — is RESOLVED.**

---

## 1. The headline

**PSR J2032+4127 is a neutron star inside Cyg OB2, and neutron stars are made by
successful explosions.** WP7 left one question dominating everything: on the
all-explode branch the association has produced ~8 supernovae; on any branch
where stars above ~40 M☉ collapse directly to black holes it has produced
**exactly zero**.

| branch | P(≥1 supernova) |
|---|---:|
| all-explode | **0.99965** |
| islands (BH above 40 M☉) | **0.000000** |

The pulsar's existence is incompatible with the second. **An observation we
already had in hand rules out the branch that would have made the entire
supernova budget vanish** — subject to one caveat that cannot be removed here
(§2c).

## 2. The pulsar is not an external object — it is beside our own star

This was found while assembling the inputs, and it sharpens the test
considerably.

**MT91 213, the pulsar's B0V companion, is this project's own anchor**
`gaia_dr3:2067835682818358400`:

| | |
|---|---|
| channel | **orphan anchor** (one of the 27), `counts_in_census = True` |
| spectral type | B0V, **17 M☉** by the B0V mass rule |
| position | 308.054458, +41.456917 — the pulsar is at 308.0546625, +41.456753 |
| subgroup (positional) | **CygOB2-A**, 0.115° from its centroid |
| Gaia astrometry | RUWE 1.07, parallax error 0.016 mas — clean |

Its absence from WP2 membership is a **selection outcome, not a data failure** —
it missed the astrometric cuts, which is exactly what the orphan-anchor channel
exists for. A reader will ask whether the pulsar's companion is in the census.
It is, and it sits in the footprint of the subgroup that produces half the
ledger's supernovae.

### 2a. The progenitor mass argument

MT91 213 is 17 M☉. The neutron star's progenitor **evolved first**, so it began
*more* massive than its companion. If it was coeval with CygOB2-A — 4.00 Myr,
turnoff **57.9 M☉** — it must have exceeded 57.9 M☉.

**That is precisely the regime where the islands prescription predicts a black
hole, not a neutron star.** The pulsar sits exactly on the contested point.

### 2b. It was born here

| | |
|---|---|
| pulsar proper motion | (−2.99, −0.74) mas/yr |
| association systemic | (−2.7067, −4.3168) mas/yr |
| peculiar transverse velocity | **27.6 km/s** |

The systemic reference is the member median — the same value the WP6 Orellana
cross-check validated against their −2.71 ± 0.02.

27.6 km/s is low, and that is what a **bound** binary requires: the system has a
decades-long eccentric orbit that a large natal kick would have disrupted. The
physical expectation and the measurement agree, and both point to the pulsar
having been born in Cyg OB2 rather than drifting in.

### 2c. What the pulsar cannot settle

Three readings produce a neutron star here, and no observable available to this
project separates them:

| | reading | consequence |
|---|---|---|
| **(a)** | stars above ~58 M☉ do make neutron stars | all-explode branch confirmed |
| **(b)** | the progenitor was less massive, hence **older** than CygOB2-A | a population WP4 does not resolve |
| **(c)** | **binary mass transfer** stripped the progenitor, so a lower-mass star died early | the effect WP7 §10 already names as its largest unmodelled physics |

(c) is not a special plea: WP6 measured multiplicity above 8 M☉ at f_bin ≈ 0.7,
so binary interaction is the norm for these stars, not the exception. **The
honest statement is that the pulsar excludes "no supernovae ever", but does not
by itself establish that 60 M☉ stars explode.**

### 2d. The age agrees too

The characteristic age is **200.7 kyr**, and `wp1_sn_markers.md` already recorded
that this is not an explosion-age measurement. Widening it for braking index
n ∈ [2, 3] and birth period P₀/P ∈ [0, 0.5] gives **151–401 kyr**.

The ledger, built without ever looking at the pulsar, gives:

- P(last SN within 401 kyr) = **0.960**
- P(last SN inside 151–401 kyr) = **0.259**

The pulsar's age is a typical draw from the ledger, not a tail event. A further
systematic remains: the measured `Ṗ` carries line-of-sight acceleration from the
decades-long orbit, which is a timing problem beyond this project's data.

## 3. γ Cygni — allowed, but the association stays unsettled

Two-sided, as pre-registered.

**Does the ledger allow it?** At 8.01 SNe/Myr — one per **125 kyr** — the
probability of a supernova within γ Cygni's 6.8–10 kyr age is **7.7%**.
Comfortably non-zero, so the ledger does not forbid it; far from predicting it.
Pre-registered range [0.03, 0.12], **X5 PASS**.

**Is it physically associated?** Our distance is 1.62 kpc; Leahy+2013 infer
1.7–2.6 kpc from H I absorption. These overlap **only at the extreme low end of
theirs**. The older ~1.5 kpc class of estimates is carried as an explicit
literature branch.

**Geometry.** The remnant sits 2.27° from the association centroid — **64 pc**
projected at 1.62 kpc. Outside the 1° footprint, but well inside the distance a
runaway covers in a few Myr at the ~34 km/s that implies. Given WP6 recovered
54.9 runaways and WP7 bounds the out-of-association supernova fraction at
≤14.6%, **a runaway progenitor is a live possibility rather than a contrivance**.

**This report does not resolve the association, and does not resolve it in our
favour.** Recorded as tension T1.

## 4. The absent remnants are weak evidence — and saying so is the result

No catalogued remnant lies within **2.27°** of the association centroid:

| remnant | separation | projected |
|---|---:|---:|
| G078.2+02.1 (γ Cygni) | 2.27° | 64 pc |
| G083.0−00.3 | 3.13° | 89 pc |
| G076.9+01.0 | 3.18° | 90 pc |

The Härer scenario predicts an invisible remnant in a low-density wind-blown
cavity, and the non-detection is consistent with it. **But it is consistent with
ordinary Poisson luck too.** At 8.01 SNe/Myr, the expected number of *currently
visible* remnants is:

| assumed visible lifetime | expected number |
|---|---:|
| 20 kyr | **0.16** |
| 100 kyr (generous, no cavity) | **0.80** |

Below one even at the generous end. **The non-detection cannot discriminate
between the cavity scenario and chance**, and claiming otherwise would be
overreach. Pre-registered as X4, **PASS**.

## 5. ²⁶Al — a consistency band, never a count

The ledger's supernovae alone sustain about **8.4 × 10⁻⁴ M☉** of ²⁶Al in steady
state (8.01 Myr⁻¹ × 1.05 Myr mean lifetime × 10⁻⁴ M☉ per event), roughly **1000×
below** the ~1 M☉ inferred for the whole Cygnus complex.

**This is the expected ordering, not a tension.** The INTEGRAL measurement is
complex-wide rather than Cyg OB2-only, and in a region this young **Wolf-Rayet
winds dominate ²⁶Al production** over supernovae. The measurement constrains a
*combination*, and the pre-registration forbade inverting it into a supernova
count. That prohibition is honoured. **G8c satisfied.**

## 6. Neighbouring associations — coarse, by design

Cyg OB1 and Cyg OB9 are older and less massive than Cyg OB2 but not negligible,
so a cavity supernova need not have come from Cyg OB2. Two of the nine
catalogued remnants in the wide box sit within 5° of the association — a
reminder that the field is not empty. No quantitative budget is claimed.

## 7. Predictions, scored

| | statement | outcome |
|---|---|---|
| **X1** | P(≥1 SN) > 0.5 on the all-explode branch | **PASS** — 0.9997 |
| **X2** | the islands branch assigns exactly zero | **PASS** — 0.000000 |
| **X3** | P(last SN within the widened pulsar age) > 0.2 | **PASS** — 0.960 |
| **X4** | expected visible remnants < 1 | **PASS** — 0.80 |
| **X5** | P(SN within γ Cygni's age) ∈ [0.03, 0.12] | **PASS** — 0.077 |

## 8. Gate

| criterion | status |
|---|---|
| **G8a** pulsar consistency resolved | **PASS** — agreement on existence, age and kinematics; the (a)/(c) degeneracy documented |
| **G8b** γ Cygni reported unsettled | **PASS** — §3, tension T1 |
| **G8c** ²⁶Al not inverted into a count | **PASS** — §5 |
| **G8d** tension list produced | **PASS** — §9 |

## 9. Tension list

| id | severity | statement | resolution |
|---|---|---|---|
| **T1** | reported, not resolvable here | γ Cygni's distance (1.7–2.6 kpc) overlaps our 1.62 kpc only at its extreme low end | left unsettled by design; the ~1.5 kpc literature class carried as a branch |
| **T2** | degeneracy, not disagreement | the pulsar cannot distinguish high-mass explodability from binary stripping | both reported; WP7 already names binary mass transfer as its largest unmodelled effect |
| **T3** | acknowledged systematic | the pulsar's `Ṗ` is contaminated by orbital acceleration | the comparison uses a systematics-widened range, never the point value |

**No tension rises to gate level.** Nothing in the chain was moved toward any
marker.

## 10. What this changes for the project

1. **WP7's worst-case branch is observationally excluded.** The "N_SN = 0"
   outcome required stars above 40 M☉ to leave no neutron stars. One sits inside
   the association, beside a star we count.
2. **The remaining explodability question is narrower.** Not "did any supernova
   happen" — that is settled — but "did the very massive stars explode, or did
   binary stripping let lower-mass ones die first". Both give a non-zero budget.
3. **Four independent external comparisons now agree** with the chain: Härer's
   supernova rate and association mass, Orellana's systemic proper motion, and
   now the pulsar's existence, age and kinematics. None was used as a
   calibration.
4. **The honest headline for a paper**: *Cyg OB2 has produced of order 8
   supernovae, the most recent probably within the last ~100 kyr, and an
   independent compact remnant confirms at least one occurred.*
