# WP11 Part B — the ²⁶Al / ⁶⁰Fe forward prediction

*Pre-registered in [wp11_isotope_prereg.json](../provenance/wp11_isotope_prereg.json)
before any isotope mass or line flux was computed. Machine-readable record:
[wp11_isotope_forecast_execution.json](../provenance/wp11_isotope_forecast_execution.json).
Brief: [wp11_bowshock_isotope_brief.md §4](../tasks/wp11_bowshock_isotope_brief.md).
Consumes the frozen `repair_v7` chain read-only; changes no published number.*

> ⚠️ **POST-HOC, AND SAYING SO.** Unlike the WP8 markers — frozen at WP1 on
> 2026-07-22, before WP5, WP6 and WP7 existed — this comparison was **not**
> frozen at WP1. It was chosen after the ledger existed and pre-registered
> before scoring. The credibility of the WP8 layer rests on its freeze and the
> two must not be blurred. The manuscript is required to say this in the same
> sentence that introduces the forecast.

**Four of five pre-registered predictions passed on their own terms. I1 passed
but is VACUOUS — it was written against a denominator that turns out to be wrong
by two orders of magnitude, and against the corrected denominator it FAILS. Both
facts are recorded below; neither is suppressed.**

---

## 1. The headline

**⁶⁰Fe splits the headline branch set exactly along α, and COSI can measure it.**

| | branches above COSI's 3σ ⁶⁰Fe sensitivity |
|---|---:|
| α = 2.0 | **18 of 18** |
| α = 2.3 | **0 of 18** |

α is the axis the WP9 verdict hinges on — supported on 18/18 α = 2.0 branches,
0/18 α = 2.3 branches — and it is the one axis Gaia DR4 will *not* settle, since
DR4 sharpens ages and 3D kinematics rather than the high-mass IMF slope. A
2-year COSI survey separates the two arms cleanly on the primary yield arm.
**The MeV line is a genuinely independent probe of the quantity this paper
cannot decide from Gaia alone.**

Two caveats bound that claim immediately and are given equal prominence in §5:
the split holds on one of three yield arms, and the yield arm is worth more than
the branch set (prediction I4).

## 2. What was computed

Per headline branch, importing the WP7 population engine so the forecast cannot
drift from the ledger:

```
M_iso(now) = < Σ over the branch's supernovae of
               y_iso(m_progenitor) × exp(−t_explosion / τ_iso) >
```

evaluated at **each supernova's own progenitor mass** — not at an ensemble
average — then converted to a line flux at the WP3/WP2 distance artifact
(1.6245 kpc), with Cyg OB2 treated as a point source (it subtends ~1°, well
inside the ~3° resolution of both SPI and COSI).

**This is deliberately not the brief's formula.** The brief (§4.2 step 2) writes
*rate × mean lifetime × yield*, which is the **steady-state limit**. Steady state
requires the supernova rate to have been roughly constant for several mean
lifetimes; Cyg OB2's first supernova was ~1.4 Myr ago, against τ = 1.05 Myr for
²⁶Al and **3.78 Myr for ⁶⁰Fe**. The association is nowhere near saturation:

| isotope | M(now) / M(steady state), headline range |
|---|---:|
| ²⁶Al | 0.44 – 0.66 |
| ⁶⁰Fe | **0.16 – 0.31** |

The brief's formula would have **overstated ⁶⁰Fe by a factor 3–6**, in the
direction that manufactures detectability. The decay-weighted sum is the same
physics without the assumption and reduces to the brief's formula in the limit
the brief assumes.

Split-half convergence at 500,000 iterations: max 0.3%, median 0.07%.

## 3. The yield branch — fixed before any flux

Per-supernova ²⁶Al and ⁶⁰Fe yields above 30 M☉ are the dominant uncertainty, and
the brief made fixing them "the whole pre-registration for Part B". Three
published arms were declared and none may be added, dropped or reweighted now.

| arm | role | source | explodability |
|---|---|---|---|
| **LC06_NL** | primary | Limongi & Chieffi 2006, ApJ 647, 483, **Table 3** | every model 11–120 M☉ explodes; Nugis & Lamers (2000) WR mass loss |
| **LC06_Langer** | low | same paper, **Table 5** | same, with the Langer (1989) WNE+WCO mass-loss rate |
| **LC18_REC** | null | Limongi & Chieffi 2018 Recommended, as applied by Falla et al. 2025 | **M > 25 M☉ collapses fully** |

SN-only ²⁶Al is *total − wind*: the wind component was released during the
star's life, not by its explosion. ⁶⁰Fe has no wind component at all — it
reaches the ISM only through the explosion, which is why it is the
supernova-specific tracer and the ledger's cleanest future observable.

### 3a. The null arm is not a low end — it is a different prediction

**The current standard compilation predicts exactly zero.** Every supernova in
this ledger came from a progenitor above 30 M☉ on every branch (WP7: for any
black-hole threshold ≤ 30 M☉ the ledger returns *exactly zero*), and LC18's
Recommended scenario collapses everything above 25 M☉ into the remnant. The
forecast verifies rather than assumes this: **0 of 315,462,579 sampled
supernovae fell at or below 25 M☉**, so the null arm is identically zero, not
merely nearly so.

That is reported, not adopted. WP8 already excluded a zero-supernova Cyg OB2 on
independent grounds — PSR J2032+4127 is a neutron star inside the association,
and neutron stars require successful explosions. The null arm's role is to say
what a COSI **non-detection** would mean: it would favour exactly the
explodability prescription the pulsar disfavours, and the tension between those
two would itself be the result.

## 4. Results

Headline set: α ∈ {2.0, 2.3}, all-explode, 36 branches. Full per-branch table:
[wp11_isotope_forecast.csv](../tables/wp11_isotope_forecast.csv); roll-up:
[wp11_isotope_summary.csv](../tables/wp11_isotope_summary.csv). α = 2.6 is
computed and carried in the table as excluded, not deleted.

**Isotope masses now (M☉), across the 36 headline branches**

| arm | ²⁶Al (SN-only) | ⁶⁰Fe |
|---|---|---|
| LC06_NL | 9.2 × 10⁻⁴ – 3.7 × 10⁻³ | **3.0 × 10⁻³ – 1.4 × 10⁻²** |
| LC06_Langer | 2.6 × 10⁻⁴ – 1.0 × 10⁻³ | 9.6 × 10⁻⁵ – 4.6 × 10⁻⁴ |
| LC18_REC | 0 | 0 |

**Line fluxes (ph cm⁻² s⁻¹), headline medians, and the instrument verdicts**

| arm | F(1809) ²⁶Al | F(1173) ⁶⁰Fe | COSI ⁶⁰Fe | SPI ⁶⁰Fe |
|---|---:|---:|---|---|
| LC06_NL | 8.3 × 10⁻⁶ | 3.6 × 10⁻⁶ | **MARGINAL** | below limit |
| LC06_Langer | 2.3 × 10⁻⁶ | 1.2 × 10⁻⁷ | BELOW_REACH | below limit |
| LC18_REC | 0 | 0 | BELOW_REACH | below limit |

Comparators, all from the literature and cited, none re-derived here:

- **SPI Cygnus ⁶⁰Fe upper limit 1.6 × 10⁻⁵ ph cm⁻² s⁻¹** — combined 1173+1332
  keV, 2σ (Martin et al. 2009, A&A 506, 703). Frozen at WP1.
- **SPI Cygnus ²⁶Al complex flux (3.9 ± 1.1) × 10⁻⁵ ph cm⁻² s⁻¹** — same paper.
  Frozen at WP1. Winds *plus* supernovae, complex-wide.
- **COSI 3σ narrow-line sensitivity 3.0 × 10⁻⁶ ph cm⁻² s⁻¹**, identical at 1173,
  1333 and 1809 keV, in 2 years of survey observations (Tomsick et al. 2023,
  PoS(ICRC2023)745, Table 1). Launch planned 2027.

### 4a. Existing INTEGRAL data already sit at the edge

Prediction I2 passed — but by a **factor of 1.06**. The most supernova-rich
headline branch predicts a combined ⁶⁰Fe flux of 1.51 × 10⁻⁵ ph cm⁻² s⁻¹
against the 1.6 × 10⁻⁵ SPI upper limit. On the primary yield arm, **the existing
INTEGRAL non-detection is already within 5% of excluding the richest α = 2.0
branches.** Martin et al. read their own non-detection as evidence that "no or
very few supernovae went off in the Cygnus complex"; this forecast makes that
reading quantitative and shows it is at the threshold of biting.

### 4b. Predicted ⁶⁰Fe/²⁶Al ratio — context, not a check

The forecast's ⁶⁰Fe(combined)/²⁶Al ratio is **0.76 – 1.00**, against the Galactic
measured 0.184 ± 0.042 (Wang et al. 2020). These are **not the same quantity**:
our ²⁶Al denominator is supernova-only while the Galactic ratio's includes wind
²⁶Al, and every progenitor here sits above 30 M☉ where LC06's ⁶⁰Fe yield is
largest. The comparison is recorded as context and was never used to set or
rescale any yield.

## 5. The predictions, scored

| id | statement | outcome |
|---|---|---|
| I1 | SN-only ²⁶Al mass ≤ 1% of the "~1 M☉" complex inventory | **PASS — but VACUOUS**, see §6 |
| I2 | ⁶⁰Fe below the frozen SPI Cygnus upper limit on every branch | **PASS** (margin 1.06×) |
| I3 | α = 2.0 median ⁶⁰Fe flux ≥ 2× the α = 2.3 median | **PASS** — ratio **2.68** |
| I4 | between-arm yield spread > within-arm branch spread | **PASS** — 30.3× vs 4.8× |
| I5 | COSI ⁶⁰Fe verdict on the primary arm is MARGINAL | **PASS** |

**I3's threshold was deliberately set below the prior estimate.** The brief and
the PROJECT_TRACE WP11 row already recorded a back-of-envelope expectation of
"~3×", so the pre-registration disclosed that partial knowledge and set the bar
at 2× — below the remembered figure — so that I3 could fail against the estimate
that motivated the work. It did not: the measured ratio is 2.68, and the
*threshold* split is sharper still (18/0).

**I4 is the honest limit on all of this.** The yield arm is worth a factor
**30.3** in the predicted ⁶⁰Fe flux; the entire 36-branch headline set is worth
**4.8**. The dominant uncertainty in this forecast is nuclear astrophysics and
Wolf-Rayet mass loss, not the census this paper measures. That spread *is* part
of the result, exactly as the brief's risk register anticipated, and it is why
§1's claim is stated as conditional on the primary arm.

## 6. Finding T4 — a WP8 number does not survive contact with WP1

**WP8 §5 states that the ledger's supernovae sustain ~8.4 × 10⁻⁴ M☉ of ²⁶Al,
"roughly 1000× below the ~1 M☉ inferred for the whole Cygnus complex".** That
denominator is wrong.

The frozen WP1 marker's own measured complex-wide 1809 keV flux, 3.9 × 10⁻⁵ ph
cm⁻² s⁻¹, corresponds at 1.6245 kpc to **8.9 × 10⁻³ M☉** of ²⁶Al — about **112×
less** than 1 M☉. (Over Martin et al.'s quoted 1.0–2.0 kpc range for the
complex: 3.4 × 10⁻³ to 1.3 × 10⁻² M☉.) The 1 M☉ figure appears to be Martin et
al.'s **Galactic** stationary mass, 1.7–2.0 M☉, used for the complex.

Restated in flux space, where no mass conversion of the measurement is required
and the error cannot propagate:

> **This ledger's supernovae alone supply 10–42% of the measured complex-wide
> 1809 keV flux.**

Supernovae are a **sub-dominant but not negligible** contributor — a stronger and
more interesting statement than WP8 made. It is still a *lower bound* on the
complex-wide signal, since Wolf-Rayet winds add ²⁶Al and are expected to
dominate at this age, so it remains **consistent with the measurement and is not
a tension with data**. WP8's qualitative conclusion — the expected ordering, not
a tension — survives; its quantitative margin does not.

**Nothing upstream was retuned.** WP11 is one-way validation: the defect is
recorded, `reports/wp8_crosschecks.md` and
`provenance/wp8_crosschecks_execution.json` are left exactly as they were, and
the manuscript quotes the flux-space comparison. I1 was scored against the
denominator named in its own pre-registration and passed; against the corrected
denominator the same test **fails**, at 10–42% versus a 1% threshold. Recording
the vacuous PASS and the corrected FAIL together is the only honest option: the
pre-registration is not amended after the fact, and the reader is not left with
a number that means nothing.

## 7. The WP8 prohibition is honoured in full

The measured complex-wide ²⁶Al flux was **never inverted into a supernova
count**. The prediction runs forward only: from the ledger's rates and
progenitor masses, through declared literature yields, to a flux.

The one inversion performed — turning the measured flux into an ²⁶Al *mass* in
§6 — is a different operation and is not what the prohibition forbids. It
requires only the decay constant and the distance; inverting to a supernova
count would additionally require a yield and an explodability assumption, which
is exactly the overreach WP8 pre-emptively banned.

## 8. What this cannot do

- It cannot measure a yield. It carries published ones as a branch, and the
  branch is worth 30× — more than everything this project measures.
- It cannot separate Cyg OB2 from the wider Cygnus complex in any existing
  measurement; SPI/COMPTEL resolution is ~3°. COSI does not fix this either.
- It cannot predict the **wind** ²⁶Al, which dominates the complex-wide signal at
  this age and is outside the ledger's scope. This is why the SN-only component
  is stated as a lower bound and never as the prediction for the measurement.
- It cannot decide between the three degenerate WP8 readings of the pulsar. It
  shows only that the null arm and the other two make observationally
  distinguishable predictions — which is what makes a COSI non-detection
  informative rather than merely disappointing.

## 9. Gate

| criterion | requirement | status |
|---|---|---|
| G11a | pre-registration exists with hashed inputs **before** any flux was computed | **MET** — `wp11_isotope_prereg.json`, written and committed before `wp11_isotope_forecast.py` was run |
| G11d | no upstream artifact modified | **MET** — verified against the `audit.py` inventory before and after |
| G11e | yield branch declared before fluxes; ²⁶Al never inverted into a count | **MET** — §3, §7 |
| G11f | manuscript integration passes `wp10_validate.py`; post-hoc disclosure present | **MET** — see [manuscript/README.md](../manuscript/README.md) |
| G11g | failed predictions recorded as failed | **MET** — I1's vacuity and corrected FAIL recorded in §6, not reinterpreted |

G11b and G11c apply to Part A (the bow-shock cross-match), which is **not
adopted for Paper 1** and was not executed. Part A remains a referee-stage
reserve, deferred to the DR4 rerun where radial velocities make the traceback 3D.

## 10. Reproduce

```bash
PYTHONPATH=scripts python3 scripts/wp11_isotope_prereg.py
WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp11_isotope_forecast.py --iterations 500000
```
