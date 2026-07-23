# Glossary and Reference — Cygnus OB2 Supernova History Project

*Companion reference to `paper1_execution_plan.md` and `method_explained.md`. Covers every abbreviation, Gaia catalogue column, variable name, statistical term and physical quantity that appears in this project.*

**How to use this document.** Sections 1–5 are the physical foundations — read these once, in order. Section 6 is a lookup table for the Gaia columns sitting in your data files. Sections 7–11 explain the machinery. Sections 12–14 are pure lookup. Section 16 walks through one real star from your own catalogue end to end; if you only read one thing, read that.

`method_explained.md` already has a short glossary at the end, oriented toward the *wind luminosity / cosmic-ray* method (that's Paper 2). This document covers Paper 1 — the supernova ledger — and the Gaia data mechanics that underpin both.

---

## 0. Project shorthand

| Term | Meaning |
|---|---|
| **WP0 … WP10** | Work Package. The ten stages of `paper1_execution_plan.md`, from requirements extraction (WP0) to manuscript (WP10). Each is meant to be self-contained and separately delegable. |
| **Gate** | A blocking validation criterion at the end of a WP. The WP is not "done" until the gate passes *or* the deviation is documented and explained. Failing a gate silently and continuing is the failure mode the whole plan is built to prevent. |
| **Artifact** | A named output file, `wpN_<name>`. Downstream WPs consume only these, never ad-hoc intermediates. |
| **Branch** | A discrete modelling choice that is *carried through the whole analysis in parallel* rather than averaged — e.g. PARSEC vs MIST isochrones, IMF slope 2.0/2.3/2.6. Results are reported per branch. Averaging over branches hides model systematics; this project forbids it. |
| **Provenance log** | The running record of every cut, threshold, catalogue version and query text, with dates. Becomes the paper's reproducibility appendix. |
| **Author+YY** | Standard astronomy citation shorthand. "Wright+15" = Wright et al., published 2015. "Berlanas+19/20" = two papers by the same lead author in 2019 and 2020. |
| **Ledger** | This project's term for the bookkeeping of stars born vs stars still alive vs stars dead — the difference being the supernova count. |

---

## 1. Units and conversions

| Unit | Meaning | Conversion |
|---|---|---|
| **deg** | degree of angle | 1 deg = 60 arcmin = 3600 arcsec |
| **arcsec (")** | arcsecond | 1/3600 deg |
| **mas** | milliarcsecond | 1/1000 arcsec. The natural unit for Gaia. |
| **mas/yr** | milliarcsec per year | unit of proper motion |
| **pc** | parsec | 3.26 light-years = 3.086 × 10¹⁶ m |
| **kpc** | kiloparsec | 1000 pc |
| **M<sub>☉</sub>** | solar mass | 1.989 × 10³⁰ kg. Written "solar masses" in the plan. |
| **L<sub>☉</sub>** | solar luminosity | 3.828 × 10²⁶ W |
| **Myr / kyr** | mega/kilo-year | 10⁶ / 10³ years |
| **erg** | energy unit (CGS) | 10⁻⁷ J. A supernova releases ~10⁵¹ erg. |
| **K** | kelvin | temperature; the Sun's surface is 5772 K, an O star is 30,000–50,000 K |
| **mag** | magnitude | logarithmic brightness — see §4 |

**The single most useful conversion in this project:**

> **distance (pc) = 1000 / parallax (mas)**

So 1 mas ↔ 1 kpc, 0.7 mas ↔ 1.43 kpc. Cygnus OB2 sits near ϖ ≈ 0.6–0.7 mas.

**Angular size ↔ physical size:** at distance *d*, an angle θ subtends a physical length *s* = θ(radians) × *d*. Practically: **1° at 1.4 kpc ≈ 24 pc**. This is why a 6° × 5.5° selection box is ~150 pc across at Cygnus — far larger than the association itself.

---

## 2. Where things are — coordinate systems

| Term | Symbol | Meaning |
|---|---|---|
| **Right ascension** | α, `ra` | Celestial longitude, 0–360°. The sky's equivalent of longitude, measured eastward from the vernal equinox. |
| **Declination** | δ, `dec` | Celestial latitude, −90° to +90°. Sky equivalent of latitude. |
| **ICRS** | — | International Celestial Reference System. The modern, non-rotating reference frame that (α, δ) are measured in. Gaia's native frame. |
| **Galactic longitude** | *l*, `l_deg` | Longitude measured *in the plane of the Milky Way*, from the Galactic centre. Cygnus OB2 is near *l* = 80°. |
| **Galactic latitude** | *b*, `b_deg` | Angle above/below the Galactic plane. *b* = 0 is the midplane. Cyg OB2 sits near *b* = +1°, i.e. essentially in the plane — which is why the field is so crowded. |
| **Ecliptic latitude** | β | Angle from the plane of Earth's orbit. Matters here only because Gaia's parallax zero point depends on it (Gaia scans the sky in a pattern tied to the ecliptic). |
| **Epoch** | `ref_epoch` | The date the positions refer to. Gaia DR3 uses **J2016.0**. Stars move, so a position is meaningless without an epoch. |
| **Cone / box search** | — | Selecting stars within a circle (cone) or rectangle (box) on the sky. Your ADQL uses an ICRS cone because the Gaia archive's geometry functions work in ICRS, then applies the exact Galactic box locally in Python. |

**Why two coordinate systems?** ICRS is what the telescope measures. Galactic coordinates are what the *physics* cares about — the Milky Way's disc, spiral arms and OB associations are organised in the Galactic frame, not the equatorial one. Your selection box (*l* ∈ [77, 83], *b* ∈ [−1.5, 4]) is defined in Galactic coordinates for that reason.

---

## 3. How far and how fast — astrometry

**Astrometry** = the measurement of stellar positions and their changes. Gaia is an astrometry mission: it measures where stars are, extremely precisely, repeatedly.

### Parallax (ϖ, `parallax`)

As Earth orbits the Sun, nearby stars appear to shift back and forth against the distant background over a year. Half the angular size of that annual wobble is the **parallax**, and it is *inversely* proportional to distance. Gaia measures it in mas.

**The trap:** distance = 1000/ϖ is only valid when the parallax is precise. When the fractional error σ<sub>ϖ</sub>/ϖ is large, inverting the parallax gives a *biased* distance — the errors are symmetric in parallax but wildly asymmetric in distance, and inverting a noisy small parallax can even give a negative or absurd distance. This is why the plan requires distance *posteriors* (Bailer-Jones-style) rather than naive inversion for faint stars.

**Fractional parallax error** f = σ<sub>ϖ</sub>/ϖ is the number to watch. f = 0.01 is superb; f = 0.2 is the conventional edge of "you may invert this"; f > 0.5 means the parallax carries almost no distance information on its own.

### Parallax zero point

Gaia's parallaxes have a small *systematic* offset — a bias, not random noise — of order −0.02 to −0.05 mas, meaning Gaia's parallaxes are on average slightly too small. It depends on the star's magnitude, colour and position on the sky.

- **`parallax_zero_point`** — the per-star correction, computed by the official `gaiadr3-zeropoint` package implementing the **Lindegren et al. (2021)** recipe.
- **`parallax_corrected` = `parallax_raw` − `parallax_zero_point`.**

In your data the median zero point is **−0.032 mas**. At ϖ ≈ 0.7 that's a 4.5% distance shift — about **65 pc at 1.4 kpc**. Not optional: it is comparable to the separation between the two Berlanas distance populations you are trying to test.

### Proper motion (μ, `pmra` / `pmdec`)

The star's real drift across the sky, in mas/yr, after the parallax wobble is removed. Two components:

- **`pmra`** = μ<sub>α*</sub> = μ<sub>α</sub>cos δ — motion in the RA direction. The cos δ factor converts an angle *of longitude* into a true angle on the sky (like how lines of longitude converge at the poles). Gaia's `pmra` already includes this factor.
- **`pmdec`** = μ<sub>δ</sub> — motion in the declination direction.

**Why it matters:** stars born together move together. A cluster shows up as a tight clump in proper-motion space even when it's spread out on the sky. This is the single most powerful membership discriminator. Cyg OB2's mean motion is roughly (−2.7, −4.5) mas/yr with an *internal* dispersion of ~0.5 mas/yr; the surrounding field spans ~10 mas/yr.

### Radial velocity (`radial_velocity`)

Motion *along the line of sight*, in km/s, from the Doppler shift of spectral lines. Together with the two proper-motion components it gives full 3D velocity.

**Availability is the problem.** Gaia measures RV only for relatively bright stars (G<sub>RVS</sub> ≲ 14). In your WP1 catalogue only **18,242 of 245,843 rows (7.4%)** have one, and Cyg OB2 members are heavily dust-obscured, so the member fraction will be far worse. This is exactly why WP6's runaway traceback is 2D and yields a *lower bound*, with 3D deferred to Gaia DR4.

### VPD — Vector Point Diagram

A scatter plot of `pmra` vs `pmdec`. Cluster members form a compact clump; field stars form a broad cloud. One of the standard membership diagnostic plots.

### Correlation columns (`ra_dec_corr`, `parallax_pmra_corr`, …)

Gaia solves for five astrometric parameters (α, δ, ϖ, μ<sub>α*</sub>, μ<sub>δ</sub>) simultaneously, so their errors are **correlated**: if the fit pushed the parallax high it may have pushed the proper motion low to compensate. These ten columns are the off-diagonal terms of the 5×5 covariance matrix, each in [−1, 1].

**Why your project cares:** WP2's Monte Carlo currently uses a *diagonal approximation* — it perturbs each parameter independently, ignoring these correlations. That is flagged in the manifest as a development shortcut and must be replaced with the full covariance before publication, or membership probabilities will be subtly wrong.

---

## 4. How bright and what colour — photometry

**Photometry** = measuring brightness through filters.

### Magnitudes

A logarithmic, **backwards** brightness scale: *smaller magnitude = brighter*. A difference of 5 magnitudes = a factor of exactly 100 in flux; 1 magnitude ≈ factor 2.512.

- **Apparent magnitude** (*m*, e.g. `phot_g_mean_mag`) — how bright the star *looks* from Earth. Depends on intrinsic brightness, distance, and dust.
- **Absolute magnitude** (*M*) — how bright it *would* look at a standard 10 pc. An intrinsic property.
- **Distance modulus** μ = *m* − *M* = 5 log₁₀(*d*/pc) − 5. This is how distance enters a CMD.

With dust, the full relation is:

> ***m* − *M* = 5 log₁₀(*d*) − 5 + *A***

where *A* is extinction (§5). Three unknowns, one equation — which is why extinction, distance and luminosity are entangled and why WP3 exists.

### Gaia's three bands

| Column | Band | Coverage |
|---|---|---|
| `phot_g_mean_mag` | **G** | Very broad white-light band, ~330–1050 nm. The most precise. |
| `phot_bp_mean_mag` | **BP** (blue photometer) | ~330–680 nm |
| `phot_rp_mean_mag` | **RP** (red photometer) | ~630–1050 nm |

### Colour

The **difference** between two magnitudes, e.g. **BP − RP**. Because magnitude is logarithmic, a colour is a flux *ratio* — and a flux ratio between a blue and a red band is a temperature indicator. Hot stars are blue (small BP−RP); cool stars are red (large BP−RP).

**But dust also reddens stars**, so an observed colour mixes temperature and dust. Separating them is the entire job of WP3.

### 2MASS J, H, K<sub>s</sub>

The Two Micron All Sky Survey measured the whole sky in three **near-infrared** bands: J (1.25 μm), H (1.65 μm), K<sub>s</sub> (2.16 μm).

**Why you need them:** dust extinction is far weaker in the infrared than the optical — roughly A<sub>K</sub> ≈ 0.1 A<sub>V</sub>. For a region like Cygnus with A<sub>V</sub> = 4–20 mag, the infrared is where the obscured stars are actually measurable. WP3 fits Gaia *plus* 2MASS photometry jointly; optical-only fits get inflated errors precisely because they cannot break the temperature/extinction degeneracy. **This cross-match is currently missing from your WP1 outputs and blocks WP3.**

---

## 5. Dust — extinction and reddening

Interstellar dust between us and a star does two things: it **dims** the star (extinction) and it **reddens** it (because dust scatters blue light more efficiently than red).

| Term | Meaning |
|---|---|
| **Extinction A<sub>λ</sub>** | Dimming in magnitudes in band λ. A<sub>V</sub> = visual band; A<sub>G</sub> = Gaia G band. A<sub>V</sub> = 5 means the star is ~100× fainter than it would be without dust. |
| **Reddening E(B−V)** | The *colour* change: E(B−V) = A<sub>B</sub> − A<sub>V</sub>. "Colour excess." |
| **R<sub>V</sub>** | The **total-to-selective extinction ratio**, R<sub>V</sub> = A<sub>V</sub>/E(B−V). Describes the *shape* of the dust law — effectively the typical grain size. Diffuse Galactic dust has R<sub>V</sub> ≈ 3.1 (your baseline); denser regions run higher. Your branches are 3.0 / 3.1 / 3.5. |
| **Extinction law** | The curve A<sub>λ</sub>/A<sub>V</sub> vs wavelength. Parameterised by R<sub>V</sub>. |
| **Differential extinction** | Extinction that varies *across* the field, star to star. |
| **3D extinction map** | A model of dust as a function of position *and distance* (Bayestar/Green; Lallement/Vergely). Used in WP3 as an independent consistency check. |

**Why this is the hardest part of Cygnus.** A<sub>V</sub> in this region ranges from ~4 to ~20 magnitudes and varies strongly on small scales. A single average extinction for the field is forbidden — it must be per star. And extinction errors propagate straight into masses (WP4) and completeness (WP5), and therefore into the supernova count. The plan's phrase "differential extinction is THE Cygnus problem" is not rhetoric.

**The degeneracy:** for hot stars, colour becomes insensitive to temperature (all O and early-B stars are about equally blue intrinsically). So a broadband fit cannot distinguish "hot star behind lots of dust" from "cooler star behind less dust." The fix is spectroscopy: a spectral type fixes the *intrinsic* colour, so the observed colour then gives extinction directly. This is why the spectroscopic anchor table matters and why WP3's caveat insists OB-star extinctions come from spectral types, not broadband fits.

---

## 6. The Gaia columns in your files

Every column present in `wp1_gaia_narrow.parquet`. `*_error` columns are the 1σ uncertainty on the corresponding quantity. `*.mask` columns are artefacts of the FITS→Parquet conversion marking missing values.

### Identity and position

| Column | Meaning |
|---|---|
| `source_id` | Unique 64-bit Gaia identifier. The primary key for all cross-matching. Encodes sky position in its high bits, so it is *not* stable across data releases — always match DR3 to DR3. |
| `ra`, `dec` | ICRS position in degrees at epoch J2016.0. |
| `ra_error`, `dec_error` | Positional uncertainty, **mas** (not degrees). |
| `l_deg`, `b_deg` | Galactic coordinates, computed locally by you with Astropy — not native Gaia columns. |
| `ref_epoch` | Epoch of the position, 2016.0 for DR3. |
| `tile_id` | Your own bookkeeping column recording which of the five RA strips a row came from. |
| `local_selection` | Your own flag recording whether the row passed the local Galactic-box cut. |

### Astrometry

| Column | Meaning |
|---|---|
| `parallax` | Raw parallax, mas. **Not yet zero-point corrected.** |
| `parallax_error` | 1σ uncertainty, mas. |
| `pmra`, `pmdec` | Proper motion components, mas/yr. `pmra` includes the cos δ factor. |
| `pmra_error`, `pmdec_error` | Uncertainties, mas/yr. |
| `ra_dec_corr` … `pmra_pmdec_corr` | The ten correlation coefficients of the 5-parameter solution (§3). |
| `radial_velocity`, `radial_velocity_error` | Line-of-sight velocity, km/s. Present for only 7.4% of your rows. |

### Photometry

| Column | Meaning |
|---|---|
| `phot_g_mean_mag` | Mean G magnitude. |
| `phot_bp_mean_mag`, `phot_rp_mean_mag` | Mean BP and RP magnitudes. |
| `phot_g_mean_flux_error` etc. | Flux uncertainties. Note Gaia publishes flux *errors* but magnitude *values*; converting requires σ<sub>mag</sub> ≈ 1.086 × σ<sub>flux</sub>/flux. |
| `phot_bp_rp_excess_factor` | (BP flux + RP flux)/G flux. See §7. |

### Quality and solution diagnostics

| Column | Meaning |
|---|---|
| `ruwe` | Renormalised Unit Weight Error — see §7. |
| `visibility_periods_used` | Number of distinct groups of observations, separated in time. |
| `astrometric_params_solved` | **31** = 5-parameter solution; **95** = 6-parameter solution (see §7). |
| `astrometric_excess_noise` | Extra scatter, in mas, needed to make the astrometric fit acceptable. Non-zero means the single-star model fits poorly. |
| `astrometric_excess_noise_sig` | Significance of that excess. Values > 2 mean the excess is real, not noise. |
| `astrometric_n_good_obs_al` / `_bad_obs_al` | Counts of good/bad along-scan observations. |
| `nu_eff_used_in_astrometry` | "Effective wavenumber" — a colour proxy (μm⁻¹) used in the astrometric solution, for **5-parameter** sources. Input to the zero-point recipe; valid range 1.1–1.9. |
| `pseudocolour` | The equivalent colour proxy for **6-parameter** sources, where colour had to be solved for rather than measured. Valid range 1.24–1.72. |

### Columns you added

| Column | Meaning |
|---|---|
| `parallax_raw` | Copy of the original `parallax`, preserved before correction. |
| `parallax_zero_point` | Per-star Lindegren correction, mas. |
| `parallax_corrected` | `parallax_raw` − `parallax_zero_point`. |
| `parallax_snr` | `parallax_corrected` / `parallax_error`. The inverse of fractional parallax error. |
| `zero_point_boundary_flag` | True if the star fell outside the recipe's validated range (G, ν<sub>eff</sub> or pseudocolour) and had to be clipped — i.e. the correction is an extrapolation. |
| `zero_point_reliable` | The negation of the above. |
| `quality_pass` | Your combined quality flag. |
| `subgroup` | Canonical WP2 kinematic subgroup (`CygOB2-A/B/C`) copied from the authoritative `tables/wp2_subgroup_labels.parquet` sidecar. Rows outside the 1,331 clean automatic-member fit, including the 61 spectroscopic quality exceptions, are explicitly `unassigned`; no second subgroup column name is valid. |

---

## 7. Gaia quality diagnostics, explained

These are the flags that decide which stars you trust. Getting them wrong biases everything downstream.

### RUWE — Renormalised Unit Weight Error

Gaia fits each star with a **single-star model**: one point source moving in a straight line plus a parallax wobble. RUWE measures how badly that model fits, normalised so that a well-behaved single star gives RUWE ≈ 1.0.

**RUWE > 1.4 conventionally flags a bad solution** — most often an unresolved binary, where two stars orbiting each other make the photocentre wobble in a way the single-star model cannot absorb. Also triggered by crowding and by partially resolved close pairs.

**The subtlety for your project:** cutting RUWE > 1.4 preferentially removes **binaries**, and massive stars are overwhelmingly binary (the multiplicity fraction for O stars is ~70%). So a hard RUWE cut is not a neutral quality filter — it is a cut that removes exactly the population you are counting. Your plan already says to exempt the spectroscopic anchors from these filters and treat them manually. Worth a sensitivity test.

### `visibility_periods_used`

The number of separate time-groups in which the star was observed. Astrometry needs observations well spread in time to disentangle parallax (annual, periodic) from proper motion (linear, cumulative). Fewer than ~8 periods and the two become degenerate. Your cut is ≥ 8. The example star in §16 has 28 — plenty.

### `astrometric_params_solved` — 31 vs 95

- **31 (5-parameter):** Gaia had a good colour measurement, so it solved for position (2), parallax (1) and proper motion (2). Uses `nu_eff_used_in_astrometry`.
- **95 (6-parameter):** the colour was unreliable (faint, crowded or very red star), so colour had to be solved for as a sixth free parameter. Uses `pseudocolour`. These solutions are **less reliable** and have a different, larger zero-point correction.

This matters because the Lindegren recipe branches on this value — which your notebook correctly handles.

### `phot_bp_rp_excess_factor`

The ratio (BP flux + RP flux) / G flux. G is measured in a narrow window; BP and RP in wider ones. If a nearby star or nebulosity contaminates the wider BP/RP windows, this ratio inflates. High values flag **unreliable colours**, which is critical in a crowded, nebulous region like Cygnus.

Your current filter checks only that the value is *present*, not that it is *sane*. The standard practice is a colour-dependent cut. Flagged in your own notebook as needing a sensitivity test.

### `astrometric_excess_noise`

Extra angular scatter (mas) needed to reconcile the observations with the model. Like RUWE, elevated values indicate binarity or a problematic fit. Use with `astrometric_excess_noise_sig` — the excess is only meaningful if the significance exceeds ~2.

---

## 8. Stars — types, physics, life and death

### Spectral types

The **OBAFGKM** sequence orders stars from hottest to coolest. Mnemonic: *Oh Be A Fine Girl/Guy, Kiss Me*.

| Type | T<sub>eff</sub> | Mass | Main-sequence lifetime | Note |
|---|---|---|---|---|
| **O** | 30,000–50,000 K | 16–100+ M<sub>☉</sub> | 3–10 Myr | The supernova progenitors. Cyg OB2 has 52. |
| **B** | 10,000–30,000 K | 2.5–16 M<sub>☉</sub> | 10–400 Myr | Those above ~8 M<sub>☉</sub> also explode. |
| **A** | 7,500–10,000 K | 1.6–2.5 M<sub>☉</sub> | ~1 Gyr | |
| **F, G, K, M** | cooler | < 1.6 M<sub>☉</sub> | ≥ 3 Gyr | The Sun is G2. Never explode. |

Each type subdivides 0–9 (O3, O9, B0…) plus a **luminosity class** in Roman numerals: V = main sequence (dwarf), III = giant, I = supergiant.

**"OB stars"** = O and early-B stars collectively — the hot, massive, short-lived population. An **OB association** is a loose, gravitationally unbound grouping of them.

### Core stellar quantities

| Term | Meaning |
|---|---|
| **T<sub>eff</sub>** | Effective temperature — the surface temperature, in K. |
| **log g** | Log of surface gravity (CGS). A proxy for how compact the star is: dwarfs have log g ≈ 4, supergiants ≈ 2. Combined with T<sub>eff</sub> it locates a star in the HR diagram. |
| **Luminosity (L)** | Total energy output, in L<sub>☉</sub>. |
| **Metallicity (Z)** | Abundance of elements heavier than helium. Assumed solar throughout this project. |
| **Main sequence** | The long, stable phase of hydrogen fusion in the core — ~90% of a star's life. |
| **Pre-main-sequence (PMS)** | The contraction phase *before* hydrogen ignition. Low-mass stars in a young cluster are still in it, sitting *above* the main sequence in the CMD. The "turn-on" they form is a second, independent age indicator (WP4 uses both). |
| **Turnoff / turnoff mass m<sub>TO</sub>(t)** | The point where stars are just now leaving the main sequence. The mass whose lifetime equals the cluster's age. **Everything more massive is already dead** — this is the clock at the heart of the ledger. |
| **Wolf-Rayet (WR) star** | An evolved massive star that has shed its hydrogen envelope, exposing the core. Enormous winds. A brief (~0.3 Myr) phase for stars born above ~25 M<sub>☉</sub>. Cyg OB2 has 3 — and their presence alone proves very massive stars formed here. |
| **Coeval** | Formed at the same time. The working approximation for a subgroup; the real spread is 1–3 Myr, which is why "star-formation duration" is a branch. |

### Diagrams

- **HRD (Hertzsprung–Russell diagram)** — luminosity vs temperature. The *theorist's* plot; axes are physical quantities.
- **CMD (colour–magnitude diagram)** — magnitude vs colour. The *observer's* version of the same thing; axes are directly measured. Converting a CMD into an HRD requires distance and extinction — hence WP3 before WP4.
- **De-reddened CMD** — a CMD after extinction correction has been applied. What WP3 delivers.

---

## 9. Populations — the IMF, isochrones and ages

### IMF — Initial Mass Function

The distribution of masses stars are *born* with. Above ~0.5 M<sub>☉</sub> it is a power law:

> **ξ(m) = k · m<sup>−α</sup>**, with α ≈ 2.35 (Salpeter) for the high-mass end.

- **α (slope)** — the *shape*, believed universal. Your branches: 2.0 / 2.3 / 2.6.
- **k (normalisation)** — the *scale*: how many stars this particular association made. Measured per subgroup in WP5.

**The core logic of the whole project:** the IMF shape is universal, so if you count the *surviving* intermediate-mass stars (2–8 M<sub>☉</sub>, which live long enough that none have died yet), you can infer how many high-mass stars were *born*. Compare that to how many are still *alive*. The difference is the number that died — i.e. the number of supernovae. This is the "ledger."

| Term | Meaning |
|---|---|
| **Calibration window** | The mass range used to fit *k* — nominally 2–8 M<sub>☉</sub>. Chosen because it is above the completeness limit but below the mass where stars start dying. |
| **Completeness** | The fraction of stars actually recovered by your pipeline, as a function of mass. Measured by **injection testing**: insert fake stars of known mass into the real data, run the full selection, count how many come back. |
| **Poisson likelihood** | The correct statistics for counting small numbers of objects in bins. Scatter equals counting noise, √N. |
| **Mass function** | The observed distribution of masses (as opposed to the IMF, which is the birth distribution). |

### Isochrones

An **isochrone** (Greek: "same time") is a model curve showing where stars of a **single age** but different masses sit on the CMD/HRD. Fit an isochrone to an observed cluster CMD and the best-fitting age is the cluster's age.

- **PARSEC** — PAdova and TRieste Stellar Evolution Code.
- **MIST** — MESA Isochrones and Stellar Tracks.

**Why both?** They disagree by 20–30% at young ages, from different treatments of convective overshooting and rotation. That disagreement is a real systematic, so the plan carries both as branches and forbids averaging. Since age is the load-bearing parameter of the entire paper, this matters more here than almost anywhere else.

**Rotating vs non-rotating models** — rotation mixes fresh fuel into the core, extending lifetimes and changing inferred ages. Wright+15 quote 2–3 Myr non-rotating, 4–5 Myr rotating for the *same stars*. When comparing your ages to theirs, compare like with like.

### Structure terms

| Term | Meaning |
|---|---|
| **Open cluster (OC)** | A gravitationally bound group of stars formed together. |
| **OB association** | A young, *unbound*, dispersing group of massive stars. Cyg OB2 is one — it is expanding, not held together. |
| **Subgroup** | A distinguishable sub-population within an association, often with its own age and distance. |
| **Substructure** | The general fact that associations are clumpy rather than smooth. Cyg OB2 is strongly substructured. |
| **Superbubble** | A giant cavity in the interstellar medium blown by the combined winds and supernovae of many massive stars. The reason this project exists. |
| **Runaway star** | A massive star ejected from its birthplace at 10–100 km/s, either by a companion's supernova kick or by a dynamical slingshot. Must be traced back and restored to the census — otherwise you count it as "dead" when it is merely absent. |

---

## 10. Supernovae and their tracers

| Term | Meaning |
|---|---|
| **CCSN (core-collapse supernova)** | The terminal explosion of a star born above ~8 M<sub>☉</sub>, when its iron core collapses. Releases ~10⁵¹ erg. |
| **SNR (supernova remnant)** | The expanding shocked shell left behind, visible in radio and X-ray for ~10⁴–10⁵ yr. **γ Cygni (G78.2+2.1)** is the one in your field: age ~7 kyr, distance 1.5–2.6 kpc, association with Cyg OB2 unsettled. |
| **Explodability** | Which masses actually explode versus quietly collapsing to a black hole. Not settled physics. Your branches: "everything above 8 M<sub>☉</sub> explodes" vs **Sukhbold-style "islands of implosion"**, where explodability is a non-monotonic function of mass. |
| **Pulsar** | A rapidly rotating, magnetised neutron star — the compact remnant of a CCSN. Its existence is *proof* a supernova happened. |
| **PSR J2032+4127** | The pulsar in Cyg OB2, in a long-period binary with the Be star MT91 213. Your single strongest independent SN marker. |
| **Characteristic age** | A pulsar age estimate from its spin period *P* and slowdown rate *Ṗ*: τ = P/(2Ṗ). Convenient but can be wrong by factors of several — comparisons must be probabilistic. |
| **P, Ṗ ("P-dot")** | Pulsar spin period and its time derivative. |
| **²⁶Al (aluminium-26)** | A radioactive isotope with a 0.72 Myr half-life, produced by massive-star winds and supernovae. Its 1.809 MeV gamma-ray line traces recent massive-star activity. Measured in Cygnus by **INTEGRAL**. Constrains winds *and* SNe combined, so it is a consistency band, not a fit. |
| **PeV / PeVatron** | 10¹⁵ eV. A "PeVatron" accelerates cosmic rays to PeV energies. The **Cygnus Cocoon** is a candidate — and the Härer+25 interpretation requires a recent supernova to power it. That assumption is what this paper tests. |
| **Cavity** | The low-density bubble carved by winds and prior SNe. A supernova exploding inside one produces a remnant that is faint or invisible in radio — which is how Härer+25 reconcile "recent SN" with "no visible SNR." |
| **N<sub>SN</sub>** | The number of past supernovae — your headline quantity. |
| **R<sub>SN</sub>(t)** | Supernova rate as a function of look-back time: explosions per Myr. |

---

## 11. Statistics and algorithms

### Clustering

| Term | Meaning |
|---|---|
| **DBSCAN** | Density-Based Spatial Clustering of Applications with Noise. Groups points that are densely packed; labels sparse points as noise. Requires no preset number of clusters. |
| **`eps` (ε)** | DBSCAN's neighbourhood radius — how close two points must be to count as neighbours. **In scaled units, not degrees.** The most consequential hyperparameter: too large and everything chains into one blob (this is what happened in your WP2 run at eps = 0.42). |
| **`min_samples`** | Minimum neighbours within ε for a point to be a "core point." Controls how sparse a real cluster may be. |
| **HDBSCAN** | Hierarchical DBSCAN. Builds a hierarchy across *all* density scales instead of fixing one ε, so it handles clusters of differing density — better suited to this region. Used by Paíz+25. |
| **RobustScaler** | Rescales each feature by its interquartile range and centres on the median, so features in different units (degrees, mas, mas/yr) become comparable. Robust to outliers, unlike standard scaling. **Consequence:** an ε in scaled units means different physical distances in each feature — a subtlety worth keeping in mind. |
| **Percolation / chaining** | The failure mode where a density algorithm links neighbour-to-neighbour across the whole field and returns one giant cluster. Diagnostic: the largest cluster contains a large fraction of the input, and its extent matches the selection box rather than the object. |
| **Feature space** | The set of coordinates clustering is performed in — here (*l*, *b*, ϖ, μ<sub>α*</sub>, μ<sub>δ</sub>). |

### Model fitting and comparison

| Term | Meaning |
|---|---|
| **GMM (Gaussian Mixture Model)** | Models a distribution as a sum of Gaussians. Used in WP2 to ask whether parallaxes are better described by one component or two. |
| **BIC (Bayesian Information Criterion)** | A model-comparison score penalising extra parameters: lower is better. ΔBIC < −10 conventionally means "strong preference." |
| **The BIC trap** | BIC's penalty grows as log *N* while the likelihood gain grows as *N*. At large *N* it will prefer more components for *any* distribution that is not exactly Gaussian — including a merely skewed one. With N ≈ 160,000 a large ΔBIC is weak evidence of physical bimodality. |
| **Extreme deconvolution (XD)** | A mixture fit that accounts for each point's *individual* measurement errors, recovering the underlying intrinsic distribution. The right tool for the two-population distance test, because heteroscedastic errors alone can manufacture apparent components. |
| **Heteroscedastic** | Having different measurement errors for different data points. Gaia parallax errors vary by an order of magnitude with magnitude, so this is very much the case here. |
| **Posterior** | In Bayesian inference, the probability distribution of a parameter *after* accounting for the data. This project treats every derived quantity as a distribution, not a number. |
| **Credible interval** | The Bayesian analogue of a confidence interval — e.g. the range containing 68% of the posterior. |
| **Monte Carlo (MC)** | Estimating uncertainty by running many random realisations through the pipeline. `N_MC` = number of draws; `SEED` fixes the random number generator for reproducibility. |
| **MC noise** | The sampling error from using finite draws. With N = 100 draws, a probability near 0.5 carries ~0.05 uncertainty — fine for development, too coarse for publication. |
| **Membership probability** | The probability a star belongs to the cluster rather than the field. Done properly it requires **both** a cluster model and a field model, with the probability being the posterior odds. A test against only a cluster locus, with no field alternative, is circular and returns near-1 for everything it was fitted on. |
| **Control field** | A patch of sky at matched Galactic latitude but *away* from the association, run through the identical pipeline. Whatever "cluster" it returns is your false-positive rate. The essential precision check. |
| **Recall vs precision** | Recall = fraction of true members recovered. Precision = fraction of claimed members that are real. **A gate testing only recall can be passed by a pipeline that labels everything a member.** |
| **Stochastic sampling** | Drawing a discrete, randomly-sampled population rather than using expectation values. Essential in WP7 because with only ~50 massive stars, integer randomness dominates. |
| **Convergence** | Posteriors stop changing when you double the number of iterations. |

---

## 12. Catalogues, surveys, archives and software

| Name | What it is |
|---|---|
| **Gaia** | ESA astrometry mission (2013–), measuring positions, parallaxes, proper motions and photometry for ~1.8 billion stars. |
| **Gaia DR2 / DR3 / DR4** | Data Release 2 (2018), 3 (2022), 4 (**due 2026-12-02**). DR3 is this project's baseline; DR4 adds epoch astrometry and many more RVs, and drives Paper 2. |
| **EDR3** | Early Data Release 3 — the astrometry of DR3, released earlier. The Lindegren zero point is defined on EDR3/DR3 astrometry. |
| **`gaiadr3.gaia_source`** | The main Gaia DR3 catalogue table you query. |
| **2MASS** | Two Micron All Sky Survey. Near-infrared J, H, K<sub>s</sub> photometry of the whole sky. |
| **`tmass_psc_xsc_best_neighbour`** | The Gaia archive's official pre-computed Gaia↔2MASS cross-match table. |
| **VizieR** | CDS's database of published astronomical catalogues. Wright+15 lives there as **J/MNRAS/449/741**. |
| **CDS** | Centre de Données astronomiques de Strasbourg. Hosts VizieR, SIMBAD, Aladin. |
| **SIMBAD** | Object-by-object database — names, types, references for individual stars. |
| **ADQL** | Astronomical Data Query Language. SQL plus geometry functions (`CONTAINS`, `CIRCLE`, `POINT`). |
| **TAP** | Table Access Protocol. The IVOA standard for running ADQL queries against a remote archive. |
| **IVOA** | International Virtual Observatory Alliance — sets these standards. |
| **pyvo** | The Python client for TAP/IVOA services. Your `scripts/gaia_download.py` uses it. |
| **Astropy** | The core Python astronomy library — coordinates, units, FITS I/O. |
| **`gaiadr3-zeropoint`** | The official Python package implementing the Lindegren parallax zero-point recipe. |
| **GOSC** | Galactic O-Star Catalog (Maíz Apellániz). Spectral types for Galactic O stars. |
| **ATNF pulsar catalogue** | The standard pulsar database — periods, derivatives, distances, proper motions. |
| **Green's SNR catalogue / SNRcat** | The two standard supernova-remnant catalogues. |
| **INTEGRAL** | ESA gamma-ray observatory; source of the ²⁶Al measurements. |
| **Bayestar / Green; Lallement / Vergely** | The two main families of 3D dust maps. |
| **FITS** | Flexible Image Transport System — the standard astronomical file format. |
| **Parquet** | Columnar data format; much faster than FITS for tabular analysis in pandas. You store both. |
| **SHA-256** | Cryptographic checksum. Your manifests record one per file so you can prove the data hasn't silently changed. |

---

## 13. Papers referred to by shorthand

| Shorthand | Full reference | Why it's in this project |
|---|---|---|
| **Härer+25** | Härer et al. 2025, A&A, arXiv:2508.21644 | Interprets the Cygnus PeV gamma-ray bubble as requiring a recent SN. **The assumption this paper tests.** |
| **Menchiari+24** | Menchiari et al. 2024, arXiv:2306.00946 | Cygnus cocoon modelling; source of stellar-side assumed values. |
| **Knödlseder+02** | Knödlseder et al. 2002, astro-ph/0206045 | Early Cygnus massive-star census and SN-rate estimate. Prior-art baseline. |
| **Martin+10** | Martin et al. 2010, arXiv:1001.1522 | Cygnus ²⁶Al modelling and SN rates. Prior-art baseline. |
| **Wright+15** | Wright, Drew & Mohr-Smith 2015, MNRAS 449, 741 | The massive-star census of Cyg OB2: 169 OB stars, 52 O-type, 3 WR. **VizieR J/MNRAS/449/741.** Your WP1 gate benchmark. |
| **Wright+16** | Wright et al. 2016, MNRAS 460, 2593 | Cygnus OB2 DANCe high-precision proper motions; found no clear expansion signature. |
| **Berlanas+19** | Berlanas et al. 2019, MNRAS 484, 1838 | Gaia DR2 substructure of Cyg OB2; the **two distance populations** (~1.35 and ~1.6 kpc). Your WP2 gate benchmark. |
| **Berlanas+20** | Berlanas et al. 2020, A&A, arXiv:2008.09917 | O-star spectroscopy; evidence for multiple star-forming bursts. |
| **Fuchs+06** | Fuchs et al. 2006, MNRAS 373, 993 | The IMF-deficit method applied to the Local Bubble. Methodological template. |
| **Zucker+22** | Zucker et al. 2022, arXiv:2201.05124 | Local Bubble formation; Methods section is the second methodological template. |
| **Lindegren+21** | Lindegren et al. 2021, A&A 649, A4 | The Gaia EDR3/DR3 parallax zero-point recipe. |
| **Sukhbold et al.** | — | Source of the "islands of implosion" explodability prescription. |
| **Hunt & Reffert 2023** | — | Large Gaia DR3 open-cluster catalogue (4105 clusters). |
| **Paíz+25** | Paíz et al. 2025, IJAA 15, 171 | Open-cluster census of the Cyg OB2 region with HDBSCAN. Adjacent, not a duplicate — contains no SN analysis. Note the venue is low-quality; use with care. |
| **Celli+24** | — | The cosmic-ray efficiency method dissected in `method_explained.md` Part II. Paper 2 material. |

---

## 14. Variables in your notebooks

| Variable | Meaning |
|---|---|
| `df` | The full WP1 catalogue as loaded, before quality filtering. |
| `analysis` | Rows passing `quality_pass` with all clustering features present. |
| `features` | The list `['l_deg','b_deg','parallax_corrected','pmra','pmdec']` — the clustering feature space. |
| `X` | The RobustScaler-transformed feature matrix fed to DBSCAN. |
| `labels` / `cluster_label` | DBSCAN's output. **−1 means noise** (unassigned); 0, 1, 2 … are cluster indices. |
| `candidate` | Stars in any cluster with ≥ 20 members. |
| `probable` / `members` | Stars surviving the membership probability cut (P > 0.05). |
| `membership_probability` | Fraction of Monte Carlo draws in which the star stayed within the accepted locus. |
| `membership_status` | Categorical binning of the above: rejected / tentative / probable / high. |
| `subgroup` | Canonical GMM subgroup assignment: `CygOB2-A`, `CygOB2-B`, `CygOB2-C`, or `unassigned`. |
| `N_MC`, `SEED` | Monte Carlo draw count (100) and random seed (20260722). |
| `centre`, `scale` | Median and MAD of the candidate set, defining the acceptance locus. **Note:** derived from the same set they are then used to test — the circularity discussed in the WP2 review. |
| `MAD` | Median Absolute Deviation — a robust alternative to standard deviation. |
| `zpt` | The `gaiadr3-zeropoint` module. |
| `quality_pass` | RUWE < 1.4 AND visibility_periods ≥ 8 AND BP/RP excess present AND positive errors AND reliable zero point. |

---

## 15. Cygnus OB2 numbers worth memorising

| Quantity | Value | Source |
|---|---|---|
| Distance | ~1.35–1.75 kpc (two populations debated) | Berlanas+19 |
| Parallax | ϖ ≈ 0.57–0.74 mas | — |
| Age | 1–7 Myr spread; 2–3 Myr non-rotating / 4–5 Myr rotating | Wright+15 |
| Total mass | ~1.6 × 10⁴ M<sub>☉</sub> | literature |
| O stars | 52 | Wright+15 |
| OB stars total | 169 | Wright+15 |
| Wolf-Rayet stars | 3 | Wright+15 |
| Extinction A<sub>V</sub> | 4–20 mag, strongly patchy | — |
| Galactic position | *l* ≈ 80°, *b* ≈ +1° | — |
| Angular scale | 1° ≈ 24 pc at 1.4 kpc | — |
| Bound? | **No** — unbound, expanding, substructured | Wright+16 |

**Sanity anchors for member counts.** Wright+15 catalogue 169 OB stars. Gaia-based member lists for this association run to order 10³–10⁴ stars down to G ≈ 19. A pipeline returning 10⁵ "members" has found the field, not the association.

---

## 16. Worked example — one real star from your catalogue

Taken from `wp2_members.parquet` (i.e. your pipeline currently calls this a member).

```
source_id                  2061169580040317952
ra, dec                    306.2172, +38.9668     deg (ICRS, epoch J2016.0)
l_deg, b_deg                77.3832,  +0.7109     deg (Galactic)
parallax                     1.00323 mas
parallax_error               0.01054 mas
pmra, pmdec                 -1.3677, -3.6856      mas/yr
phot_g_mean_mag             10.3654
phot_bp_mean_mag            10.7713
phot_rp_mean_mag             9.7461
ruwe                         0.8815
visibility_periods_used     28
astrometric_params_solved   31
radial_velocity             (none)
nu_eff_used_in_astrometry    1.5182
```

**Reading it line by line:**

1. **Solution type.** `astrometric_params_solved = 31` → a 5-parameter solution, so Gaia had a good colour and used `nu_eff = 1.518` (comfortably inside the valid 1.1–1.9 range). This is the more reliable of the two solution types.

2. **Astrometric quality.** RUWE = 0.88 — below 1.0, an excellent single-star fit, no sign of an unresolved companion. 28 visibility periods, far above the ≥ 8 requirement. The astrometry is trustworthy.

3. **Parallax and distance.** σ<sub>ϖ</sub>/ϖ = 0.01054/1.00323 = **1.05%** — superb. Applying the median zero point (−0.032 mas): ϖ<sub>corrected</sub> ≈ 1.035 mas, so

   > *d* = 1000/1.035 ≈ **966 pc**

   and the distance modulus is μ = 5 log₁₀(966) − 5 ≈ **9.9 mag**.

4. **Colour.** BP − RP = 10.7713 − 9.7461 = **1.025 mag**. Fairly red — but for a star behind Cygnus dust, that could be an intrinsically hot star that has been reddened, or an intrinsically cool one. Without 2MASS photometry or a spectral type you cannot tell. That ambiguity is precisely what WP3 exists to resolve.

5. **No radial velocity.** At G = 10.4 this star is bright enough that you might expect one, but it falls outside Gaia's RVS coverage. Typical for this field.

6. **The punchline.** Cygnus OB2 sits at 1.35–1.75 kpc, i.e. ϖ ≈ 0.57–0.74 mas. This star has ϖ ≈ 1.035 mas → 966 pc. It is roughly **400–800 pc in front of the association**, and its parallax is precise enough (1%) that this is not a measurement error — it is a genuine foreground star.

   Your current `wp2_members.parquet` nevertheless lists it as a member. That single row is the contamination problem in miniature: it demonstrates why the WP2 gate needs a **precision** criterion (control fields, member-count sanity, spatial extent) and not only the recall criterion currently written into the plan.

---

## Further reading

- **Gaia DR3 documentation** — https://gea.esac.esa.int/archive/documentation/GDR3/ — the authoritative description of every column.
- **Gaia data model, `gaia_source`** — the per-column reference; worth bookmarking alongside §6 above.
- **Lindegren et al. 2021, A&A 649, A4** — the parallax zero point, and why it is not a single number.
- **Bailer-Jones et al. 2021, AJ 161, 147** — how to turn parallaxes into distances properly.
- **`method_explained.md`, Part I** — this project's own derivation of the ledger method, from zero.
