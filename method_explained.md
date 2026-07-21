# The Feedback-Ledger Method, Explained From Zero

## A complete step-by-step guide to the Gaia-based measurement of cosmic-ray acceleration efficiency in star clusters, and a detailed reading guide to Celli et al. (2024)

*Prepared for Vadym — July 2026. Companion to the Cygnus OB2 action plan (v3). Every step includes the mechanics (what you do), the rationale (why you do it), and the physics (why nature lets you do it). Analogies are marked with a → symbol.*

---

# PART I — THE METHOD

## 1. The big picture

Gamma-ray observatories (Fermi-LAT, H.E.S.S., HAWC, LHAASO) have detected extended emission around roughly a dozen massive star clusters and the superbubbles they inflate. This emission is believed to come from cosmic rays (CRs) — protons and nuclei accelerated to relativistic energies — colliding with the surrounding gas. The central quantity of the field is the **acceleration efficiency**:

**eta_CR =~ 3 × L_gamma / L_w**

where L_gamma is the gamma-ray luminosity (measured by the gamma-ray instruments) and L_w is the total mechanical power of the stellar winds of the cluster (the energy source available for acceleration). The factor 3 accounts for the fraction of a proton's energy that ends up in gamma-rays through pion production, under the assumption that protons lose essentially all their energy in the region (the "calorimetric" assumption).

The numerator, L_gamma, is measured. The denominator, L_w, is currently *assumed* — taken from heterogeneous literature with factor 2–10 uncertainties. Consequently, eta_CR is known only to order of magnitude (~1%), and the question "is it universal across systems, and if not, what controls it?" — which discriminates between acceleration mechanisms — cannot be answered.

The method described here measures the denominator, star by star, from Gaia data, together with a second quantity nobody currently provides: the **supernova history** of each system, which decides whether the gamma-ray emission can be attributed to winds at all.

**The whole chain in one breath:** stellar lifetimes fall steeply with mass (physics) → below 8 solar masses nobody in a young cluster has died yet and Gaia sees everyone (window choice) → counting stars there fixes the one free parameter of a universal birth distribution (the IMF normalization k) → k plus the cluster age splits the born massive stars into dead ones (giving the supernova count and timeline) and living ones (whose measured winds sum to L_w) → a Monte Carlo turns every arrow into a probability distribution → eta_CR with honest error bars.

Each arrow is now unpacked completely.

---

## 2. Step 0 — Membership: deciding who belongs to the cluster

### What you do

You cannot simply query "all stars within this patch of sky." A sky-region query returns mostly *strangers*: foreground stars between us and the cluster, and background stars far behind it. In the Galactic plane (where Cygnus OB2 sits), strangers outnumber true members enormously — the line of sight passes through kiloparsecs of star-filled disk. The membership step filters the query result down to the stars that physically belong to the association. Three filters, applied in sequence:

1. **Parallax (distance) window.** Parallax is the tiny apparent annual wobble of a star's position caused by Gaia orbiting the Sun; its size is inversely proportional to distance (distance in parsecs = 1/parallax in arcseconds). Cyg OB2 sits at roughly 1.4–1.7 kpc, i.e. parallax around 0.6–0.7 milliarcseconds. You keep stars with parallax in a *generous* window (say 0.45–0.9 mas), wide on purpose: parallax errors at this distance are significant, and a tight cut would silently discard genuine members with noisy measurements.

2. **Proper-motion clustering.** Proper motion is a star's angular drift across the sky in milliarcseconds per year. Stars born together move together — they inherited the velocity of one parent cloud. Members of Cyg OB2 therefore clump tightly in proper-motion space, while the strangers (at all distances, on all orbits) smear broadly. You run a clustering algorithm (HDBSCAN — a density-based method that finds clumps of points without being told how many clumps to expect or what shape they have) on the combined space of sky position, parallax, and the two proper-motion components. The dense clump *is* the association.

3. **Quality cuts.** RUWE (renormalized unit weight error) is Gaia's goodness-of-astrometric-fit statistic; RUWE > 1.4 usually flags unresolved binary stars or contaminated measurements, whose parallaxes and proper motions cannot be trusted. You also require the parallax to be measured at better than ~3 sigma.

→ *Analogy: a photo of a crowded square during a marathon. You want the runners, not the bystanders. Distance cut = focus depth (ignore people much closer or farther than the course). Proper-motion clustering = take two photos a minute apart and keep only the people who all moved the same direction at the same speed. The bystanders wander randomly; the runners move coherently.*

### Two technical details that matter

**The parallax zero-point.** Gaia parallaxes carry a small systematic offset (on average about −17 microarcseconds, varying with star color and brightness). At 1.5 kpc a star's parallax is ~670 microarcseconds, so ignoring the offset biases all distances by ~2.5% coherently — a distortion of the whole 3D structure. An official correction recipe (Lindegren et al.) exists; apply it before anything else.

**Probabilities, not verdicts.** Every star's parallax and proper motion have error bars. Instead of a binary member/non-member decision, you draw many random realizations of each star's astrometry from its error distribution, re-run the clustering, and record the *fraction* of realizations in which the star lands in the cluster. A star that makes it 60% of the time enters all later sums with weight 0.6. This way measurement uncertainty flows honestly into every downstream number instead of being hidden by a sharp cut.

### Why this step is half the project

Every subsequent quantity — ages, masses, supernova counts, wind luminosity — is computed *from the member list*. Contamination (strangers counted as members) inflates star counts and corrupts the mass function; incompleteness (members missed) deflates them. The published member lists for Cyg OB2 (Wright et al. 2015 census; Berlanas et al. 2019 spectroscopic work) exist precisely so you can validate: if your selection recovers >=80% of the spectroscopically confirmed members, your machinery works; if not, stop and debug. This deliberate reproduction of a known result before producing new ones is the "validation gate" — the pipeline equivalent of calibrating an instrument on a standard source.

### The Cygnus-specific complication

Cygnus OB2 is not a neat bound cluster; it is an *association* — spatially extended, kinematically loose, with internal subgroups of different ages, and (per Berlanas et al.) at least two stellar populations at different distances (~1.35 and ~1.6 kpc) overlapping on the sky. The clustering will fragment it. That fragmentation is *physical* — keep the subgroups separate, because they have different ages, and age is the most sensitive input of everything that follows.

---

## 3. Step 1 — Extinction and the color–magnitude diagram

### The physics of extinction

Interstellar space contains dust. Dust grains absorb and scatter starlight, making stars look **dimmer** (extinction, denoted A_V in the visual band or A_G in the Gaia band) and **redder** (because blue light scatters more efficiently than red — the same reason sunsets are red). A hot blue massive star behind a dust cloud can masquerade as a cooler, fainter star. If you do not correct for this, every mass and age you derive is wrong.

Cygnus is one of the worst directions in the sky for this: the line of sight runs *along* a spiral arm, and extinction reaches A_V ~ 5–10 magnitudes, varying strongly from star to star across arcminutes. Consequently:

- **A single average extinction for the whole region is forbidden.** You must estimate extinction *per star*, by fitting each star's brightness in several bands simultaneously (Gaia's G, BP, RP plus the infrared J, H, Ks bands from the 2MASS survey — infrared is essential because dust dims it far less; the Gaia archive provides a pre-computed cross-match to 2MASS).
- **The extinction law itself is a systematic.** The ratio of total to selective extinction (R_V, normally 3.1) describes the grain-size mixture; it is anomalous in parts of Cygnus. Repeating the analysis with R_V between 3.0 and 3.5 and carrying the difference as an error term is part of honest bookkeeping.

### The CMD

With extinction corrected, you plot the **color–magnitude diagram (CMD)**: each member star is one point, with color (a proxy for surface temperature: blue = hot, red = cool) on the horizontal axis and intrinsic brightness on the vertical. This is the observational version of the Hertzsprung–Russell diagram, and it is the canvas on which the next step operates. A healthy young cluster's CMD shows a well-defined diagonal band (the main sequence — stars fusing hydrogen in their cores) with a bend at the top.

---

## 4. Step 2 — Isochrones: reading the cluster's age off the CMD

### The physics: massive stars live fast

A star's engine burns hydrogen at a rate that rises very steeply with mass: luminosity scales roughly as mass to the power 3.5. Fuel scales only linearly with mass. Lifetime =~ fuel/burn-rate therefore *falls* steeply with mass, roughly as mass to the power −2.5. Concrete numbers worth memorizing:

| Birth mass (solar masses) | Main-sequence lifetime |
|---|---|
| 1 (Sun) | ~10,000 Myr |
| 2 | ~1,000 Myr |
| 8 | ~40 Myr |
| 20 | ~9 Myr |
| 40 | ~5 Myr |
| 60 | ~4 Myr |
| 100 | ~3 Myr |

→ *Analogy: candles of different thickness lit at the same moment — except backwards. In a cluster, the THICKEST candles burn out first. A cluster is a shelf of candles all lit simultaneously; glancing at which candles are still burning tells you how long ago the match was struck.*

### What an isochrone is

An **isochrone** (Greek: "equal time") is a theoretical curve computed from stellar-evolution models: it marks where stars of *one common age* but *all different masses* fall on the CMD. A 3 Myr isochrone looks different from a 10 Myr isochrone, chiefly at the top: as the cluster ages, its most massive stars leave the main sequence (swelling, cooling, then dying), so the top of the observed band peels away from the model line at a location called the **main-sequence turnoff**. The turnoff position is effectively a clock hand.

### What you do

Overlay a grid of model isochrones (ages 1–10 Myr, in fine steps) on the de-reddened CMD of each subgroup, and find the age whose isochrone best matches, via a proper likelihood over all stars (not by eye). Two independent families of stellar models exist — PARSEC and MIST — which disagree at the ~20–30% level for very young ages. Fit with both, carry both results forever, never average them silently: the spread *is* the model systematic.

Two things fall out of the fit:

1. **The subgroup age t** — the master clock for the supernova ledger.
2. **A mass for every member star** — each star's position along the isochrone corresponds to a unique mass, so you read masses off by interpolation.

### Pitfalls (each becomes an error term)

- **Unresolved binaries** sit ~0.75 magnitudes above the single-star sequence (two stars' light in one point), mimicking younger/heavier stars.
- **Differential extinction residuals** scatter stars off the sequence, blurring the turnoff.
- **The brightest stars are the worst measured** — very bright stars saturate Gaia's detectors, and broadband colors barely distinguish a 35,000 K star from a 45,000 K one (both look "very blue"; the difference is in the ultraviolet, which Gaia doesn't see). For the ~50 brightest members, override the photometric estimates with spectroscopic temperatures and luminosities from the literature (spectral types from the Wright/Berlanas censuses convert to calibrated temperatures). This photometric+spectroscopic hybrid is standard practice; documenting which star got which treatment is what makes it defensible.

---

## 5. Step 3 — The initial mass function: the universal birth ratio

### What it is

When a molecular cloud collapses and fragments into stars, the birth masses follow a statistical distribution called the **initial mass function (IMF)**: xi(m) = dN/dm, the number of stars born per unit interval of mass. Empirically, above ~0.5 solar masses it is a power law:

**xi(m) = k × m^(−2.3)**

(the Kroupa/Salpeter slope). Nature makes many featherweights and few heavyweights, in a specific steep ratio: for every star born at 20 solar masses, roughly a couple of hundred are born near 1 solar mass.

### Why we are allowed to use it

The decisive *empirical* fact is that this shape appears to be **universal**: measured across clusters, environments, and metallicities, the same slope keeps emerging (its physical origin lies in the turbulent fragmentation of clouds). Universality is an assumption and must be stated as one — it is debated at the extremes (very massive clusters, starburst environments) — but it is the best-tested assumption in this entire chain, far more solid than any wind prescription.

Because the *shape* is fixed, an individual cluster's entire birth population is described by **one number**: the normalization k — "how big was this star-formation event." The whole method reduces to measuring k for your cluster, then reading predictions off the fixed shape.

→ *Analogy: a bakery with a fixed recipe book — for every wedding cake it bakes, it always makes 30 loaves and 200 rolls. If you walk in at closing time and count 100 rolls left (rolls keep — nobody "eats" them in a day), you know the day's production run: half a standard batch. And therefore you know how many wedding cakes were made this morning — even though every cake is already gone. The rolls are witnesses to the cakes.*

---

## 6. Step 4 — Why the calibration window is 2–8 solar masses

To measure k you must count stars in a mass range where **the stars you count today equal the stars that were born**. Two conditions must hold simultaneously:

**Condition 1: Nothing has died there ("immortal").** From the lifetime table: an 8 solar-mass star lives ~40 Myr. Your cluster is a few Myr old. Therefore every star ever born below 8 solar masses is still alive — the count is complete *in time*.

**Condition 2: Gaia sees all of them ("visible").** At Cygnus distance and extinction, stars below ~1–2 solar masses start falling below detection limits or acquiring garbage astrometry — the count there is incomplete *in brightness*. Above ~2 solar masses, Gaia plus 2MASS catch essentially everyone (verifiable, see completeness below).

The window [2, 8] solar masses is where both conditions hold: **immortal and visible**. Counting there is counting births directly. Below 2: alive but partly invisible. Above the turnoff: visible in principle but partly dead. The window is where the bookkeeping is exact.

---

## 7. Step 5 — Fitting the normalization k

### Mechanics, in full

1. Take all members with masses in [2, 8] solar masses (masses from Step 2).
2. Bin them in mass — for example six logarithmic bins. Each star contributes its *membership probability*, not 1 (the 60% member counts 0.6).
3. The model prediction for each bin [m1, m2] is k × integral of m^(−2.3) from m1 to m2.
4. Fit k by maximizing a **Poisson likelihood** — Poisson, not least-squares, because bins contain small counts, and the natural scatter of a count of 12 is sqrt-of-12-like counting noise, not Gaussian noise with some independent error bar. Least-squares on small counts biases the fit.
5. Optionally let the slope vary as well (2.3 ± 0.3) — better: keep it fixed at 2.3 for the headline result and repeat with 2.0 and 2.6 as bracketing systematics.

### The completeness correction

Before trusting the counts, ask: of the stars of mass m that are truly there, what fraction did my selection actually catch? Answer by **injection testing**: sprinkle synthetic stars of known masses into the real catalog (with realistic extinction for their sky position), run your entire selection pipeline on them, and record the recovery fraction as a function of mass. Divide the model prediction by that fraction. If recovery is ~100% across [2, 8] — the design goal of choosing that window — the correction is cosmetic; but you must *demonstrate* that, not assume it.

### Toy example (numbers for intuition)

Suppose you count 300 members in [2, 8] solar masses. Fitting k so that the integral over [2, 8] equals 300 fixes k. That same k then predicts (for slope 2.3, ceiling ~120 solar masses; exact values shift with details — the mechanics is the point):

- born above 20 solar masses: ~20–25 stars
- born above 40 solar masses: ~7–9 stars
- born above 60 solar masses: ~4 stars

You never observed most of these heavyweights — many are dead — yet their birth numbers are pinned by the 300 lightweights, because both groups were drawn from the same universal distribution at the same event. **The lightweight survivors are demographic witnesses to their dead siblings.**

---

## 8. Step 6 — The dead branch: counting supernovae

### The turnoff mass as executioner's ledger

Define **m_TO(t)**: the mass whose lifetime exactly equals the cluster age t (from the same stellar models as the isochrones — consistency matters). Every star born heavier than m_TO is dead; every star born lighter still shines. For t = 5 Myr, m_TO =~ 40 solar masses.

Stars born above ~8 solar masses end their lives as **core-collapse supernovae** — the core runs out of fusable fuel, collapses to a neutron star or black hole, and the envelope explodes with ~10^51 erg of kinetic energy. Stars below 8 die quietly as white dwarfs. In a young cluster m_TO is far above 8, so the *entire* dead branch consists of supernovae:

**N_SN = k × integral of m^(−2.3) from m_TO(t) to m_max**

With the toy numbers and t = 5 Myr: the ~7–9 stars born above 40 solar masses are your supernova count.

### Corrections that must be applied

1. **Escaped runaways.** 10–20% of O stars get ejected from their birthplace at high velocity — either flung out by the supernova explosion of a binary companion, or slingshotted by close gravitational encounters in the dense cluster core. A runaway is *missing from the cluster but not dead*; naively it would be counted as a supernova. Fix: search a wide field around the association for massive stars whose proper motions, traced backward in time, converge on the cluster within its lifetime — then add them back to the living census. With Gaia DR3 the traceback uses mostly 2D sky motion (a lower bound on recoveries); DR4's radial velocities and epoch astrometry upgrade it to full 3D.
2. **Not every massive star explodes.** Some collapse directly into black holes with no visible supernova (the "failed supernova" channel; plausibly 15–25% of core collapses, mass-dependent, uncertain). Quote both numbers: "stars that died" (robust) and "likely luminous supernovae" (with an explodability prior).
3. **Binary mass transfer.** Interacting binaries exchange mass, which shifts effective lifetimes and can even produce supernovae from stars born slightly below 8 solar masses (and prevent some above). On a ledger of single digits, treat as a ±1 systematic.

---

## 9. Step 7 — The supernova *history*: why mass is a clock

### The common confusion, addressed head-on

"The IMF deficit tells how many exploded — but how can we possibly tell *when*?" The missing link: **the same lifetime law that tells you who died tells you when each of them died, because death date is a function of mass.**

All the cluster's stars were born together, t Myr ago (coeval formation — an assumption, see caveats). Lifetime tau(m) is a known, steeply decreasing function. So a star of birth mass m did not die at some unknown moment — it died at exactly one moment: tau(m) after birth, which is (t − tau(m)) before today. Heavier means earlier. The dead stars exploded in strict mass order, heaviest first.

→ *Analogy: fuses of different lengths, all lit at the same instant. Even if you arrived after some have burned out, the set of stubs plus the known burn rates reconstructs the full firing schedule: which fuse popped when. Mass is fuse length; the lighting moment is the cluster's birth; the stubs are the missing stars.*

### Worked example

Subgroup age t = 6 Myr, so m_TO =~ 32 solar masses — everything born heavier is dead. Say the normalized IMF implies three stars born above the turnoff: roughly one near 90, one near 50, one near 35 solar masses. Then:

- the 90 solar-mass star (lifetime ~3 Myr) exploded **3 Myr ago**;
- the 50 solar-mass star (lifetime ~4.3 Myr) exploded **1.7 Myr ago**;
- the 35 solar-mass star (lifetime ~5.7 Myr) exploded **0.3 Myr ago**.

That ordered list *is* the supernova history: first explosion 3 Myr ago, roughly one every ~1.3 Myr since, most recent 0.3 Myr ago. Formally, the explosion rate as a function of look-back time is R(tau) = xi(m(tau)) × |dm/dtau|: the IMF says how many stars occupy each mass sliver; the lifetime function converts mass slivers into time slivers.

### The two derived numbers that matter most

1. **Time of the first supernova** = t − tau(m_max): when the bubble began to be supernova-driven. Directly comparable to the bubble's dynamical age estimated from its size and expansion speed — an independent cross-check between stellar and ISM clocks.
2. **Time since the last supernova**, set by the mass just above today's turnoff. This is the *classifier* for the gamma-ray interpretation: if the last explosion was ~50 kyr ago, the observed gamma emission may well be powered by a young supernova remnant rather than by winds, and eta_CR computed against L_w alone is meaningless; if no supernova has ever occurred (very young system), winds are the only possible engine and eta_CR is clean. (This is precisely the live issue in Cygnus: the leading model postulates an energetic supernova ~50 kyr ago; the ledger assigns that postulate a probability.)

### Caveats — each one a Monte Carlo branch

1. **You never know the actual masses of the dead stars** — only the distribution they were drawn from. So you obtain not *the* history but an ensemble of possible histories from repeated random IMF draws consistent with the surviving 2–8 count. "Last SN 0.3 Myr ago" becomes "probability 60% that the last SN was within the past 0.5 Myr." For classification purposes, the probabilistic statement is exactly what is needed.
2. **Coeval formation is approximate** — real clusters form stars over ~1–3 Myr, smearing every date by that much. Working at subgroup level (each subgroup much closer to single-age) absorbs most of this.
3. **Binaries reshuffle lifetimes** by ~1 Myr — another smearing term.
4. **The whole timeline rides on the age**: if t is wrong by 1 Myr, every date shifts by 1 Myr *and* the turnoff mass moves, changing N_SN itself. Between ages 2 and 4 Myr, m_TO slides from ~120 to ~55 solar masses — i.e., between "no supernovae yet" and "several." Age is the load-bearing input of the ledger; always present N_SN and the timeline *as functions of age*, not only at the best-fit age.

**Summary sentence: "how many were born vs how many remain" is the count; the history is the same calculation read along the time axis instead of the mass axis — no extra data needed, because for coeval stars, mass is a clock.**

---

## 10. Step 8 — The living branch: wind luminosity, star by star

### Critical design point first

For the **living** massive stars you do *not* use the IMF — you use the actual stars. They are the brightest objects in the field; Gaia, 2MASS, and decades of targeted spectroscopy see them individually. The IMF's role for the living branch is only to *audit* the census (Step 9). This is the fundamental difference from the statistical approach (Part II): where the power resides, you measure; the IMF only fills bookkeeping gaps.

### The physics of line-driven winds

A massive star's atmosphere is so luminous that light itself blows off the surface layers. Photons streaming outward scatter in the millions of absorption lines of metal ions (iron, carbon, nitrogen...) in the atmosphere; each scattering transfers a photon's momentum to the gas. The cumulative push launches a continuous outflow — the **stellar wind** — described by two parameters:

- **Mass-loss rate (M-dot):** how much stellar material leaves per year. Because the driving is radiative, M-dot climbs steeply with luminosity — roughly as L to the power 1.8–2.2. An O-type star loses ~10^−7 to 10^−6 solar masses per year; the Sun, for comparison, ~10^−14.
- **Terminal velocity (v_inf):** the wind's final coasting speed, scaling with the escape velocity from the stellar surface (=~2.6 × v_esc for hot stars) — typically 2000–3000 km/s for O stars.

The wind's mechanical power — kinetic energy exported per second — is **(1/2) × M-dot × v_inf²** per star, and the cluster total is the sum:

**L_w = sum over living massive members of (1/2) M-dot_i v_inf,i²**

### Why a handful of stars dominates (feel for the numbers)

- O5V star (M-dot ~ 10^−6, v_inf ~ 2800 km/s): ~2×10^36 erg/s
- B2V star: several orders of magnitude less
- One **Wolf-Rayet star** (an evolved massive star that has shed its hydrogen envelope; M-dot ~ 10^−5, v_inf up to 5000 km/s): 10^37–10^38 erg/s — *comparable to all O stars of the cluster combined*.

Because of this steepness, L_w is effectively the property of the ~10 most extreme individual objects. This is simultaneously (a) why the star-by-star treatment is mandatory for per-system physics, (b) why Wolf-Rayet stars get special treatment — their parameters come from dedicated WR spectroscopy (they cannot be characterized by broadband photometry at all), and (c) why the statistical/ensemble approach of Part II fails for individual systems: the top of the IMF is where small-number randomness is largest, and it is exactly where the power lives.

### Prescriptions: converting a star into its wind

You cannot measure M-dot for every star directly; you compute it from the star's basic parameters (luminosity, mass, temperature, metallicity) using a **prescription** — a formula fitted to wind theory and calibrated on observed stars. Two generations coexist:

- **Vink et al. (2001):** the community standard for two decades, based on Monte Carlo radiative-transfer models.
- **Björklund et al. (2021):** modern theoretical models, systematically ~2–3× *lower*.

The gap between them reflects the unresolved **clumping problem**: if winds are clumpy rather than smooth (they are), then classical observational M-dot estimates are inflated, and the true rates are lower. Nobody knows the final answer, so the honest procedure is to run *both* prescriptions as parallel branches and report the pair. The spread is not your failure — it is the current state of wind physics, quantified, and propagated for the first time into cosmic-ray efficiency estimates.

### Assembly

For each living member above ~8 solar masses: take temperature and luminosity (photometric from Step 2 for most; spectroscopic literature values overriding for the top ~50 stars), apply both prescriptions to get M-dot, get v_inf from the escape-velocity scaling (or spectroscopic values where published), compute the star's (1/2) M-dot v_inf², multiply by membership probability, sum. Add the Wolf-Rayet stars from dedicated literature. Apply the completeness correction from the audit (Step 9) for extinction-hidden members. All inside the Monte Carlo.

---

## 11. Step 9 — The closure test: making the method falsifiable

The fitted k makes one more prediction not yet used: how many massive stars should be **alive right now** — k × integral from 8 to m_TO(t). Compare with the directly observed count of living massive members. Three possible outcomes:

1. **Agreement** → membership is complete, the IMF assumption holds in this cluster, and the L_w sum rests on a verified census.
2. **Observed < predicted** → either extinction hides members (quantify → this *is* the completeness correction applied to L_w), or members escaped as runaways (→ the traceback correction), or both.
3. **Persistent disagreement after corrections** → the IMF assumption itself is failing in this system; flag honestly rather than paper over.

This closure is what elevates the construction from curve-fitting to physics: the intermediate-mass stars *calibrate*, the massive stars are *predicted*, and the prediction is *checked* against the very stars whose winds you then sum. Each branch of the calculation constrains the others.

---

## 12. Step 10 — The Monte Carlo: why, exactly

"Monte Carlo" = estimating uncertainty by simulating thousands of random realizations instead of propagating error formulas. Three distinct reasons force it here — no analytic formula survives this chain:

1. **Measurement errors.** Every star's parallax, proper motion, extinction, temperature, and mass has an error bar, and they feed a nonlinear pipeline (clustering → CMD → isochrone → integrals). Redraw all measurements from their error distributions thousands of times, re-run everything, and the spread of outcomes *is* the statistical error on k, N_SN, L_w.
2. **IMF stochasticity — the deep one.** The IMF is a probability distribution; a single cluster is one random *draw* from it. "The IMF predicts 5 stars above 35 solar masses" means 5 *on average over many clusters*; this particular cluster may have drawn 2 or 9 — Poisson-like randomness on single digits. Simulate many random populations consistent with the observed 2–8 census; the spread of dead-star counts across simulations is the honest, typically **asymmetric** error on N_SN ("5, with 68% interval 3–8" — never a clean ±1). The same simulations automatically produce the ensemble of supernova timelines of Step 7.
3. **Discrete model choices.** Vink vs Björklund; PARSEC vs MIST; IMF slope 2.0/2.3/2.6; R_V 3.0/3.5. These are not statistical errors — they are forks in the road. Carry them as labeled parallel branches; report headline numbers per branch. Averaging over model choices hides exactly the information a referee (or collaborator) needs.

→ *Analogy: a court reconstructing an event from witnesses. Instead of one "best guess" narrative, the court runs every scenario consistent with all testimony and reports which conclusions hold in 95% of scenarios. The 2–8 solar-mass stars are the witnesses; the scenarios are the simulated clusters; the verdict on N_SN comes with a confidence, not a pretense of certainty.*

The deliverable of the entire pipeline is not "eta_CR = 0.9%" but "eta_CR = 0.9%, 68% interval [0.5, 1.6], prescription systematic ×2.3, calorimetric assumption stated" — and the interval is the publishable part, because the whole scientific point is replacing unquantified denominators with quantified ones.

---

## 13. Step 11 — Closing the loop: eta_CR

1. Take the gamma-ray luminosity from the literature: for Cygnus, the Fermi-LAT cocoon spectrum and the LHAASO bubble. **Rescale to your Gaia distance** — published luminosities assumed some distance d0, and L scales as distance squared; Gaia has shifted cluster distances by tens of percent (Westerlund 1 moved ~30%), so even this trivial step changes eta_CR materially.
2. Compute **eta_CR = 3 L_gamma / L_w** with full error propagation (L_gamma errors from the gamma catalogs; L_w posterior from your Monte Carlo; distance uncertainty correlating both).
3. State the calorimetric assumption explicitly: the factor 3 presumes accelerated protons lose essentially all their energy to pion production within the region. Where confinement is weak and protons escape, the true eta_CR is *higher* than computed — the result is then a lower limit. (Whether calorimetry holds per system is itself a theory question — precisely the discussion to hand to the theory-side collaborators.)
4. Classify the system with the supernova ledger *before* interpreting eta_CR: wind-only (no SN yet — clean), SN-contaminated (recent SN plausible — eta_CR against winds alone is ill-defined), or post-SN (SNe long ago, remnants faded — winds again dominant but history matters for the bubble's CR reservoir).

**Data availability, stated honestly:** the GeV side (Fermi-LAT) is fully public — photon-level data, catalogs, diffuse models — so both luminosities and *upper limits for non-detected systems* can be computed independently. The TeV side (H.E.S.S., HAWC, LHAASO) publishes catalogs and spectra but not event data — TeV enters as published values only, and no custom TeV upper limits are possible. The homogeneous, self-computed layer of any multi-system study therefore lives at GeV.

---

## 14. The corrected one-paragraph pipeline

Query Gaia around Cygnus OB2 (the association — the engine; not the cocoon — the emission region) → filter to actual members by parallax window plus common proper motion, with per-star membership probabilities → extinction-correct each star individually and build the CMD per subgroup → isochrone fits give subgroup ages and per-star masses → the IMF, normalized on the 2–8 solar-mass survivors (immortal and visible), gives the born massive-star count, which splits by the age into the dead branch (N_SN and the explosion timeline, since mass is a clock) and the predicted living branch (audited against the directly observed massive stars, whose measured winds — both prescriptions — sum to L_w) → Monte Carlo over measurements, IMF draws, and model branches turns every number into a distribution → classify the system by time-since-last-SN, then eta_CR = 3 L_gamma / L_w with the literature gamma flux rescaled to the Gaia distance.

---

## 15. The method's lineage (know it, cite it)

The IMF-deficit supernova counting is an *established* technique with a 25-year pedigree — a feature, not a weakness: referees cannot attack the method's validity, only the inputs, which is where DR4 makes you strongest.

- **Maíz-Apellániz (2001)** — the founding application: evolutionary synthesis of the Sco-Cen subgroups → ~20 supernovae in 10–12 Myr, plus Hipparcos traceback showing the groups crossed the Local Bubble's position.
- **Fuchs et al. (2006)** — the canonical IMF-deficit form: fit the IMF to surviving members, count the missing high-mass tail → 14–20 supernovae for the Local Bubble.
- **Breitschwerdt et al. (2016)** and the iron-60 literature — consumed those counts for radioisotope timing of nearby supernovae.
- **Voss, Diehl, Knödlseder, Martin (~2009–2012)** — general population-synthesis framework predicting supernova histories, aluminum-26, iron-60, and energy budgets from massive-star censuses; applied to Orion and to **Cygnus** — the closest thing to a prior supernova ledger for the prototype system, validated against the INTEGRAL aluminum-26 gamma-line map.
- **Zucker et al. (2022)** — the modern Gaia-era rerun for the Local Bubble (~15 supernovae) with proper cluster memberships.

**What is genuinely new in this project:** (1) DR4-grade inputs — memberships, per-star masses, subgroup ages, runaway tracebacks with epoch astrometry — changing the error budget qualitatively; (2) homogeneous application across *all* gamma-detected systems instead of one bubble; (3) the forward couplings nobody has built: supernova timeline → SNR occupancy, and supernova classification + measured L_w → per-system eta_CR. Familiar ingredients, one new combination — the safest kind of novelty: easy to defend, hard to dismiss, slow to steal because the value is in execution.

---

# PART II — READING CELLI ET AL. (2024) LIKE A REFEREE

*"Mass and wind luminosity of young Galactic open clusters in Gaia DR2", Celli, Specovius, Menchiari, Mitchell & Morlino, A&A 686, A118 (2024); arXiv:2311.09089. This is the closest existing work to the L_w part of the method above — and the clearest published demonstration of why the gamma-detected systems need the star-by-star treatment instead.*

## 1. What the paper actually does (mechanics)

Sample: start from the Gaia DR2 open-cluster catalog built by neural-network classification (~2017 clusters); keep those younger than 30 Myr → 387 usable clusters. Median age 14.5 Myr, median distance 2.3 kpc, median **58 observed member stars per cluster**.

For each cluster, the pipeline never fits individual massive stars. Instead:

1. Take the G-band magnitude distribution of the catalog members; kernel-fit it to extract two numbers — the brightest star's magnitude and the *mode* of the distribution (taken as the completeness edge).
2. Convert both magnitudes to bolometric (whole-spectrum) luminosities using temperature-dependent corrections, and correct for extinction using catalog values.
3. Convert the two luminosities into two masses via a mass–luminosity relation → the observed mass window [M*_min, M*_max].
4. Count the stars N* in that window; normalize a fixed four-segment IMF so its integral over the window equals N*.
5. Integrate the normalized IMF from 0.08 solar masses up to a maximum mass set jointly by the cluster age (older → lower ceiling, since heavier stars have died) and by an iterative cluster-mass/heaviest-star relation → **total cluster mass**, explicitly declared a *lower limit*.
6. Wind luminosity: integrate mass-loss-rate and wind-speed formulas *over the IMF* (an ensemble average, not a sum over actual stars), using the Nieuwenhuijzen & de Jager (1990) mass-loss prescription and an escape-velocity wind-speed scaling; then *add* observed Wolf-Rayet stars on top with one average mass-loss rate for all WN types and one for WC/WO types.

Note the philosophical difference from Part I: their unit of analysis is the *ensemble-average cluster of a given mass*; ours is *this particular realized cluster*. Both are valid — for different questions (see Section 10).

## 2. Figure 1 — the sample, and the "58 stars" reality check

Three panels: distributions of cluster age, distance, and number of member stars. Extract three numbers: median age 14.5 Myr, median distance 2.3 kpc, median 58 members. Internalize the last one: the typical cluster's entire mass estimate rests on ~58 counted stars and two magnitudes read off a histogram. That is population-statistics territory — perfectly sensible for 387 clusters at once, and structurally incapable of resolving any single cluster's individuality. Cygnus OB2 has *thousands* of measurable members; the regimes do not overlap. Also note panel (b): the flat distance distribution out to ~3 kpc is their completeness claim — the volume within which any statistical DR4 update would compete.

## 3. Figures 5 and 6 — the Wolf-Rayet self-audit (the paper's biggest money plot)

**The logic.** Wolf-Rayet stars are the evolved descendants of stars born above ~25 solar masses, in a phase lasting only ~0.25–0.4 Myr. Given a cluster's mass and age, the same IMF machinery that produced the mass estimate *predicts* how many WR stars the cluster should host right now. Figure 5 shows the prediction: WR numbers peak for cluster ages ~3–5 Myr and vanish after ~6.5 Myr (when even 25 solar-mass stars have died).

**The result.** Applying this to all 387 clusters with the catalog masses and ages predicts **zero** Wolf-Rayet stars in the entire sample. Observed, from the Galactic WR catalog: **49**, in 14 clusters (24 in Westerlund 1 alone). A more dramatic self-refutation is hard to design: the method's own extrapolation machinery, tested against the one directly countable product of the IMF's top end, misses by 49 to 0.

**The repair (Figure 6).** They ask: by what factor chi must all cluster masses be multiplied to reproduce the 49 observed WRs? Answer, depending on assumptions: chi =~ 4.8–5.7 (conservative WR physics, "Case A"), chi =~ 2.7 (generous, "Case B"), chi =~ 1.2 only if all cluster ages are forced artificially young ("Case C" — which their own statistical test disfavors). So the *published headline masses carry a known ×2.7–5.7 correction of uncertain size* — applied globally, as one number for all clusters.

**Why this is your exhibit A.** The top of the IMF is where the WR stars live — and it is also exactly where the wind power lives (Part I, Step 8). A method that misses the WR count by this margin misses L_w by a related margin, cluster by cluster. In conversation, the move is never "your method fails" but: *"your Figure 6 is what convinced me the gamma-detected systems need the star-by-star census — do you agree?"* — her own figure, her own conclusion, your project as its logical continuation.

## 4. Figure 4 and the age problem

They compare with an independent DR2 analysis (Almeida et al. 2023; isochrone fitting with synthetic clusters; 102 clusters in common). Distances agree within 2 kpc and diverge beyond. **Ages disagree systematically**: the catalog ages they inherited are underestimated for the youngest clusters (below ~6 Myr) and overestimated for older ones — and the WR analysis independently demands age revisions (some clusters with observed WR stars have catalog ages at which WR stars are impossible). Figure 4 itself shows the mass comparison (Almeida ~40% higher on average); the age discussion is in the surrounding text — notably, there is *no dedicated age figure*, which quietly tells you ages were inherited from the catalog, never derived.

**Why it matters to you:** age is the most sensitive input of the supernova ledger (Part I, Step 9, caveat 4). Their weakest parameter is your most critical one — which is why Phase 2 of your pipeline (own isochrone fits, per subgroup, two model families) is a wall, not wallpaper.

## 5. The wind prescription (Section 6 of the paper — read the equations, no figure)

Mass-loss rate from **Nieuwenhuijzen & de Jager (1990)** — a smooth global fit across the whole HR diagram, predating the modern line-driven-wind era; wind speed from escape velocity times a temperature-dependent factor between 1.0 and 2.65; cluster wind luminosity from a momentum-conservation average over the IMF. Wolf-Rayet stars added on top with a single average mass-loss rate (10^−4.9 solar masses/yr) for every WN star and fixed wind speeds (1600 km/s WN; 2300 km/s WC/WO).

Two observations. First, modern prescriptions (Vink 2001 vs Björklund 2021) differ from each other by ×2–3 and from NdJ 1990 by comparable factors — a systematic the paper cannot see because it runs one prescription. Your bracketing design quantifies what their single choice hides. Second, for a massive cluster the WR add-on can dominate the total: Westerlund 1's L_w is substantially "24 × an assumed constant." A per-star WR treatment from dedicated spectroscopy replaces that constant with measurements.

## 6. Figure 7 — the deliverable (and what any current correlation study would be built on)

Two panels: distributions of cluster mass-loss rate and wind luminosity across the 387 clusters. Medians: 1.3×10^−8 solar masses/yr and **L_w = 3×10^34 erg/s**, with a tail reaching ~3×10^38; with WR stars added, Danks 2 exceeds 10^39. This is the paper's product — the only homogeneous L_w compilation in existence — and it feeds their companion paper's gamma-ray/neutrino predictions. Anyone attempting an L_gamma vs L_w correlation today would draw the x-axis from this figure — with the ×2.7–5.7 chi-factor and the per-object unpredictability of Table 1 baked in.

## 7. Binaries and membership (Section 5 — two paragraphs, no figure, read verbatim)

Two self-declared factor-~2 systematics: (a) binaries ignored — if most stars are paired, star counts undercount mass by up to ×2; (b) the member lists come from one neural-network catalog, and alternative membership algorithms find up to ×2 more members. Both push masses and L_w *upward*, consistent with the WR verdict. These paragraphs are effectively the introduction of your paper, written by them: the limitations listed are precisely what a star-by-star DR4 census removes (probabilistic astrometric membership instead of one classifier; per-star binarity flags via RUWE and DR4 non-single-star tables).

## 8. Table 1 — the failure localized (your second money plot)

The six most massive clusters in their sample, against dedicated single-cluster literature studies:

| Cluster | Their mass (Msun) | Literature mass | Their age (Myr) | Literature age |
|---|---|---|---|---|
| Westerlund 1 | 22,227 | ~49,000 | 7.9 | 4.0 ± 0.5 |
| Danks 2 | 7,812 | 3,000 ± 800 | 1.0 | ~3 |
| NGC 3603 | 4,808 | 13,000 ± 3,000 | 1.0 | 1–2 |
| Danks 1 | 4,358 | 8,000 ± 1,500 | 1.0 | ~1.5 |
| NGC 6231 | 3,199 | 3,750 ± 450 | 13.8 | 2–7 |
| Westerlund 2 | 2,172 | 3,600 | 4.0 | 2.5 |

Three escalating lessons:

1. **These names are the gamma-ray sample.** Westerlund 1, Westerlund 2, NGC 3603, Danks 1 and 2 are on the list of firmly gamma-detected clusters. The table is the pipeline being tested on exactly the objects the efficiency question needs — not an abstract benchmark.
2. **The errors are large where they hurt.** Masses off by ×2–3; ages worse — Westerlund 1's age wrong by a factor 2 (7.9 vs 4.0 Myr) on the single most important TeV cluster in the sky. An age wrong by ×2 does not perturb the supernova ledger; it *changes its answer qualitatively* (turnoff mass, N_SN, time-since-last-SN all move together).
3. **The errors go in both directions — the decisive point.** Westerlund 1 and NGC 3603 are underestimated; Danks 2 is *over*estimated by ×2.6 with its age wrong the opposite way. Figure 6 showed a global *bias* (fixable by a global chi factor — which is what they apply). Table 1 shows per-object errors **unpredictable in sign and size** — and no global correction can fix unpredictability. For a population sum, tolerable; for per-system efficiencies, fatal. Table 1 is the one-table proof that the universality-of-eta_CR question *cannot* be answered with the statistical method — through no fault of the paper's execution; it is a regime mismatch, not a mistake.

## 9. The companion consumer

"Probing Stellar Clusters from Gaia DR2 as Galactic PeVatrons I" (arXiv:2403.16650) feeds these masses and wind luminosities into the wind-termination-shock acceleration model to predict gamma-ray and neutrino fluxes for the cluster population. Read it with two questions: which named clusters are predicted detectable (and does the list match the actually-detected ones — mismatches are diagnostic), and how do predicted fluxes scale with the input L_w (that scaling is the lever your improved values would move). Strategic meaning: her group operates a prediction machine whose fuel is exactly the quantity the star-by-star pipeline produces better. You are not building a rival machine; you are refining its fuel.

## 10. The regime distinction, condensed

| | Celli et al. (2024) | This project |
|---|---|---|
| Question | Population-level: total wind power of the Galactic cluster population | Per-system: eta_CR of individual gamma-detected systems |
| Unit of analysis | Ensemble-average cluster of given mass | This particular realized cluster |
| Massive stars | IMF-integrated expectation values | Actual observed stars, individually |
| Ages | Inherited from DR2 catalog | Fitted per subgroup, two model families |
| Wind physics | One 1990 prescription | Vink 2001 vs Björklund 2021 bracket |
| WR stars | One average value per type | Per-star from dedicated spectroscopy |
| Supernova history | Absent | Central deliverable |
| Errors | Global chi factor ×2.7–5.7 | Per-quantity posteriors via Monte Carlo |
| Right for | 387 clusters, population totals | 10–15 systems, efficiency scatter |
| Sample includes Cyg OB2? | No — associations excluded by construction | Yes — the prototype |

**One sentence to keep: statistical where the numbers are, star-by-star where the power is.**

## 11. Money plots, ranked for the conversation

1. **Figure 6** (the WR chi-factor) — her own proof that the statistical method misses the top of the IMF by ×3–6, which is where the wind power lives.
2. **Table 1** (massive-cluster comparison) — the failure localized to exactly the gamma-relevant systems, in both directions, hence uncorrectable globally.
3. **Figure 7b** (the L_w distribution) — the deliverable your pipeline supersedes; the x-axis of any premature correlation study.
4. **Figure 4** (Almeida comparison) — the age unreliability: your supernova-ledger risk, documented by them.
5. **Figure 1** (sample properties) — the 58-stars-per-cluster reality check separating her regime from yours.

Internalize Figures 6 and Table 1 to the point of fluent discussion, and you will read as someone who studied the paper more carefully than most referees — the single strongest impression a junior researcher can make on a potential collaborator.

---

# APPENDIX — GLOSSARY

- **Parallax** — apparent annual positional wobble of a star due to the observer's orbit; inverse of distance. 1 milliarcsecond (mas) <-> 1 kpc.
- **Proper motion** — angular drift of a star across the sky, mas/yr.
- **Radial velocity (RV)** — line-of-sight velocity from the Doppler shift of spectral lines.
- **RUWE** — Gaia's astrometric goodness-of-fit; >1.4 flags unreliable solutions (often unresolved binaries).
- **CMD (color–magnitude diagram)** — stellar color (temperature proxy) vs brightness; the observational Hertzsprung–Russell diagram.
- **Isochrone** — model curve of where stars of one common age but different masses lie on the CMD.
- **Main-sequence turnoff** — CMD location where stars are currently leaving the main sequence; a clock hand for cluster age.
- **Turnoff mass m_TO(t)** — the mass whose lifetime equals the cluster age; all heavier stars are dead.
- **IMF (initial mass function)** — distribution of stellar birth masses, xi(m) = k m^−2.3 above 0.5 solar masses; shape universal, normalization k per cluster.
- **CCSN (core-collapse supernova)** — terminal explosion (~10^51 erg) of a star born above ~8 solar masses.
- **SNR (supernova remnant)** — the expanding shocked shell left behind, radio/X-ray visible for ~10^4–10^5 years.
- **Wolf-Rayet (WR) star** — evolved massive star that shed its hydrogen envelope; extreme wind (M-dot ~10^−5 solar masses/yr, up to 5000 km/s); ~0.25–0.4 Myr phase of stars born above ~25 solar masses.
- **M-dot, v_inf** — mass-loss rate and terminal wind velocity; single-star wind power = (1/2) M-dot v_inf².
- **L_w** — cluster wind luminosity: sum of member wind powers.
- **Prescription** — calibrated formula giving M-dot from stellar parameters (Vink 2001; Björklund 2021; Nieuwenhuijzen & de Jager 1990).
- **Clumping** — small-scale wind inhomogeneity; inflates classical M-dot estimates by ×2–3; the main wind-physics uncertainty.
- **Extinction (A_V, A_G)** — dust dimming; reddening is the associated color change; R_V parameterizes the dust law.
- **eta_CR** — CR acceleration efficiency, =~ 3 L_gamma / L_w under the calorimetric assumption.
- **Calorimetric assumption** — accelerated protons lose essentially all energy to pion production locally; if violated, computed eta_CR is a lower limit.
- **HDBSCAN** — density-based clustering algorithm; finds clumps without preset number or shape.
- **Poisson likelihood** — correct statistical treatment for small counts (bin scatter = counting noise).
- **Monte Carlo** — uncertainty estimation by simulating many random realizations through the full pipeline.
- **ADQL / TAP** — the SQL dialect and web protocol of the Gaia archive; queries run server-side.
- **Runaway star** — massive star ejected from its birthplace (companion's supernova kick or dynamical slingshot); must be traced back and restored to the census.
- **Coeval** — formed at the same time; the working approximation for cluster subgroups (real spread ~1–3 Myr).
