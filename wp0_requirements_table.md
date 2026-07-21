# WP0 — Requirements and evidence table

**Study:** Gaia SNR History in Cygnus / Paper 1  
**WP0 completed:** 2026-07-20  
**Scope:** translate the scientific question into auditable requirements, extract the stellar-population assumptions used by the main prior-art papers, and establish whether the planned Gaia-based analysis is already published.

## WP0 gate result

**Requirements table:** complete.  
**Prior-art assumptions:** extracted and cross-checked against the cited papers.  
**Duplicate search:** no direct duplicate identified as of 2026-07-20.  
**Action:** proceed to WP1, while repeating the ADS/arXiv search immediately before manuscript submission.

## Source register and citation correction

| ID | Reference | Role in WP0 |
|---|---|---|
| H25 | Hä​rer et al. 2025, [arXiv:2508.21644](https://arxiv.org/abs/2508.21644), [A&A DOI](https://doi.org/10.1051/0004-6361/202555531) | Direct gamma-ray/SNR hypothesis for Cygnus OB2; required physical plausibility test |
| M24 | Menchiari et al. 2024, [arXiv:2402.07784](https://arxiv.org/abs/2402.07784), [A&A 686 A242](https://doi.org/10.1051/0004-6361/202348817) | Recent wind/CR model and explicit expected-SN numbers |
| K02 | Knödlseder et al. 2002, [astro-ph/0206045](https://arxiv.org/abs/astro-ph/0206045) | Population-synthesis and radioactive-isotope prior |
| Ma10 | Martin et al. 2010, [arXiv:1001.1522](https://arxiv.org/abs/1001.1522) | Cygnus-wide population synthesis and recent-SN-rate prior |
| W15 | Wright et al. 2015, [arXiv:1502.05718](https://arxiv.org/abs/1502.05718) | Modern Cyg OB2 census, IMF and mass normalization |
| B19 | Berlanas et al. 2019, [MNRAS 484, 1838](https://academic.oup.com/mnras/article/484/2/1838/5289922) | Gaia DR2 distance substructure and membership context |
| B20 | Berlanas et al. 2020, [arXiv:2008.09917](https://arxiv.org/abs/2008.09917), [A&A 642 A168](https://doi.org/10.1051/0004-6361/202038642) | Spectroscopic O-star parameters, age bursts, binaries and runaways |
| F06 | Fuchs et al. 2006, [astro-ph/0609227](https://arxiv.org/abs/astro-ph/0609227), [MNRAS DOI](https://doi.org/10.1111/j.1365-2966.2006.11044.x) | Reusable IMF-deficit / missing-high-mass-star counting template |
| Z22 | Zucker et al. 2022, [arXiv:2201.05124](https://arxiv.org/abs/2201.05124) | Independent superbubble-momentum cross-check on an inferred SN count |

The execution plan lists Menchiari as `arXiv:2306.00946`; that identifier is not the cited 2024 Cyg OB2 wind/CR paper. The paper used here is `arXiv:2402.07784`.

## A. Hörer et al. requirement extraction

The H25 paper is a gamma-ray interpretation, not a stellar-demography ledger. Its result is therefore a **hypothesis to test**, not an input count to adopt.

| Requirement to carry into Paper 1 | Evidence extracted from H25 | Operational test / WP9 consequence |
|---|---|---|
| Explosion energy must be physically plausible | The preferred model requires a powerful event. H25’s discussion gives an order of `3 × 10^51 erg` for 10% proton-acceleration efficiency; its hydrodynamic setup uses a higher, type-Ic-like `5 × 10^51 erg` injection. | Report an energy branch/range and compare it with the inferred progenitor mass/type and explodability branch. Do not silently equate the simulation energy with a measured SN energy. |
| Progenitor type must be consistent with the population | H25 argues that, at 3–5 Myr, type-Ic explosions are the most likely among its modeled outcomes and that their energies may exceed the canonical `10^51 erg`. | Condition the plausibility statement on the progenitor channel. Mark the result as weak/conditional if only an ordinary core-collapse event is supported. |
| The event age must match the gamma-ray morphology and the surviving population | Preferred age is about 50 kyr; the lower bound is roughly 10–20 kyr because an older remnant should have faded. H25 also allows an event within the last few hundred kyr in its conclusion. | Produce a posterior for the last-SN lookback time and compare it with a 10–20 kyr lower bound, a ~50 kyr preferred value, and the broader “few hundred kyr” allowance. |
| The location must be in the relevant low-density cavity/core | H25 places its preferred event in a low-density, roughly 150 pc superbubble associated with Cyg OB2; it distinguishes this from the offset Gamma Cygni SNR. | Use Gaia-defined subgroup positions and a cavity/core geometry test. State whether the inferred progenitor belongs to the central association or to a line-of-sight substructure. |
| The hypothesis must be population-plausible | Cyg OB2 contains more than 70 O stars and 3 WR stars in H25’s adopted context, with heterogeneous ages around 3–7 Myr and distance groups near 1.3, 1.5 and 1.7 kpc. | Recompute the expected progenitor pool from the Gaia DR3 membership ledger rather than importing the H25 census. Carry distance-mixture and age-mixture uncertainty into the final probability. |

**Verbatim evidence from H25** (Conclusion, §8): “a powerful supernova, which exploded about 50 kyr ago in the stellar association Cygnus OB2 is the best candidate”. This is the paper’s preferred interpretation, not a direct detection of a remnant or progenitor.

## B. Stellar-population assumptions to audit

| Quantity | Prior-art values / assumptions | Consequence for the planned analysis |
|---|---|---|
| Distance | Hanson: about `1400 ± 80 pc`; B19 resolves a foreground group near `1350 pc` and a main group near `1755 pc`, with the foreground component about 19% of the modeled OB sample; Wright/Martin use values near 1.6 kpc. | Distance must be a per-star or probabilistic subgroup variable. A single 1.4–1.6 kpc distance can bias luminosity, mass, and therefore the inferred number of dead massive stars. |
| Present-day O-star count | Historical estimates range from about 75–120 O stars depending on sky area, selection, and whether a value is an upper limit. Wright’s core census has 52 O stars among 169 primary OB systems; B20’s spectroscopic compilation contains 78 O stars. | Define the Paper 1 sample explicitly and report completeness. Do not compare unlike counts without matching the footprint, membership probability, multiplicity treatment, and mass range. |
| Total stellar mass | W15 gives `16,500^{+3,800}_{-2,800} M☉` for Cyg OB2. M24 adopts the same scale in its mock populations. | Use this as a prior/check, not as the sole normalization. Compare the IMF-normalized mass from the Gaia sample with the literature value and quantify the effect of excluded low-mass members. |
| Age | Literature spans roughly 2–3 Myr for non-rotating isochrones, 4–5 Myr for rotating interpretations, and an extended/continuous formation interval of about 1–7 Myr. B20 identifies at least two bursts near 3 and 5 Myr. | Age is load-bearing because it sets the turnoff mass and the SN clock. Require PARSEC/MIST and age-duration branches, plus spectroscopic anchors. |
| O-star / WR census and evolution | K02 uses a 1–4 Myr Cyg OB2 age, about 120 O stars as an upper limit, 3 WR stars, and a high-mass interval extending to roughly 120 M☉. Its synthesis places first SNe around 4 Myr after a coeval birth. | Use WR/O-star counts as validation diagnostics, not as a deterministic SN count. Test whether the adopted age branch can produce the observed evolved population. |
| IMF slope | K02 reports a high-mass logarithmic slope near `Γ = −1.1 ± 0.3`; W15 reports an IMF close to the Salpeter-like slope in its convention. F06 fits `Γ = −1.1 ± 0.1`, equivalent to `dN/dM ∝ M^{-2.1}`. | Implement the planned `α = 2.0, 2.3, 2.6` branches and make the slope convention explicit. Never compare `Γ` and `α` without converting definitions. |
| IMF calibration interval | F06 normalizes its IMF with surviving A0–B0 stars, approximately `2.6–8.2 M☉`, then integrates above the turnoff. Martin instead samples until the observed high-mass criterion is met. | Paper 1 should expose the calibration interval, completeness correction, and likelihood. The calibration window is part of the scientific result, not a hidden implementation detail. |
| Binary/multiple fraction | B20 finds 29 of 78 O stars in known binaries/multiples, while the true fraction is likely higher; it cites a roughly 45–55% spectroscopic-binary expectation. Martin notes that unresolved multiples affect the inferred IMF. | Include the planned binary sensitivity branch. State whether masses/counts refer to systems or individual stars and how mergers/stripped stars alter the dead-star ledger. |
| Runaways / ejections | B20 identifies 9 O stars with peculiar proper motions or possible ejection histories. Martin links a nearby young pulsar/runaway picture to a possible first SN. | Include a runaway correction or sensitivity test and do not assign every missing high-mass star to an in-situ SN. |
| Failed SNe / explodability | The prior papers do not provide a Gaia-calibrated failed-SN posterior for Cyg OB2. H25’s gamma interpretation needs an energetic successful event, while population synthesis counts can include progenitors that do not leave an observable SN. | Retain the planned all-explode versus Sukhbold-like failed-island branches and report successful explosions separately from all deaths. |

## C. Predicted SN rates and recent-event counts

| Source | Quantity extracted | Interpretation |
|---|---|---|
| K02 | Coeval models place the first Cyg OB2 SNe around 4 Myr after birth; the paper emphasizes low recent SN activity through the weak/absent `60Fe` prediction. | A broad evolutionary prior, not a direct recent-SN count. Its age constraint is strongly model-dependent. |
| Martin 2010 | For the **whole Cygnus complex**, the simulation predicts typically 10–20 SNe over the last Myr. The Cyg OB2 row is normalized to 120 stars in 20–120 M☉ at 1584 pc and 2.5 Myr. | Do not quote 10–20 as a Cyg OB2 count. It is a complex-wide output from an instantaneous-burst, Monte Carlo population model. |
| Menchiari 2024 | Their mock Cyg OB2 populations give `7 ± 2.5` SNe at 3 Myr and `26 ± 5` at 5 Myr, with a total mass scale near `16,500 M☉`. | Useful recent benchmark, but tied to a wind/CR model, adopted ages, and mock populations rather than a Gaia DR3 member-by-member census. |
| Fuchs 2006 | An IMF-deficit calculation gives 8–12 SNe in UCL and 6–8 in LCC over the relevant lookback intervals, 14–20 total, using surviving B/A stars and single-star lifetimes. | This is the closest published counting template, but it targets Sco-Cen associations and is deterministic/coeval rather than Gaia-probabilistic. |
| Zucker 2022 | Superbubble momentum independently requires `15^{+11}_{-7}` SNe for UCL/LCC, broadly consistent with the F06 IMF estimate. | Use as a methodological cross-check: dynamics can validate an IMF-based total, but it does not replace the Paper 1 stellar ledger. |

Martin’s relevant evidence is: “Our simulation predicts that a certain number of SNe, typically 10-20 over the last Myr, exploded in the Cygnus complex” (§3.4). The geographic qualifier is essential.

## D. Requirements for the Paper 1 result

The WP0 audit makes the following requirements explicit:

1. Define the Gaia DR3 membership sample, footprint, distance-mixture treatment, and system-versus-star counting convention.
2. Infer per-star extinction and luminosity/mass; do not use a regional mean extinction as the primary measurement.
3. Use spectroscopic O-star parameters and the WR/O-star census as anchors and validation checks.
4. Fit/normalize the IMF in a declared surviving-star interval, with completeness and Poisson/stochastic uncertainty visible.
5. Convert the inferred turnoff/dead-star population to a successful-SN history using explicit lifetime, binary, runaway, and failed-SN branches.
6. Report both an all-deaths ledger and a successful-explosion ledger; these answer different questions.
7. Test the H25 hypothesis through a joint plausibility statement involving age, energy, progenitor channel, subgroup location, and last-SN lookback time.
8. Keep the following branches in the reported sensitivity envelope: PARSEC/MIST; `α = 2.0/2.3/2.6`; `R_V = 3.1/3.0/3.5`; all-explode/Sukhbold-like failed islands; and SF duration 0/1/2 Myr.

## E. WP0 acceptance checklist

- [x] H25 energy, progenitor, age, location, and population-plausibility requirements extracted.
- [x] Menchiari, Knödlseder, Martin, Wright, Berlanas, Fuchs, and Zucker assumptions recorded.
- [x] The Menchiari arXiv identifier in the execution plan corrected for this audit.
- [x] Recent-SN values distinguished between Cyg OB2 and the whole Cygnus complex.
- [x] IMF-deficit implementation choices identified and mapped to Paper 1 decisions.
- [x] Duplicate search expanded to 2024–2026 ADS-style citation/title searches, arXiv, ICRC and TeVPA materials, plus all-time exact/near-exact keyword searches.
- [x] No direct duplicate identified as of 2026-07-20; residual search limitations documented in [wp0_dedup_report.md](wp0_dedup_report.md).
