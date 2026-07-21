# WP0 — Prior-art and IMF-deficit summary

**Completed:** 2026-07-20  
**Purpose:** establish what has already been done, extract reusable methods, and define the exact novelty boundary for Paper 1.

## 1. Prior-art map

| Paper | Data / target | Main method | SN-history result relevant to Paper 1 | What it does not provide |
|---|---|---|---|---|
| Knödlseder et al. 2002, [astro-ph/0206045](https://arxiv.org/abs/astro-ph/0206045) | 2MASS-era Cyg OB2 census, radioactive isotopes, evolutionary synthesis | Bayesian population synthesis with Geneva tracks, coeval/continuous formation alternatives, isotope yields | Cyg OB2 is modeled as young (about 1–4 Myr), with first SNe near 4 Myr in the coeval picture; weak/absent `60Fe` indicates little recent activity | No Gaia membership, per-star extinction, probabilistic distance mixture, or explicit last-SN posterior |
| Martin et al. 2010, [arXiv:1001.1522](https://arxiv.org/abs/1001.1522) | Cygnus complex; Cyg OB2 row uses 120 stars in 20–120 M☉, 1584 pc, 2.5 Myr | Monte Carlo IMF sampling until the observed population criterion is met; instantaneous burst; stellar evolution in 10 kyr steps | Typically 10–20 SNe over the last Myr for the **whole Cygnus complex**; links the young pulsar/runaway picture to possible past activity | The result is not a Cyg OB2-only Gaia ledger and does not model the planned per-star extinction/distance/membership uncertainty |
| Wright et al. 2015, [arXiv:1502.05718](https://arxiv.org/abs/1502.05718) | Cyg OB2 OB census and IMF | Spectroscopic/photometric census with IMF and total-mass inference | 169 primary OB systems, 52 O stars, 3 WR stars; total mass about `16,500^{+3,800}_{-2,800} M☉`; continuous SF broadly 1–7 Myr | It is the key census baseline but not a Gaia DR3 probabilistic dead-star/SN-history analysis |
| Berlanas et al. 2019, [MNRAS 484, 1838](https://academic.oup.com/mnras/article/484/2/1838/5289922) | Gaia DR2 astrometry for known OB members | Unbinned parallax mixture modeling with MCMC and BIC | Main group near 1760 pc and foreground group near 1350 pc; foreground component about 19% of the modeled sample and about 10% of the cited mass scale | Does not convert the distance-resolved membership structure into an IMF-normalized SN history |
| Berlanas et al. 2020, [arXiv:2008.09917](https://arxiv.org/abs/2008.09917) | 78 O stars with spectroscopy and Gaia DR2 | Stellar parameters, ages, multiplicity and proper-motion analysis | At least two formation episodes near 3 and 5 Myr; 9 O stars have peculiar motions / possible ejection histories; 29/78 known binary or multiple | No stochastic SN ledger or joint last-SN test of the H25 gamma-ray hypothesis |
| Fuchs et al. 2006, [astro-ph/0609227](https://arxiv.org/abs/astro-ph/0609227) | Sco-Cen UCL/LCC surviving A0–B0 stars | Normalize an IMF with surviving stars, infer the turnoff mass from lifetimes, integrate the missing high-mass population | 8–12 SNe in UCL plus 6–8 in LCC, 14–20 total; derives a time-dependent SN rate from the IMF and lifetimes | Not Cygnus; deterministic/coeval; no Gaia probabilistic membership, binaries, failed-SN branch, or per-star extinction model |
| Zucker et al. 2022, [arXiv:2201.05124](https://arxiv.org/abs/2201.05124) | UCL/LCC superbubble shell | Analytic superbubble dynamics and momentum budget | `15^{+11}_{-7}` SNe required by the shell momentum, consistent with the F06 14–20 estimate | A dynamical cross-check, not an IMF census or Gaia-based demographic reconstruction |
| Menchiari et al. 2024, [arXiv:2402.07784](https://arxiv.org/abs/2402.07784) | Cyg OB2 wind/CR environment and mock populations | Wind prescriptions plus mock stellar populations | Expected SNe `7 ± 2.5` at 3 Myr and `26 ± 5` at 5 Myr; adopted total mass near `16,500 M☉` | No Gaia DR3 member-level mass/extinction ledger; values depend on adopted age and mock-population assumptions |
| Hä​rer et al. 2025, [arXiv:2508.21644](https://arxiv.org/abs/2508.21644) | Cygnus gamma rays and molecular environment | 3D gas/hydrodynamic gamma-ray modeling | Preferred interpretation is a powerful Cyg OB2 SN about 50 kyr ago; lower age scale about 10–20 kyr; energy requirement is several `10^51 erg` for the assumed efficiency | Tests the gamma-ray morphology/spectrum, not whether the Cyg OB2 stellar population statistically produced the event |

K02’s compact warning is useful context: “the small number of recent supernova events suggests only little 60Fe production” (§6.1). It should be treated as an isotope constraint with large model dependence, not as a direct count.

## 2. The reusable IMF-deficit template

The transferable logic appears most clearly in F06 and has related forms in K02 and Martin:

1. **Select a surviving calibration population.** Use a mass interval below the current turnoff where the census is expected to be reasonably complete. F06 uses approximately `2.6–8.2 M☉` A0–B0 stars.
2. **Specify the IMF convention.** F06 writes the IMF in a `Γ` convention; Paper 1 should use `dN/dM ∝ M^{-α}` and report the conversion explicitly. The planned branches are `α = 2.0, 2.3, 2.6`.
3. **Normalize with the observed survivors.** Estimate the normalization from the calibration stars, including membership and completeness uncertainty. F06 uses a deterministic normalization; Martin samples the IMF until the observed high-mass criterion is met.
4. **Infer the turnoff mass from an age/lifetime relation.** The dead-star boundary is not a fixed mass. It changes with age, rotation/isochrone family, metallicity, and the assumed duration of star formation.
5. **Integrate the missing high-mass population.** The basic count is proportional to `∫_{m_TO}^{m_max} ξ(m) dm`, but Paper 1 must separate all stellar deaths from successful explosions.
6. **Map mass to lookback time.** Use `τ(m)` to construct a rate/history, convolving with measurement, age-duration, and stellar-evolution uncertainty. A coeval burst produces a sharp mapping; an extended burst broadens it.
7. **Add physical corrections.** Binary evolution and mergers can create stripped stars or change lifetimes; runaways can move progenitors out of the present association; failed explosions reduce the successful-SN count; incompleteness and unresolved multiplicity affect the normalization.
8. **Propagate stochasticity.** A high-mass tail contains few objects, so a Poisson or posterior-predictive treatment is essential. A single plug-in count is not sufficient for the H25 plausibility test.

The F06 evidence can be summarized in one short excerpt: “The expected number of supernovae, i.e. the number of ‘missing’ stars” (§3.2). Paper 1 should retain that intuitive interpretation while replacing the deterministic census with a Gaia-era probabilistic ledger.

## 3. Method comparison against the planned Paper 1

| Planned Paper 1 element | Closest precedent | Novelty status |
|---|---|---|
| Gaia DR3 membership with probabilistic subgroup assignment | B19 Gaia DR2 distance mixture | Incremental data release plus a new downstream SN-history use; not a duplicate by itself |
| Per-star extinction rather than a regional mean | B20 spectroscopy and earlier photometric work | Methodologically important; must be shown to change mass/turnoff inference |
| PARSEC/MIST age and stellar-mass branches | K02/Martin evolutionary synthesis; B20 age analysis | Robustness design, not standalone novelty |
| IMF normalization from surviving members with completeness and Poisson uncertainty | F06 missing-star method; Martin Monte Carlo sampling | Directly related prior art; novelty is its Gaia DR3, subgroup-resolved implementation and uncertainty ledger |
| Last-SN lookback posterior | F06 time-dependent rate; H25 fixed/preferred gamma-ray event age | A potentially novel synthesis if it is derived from the stellar population and explicitly compared with H25 |
| Successful versus failed SN history | Prior studies discuss evolution and yields but do not provide the planned Cyg OB2 branch structure | Strong differentiator if implemented transparently and not overclaimed |
| Joint age–energy–progenitor–location plausibility | H25 provides the gamma-ray side; W15/B20 provide stellar context | Core scientific contribution: test whether the proposed event is demographically plausible |

## 4. Precise novelty boundary

The defensible claim is:

> Paper 1 provides a Gaia DR3, distance-subgroup-aware, probabilistic reconstruction of the Cyg OB2 massive-star population and converts it into a stochastic, branch-resolved history of stellar deaths and successful SNe, including a posterior for the last event, then tests the demographic plausibility of the H25 gamma-ray hypothesis.

Avoid these broader claims:

- “First evidence for a past supernova in Cyg OB2.” Wright 2015 already discusses a young pulsar/runaway interpretation.
- “First SN-rate estimate for Cygnus.” Martin 2010 already gives a Cygnus-complex rate/count.
- “First Gaia study of Cyg OB2.” B19 and B20 are Gaia DR2 studies.
- “The H25 event is confirmed.” H25 is a model interpretation; Paper 1 can test demographic plausibility, not establish the gamma-ray association alone.

## 5. Prior-art conclusions for WP1

The literature supports the proposed project because the pieces exist separately but not as the same analysis: Gaia distance/membership work, spectroscopic O-star work, population synthesis, IMF-deficit counting, and gamma-ray modeling are distributed across different papers and datasets. The remaining scientific gap is the calibrated connection between a Gaia-era Cyg OB2 member ledger and the probability distribution of recent successful SNe, especially after accounting for distance substructure, per-star extinction, age duration, binaries, runaways, and failed explosions.

The closest 2024–2026 data-level overlap found in the expanded search is Paíz et al. (2025), [A New Open Cluster Census in the Region of Cygnus OB2 with Gaia Data](https://doi.org/10.4236/ijaa.2025.152013). It uses Gaia DR3 and HDBSCAN to identify open clusters and astrometric members in the wider region, with age and reddening estimates, but it does not perform the massive-star IMF-deficit calculation or reconstruct an SN history. Paper 1 should cite it and explain the different sample definition and scientific endpoint.
