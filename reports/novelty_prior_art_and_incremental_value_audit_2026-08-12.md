# Cygnus OB2 supernova-history study: prior art, non-novel results, and a path to demonstrable incremental value

**Date:** 2026-08-12  
**Status:** candid internal scientific audit, not manuscript prose  
**Scope:** the massive-star census, inferred stellar-death history, conversion of deaths into successful supernovae, and the value of resolving Cygnus OB2 into subgroups.

## Executive verdict

The project is technically substantial, but technical complexity is not itself scientific novelty. Several statements currently emphasized in the manuscript are already known, obvious consequences of known physics, or direct repetitions of earlier calculations:

- a roughly 4 Myr stellar population can have lost only very massive stars;
- the number of past deaths is extremely sensitive to an assumed age near the first-supernova boundary;
- an IMF normalized by surviving stars can be integrated above the turnoff to estimate missing/dead massive stars;
- not every massive-star core collapse necessarily produces a successful luminous supernova;
- Cygnus OB2 may already have experienced a supernova;
- Cygnus OB2 has age structure and line-of-sight/spatial substructure;
- mock-population estimates of the number of past Cygnus OB2 supernovae already exist.

The strongest remaining candidate contribution is narrower:

> **Because the death-count function is strongly nonlinear at ages of approximately 2--5 Myr, replacing the resolved star-formation history of Cygnus OB2 by one association-wide age may bias both the inferred number and timing of stellar deaths. A subgroup-resolved Gaia analysis could identify the responsible population and quantify when a simple fixed-age calculation is adequate.**

That statement is currently **plausible but not yet demonstrated**. The existing comparison against a forced common age of 4 Myr is suggestive, but the required fair comparison against the best-fitting one-population model has not been performed. Until that ablation is completed, the novelty should be described as a hypothesis under test rather than a result.

## 1. What the main prior studies already did

| Study | What it already established or calculated | What overlaps this project | What it did not do |
|---|---|---|---|
| [Wright, Drew & Mohr-Smith (2015)](https://arxiv.org/abs/1502.05718), *The Massive Star Population of Cygnus OB2* | Compiled 169 primary OB systems, inferred stellar masses and ages, assessed the star-formation history and high-mass mass function, estimated a total mass of about `16,500 (+3,800/-2,800) Msun`, and argued from the age distribution, mass-function deficit, a pulsar, and a runaway O star that Cygnus OB2 may already have experienced its first supernova. It explicitly discusses the loss of the most massive stars from an older population. | Massive-star census, IMF, age distribution, evolved/missing high-mass stars, total mass, and evidence for a previous supernova. | No Gaia DR3 probabilistic subgroup-to-death ledger or posterior distribution for the number/time of deaths. |
| [Berlanas et al. (2019)](https://arxiv.org/abs/1901.02959), *Disentangling the spatial substructure of Cygnus OB2 from Gaia DR2* | Used Gaia DR2 parallaxes and a mixture model to identify a main group near 1760 pc and a foreground group near 1350 pc, with individual membership probabilities. | The premise that Cygnus OB2 is not one homogeneous distance population and that Gaia should be used probabilistically. | No IMF-normalized stellar-death or supernova history. |
| [Berlanas et al. (2020)](https://arxiv.org/abs/2008.09917), *Spectroscopic characterization ... Evidence of multiple star-forming bursts* | Analysed 78 O stars, derived stellar parameters and masses for stars with reliable Gaia astrometry, found at least two formation episodes near 3 and 5 Myr, and discussed peculiar motions and possible past supernovae. | Non-coeval star formation, subgroup/episode-dependent ages, spectroscopy, multiplicity, runaways, and possible past explosions. | No Gaia DR3 low/intermediate-mass normalization connected to a stochastic death-time ledger. |
| [Knödlseder et al. (2002)](https://arxiv.org/abs/astro-ph/0206045), *Gamma-ray line emission ... The Cygnus region* | Performed Bayesian evolutionary population synthesis for Cygnus OB associations and clusters, propagating age, distance, and richness uncertainties; predicted isotope and massive-star feedback histories and concluded that few recent supernovae imply little `60Fe`. | Age-dependent massive-star evolution and supernova activity in Cygnus with uncertainty propagation. | Based on heterogeneous pre-Gaia populations; not a Gaia DR3, Cygnus-OB2-only subgroup ledger. |
| [Martin et al. (2010)](https://arxiv.org/abs/1001.1522), *Predicted gamma-ray line emission from the Cygnus complex* | Used Monte Carlo IMF population synthesis and time-dependent stellar evolution for the Cygnus complex. It reports approximately 10--20 supernovae over the last Myr for the **whole complex**, not Cygnus OB2 alone. | Monte Carlo population synthesis, time-dependent SN histories, and sensitivity to population inputs. | Not a Gaia-selected Cygnus-OB2-only reconstruction and not subgroup resolved. |
| [Fuchs et al. (2006)](https://arxiv.org/abs/astro-ph/0609227), *The search for the origin of the Local Bubble redivivus* | Normalized an IMF using surviving stars, used ages/lifetimes to determine the missing high-mass population, estimated 14--20 past supernovae, and reconstructed the association trajectories. | The central "missing stars above the turnoff" or IMF-deficit logic. This is a direct methodological precedent. | Different associations and older data; no Cygnus-specific Gaia analysis. |
| [Menchiari et al. (2024)](https://arxiv.org/abs/2402.07784), *Cygnus OB2 as a test case for particle acceleration in young massive star clusters* | Generated 1000 mock populations using the Wright et al. total mass and an assumed IMF, removed stars whose main-sequence lifetime was shorter than the cluster age, and obtained `7 +/- 2.5` supernovae at about 3 Myr and `26 +/- 5` at about 5 Myr. | A Cygnus-OB2-specific mock-population estimate of the number of stars already dead/exploded. Its computational skeleton is close to the present death-count ledger. | It did not infer the member-level normalization, extinction, ages, or subgroups from Gaia DR3 and did not reconstruct a subgroup-resolved event-time distribution. Its main scientific target was the wind/particle-acceleration model. |
| [Sukhbold et al. (2016)](https://arxiv.org/abs/1510.04643) and [Ertl et al. (2016)](https://arxiv.org/abs/1503.07522) | Showed that explodability depends on presupernova structure and is not a simple monotonic function of ZAMS mass; some collapsing progenitors fail and form black holes. | The distinction between stellar death/core collapse and a successful supernova. | They do not apply that physics to the Gaia-inferred Cygnus OB2 population. |
| [Zapartas et al. (2017)](https://arxiv.org/abs/1701.07032) | Demonstrated with binary population synthesis that interactions alter the number and delay-time distribution of core-collapse supernovae, including late events. | Shows why a coeval single-star death ledger cannot be treated as a complete physical supernova history. | No Cygnus OB2 application. |
| [Härer et al. (2025)](https://arxiv.org/abs/2508.21644) | Modelled the Cygnus gamma-ray emission and favoured, under its assumptions, a powerful recent supernova interpretation. | Supplies an external high-energy event hypothesis that a stellar-demographic history can test. | Does not derive the event probability from a Gaia stellar census. |

The project's own prior-art sweep found no paper combining all of Gaia-era membership, extinction/mass inference, subgroup-resolved IMF normalization, a stochastic past-event history, and a Cygnus high-energy comparison as of 2026-07-30 ([local sweep](wp0_dedup_resweep_2026-07-30.md)). This means that no exact duplicate was found. It does **not** by itself prove that combining those pieces changes a scientific conclusion or is worth a full journal paper.

## 2. The Menchiari comparison: what is true and what it means

The manuscript's sentence about Menchiari et al. is factually supported. In Appendix B, Menchiari et al. adopt the Wright et al. mass, simulate mock stellar populations, remove stars with main-sequence lifetimes shorter than the assumed age, and report:

- approximately `7 +/- 2.5` supernovae at 3 Myr;
- approximately `26 +/- 5` supernovae at 5 Myr.

The ratio is `26/7 = 3.7`, so "approximately a factor of four across a 2 Myr age change" is accurate. The result is not surprising: this age interval lies near the onset of the supernova era, where the turnoff moves rapidly through the sparse upper IMF.

This prior work means that the following cannot be claimed as new:

- predicting how many Cygnus OB2 massive stars have already died from an assumed mass, IMF, age, and lifetime relation;
- demonstrating that the inferred count changes strongly between 3 and 5 Myr;
- obtaining a count of order several to several tens.

The current baseline gives `N_SN | all-explode = 8.43`, but numerical proximity to Menchiari's `7 +/- 2.5` is not yet evidence that the complicated model is superior. The values are not directly like-for-like: the current analysis uses different normalization conventions, mixed subgroup ages, different evolutionary grids, and explicit branches. A controlled reproduction of Menchiari's assumptions is needed before interpreting agreement or disagreement.

## 3. Results that are not scientific novelty

### 3.1 Only very massive stars could have died

For a young population, the statement follows immediately from the standard age--lifetime relation:

```text
young age -> short available lifetime -> only high-mass stars can be dead
```

The precise ledger floors are:

- `33.9 Msun` on the most permissive 2 Myr formation-window branch;
- `45.4 Msun` as the minimum over coeval branches;
- `52.1 Msun` on the baseline branch;
- central baseline turnoffs of approximately `58` and `56 Msun` for A and B.

The global `33.9 Msun` floor is not measured from a vanished star. It is produced by the PARSEC, `R_V=3.0`, 2 Myr formation-window branch when the oldest sampled realization of subgroup B reaches about 5.619 Myr. The turnoff relation then gives `33.921 Msun`; the smallest sampled dead progenitor is `33.934 Msun` ([branch table](../tables/wp7_alpha_headline_branch_sets.csv)).

This is useful bookkeeping and a necessary boundary condition, but it is not a discovery. The manuscript should say that the high progenitor masses are an **expected consequence of the inferred young ages**, not that the Gaia analysis discovered a new very-massive death regime.

### 3.2 Stellar deaths are not identical to successful explosions

This physical distinction is well known. Sukhbold et al., Ertl et al., and a large explodability literature already establish that final fate is not a simple function of initial mass. The present project can quantify how an assumed explodability prescription changes the **Cygnus OB2** ledger, but it did not discover the distinction.

The current hard-cutoff branches are intentionally simple threshold experiments, not a detailed physical explodability map. Therefore:

- "the all-explode death count is conditional on explodability" is correct;
- "we discover that death and explosion differ" is not;
- "a hard cutoff at or below 30 Msun makes the ledger zero" describes the adopted model experiment, not a prediction that every star above 30 Msun fails in nature.

### 3.3 Age sensitivity

Age sensitivity is real and severe, but it is already shown by Menchiari et al. and is inherent to the turnoff method. The current common-age scan gives zero below about 2.75 Myr, `12.73` deaths at 4 Myr, and `38.20` at 6 Myr ([age scan](../tables/wp7_age_sensitivity.csv)). This is valuable disclosure of model dependence, not standalone novelty.

### 3.4 Monte Carlo sampling and a branch grid

Monte Carlo propagation, two evolutionary families, several IMF slopes, extinction-law branches, and star-formation-duration branches improve reliability and transparency. They are methodological quality unless they reveal a conclusion that a simpler model misses. A larger uncertainty table is not automatically a new astrophysical result.

### 3.5 Evidence that Cygnus OB2 already experienced a supernova

Wright et al. already made this case using the age/mass-function deficit, a pulsar, and a runaway O star. The present work may update the probability or test consistency under Gaia-defined membership, but it cannot claim first evidence for a past Cygnus OB2 supernova.

## 4. What the present analysis does provide now

The following are genuine outputs, although not all are yet demonstrated as material improvements:

- a Gaia DR3-selected, probabilistic stellar census divided into three kinematic subgroups;
- per-star extinction and stellar-mass inference rather than one regional correction;
- subgroup ages and separate IMF normalizations;
- a stochastic ledger of deaths and lookback times;
- explicit branch dependence rather than one assumed age/IMF/track;
- a conditional probability for a recent event;
- a reproducible chain from member catalogue to manuscript tables.

The baseline, under the permissive all-explode prescription, is:

| Quantity | Baseline result |
|---|---:|
| A age | 4.00 Myr |
| B age | 4.09 Myr |
| C age | 2.52 Myr |
| A deaths | 4.17 |
| B deaths | 4.26 |
| C deaths | 0.00 |
| Total deaths, labelled `N_SN` under all-explode | 8.43 |
| `P(last event < 100 kyr)` | 0.552 |
| Headline branch range | 5.63--28.7 |

These numbers are traceable to [the ledger](../tables/wp7_ledger.csv) and [generated manuscript numbers](../manuscript/numbers.tex). They are conditional predictions, not observed historical events.

## 5. The strongest candidate novelty: quantify aggregation bias

The scientifically interesting possibility is not that the population is young. It is that averaging a multi-age population can give the wrong death history because

```text
N_death(mean age) != sum over subgroups of N_death(subgroup age),
```

especially close to the first-death boundary.

The existing diagnostic is promising:

| Model | Mean deaths | `P(last < 100 kyr)` |
|---|---:|---:|
| Every subgroup forced to 4.0 Myr | 12.73 | 0.727 |
| Subgroup-resolved baseline | 8.43 | 0.552 |
| Difference | 4.30 (51% of resolved value) | 0.175 |

At a common 4 Myr, the expected subgroup contributions are approximately A `4.18`, B `3.97`, and C `4.57`. In the resolved baseline they are A `4.17`, B `4.26`, and C `0.00`. Almost the full difference is therefore caused by assigning the young C population the age of A/B.

This is **not yet a fair proof**. Four Myr was imposed rather than inferred by fitting the whole association as one population. A different single fitted age may reproduce the total count, even if it does not reproduce subgroup attribution or timing. In addition, C's contribution is branch dependent: it ranges from zero to about 7.3 deaths over retained all-explode branches. "C has produced no deaths" is therefore not robust; "C lies near a model-dependent first-death boundary and drives leverage" is the safer hypothesis.

## 6. Turning "potentially yes" into a rigid "yes"

### 6.1 Required model hierarchy

Run three matched analyses:

1. **`M_prior`: Menchiari-style model.** Fixed total mass, assumed common age, adopted IMF, and the same lifetime prescription used by Menchiari as closely as possible. First reproduce their 3 and 5 Myr results.
2. **`M0`: best simple model.** Analyse the same Gaia-selected stars as one population. Refit one common age and one global normalization; do not choose the common age to match the resolved answer.
3. **`M1`: subgroup model.** Use the current A/B/C memberships, ages, and normalizations.

Hold all other ingredients fixed between `M0` and `M1`: input catalogue, extinction machinery, isochrones, IMF family, upper-mass limit, formation-duration definition, multiplicity convention, and explodability prescription. Otherwise any difference cannot be attributed specifically to subgroup resolution.

### 6.2 Recommended pre-registered materiality tests

The exact tolerances should be finalized before examining the new `M0` result. A defensible starting set is:

| Question | Recommended pass criterion for a material improvement |
|---|---|
| Do the data support separate subgroup ages? | Held-out predictive improvement of `M1` over `M0`, for example `Delta ELPD > 2` standard errors, plus acceptable posterior predictive checks. |
| Does subgroup resolution materially change the total death count? | The paired 95% interval for `N_death(M1)-N_death(M0)` excludes a practical-equivalence interval of `+/-10%`, and the absolute change is at least about 3 deaths (comparable to the baseline stochastic half-width). |
| Does it change the recent-event inference? | `|Delta P(last < 100 kyr)| >= 0.10`, with uncertainty excluding a negligible difference. |
| Does it change the time history rather than only the integral? | A predeclared distance between the two `R_death(t)`/last-event distributions exceeds an injection-calibrated threshold, with a physically meaningful shift in an event-time quantile (candidate scale: 0.1 Myr). |
| Is subgroup attribution informative? | At least one subgroup changes by a predeclared material amount (candidate: 2 deaths) or crosses a predeclared contribution boundary, and this cannot be reproduced by the pooled model. |
| Is the conclusion robust to model branches? | The direction and material size of the effect survive in at least 80% of retained branches and under both PARSEC and MIST; exceptions are reported, not averaged away. |
| Is the extra complexity calibrated? | In mock/injection populations with known multi-age histories, `M1` reduces bias or restores nominal interval coverage relative to `M0`; in truly coeval injections it does not manufacture a false subgroup advantage. |

A rigid "yes" requires both **statistical evidence** that the extra structure is supported and **scientific materiality** in an output that matters. A tiny improvement in fit with unchanged death history is not enough; a large change produced by an unsupported overfit is also not enough.

### 6.3 Variance and leverage decomposition

To show how the work helps other researchers, decompose uncertainty in `N_death`, the last-event time, and recent-event probability into contributions from:

- subgroup resolution versus one population;
- age uncertainty;
- extinction law;
- isochrone family;
- IMF slope and normalization;
- formation duration;
- stochastic upper-IMF sampling;
- explodability;
- unmodelled binaries.

The practical result should be a rule such as:

> A one-age approximation is adequate when the intrinsic age spread is below `X` and every subgroup lies safely on the same side of the first-death boundary; it becomes biased by more than `Y%` when a young component containing at least fraction `f` of the IMF normalization is mixed with an older component.

That rule would be portable to other young associations and is more useful than one additional Cygnus-specific number.

### 6.4 External checks must remain independent

Use the pulsar, runaway stars, isotope limits, and high-energy event hypothesis only after the stellar-demographic models are fixed. They can test whether the predicted timing is plausible, but should not be used to tune the ages or select a branch and then be presented as validation.

## 7. Other possible contributions, with requirements

| Candidate contribution | Current status | What would make it strong |
|---|---|---|
| First Gaia DR3 subgroup-resolved Cygnus OB2 death-time distribution | No exact duplicate found in the 2026-07-30 sweep, but uniqueness alone is insufficient. | Demonstrate calibration, release the posterior/data product, and show a conclusion that differs materially from `M0`/Menchiari. |
| Identify CygOB2-C as the leverage point | Supported qualitatively; exact contribution is branch dependent (`0` to about `7.3`). | Show a robust variance/leverage decomposition and specify which new age/extinction observation would reduce the SN-history uncertainty most. |
| Last-event probability | Potentially a useful new Cygnus-OB2-specific product. | Demonstrate that subgroup resolution changes/calibrates it, report full conditionality, and compare independently to external markers without tuning. |
| Empirical test of the fixed-age estimator | Not yet performed. | Reproduce Menchiari under matched assumptions, then map the domain where the simple estimator is adequate or biased using controlled injections. |
| Death-to-explosion uncertainty for Cygnus OB2 | Application of known theory, not new physics. | Replace or supplement hard cutoffs with published explodability prescriptions and clearly separate data-constrained death history from theory-conditioned explosion history. |
| Reusable Gaia-to-feedback workflow | Strong technical product but not automatically a strong astrophysical paper. | Validate on injections and preferably one comparison association; document inputs/outputs and release a portable implementation. |
| Test of a recent gamma-ray-producing supernova | Conditional and not unique; several source scenarios exist in Cygnus. | Predefine the external event window/energy/location requirements, propagate demographic uncertainties, and phrase the result as compatibility/exclusion rather than confirmation. |

## 8. Manuscript consequences

### Claims to remove or demote

- Do not present "only very massive stars have died" as a discovery. Present it as the expected age--lifetime boundary.
- Do not claim discovery that stellar deaths and successful explosions differ. Cite the established explodability literature and call this a model dependency.
- Do not imply that Menchiari merely made an uncontrolled rough estimate while this work uniquely predicts a death count. Their Appendix B already performs a Cygnus-OB2-specific mock-population death/SN calculation.
- Do not call `33.9 Msun` an observed progenitor floor. It is the lowest modelled turnoff over the widest formation-duration branch.
- Do not call the all-explode ledger an observed supernova history. Prefer `stellar-death history` or `N_death`; write `N_SN | all explode` only when the conditioning is explicit.
- Do not make technical complexity, number of branches, or Monte Carlo size the novelty claim.

### Defensible wording now

> We infer a branch-resolved, single-star stellar-death history from a Gaia DR3 subgroup census. The resulting supernova count remains conditional on explodability and binary evolution. A comparison with fixed-age population estimates shows that age assumptions dominate the answer; whether subgroup resolution yields a material improvement over the best one-population fit is tested explicitly rather than assumed.

### Stronger wording only after the ablation passes

> A one-age treatment of Cygnus OB2 biases the inferred death count and recent-event probability because its young subgroup lies below the first-death boundary while the older subgroups lie above it. The subgroup-resolved model is predictively favoured, reduces calibrated reconstruction bias, and changes the inferred history by `[measured effect]` across `[robust branch fraction]`.

## 9. Publication decision tree

1. **If `M1` is predictively supported and materially changes the history:** the project has a clear astrophysical result suitable for a full A&A/MNRAS-style paper, subject to the remaining model limitations.
2. **If `M1` improves the stellar fit but the death history is equivalent to `M0`:** subgroup science may still be publishable, but the supernova-history novelty is weak. The paper should focus on the resolved stellar population and treat the SN ledger as an application.
3. **If `M0` and `M1` are equivalent in both fit and outputs:** the complex analysis validates a simpler estimator. That is useful as a methods/benchmark result only if the validity domain is mapped across realistic simulations or multiple associations; it is unlikely to support the current strong narrative by itself.
4. **If `M1` overfits or fails calibration:** do not use the subgroup-resolved SN history as the headline result. Preserve the negative finding and simplify the paper.

The work already done is not foolish even if a simpler count agrees. It has located the assumptions and made a decisive validation experiment possible. But the publication case must come from the result of that experiment, not from the amount of work invested.

## 10. Internal symposium talk

The requested title, **"Tracing the Massive-Star and Supernova History of Cygnus OB2,"** remains accurate and does not commit to a novelty claim. No abstract was submitted, so there is no scientific reason to withdraw the request. A defensible ongoing-work talk can be organized around:

1. what earlier fixed-mass/fixed-age calculations already established;
2. the Gaia DR3 subgroup and age reconstruction;
3. why the nonlinear first-death boundary creates possible aggregation bias;
4. the existing 4 Myr diagnostic;
5. the pending `M0` versus `M1` ablation and its predeclared pass/fail criteria;
6. what will be concluded if the simple model wins.

This is scientifically honest and likely more interesting than presenting a polished but overstated novelty claim.

## 11. How to avoid this problem in future projects

Before undertaking the expensive analysis:

1. **Write the closest-prior-work table first.** Identify the paper with the same target, the paper with the same method, and the paper with the same claimed result. A project is not safe merely because no single paper contains all three.
2. **Reproduce the strongest simple baseline.** Implement the prior-art calculation before the complex pipeline. If the proposed method cannot state in advance what it may change, stop and redesign.
3. **Write one falsifiable incremental-value sentence.** Use the form: "Relative to baseline `B`, method/data `X` changes scientifically relevant quantity `Y` by at least `T`, because `mechanism`."
4. **Pre-register a practical-equivalence region.** Decide how small a difference would mean that the simple model is adequate.
5. **Use staged stop/go gates.** Pilot membership -> age separation -> pooled-versus-resolved output -> full systematic grid. Do not build the final machinery before the central effect survives a cheap pilot.
6. **Separate categories of novelty.** New data, new method, new result, and new implication are different. "First combination" and "more rigorous" are weak unless they change inference or enable a reusable capability.
7. **Run a hostile-referee test early.** Ask: "Could the headline conclusion be obtained from one standard lifetime lookup and one IMF integral?" If yes, that conclusion cannot carry the paper.
8. **Repeat the literature search before submission.** The local sweep already states that its 2026-07-30 result must be refreshed immediately before posting.

## 12. Immediate work plan

- [ ] Reproduce Menchiari's 3 and 5 Myr counts under matched mass, IMF, lifetime, and upper-mass assumptions.
- [ ] Implement and pre-register the pooled one-population model `M0`.
- [ ] Freeze the materiality thresholds before reading the `M0` versus `M1` outcome.
- [ ] Compare predictive fit, total deaths, subgroup attribution, event-time history, and recent-event probability with paired uncertainties.
- [ ] Run coeval and multi-age injection tests to measure bias and coverage.
- [ ] Decompose which subgroup and which uncertainty source controls each headline output.
- [ ] Decide the manuscript narrative from the pass/fail result; do not decide it in advance.
- [ ] Redraw the sensitivity figure to show progenitor/death-mass distributions or branchwise mass floors if that result remains in the paper.
- [ ] Refresh the strict prior-art search immediately before submission.

## Bottom line

The current study should not claim novelty from the facts that a young cluster loses only its most massive stars, that death counts depend strongly on age, or that core collapse and a successful explosion are different. Those are established.

The viable scientific question is whether a Gaia-resolved multi-age population produces a materially different and better-calibrated death history than the best possible one-age population model. The existing numbers suggest that it might, mainly because the young C subgroup sits near the first-death boundary. The controlled pooled-versus-resolved analysis is what can turn that possibility into a rigid yes—or demonstrate honestly that the simpler calculation is sufficient.
