# Manuscript reframing and scientific-revision brief

## Purpose

Use this document as the complete brief for revising the Cygnus OB2 Paper 1
manuscript. The current manuscript is in `manuscript/main.tex`; its numerical
macros are generated in `manuscript/numbers.tex`. The project plan is
`paper1_execution_plan.md`, and the analysis history is `PROJECT_TRACE.md`.

The revision must preserve the central scientific achievement: a
branch-resolved reconstruction of the massive-star and supernova history of
Cygnus OB2 from a homogeneous Gaia DR3-based census. The paper should not be
reduced to a generic astrometry-method article. Its strongest contribution is
that it makes the supernova history measurable while showing, quantitatively,
which conclusions are controlled by the data and which depend on age, IMF,
explodability, extinction, membership, binaries, and runaway corrections.

The Cygnus cocoon application and the COSI forecast should remain, but as
applications of the stellar-population result rather than the foundation on
which the paper's importance rests.

## Non-negotiable working rules

1. Do not overwrite or silently replace any upstream WP artifact. Preserve the
   existing repair/version history and failed products.
2. Do not change a scientific threshold, gate, branch definition, or headline
   branch set merely to improve the conclusion.
3. Do not type new measured numbers directly into `manuscript/main.tex`. Add
   derived quantities to `scripts/wp10_numbers.py`, resolve their inputs through
   `scripts/wp10_inputs.py`, regenerate `manuscript/numbers.tex`, and rerun
   `scripts/wp10_validate.py`.
4. Clearly distinguish:
   - observed catalogue facts;
   - fitted/model-dependent quantities;
   - conditional scenario calculations;
   - formal probability statements;
   - bounds and sensitivity tests.
5. A failed pre-registered prediction or gate must remain labelled as failed.
   It may be explained, but not reclassified after seeing the result.
6. Use final `repair_v7` products unless a new, explicitly versioned analysis is
   required by this brief. Do not consume an unversioned or historical WP5
   product by accident.
7. Before changing the manuscript, inspect and preserve the current dirty
   worktree. Do not stage or commit unless separately authorized.
8. Every substantive literature claim must have a primary source. Do not add a
   scientific statement that cannot be sourced.

## Recommended scientific framing

### Central question

Replace the implicit question

> What is the probability that Cygnus OB2 supplied the supernova required by a
> particular gamma-ray cocoon model?

with the broader but still concrete question

> What supernova history is implied by the present-day massive-star census of
> Cygnus OB2, which parts of that reconstruction are robust to the carried model
> choices, and which physical uncertainties prevent a unique history?

### Central contribution

The paper should present:

1. A probabilistic Gaia DR3 membership and kinematic decomposition of Cygnus
   OB2, with spectroscopic anchors retained explicitly.
2. A spatially resolved extinction and age analysis coupled consistently to
   the mass-function inference.
3. A response-aware normalization of the 2--8 solar-mass census and an
   out-of-sample closure test above 8 solar masses.
4. A branch-resolved supernova ledger that never hides discrete model choices
   inside a single marginalized number.
5. A diagnosis of the dominant physical uncertainties: subgroup ages,
   high-mass IMF slope, very-massive-star explodability, subgroup-dependent
   closure, and incomplete 3D runaway information.
6. Independent checks from the pulsar, gamma Cygni, remnant visibility, and
   radioactive isotopes.
7. Two forward-looking tests: Gaia DR4 for kinematics/ages and COSI for the
   isotope-yield/IMF combination.

### Recommended one-sentence paper claim

> We reconstruct a branch-resolved supernova history for Cygnus OB2 from its
> Gaia DR3 massive-star census and show that the existence of past explosions is
> secure, whereas their number, recency, and relevance to the Cygnus cocoon are
> controlled primarily by the high-mass IMF, subgroup ages, and the uncertain
> fate of stars above roughly 30 solar masses.

### Claim hierarchy

The manuscript should explicitly organize conclusions into three levels.

#### Observationally or internally well-supported

- The selected member sample contains three stable kinematic subgroups.
- Within the parallax-selected sample, the data do not require two distance
  components.
- Differential extinction is spatially coherent and cannot be represented by a
  single regional value.
- The adopted age solution has two older subgroups and one younger subgroup,
  with a much wider systematic age envelope than the three central values alone
  imply.
- At least one successful explosion occurred, independently established by
  PSR J2032+4127 if its association is accepted.
- The model-branch spread in the inferred supernova budget exceeds the
  stochastic interval of any one branch.

#### Model-dependent but quantitatively constrained

- Baseline all-explode ledger: mean `N_SN = 8.43`, median 8, 68 percent interval
  `[5, 11]`.
- Headline all-explode range: `5.63--28.74`, depending principally on IMF slope.
- Baseline probability of a most recent event within 100 kyr: about 0.552.
- Association-wide closure is near a Salpeter-like slope, but subgroup closure
  is heterogeneous and Cyg OB2-C remains anomalous.
- All already-dead stars in the adopted single-star age/turnoff calculation are
  very massive; therefore the ledger is a probe of very-massive-star
  explodability rather than the canonical 15--25 solar-mass regime.

#### Conditional applications, not direct measurements

- Availability of an event satisfying the cocoon model's age/type/location
  requirements.
- The mapping from ZAMS mass to stripped-envelope SN type.
- The in-situ fraction inferred from an incomplete two-dimensional runaway
  census.
- COSI detectability of 60Fe on a chosen yield arm.
- The contribution of Cygnus OB2 rather than neighbouring Cygnus populations to
  an event somewhere in the cavity.

## Mandatory scientific corrections

### 1. Correct the membership-method description

Current problem:

- `manuscript/main.tex` says membership was determined in
  `(l, b, pmra, pmdec)` and that parallax was excluded.
- The final membership mixture in
  `provenance/wp2_membership_manifest.json` used
  `(relative_l, relative_b, parallax_corrected, pmra, pmdec)`.
- Only the subsequent subgroup clustering in
  `scripts/wp2_derive_subgroups.py` used `(l, b, pmra, pmdec)` without parallax.

Required fix:

- Rewrite the membership subsection so association selection and subgroup
  decomposition are described as two distinct operations.
- Correct the Figure 1 caption.
- Retain and strengthen the existing warning that the non-confirmation of the
  Berlanas two-distance structure is sample-conditional because membership used
  parallax.
- Use the exact claim: there is no evidence for two distances among the selected
  members. Do not claim that the foreground population does not exist.

### 2. Explain the WP2 recall gate honestly

The recorded gate passes at `189/229 = 82.53 percent`, but this consists of:

- 128 automatically recovered stars among 168 quality-analyzable literature
  members;
- 61 manually retained spectroscopic quality exceptions.

The automatic recovery fraction is `128/168 = 76.2 percent`. Manual treatment
was allowed and required by the plan, so this is not a protocol violation. It is
not, however, an independent 82.5 percent machine-classifier recall.

Required fix:

- State automatic and manual recovery separately.
- Explain that the control fields were used to calibrate the mixture prior and
  are therefore not a wholly independent validation set.
- Avoid describing the measured control yield as a fully external false-positive
  rate without that qualification.

### 3. Correct the WP5 residual-failure description

Final `repair_v7` recount:

- 40 of 54 subgroup fit cells pass; 14 fail.
- Failures by subgroup: A 4, B 7, C 3.
- Failures by family: PARSEC 4, MIST 10.
- Failures by R_V: 4 at 3.0, 5 at 3.1, 5 at 3.5.
- Failures by alpha: 7 at 2.0, 2 at 2.3, 5 at 2.6.

The statement that failures concentrate in the `R_V=3.5` and `alpha=2.6`
corners is false for the final table.

Required fix:

- Replace that statement with the actual distribution.
- Distinguish clearly between:
  - baseline branch accepted for all three subgroups;
  - 40/54 subgroup cells passing;
  - the policy decision to retain failed non-baseline branches as explicit
    sensitivity cases.
- Do not say or imply that the full grid passed the WP5 residual gate.

### 4. Add an all-subgroup-pass sensitivity analysis

Within the 36 headline branches:

- 27/36 subgroup cells pass individually.
- Only 5/12 family--R_V--alpha combinations pass in all three subgroups.
- With the three star-formation durations, 15/36 headline ledger branches are
  based entirely on passing subgroup fits.

The qualitative alpha split survives when only all-subgroup-pass combinations
are retained:

- three valid alpha=2.0 duration branches remain above the 0.5 threshold;
- twelve valid alpha=2.3 duration branches remain below it.

Required fix:

- Generate a small table and/or overplot showing the full headline set and the
  strict all-subgroup-pass subset.
- Present this as a robustness result: the alpha split survives, but the
  branch representation is imbalanced and the 18/18 versus 0/18 language uses
  branches that include failed WP5 cells.

### 5. Reframe census closure

Current association-wide results near alpha=2.3 are approximately:

- A: 0.95;
- B: 1.07;
- C: 1.39;
- association grid median: about 1.07;
- closing-alpha median: about 2.25.

Required fix:

- Replace the heading `The census closes at Salpeter` with wording such as
  `Association-wide closure is near Salpeter, but the subgroups differ`.
- Do not claim that a single power law is validated uniformly from 2 to 120
  solar masses.
- Present C's closure excess, shallower preferred slope, family-dependent age,
  and possible distance contamination as a connected scientific finding.
- State that a global alpha forces the same slope on all subgroups. Add or at
  least discuss a mixed/subgroup-dependent slope sensitivity, because A and B
  prefer approximately 2.3 while C tends toward approximately 2.0.

### 6. Report the WP7 formal gate correctly

Prediction L6 and gate G7a failed the pre-registered relative convergence
criterion. The stored execution record has `G7a_converged: false`.

At the same time, substantive convergence of material cells is good; the formal
failure was driven by an ill-specified relative criterion around small or zero
quantities.

Required fix:

- Say explicitly: the formal convergence gate failed as registered.
- Then report the absolute/material-cell drift demonstrating numerical
  stability.
- Do not rewrite the failed gate as a pass.

### 7. Rename the explodability branch

The implemented `islands` branch is not a Sukhbold/Ertl non-monotonic
explodability map. It is a hard cutoff: no SN above 25 solar masses.

Required fix:

- Rename it throughout to `extreme 25-solar-mass direct-collapse cutoff` or
  equivalent.
- Present the black-hole-threshold scan as the primary explodability result.
- Explain that the physical conclusion is conditional on the fate of stars
  above approximately 30--50 solar masses.
- Do not imply that a detailed high-mass explodability prescription was
  implemented.

### 8. Narrow the pulsar claim

PSR J2032+4127 establishes that at least one successful explosion occurred if
its association with Cygnus OB2 is accepted. It excludes a zero-supernova model.
It does not establish that one of the currently inferred >34--52 solar-mass
single-star progenitors exploded, because an older population or a
binary-stripped lower-mass progenitor remains possible.

Required fix:

- Replace `the pulsar settles that branch` in the abstract and conclusions.
- Use: `the pulsar independently excludes an entirely explosion-free history`.
- Retain the three alternative interpretations already acknowledged in the
  Results.

### 9. Treat C3 as a model mapping, not a measured probability

The code sets `C3=1` because every simulated progenitor exceeds a fixed
30-solar-mass stripping threshold. This does not demonstrate that every event
would observationally be Type Ib/c.

Required fix:

- Write `C3=1 under the adopted deterministic mass-to-type mapping`.
- Cite primary stellar/binary-evolution literature showing the dependence of
  envelope stripping on winds, metallicity, rotation and binarity.
- Preferably add alternative subtype mappings or a bounded sensitivity branch.
- If no new subtype model is added, do not call C3 a measured probability and do
  not state that all explosions `should have been` Type Ib/c without a clear
  conditional qualifier.

### 10. Replace P_verdict with a defensible quantity

`P_verdict = C1 * C3 * C4` is computed correctly from the declared formula, but
its interpretation is not defensible as a calibrated probability:

- C4=0.854 is an upper bound, not a sampled probability.
- The runaway sample is two-dimensional and footprint-limited.
- A living massive-star runaway fraction need not equal the escape fraction of
  already-dead, more massive progenitors.
- C3 is a deterministic model mapping.
- Independence of C1, C3 and C4 is asserted rather than demonstrated.
- The quantitative neighbouring-association contribution required by WP8.5 and
  WP9 is absent.

Required fix:

- Preferred: rename the result a `conditional scenario-availability score` or
  `optimistic upper-bound scenario probability`.
- Show the result as a function of C4 rather than fixing C4 at 0.854.
- Correct the threshold statement:
  - C4 below approximately 0.720 moves at least one alpha=2.0 branch below 0.5;
  - C4 below approximately 0.581 moves every alpha=2.0 branch below 0.5.
- Do not say `below 0.60 moves any branch`; that should say `all branches`.
- If the paper retains a formal probability interpretation, construct a
  probabilistic escape model and a probabilistic subtype model, propagate them
  jointly, and justify dependencies.

### 11. Complete or explicitly remove WP8.5 from the claimed probability

The original plan required a coarse literature-based SN budget for neighbouring
associations such as Cyg OB1 and Cyg OB9 to bound the probability that a cavity
event originated outside Cyg OB2. The WP8 report explicitly says no quantitative
budget was claimed.

Required choice:

- Either complete a sourced, explicitly coarse alternative-source budget and
  propagate it into the cocoon application;
- or state that the analysis estimates only the availability of an event from
  Cygnus OB2 and cannot quantify whether Cygnus OB2 was the source of an event
  elsewhere in the larger cavity.

Do not retain language suggesting that the originally planned combined
alternative-source probability was calculated when it was not.

### 12. Qualify the binary `bound`

The +/-30 percent binary calculation is a deliberate scenario bracket, not a
formal statistical bound on binary evolution at 4 Myr.

Required fix:

- Call it a conservative sensitivity bracket.
- Avoid claiming the unmodelled binary systematic has been bounded in a strict
  probabilistic sense.
- Preserve the important result that the tested bracket is smaller than the
  carried model-branch spread.

### 13. Correct WP11 provenance language

The WP11 report claims its preregistration had hashed inputs and was written and
committed before execution. In the current worktree:

- the preregistration and WP11 files are untracked;
- the preregistration contains no input hashes;
- its timestamp precedes execution by about 16 minutes.

Required fix:

- Do not claim it was committed or immutable before execution.
- Describe it accurately as a timestamped local specification written before
  the scoring execution, with prior qualitative expectations explicitly
  disclosed.
- Do not attempt to retroactively manufacture stronger preregistration status.
- Keep the manuscript's current disclosure that WP11 was selected after the
  ledger existed.

### 14. Keep the corrected 26Al result

WP8's old approximately 1-solar-mass denominator was wrong by about two orders
of magnitude. The frozen measured flux and distance imply approximately
`8.86e-3 solar masses` of 26Al, subject to the stated line-yield conversion.

The manuscript's current flux-space statement that SN ejecta provide roughly
10--42 percent of the measured complex-wide flux is the correct downstream
presentation.

Required fix:

- Preserve the WP8 historical error as a withdrawn result.
- Use the corrected flux-space comparison in the paper.
- Do not describe the SN component as the whole complex measurement; winds and
  neighbouring populations contribute.

### 15. Qualify the COSI forecast

The primary yield arm gives a clean split: 18/18 alpha=2.0 branches above the
adopted COSI sensitivity and 0/18 alpha=2.3 branches above it. However:

- the split is conditional on one yield arm;
- the between-yield-arm spread is about 30-fold;
- the within-primary-arm headline branch spread is about 4.8-fold;
- a standard high-mass collapse prescription gives zero isotope yield for this
  ledger.

Required fix:

- Lead with joint IMF/yield discrimination, not with an unconditional promise
  that COSI will determine alpha.
- Preferred wording: `A COSI measurement or upper limit would constrain the
  combination of high-mass IMF and very-massive-star yield/explodability models`.
- State the 18/18 versus 0/18 result immediately afterward as the prediction on
  the primary yield arm.
- Retain the verified 3-sigma, two-year narrow-line sensitivity of
  `3e-6 ph cm^-2 s^-1` for the 60Fe lines, citing the primary COSI mission paper.

## Recommended manuscript structure

### Title and subtitle

Keep the current main title or use one of these:

1. `The supernova history of Cygnus OB2 from Gaia DR3`
   - Subtitle: `A branch-resolved census of massive-star deaths`
2. `Reconstructing the supernova history of Cygnus OB2`
   - Subtitle: `What Gaia measures and what stellar physics still controls`
3. `The massive-star and supernova history of Cygnus OB2`
   - Subtitle: `A probabilistic census with explicit model branches`

Avoid a title led by the cocoon probability or COSI prediction.

### Abstract

Rebuild the abstract around this order:

1. Context: supernova histories of young associations are needed for feedback,
   remnants, isotopes and high-energy interpretations, but are difficult because
   they combine membership, extinction, age, completeness and high-mass stellar
   evolution.
2. Aim: reconstruct the Cygnus OB2 history and identify which conclusions are
   robust versus branch-dependent.
3. Method: Gaia+2MASS+spectroscopic anchors, two-stage membership/subgroups,
   spatial extinction, age/mass inference, response-aware IMF normalization,
   closure test, stochastic ledger, explicit branches.
4. Results: three subgroups; baseline and branch-range SN counts; recent-event
   probability as a branch-dependent quantity; high-mass explodability
   conditional; heterogeneous closure.
5. Validation: pulsar excludes an explosion-free history; other markers are
   consistent but not individually decisive.
6. Application: cocoon availability is conditional and IMF-sensitive, not a
   stand-alone calibrated probability.
7. Forecast: COSI can constrain the joint IMF/yield problem; Gaia DR4 improves
   ages and 3D runaways but not the high-mass IMF by itself.

Remove or rewrite:

- `the pulsar settles that branch`;
- unconditional `P = ...` language for the cocoon scenario;
- any implication that every marker, including WP11, was frozen before the
  ledger existed.

### Introduction

Expand the motivation beyond the cocoon while retaining it as a concrete use
case.

Suggested sequence:

1. Why young-association supernova histories matter: feedback, superbubbles,
   compact remnants, radionuclides and cosmic-ray acceleration.
2. Why the measurement is hard: present-day census versus initial population;
   extinction, membership, incompleteness, stochastic sampling, age and
   explodability.
3. Why Cygnus OB2 is the appropriate test case: rich massive-star population,
   differential extinction, substructure, young pulsar, gamma-ray environment.
4. Limitations of previous population-synthesis estimates: heterogeneous census
   and assumed global age/mass.
5. This paper's methodological advance: the whole chain is internally linked,
   branch-resolved and provenance-controlled.
6. The cocoon and isotope applications as demonstrations of what the ledger can
   and cannot decide.

### Data section

Add compact tables for:

- Gaia narrow and wide samples;
- 2MASS match counts;
- spectroscopic anchor sources and evidence rows;
- quality-analyzable versus quality-exempt anchors;
- external marker catalogues and whether they were frozen at WP1 or added
  post-hoc.

Explicitly distinguish source identities from literature-evidence rows.

### Methods section: expand substantially

#### Membership selection

Include:

- selection box and Gaia quality criteria;
- parallax zero-point treatment;
- DBSCAN's limited role as a training seed;
- final cluster-versus-field mixture and all five features;
- full Gaia covariance Monte Carlo;
- control-prior calibration;
- soft probabilities and manual anchor exceptions;
- separate automatic/manual recall numbers;
- limitations from selection circularity.

#### Kinematic subgroups

Include:

- clean automatic-member sample;
- four-dimensional `(l,b,pmra,pmdec)` feature space;
- deliberate exclusion of parallax here, not at membership stage;
- seed-stability criterion and null/bootstrap comparison;
- subgroup-label uncertainty and unassigned anchors.

#### Distance test

Include:

- one- versus two-component comparison;
- depth estimate and delta BIC;
- exact scope of the non-confirmation;
- why a parallax-blind membership rerun is deferred.

#### Extinction

Include:

- anchor derivation from fixed-temperature intrinsic colours;
- six-band photometric treatment and missing-band masking;
- variogram/simple-kriging prior mean and uncertainty;
- why the earlier nearest-anchor prior was repaired;
- B's sparse local anchor support;
- absolute approximately 0.5-mag anchor/photometry scale discrepancy;
- downstream effect on ages and normalization.

#### Ages and masses

Include:

- PARSEC and MIST branches;
- upper-MS and counts-based age indicators;
- measurable/excluded/railed branch counts;
- native isochrone-node treatment;
- age--normalization covariance and joint age-k inference;
- B's posterior railing and one-sided lower-bound consequence for N_SN;
- C's family-dependent boundary at the 120-solar-mass IMF ceiling.

#### Completeness and IMF normalization

Include:

- end-to-end catalogue-level injection chain;
- no absolute 95-percent completeness edge;
- response-corrected Poisson likelihood rather than a simple completeness cut;
- 2--8 solar-mass calibration window;
- residual gate definitions;
- truth-age marginalization repair;
- mass-dependent binary-fraction repair_v7;
- baseline pass versus partial grid failure;
- why failed branches are retained as sensitivity cases.

#### Massive-star closure

Include:

- observed and predicted sides of the >8-solar-mass comparison;
- response integration below the observed threshold to capture up-scatter;
- members, orphan anchors and false-positive-corrected runaways;
- association and subgroup closure values;
- distinction between an aggregate near-Salpeter result and subgroup
  heterogeneity;
- limits of the multiplicity test.

#### Runaways

Include:

- peculiar rather than absolute proper motions;
- control-field false-positive correction;
- two-dimensional traceback;
- finite wide-box velocity/time ceiling;
- why the recovered runaway count is a lower bound;
- why it constrains location rather than the number of deaths;
- why applying the living runaway fraction to dead progenitors is an assumption.

#### Supernova ledger

Include:

- paired age-k posterior draws;
- stochastic Poisson population sampling;
- IMF truncation and stellar lifetime inversion;
- star-formation-duration convention;
- definition of death versus successful explosion;
- all-explode branch and hard direct-collapse threshold scan;
- baseline branch definition;
- 2 million iterations and the formal convergence-gate failure;
- separation of stochastic intervals from discrete branch spread.

#### External validation and applications

Separate three categories:

1. Frozen independent markers: pulsar, gamma Cygni, SNR catalogues, original
   26Al measurement.
2. Coarse environmental context: neighbouring populations, currently lacking a
   quantitative budget.
3. Post-hoc forward prediction: WP11 isotope forecast, selected after the ledger
   but specified before scoring.

Do not mix these evidential categories.

## Recommended Results structure

1. `Membership, contamination and kinematic substructure`
2. `Differential extinction and subgroup age structure`
3. `Response-aware mass functions and fit-quality map`
4. `The living massive-star census and subgroup-dependent closure`
5. `A branch-resolved supernova history`
6. `What is robust about the explosion history`
7. `Independent evidence for past explosions`

The Results should culminate in the ledger and its uncertainty structure. The
cocoon evaluation should remain in the Discussion.

## Recommended Discussion structure

1. `A measured history, not a single supernova number`
   - baseline versus branch range;
   - statistical, stochastic and model uncertainties separated.
2. `The very-massive-star explodability problem`
   - threshold scan;
   - pulsar evidence;
   - limitations of hard-cutoff prescriptions.
3. `Does one IMF describe all three subgroups?`
   - A/B versus C closure and slope preferences;
   - distance-contamination hypothesis;
   - mixed-slope sensitivity.
4. `Application to the Cygnus cocoon`
   - conditional scenario score;
   - C1/C3/C4 definitions and limitations;
   - C4 sensitivity;
   - missing neighbouring-association probability;
   - preferred versus permissive age windows.
5. `Radioactive-isotope consistency and forecast`
   - corrected 26Al flux comparison;
   - 60Fe forecast;
   - yield-arm dominance;
   - COSI as joint IMF/yield constraint.
6. `What Gaia DR4 will and will not resolve`
   - 3D runaway traceback;
   - subgroup ages and distance contamination;
   - high-mass IMF still limited by census size/modeling.
7. `Limitations and reproducibility`
   - environment/data release;
   - formal failed gates;
   - deferred parallax-blind membership and neighbour budget.

## Appendices to add or expand

### Appendix A: Full branch and gate table

The current Appendix A is empty and must be populated. Generate, do not type, a
table containing at least:

- family;
- R_V;
- alpha;
- star-formation duration;
- subgroup WP5 residual-gate status;
- all-subgroup-pass flag;
- subgroup and association closure values;
- N_SN central value and interval;
- recent-event probability;
- minimum progenitor/turnoff mass;
- conditional scenario score;
- notes for failed or excluded branches.

If too large for print, place a compact table in the paper and publish the full
machine-readable table electronically.

### Appendix B: Membership validation and distance-selection limitation

Add:

- automatic versus manual recall accounting;
- per-control-field yields;
- prior calibration procedure;
- subgroup seed-stability summary;
- one/two-distance comparison;
- an explicit diagram or explanation of why the test is conditional on the
  parallax-selected sample.

### Appendix C: Extinction, age and mass-function diagnostics

Add:

- anchor distribution by subgroup;
- variogram/kriging summary;
- anchor-versus-photometry offset;
- excluded and grid-railed age branches;
- B's posterior boundary;
- complete WP5 residual map and gate definitions;
- historical repair summary only where needed to understand the final method.

### Appendix D: Ledger algorithm and convergence

Add:

- compact pseudocode;
- branch definitions;
- lifetime/turnoff consistency;
- convergence-doubling results;
- explicit statement that L6/G7a failed formally;
- material-cell absolute drift showing practical stability;
- threshold-scan definition.

### Appendix E: External checks and conditional cocoon calculation

Add:

- preregistered prediction table with PASS/FAIL unchanged;
- pulsar assumptions and alternative interpretations;
- gamma Cygni and remnant-visibility calculations;
- definitions and status of C1, C3, C4;
- C4 sensitivity table;
- neighbouring-association limitation;
- distinction between frozen WP8 evidence and post-hoc WP11 forecast.

### Appendix F: Provenance and reproducibility

Expand the current short paragraph into a usable reproducibility map:

- versioned final input list;
- script-to-output map;
- gate records;
- superseded artifacts and why they are forbidden;
- exact Conda environment export;
- data acquisition/retrieval instructions;
- hashes and archival location/DOI when available;
- command sequence for numbers, figures, validator and LaTeX build;
- current known limitations.

Do not claim every historical manifest hash matches current files. Final WP10
authorized inputs match, while several older WP2--WP4 manifests refer to
superseded states.

### Appendix G: Withdrawn and failed results

Keep and expand the current idea. Include:

- absolute-versus-peculiar proper-motion runaway result;
- closure lower-integration-bound defect;
- earlier 26Al denominator error;
- failed WP7 L1 and L6 predictions;
- failed WP6/WP11 predictions as recorded;
- what changed downstream and what was preserved.

## Figures to add or redesign

Figures should explain the inference chain, not only display its last numbers.
Use readable labels at final A&A column width and avoid multi-panel plots whose
text becomes illegible.

### Figure 1: Analysis flow and evidence hierarchy

New figure. A compact flow diagram:

`Gaia+2MASS+spectroscopy -> membership -> subgroups -> extinction -> ages/masses
-> IMF normalization -> massive-star closure -> SN ledger -> external checks ->
conditional applications`.

Mark:

- observed inputs;
- inferred quantities;
- branch axes;
- independent validation inputs;
- post-hoc forecast;
- where failed gates or one-sided bounds enter.

This will make the paper's methodological contribution immediately visible.

### Figure 2: Membership and subgroups

Retain the current sky/proper-motion concept but revise it to show clearly:

- all candidate members by probability;
- the three subgroup labels;
- manual quality-exempt anchors with a distinct symbol;
- literature samples;
- which dimensions enter membership versus subgroup clustering.

Do not write in the caption that parallax was excluded from membership.

### Figure 3: Control fields and membership calibration

Redesign the current control-field figure. The present panels are too small.

Preferred design:

- one readable comparison of target and three controls using identical axes;
- a bar or dot plot of P>0.5 yields with area normalization;
- automatic/manual literature recovery alongside the controls;
- note that controls calibrated the mixture prior.

Move detailed sky panels to an appendix if necessary.

### Figure 4: De-reddened CMD/HRD and age evidence

Mandatory new figure. This was in the original figure plan and is currently
missing despite age being load-bearing.

Suggested panels:

- extinction-corrected CMD or HRD for A, B and C;
- PARSEC and MIST isochrones at adopted ages;
- spectroscopic anchors highlighted;
- stars used in upper-MS/counts fits versus excluded points;
- age posterior per subgroup and family;
- B's top-node railing visible;
- the full 2.25--5.67 Myr systematic envelope.

The figure must distinguish observational scatter from model-family spread.

### Figure 5: Extinction map and calibration support

New or appendix figure:

- spatial A_V map;
- anchor locations;
- local kriging support/uncertainty;
- B's sparse-anchor region;
- comparison of anchor and prior-free photometric A_V;
- optional variogram panel.

This figure explains why the repair changed B and why A/C stayed stable.

### Figure 6: Mass-function fits and residual-gate map

Expand the current mass-function figure:

- observed response-weighted mass functions and fitted models by subgroup;
- baseline residuals;
- a 3 x 18 or equivalent heat map of all 54 residual-gate verdicts;
- visual distinction between baseline pass, passing sensitivity cells, and
  retained failed cells;
- mark the strict all-subgroup-pass headline subset.

This directly prevents readers from interpreting 40/54 as a fully validated
grid.

### Figure 7: Massive-star closure by subgroup and alpha

New figure:

- predicted/observed >8-solar-mass ratio versus alpha;
- separate A, B, C and association curves;
- unity line;
- closing alpha for each subgroup;
- PARSEC/MIST and R_V variation as bands or points;
- C's excess made visually obvious.

This may be the paper's most scientifically interesting diagnostic figure.

### Figure 8: Branch-resolved supernova history

Retain the R_SN(t) concept, but consider:

- baseline subgroup stack;
- band or side panel showing the headline branch envelope;
- marker windows for the pulsar, 100 kyr question and gamma Cygni;
- clear separation between baseline stochastic uncertainty and model-branch
  spread.

### Figure 9: Age and explodability sensitivity

Combine or coordinate:

- N_SN and recent-event probability versus assumed age;
- N_SN versus direct-collapse threshold;
- IMF upper-limit/turnoff boundary for C;
- mark the adopted subgroup ages and their allowed envelope.

Rename all `islands` labels to the actual hard-cutoff prescription.

### Figure 10: Conditional cocoon score

Replace the current simple branch verdict figure with a more honest diagnostic:

- scenario score versus C4 for every headline branch;
- alpha=2.0 and 2.3 separated;
- current C4 upper-bound value marked;
- C4=0.720 and 0.581 threshold points marked;
- all-subgroup-WP5-pass branches emphasized;
- failed-WP5 branches de-emphasized or open-symbol;
- optional panels for preferred and permissive age windows.

The y-axis should not be labelled as an unconditional probability unless a new
probabilistic C3/C4 model is implemented.

### Figure 11: 26Al/60Fe forecast and yield uncertainty

New figure or appendix figure:

- predicted 60Fe flux by alpha and yield arm;
- COSI sensitivity line;
- current SPI upper limit;
- predicted 26Al SN fraction of the measured complex flux;
- visual comparison of between-yield-arm and within-census-branch spreads.

The figure should show immediately that yield uncertainty dominates, while the
clean alpha split exists within the primary arm.

## Tables to add

1. Data and validation inventory: rows, matches, anchors, controls and frozen
   external markers.
2. Subgroup summary: member count, A_V, age by indicator/family, turnoff,
   normalization and closure.
3. Baseline ledger plus headline range: N_SN, recent-event probability, first
   and last-event summaries, minimum progenitor mass.
4. Gate summary: every WP gate with PASS/FAIL/PARTIAL and the exact consequence.
5. Conditional cocoon terms: definition, numerical treatment, whether it is a
   probability/assumption/bound, and dominant limitation.
6. Isotope forecast: yield arms, branch range, instrument comparison and status.

## Conclusions: recommended structure

The conclusions should contain concise, scoped statements in this order:

1. Three stable kinematic subgroups are found within the selected Cygnus OB2
   member sample; the distance result is conditional on parallax-based
   membership.
2. Differential extinction and subgroup age structure materially change the
   high-mass history; B remains a one-sided age/SN lower bound.
3. The baseline supernova ledger is approximately eight events, but the honest
   result is the branch range, not the baseline alone.
4. All inferred deaths occupy the uncertain very-high-mass explodability
   regime.
5. The pulsar independently excludes an entirely explosion-free history.
6. Association-wide closure is near Salpeter-like, but subgroup heterogeneity,
   especially C, prevents claiming a single universal slope without
   qualification.
7. The cocoon application is conditional and dominated by IMF/type/location
   assumptions; it is not currently a fully calibrated probability.
8. Gaia DR4 will improve ages, distance contamination and 3D runaways; it will
   not by itself determine the high-mass IMF.
9. COSI can constrain the combined IMF/yield/explodability problem, with a clean
   alpha split on one specified yield arm.

## Reproducibility work required before submission

1. Populate Appendix A; the current empty heading means WP10's original
   external-reader gate is not met.
2. Add a data/code availability statement.
3. Export and pin the `cygob2-gaia` environment, including non-Python tools used
   to compile the manuscript.
4. Provide retrieval instructions and an archival location for the ignored
   approximately 8.4 GB data tree, with hashes.
5. Update `manuscript/README.md` and the header of `main.tex`: `aa.cls` is now
   present and the manuscript has been compiled.
6. Update `PROJECT_TRACE.md` so its title/status reflects WP10 and WP11
   accurately.
7. Make clear that the original plan contains WP0--WP10. WP11 is an additional
   post-hoc extension; only Part B was executed.
8. Label old manifests as historical/superseded. The final WP10 authorized
   inputs currently match, but some older WP2--WP4 hashes no longer describe
   current repaired files.
9. Ensure WP8/WP9/WP11 records carry complete input hashes in future executions.
   Do not claim hashes existed before a run when they did not.
10. Decide whether notebooks are a release requirement. If they are, create
    dedicated, commented notebooks or literate reports for WP4 and WP6--WP11;
    do not claim the current project has one executed notebook per step.
11. Re-run the WP0 literature/deduplication sweep immediately before submission.
12. Compile from a clean release checkout and archive the PDF build command,
    tool versions, log, and resulting SHA-256.

## Required validation after revision

At minimum:

1. Rerun `scripts/wp10_inputs.py`.
2. Rerun `scripts/wp10_numbers.py`.
3. Regenerate all affected figures.
4. Rerun `scripts/wp10_validate.py`; all seven checks must pass.
5. Compile with `latexmk -pdf main.tex` from `manuscript/`.
6. Confirm no undefined references/citations, no missing figures, and no
   materially overfull boxes.
7. Render every PDF page and visually inspect final-column label sizes,
   multi-panel readability, whitespace, appendix tables and figure captions.
8. Re-hash every authorized input and output.
9. Verify every number in the abstract, conclusions, tables and figure captions
   against its generated macro/source artifact.
10. Produce a final gate table distinguishing:
    - formal gate status;
    - scientific interpretation;
    - environment/build blockers;
    - remaining deferred work.

## Scope of necessary reanalysis

Do not restart WP1--WP7 wholesale. The core numerical chain is largely intact.
The minimum new quantitative work is concentrated in:

1. WP5 all-subgroup-pass sensitivity extraction.
2. WP6 subgroup-closure presentation and optional mixed-slope sensitivity.
3. WP9 C4 sensitivity and renaming/reinterpretation of P_verdict.
4. A probabilistic or explicitly conditional C3 treatment.
5. WP8.5 neighbouring-association budget, if the paper wants a quantitative
   cavity-source probability.
6. Figure/table generation and reproducibility packaging.

Any new analysis must be pre-specified, versioned, and must consume the frozen
accepted chain read-only unless its purpose explicitly requires a new branch.

## Definition of done for the revising agent

The revision is complete only when:

- the manuscript describes the implemented membership and subgroup methods
  correctly;
- baseline acceptance is not confused with full-grid WP5 validity;
- failed formal gates remain visible;
- the cocoon quantity is no longer misrepresented as an unconditional measured
  probability;
- C3 and C4 are labelled according to their actual epistemic status;
- the neighbouring-association limitation is resolved or explicitly scoped out;
- the age evidence, closure heterogeneity, residual-gate landscape and C4
  sensitivity have readable figures;
- Appendix A contains the branch/gate table;
- provenance and reproducibility instructions are sufficient for an external
  reader to trace the final numbers;
- WP10 validation and a clean LaTeX build pass with recorded artifacts;
- the final handoff lists every changed file, new analysis, validation command,
  PASS/FAIL result and remaining scientific limitation.

