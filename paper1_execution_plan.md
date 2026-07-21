# Paper 1 Execution Plan — "The Supernova History of Cygnus OB2 from Gaia DR3"

*Version 1.0, July 2026. Implementation-free master plan: principal steps only, each work package (WP) self-contained and delegable to a separate agent/session. Code design decisions are deliberately left open. Companion documents: `cygob2_zero_to_hero_action_plan_v3.md` (context and collaboration strategy) and `method_explained.md` (physics and rationale of every step — read the relevant section of that document before executing a WP here).*

---

## 0. Mission statement

**Objective.** Produce a probabilistic supernova (SN) history of the Cygnus OB2 association from Gaia DR3 + literature spectroscopy: the number of past core-collapse SNe, their timeline, and the probability that a SN occurred within the last ~50–100 kyr — the assumption underlying the Härer et al. 2025 (A&A, arXiv:2508.21644) interpretation of the Cygnus PeV gamma-ray bubble. Validate the ledger against independent SN markers (PSR J2032+4127, γ Cygni SNR, INTEGRAL 26Al).

**Primary deliverable.** A submittable manuscript (A&A Letters if the verdict is sharp; regular A&A article if constraints are soft) + a reproducible pipeline that will be rerun on Gaia DR4 (release 2026-12-02) for Paper 2.

**Explicitly OUT of scope for Paper 1:** wind luminosity L_w (Paper 2), gamma-ray data analysis (Paper 2), any other cluster (Paper 3), DR4 data (rerun later).

**Headline result formats (decide at WP9, not before):**
- Sharp: "P(≥1 SN in Cyg OB2 within the last 100 kyr) = X% [interval]" + verdict on Härer assumption.
- Soft: "N_SN(age) constraints; the last-SN posterior; consistency requirements imposed by the pulsar."
Both are publishable; never promise the sharp one externally before WP9.

---

## 1. Global conventions (binding for every WP)

1. **Data releases:** Gaia DR3 (`gaiadr3.gaia_source` + official cross-match tables). No DR2 inputs except for validation comparisons. DR4 rerun is a separate future task.
2. **Parallax zero-point:** always apply the Lindegren et al. (2021) correction before any distance use. Record corrected and raw values.
3. **Errors philosophy:** every derived quantity is a distribution, not a number. Three error classes, kept separate and labeled: (a) statistical (Monte Carlo over measurement errors and membership probabilities), (b) stochastic (random IMF sampling), (c) model branches (discrete choices — never averaged, always reported per branch).
4. **Mandatory model branches (minimum set):** isochrones {PARSEC, MIST}; IMF high-mass slope {2.0, 2.3, 2.6}; extinction law R_V {3.1 baseline; 3.0, 3.5 checks}; explodability {all m>8 explode; Sukhbold-like islands of implosion}; star-formation duration per subgroup {0, 1, 2 Myr}.
5. **Provenance log:** every WP maintains a running log of every cut, threshold, catalog version, and query text, with dates. This log becomes the reproducibility appendix. No exceptions — the paper's core claim is homogeneity and auditability.
6. **Artifact naming:** each WP outputs files named `wpN_<name>` with a short README block (content, provenance, upstream artifacts used). Downstream WPs consume only these named artifacts, never ad-hoc intermediate data.
7. **Validation gates are blocking:** a WP is "done" only when its gate criteria are met or the deviation is documented and explained. Unexplained disagreement never propagates downstream.
8. **Definition of the target:** "Cygnus OB2" = the association's subgroups as recovered in WP2, including the Berlanas two-distance-population structure. The wider Cygnus X populations (Cyg OB9, OB1, etc.) are handled only in WP8/WP9 as alternative SN sources, not censused in full.

---

## WP0 — Requirements extraction and dedup sweep

**Objective:** turn the target papers into a requirements table and confirm nobody has already published this analysis.

**Inputs:** Härer+25 (arXiv:2508.21644); Menchiari+24 (arXiv:2306.00946); Knödlseder+02 (astro-ph/0206045); Martin+10 (arXiv:1001.1522); Wright+15 (Cyg OB2 census); Berlanas+19 (MNRAS 484, 1838) and Berlanas+20 (A&A, arXiv:2008.09917); Fuchs+06 (MNRAS 373, 993); Zucker+22 (arXiv:2201.05124, Methods/Supplement); ADS.

**Procedure:**
1. From Härer+25, extract the exact SN requirement: explosion energy, progenitor type, age window, assumed location (cavity vs association core), and any plausibility argument they attach. Record verbatim quotes with section numbers. This defines the statement WP9 must evaluate.
2. From Menchiari+24 and Härer+25, tabulate every stellar-side number they assume (distance, association mass, age, O-star count) with quoted uncertainties → the "assumptions table."
3. From Knödlseder+02 and Martin+10, extract their predicted SN rates / recent-SN counts for Cygnus and the census inputs they used → the prior-art baseline your result will be compared against.
4. From Fuchs+06 and Zucker+22 Methods, write a one-page summary of the IMF-deficit implementation choices (mass ranges, IMF forms, corrections applied) → the methodological template.
5. Dedup sweep: ADS search for 2023–2026 papers citing (a) Martin+10, (b) Berlanas+20, (c) Härer+25; scan titles/abstracts for any Gaia-based SN-history or demographic test of the Cygnus SN hypothesis. Also check arXiv listings for the last 6 months (astro-ph.HE + astro-ph.GA, keyword "Cygnus").
6. Decision point: if a direct duplicate exists → escalate to user (options: differentiate, join, retarget to Carina/Wd1). If not → proceed.

**Outputs:** `wp0_requirements_table.md`, `wp0_prior_art_summary.md`, `wp0_dedup_report.md`.
**Gate:** requirements table complete with verbatim quotes; dedup verdict explicit.
**Caveats:** conference proceedings and theses can hide in-progress duplicates — check recent ICRC/TeVPA proceedings too. Strength: this WP costs days and de-risks months.

---

## WP1 — Data acquisition

**Objective:** assemble every input dataset, versioned and frozen.

**Inputs/sources:**
1. **Gaia DR3 region query:** sky box l ∈ [77°, 83°], b ∈ [−1.5°, +4°] (generous; covers Cyg OB2 + margins for runaway traceback in WP6 use a wider box: l ∈ [72°, 88°], b ∈ [−5°, +8°] as a second query). Parallax ∈ [0.35, 1.1] mas (wide on purpose), G < 19. Retrieve: astrometry + errors + correlations, photometry (G, BP, RP), RUWE, radial velocities where present, `phot_bp_rp_excess_factor`, astrometric quality flags.
2. **2MASS J/H/Ks** via the official Gaia archive pre-computed cross-match (`tmass_psc_xsc_best_neighbour` + join table). Record match rates.
3. **Spectroscopic anchors:** Wright+15 massive-star census (spectral types, ~169 OB stars); Berlanas+19/20 member lists and O-star spectroscopy (Teff, log g where given); the Galactic O-Star Catalog (GOSC, Maíz Apellániz) entries in the field; Galactic Wolf-Rayet catalog (Crowther, v. latest — 3 WRs attributed to Cyg OB2, record which).
4. **SN markers:** ATNF pulsar catalogue entries within the wide box (record P, P-dot, characteristic age, proper motion if any, distance estimates; PSR J2032+4127 is the anchor); Green's SNR catalogue + SNRcat entries in the region (γ Cygni G78.2+2.1 parameters: age estimates, distance estimates, references); INTEGRAL 26Al measurements of the Cygnus region (from Martin+09/10 papers — values with errors).
5. **Extinction reference:** a 3D extinction map covering 1–2 kpc in Cygnus (e.g., Green/Bayestar or Lallement/Vergely versions; pick one as baseline, one as check) — used for consistency checks of per-star extinctions, not as substitute.

**Outputs:** `wp1_gaia_narrow.cat`, `wp1_gaia_wide.cat`, `wp1_2mass_join.cat`, `wp1_spectroscopic_anchors.cat` (one merged table: source_id ↔ spectral type ↔ Teff/logg ↔ provenance), `wp1_sn_markers.md`, `wp1_extinction_refs.md`, plus the query texts in the provenance log.
**Gate:** row counts sane (narrow box: order 10^5–10^6); ≥90% of Wright+15 stars matched to Gaia source_ids (unmatched ones individually listed with reasons — saturation is expected for the very brightest).
**Caveats:** the brightest O stars may have poor or missing Gaia astrometry (saturation) — they must not silently drop out; carry them via the spectroscopic anchor table. Weakness: any error here poisons everything; hence the freeze-and-version discipline.

---

## WP2 — Membership and substructure

**Objective:** probabilistic member list of Cyg OB2 and its subgroups.

**Procedure:**
1. Quality filtering on the narrow catalog (RUWE, parallax S/N, BP/RP excess factor for photometric reliability). Log the fraction removed; do NOT apply these filters to the spectroscopic-anchor stars (they get manual treatment).
2. Zero-point-correct parallaxes.
3. Cluster in (l, b, parallax, pmRA, pmDec) space with a density-based algorithm; scan hyperparameters for stability (a result that appears only for one narrow hyperparameter choice is not a result).
4. Monte Carlo membership: resample astrometry from per-star covariance N times; membership probability = fraction of runs in cluster. Keep P > 0.05 stars with their probabilities (soft threshold; hard cuts only in sensitivity checks).
5. Substructure: identify subgroups; explicitly test the Berlanas two-population (~1.35/1.6 kpc) hypothesis — fit distance distribution of members with 1 vs 2 components, report preference.
6. Assign each spectroscopic-anchor star to a subgroup (or "unassigned" with reason).

**Outputs:** `wp2_members.cat` (source_id, membership P, subgroup label, corrected parallax), `wp2_subgroups.md` (per-subgroup: N, distance posterior, sky footprint), diagnostic plots (sky, PM plane, parallax histograms).
**Gate (blocking):** recover ≥80% of Berlanas+19 spectroscopic members with P > 0.5; every missed one individually explained (bad RUWE, saturation, disputed membership...). Compare subgroup structure qualitatively with Berlanas+19 and Wright+15 maps.
**Caveats/weaknesses:** Cygnus sits along a spiral arm — contamination by unrelated young stars at similar distance and PM is the dominant risk; quantify it via the PM-offset control fields (run the same pipeline on 2–3 nearby control boxes at same |b| — the "cluster" yield there estimates false-positive rate). The association is loose: expect fragmentation; do not force one blob. Strength: this is the best-validated step — two independent published lists to check against.

---

## WP3 — Per-star extinction and de-reddened CMD

**Objective:** extinction A_G (or A_V) per member star; extinction-corrected CMDs per subgroup.

**Procedure:**
1. For each member: fit Gaia G/BP/RP + 2MASS J/H/Ks photometry against reddened synthetic photometry (from the same isochrone families used in WP4 — consistency requirement) with extinction as free parameter. Baseline R_V = 3.1.
2. Where 2MASS is missing, flag the star; optical-only extinctions get inflated errors.
3. Sanity check the resulting extinction map: (a) spatial coherence (extinction should correlate across neighboring stars), (b) comparison against the 3D extinction reference maps at member distances, (c) the known A_V range for Cyg OB2 (~4–8 mag in the core; strong patchiness is expected and physical).
4. Repeat with R_V = 3.0 and 3.5 → branch outputs.
5. Produce per-subgroup de-reddened CMDs (absolute magnitude using each star's distance from subgroup posterior, not individual inverted parallax for faint stars).

**Outputs:** `wp3_extinction.cat` (per star, per R_V branch), `wp3_cmds/` (per subgroup, per branch), extinction map figure.
**Gate:** extinction distribution consistent with literature range; no unphysical negative extinctions beyond noise; spatial coherence demonstrated.
**Caveats:** differential extinction is THE Cygnus problem — a single regional A_V is forbidden. Degeneracy between temperature and extinction in broadband fits is real: for hot stars colors saturate, so extinction for OB stars is better derived from spectral-type-based intrinsic colors (use the anchor table). Weakness: extinction errors leak directly into masses (WP4) and completeness (WP5) — propagate, don't hide.

---

## WP4 — Subgroup ages and per-star masses

**Objective:** age posterior per subgroup; mass per member star.

**Procedure:**
1. Isochrone grids: PARSEC and MIST, solar metallicity, ages 1–10 Myr (log-spaced fine grid).
2. Per subgroup, fit the de-reddened CMD by likelihood over member stars (weight = membership P). Use both the upper main sequence AND the pre-main-sequence turn-on if detectable — two independent age indicators; report both.
3. Explicit binary treatment: include an unresolved-binary population in the CMD likelihood (fraction ~0.3–0.5 branch) rather than sigma-clipping the binary sequence away.
4. For the ~50 brightest/anchored stars: override photometric Teff/L with spectroscopic values from the anchor table (spectral type → calibration). Per-star provenance recorded (photometric vs spectroscopic).
5. Read per-star masses off the best isochrone per branch; for anchored stars, masses from spectroscopic HRD position.
6. Sensitivity: refit ages with the R_V branches from WP3; with distance fixed at both Berlanas populations.
7. Cross-check ages against literature: Wright+15 (1–7 Myr spread), Berlanas+20 burst structure (their multiple star-forming bursts should map onto your subgroups — compare explicitly).

**Outputs:** `wp4_ages.md` (subgroup × branch age posteriors, both indicators), `wp4_masses.cat` (per star: mass, method, branch values), HRD/CMD overlay figures.
**Gate:** recover an age *spread* consistent with the literature (failure to see subgroup age differences = red flag); anchored-star spectroscopic HRD consistent with chosen isochrones within quoted model differences.
**Caveats:** age is the load-bearing parameter of the whole paper — this WP earns the most careful treatment. PARSEC vs MIST young-age disagreement (~20–30%) must be carried, never averaged. Weakness: if age posteriors exceed ~±1.5 Myr per subgroup, the final verdict will be soft (see WP9 contingency). Strength: Cygnus OB2's age structure is well-studied; you are refining, not pioneering, so anomalies are detectable.

---

## WP5 — IMF normalization and completeness

**Objective:** the birth-population normalization k per subgroup, with completeness correction.

**Procedure:**
1. Completeness via injection testing: insert synthetic stars (mass → magnitudes via WP4 isochrones + WP3 extinction distribution) into the real field; run the WP2 selection; recovery fraction vs mass, per subgroup footprint.
2. Choose the calibration window: nominal 2–8 solar masses; verify with the completeness curve that the lower edge is ≥95% complete — raise it if not (e.g. 3–8 at Cygnus distance/extinction; document the choice, it is a legitimate per-field decision).
3. Fit IMF normalization per subgroup by Poisson likelihood over mass bins (members weighted by P; model = branch IMF slope × completeness).
4. Branches: slope {2.0, 2.3, 2.6} × isochrone family (masses shift with family).
5. Consistency: the summed normalization implies a total association mass — compare with literature values (~1.6×10^4 solar masses scale) as an order-of-magnitude sanity check.

**Outputs:** `wp5_completeness_curves`, `wp5_imf_norm.md` (k posterior per subgroup per branch), mass-function figure with fit overlay.
**Gate:** fit residuals consistent with Poisson scatter (no systematic bin trends = no hidden completeness slope); implied association mass within factor ~2 of literature.
**Caveats:** the window must contain enough stars per subgroup (target ≥50) — if a subgroup is too poor, merge it with its nearest kin for IMF purposes and document. Weakness: completeness modeling inherits WP3 extinction errors — bracket by re-running with the extinction branches. Strength: this step's assumptions (IMF universality) are the field's best-tested; referees rarely attack a properly-Poisson IMF normalization.

---

## WP6 — Massive-star census closure and runaway correction

**Objective:** verify the census of LIVING massive stars against the IMF prediction; recover escaped members.

**Procedure:**
1. Closure test: per subgroup, IMF-predicted count of living stars above 8 solar masses (from k, age, turnoff) vs directly observed count (members + anchors). Report the ratio with uncertainties.
2. If observed < predicted: attribute between (a) extinction-hidden stars (use completeness curve at high mass — should be ~1 there, but check the heavily embedded sightlines), (b) escaped runaways, (c) genuine IMF deviation. The attribution feeds WP7 corrections.
3. Runaway recovery: in the WIDE box catalog, select OB candidates (color/magnitude + parallax compatible with Cygnus distance); trace proper motions backward over subgroup age windows; flag stars whose paths intersect the association footprint with plausible ejection velocities (10–100 km/s). Without RVs this is 2D — treat recovered count as a lower bound; note the DR4 upgrade path.
4. Cross-match recovered runaway candidates against literature runaway catalogs and the anchor table (known Cygnus runaways exist — e.g., discussions around bow-shock stars in the region).
5. Update the living-star ledger: members + anchors + runaways (with provenance flags).

**Outputs:** `wp6_closure.md` (predicted vs observed per subgroup), `wp6_runaways.cat` (candidates, traceback parameters, confidence class), updated `wp6_massive_census.cat`.
**Gate:** closure ratio explained within quantified contributions; runaway search reproduces at least the well-known literature candidates in the field.
**Caveats:** 2D traceback produces false positives (chance alignments) — quantify with the control-field trick (trace back control-field OB stars to the association; the intersection rate = false-positive rate). Weakness: runaways ejected >3 Myr ago with modest velocities have left the wide box or blended into the field — an irreducible DR3 incompleteness; state it. Strength: this WP is what makes N_SN defensible — reviewers of Fuchs-style analyses always ask about runaways.

---

## WP7 — The supernova ledger

**Objective:** N_SN posterior, explosion timeline R_SN(t), and time-since-last-SN posterior; all per branch, all as functions of age.

**Procedure:**
1. Monte Carlo engine: for each branch and each iteration — draw subgroup age from WP4 posterior; draw star-formation duration (branch); draw a full stellar population from the IMF conditioned on reproducing the observed calibration-window counts (WP5) — i.e., discrete stochastic sampling, not expectation values; apply lifetimes τ(m) from the matching stellar-model family; determine which drawn stars have died by now.
2. Corrections per iteration: subtract recovered runaways from the "missing = dead" bookkeeping (they are alive); apply the explodability branch (all-explode vs islands-of-implosion) to convert "died" into "exploded as SN."
3. Record per iteration: N_SN, the list of explosion epochs, time since the most recent explosion.
4. Aggregate: posteriors of N_SN and last-SN time per subgroup and for the association total; R_SN(t) curves (explosions per Myr vs look-back time) with credible bands.
5. Mandatory presentation: N_SN and P(last SN < 100 kyr) as explicit FUNCTIONS of assumed age, in addition to marginalized values — this is the honesty plot for the age sensitivity.
6. Sanity anchors: compare N_SN total against Knödlseder+02/Martin+10 predictions (WP0 baseline) — explain any factor-level difference by input differences (census size, ages), not hand-waving.

**Outputs:** `wp7_ledger.md` (all posteriors, per branch), `wp7_rsn_curves`, the age-sensitivity figure, MC convergence diagnostics.
**Gate:** MC converged (posteriors stable under doubling iterations); branch spread documented; Knödlseder/Martin comparison written.
**Caveats:** small-number stochasticity dominates — expect asymmetric posteriors; never quote means without intervals. Binary mass transfer (±1-level effect) acknowledged as unmodeled or as a crude extra branch. Weakness: everything inherits WP4 ages. Strength: this WP is pure computation on frozen inputs — fully reproducible, fully delegable.

---

## WP8 — External cross-checks (the validation layer)

**Objective:** confront the ledger with every independent SN marker in the field.

**Checks, in order of constraining power:**
1. **PSR J2032+4127 (the anchor).** The pulsar is in Cyg OB2 with a Be-star companion (MT91 213). Its existence requires ≥1 SN. Compare its characteristic age (with the standard caveat that characteristic ages can be off by factors) against the ledger's last-SN posterior; check the ledger assigns non-negligible probability to ≥1 SN within the pulsar's plausible age. A ledger that makes the pulsar improbable indicates a problem (missed older subgroup, IMF issue) — treat as gate-level tension requiring resolution or honest discussion.
2. **γ Cygni SNR (G78.2+2.1).** Age ~7 kyr, distance estimates 1.5–2.6 kpc — association with Cyg OB2 unsettled. Two-sided check: does the ledger comfortably allow a ~7 kyr SN? Conversely, use your Gaia-based association distance to comment on the plausibility of physical association. Do not overclaim; the distance overlap is genuinely ambiguous.
3. **INTEGRAL 26Al.** 26Al (half-life 0.7 Myr) traces recent massive-star ejecta + SNe. Compare the ledger-implied recent activity qualitatively with the measured Cygnus 26Al flux and Martin+10's modeling. This is a consistency band, not a fit — 26Al also comes from WR winds, so it constrains the combination.
4. **Absence of radio SNRs in the cavity.** The Härer scenario predicts an invisible remnant (low-density cavity). Check Green/SNRcat for any candidate; frame the non-detection as *consistent* with the cavity argument, and quantify what the ledger + standard SNR visibility would predict for detectable remnants — this is the bridge to the Vink-flavored discussion.
5. **Neighboring associations (Cyg OB9, OB1, etc.).** The Härer SN needs to be in the cavity, not necessarily in Cyg OB2 proper. Produce a coarse (literature-based, not full-pipeline) estimate of the SN budget of the neighboring populations to bound the alternative-source probability. Label clearly as coarse.

**Outputs:** `wp8_crosschecks.md` (one section per check: inputs, comparison, verdict), tension list if any.
**Gate:** pulsar consistency resolved (agreement or documented tension with candidate explanations).
**Caveats:** characteristic pulsar ages and SNR age/distance estimates carry factor-level uncertainties — comparisons are probabilistic, not binary. Strength: this layer is what elevates the paper from an IMF exercise to a validated measurement; it is also the most referee-pleasing section.

---

## WP9 — The verdict and result framing

**Objective:** compute the exact statement the paper leads with; choose the framing.

**Procedure:**
1. Evaluate the WP0 requirement: P(≥1 SN satisfying Härer's energy/age/location window), marginalized and per branch. Include the WP8.5 alternative-source bound in a combined "SN available in the cavity" probability.
2. Sensitivity table: verdict vs {age treatment, IMF slope, explodability, isochrone family}. Identify the dominant driver.
3. Framing decision rule: if the marginalized verdict is stable across branches and its 68% interval does not straddle "plausible/implausible" — Letter framing (lead with the verdict). Otherwise — regular-article framing (lead with the ledger + constraints + DR4 forecast). Compute what DR4 improves (age precision, runaway 3D traceback) and quote the expected sharpening.
4. Pre-submission disclosure per the collaboration protocol: Vink first (draft), Brian second (verdict + draft), before arXiv.

**Outputs:** `wp9_verdict.md` (the statement, intervals, sensitivity table, framing decision + justification).
**Gate:** the verdict sentence exists, is branch-annotated, and survived a devil's-advocate pass (write the strongest referee objection and the response).

---

## WP10 — Manuscript assembly

**Objective:** the submittable manuscript.

**Structure (Letter variant):** Title/abstract (verdict-led) → Introduction (Cygnus PeV context; Härer assumption; prior SN budgets Knödlseder/Martin; the gap) → Data & Method (compressed; points to appendix) → Results (ledger + verdict + cross-checks) → Discussion (implications for cocoon/bubble interpretation; SNR visibility in cavity; DR4 outlook) → Appendix (full method, branches, provenance, validation gates).
**Structure (regular variant):** same content, method promoted to main text, verdict as one of several results.

**Figure plan (draft; ~6 figures):** (1) membership/substructure sky+PM figure; (2) de-reddened CMD with isochrone fits per subgroup; (3) mass function + IMF fit; (4) R_SN(t) with credible bands; (5) the age-sensitivity honesty plot; (6) verdict summary vs branches (+ cross-check markers overlaid on the timeline).

**Checklist:** provenance appendix compiled from WP logs; data/code availability statement (pipeline public on release — recommended); acknowledgments (Gaia/DPAC boilerplate, 2MASS, ADS); language pass; Vink iteration; Brian disclosure; arXiv posting simultaneous with journal submission.

**Outputs:** manuscript + figures + reproducibility package.
**Gate:** an external reader (Vink) can follow every number from query to verdict using the appendix alone.

---

## Dependency graph

WP0 → WP1 → WP2 → WP3 → WP4 → WP5 → WP7
WP4/WP5 → WP6 → WP7 → WP8 → WP9 → WP10
(WP8 items 1–2 data acquisition can start at WP1; WP0 dedup repeats cheaply mid-project.)

## Risk register (top 5)

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Age posteriors too wide → soft verdict | Medium | Reframe, not fail | WP4 dual indicators + spectroscopic anchors; regular-article contingency built into WP9 |
| Contamination in membership inflates IMF norm | Medium | Biases N_SN high | Control-field false-positive quantification (WP2); soft probabilities everywhere |
| Someone publishes the same test mid-project | Low-Medium | Scoop | WP0 dedup, early Brian conversation, speed; retarget option (Carina/Wd1) |
| Pulsar-ledger tension unresolvable | Low | Delays; but is itself a result | WP8 treats tension as reportable finding |
| Scope creep toward Paper 2 (L_w, gamma) | High | Schedule slip | Mission statement's OUT-of-scope list is binding; log temptations for Paper 2 |

## Suggested execution order for delegation

Each WP is a self-contained brief for a sub-task/agent session: hand it this file + `method_explained.md` + the upstream artifacts listed in its Inputs. Natural session-sized units: WP0; WP1; WP2; WP3+WP4 (coupled via isochrones); WP5+WP6; WP7; WP8; WP9+WP10.
