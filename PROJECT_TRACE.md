# Project trace — WP0 to WP9

*Master index for "The Supernova History of Cygnus OB2 from Gaia DR3" (Paper 1).
Created 2026-07-27. Every artifact produced so far is linked from here, with its
gate verdict and the numbers that downstream work depends on. This is a
navigation and audit document — it duplicates no content, it points at it.*

**Governing documents:** [paper1_execution_plan.md](paper1_execution_plan.md)
(what to do) · [method_explained.md](method_explained.md) (why) ·
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) (how every number was chosen) ·
[GLOSSARY.md](GLOSSARY.md) (nomenclature) · [AUDIT.txt](AUDIT.txt) (generated
file inventory with checksums, produced by [audit.py](audit.py))

---

## 1. Status board

| WP | Objective | Gate verdict | Downstream authorized |
|---|---|---|---|
| WP0 | Requirements + dedup | **PASS** — no direct duplicate; one Gaia DR3 partial overlap to cite | yes |
| WP1 | Data acquisition | **PASS** | yes |
| WP2 | Membership + substructure | **PASS** — 8/8 criteria, after one rejected run | yes |
| WP3 | Per-star extinction | **PASS**, then repaired (repair_v1) | yes |
| WP4 | Ages + masses | **SATISFIED WITH CONSTRAINTS**, then repaired (repair_v1) | yes, with constraints |
| WP5 | IMF normalization + completeness | **PASS — accepted 2026-07-28 (`repair_v6`).** Gate G3 clears under the **strict** per-branch reading with `A_or_C_regressions: []`; all three subgroups pass the baseline residual gate under both trend statistics; branch grid **40/54**; association mass 29,122 M☉ (**definition**: multiplicity-adjusted stellar mass over 0.08–120 M☉; 29,246 on the `repair_v7` chain — [B1 reconciliation](reports/wp5_association_mass_reconciliation.md)) | **yes** |
| WP6 | Census closure + runaways | **COMPLETE — BOTH GATE CRITERIA MET.** The census very nearly closes at Salpeter: on the current `repair_v7` chain observed massive stars exceed prediction by **6.7%** (α=2.3, grid median **1.067**) and the closing slope is **α ≈ 2.25**, with **18/18** cells inside the carried grid — up from 1.105 / α 2.230 / 17-18 under `repair_v6`, which remains the accepted-and-preserved baseline ([CUTS §18.8](CUTS_AND_THRESHOLDS.md)). Three defects were found and fixed after the first run: **issue #17** (the forward integral truncated at the 8 M☉ threshold the response smears across) absorbed ~3/4 of the excess; **issue #16** (absolute instead of peculiar proper motions) invalidated the traceback; **issue #15** (massive-star multiplicity) was measured and accounts for only **3.7%**. Runaways **119 raw → 54.9** corrected; living ledger **380.6** above 8 M☉. External gate **PASSED** — BD+43 3654 at p = 1.000, 38.8 km/s, 1.36 Myr vs literature ~40 km/s, 1.6 Myr. Three independent literature cross-checks agree (§10b). Withdrawn, not to be quoted: 45% excess, grid median 1.444, closing α 2.070, 260/109 runaways, 471.9 ledger. **Open is a science question, not a method one: CygOB2-C's residual 1.405** | see [wp6_closure.md](wp6_closure.md) |
| WP7 | The supernova ledger | **COMPLETE — 3 of 4 gate criteria met.** Baseline (PARSEC, R_V=3.1, α=2.3, coeval, all-explode): **N_SN = 8.43** for the association (median 8, 68% [5,11]), A **4.17** · B **4.26** · C **0**; P(≥1 SN) = 0.9997; **P(last SN < 100 kyr) = 0.552**; median time since last SN **0.086 Myr**. **The number is not 8.43** — across 54 all-explode branches it spans **1.93–28.74**, a factor of 14.9, driven mainly by α. **Headline branch set adopted 2026-07-30 (item D1): α ∈ {2.0, 2.3}, 36 branches, N_SN 5.63–28.74 (factor 5.1); α = 2.6 reported in the sensitivity table, not deleted; census closure (E2) deliberately not spent** — [wp7_alpha_headline_adoption.md](reports/wp7_alpha_headline_adoption.md). Two findings dominate: (i) **the whole budget lives above 30 M☉** — for any black-hole threshold ≤ **30** M☉ the ledger returns *exactly zero* on every branch, and at 40 M☉ on every coeval branch (**corrected 2026-07-30, item D1**: the previously stated "≤ 40 M☉ on every branch" is false on the 1–2 Myr formation-window branches, where the lowest turnoff falls to 33.9 M☉ and ≤1.8% of N_SN sits below 40; **L2 itself is unaffected**, having been written against the 25 M☉ islands cut — [wp7_alpha_headline_adoption.md §4](reports/wp7_alpha_headline_adoption.md)), so the result is entirely conditional on whether very massive stars explode; (ii) **CygOB2-C sits on the boundary** where the turnoff crosses the 120 M☉ IMF ceiling, contributing 0 to ~7 SNe depending on family, R_V *and* formation duration. **L1 and L6 FAILED** and are recorded as failed, not reinterpreted. **The execution plan's 'missing = dead' route was found not well-posed at design time** and the turnoff route was fixed as the measurement before any Monte Carlo ran; runaways bound the out-of-association fraction (≤14.6%) rather than entering N_SN | see [wp7_ledger.md](reports/wp7_ledger.md) |
| WP8 | External cross-checks | **COMPLETE — ALL FOUR GATE CRITERIA MET, all five predictions PASS.** **PSR J2032+4127 rules out WP7's zero-supernova branch**: a neutron star requires a successful explosion, and the ledger gives P(≥1 SN) = **0.9997** all-explode against **exactly 0** on the islands branch. **The pulsar's companion is our own star** — MT91 213 is anchor `gaia_dr3:2067835682818358400`, a B0V 17 M☉ orphan anchor positionally in CygOB2-A, 0.115° from its centroid, already counted in the living census. Its progenitor evolved first, so if coeval with A it exceeded A's 57.9 M☉ turnoff — exactly where the islands prescription predicts a black hole. Age agrees: characteristic 200.7 kyr, widened to **151–401 kyr**, against ledger P(last SN within 401 kyr) = **0.960**. Kinematics agree: peculiar transverse velocity **27.6 km/s**, low as a bound decades-long binary requires. **Three readings remain degenerate** (high-mass explodability / older population / binary stripping) and all give a non-zero budget. γ Cygni **allowed but unsettled** (P = 7.7% within 10 kyr; our 1.62 kpc vs Leahy 1.7–2.6 kpc). Absent remnants are **weak evidence** — expected visible number 0.80 even at a generous 100 kyr lifetime with no cavity. ²⁶Al kept a consistency band, never inverted. Tensions T1–T3, none gate-level | see [wp8_crosschecks.md](reports/wp8_crosschecks.md) |
| WP9 | The verdict and framing | **COMPLETE — gate met, all three predictions PASS.** **P_verdict = 0.323–0.736 (median 0.533)** over 36 headline branches, decomposing as C1 age **0.379–0.861**, C3 progenitor type **1.000** (every SN came from a **>34 M☉** star on every branch and a >52 M☉ star on the coeval ones — above the ~30 M☉ stripping threshold everywhere, so all were type Ib/c, exactly Härer's required channel; the flat ">52 M☉, far above 30" phrasing was **qualified 2026-07-30, item D1** — C3 = 1.000 is unaffected as computed), C4 in-situ **0.854**. The energy condition is **not** multiplied in and is reported as a stated conditional. **The verdict hinges on α alone**: supported on **18/18** α = 2.0 branches (0.593–0.735) and on **0/18** α = 2.3 branches (0.323–0.474), with the 0.5 boundary falling exactly between them (axis spreads: α 0.264, R_V 0.065, family 0.045, δ 0.043). **Framing rule applied mechanically → INCONCLUSIVE → regular article.** Under Härer's own permissive "few hundred kyr" allowance every branch clears 0.5 (**0.727–0.854**), so the scenario is *supported* against their permissive reading and *unresolved* against their preferred ~50 kyr; both windows were pre-registered before either was computed. α = 2.6 reported as excluded (P = 0.142–0.250), not deleted. Six-objection devil's-advocate pass, three conceded | see [wp9_verdict.md](reports/wp9_verdict.md) |
| WP10 | Manuscript assembly | **DRAFTED 2026-07-30 — validation gate PASSES, NOT YET COMPILED.** Regular A&A article per WP9's mechanical framing decision. [manuscript/main.tex](manuscript/main.tex) (~5,250 words) with **no hand-typed numbers**: all 117 quoted quantities are macros in the generated [numbers.tex](manuscript/numbers.tex), read from versioned artifacts resolved through [wp10_inputs.py](scripts/wp10_inputs.py), which refuses the 11 superseded products by name (item S3). Six figures; 3–6 newly built by [wp10_figures.py](scripts/wp10_figures.py). [wp10_validate.py](scripts/wp10_validate.py) passes all seven checks (V1 undefined macros, V2 dead macros, V3 refs, V4 citations, V5 figures, V6 forbidden inputs, V7 bare numbers). **No LaTeX toolchain in this environment, so the manuscript has not been compiled** and `aa.cls` is not vendored — both are pre-submission items, with the author list, Appendix A's branch table, and the final dedup re-sweep | [manuscript/README.md](manuscript/README.md) · [wp10_validation.json](provenance/wp10_validation.json) |
| WP11 | Post-hoc validation extensions (IR bow shocks + isotope forecast) | **PART B COMPLETE 2026-07-31 — gate met, I2–I5 PASS, I1 passes but is VACUOUS and its corrected form FAILS (recorded, not reinterpreted). PART A still not adopted** (referee-stage reserve / deferred to the DR4 rerun where 3D velocities double its value). **The punchline held and sharpened: ⁶⁰Fe splits the headline set exactly along α — 18/18 α = 2.0 branches clear COSI's 3σ 2-yr narrow-line sensitivity (3.0e-6 ph cm⁻² s⁻¹), 0/18 α = 2.3 branches do**, a factor 2.68 between medians against the pre-registered 2× bar (the "~3×" guess was disclosed as partial prior knowledge and the threshold set *below* it so I3 could fail). COSI therefore discriminates the one axis the WP9 verdict hinges on and DR4 will not fix. Primary arm predicts ²⁶Al 0.92–3.7 ×10⁻³ M☉ (SN-only) and **⁶⁰Fe 3.0–14.3 ×10⁻³ M☉**. Three published yield arms fixed before any flux: **LC06 Table 3** (primary, all masses explode), **LC06 Table 5** (Langer 1989 WR mass loss, low), **LC18 Recommended** (null — collapses everything >25 M☉, and **0 of 315M sampled SNe fall below 25 M☉**, so it predicts *identically zero*; reported not adopted, since WP8's pulsar excludes a zero-SN Cyg OB2). **Two honest bounds given equal weight: the yield arm is worth 30× against 4.8× for the whole branch set** (I4), and the existing INTEGRAL/SPI Cygnus ⁶⁰Fe limit is already within **5%** of excluding the richest branch (I2 passes by 1.06×). The brief's steady-state formula was **not** used — Cyg OB2 has accumulated only 0.16–0.31 of saturated ⁶⁰Fe (first SN ~1.4 Myr ago vs τ = 3.78 Myr), so rate×τ would have overstated it 3–6×; a decay-weighted, mass-resolved sum was pre-registered instead. **Finding T4: WP8 §5's "~1 M☉ for the whole Cygnus complex" is ~100× too large** (the frozen WP1 flux gives 8.9e-3 M☉ at 1.62 kpc; it looks like Martin+10's *Galactic* stationary mass) — so SN-only ²⁶Al is **10–42%** of the measured complex flux, not 0.1%; WP8's qualitative conclusion survives, its margin does not, and **nothing upstream was retuned**. ²⁶Al never inverted into a SN count. Rejected for Paper 1, with reasons recorded: diffuse WISE/MSX cavity morphology (Paper 2), hidden-OB search (closure test already answers it), ²⁶Al map re-derivation (resolution) | [wp11_isotope_forecast.md](reports/wp11_isotope_forecast.md) · [brief](tasks/wp11_bowshock_isotope_brief.md) |

**Current status: WP0–WP9 complete, WP10 drafted, WP11 Part B complete.** The
chain runs on `repair_v7`; WP5 is accepted, WP6/WP7/WP8/WP9 are complete, and
the framing decision (regular A&A article) was made mechanically at WP9. **WP11
Part B was executed on 2026-07-31** — the ²⁶Al/⁶⁰Fe forward prediction, one
Discussion subsection (`sec:cosi`) plus one sentence in the ²⁶Al cross-check
paragraph, `wp10_validate.py` re-run and passing all seven checks with 138
macros. It changed no published number: the only files touched outside WP11's
own are `scripts/wp10_inputs.py`, `scripts/wp10_numbers.py` and the manuscript,
verified against the `audit.py` inventory before and after. WP11 Part A
(bow-shock cross-match) remains **not adopted** for Paper 1: a ready-to-execute
reserve for the referee stage, otherwise deferred to the DR4 rerun.

**Two things the WP11 execution changed about how the paper reads.** First, the
paper now carries a **falsifiable forward prediction** rather than only a
retrospective verdict, and the prediction discriminates α — the single axis WP9
showed the verdict hinges on and the DR4 outlook admits DR4 cannot settle.
Second, **finding T4 corrects a published WP8 sentence**: supernovae supply
10–42% of the measured complex-wide ²⁶Al flux, not the ~0.1% WP8 §5 implies. The
²⁶Al cross-check is still a consistency band and still not a tension with data,
but it is a much tighter band than WP8 claimed. Per the one-way-validation rule,
WP8's own report and execution JSON were left untouched; the correction lives in
WP11 and the manuscript quotes the flux-space form, which is immune to the
error.

**The pre-WP10 brief was executed on 2026-07-30**
([pre_wp10_assessment_brief.md](tasks/pre_wp10_assessment_brief.md)):

| item | outcome | report |
|---|---|---|
| **B1** association mass | three nested integrals of one `k`, differing by 1.67×; the "agree, 5%" cross-check compared mismatched definitions and is **withdrawn as stated** | [wp5_association_mass_reconciliation.md](reports/wp5_association_mass_reconciliation.md) |
| **B2** CygOB2-B's age | both real: 3.548 is the CMD prior, **4.09 the counts-based age the ledger uses**; worth **1.92×** in B's SN count. B's fitted age **rails** against its prior grid → one-sided lower bound | [wp4_wp5_age_reconciliation.md](reports/wp4_wp5_age_reconciliation.md) |
| **B3** stale status blocks | §1, §7, §8, §9 and §10 marked historical; WP5's header now reads ACCEPTED | this file |
| **D1** drop α = 2.6 | pre-registered, then adopted. **D1-P1 FAILED** and caught a real over-generalization (see below); P2–P4 passed | [wp7_alpha_headline_adoption.md](reports/wp7_alpha_headline_adoption.md) |
| **T3** binary mass transfer | bounded at **±30%**, ~5× smaller than the branch spread, from a digitized BPASS comparison and Zapartas+2017's delay times | [wp7_binary_bound.md](reports/wp7_binary_bound.md) |
| **S1** dedup re-sweep | **no duplicate**; 4 new papers to cite, incl. a *microquasar* PeVatron candidate that widens the introduction's framing | [wp0_dedup_resweep_2026-07-30.md](reports/wp0_dedup_resweep_2026-07-30.md) |
| **S2** WP7 §5 arithmetic | closed: 322.4 + **58.2** = 380.6 (binned) vs 322.4 + 54.9 = 377.3 (unclipped); the two conventions must not be mixed inside one sum | [wp7_ledger.md §5](reports/wp7_ledger.md) |
| **S3** versioned tables only | enforced mechanically, not by discipline | [wp10_inputs.py](scripts/wp10_inputs.py) |
| **S4** obligations O1/O2/O3 | discharged, with the exact manuscript wording and location fixed for each | [wp3_obligations_discharge.md](reports/wp3_obligations_discharge.md) |

**Two corrections the brief itself needed.** (i) The retained 36-branch set's
ensemble median is **13.29, not ≈9** — 8.79 is the *full* 54-branch median, so
dropping α = 2.6 *raises* the centre by ~50% rather than only trimming a low
tail; the manuscript therefore leads with the baseline branch value, never an
ensemble median. (ii) The "entire SN budget above ~52 M☉" claim holds on the
baseline and coeval branches but not on the 1–2 Myr formation-window ones; the
defensible floor is **30 M☉** (see the WP7 row).

**Note dated 2026-07-28 (evening) — the CygOB2-B problem is
solved.** It had two independent causes, each found, separately pre-registered,
and confirmed:

1. **The injection truth model used a single age** (the upper-MS MAP), and at
   B's fitted age the PMS/Henyey isochrone fold sat on the 3.17–4.0 M☉ bin
   boundary, manufacturing a localized bump. Fixed in `repair_v4` by truth-side
   joint age–k marginalization.
2. **B's extinction was biased high by the anchor spatial prior.** The prior's
   width grew with distance from the calibrating anchors, correctly declaring
   them uninformative, while its *mean* still took those same anchors at full
   weight — the two moments contradicted each other. B's nearest anchors sit
   0.377° away against 0.089° for A and 0.139° for C. Fixed in `repair_v5` by
   replacing the neighbour median with a simple-kriging mean using the
   *already-fitted* variogram; no new parameter, and repair_v3 behaviour is
   bit-preserved under the old mode.

Fix 2 was pre-registered before it was run
([wp3_kriging_prior_prereg.json](provenance/wp3_kriging_prior_prereg.json)) and
**all four predictions were confirmed**
([wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json)):
B's A_V fell 0.359 mag; B's age rose 2.818 → **3.548 Myr** with a 68% interval
[3.389, 3.991] that reaches the age at which WP5's mass function independently
passed; A and C moved by **0.000 Myr**; and B's residual tilt flattened, the
replacement trend statistic going |T| 2.40 → **1.31**.

B's baseline residuals across the three versions:

| version | bin 0 | bin 1 | bin 2 | bin 3 | bin 4 | bin 5 | T | gate |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| repair_v3 | +1.45 | +2.41 | **+3.62** | −0.36 | −1.30 | −0.92 | −3.22 | fail |
| repair_v4 | +1.20 | +2.08 | **+2.51** | −0.24 | −0.79 | −0.54 | −2.40 | fail |
| **repair_v5** | +0.87 | +2.27 | **+0.36** | −0.24 | −0.03 | +0.17 | **−1.31** | **pass** |

**`repair_v5` is the first version whose baseline passes for all three
subgroups under *both* the incumbent and the better-calibrated replacement
trend statistic.** Branch grid **26 → 38 of 54** (37 under the replacement).

> ⚠️ **HISTORICAL — superseded 2026-07-28.** The paragraph below records the
> state of WP5 at `repair_v5`. It is **not** current status. WP5 was **accepted**
> at `repair_v6` (gate G3 passes with `A_or_C_regressions: []` under the strict
> per-branch reading; `downstream_wp6_authorized: true`) and **re-passed** under
> `repair_v7`. Issue #11 is closed, issue #13 is closed. Current verdict is the
> §1 status board. Kept for the audit trail, per the project's
> nothing-is-overwritten convention.

**~~WP5 is nevertheless still not accepted.~~** Under the strict per-branch reading
of gate G3 adopted on 2026-07-27, three A/C non-baseline cells regress against
repair_v4: two sit at trend p = 0.040 and 0.048, inside the indeterminate band
around the 0.05 threshold that **issue #11** identified and that criterion R3
has not yet resolved; the third, CygOB2-C MIST R_V = 3.5 α = 2.0, is a genuine
max-residual failure at 3.25. `downstream_wp6_authorized` stays `false`. What
blocks WP5 is now a threshold-methodology question, not an astrophysical one.

---

## 2. WP0 — Requirements extraction and dedup

| | |
|---|---|
| **Verdict** | No direct duplicate. One relevant Gaia DR3 partial overlap identified and must be cited. |
| **Deliverables** | [wp0_requirements_table.md](wp0_requirements_table.md) · [wp0_prior_art_summary.md](wp0_prior_art_summary.md) · [wp0_dedup_report.md](wp0_dedup_report.md) |

The requirements table defines the statement WP9 must eventually evaluate. Dedup
covered ADS 2023–2026 citers of Martin+10 / Berlanas+20 / Härer+25 plus ICRC and
TeVPA 2025 materials.

---

## 3. WP1 — Data acquisition

| | |
|---|---|
| **Verdict** | **PASS** · WP3 input gate READY |
| **Report** | [wp1_completion_report.md](wp1_completion_report.md) |
| **Provenance** | [wp1_provenance.md](provenance/wp1_provenance.md) · [wp1_manifest.json](provenance/wp1_manifest.json) |

**Headline numbers**

| quantity | value |
|---|---:|
| Gaia DR3 narrow-box unique sources | 245,843 |
| Gaia DR3 wide-box unique sources | 3,133,326 |
| Official 2MASS PSC matches | 228,654 |
| Complete J/H/Ks with finite errors | 203,646 |
| Wright+15 → Gaia match rate | **165 / 165 (100%)** |
| Spectroscopic anchors (unique Gaia rows) | 252 (540 literature-evidence rows) |
| ATNF pulsars / Green SNRs in wide box | 80 / 9 |

**Scripts** — [gaia_download.py](scripts/gaia_download.py) ·
[wp1_import_gaia_votable.py](scripts/wp1_import_gaia_votable.py) ·
[wp1_combine_2mass_tiles.py](scripts/wp1_combine_2mass_tiles.py) ·
[wp1_finalize_2mass.py](scripts/wp1_finalize_2mass.py) ·
[wp1_finalize_wide.py](scripts/wp1_finalize_wide.py) ·
[wp1_build_anchors.py](scripts/wp1_build_anchors.py) ·
[wp1_refreeze_anchor_uniqueness.py](scripts/wp1_refreeze_anchor_uniqueness.py) ·
[wp1_freeze_extinction_refs.py](scripts/wp1_freeze_extinction_refs.py) ·
[wp1_freeze_sn_markers.py](scripts/wp1_freeze_sn_markers.py) ·
[wp1_validate.py](scripts/wp1_validate.py)

**Data** — `data/processed/wp1_gaia_narrow.parquet` ·
`wp1_gaia_wide.parquet` · `wp1_2mass_join.parquet` ·
`wp1_spectroscopic_anchors.parquet` · `wp1_wright15_gaia_crossmatch.parquet` ·
`wp1_atnf_pulsars_wide.parquet` · `wp1_green_snrs_wide.parquet` ·
`wp1_chandra_snrcat_wide.parquet`

**Supporting docs** — [wp1_wright15_crossmatch.md](wp1_wright15_crossmatch.md) ·
[wp1_sn_markers.md](wp1_sn_markers.md) ·
[wp1_extinction_refs.md](wp1_extinction_refs.md) ·
[table1_wp1_inventory.md](tables/table1_wp1_inventory.md)

**Notebooks** — [gaia_pilot_query.ipynb](notebooks/gaia_pilot_query.ipynb) ·
[wp1_step1_build_narrow_catalogue.ipynb](notebooks/wp1_step1_build_narrow_catalogue.ipynb)

**Per-step execution records** — `provenance/wp1_*_execution.json` (one per tile
and per step), `provenance/wp1_*_job.json` (archive job IDs),
`provenance/wp1_gaia_narrow_validation.json`,
`provenance/wp1_gaia_wide_validation.json`

---

## 4. WP2 — Membership and substructure

| | |
|---|---|
| **Verdict** | **PASS** — all 8 gate criteria |
| **Deliverable** | [wp2_subgroups.md](wp2_subgroups.md) |
| **Gate table** | [table2_wp2_gate.md](tables/table2_wp2_gate.md) |
| **Provenance** | [wp2_provenance.md](provenance/wp2_provenance.md) · [wp2_membership_manifest.json](provenance/wp2_membership_manifest.json) |

**Gate results**

| criterion | threshold | achieved |
|---|---|---:|
| Berlanas+19 recall at P > 0.5 | ≥ 0.80 | 0.825 |
| Control / target yield ratio | ≤ 0.10 | 0.0366 |
| Member count | 10²–10⁴ | 1,392 |
| Central-90% l span | < 4.8° | 1.07 |
| Central-90% b span | < 4.4° | 0.89 |
| Convex-hull area | < 16.5 deg² | 2.08 |
| Largest seed / analysis sample | ≤ 0.10 | 0.0080 |
| Published-structure comparison | documented | yes |

**Sample structure** — `data/processed/wp2_members.parquet` holds 2,112 rows
(everything above P > 0.05). Of these 1,392 have P > 0.5; 1,331 carry an A/B/C
subgroup label; 61 are anchor-quality-exempt spectroscopic members; 781 are
`unassigned` (median P = 0.25) and are carried for provenance only.

| subgroup | N labelled | median P | mean P |
|---|---:|---:|---:|
| CygOB2-A | 476 | 0.97 | 0.89 |
| CygOB2-B | 426 | 0.91 | 0.84 |
| CygOB2-C | 429 | 0.91 | 0.85 |
| *unassigned* | 781 | 0.25 | 0.30 |

**Method notes** — clustering in (l, b, μα*, μδ); **parallax deliberately
excluded** (one population at 1.62 kpc, intrinsic depth 45 pc, exhausted at DR3
precision). Gaussian mixtures k = 2…8 over 50 deterministic seeds; acceptance by
**seed stability (ARI), not BIC**. Only k = 3 is stable. HDBSCAN scanned for
completeness, not used for labels. The Berlanas two-distance hypothesis is
**not confirmed** — see [wp2_subgroups.md](wp2_subgroups.md) §6.

**Scripts** — [wp2_membership_pipeline.py](scripts/wp2_membership_pipeline.py) ·
[wp2_finalize_membership.py](scripts/wp2_finalize_membership.py) ·
[wp2_derive_subgroups.py](scripts/wp2_derive_subgroups.py) ·
[wp2_distance_population_test.py](scripts/wp2_distance_population_test.py) ·
[wp2_finalize_audits.py](scripts/wp2_finalize_audits.py)

**Data** — `wp2_members.parquet` · `wp2_control_members.parquet` ·
`wp2_anchor_assignments.parquet` ·
`wp2_distance_component_probabilities.parquet` ·
[wp2_subgroup_labels.parquet](tables/wp2_subgroup_labels.parquet)

**Figures** — [wp2_membership_sky_pm.png](figures/wp2/wp2_membership_sky_pm.png) ·
[wp2_subgroups_sky.png](figures/wp2/wp2_subgroups_sky.png) ·
[wp2_subgroups_vpd.png](figures/wp2/wp2_subgroups_vpd.png) ·
[wp2_subgroups_extinction.png](figures/wp2/wp2_subgroups_extinction.png) ·
[wp2_distance_population_test.png](figures/wp2/wp2_distance_population_test.png) ·
[kdistance_min15.png](figures/wp2/kdistance_min15.png) ·
[target_kdistance_min15.png](figures/wp2/target_kdistance_min15.png)

**Paper figures** — [fig1_membership_literature.pdf](figures/paper/fig1_membership_literature.pdf) ·
[fig2_control_fields.pdf](figures/paper/fig2_control_fields.pdf), built by
[make_paper_figures_wp1_wp2.py](scripts/make_paper_figures_wp1_wp2.py)

**Notebook** — [wp2_membership_and_substructure.ipynb](notebooks/wp2_membership_and_substructure.ipynb)

**Rejected / superseded** — the 2026-07-22 baseline run is preserved, not
deleted: `data/processed/wp2_members_failed_20260722.parquet`,
[wp2_membership_manifest_failed_20260722.json](provenance/wp2_membership_manifest_failed_20260722.json),
[wp2_membership_attempt_deconvolved_failure.json](provenance/wp2_membership_attempt_deconvolved_failure.json),
[wp2_distance_population_attempt_failure.json](provenance/wp2_distance_population_attempt_failure.json).
Hyperparameter scans: `wp2_dbscan_scan.csv`, `wp2_hdbscan_scan.csv`,
`wp2_hdbscan_leaf_scan.csv`, `wp2_hdbscan_subgroup_scan.csv`,
`wp2_gmm_seed_stability.csv`, `wp2_mixture_prior_control_scan.csv`,
`wp2_literature_footprint_density_scan.csv`, all under `provenance/`.

**Literature recovery** — [table3_literature_recovery.md](tables/table3_literature_recovery.md) ·
[wp2_berlanas_recovery_audit.csv](provenance/wp2_berlanas_recovery_audit.csv) ·
[wp2_paiz_crossmatch.csv](provenance/wp2_paiz_crossmatch.csv)

---

## 5. WP3 — Per-star extinction and de-reddened CMD

| | |
|---|---|
| **Verdict** | **PASS** · WP4 input gate READY |
| **Report** | [wp3_completion_report.md](wp3_completion_report.md) |
| **Provenance** | [wp3_provenance.md](provenance/wp3_provenance.md) · [wp3_manifest.json](provenance/wp3_manifest.json) |

**Headline numbers** — distance posterior 1.6245 ± 0.045 kpc (μ = 11.05).
A_V median 5.85, span ~2–10.6 mag, core 4–8 populated. Spatial coherence
ρ = 0.66 at 18σ (Vergely cube ρ = 0.29). Two marginally negative extinctions,
none beyond 2σ.

**Method** — anchors get A_V from fixed-Teff intrinsic colours (avoids the
Teff/A_V degeneracy for hot stars); all other members get a χ²-posterior over
PARSEC+MIST templates with μ fixed. Isochrones are **shared with WP4** by
construction. R_V branches 3.0 / 3.1 / 3.5 carried in parallel from the start.

**Scripts** — [wp3_fetch_parsec.py](scripts/wp3_fetch_parsec.py) ·
[wp3_fetch_mist.py](scripts/wp3_fetch_mist.py) ·
[wp3_extinction_law.py](scripts/wp3_extinction_law.py) ·
[wp3_assemble_photometry.py](scripts/wp3_assemble_photometry.py) ·
[wp3_anchor_extinction.py](scripts/wp3_anchor_extinction.py) ·
[wp3_broadband_extinction.py](scripts/wp3_broadband_extinction.py) ·
[wp3_build_catalogue.py](scripts/wp3_build_catalogue.py) ·
[wp3_cube_comparison.py](scripts/wp3_cube_comparison.py) ·
[wp3_figures_and_coherence.py](scripts/wp3_figures_and_coherence.py) ·
[wp3_common.py](scripts/wp3_common.py)

**Data** — `wp3_extinction.parquet` · `wp3_member_photometry.parquet` ·
`wp3_anchor_extinction.parquet` · `wp3_broadband_extinction.parquet` ·
`wp3_cube_comparison.parquet` · `wp3_isochrones_parsec.parquet` ·
`wp3_isochrones_mist.parquet` (see
[wp3_isochrones_README.md](provenance/wp3_isochrones_README.md))

**Tables** — [wp3_extinction_rv30.cat](tables/wp3_extinction_rv30.cat) ·
[wp3_extinction_rv31.cat](tables/wp3_extinction_rv31.cat) ·
[wp3_extinction_rv35.cat](tables/wp3_extinction_rv35.cat)

**Figures** — [wp3_dereddened_cmd.png](figures/wp3/wp3_dereddened_cmd.png) ·
[wp3_extinction_map.png](figures/wp3/wp3_extinction_map.png) ·
[wp3_av_distribution.png](figures/wp3/wp3_av_distribution.png) ·
[wp3_wp5_truncation.png](figures/wp3/wp3_wp5_truncation.png)

**Note for WP5** — [wp3_provenance.md](provenance/wp3_provenance.md) Step 8
anticipates the WP5 truncation bias: the G < 19 limit removes 2 M☉ stars above
A_V ≈ 6.5, 3 M☉ above ≈ 8.7, 5 M☉ above ≈ 10.1 (4 Myr PARSEC). Completeness must
therefore be a function of both mass and position.

**Extinction method as adopted (2026-07-28).** The working extinction is
`repair_v5`: the anchor spatial prior's **mean** is a simple-kriging estimate
over the 8 nearest anchors using the repair_v3 variogram, and its **width** is
the unconditional variogram sigma. See
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §15 and
[WP3_KRIGED_PRIOR_repair_v5.md](reports/WP3_KRIGED_PRIOR_repair_v5.md).
Median A_V by subgroup: A 5.501, B **6.031** (was 6.391), C 6.435.
Carried systematic: broadband photometry sits ~0.5 mag below spectroscopically
calibrated anchors at matched position, so the absolute extinction and mass
scales rest on the anchor calibration (issue #1d, obligation O2).

---

## 6. WP4 — Subgroup ages and per-star masses

| | |
|---|---|
| **Verdict** | **SATISFIED** — one association-wide PASS with a subgroup-B limitation and one documented indicator disagreement · WP5/WP7 input gate READY WITH CONSTRAINTS |
| **Report** | [wp4_completion_report.md](wp4_completion_report.md) · [wp4_ages.md](wp4_ages.md) |
| **Provenance** | [wp4_provenance.md](provenance/wp4_provenance.md) · [wp4_manifest.json](provenance/wp4_manifest.json) |

**Headline result** — the three kinematic subgroups remain consistent with a
coeval upper main sequence. Retained upper-MS branch MAPs span
**3.16–4.50 Myr**; the honest envelope across both retained indicators is
**2.25–5.67 Myr**. ≈3.5–4.5 Myr is a central upper-MS summary, *not* the full
envelope. The A/B/C age ordering is unresolved, so the star-formation-duration
branch (0/1/2 Myr) remains mandatory for WP7.

| audit quantity | value |
|---|---|
| Age branches measurable | 104 / 132 (28 excluded, 9 grid-railed) |
| Mass branches retained | 6 / 6 · `mass_baseline` = PARSEC R_V=3.1, reporting only |
| HRD anchors A/B/C/unassigned | 60 / 5 / 43 / 48 |
| WP6 proper-motion hand-off | 7 candidates, none removed from membership |

**Scripts** — [wp4_fit_ages.py](scripts/wp4_fit_ages.py) ·
[wp4_masses.py](scripts/wp4_masses.py) ·
[wp4_anchors_hrd.py](scripts/wp4_anchors_hrd.py) ·
[wp4_clump.py](scripts/wp4_clump.py) ·
[wp4_figures.py](scripts/wp4_figures.py) ·
[wp4_closure_audit.py](scripts/wp4_closure_audit.py) ·
[wp4_report.py](scripts/wp4_report.py) ·
[wp4_schema_repair.py](scripts/wp4_schema_repair.py) ·
[wp4_common.py](scripts/wp4_common.py)

**Data** — `wp4_age_posteriors.parquet` · `wp4_masses.parquet` ·
`wp4_anchor_hrd.parquet` · `wp4_clump.parquet` · `wp4_posterior_curves.npz`

**Tables** — [wp4_ages_table.md](tables/wp4_ages_table.md) ·
[wp4_ages_envelope.md](tables/wp4_ages_envelope.md) ·
[wp4_ages_summary.csv](tables/wp4_ages_summary.csv) ·
[wp4_masses.cat](tables/wp4_masses.cat)

**Figures** — [wp4_age_posteriors.png](figures/wp4/wp4_age_posteriors.png) ·
[wp4_age_summary.png](figures/wp4/wp4_age_summary.png) ·
[wp4_cmd_subgroups.png](figures/wp4/wp4_cmd_subgroups.png) ·
[wp4_hrd_anchors.png](figures/wp4/wp4_hrd_anchors.png)

**Audits** — [wp4_closure_audit.json](provenance/wp4_closure_audit.json) ·
[wp4_pm_outliers.csv](provenance/wp4_pm_outliers.csv) ·
[task_B_wp4_closure_1.md](tasks/task_B_wp4_closure_1.md)

> **Read this first (2026-07-30, item B2).** Two age sets circulate and both are
> real. The **upper-MS MAPs below (A 3.981, B 3.548, C 2.512) are photometric
> CMD ages and are the PRIOR, not the adopted values.** What WP7 and WP9
> actually consume is the **WP5 fitted truth-age posterior** — the counts-based
> age — whose baseline means are **A 4.00, B 4.09, C 2.52 Myr** on `repair_v7`
> (A 4.00, B 4.07, C 2.52 on `repair_v6`, which is what obligation O1 quotes).
> Only B differs materially, by +0.54 Myr, and that difference is worth a factor
> **1.92** in B's supernova count. The paper quotes the counts-based set, with
> the upper-MS MAPs labelled as the prior.
> Full reconciliation, including the finding that B's fitted age **rails against
> the top of its own prior node support** (53% of the posterior on the topmost
> of nine nodes, so B's N_SN is a one-sided lower bound):
> [wp4_wp5_age_reconciliation.md](reports/wp4_wp5_age_reconciliation.md) ·
> [wp4_wp5_age_reconciliation_execution.json](provenance/wp4_wp5_age_reconciliation_execution.json).

**Ages as adopted (repair_v5, 2026-07-28).** Upper-MS MAP ages on the baseline
branch: **A 3.981, B 3.548, C 2.512 Myr**. CygOB2-B moved +0.730 Myr from
repair_v3's 2.818 when its extinction was corrected; A and C did not move at
all. The total spread is unchanged at 1.469 Myr, but the association now reads
as **two older subgroups (A, B) plus one younger (C)** rather than one older
plus two younger — a revised star-formation history that WP7 inherits
(obligation O1). B's PMS indicator remains unmeasurable (n = 4, grid-railed),
so issue #9's indicator disagreement is documented rather than resolved.

---

## 7. WP5 — IMF normalization and completeness · **ACCEPTED** (repair_v6, re-passed under repair_v7)

> ⚠️ **This section below the header is a HISTORICAL record of the blocked
> era (2026-07-23 → 2026-07-28) and is not current status.** WP5 was accepted on
> 2026-07-28 at `repair_v6` and re-passed under `repair_v7` on 2026-07-29 — see
> the §1 status board, [§8 repair_v6](#repair_v6--truth-age-isochrone-interpolation-2026-07-28--accepted-working-version)
> and [wp5_completion_report_repair_v6.md](wp5_completion_report_repair_v6.md).
> Current baseline: all three subgroups pass under both trend statistics,
> branch grid 40/54, association mass 29,122 M☉ (29,246 on `repair_v7`).
> The blocked-era text is preserved unedited below because every downstream
> repair argument refers to it. Marked historical 2026-07-30, item B3.

| | |
|---|---|
| **Verdict** | **ACCEPTED 2026-07-28 (`repair_v6`), re-passed under `repair_v7`.** *(Historical: **BLOCKED** — mass-function residual gate fails)* |
| **Report (current)** | [wp5_completion_report_repair_v6.md](wp5_completion_report_repair_v6.md) · [wp5_imf_norm_repair_v6.md](wp5_imf_norm_repair_v6.md) · normalization consumed downstream: `data/processed/wp5_imf_normalization_repair_v7.parquet` |
| **Report (historical, pre-repair 0/54)** | [wp5_completion_report.md](wp5_completion_report.md) · [wp5_imf_norm.md](wp5_imf_norm.md) (copy at [tables/wp5_imf_norm.md](tables/wp5_imf_norm.md)) — **frozen record, issue #2; must not reach the manuscript (item S3)** |
| **Provenance** | [wp5_provenance.md](provenance/wp5_provenance.md) · [wp5_manifest.json](provenance/wp5_manifest.json) · [wp5_validation.json](provenance/wp5_validation.json) |

**Gate definition** ([wp5_fit_imf.py:297-299](scripts/wp5_fit_imf.py#L297-L299)) —
a branch passes only if `chi_p ≥ 0.01` **and** `trend_p ≥ 0.05` **and**
`max_abs_pearson_residual ≤ 3.0`, for all three subgroups.

**State at repair_v5** *(historical — the accepted baseline is `repair_v6`, §8)* — baseline (PARSEC, R_V = 3.1, α = 2.3):

| subgroup | χ² p | trend p (incumbent) | trend p (replacement) | max abs residual | pass |
|---|---:|---:|---:|---:|:--:|
| CygOB2-A | 0.448 | 0.468 | 0.590 | 1.34 | yes |
| CygOB2-B | 0.293 | 0.111 | 0.186 | **2.27** (was 3.62) | **yes** |
| CygOB2-C | 0.350 | 0.957 | 0.858 | 1.56 | yes |

This is the **first version to pass the baseline under both trend statistics**.

Branch grid across versions: **0/54** (original) → **31/54** (repair_v1) →
**33/54** (repair_v2) → **26/54** (repair_v3) → **29/54** (repair_v4) → **38/54** (repair_v5). The
repair_v3 dip is real and explained in §8 — the corrected extinction prior
unmasked the CygOB2-B residual rather than causing it, and repair_v4 then
removed it at the source.

Association mass at baseline: **29,185 [28,116, 30,275] M☉** (repair_v5), within a
factor two of the 16,500 M☉ literature scale.

**~~WP5 is still not accepted.~~** *(historical, true at `repair_v5`; accepted at
`repair_v6` — issues #11 and #13 both closed.)* Under the incumbent trend statistic the block is
issue #11 (Monte-Carlo instability) via the strict reading of gate G3; under
the pre-declared replacement statistic B's baseline fails outright on a
surviving mass-dependent tilt (issue #12, trend p = 0.017). The table above
therefore reads "pass" only against a trend test now known to be underpowered
against exactly this residual shape. Gate record:
[wp5_repair_v4_gate.json](provenance/wp5_repair_v4_gate.json); verdict report:
[wp5_completion_report_repair_v4.md](wp5_completion_report_repair_v4.md).

324,000 catalogue injections. Association mass 30,696 M☉ (within a factor 2 of
the 16,500 M☉ literature scale). No subgroup reaches an absolute 95%
completeness edge; bright-plateau completeness is ~0.79–0.80.

**Blocker 1 — top-bin parent truncation: FIXED in repair_v2.** The forward
response integrated the parent mass distribution only to 8 M☉ while the mass
estimator has 22–25% scatter, so the ~300 living stars above 8 M☉ that scatter
down into the top bin had no term in the model. Diagnosis:
[WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md).
Execution and outcome:
[WP5_PARENT_RANGE_FIX_repair_v2.md](reports/WP5_PARENT_RANGE_FIX_repair_v2.md).
Top-bin residuals fell in all three subgroups (A −1.16, B −0.82, C −1.37) and
CygOB2-C now passes on every criterion.

**Blocker 2 — CygOB2-B mid-window excess: ~~OPEN~~ CLOSED 2026-07-28** (issue
#1c; both causes found — the single-age injection truth model on the PMS/Henyey
fold, fixed in `repair_v4`, and the anchor-prior extinction-scale error, fixed
in `repair_v5`; B's bump-bin residual 3.62 → 2.51 → 0.36). *Text below is the
diagnosis as it stood while the blocker was open.* Pre-existing (it was already
B's worst bin at +2.75 under repair_v1) and independent of the truncation; the
~5% drop in k from the fix pushed it over the bound. It is **not** a lower-edge
artifact — raising the calibration lower edge to 2.5 / 3.0 / 3.5 M☉ only
relocates the failure (33 → 30 → 28 → 33 branches passing). It **is** strongly
and monotonically R_V-dependent in both isochrone families (B bin-2 residual at
α = 2.3: PARSEC 1.72 / 3.02 / 4.20 and MIST 0.55 / 1.53 / 4.26 at
R_V = 3.0 / 3.1 / 3.5), which is an extinction-law signature rather than an IMF
or membership one. Likely mechanism: CygOB2-B has the narrowest A_V distribution
of the three subgroups (half-central-68% spread 0.196 mag vs 0.917 for A) and
the highest median (6.56), so an extinction-law error displaces B's mass
function rigidly instead of smearing it. Next diagnostic is the WP3 anchor
spatial prior — see the report's "Recommended next diagnostic".

**Scripts** — [wp5_common.py](scripts/wp5_common.py) ·
[wp5_injections.py](scripts/wp5_injections.py) ·
[wp5_fit_imf.py](scripts/wp5_fit_imf.py) ·
[wp5_report.py](scripts/wp5_report.py) ·
[wp5_finalize.py](scripts/wp5_finalize.py) ·
[wp5_make_notebook.py](scripts/wp5_make_notebook.py)

**Data** — `wp5_injection_response.parquet` ·
`wp5_completeness_curves.parquet` · `wp5_imf_normalization.parquet` ·
`wp5_mass_function_bins.parquet` · `wp5_imf_posterior_draws.npz` ·
`wp5_association_mass.parquet`

**Tables** — [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) **(stale — pre-repair,
0/54 passing; must be regenerated)** ·
[wp5_association_mass.csv](tables/wp5_association_mass.csv)

**Figures** — [wp5_completeness_curves.png](figures/wp5/wp5_completeness_curves.png) ·
[wp5_mass_function.png](figures/wp5/wp5_mass_function.png) ·
[wp5_association_mass.png](figures/wp5/wp5_association_mass.png)

**Notebook** — [wp5_imf_normalization_and_completeness.ipynb](notebooks/wp5_imf_normalization_and_completeness.ipynb)

---

## 8. Repair history

Repairs are versioned side-by-side with the frozen originals, never overwriting
them. `repair_v1` covers WP3, WP4 and WP5 as one coupled attempt.

| | |
|---|---|
| **Brief** | [wp3_extinction_repair_brief.md](tasks/wp3_extinction_repair_brief.md) |
| **Handoff brief (discharged)** | [wp5_cygob2b_massfunction_brief.md](tasks/wp5_cygob2b_massfunction_brief.md) — written for a fresh agent while #1c was the blocker; **#1c closed 2026-07-28** |
| **Report** | [WP3_WP5_REPAIR_REPORT_repair_v1.md](reports/WP3_WP5_REPAIR_REPORT_repair_v1.md) *(frozen; diagnosis superseded)* |
| **Correction** | [WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md) *(2026-07-27)* |
| **Gate record** | [wp3_repair_gate.json](provenance/wp3_repair_gate.json) — `accepted: false`, `downstream_wp6_authorized: false` |
| **Manifest** | [wp3_wp5_repair_manifest.json](provenance/wp3_wp5_repair_manifest.json) |
| **Notebook** | [wp3_extinction_repair_and_wp5_regate.ipynb](notebooks/wp3_extinction_repair_and_wp5_regate.ipynb) |

**What repair_v1 changed**

- **WP3** — 0.03-mag error floor in every band; eight-anchor spatial prior with
  measured widths 0.452 / 0.453 / 0.475 mag; full gridded extinction posterior;
  hidden template-branch width calibrated on the 149 spectroscopic anchors.
- **WP4** — ages refitted from the repaired extinction catalogue; masses are
  posterior samples over the full extinction distribution, all six bands, the
  binary branch and the subgroup age posterior; six-band model width 0.38 mag;
  log-uniform mass measure, **not** an IMF prior; spectroscopic-HRD overrides
  unchanged.
- **WP5** — injections re-run through the actual repaired WP3/WP4 estimators;
  baseline response resolution raised from 16 to 64 mass draws.

**Effect** — the branch grid went from **0 / 54** passing to **31 / 54**.
The baseline still fails on CygOB2-C.

### repair_v2 — WP5 parent-range fix (2026-07-27)

| | |
|---|---|
| **Diagnosis** | [WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md](reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md) |
| **Report** | [WP5_PARENT_RANGE_FIX_repair_v2.md](reports/WP5_PARENT_RANGE_FIX_repair_v2.md) |
| **Provenance** | [wp5_injections_repair_execution_repair_v2.json](provenance/wp5_injections_repair_execution_repair_v2.json) · [wp5_imf_fit_execution_repair_v2.json](provenance/wp5_imf_fit_execution_repair_v2.json) |
| **Scope** | WP5 only. Consumes the `repair_v1` WP3/WP4 upstream unchanged; every `repair_v1` artifact is preserved byte-identical. |

**What changed** — `MASS_GRID` extended 0.5–8.0 → **0.5–18.0 M☉** (45 points),
new `PARENT_MASS_HI = 18.0`; the forward-response clip in `wp5_fit_imf.py` moved
from `CALIBRATION_HI` to `MASS_GRID.max()`. `CALIBRATION_HI` stays 8.0 — the
*observed* window is unchanged. New `--output-version` / `--wp5-version` flags
keep WP5-only re-runs from overwriting earlier versions.
[CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) gained §7.2 (window edge vs
parent range) and a §7.1 amendment.

**Effect** — 324,000 injections (was 223,200); branch grid **31 → 33 / 54**;
top-bin residuals fell in all three subgroups; CygOB2-C cleared. Baseline still
blocked, now on CygOB2-B bin 2.

**Tables** — [wp5_imf_norm_repair_v2.csv](tables/wp5_imf_norm_repair_v2.csv) ·
[wp5_baseline_residuals_repair_v2.csv](tables/wp5_baseline_residuals_repair_v2.csv) ·
[wp5_completeness_baseline_repair_v2.csv](tables/wp5_completeness_baseline_repair_v2.csv)

### repair_v3 — WP3 variogram anchor prior (2026-07-27) · **negative result**

| | |
|---|---|
| **Report** | [WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md](reports/WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md) |
| **Chain runner** | [run_repair_v3_chain.sh](scripts/run_repair_v3_chain.sh) — WP3 → WP4 ages → WP4 masses → WP5 injections → WP5 fit |
| **Diagnostic** | [wp3_anchor_prior_diagnostic.py](scripts/wp3_anchor_prior_diagnostic.py) → [wp3_anchor_prior_diagnostic_execution.json](provenance/wp3_anchor_prior_diagnostic_execution.json) |
| **Lower-edge scan** | [wp5_lower_edge_scan.py](scripts/wp5_lower_edge_scan.py) → [wp5_lower_edge_scan_execution.json](provenance/wp5_lower_edge_scan_execution.json) |
| **Scope** | Full WP3→WP5 chain. repair_v1 and repair_v2 artifacts untouched. |

**Diagnosis (confirmed).** The repair_v1 prior gives each member the median A_V
of its 8 nearest spectroscopic anchors with one global width per R_V, calibrated
by leave-one-out at the anchor density (0.071°). CygOB2-B's members sit 0.377°
from their anchors, where the fitted variogram (nugget 0.000, sill 1.228 mag,
range 0.853°) says the required width is 1.052 mag — **understated 2.32×**. The
prior collapsed B's differential extinction **9.5×** (1.860 → 0.196) versus 2.0×
for A, and the independent Vergely/Dharmawardena cubes rank B's patchiness at or
above A's.

**Fix.** New `ANCHOR_PRIOR_MODE` in
[wp3_repair_common.py](scripts/wp3_repair_common.py): per-star prior width read
off the fitted variogram at each star's own 8th-anchor separation. Both
`WP_REPAIR_VERSION` and `WP3_ANCHOR_PRIOR_MODE` are environment variables
defaulting to repair_v1 behaviour, so the chain re-runs at any version with no
script edits and repair_v1 stays exactly reproducible. Extinction outcome was as
intended: B's spread restored 3.1× (0.196 → 0.603) while A and C moved −7%/+4%.

**Outcome — the fix did not work.** Branch grid **33 → 26 / 54**; CygOB2-B
baseline residual **3.02 → 3.62**, χ² p 0.036 → 0.0002. B's residuals sharpened
into a monotone tilt (+1.45, +2.41, +3.62, −0.36, −1.30, −0.92): its mass
function is steeper than the model across the entire window. The over-tight
repair_v1 prior had been partly *masking* this.

**Lower-edge scan also negative.** Raising B's calibration lower edge makes B
monotonically worse (3.62 → 5.11 → 5.23 at edges 2.0 / 2.5 / 2.75), which is the
opposite of a completeness-ramp artifact.

**Version recommendation (SUPERSEDED 2026-07-28 by repair_v5, which keeps the
variogram width but fixes the prior mean — see below).** Keep repair_v3 as the working extinction on the
merits of the variogram argument, *not* on gate score; picking repair_v2 because
it scores better would be selecting the extinction model by its residuals.
Caveat carried: with a ~1.05 mag prior width, B's A_V is now set almost entirely
by broadband photometry, reintroducing the Teff/A_V degeneracy the anchors were
meant to break. This is an explicit science decision, not a default.

**Repair scripts** — [wp3_extinction_repair.py](scripts/wp3_extinction_repair.py) ·
[wp3_repair_common.py](scripts/wp3_repair_common.py) ·
[wp3_repair_make_notebook.py](scripts/wp3_repair_make_notebook.py) ·
[wp3_wp5_repair_finalize.py](scripts/wp3_wp5_repair_finalize.py) ·
[wp4_mass_posteriors_repair.py](scripts/wp4_mass_posteriors_repair.py) ·
[wp4_repair_common.py](scripts/wp4_repair_common.py) ·
[wp5_injections_repair.py](scripts/wp5_injections_repair.py)

**Repair artifacts** — `data/processed/*_repair_v1.parquet` / `.npz` ·
[wp3_extinction_mass_invariant_repair_v1.csv](tables/wp3_extinction_mass_invariant_repair_v1.csv) ·
[wp5_baseline_residuals_repair_v1.csv](tables/wp5_baseline_residuals_repair_v1.csv) ·
[wp5_completeness_baseline_repair_v1.csv](tables/wp5_completeness_baseline_repair_v1.csv) ·
[wp3_wp5_repair_gates_repair_v1.png](figures/wp3_repair/wp3_wp5_repair_gates_repair_v1.png) ·
execution records `provenance/wp*_repair*_execution.json`

### Steps 1–2 of the #1c gated plan (2026-07-27, evening) · **cause identified**

| | |
|---|---|
| **Plan** | [wp5_cygob2b_age_caustic_fix_brief.md](tasks/wp5_cygob2b_age_caustic_fix_brief.md) (status header updated with results) |
| **Report** | [WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md) |
| **Gate records** | [wp5_stable_label_refit_execution.json](provenance/wp5_stable_label_refit_execution.json) (G1) · [wp5_age_scan_execution.json](provenance/wp5_age_scan_execution.json) (G2) · [wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json) (step-3 forecast) |
| **Scope** | Diagnostics only; every stored repair_v3 artifact preserved byte-identical. Snapshots `data/processed/wp5_age_scan_B_*.parquet` are new files, hashed in the scan JSON. |

**Step 1 (G1 PASSED, issue #6 closed).** Per-star label confusion over the 50
frozen GMM seeds ([wp2_label_stability.py](scripts/wp2_label_stability.py)):
1327/1331 members ≥90% seed-stable, cross-assignment ≤ 0.16%, no mass trend in
B (ρ = 0.011, p = 0.84). Stable-label refit of the full 54-branch grid moves
B's baseline bin-2 residual by **0.009** (limit 0.5); grid stays 26/54.
Contamination by mislabeling excluded.

**Step 2 (G2 PASSED).** B-only baseline injections with the truth age forced
to eight native isochrone ages, recovery side untouched, official `fit_one`
throughout ([wp5_age_conditional_scan.py](scripts/wp5_age_conditional_scan.py)):
bin-2 residual strictly monotone in truth age — **5.49 / 4.47 / 3.77 / 2.85 /
2.29 / 1.81 / 1.22 / 0.74 at 2.24 / 2.51 / 2.82 / 3.16 / 3.55 / 3.98 / 4.47 /
5.01 Myr** — with the entire tilt flattening in step; the control node (2.82 =
production truth) reproduces the stored failure (MC floor ~0.15). Nodes
**3.162 Myr (inside B's 68% interval, 2/9 posterior weight) and 3.548 Myr
(1/9) pass all three gate statistics** → the anomaly is the single-age
injection truth model sitting on the PMS/Henyey fold, displaced by an age
error, exactly finding F2/F3. B's PMS MAP itself is unmeasurable (n = 2,
grid-railed at 10 Myr); the older probes stand in for it.

**Plain truth-age marginalization is insufficient**: both posterior mixtures
fail (bin-2 3.40 / 2.89 with 1 Myr SF spread) because B's WP4 posterior is
bottom-heavy (σ_lo ≈ 0.02 vs σ_hi ≈ 0.43 Myr; 6/9 nodes on the MAP). The
brief-§4-authorized **joint age–k fit (WP4 posterior as prior, Jeffreys k
integrated analytically) passes under both prior variants** — max|res| 2.48
(χ²p 0.014) and 2.06 (χ²p 0.063), with B's joint age posterior concentrating
at 3.5–4.0 Myr, toward the PMS-indicator direction of issue #9. **repair_v4
must therefore implement the joint age–k fit uniformly for all three
subgroups**, with the anti-tuning rule intact (no gate statistic enters the
node weights).

**Steps 3a–3b (2026-07-27) — machinery implemented and validated.**
[wp5_joint_age_fit.py](scripts/wp5_joint_age_fit.py) implements the joint fit;
`inject_curve` gained an optional `truth_age_override` whose default preserves
the frozen repair_v1–v3 behaviour exactly. The node rule adds **no new
parameter** — it is the recovery side's own `age_posterior_nodes`
discretization snapped to native isochrone ages, needing **36 node responses**
for the whole 54-branch grid. Validation
([wp5_joint_fit_baseline_check_execution.json](provenance/wp5_joint_fit_baseline_check_execution.json),
zero new injections): **single-node equivalence is exact** — where a subgroup
spans one native age (A and C on the baseline) `fit_joint` reproduces the
official `fit_one` bit-for-bit, max k-draw difference **0.0**, so the new
machinery provably cannot move such a subgroup; and the **baseline branch now
passes for all three subgroups** — B's bin-2 residual 3.62 → **2.51**
(χ²p 0.0002 → 0.023) with k moving only −3.6%, A 1.36 and C 1.53 unchanged to
within Monte-Carlo noise (≤0.012 residual, ≤0.16% in k, from RNG stream
position). B's fitted truth-age posterior is 2.818: 0.02 · 3.162: 0.36 ·
3.548: 0.62.

**Step 3c is written and ready but NOT yet run** —
[run_repair_v4_chain.sh](scripts/run_repair_v4_chain.sh) →
[wp5_injections_agenodes.py](scripts/wp5_injections_agenodes.py) (33 node
injections; 3 of 36 reused from the G2 scan; measured 49 s per node) then
[wp5_fit_imf_joint.py](scripts/wp5_fit_imf_joint.py) (54-branch fit,
repair_v4 products, gate G3 with an explicit A/C per-branch regression check).
repair_v4 is **WP5-only**: WP3/WP4 repair_v3 products are consumed unchanged,
as repair_v2 consumed repair_v1.

### repair_v4 — WP5 truth-side joint age–k fit (2026-07-27) · **baseline passes**

| | |
|---|---|
| **Chain runner** | [run_repair_v4_chain.sh](scripts/run_repair_v4_chain.sh) — node injections → 54-branch joint fit |
| **Scripts** | [wp5_joint_age_fit.py](scripts/wp5_joint_age_fit.py) · [wp5_injections_agenodes.py](scripts/wp5_injections_agenodes.py) · [wp5_fit_imf_joint.py](scripts/wp5_fit_imf_joint.py) |
| **Provenance** | [wp5_injections_agenodes_execution_repair_v4.json](provenance/wp5_injections_agenodes_execution_repair_v4.json) · [wp5_imf_fit_execution_repair_v4.json](provenance/wp5_imf_fit_execution_repair_v4.json) · [wp5_report_execution_repair_v4.json](provenance/wp5_report_execution_repair_v4.json) |
| **Report** | [wp5_imf_norm_repair_v4.md](wp5_imf_norm_repair_v4.md) (data-driven; `wp5_report.py --wp5-version`) |
| **Scope** | **WP5 only.** WP3 extinction and WP4 ages/masses are the repair_v3 products, consumed unchanged. Every repair_v1–v3 artifact preserved byte-identical. |

**The fix.** The injection truth age is no longer a single upper-MS MAP. The
forward response is the WP4-posterior mixture over truth-age nodes, with node
weights updated by the Poisson likelihood of the observed counts (Jeffreys k
integrated analytically). Node rule adds **no new parameter** — it is the
recovery side's own `age_posterior_nodes` discretization snapped to native
isochrone ages. Identical machinery for all three subgroups; no gate threshold
moved; no per-subgroup age chosen. With a single node the estimator reduces to
the unmodified `fit_one` **bit-for-bit**.

**Result — the baseline branch passes for all three subgroups:**

| subgroup | max abs residual v3 → v4 | χ² p v3 → v4 | k v3 → v4 | pass |
|---|---:|---:|---:|:--:|
| CygOB2-A | 1.37 → 1.45 | 0.350 → 0.380 | 1724 → 1708 | yes |
| **CygOB2-B** | **3.62 → 2.51** | **0.0002 → 0.023** | 1925 → 1855 | **yes (was no)** |
| CygOB2-C | 1.52 → 1.52 | 0.341 → 0.492 | 1884 → 1873 | yes |

Branch grid **26 → 29 / 54**. Association mass 30,155 [28,987, 31,312] M☉,
within a factor two of the 16,500 M☉ literature scale. B's fitted truth-age
posterior mean is 3.40 Myr on the baseline (3.16–3.96 across branches), above
its 2.82 Myr upper-MS MAP and toward the PMS-indicator direction of issue #9.

**Open verdict question — the "A and C must not regress" clause of gate G3.**
Two CygOB2-C cells (MIST, R_V = 3.1, α = 2.0 and 2.3) flip pass → fail, so the
per-branch reading of that clause evaluates **false** while the baseline
reading passes. The flips are **not** caused by the fix, and this is measured,
not asserted ([wp5_trend_stability_check_execution.json](provenance/wp5_trend_stability_check_execution.json)):

- Of the 7 pass/fail flips between v3 and v4, **3 are the model change** (all
  fail → pass, driven by large max-residual improvements: B baseline α = 2.3
  3.63 → 2.51, α = 2.6 3.77 → 2.66, C MIST R_V = 3.5 α = 2.0 3.03 → 2.27).
- The other **4 are single-truth-age-node cells where the estimator is provably
  identical to repair_v3** — 2 flip each way (A MIST R_V = 3.1 α = 2.0 and
  A MIST R_V = 3.5 α = 2.6 gain; C MIST R_V = 3.1 α = 2.0 and α = 2.3 lose).
  In all four the χ² and max residual barely move; only the rank-based trend
  statistic crosses 0.05.
- **Paired test:** fitting the *identical* estimator to the two independent
  Monte-Carlo realizations of the *same* model (the repair_v3 injection
  response and the repair_v4 node response, generated at the same truth age)
  gives **opposite gate verdicts** in all four cells, and `trend_p` is
  completely insensitive to the fit RNG seed. The flips are injection
  Monte-Carlo noise.
- **Mechanism:** with 6 bins the two-sided Spearman p-value lives on a coarse
  lattice (ρ = 0.771 → p = 0.072; 0.829 → 0.042; 0.943 → 0.005), so a small
  realization change in a *well-fitting* residual vector can cross p = 0.05
  while χ² and max residual are unchanged or better. The instability is
  unbiased — it moved two cells each way.

This is a genuine weakness of the trend diagnostic, not of repair_v4. Resolving
it is a science decision recorded in §10; **no threshold was moved and
`downstream_wp6_authorized` has not been flipped.**

### repair_v5 — kriged anchor prior (2026-07-28) · **ADOPTED WORKING VERSION**

| | |
|---|---|
| **Decision record** | [wp3_kriging_adoption.json](provenance/wp3_kriging_adoption.json) |
| **Report** | [WP3_KRIGED_PRIOR_repair_v5.md](reports/WP3_KRIGED_PRIOR_repair_v5.md) |
| **Method spec** | [CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §15 |
| **Pre-registration / outcome** | [wp3_kriging_prior_prereg.json](provenance/wp3_kriging_prior_prereg.json) · [wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json) |
| **Downstream impact** | [wp3_kriging_downstream_impact.json](provenance/wp3_kriging_downstream_impact.json) |
| **Gate record** | [wp5_repair_v5_gate.json](provenance/wp5_repair_v5_gate.json) · report [wp5_completion_report_repair_v5.md](wp5_completion_report_repair_v5.md) |
| **Chain runner** | [run_repair_v5_chain.sh](scripts/run_repair_v5_chain.sh) — full WP3 → WP5 |
| **Scope** | Full chain. repair_v1–v4 artifacts preserved byte-identical; `WP3_ANCHOR_PRIOR_MODE=variogram\|global` still reachable. |

**The change.** The anchor spatial prior's **mean** is now a simple-kriging
estimate using the variogram *already fitted in repair_v3*, replacing the plain
median of the 8 nearest anchors. The **width is unchanged**. No new parameter.
The prior previously widened its uncertainty toward the sill for distant-anchor
stars — declaring those anchors uninformative — while still centring on them at
full weight; the two moments contradicted each other.

**Why the prior and not the photometry was at fault** — four independent lines:
the misfit is a *pure* A_V offset (per-band residual ÷ extinction coefficient
constant to <0.1 mag); optical-only and near-IR-only A_V agree to 0.04–0.07 mag
with B *not* the outlier; the prior-vs-photometry disagreement grows with anchor
distance **only in B** (ρ = +0.42, p = 1×10⁻¹⁸, versus −0.13 for A and −0.08
n.s. for C); and both 3D dust cubes rank B *least* extinguished while the
uncorrected prior ranked it *most*.

**Pre-registered, all four predictions confirmed:** B's A_V −0.359 mag; B's age
2.818 → **3.548 Myr** (68% [3.389, 3.991], reaching the age at which WP5's mass
function independently passed); A and C ages moved **0.000 Myr**; B's tilt
statistic |T| 2.40 → **1.31**.

**Result.** B's bump-bin residual **3.62 → 2.51 → 0.36** across v3/v4/v5 and the
high-mass deficit −1.30/−0.92 → −0.03/+0.17. Branch grid **26 → 29 → 38 of 54**.
`repair_v5` is the **first version whose baseline passes for all three subgroups
under both trend statistics**. Together with repair_v4 this closes finding F1 of
the fix brief: the *bump* was the isochrone fold under a single-age truth model,
the *mass-scale displacement* was this extinction error.

**Still not accepted.** Three A/C non-baseline cells regress (A PARSEC R_V=3.1
α=2.0 at trend p 0.040; A MIST R_V=3.5 α=2.6 at 0.048; C MIST R_V=3.5 α=2.0 at
max|r| 3.25). Two sit inside issue #11's indeterminate band.
`downstream_wp6_authorized` stays **false**.

**Downstream effects.** Association mass 30,155 → 29,185 M☉ (**−3.2%**, inside
the ±18% branch spread, so N_SN moves ~3%). k: A +2.2%, **B −12.2%**, C +1.0%.
The consequential change is structural, not numerical: the age spread is
unchanged (1.469 Myr) but B moves from near C to near A, so the association
reads as **two older subgroups plus one younger** rather than one older plus
two younger — a different star-formation history, which WP7 inherits.
Anchor selection bias was tested and is **absent** (anchors read +0.500 mag
*higher* than co-located members' prior-free photometry, p = 4.3×10⁻¹⁶).

**Obligations accepted with adoption:** O1 report the revised SF history as a
result; O2 carry the anchor absolute-scale systematic (photometry sits ~0.5 mag
below anchors, so the mass scale and N_SN inherit the anchor calibration);
O3 carry B's calibration asymmetry (4 anchors vs 59 and 42; Gaia XP spectra the
natural remedy); ~~O4 resolve issue #11's R3 before WP5 acceptance~~ —
**O4 discharged 2026-07-28**, see
[WP5_TREND_STATISTIC_R3_RESOLUTION.md](reports/WP5_TREND_STATISTIC_R3_RESOLUTION.md).

**Correction to the "three regressing cells" line above (2026-07-28).** Issue
#11's resolution re-classified them. The two CygOB2-A cells are **indeterminate
in both `repair_v4` and `repair_v5`** — their verdict is a coin flip on the
injection realization (π 0.46 → 0.30 and 0.79 → 0.78) — and are discounted under
[CUTS §14.7](CUTS_AND_THRESHOLDS.md). **One determinate regression survives:**
CygOB2-C MIST R_V = 3.5 α = 2.0, determinate pass in v4 (π = 1.000), and it is
issue #13's node-snapping artifact rather than a physical failure.

### repair_v6 — truth-age isochrone interpolation (2026-07-28) · **ACCEPTED WORKING VERSION**

| | |
|---|---|
| **Pre-registration** | [wp5_node_interpolation_prereg.json](provenance/wp5_node_interpolation_prereg.json) — written before any injection |
| **Mechanism evidence** | [wp5_node_rule_continuity_execution.json](provenance/wp5_node_rule_continuity_execution.json) · [wp5_node_rule_continuity.py](scripts/wp5_node_rule_continuity.py) |
| **Report** | [WP5_NODE_INTERPOLATION_repair_v6.md](reports/WP5_NODE_INTERPOLATION_repair_v6.md) |
| **Outcome** | [wp5_node_interpolation_outcome.json](provenance/wp5_node_interpolation_outcome.json) — all four predictions confirmed |
| **Gate record** | [wp5_repair_v6_gate.json](provenance/wp5_repair_v6_gate.json) · report [wp5_completion_report_repair_v6.md](wp5_completion_report_repair_v6.md) |
| **Chain runner** | [run_repair_v6_chain.sh](scripts/run_repair_v6_chain.sh) — WP5 only |
| **Backward compatibility** | [wp5_v6_backward_compatibility.json](provenance/wp5_v6_backward_compatibility.json) — all four checks pass |
| **Scope** | **WP5-only.** WP3 extinction and WP4 ages and masses consumed from `repair_v5` unchanged, exactly as `repair_v4` consumed `repair_v3`. repair_v1–v5 artifacts byte-identical. |

**The change.** The nine WP4 age-posterior quantiles are no longer snapped to
native isochrone ages. Each keeps weight 1/9 where the posterior puts it, and
`wp5_common.load_isochrone_between_ages` builds the truth isochrone there using
the **same bracketing and linear age blend the recovery side already used**
(`wp4_repair_common._interpolate_age_sequence`). The two sides had been
inconsistent — the recovery side interpolated while the truth side snapped. No
new parameter: the node count is still `N_AGE_NODES = 9`.

**Why.** Measured, not inferred. Translating a branch's WP4 posterior by
**0.0005 Myr** moved the snapped node distribution by up to **0.0546 Myr** in
1-Wasserstein distance — a **109× amplification** — and the node count varied
between 1 and 4 within a single branch under a ±0.1 Myr translation. The
unsnapped rule is **Lipschitz-1**: it moves by exactly the shift that caused it,
the tightest continuity a rule can have. The interpolated loader was verified to
reproduce the native isochrone to **0.00e+00 mag** at every native age, so at a
node that happens to land on the grid nothing changes at all.

**Cost.** 37 → **162** node injections. No reuse is possible because every node
now sits at a new age.

**Pre-registered predictions — all four confirmed.** P1 continuity (established
before the run, Lipschitz-1); **P2** the surviving regression clears, max|r|
**3.25 → 2.94** — explicitly *not* the adoption criterion, since it is the cell
that motivated the work; **P3** the baseline passes all three subgroups under
both statistics with B's |T| = **1.43**; **P4** the change is local — worst
truth-age posterior mean shift **0.118 Myr** (limit 0.15, passing but not by a
wide margin), association mass **−0.22%**, grid +3. **Adoption rested on P3 and
P4, which bound the size of the change without prescribing its direction.**

**Result — WP5 accepted.** Gate G3 passes with `A_or_C_regressions: []` under
the **strict** per-branch reading; §14.7's refinement discounted nothing and was
not needed (`readings_agree: true`). `accepted: true`,
`blocking_reason: null`, **`downstream_wp6_authorized: true`**. Branch grid
38 → **40 of 54**. The baseline barely moved — a nine-node mixture reproduced
what a two-node mixture found, so repair_v5's physics was not an artifact of the
coarse approximation. k moved −1.1% / +0.6% / −0.3% (A/B/C), association mass
29,185 → **29,122 M☉**. Truth-age nodes per baseline branch 2/4/1 → **9/9/9**.

**Unanticipated benefit.** The indeterminate fraction of the branch grid — 41–44%
in every earlier version, which had triggered §14.7(5) as an underpowered-grid
finding — falls to **17%** (9 of 54), below the 25% limit. Mechanism, stated
rather than treated as luck: a nine-node mixture averages the forward response
over nine independent injection realizations instead of two, cutting the
response's Dirichlet uncertainty by roughly √(9/2). The finding stands for
repair_v3–v5 and **does not apply to the accepted version**.

**Backward compatibility proven** ([wp5_v6_backward_compatibility.json](provenance/wp5_v6_backward_compatibility.json)):
the default node rule reproduces the exact repair_v4/v5 node sets on disk; the
interpolating loader equals the native table to **0.00e+00 mag** at all 35
native ages; the vectorized Jeffreys k is bitwise identical; the interpolation
flag is off for every pre-v6 version.

---

## 9. Open issues register

| # | Issue | Severity | Owner WP | Status |
|---|---|---|---|---|
| 1 | Forward response truncated parent mass distribution at 8 M☉; top bin under-predicted in all three subgroups | ~~blocking~~ | WP5 | **CLOSED 2026-07-27** by repair_v2; top-bin residuals fell in all three subgroups |
| 1b | WP3 anchor spatial prior over-constrains CygOB2-B: width calibrated at 0.071° applied at 0.377°, understated 2.32×, collapsing B's A_V spread 9.5× | high | WP3 | **CLOSED 2026-07-27** by repair_v3 variogram prior. Real bug, correctly fixed — but *not* the cause of #1c |
| 1c | **CLOSED 2026-07-28.** Both causes found and fixed: the single-age injection truth model on the PMS/Henyey fold (`repair_v4`) and the anchor-prior extinction-scale error (`repair_v5`). B's baseline worst residual went **3.62 → 2.51 → 0.36** in the bump bin, the tilt statistic \|T\| **3.22 → 2.40 → 1.31**, and the grid **26 → 29 → 38** of 54. `repair_v5` is the first version whose baseline passes for all three subgroups under both trend statistics. WP5 acceptance is now held by issue #11, not by #1c. Earlier status: **PARTIALLY SOLVED — the bump is fixed, a tilt remains.** ~~CAUSE FOUND AND FIXED~~ (that claim, made 2026-07-27, was **corrected on 2026-07-28** and was too strong). One real cause was identified and removed: the single-age injection truth model sitting on the PMS/Henyey isochrone fold. The fix (truth-side joint age–k marginalization, `repair_v4`) takes CygOB2-B's baseline residual **3.62 → 2.51**, χ² p 0.0002 → 0.023, grid 26 → 29/54. **But a properly powered trend test shows B's residuals still drift systematically with mass** (+1.20, +2.08, +2.51, −0.24, −0.79, −0.54; replacement trend p = **0.017**, T = −2.40, versus A −0.39 and C +0.19 which are flat). The incumbent rank test scored this 0.111 only because the sequence is not *monotone*. This matches finding F1 of the fix brief, which decomposed B's anomaly into **a localized bump plus a global mass-scale displacement** — the age fix removed the bump; the mass-scale displacement was never addressed. Next hypothesis is therefore well-posed and already written down: issue #1d (B's extinction is set by broadband photometry alone, so a scale error displaces B's masses coherently) and step 4b of the gated plan. Original statement: CygOB2-B's observed mass function is steeper than the model across the whole 2–8 M☉ window (bin-2 residual 3.62; monotone tilt +1.45/+2.41/+3.62/−0.36/−1.30/−0.92). Robust to every instrumental explanation tested. **2026-07-27 shape diagnostic** ([wp5_bump_shape_diagnostic_execution.json](provenance/wp5_bump_shape_diagnostic_execution.json)): no single IMF slope explains B (bin-2 excess survives all α at R_V ≥ 3.1); the bump sits on the PMS/Henyey isochrone fold (2.77–3.22 M☉ at B's MAP age); injections generate truth at a single UMS MAP age and are structurally blind to an age error; contamination signatures absent in B's bin-2 stars | ~~blocking~~ | WP5←WP4 | **CLOSED 2026-07-28** (row header above). Historical trail: **cause identified 2026-07-27** — steps 1–2 of the gated plan executed ([WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md)): **G1 passed** (labels exonerated, shift 0.009; contamination excluded), **G2 passed** — B's residual is strictly monotone in the injection truth age (bin-2: 5.49 → 0.74 over 2.24 → 5.01 Myr) and ages inside B's WP4 posterior (3.162, 3.548 Myr) pass the full gate ([wp5_age_scan_execution.json](provenance/wp5_age_scan_execution.json)). Plain truth-age marginalization fails (posterior bottom-heavy); the §4-authorized **joint age–k fit passes both prior variants** ([wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json)). **Steps 3a–3c executed** — `repair_v4` chain complete ([wp5_repair_v4_gate.json](provenance/wp5_repair_v4_gate.json), [wp5_completion_report_repair_v4.md](wp5_completion_report_repair_v4.md)); `accepted: false`, `downstream_wp6_authorized: false` under the adopted strict G3 reading. The blocker is now **issue #11**, not the CygOB2-B mass function |
| 1d | **FIXED 2026-07-28 in `repair_v5`** by the kriged prior mean (`ANCHOR_PRIOR_MODE=kriging`): the prior now shrinks toward the field mean as anchor information decays, exactly as the fitted variogram says it should. Kriging weights sum to 1.000 for A and 0.992 for C but only **0.772 for B**, so the correction is intrinsically B-specific and no per-subgroup choice was made. B's A_V fell 0.359 mag and its age rose to 3.548 Myr while A and C did not move. The underlying observational limitation stands as a carried caveat — B still has only 4 spectroscopic anchors against A's 59 and C's 42, so its extinction rests on weaker calibration than its siblings', and Gaia XP spectra remain the natural way to remove that asymmetry. Original statement: with the variogram prior, CygOB2-B's A_V is set almost entirely by broadband photometry, reintroducing the Teff/A_V degeneracy the anchors were meant to break. No DR3 spectroscopic anchors lie near B's footprint | **high** (was medium) | WP3 | **CONFIRMED WITH A MEASUREMENT 2026-07-28** ([wp3_nir_extinction_crosscheck_execution.json](provenance/wp3_nir_extinction_crosscheck_execution.json), [wp3_nir_extinction_crosscheck.py](scripts/wp3_nir_extinction_crosscheck.py)). Anchor counts inside the member sample are **A 59, C 42, B 4** — the anchor famine is real and B-specific. And B's WP3 A_V sits **+0.385 mag above** what its near-IR colour excess implies, relative to A and C (Mann-Whitney p = 2.4×10⁻¹², n=334 vs 684 in the 2–8 M☉ window). The offset keeps its sign and significance at *every* candidate age for B (+0.385 at 2.82 Myr → +0.265 at 4.47 Myr), so it is not an artifact of B's disputed age. **Mechanism corrected 2026-07-28** ([wp3_band_tension_diagnostic_execution.json](provenance/wp3_band_tension_diagnostic_execution.json)): my first reading — that the optical bands were outvoting the near-IR ones — is **wrong and is retracted**. Inside the fit, the A_V implied by G/BP/RP alone and by J/H/Ks alone agree to 0.04–0.07 mag in every subgroup, and B is *not* the outlier (B +0.044 vs A +0.074). The misfit is instead a **pure A_V offset**: per-band residuals are proportional to the extinction coefficient (residual/k constant to <0.1 mag across all six bands), meaning the model has the colours right and only the extinction wrong. Its source is the **anchor spatial prior**, whose mean sits above what the photometry alone prefers by **+0.48 mag (A), +0.71 mag (B), −0.10 mag (C)**, pulling the adopted A_V up by +0.415 / +0.485 / −0.080 respectively. B's eighth-nearest anchor is **0.374°** away against 0.089° for A and 0.139° for C, so B's prior is the least locally informed. **Caveat, stated plainly:** A is pulled almost as hard as B and still passes the WP5 gate, so the prior pull is necessary-but-not-sufficient and does not by itself establish causation for B |
| 2 | [wp5_imf_norm.csv](tables/wp5_imf_norm.csv) / `.md` are the pre-repair run (0/54) and contradict the current fit. `wp5_report.py` is hardwired to the unversioned frozen products and cannot regenerate from a repair version | high | WP5 | **CLOSED 2026-07-27** — `wp5_report.py` gained `--wp5-version`, and the versioned path is fully **data-driven**: every pass/fail statement is read from the run's own gate record instead of the hardcoded blocked-gate prose. Versioned products exist for repair_v3 and repair_v4 ([wp5_imf_norm_repair_v4.md](wp5_imf_norm_repair_v4.md), [tables/wp5_imf_norm_repair_v4.csv](tables/wp5_imf_norm_repair_v4.csv), matching figures). The unversioned files are left untouched as the frozen historical record of the 0/54 run |
| 3 | **CLOSED 2026-07-28 — measured, corrected, and the correction turned out to differ from the flagged estimate in both size and sign.** The defect was real: recovery is flat in mass out to 18 M☉ (plateau **0.787** grid median), so the loss is the WP2 quality filter and does not vanish at the bright end. But the naive repair — dividing the observed count by that plateau — is **forbidden**, because the response also scatters mass estimates across the 8 M☉ threshold. The correct forward estimator is `k · ∫ dM M^(−α) R(observed above threshold \| M)`, the same construction WP5 already uses inside 2–8 M☉ extended to one open-ended bin. Its **effective completeness is 0.871** (grid median), range 0.808–**1.041**; in **6 of 54 cells it exceeds 1**, so net up-scatter beats recovery loss and the scalar divisor has the *wrong sign*. Spurious deficit if uncorrected: **12.9%**, not the ~17–20% originally flagged. A per-subgroup systematic also emerged and is now binding: **CygOB2-B sits at 0.969 against A 0.840 and C 0.852**, because B's wider mass posteriors (obligation O3) scatter more stars up — a single association-wide completeness would under-count B's explosions. Specification is binding in [CUTS_AND_THRESHOLDS.md §16](CUTS_AND_THRESHOLDS.md); must be re-run on whichever WP5 version is accepted. Original statement: bright-mass completeness plateaus at ~0.83; WP6 step 2(a) assumes ~1.0 and will manufacture a false ~17% deficit → biases N_SN high | ~~high~~ | WP6 | **CLOSED** — [wp6_bright_completeness_execution.json](provenance/wp6_bright_completeness_execution.json) · [wp6_bright_completeness.py](scripts/wp6_bright_completeness.py) |
| 4 | No subgroup reaches an absolute 95% completeness edge; `CUTS_AND_THRESHOLDS.md` §7.1 procedure is unachievable on this field | medium | WP5 | **CLOSED 2026-07-27** — formal supersession recorded in §7.1; `corrected_no_absolute95_edge` fallback in force |
| 4b | Plan WP5 step 2 / §7.1 prescribe raising the lower edge when completeness fails. Measured twice — global edges on repair_v2 (30/28/33 of 54) and per-subgroup completeness-driven edges on repair_v3 (25/28/29 of 54). Neither works; raising B's edge makes B monotonically worse. The prescribed remedy does not apply on this field | medium | WP5 | measured and documented in [wp5_lower_edge_scan_execution.json](provenance/wp5_lower_edge_scan_execution.json); do not retry blind |
| 5 | Branches failing on the top bin now confined to the R_V = 3.5 and α = 2.6 corners; plan §1.4 forbids silently dropping branches | medium | WP5 | **CLOSED 2026-07-27** — retention policy adopted as [CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §13, written *before* the repair_v4 grid was evaluated so it could not be tailored to the outcome: no branch is ever dropped; the baseline must pass; non-baseline failures are carried as a systematic statement and reported by their concentration; downstream consumes the ensemble, never the best-scoring branch |
| 6 | Subgroup-label (A vs B vs C) uncertainty never quantified — `membership_probability` measures cluster-vs-field only | medium | WP2 | **CLOSED 2026-07-27** — quantified per star over the 50 frozen GMM seeds ([wp2_label_stability.py](scripts/wp2_label_stability.py) → [wp2_label_stability_execution.json](provenance/wp2_label_stability_execution.json), per-star CSV alongside): 1327/1331 clean members are ≥90% seed-stable, seed-averaged cross-assignment ≤ 0.16%, and B shows **no mass-dependent instability** (window Spearman ρ = 0.011, p = 0.84; 341/342 window stars stable). The stable-label WP5 refit ([wp5_stable_label_refit_execution.json](provenance/wp5_stable_label_refit_execution.json)) moves B's baseline bin-2 residual by **0.009** (G1 limit 0.5) and leaves the grid at 26/54 — **labels exonerated; gate G1 of the #1c plan passed**, contamination-by-mislabeling excluded as the bump cause |
| 7 | A/B/C age ordering unresolved → star-formation-duration branch (0/1/2 Myr) stays mandatory | medium | WP4→WP7 | carried by design |
| 8 | **Berlanas+19's two-distance population is not confirmed *within the member sample*, and our test cannot settle it more strongly than that.** A latent-distance extreme-deconvolution fit gives a single population at **1.6245 kpc, depth 45.4 pc**; the two-component fit collapses to components **31 pc apart** (ΔBIC +12.17 favours one), the two components have **identical proper motions** (KS p = 0.77), three control fields agree, and **no member lies closer than 1350 pc**. **But WP2's classifier uses `FEATURES = ['parallax_corrected', 'pmra', 'pmdec']` — parallax is a clustering feature**, so a 1350 pc foreground group would have been removed *before* the test ran. Testing survivors of a parallax-based selection for a parallax split is partly circular. Correct scope: *within the WP2 member sample there is no evidence of two distances* — **not** evidence against Berlanas's 19-star foreground group. Corroborating context: Orellana+21 independently find foreground groups of 179 and 188 stars at ~1280 pc plus UCB585 at ~1460 pc, all separable by proper motion, so all three analyses agree foreground structure exists and is PM-separable; the disagreement is about membership, not existence | low for WP2, **open for the CygOB2-C residual** | WP2→WP6 | documented — [berlanas_2019_two_distance.md](cross_checks/berlanas_2019_two_distance.md) · [wp6_external_crosschecks_execution.json](provenance/wp6_external_crosschecks_execution.json). Two follow-ups carried: (1) a **parallax-blind membership test** is the only way to break the circularity; (2) **distance contamination inflates masses and therefore closure ratios** — the direction of CygOB2-C's 1.405 — and is covered by **no** WP6 alternative (A4 tested membership weighting, not distance) |
| 9 | Upper-MS vs PMS indicator disagreement in subgroup B (N = 19 tail) | low | WP4 | documented; **2026-07-27**: the #1c joint age–k forecast independently pulls B's 2–4 M☉ truth age to 3.5–4.0 Myr — the same direction as the PMS indicator ([wp5_age_joint_fit_diagnostic_execution.json](provenance/wp5_age_joint_fit_diagnostic_execution.json)) |
| 10 | Wide-box definition inconsistency noted in `CUTS_AND_THRESHOLDS.md` §4.2 | low | WP1→WP6 | resolve before runaway traceback |
| 11 | **CLOSED 2026-07-28 — the criterion was measuring the wrong thing.** Criterion R3 required both realizations of an identical model to return the same verdict. That is **unsatisfiable by any statistic**: a cell with pass-probability π disagrees across independent realizations with probability 2π(1−π), so only a degenerate test could comply. R3 is **withdrawn** and replaced by R3a/R3b/R3c ([CUTS §14.6–14.7](CUTS_AND_THRESHOLDS.md), pre-declared before measurement); R1, R2, R4 stand. **R3a PASS** — the incumbent has no achievable p between **0.0416 and 0.0724** (gap 0.031 straddling the threshold), the replacement's gap is 0.00005. **R3b PASS, and it independently convicts the incumbent** — using each cell's own injection uncertainty (Dirichlet response posterior, 400 replicates, no new parameter) the implied flip count is 1.43 [0,3], and the incumbent flipped **4 of 4**, more often than its own noise can explain; the replacement flipped 1 against 0.74 expected. **R3c: implemented, and it triggers §14.7(5)** — **41–44% of the 54-cell grid is indeterminate** (0.05 < π < 0.95) in *every* version, so the branch grid is underpowered at N_INJECT_PER_MASS = 400 and per-cell verdicts are not independent evidence. Determinate failures nevertheless fall **15 → 11 → 5** across v3/v4/v5, so repair_v5's improvement is not an artifact. **Consequence for G3:** the two A cells at trend p 0.040 and 0.048 are indeterminate in *both* v4 and v5 and are discounted; **one regression survives — CygOB2-C MIST R_V=3.5 α=2.0**, determinate pass in v4 (π = 1.000) and indeterminate in v5 (π = 0.900). That cell is issue #13, exactly as forecast, and `repair_v6` cleared it. **Postscript:** the underpowered-grid finding this raised (§14.7(5)) does **not** apply to the accepted version — `repair_v6`'s indeterminate fraction is **17%**, against 41–44% for v3–v5, because a nine-node truth-age mixture averages the response over nine injection realizations instead of two. Prior status: **NOW THE SOLE WP5 BLOCKER**. Original statement: **The WP5 residual-trend gate statistic is Monte-Carlo unstable.** With 6 bins the two-sided Spearman p-value lives on a coarse lattice (ρ = 0.771 → p = 0.072; 0.829 → 0.042; 0.943 → 0.005), so an independent injection realization of an *identical* model can move a well-fitting cell across p = 0.05 while χ² and max residual are unchanged or better. Measured directly: 4 of the 7 repair_v3→v4 gate flips are single-truth-age-node cells whose estimator is provably identical, and a paired refit on the two realizations gives opposite verdicts (2 flips each way, so the instability is unbiased). `trend_p` is insensitive to the fit RNG seed — the driver is injection noise | ~~blocking~~ | WP5 | **RESOLVED** — [WP5_TREND_STATISTIC_R3_RESOLUTION.md](reports/WP5_TREND_STATISTIC_R3_RESOLUTION.md) · [wp5_verdict_stability_execution.json](provenance/wp5_verdict_stability_execution.json) · [wp5_verdict_stability.py](scripts/wp5_verdict_stability.py) · prior evidence: [wp5_trend_stability_check_execution.json](provenance/wp5_trend_stability_check_execution.json) · [wp5_trend_stability_check.py](scripts/wp5_trend_stability_check.py). Fix by **replacing the diagnostic, never by moving the threshold** ([CUTS_AND_THRESHOLDS.md](CUTS_AND_THRESHOLDS.md) §6.4): done and evaluated — see the status column. Raising `N_INJECT_PER_MASS` reduces the noise but not the lattice coarseness |
| 13 | **CLOSED 2026-07-28 by `repair_v6`.** The fix was pre-registered and all four predictions confirmed: the blocking cell CygOB2-C MIST R_V=3.5 α=2.0 cleared (max\|r\| **3.25 → 2.94**), the baseline held under both statistics (B's \|T\| = 1.43), and the change stayed local (worst truth-age posterior mean shift 0.118 Myr against a 0.15 limit, association mass −0.22%, grid 37 → 40). Gate G3 then passed with `A_or_C_regressions: []` and **WP5 was accepted**. Report: [WP5_NODE_INTERPOLATION_repair_v6.md](reports/WP5_NODE_INTERPOLATION_repair_v6.md). Prior status: **The truth-age node set changes discontinuously under age-grid snapping.** The repair_v4 node rule snaps the nine WP4 posterior quantiles to native isochrone ages, whose grid is coarse (0.05 dex ≈ 12%), so an arbitrarily small change in the age posterior can add or delete a whole node carrying tens of percent of the weight. Measured: CygOB2-C MIST R_V=3.5 α=2.0 kept an identical fitted age (2.002 → 2.002 Myr) and moved only −0.059 mag in A_V, but its 68% upper edge fell 2.084 → 2.062 Myr, which deleted the 2.248 Myr node holding **41.6%** of the weight — and bin 0's residual jumped +2.27 → **+3.25**. This is the same class of defect as issue #11: a coarse discretization turning a continuous input into a discontinuous output | ~~blocking~~ | WP5 | **DIAGNOSED, FIXED, PRE-REGISTERED, RUNNING.** The discontinuity is now measured rather than inferred ([wp5_node_rule_continuity_execution.json](provenance/wp5_node_rule_continuity_execution.json), [wp5_node_rule_continuity.py](scripts/wp5_node_rule_continuity.py)): translating the WP4 posterior by **0.0005 Myr** moves the snapped node distribution by up to **0.0546 Myr** in 1-Wasserstein distance — a **109× amplification** — and the node count varies between 1 and 4 within a single branch under a ±0.1 Myr translation. Dropping the snap makes the rule **Lipschitz-1**: it moves in W1 by exactly the shift that caused it, the tightest continuity a rule can have. Fix implemented as `repair_v6`: `wp5_common.load_isochrone_between_ages` builds the truth isochrone by the same bracketing and linear age blend the recovery side already uses (verified to reproduce the native table to **0.00e+00 mag** at every native age), and `truth_age_nodes(snap=False)` keeps the nine quantiles where the posterior puts them. No new parameter; repair_v1–v5 are bit-preserved behind default arguments. Predictions **pre-registered before the run** in [wp5_node_interpolation_prereg.json](provenance/wp5_node_interpolation_prereg.json) — P2 the regressing cell clears (explicitly **not** the adoption criterion), P3 the baseline still passes under both statistics, P4 the change is local (posterior mean age < 0.15 Myr, mass ±5%, grid ±6 cells). 162 node injections running (37 → 162; no reuse possible, every node is at a new age). Confirmed by issue #11's independent analysis: this cell is the one surviving G3 regression, determinate pass in v4 (π = 1.000) and indeterminate in v5 (π = 0.900). Prior status: open — [wp5_node_snapping_discontinuity.json](provenance/wp5_node_snapping_discontinuity.json). Candidate fix: interpolate the truth isochrone between native ages instead of snapping. Note the two sides are currently **inconsistent** — `wp4_repair_common._interpolate_age_sequence` already interpolates on the *recovery* side while the truth side snaps. **Re-specifying R3 will not clear this cell**, so #11 and #13 are two separate prerequisites for a fair G3 re-evaluation |
| 14 | **FOUND AND FIXED 2026-07-28, while designing WP6 step 0a.** **The isochrone turnoff was a step function of age.** The natural way to read it — `load_isochrone_between_ages(family, age)['Mini'].max()` — is wrong: that loader takes the *intersection* of the two bracketing native tables' mass ranges, so its ceiling is always the **older** bracket's turnoff and is constant across each native interval. Measured: **48.0 M☉ at both 4.00 and 4.20 Myr** (PARSEC), **69.0 M☉ at both 3.00 and 3.50 Myr** (MIST). This matters because **N_SN is essentially the IMF integral above the turnoff**, so a stepped turnoff makes the paper's headline number a step function of age — the same defect class as issue #13, now in the quantity the paper is about. It also manufactured a spurious **22% PARSEC/MIST disagreement** at 4 Myr (48.0 vs 58.6) that **vanishes once interpolated** (57.9 vs 58.7). Fixed by interpolating log-log between the native tables' own maxima, using only tabulated values, plus a running minimum to remove two unphysical inversions from table sampling (PARSEC 331.4 → 336.0 M☉ at 1.995 → 2.239 Myr; MIST 69.0 → 74.9 at 3.181 → 3.571). Caveat recorded: below 3.16 Myr (PARSEC) and 2.52 Myr (MIST) the tabulated maximum is the **table's own ceiling** (300.0 and 210.2 M☉), not an evolutionary turnoff — both far above the 120 M☉ IMF limit, so the value is capped and never used. CygOB2-C lives entirely in that regime, which is another way of saying C has lost no stars yet. **WP5 is unaffected** — it injects nothing above 18 M☉, far below every turnoff, so no accepted result depends on this | ~~blocking for WP6/WP7~~ | WP6→WP7 | **FIXED** — `turnoff_mass` in [wp6_mass_extension_decision.py](scripts/wp6_mass_extension_decision.py); recorded in [wp6_mass_extension_decision.json](provenance/wp6_mass_extension_decision.json) |
| 15 | **The injection truth model under-states massive-star multiplicity, and it is the leading non-IMF explanation of WP6's closure excess.** `wp5_common.F_BINARY = 0.40` is applied at **every mass**, with `q ~ Uniform(0.1, 1)`. Measured values for O stars are far higher: Sana+2012 find a close-binary fraction **f_b ≈ 0.7 for 15–60 M☉** with a near-flat mass-ratio exponent (κ = −0.1 ± 0.6, consistent with the uniform q already used); Duchêne & Kraus 2013 give multiplicity >90% for O types against ~50% for solar type; and **for Cyg OB2 itself** Caballero-Nieves+2020 resolve companions around **47% of 74 O/early-B stars** at 0.08–10 arcsec, with 48 of 74 multiple once spectroscopic binaries are included. Direction of the bias: if real massive stars carry more unresolved companions than the truth model assumes, they are brighter than the model predicts at fixed true mass, their inferred masses are biased **high**, and more cross 8 M☉ — an apparent excess with no IMF change. Order of magnitude: f_bin 0.4 → 0.7 with flat q implies a ~20–30% mass overestimate near threshold, inflating the count above 8 M☉ by ~1.25^1.3 ≈ **1.3×** against an observed excess of 1.45×. It is also mass-dependent, matching the observed ordering A 1.09 < B 1.45 < C 1.71 | **high — gates the WP6 IMF interpretation** | WP5→WP6→WP7 | open — test pre-registered in [wp6_multiplicity_prereg.json](provenance/wp6_multiplicity_prereg.json), specification in [CUTS §18](CUTS_AND_THRESHOLDS.md); paired control/treatment injections running ([wp6_multiplicity_injections.py](scripts/wp6_multiplicity_injections.py)), scoring by [wp6_multiplicity_closure.py](scripts/wp6_multiplicity_closure.py). **V4 passed before anything ran** — with `truth_binary_fraction=None` the regenerated response is byte-identical to the published one, so the code change cannot have perturbed an accepted artifact. **Sharpened 2026-07-28 — this is also a defect in the audit's own master table.** CUTS row 20 carries the binary fraction as a Class E branch spanning **0.3–0.5** and marks it *correct — keep*. That is wrong at the top end twice over: the measured O-star value of **0.70 lies outside the carried range entirely**, and WP5 does not carry the branch at all — WP4 fits ages at f_bin ∈ {0.3, 0.4, 0.5} but the injections select `f_bin == F_BINARY` and inject a **single constant 0.40 at every mass**. The range was specified from field/low-mass literature and then applied to the whole mass axis. Neither point touches WP5's accepted result, whose calibration window stops at 8 M☉ where 0.40 is defensible — but WP6's closure window starts exactly where the specification stops being right. Row 20 corrected to *Class E, carried range under-specified above ~10 M☉* ([CUTS §18.0](CUTS_AND_THRESHOLDS.md)). **CLOSED 2026-07-28 — measured, and it is NOT the explanation.** 162 paired nodes, 3.3 h, realized f_bin control **0.396 ± 0.005** vs treatment **0.611 ± 0.011**. The grid median at α = 2.3 moves **1.103 → 1.099**: multiplicity above 8 M☉ absorbs **3.7%** of the excess, against the ≥50% M3 required. **M1 PASS** (all three subgroups fell; 48 of 54 cells). **M2 FAIL** — CygOB2-C's reduction is *smaller* than B's by ~2.7σ, the opposite of the predicted turnoff ordering, so the pre-registered consequence applies: the excess and the multiplicity effect have different mass dependences. **M3 PASS literally, FAIL on the governing relative reading** — its 1.222 threshold was derived arithmetically from the 1.444 baseline that **issue #17 withdrew**, and the corrected control arm already sits at 1.103, below the threshold *before the treatment arm changes anything*. The pre-registration was **not amended**; both readings are recorded and the relative form — what M3's own sentence states — governs. **Decision rule applied: the 'if M3 fails' branch.** The IMF reading is reported with the measured 3.7% correction applied and the remainder carried as a systematic; **no repair_v7 is triggered**; α = 2.6 stays disfavoured. **The paired design validated itself** — the control arm reproduces the published grid median to **0.002** (1.103 vs 1.105) despite a different RNG realization and a spliced sub-8 segment, so realization noise is about half the effect being measured. Consistency: the effect grows with α (0.0011 / 0.0022 / 0.0042), the same low-mass-weighting mechanism as issue #17's F2. **BINDING SCOPE LIMIT:** this tests multiplicity **above 8 M☉ only** — the range where the mechanism has *least* room to act, since a star already above the threshold cannot be scattered into the census by being brightened. The 4–8 M☉ up-scatter channel is held at f_bin = 0.40 in both arms. Permitted claim: *multiplicity above 8 M☉ does not explain the excess*. **Not permitted:** *multiplicity does not explain the excess*. Testing the sub-8 channel perturbs the accepted WP5 calibration and needs a `repair_v7` chain re-run — carried as an open recommendation. Evidence: [wp6_multiplicity_closure_execution.json](provenance/wp6_multiplicity_closure_execution.json) · [wp6_multiplicity_injections_execution.json](provenance/wp6_multiplicity_injections_execution.json) · [CUTS §18.6](CUTS_AND_THRESHOLDS.md) |
| 16 | **The runaway traceback used absolute proper motions instead of peculiar ones.** Gaia proper motions are dominated by the association's own bulk motion plus Galactic rotation — Cyg OB2's systemic motion is **(−2.707, −4.317) mas/yr**, larger than a typical ejection signature — so tracing back absolute motions measures the common drift rather than the ejection. **Caught by the literature cross-match**, which is exactly what that gate exists for: the canonical ejected star **BD+43 3654** (O4If, ~70 M☉, Comerón & Pasquali 2007, Gaia DR3 2069819545390584192) sits 2.78° out, is present in the candidate list, and scored recovery probability **0.000**. With the systemic motion subtracted its peculiar motion is (+0.112, +5.046) mas/yr = **38.8 km/s**, giving 2.24° of travel over its 1.6 Myr age against the 2.78° observed. The first runaway result (260 raw, 109 corrected) is **withdrawn** | ~~blocking~~ | WP6 | **FIXED AND EXTERNALLY CONFIRMED 2026-07-28.** The traceback now subtracts the systemic motion, rotated to galactic at each star's own position; controls receive identical treatment. Re-running gives **119 raw recovered** (was 260), a **1.92% effective chance rate** measured per separation bin (was 4.5%), and **54.9 false-positive-corrected** (was 109.2). BD+43 3654 now scores **recovery probability 1.000** at **38.8 km/s** and **1.36 Myr** lookback, against a literature ~40 km/s and 1.6 Myr — agreement in both quantities from an estimator that knows nothing about the published values. The ledger's binned total and the aggregate correction now agree to 6% (58.2 vs 54.9); under the withdrawn version they disagreed by 37% (149 vs 109), which is weak independent support for the fix. Gate record: [wp6_runaway_crossmatch_execution.json](provenance/wp6_runaway_crossmatch_execution.json) · [wp6_runaway_crossmatch.py](scripts/wp6_runaway_crossmatch.py) · [wp6_runaway_crossmatch.csv](tables/wp6_runaway_crossmatch.csv). **Lesson recorded:** no internal check could have found this — it was visible only against a star whose answer is known independently, which is the whole argument for gating on external data |
| 17 | **The closure integral truncated at the very threshold the response smears across — and CUTS §16.2 already said so.** `wp6_closure_test.py` integrated the forward prediction from **8 M☉**, while the observed side counts every member's `P(M > 8)` **whatever that member's true mass is**. Stars with true mass below 8 M☉ whose *estimated* mass lands above it were counted as observed and never predicted, so the two sides were not the same quantity. The omitted up-scatter is large: R(estimated > 8 \| true M) is **0.23 at 7 M☉, 0.09 at 6, 0.03 at 5**. **The binding specification was already correct** — [CUTS §16.2](CUTS_AND_THRESHOLDS.md) writes the integral with *no lower limit*, and justifies forbidding a scalar completeness precisely because "the response also scatters mass estimates across the threshold", even reporting six cells where net up-scatter exceeds the recovery loss. The implementation truncated the effect the specification exists to capture. **Same defect class as issue #3, committed a second time at the other end of the same integral.** Why it survived: the *upper* limit absorbed all the scrutiny (it is physical, it drove the mass extension, it produced issue #14), the lower limit read as a restatement of the 8 M☉ supernova threshold rather than an integration bound, and the bias has the **same sign** as the excess WP6 was already reporting, so nothing downstream contradicted it. Found by probing a single node while the issue #15 injections ran | ~~blocking for WP6/WP7~~ | WP6→WP7 | **FIXED, PRE-REGISTERED, SCORED, ADOPTED.** Predictions declared before the grid run in [wp6_closure_floor_prereg.json](provenance/wp6_closure_floor_prereg.json); scored in [wp6_closure_floor_score_execution.json](provenance/wp6_closure_floor_score_execution.json). Floor **4.0 M☉**, chosen on convergence of the integral and not on its effect. **Regression check: at floor 8.0 the refactored estimator reproduces the published ratios exactly**, so the bound is the only behaviour change. **F1 PASS** — all 54 cells fell. **F2 PASS** — the reduction is larger at steeper α, as the IMF-weighting argument requires. **F3 FAIL** — the pre-declared 5-point spread limit was exceeded at 6.0 points; recorded as failed, not reinterpreted. What the data show is that F3's *substantive* claim holds: across all 54 cells the mean reduction is **A −18.3%, C −18.5%** despite turnoffs of 59.7 and 120.0 M☉, so the effect is flatly uncorrelated with the turnoff; the spread comes from **B** alone, whose effective completeness [CUTS §16.4](CUTS_AND_THRESHOLDS.md) already documents as differing from A's and C's by ~0.12. The threshold was set too tight for a quantity whose per-subgroup spread was already on the record. **F4 PASS** — grid median **1.444 → 1.105**. **Consequences: the published closure ratios 1.087 / 1.448 / 1.706 and the closing slope α = 2.070 are WITHDRAWN.** Corrected: **A 0.894, B 1.106, C 1.405**, closing **α = 2.230** [2.03–2.31], cells consistent with unity **2 → 6 of 54**, and **17/18** cells now inside the carried branch grid (was 12/18). The "shallower than Salpeter" reading is largely an artefact of this bug. Convergence is slower at grid level than the probe node suggested (1.4% from 4.0→3.0, a further 1.0% to 2.0); the floor was **not** moved in response — it was pre-declared — and the residual ~2.4% is carried as a systematic |
| 18 | **WP8 §5's ²⁶Al comparison denominator is wrong by ~100×.** It compares the ledger's SN-only ²⁶Al mass against "the ~1 M☉ inferred for the whole Cygnus complex" and reports the ordering as roughly 1000×. The **frozen WP1 marker's own measured complex-wide 1809 keV flux**, (3.9 ± 1.1) × 10⁻⁵ ph cm⁻² s⁻¹, corresponds at 1.6245 kpc to only **8.9 × 10⁻³ M☉** of ²⁶Al (3.4 × 10⁻³ – 1.3 × 10⁻² over Martin+10's quoted 1.0–2.0 kpc for the complex). The 1 M☉ figure appears to be Martin+10's **Galactic** stationary mass, 1.7–2.0 M☉, used for the complex. Found by WP11 Part B, which needed the same comparison and could not reproduce the margin | medium — corrects a published number, changes no ledger value | WP8→WP11 | **RECORDED, NOT REPAIRED — deliberately.** WP11 is one-way validation: `reports/wp8_crosschecks.md` and `provenance/wp8_crosschecks_execution.json` are left exactly as they were. The correct statement, in flux space where no mass conversion of the measurement is needed and the error cannot propagate, is that **this ledger's supernovae supply 10–42% of the measured complex-wide 1809 keV flux** — sub-dominant but not negligible, still a lower bound on the measurement since winds add to it, so **WP8's qualitative conclusion survives and its quantitative margin does not**. The manuscript quotes the flux-space form. Pre-registered prediction I1 was scored against the denominator its own pre-registration named (PASS) and is flagged **VACUOUS**, with the corrected re-score (**FAIL**) recorded beside it rather than either being suppressed. See [wp11_isotope_forecast.md §6](reports/wp11_isotope_forecast.md) |
| 12 | **CLOSED 2026-07-28 — cause confirmed and fixed in `repair_v5`.** The extinction-scale error was real and was the cause: correcting the anchor prior's mean (kriging) flattened B's tilt, \|T\| 2.40 → **1.31** (trend p 0.017 → 0.186), with bin-2 residual 2.51 → **0.36** and the high-mass deficit −0.79/−0.54 → −0.03/+0.17. Pre-registered and confirmed ([wp3_kriging_prior_outcome.json](provenance/wp3_kriging_prior_outcome.json)). Original statement: **an independently measured extinction-scale error in B.** The near-IR cross-check (issue #1d) finds B's optical A_V high by +0.385 mag relative to A and C, which rescales B's whole mass ladder by **12.6%**; by contrast the disputed age shift moves assigned mass by only 0.4%, because it mainly relocates the PMS/Henyey fold rather than rescaling masses. That is precisely finding F1's decomposition of B's anomaly into *a localized bump plus a global mass-scale displacement*: `repair_v4` removed the bump (age/fold), and the displacement is the extinction offset. The near-IR reddening-law **shape** is the same for all three subgroups (E(J−H)/E(H−Ks) = 1.21–1.28), so B does not need its own R_V — this disfavours step 4b and favours a scale error. **Not yet demonstrated:** that correcting the extinction actually flattens the tilt; the injections pass through the same biased estimator so the bias partly cancels, and the sign of the net effect must be established by running the chain, not by argument. Original statement: **CygOB2-B's mass function still has a mass-dependent tilt after the age fix.** Under the better-powered replacement trend statistic, B's `repair_v4` baseline residuals (+1.20, +2.08, +2.51, −0.24, −0.79, −0.54) give trend p = **0.017** (T = −2.40), against A −0.39 and C +0.19 which are flat — so B's baseline **fails** the gate under that statistic even at `repair_v4`. The incumbent rank test scored it 0.111 only because the sequence rises then falls instead of increasing monotonically. This is the *global mass-scale displacement* half of finding F1 in [wp5_cygob2b_age_caustic_fix_brief.md](tasks/wp5_cygob2b_age_caustic_fix_brief.md): the age fix removed the localized bump, but the displacement was never addressed | ~~blocking~~ | WP5←WP3 | **CLOSED 2026-07-28 by `repair_v5`** (row header above). Text below is the state while it was open — leading hypothesis was issue #1d: with no spectroscopic anchors near B, its A_V is set by broadband photometry alone, so an extinction-scale error moves B's whole mass scale coherently. Step 4b of the gated plan (continuous R_V / extinction-scale scan for B, to be justified by **independent extinction evidence, never gate score**) is the pre-written next experiment |

---

## 10. Next actions

*Updated 2026-07-28 after WP6 ran. WP5 is accepted; WP6's closure test,
attribution, runaway search and ledger are complete.*

1. ~~**The one open WP6 item: the literature runaway cross-match.**~~
   **DONE 2026-07-28 — and it caught issue #16.** *(Historical text: the plan's
   gate requires recovering known Cygnus runaways (bow-shock stars). The anchor
   table carries no runaway flag and this environment has no network access, so
   it is not done. Until it is, the 109 recovered runaways are a statistical
   result with no per-star confirmation.)* The cross-match was run:
   BD+43 3654 recovers at p = 1.000, 38.8 km/s, 1.36 Myr against a literature
   ~40 km/s and 1.6 Myr, and it was this check that exposed the absolute- vs
   peculiar-proper-motion defect. The 109-runaway figure quoted above is the
   **withdrawn** one; the corrected count is 119 raw → **54.9** false-positive
   corrected. See issue #16 and
   [wp6_runaway_crossmatch_execution.json](provenance/wp6_runaway_crossmatch_execution.json).
2. **WP7 must not refit α.** The closure test is the analysis's only
   out-of-sample check — `k` is fitted from 2–8 M☉ counts and the >8 M☉ census
   never enters the WP5 likelihood. Carry **α = 2.6 as disfavoured** by the
   census and quote N_SN per branch. Converting the test into a fitted slope
   would destroy it.
3. **No extra α branch is needed any more.** Before issue #17 the closing slope
   was 1.92–2.18 and 6 of 18 cells extrapolated below α = 2.0, which argued for
   a dedicated α ≈ 2.0–2.1 run. With the integration floor corrected the closing
   slope is **2.03–2.31 and 17 of 18 cells sit inside the carried grid**, so the
   answer is measured rather than extrapolated. That recommendation is
   withdrawn — a worked example of a bug fix removing a planned experiment
   rather than adding one.
4. **WP7 inputs, ready now**: [wp6_massive_census.cat](tables/wp6_massive_census.cat)
   (380.6 living stars above 8 M☉ across three provenance-flagged channels), the
   per-subgroup turnoffs, and the runaway correction with its two totals (54.9
   unclipped, 58.2 binned) to carry as a systematic. These supersede the
   withdrawn absolute-PM figures (471.9 / 109 / 149) — see issue #16.
5. **Runaway count is a lower bound, twice over** — the ±8° box caps recovery at
   **44 km/s over 5 Myr** (issue #10, now quantified at the WP2 distance of
   1.62 kpc, superseding CUTS §4.2's 1.4 kpc figures), and the traceback is 2D
   so line-of-sight ejections are missed.
6. **Discharge the remaining WP3 adoption obligations** from
   [wp3_kriging_adoption.json](provenance/wp3_kriging_adoption.json): **O1**
   report the revised star-formation history (A 4.00, B 4.07, C 2.52 Myr — two
   older subgroups plus one younger); **O2** the anchor absolute-scale
   systematic; **O3** CygOB2-B's 4-anchor calibration asymmetry. O4 is
   discharged.
   *Definition attached 2026-07-30 (item B2): the O1 ages are the **WP5 fitted
   truth-age posterior means** (counts-based), not the WP4 upper-MS MAPs of
   3.981 / 3.548 / 2.512 Myr, which are their prior. On the accepted
   `repair_v7` chain they read **A 4.00, B 4.09, C 2.52**. See
   [wp4_wp5_age_reconciliation.md](reports/wp4_wp5_age_reconciliation.md).
   **O1, O2 and O3 are discharged in the manuscript** — see
   [wp3_obligations_discharge.md](reports/wp3_obligations_discharge.md).*
7. **Binary mass transfer is now bounded, not merely acknowledged (item T3,
2026-07-30).** Two independent lines. *Empirical*: Härer's Fig. 2 is a BPASS
calculation **with** binaries for a Cyg OB2-like population; digitized exactly
from the local PDF (its step edges land on BPASS's native 0.1-dex log-age grid),
it gives 6.00 SNe/Myr per 10⁴ M☉ at 4 Myr against our single-star 9.28–9.45 at
matched α, age and mass — a **ratio of 0.72 after correcting for BPASS's wider
IMF range**, i.e. the binary-inclusive rate sits *below* ours. *Theoretical*:
Zapartas+2017 raise the integrated CCSN count by 14⁺¹⁵₋₁₄% but attribute it to
50–200 Myr events from 4–8 M☉ binaries, and find the early delay-time
distributions "remarkably similar", diverging only at ~20 Myr — a channel that
cannot have operated in a 4 Myr association. **Adopted bracket ±30% on N_SN**,
run through WP7's engine: baseline **5.90 / 8.43 / 10.96**, a span of 5.1
supernovae against the headline branch span of 23.1. **The unmodelled systematic
is ~5× smaller than the branch spread already reported.** Discussion-section
bound, never marginalized in —
[wp7_binary_bound.md](reports/wp7_binary_bound.md) ·
[wp7_binary_bound_execution.json](provenance/wp7_binary_bound_execution.json).

**Carried methodological finding for the systematics section.** At
   `N_INJECT_PER_MASS = 400` a single-node truth model leaves 41–44% of the
   54-cell grid with an indeterminate verdict; `repair_v6`'s nine-node mixture
   cuts this to **17%**. Future branch-grid comparisons must report the
   indeterminate set.

## 10b. External cross-checks

Published measurements the chain has been compared against, in
[cross_checks/](cross_checks/). Machine-readable record with input hashes:
[wp6_external_crosschecks_execution.json](provenance/wp6_external_crosschecks_execution.json).

| source | quantity | ours | theirs | verdict |
|---|---|---|---|---|
| [Härer+2025](cross_checks/harer_2025_supernova_rate.md) (A&A 703, A111) | SN rate | 0.78 ev/100 kyr (α=2.3)<br>1.92 (α=2.0) | 0.25–2.0 range;<br>"several per Myr" | **agree** |
| [Härer+2025](cross_checks/harer_2025_supernova_rate.md) / Wright+2015 | association mass (**like-for-like**: primaries, 0.08–120 M☉) | **2.42 × 10⁴ M☉** | 1.65 × 10⁴ M☉ | **agree at the factor level, 1.47×** |
| [Orellana+2021](cross_checks/orellana_2021_astrometry.md) (MNRAS 502, 6080) | systemic μ_α* | **−2.7067** | −2.71 ± 0.02 | **agree, 0.003 mas/yr** |
| [Orellana+2021](cross_checks/orellana_2021_astrometry.md) | distance | 1616–1629 pc | 1669 ± 6 pc | 3.3% offset → **+2.4% systematic** |
| [Berlanas+2019](cross_checks/berlanas_2019_two_distance.md) (MNRAS 484, 1838) | two-distance split | one population, 45 pc depth | 1350 + 1760 pc | **not confirmed — test cannot settle it** |

**These are validations, never calibrations.** No value in the chain was tuned
toward any of them; the Härer comparison was not made until after WP6 closed. A
disagreement becomes an issue in §9, not a reason to move a number.

**Mass row corrected 2026-07-30 (item B1).** The row previously read
"1.74 × 10⁴ M☉ … agree, 5%". That compared our 0.5–120 M☉ **primary** integral
against a literature mass drawn over the **whole** IMF, so it paired mismatched
definitions and the 5% was a cancellation. Three nested integrals of the same
`k` were all being called "the association mass": primaries 0.5–120 M☉
(1.75 × 10⁴), primaries 0.08–120 M☉ (**2.42 × 10⁴**, the like-for-like one), and
that plus unresolved companions (2.92 × 10⁴, WP5's headline). Decomposition,
reproduction check and the Wright+15 definitional argument:
[wp5_association_mass_reconciliation.md](reports/wp5_association_mass_reconciliation.md)
· [wp5_association_mass_reconciliation_execution.json](provenance/wp5_association_mass_reconciliation_execution.json).
**No fit was re-run and no stored number moved**; the *rate* half of the Härer
cross-check — the load-bearing half, since it never integrates below the
turnoff — is untouched.

Three points worth carrying into WP7:

1. **The Härer agreement is the closest thing to an end-to-end check the project
   has.** Rate and mass depend on `k` differently — mass integrates the whole
   IMF, rate samples it only at the turnoff — so a broken extrapolation above
   8 M☉ would generally break them inconsistently. Both agree.
2. **Our N_SN is a lower bound relative to BPASS.** Their Fig. 2 is dominated by
   type Ic at 3–5 Myr, the stripped-envelope channel, which requires a binary.
   Our turnoff counting is single-star. Same physics as issue #15, arriving
   independently.
3. **The distance offset is not adopted**, but ~2.4% on every closure ratio
   belongs in the systematics budget — it helps CygOB2-A and hurts CygOB2-C.

## 10c. Open recommendations (not executed)

| item | why it is open | cost | doc |
|---|---|---|---|
| **`repair_v7`** — mass-dependent f_bin below 8 M☉ | **DONE 2026-07-29.** Justified by the discriminator (**D2 = +9.87% ± 1.58%**, five times the 2% threshold, against **D1 = +0.16%**), then executed over the full mass range — WP5 base grid and WP6 mass extension, 162 nodes each. **G3 re-passed**, `A_or_C_regressions: []`, 40/54 branches, V1 still 5/5 byte-identical. `k` **+0.54% ± 0.25%**, matching the discriminator's out-of-sample **+0.54% ± 0.19%** exactly. Closure grid median **1.105 → 1.067**, closing α **2.230 → 2.254**, cells inside the branch grid **17/18 → 18/18**, cells consistent with unity **6/54 → 9/54**. Issue #15's 0.0038 correction is now **internal and must not be applied again** ([CUTS §18.9](CUTS_AND_THRESHOLDS.md)) | ~6 h actual (WP3/WP4 stages were unnecessary) | [repair_v7_recommendation.md](tasks/repair_v7_recommendation.md) · [wp5_fbin_discriminator_execution.json](provenance/wp5_fbin_discriminator_execution.json) · [wp6_closure_test_execution_repair_v7.json](provenance/wp6_closure_test_execution_repair_v7.json) · [wp6_multiplicity_subsumption.json](provenance/wp6_multiplicity_subsumption.json) |
| **parallax-blind membership test** | issue #8's verdict is partly circular: parallax is a WP2 clustering feature | moderate | [berlanas_2019_two_distance.md](cross_checks/berlanas_2019_two_distance.md) |
| **distance-contamination test for CygOB2-C** | inflates masses and therefore the closure ratio — C's direction; covered by no WP6 alternative | moderate | as above |

**~~Recommendation on `repair_v7`: defer.~~ OVERTURNED 2026-07-29 by its own
test.** The recommendation was to defer, reasoning that the sub-8 M☉ effect
should be *smaller* than the 0.4% measured above 8 M☉ because the f_bin change
there is smaller. That reasoning was wrong: it treated the mechanism as a
smooth response shift when it is a **threshold** effect. Inside the calibration
window most stars are recovered either way, so recovery is insensitive
(**D1 = +0.16%**); at the 8 M☉ boundary a small brightness shift converts
directly into a large crossing probability (**D2 = +9.87%**). The 25× ratio
between them is exactly what prediction G2 anticipated.

**`repair_v7` is justified and should be scheduled.** Three things follow:

1. **The 15-minute test paid for itself.** It cost ~3 minutes of compute and
   reversed an 8-hour decision that would otherwise have been made on intuition.
   Pre-registering a cheap discriminator before an expensive re-run is now a
   pattern worth repeating.
2. **Do not take the shortcut D1 appears to offer.** D1 ≈ 0 suggests only WP6
   needs re-running, but D1 measures the *recovery fraction*, not *mass
   migration within the window*, and the WP5 likelihood consumes the full
   response matrix. D2 is direct evidence that mass estimates shift enough to
   move stars across a threshold; nothing makes that stop at 8 M☉. Scope is the
   full chain.
3. **It does not block WP7 — now measured, not assumed**
   ([wp5_fbin_k_impact_execution.json](provenance/wp5_fbin_k_impact_execution.json)).
   D1 was a recovery fraction and did not capture mass migration, so the
   k-relevant quantity was computed directly: P(estimated mass lands inside the
   2–8 M☉ window | true M), IMF-weighted. Result: **k and therefore N_SN shift
   by +0.54% ± 0.19%**. Against a Class E branch spread on N_SN of **2.33 to
   27.42, a factor of 11.8**, the correction is ~2000× smaller than the
   uncertainty WP7 must report regardless. **WP4 ages, per-star masses, the
   runaway correction and the lifetimes are bit-identical** — the injection
   truth model plays no part in fitting real stars. Of WP7's four inputs, only
   `k` moves, by 0.54%.

   This is categorically different from issues #16 and #17, which were defects
   of **unknown** size that turned out to be large. Proceeding on a bounded,
   quantified perturbation is legitimate; proceeding on an unknown one is what
   produced the two withdrawn results. WP7 results remain **provisional** until
   re-checked against `repair_v7`, which is cheap because WP7 is pure
   computation on frozen inputs.
4. **`repair_v7` is a WP6 fix, not a WP7 fix.** It moves WP6's closure ratio
   1.099 → ~1.074, removing ~23% of the 9.9% still unexplained. **The same
   physical error is 0.54% on a window-integrated normalization and 9.87% on a
   threshold-crossing probability** — an 18× asymmetry produced purely by one
   quantity being an average and the other a boundary.
5. **Cost estimate corrected: ~6 h, not 8.** The f_bin change lives entirely in
   the WP5 injection truth model, so the WP3 and WP4 stages are unnecessary. The
   `repair_v7` label is a misnomer for what is really a WP5/WP6 re-run.

**The CygOB2-C argument still stands** and was not the part that was wrong: a
uniform truth-model change moves all subgroups the same way, so `repair_v7`
cannot explain C's *direction* disagreement with A. That remains the open
science question, and `repair_v7` is not its answer.

## 10d. WP11 Part B — the ²⁶Al / ⁶⁰Fe forward prediction (2026-07-31)

| | |
|---|---|
| **Verdict** | **COMPLETE** — gate met; I2–I5 PASS; **I1 PASS but vacuous, corrected form FAILS** (recorded, not reinterpreted) |
| **Report** | [wp11_isotope_forecast.md](reports/wp11_isotope_forecast.md) |
| **Pre-registration** | [wp11_isotope_prereg.json](provenance/wp11_isotope_prereg.json) — written and run *before* any flux existed |
| **Execution record** | [wp11_isotope_forecast_execution.json](provenance/wp11_isotope_forecast_execution.json) |
| **Scripts** | [wp11_isotope_prereg.py](scripts/wp11_isotope_prereg.py) · [wp11_isotope_forecast.py](scripts/wp11_isotope_forecast.py) |
| **Tables** | [wp11_isotope_forecast.csv](tables/wp11_isotope_forecast.csv) (per branch × arm) · [wp11_isotope_summary.csv](tables/wp11_isotope_summary.csv) |
| **Brief** | [wp11_bowshock_isotope_brief.md §4](tasks/wp11_bowshock_isotope_brief.md) |

**Headline** — ⁶⁰Fe splits the headline branch set exactly along α: **18/18
α = 2.0 branches clear COSI's 3σ two-year narrow-line sensitivity, 0/18 α = 2.3
branches do.** COSI measures the axis the verdict hinges on and DR4 cannot fix.

**Two departures from the brief, both deliberate and both pre-registered.**

1. **The brief's estimator was not used.** Its `rate × mean lifetime × yield` is
   the steady-state limit, and Cyg OB2 is not in steady state — only 0.16–0.31
   of the saturated ⁶⁰Fe mass has accumulated (first SN ~1.4 Myr ago against
   τ = 3.78 Myr), so it would have overstated ⁶⁰Fe by 3–6×, in the direction
   that manufactures detectability. Replaced by a decay-weighted, mass-resolved
   sum that reduces to the brief's formula in the limit the brief assumes.
2. **The yield branch is three named published arms, not a ± factor.** The brief
   asked for "a declared literature range". The two non-null arms come from the
   *same paper and differ only in the WR mass-loss prescription*, so the 30×
   spread is Limongi & Chieffi's own stated systematic rather than one this
   project invented. The third arm is LC18-Recommended, which predicts
   **identically zero** because it collapses everything above 25 M☉ — verified,
   not assumed: 0 of 315M sampled supernovae fall below that mass.

**Finding T4 — a WP8 number does not survive contact with WP1.** WP8 §5 compares
against "~1 M☉ for the whole Cygnus complex"; the frozen WP1 flux gives
**8.9 × 10⁻³ M☉** at 1.62 kpc, ~100× less (it appears to be Martin+10's
*Galactic* stationary mass). Supernovae therefore supply **10–42%** of the
measured complex-wide ²⁶Al flux, not ~0.1%. WP8's conclusion — expected
ordering, not a tension — survives; its margin does not. **Nothing upstream was
retuned** (one-way validation): `reports/wp8_crosschecks.md` and
`provenance/wp8_crosschecks_execution.json` are untouched, and the manuscript
quotes the flux-space comparison, which needs no mass conversion and so cannot
inherit the error. Added to the open-issues register as the one WP11-discovered
defect.

**Gate G11d verified mechanically.** Diffing the `audit.py` inventory before and
after, the only pre-existing files changed are `scripts/wp10_inputs.py`,
`scripts/wp10_numbers.py` and their regenerated provenance — i.e. exactly the
manuscript machinery the brief §5 requires touching. No WP1–WP9 data product,
table, report or provenance record moved.

## 11. Conventions used throughout

- **Nothing is overwritten.** Failed runs and superseded artifacts are kept with
  `_failed_<date>` or `_repair_v<n>` suffixes and stay referenced from the
  provenance logs.
- **Every gate verdict is a JSON record**, not just prose — `provenance/*_gate.json`,
  `*_validation.json`, `*_manifest.json` carry sha256 digests of every input and
  output.
- **Model branches are carried, never averaged** — isochrones {PARSEC, MIST} ×
  R_V {3.0, 3.1, 3.5} × IMF slope {2.0, 2.3, 2.6}, per plan §1.4.
- **`P` is overloaded** — membership probability (WP2, per-star, [0,1]),
  p-values (gate statistics, low = bad), response probabilities (WP5 injections),
  and pulsar spin period (WP1/WP8). See [GLOSSARY.md](GLOSSARY.md).
