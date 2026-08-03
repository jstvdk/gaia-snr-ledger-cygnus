# WP11 — Post-hoc validation extensions: IR bow-shock runaway confirmation and radioisotope forecast

*Status: **PART B ADOPTED — execute before submission** (PI decision,
2026-07-30). **Part A is NOT adopted for Paper 1**: it is held as a
ready-to-execute reserve for the referee stage (executable in ~1 day
mid-revision if a referee asks about runaway confirmation) and is otherwise
deferred to the DR4 rerun, where radial velocities make the traceback 3D and
roughly double the orientation test's value. An executing agent should
therefore do §4 and §5–§7 as they apply to Part B, and leave §3 untouched.*

*Why Part B was adopted: beyond falsifiability, it opens a citation channel
into the MeV-astronomy community — a per-branch ⁶⁰Fe/²⁶Al line-flux prediction
for the Cygnus region is directly usable by COSI target-selection and
sensitivity studies, and either COSI outcome serves the paper (a detection
validates the census method; an upper limit below the α = 2.0 band disfavours
exactly the branch the verdict hinges on).*

*Created 2026-07-30 on PI decision, after triage of a supervisor suggestion
("radioactive-element maps, WISE/MSX infrared emission"). Two of the four
candidate additions survived triage; the other two are **rejected for Paper 1**
with reasons recorded in §2. This brief is self-contained for a fresh agent: it
lists inputs, procedure, pre-registration requirements, gates and manuscript
integration. It consumes the frozen `repair_v7` chain and MUST NOT modify
anything upstream.*

**Governing documents:** [paper1_execution_plan.md](../paper1_execution_plan.md)
(conventions §1 are binding here) · [PROJECT_TRACE.md](../PROJECT_TRACE.md) ·
[CUTS_AND_THRESHOLDS.md](../CUTS_AND_THRESHOLDS.md) ·
manuscript rules in [manuscript/README.md](../manuscript/README.md)
(**never type a number into `main.tex`** — every new quantity becomes a macro in
`scripts/wp10_numbers.py`).

---

## 1. Objective

Two bounded additions to the paper's validation layer:

- **Part A — IR bow-shock cross-match.** Convert the runaway census from a
  statistical result (119 raw → 54.9 corrected, one per-star confirmation) into
  a per-star-confirmed one, by cross-matching the WP6 candidates against
  published infrared bow-shock catalogues. Cygnus X is the best-surveyed
  bow-shock field on the sky; the original WP6 gate asked for exactly this
  cross-match and it was recorded as not done for lack of network access
  ([PROJECT_TRACE §10 item 1](../PROJECT_TRACE.md)).
- **Part B — radioisotope forward prediction.** Convert the ²⁶Al consistency
  band into a falsifiable forecast: compute the ledger's predicted SN-only
  ²⁶Al and ⁶⁰Fe masses and line fluxes per headline branch, and compare with
  INTEGRAL/SPI sensitivity and the expected COSI (launch ~2027) line
  sensitivity. ⁶⁰Fe is the supernova-specific tracer (WR winds do not produce
  it), so this is the ledger's cleanest future observable.

Neither part changes any published number. Both are **one-way validation**: a
tension is reported, nothing upstream is retuned.

## 2. Rejected scope — binding, do not re-open for Paper 1

| candidate | verdict | reason |
|---|---|---|
| Diffuse WISE/MSX emission / cavity morphology | **Paper 2** | ISM/γ-ray side; out of scope per the mission statement; would undo the 2026-07-30 measurement-first reframe |
| WISE/MSX search for extinction-hidden OB stars | **rejected** | The census closure test already answers it: observed ≥ predicted (6.7% excess at Salpeter), so no large hidden population exists inside the footprint. Cover with one sentence citing the closure result |
| Re-deriving a ²⁶Al map | **rejected** | SPI/COMPTEL ~3° resolution cannot separate Cyg OB2 from the complex; the WP8 band already uses the flux correctly; inversion is forbidden by the WP8 pre-registration |

## 3. Part A — IR bow-shock cross-match

### 3.1 Inputs

| input | where | note |
|---|---|---|
| WP6 runaway candidates | `tables/wp6_runaways.cat` (repair_v7 chain) | 119 raw candidates with traceback parameters and peculiar-PM vectors |
| Chance-rate machinery | WP6 control-field pattern, [wp6_runaway_crossmatch.py](../scripts/wp6_runaway_crossmatch.py) | reuse, do not reinvent |
| Kobulnicky+2016 catalogue | VizieR `J/ApJS/227/18` (ApJS 227, 18) | ~709 IR bow-shock candidates, 24/22 µm, tabulates position angles |
| E-BOSS I + II | Peri+2012 (A&A 538, A108) and Peri+2015 (A&A 578, A45), VizieR | independent, partially overlapping |
| Orientation caveat literature | Kobulnicky+2019 | bow-shock apex vs. velocity misalignment statistics — needed for §3.4 framing |

**Network dependency:** catalogue acquisition needs VizieR access. If the
execution environment has none, stop and hand the exact VizieR identifiers
above to the PI to fetch; do not substitute from memory.

### 3.2 Procedure

1. **Freeze the catalogues** under WP1 conventions: store under
   `data/processed/wp11_*`, write a manifest with sha256 digests and an
   execution JSON under `provenance/`.
2. **Pre-register before any match is scored** —
   `provenance/wp11_bowshock_prereg.json`, following the WP7/WP8/WP9 pattern.
   Required predictions (executing agent freezes exact thresholds before
   running; suggested content):
   - **BS1 (positive control):** BD+43 3654 appears in the catalogues and
     matches its WP6 candidate entry. If BS1 fails the pipeline is broken —
     fix and re-run; BS1 is a control, **not** a discovery, and must be labeled
     so in the paper.
   - **BS2 (headline):** the number of bow-shock matches among the 119
     candidates exceeds the chance-alignment expectation, with the chance rate
     measured per separation bin from control stars exactly as WP6 measured its
     false-positive rate.
   - **BS3 (secondary):** among matches with a tabulated apex position angle,
     the apex direction is closer to the candidate's peculiar-PM direction than
     uniform (pre-declare the circular statistic and threshold). Kobulnicky+2019
     show apex–velocity misalignment is common, so BS3 failing is reportable,
     not gate-level.
3. **Cross-match** by position (declare the radius from the catalogues' own
   astrometric quality before matching; account for the epoch difference with
   the candidates' proper motions if the radius makes it non-negligible).
4. **Score the pre-registered predictions.** Failures are recorded as failed,
   never reinterpreted — project rule.
5. **Report** `reports/wp11_bowshock_crossmatch.md`, matches table
   `tables/wp11_bowshock_matches.csv`, execution JSON with input hashes.

### 3.3 What the result feeds

- The census/runaway section of the manuscript (per-star confirmations
  alongside BD+43 3654).
- The **C4 in-situ bound** discussion: WP9 devil's-advocate objection 2
  conceded that C4 is an optimistic upper bound; per-star confirmations with
  directions are the strongest DR3-era evidence available on it. Do **not**
  recompute C4 — the bound stands; the matches are corroboration, and any
  tension is reported as tension.

### 3.4 Binding one-sidedness caveat

Bow shocks form only where the ISM density and the stellar velocity suffice;
the majority of true runaways show none. **A non-match is therefore not
evidence against runaway status**, and no candidate may be removed or
down-weighted for lacking a bow shock. The test confirms; it never refutes.

## 4. Part B — ²⁶Al / ⁶⁰Fe forward prediction

### 4.1 Inputs

| input | where | note |
|---|---|---|
| Per-branch SN rates and epochs | `tables/wp7_ledger.csv` (headline set, 36 branches) | frozen |
| Distance | the WP3 value (1.62 kpc) via its versioned artifact | use the artifact, not the number |
| Baseline SN-only ²⁶Al estimate | [wp8_crosschecks.md §5](../reports/wp8_crosschecks.md) | extend from baseline-only to per-branch |
| Per-SN yields | Limongi & Chieffi 2018, Sukhbold+2016 (²⁶Al, ⁶⁰Fe) | see yield branch below |
| Instrument sensitivities | SPI (Wang+2020 ⁶⁰Fe); COSI expected narrow-line sensitivity (Tomsick+2023) | literature values, cited |

### 4.2 Procedure

1. **Fix the yield branch before computing any flux** (this is the whole
   pre-registration for Part B, recorded in
   `provenance/wp11_isotope_prereg.json`): per-SN ⁶⁰Fe and ²⁶Al yields for
   >30 M☉ progenitors are strongly model-dependent — carry a declared
   literature range as a branch, never a single value, and do not adjust it
   after seeing the fluxes.
2. Per headline branch: steady-state isotope mass = rate × mean lifetime ×
   yield; convert to 1.809 MeV (²⁶Al) and 1.173/1.333 MeV (⁶⁰Fe) line fluxes
   at the WP3 distance.
3. Compare with SPI limits and COSI expected sensitivity; state the
   detectability verdict per branch, including the honest outcome if that
   verdict is "below reach" — a stated non-detectability is still a forecast.
4. **The WP8 prohibition stands:** the measured complex-wide ²⁶Al flux is never
   inverted into a supernova count. The prediction runs forward only, and the
   SN-only component must be stated as a lower bound on the complex-wide signal
   (WR winds dominate ²⁶Al at this age; they contribute no ⁶⁰Fe).
5. **Report** `reports/wp11_isotope_forecast.md` + execution JSON.

## 5. Manuscript integration (both parts)

1. Every quoted quantity enters as a macro via `scripts/wp10_numbers.py`,
   reading the WP11 artifacts through `scripts/wp10_inputs.py` (register them
   as authorized inputs). **No hand-typed numbers.**
2. Part A lands in the census-closure/runaway results text and the systematics
   paragraph on the in-situ bound; Part B lands as one Discussion paragraph
   adjacent to the DR4 outlook (a second forward-looking instrument, COSI,
   beside DR4) and one sentence in the ²⁶Al cross-check paragraph.
3. **Mandatory disclosure:** unlike the WP8 markers, the WP11 comparisons were
   **not frozen at WP1** — they are post-hoc additions, pre-registered before
   scoring but chosen after the ledger existed. The manuscript must say this in
   the same sentence that introduces them. The credibility of the WP8 layer
   rests on the freeze; do not blur the two.
4. Re-run `scripts/wp10_validate.py`; all seven checks must pass.

## 6. Gate (conjunctive)

| criterion | requirement |
|---|---|
| G11a | Both pre-registrations exist with sha256-hashed inputs **before** any score/flux was computed |
| G11b | BS1 positive control passes (or the failure is diagnosed and the fix documented before re-run) |
| G11c | Chance-alignment control measured, not assumed, per the WP6 pattern |
| G11d | No upstream artifact modified — verify with the [audit.py](../audit.py) inventory before and after |
| G11e | Yield branch declared before fluxes; ²⁶Al never inverted |
| G11f | Manuscript integration passes `wp10_validate.py`; post-hoc disclosure sentence present |
| G11g | Failed predictions recorded as failed |

## 7. Risks and cost

| risk | mitigation |
|---|---|
| Crowded-field chance alignments inflate matches | BS2's measured chance rate; separation-binned, as WP6 |
| Apex–velocity misalignment reads as failure | BS3 declared secondary; Kobulnicky+2019 cited in the prereg |
| Yield model spread swamps the ⁶⁰Fe forecast | that spread *is* the result — report the branch range |
| Scope creep back toward diffuse-emission analysis | §2 is binding; log temptations for Paper 2 |
| No network access | stop at §3.1; hand VizieR IDs to the PI |

**Cost:** Part A ~0.5–1 day (dominated by catalogue acquisition and the
controls); Part B ~0.5 day (pure computation on frozen inputs). Neither blocks
the pre-submission items in [manuscript/README.md](../manuscript/README.md);
if time is short, Part B alone is still worth having, and Part A without Part B
is too — they are independent.
