# Issue #19 — completion report: repair_v8

*2026-10-01. Brief: [tasks/issue19_stale_wp4_inputs_brief.md](../tasks/issue19_stale_wp4_inputs_brief.md).
Pre-registration (written before any rerun): [provenance/issue19_repair_v8_prereg.json](../provenance/issue19_repair_v8_prereg.json).
Integrity: [issue19_repair_v8_integrity.json](../provenance/issue19_repair_v8_integrity.json).
Outcome: [issue19_repair_v8_outcome.json](../provenance/issue19_repair_v8_outcome.json) ·
[tables/issue19_before_after.csv](../tables/issue19_before_after.csv).
Chain: [scripts/run_repair_v8_chain.sh](../scripts/run_repair_v8_chain.sh).*

## 1. Verdict

**repair_v8 is adopted.** All four integrity checks passed, and that is the only
thing adoption depended on. **4 of 6 predictions passed; P1 and P4 failed, both
only at R_V = 3.5.** The headline count moves from 8.43 to 8.36. The closure
narrative changes for subgroup A: its baseline census now closes at 1.006,
against 0.865 before. WP5 gate G3 **fails** on repair_v8, by its no-regression
clause; this is reported, not repaired.

## 2. What was done

| step | what | product |
|---|---|---|
| 1 | Pre-registration: P1–P6, I1–I4, adoption rule, decisions D1–D3; 12 inputs, 648 reused node files and 49 comparison artifacts hashed | `provenance/issue19_repair_v8_prereg.json` |
| 2 | `wp4_anchors_hrd.py` requires explicit `--age-version`/`--extinction-version`/`--output-version` and refuses the unversioned posterior. One mass per family × R_V, each at that branch's repair_v5 age and its own M_G0. Present-day `Mass` is kept (D1). | `wp4_anchor_hrd_repair_v8.parquet` |
| 3 | `wp4_mass_posteriors_repair.py` reads the versioned anchors (`--anchor-version`) and uses each branch's anchor mass. Label `photometric_posterior_repair_v1` → `photometric_posterior`. | `wp4_mass_posteriors_repair_v8.parquet` + samples |
| 4 | WP5 joint fit (`--mass-version` / `--age-version` / `--response-version` split out of `--upstream-version`). The repair_v7 node responses are reused: the injection code never opens the mass posteriors. | `wp5_*_repair_v8` |
| 5–7 | WP6 census/closure/attribution/living ledger, WP7 ledger (convergence scan to 2×10⁶), α-headline sets, binary bound, WP8, WP9, WP11 (5×10⁵, as published), reconciliations, α-plausibility, WP12 | `_repair_v8` siblings via `scripts/chain.py` |
| 8 | `chain.py` (`ADOPTED = repair_v8`); `wp10_inputs` resolves through it and forbids the unversioned WP4 posterior, anchor file and two pre-repair WP4 tables, plus every superseded repair_v7 product; `\ageEnvLo`/`\ageEnvHi`/`\ageEnvRows` macros; 2.25/5.67 removed from the V7 whitelist | numbers.tex, 216 macros |
| 9 | fig04/fig09 envelope computed and relabelled "retained upper-MS envelope" | `figures/paper/` |
| 10 | Text: main.tex (envelope ×3, PMS paragraph, α-2.6 sentence), WP4 report banners, reconciliation §5, O1, PROJECT_TRACE WP4 headline + issue #9, brief references | |

Nothing at repair_v7 or earlier was overwritten (I4). WP12 on repair_v8 checks
its inputs against a **new** hash record,
[`wp12_revision_prereg_repair_v8.json`](../provenance/wp12_revision_prereg_repair_v8.json),
which carries R1–R7 verbatim from the original record. The original was not
edited.

## 3. Integrity checks

| | check | result |
|---|---|---|
| I1 | 1,242 photometric rows and their sample cube identical to repair_v5 on all six branches | **PASS** |
| I2 | R_V = 3.1 anchor masses reproduce the scope diagnostic | **PASS**. The first run failed at 1.8×10⁻¹⁵ M☉, a CSV round-trip artifact; tolerance set to 10⁻¹² and recorded. |
| I3 | Replay (repair_v5 masses + repair_v7 responses through the modified fit) reproduces repair_v7: all 54 cells, all draws | **PASS**. The first run reported FAIL from a NaN-comparison bug in the check (an all-NaN column); no value differed. |
| I4 | All 697 hashed repair_v7 artifacts unchanged after the full chain | **PASS** |

## 4. Predictions

| | prediction | measured | outcome |
|---|---|---|---|
| P1 | A-PARSEC window loses 8.0 ± 0.5 stars on every PARSEC branch | 8.0 / 8.0 / **11.0** at R_V 3.0 / 3.1 / 3.5 | **FAIL** |
| P2 | MIST, and B and C everywhere: k and N_death within 1 % | worst k −0.49 % (C, MIST, R_V 3.5); 0 of 162 ledger rows outside | PASS |
| P3 | Baseline N_death ∈ [8.20, 8.43] | 8.360 | PASS |
| P4 | A's 9 PARSEC closure ratios rise 8–20 % | 16.2–16.7 % at R_V 3.0/3.1; **20.7–20.9 %** at R_V 3.5 (6/9 in range; median 16.7 %) | **FAIL** |
| P5 | Association closing slope ∈ [2.10, 2.30] | 2.186 | PASS |
| P6 | repair_v5 envelope 2.00–4.01 Myr, no PMS row | 2.00–4.01, 0 PMS, 66 rows | PASS |

**What the failures mean.** The disclosed prior knowledge said that the
R_V = 3.0/3.5 parts of P1 were genuine forecasts. At R_V = 3.0 the forecast
held. At R_V = 3.5 the per-R_V magnitudes (F3(b)) push three more A anchors
(Σp = 3.0) over 8 M☉, so the window loses 11, not 8, and the closure rise
overshoots the 20 % ceiling by about 1 point. The defect's effect was
understood in direction and in its R_V = 3.0/3.1 size. Its R_V = 3.5 size was
underestimated because the forecast was built only from the R_V = 3.1
diagnostic.

## 5. Before / after (baseline: PARSEC, R_V 3.1, α 2.3, coeval, all-explode)

| quantity | repair_v7 | repair_v8 |
|---|---:|---:|
| N_death, association | 8.43 | **8.36** |
| N_death A / B / C | 4.17 / 4.26 / 0.00 | **4.10** / 4.26 / 0.00 |
| headline range, 36 branches | 5.63–28.74 | **5.57–28.74** |
| all 54 branches | 1.93–28.74 (median 8.79) | 1.94–28.74 (median 8.72) |
| retained-36 median | 13.29 | 13.16 |
| P(≥1 event) / P(last < 100 kyr) | 0.9997 / 0.552 | 0.9996 / 0.547 |
| k_A (baseline) | 1,734 | 1,692 (−2.4 %) |
| counts-based age A / B / C (Myr) | 4.00 / 4.09 / 2.52 | 4.01 / 4.09 / 2.52 |
| association stellar mass (M☉) | 29,246 (29,122 at v6) | 29,014 |
| mass vs Wright+15 like-for-like | 1.47× | 1.45× |
| A census above 8 M☉, probabilistic / thresholded | 59.4 / 57 | **67.4 / 65** |
| A closure ratio, baseline | 0.865 | **1.006** |
| closure grid median at 2.3 (18 cells) | 1.067 | 1.073 |
| star-weighted association closure | 1.154 | 1.171 |
| closing slopes A / B / C / association | 2.34 / 2.25 / 2.06 / 2.20 | **2.29** / 2.25 / 2.05 / 2.19 |
| living ledger above 8 M☉ | 380.6 | **388.6** |
| WP5 cells passing | 40/54 | 40/54 (different four) |
| WP5 gate G3 | PASS | **FAIL** (see §6) |
| headline cells passing / combos / all-pass branches | 27/36 · 5/12 · 15/36 | **25/36** · 5/12 · 15/36 |
| S (scenario score) headline | 0.323–0.736 | 0.323–0.733 |
| WP9 verdict | INCONCLUSIVE | INCONCLUSIVE |
| α split (S > 0.5) | 18/18 vs 0/18 | 18/18 vs 0/18 |

Every progenitor is still above ~34 M☉. Full table:
[issue19_before_after.csv](../tables/issue19_before_after.csv).

## 6. Things that changed verdict or wording, and need a human

1. **WP5 gate G3 fails on repair_v8.** The baseline still passes in all three
   subgroups, the mass is within a factor two, and the grid is still 40/54. The
   no-A/C-regression clause fails: A-PARSEC α = 2.0 at R_V 3.0 and 3.1 go from
   pass to fail on the residual trend (p = 0.0048). Two α = 2.6 cells go the
   other way (A-PARSEC and A-MIST at R_V 3.5). Removing 8 stars from the top of
   A's window steepens the preferred slope, which is physical. Per the adoption
   rule this does not revert to repair_v7 and no threshold moved. **Open as an
   issue: is G3's regression clause, written for model changes, meant to apply
   to a data-defect fix?**
2. **α-plausibility wording.** On repair_v8 the calibration-window median χ²
   is α 2.0 = 10.4, 2.3 = 7.00, 2.6 = 10.06. The paper said α = 2.6 "has the
   worst median"; that is no longer true. The pre-registered criterion (A1)
   compares 2.3 against 2.6 only, and it still passes. The sentence now quotes
   all three medians. The headline set (α ∈ {2.0, 2.3}) is unchanged, as
   registered, but it now drops a slope that fits the window *better* than one
   it keeps. **That is a decision for you, not for this rerun.**
3. **A's closure.** A now closes at Salpeter on the baseline branch (1.006; grid
   median 1.02), with closing slope 2.29. "A and B want Salpeter, C wants
   shallower" still holds (2.29 / 2.25 / 2.05), and more cleanly than before.
4. **α-headline adoption D1-P4.** Its pre-registered threshold quotes 8.43,
   which fails literally on repair_v8. The chain reading (the restriction does
   not move the baseline) passes and governs on non-legacy chains. Both are
   recorded in `wp7_alpha_headline_adoption_outcome_repair_v8.json`.
5. **L6 convergence letter** flips FAIL → PASS at 2×10⁶ (worst drift 0.00998
   against 0.01). This is Monte-Carlo noise in near-zero cells, not an
   improvement; material-cell drift is 0.17 %, as before.

## 7. Deviations and errors during execution, recorded

- **WP11 iteration count.** The chain first ran WP11 at the script default,
  2×10⁵, instead of the published 5×10⁵, and the first repair_v8 WP12 hash
  record pinned that product. WP11 was re-run at 5×10⁵. The first record was
  **renamed, not deleted**, to `wp12_revision_prereg_repair_v8_void_wp11_at_200k.json`
  and is referenced from its replacement. WP12 was re-run against the
  replacement.
- I2 and I3 each failed once because of defects in the checks themselves (§3).
  Both are recorded in the integrity file.

## 8. Not done (left open)

- **F3(a)**, Mass → Mini. Deferred by D1 to a separate pre-registration.
- **fig04 anchors.** The brief's "draw the 61 unlabelled anchors, or adopt
  `draft_fig04_cmd_redesign.py`" is a presentation choice for the user.
- **Slides.** `slides/cygob2_supernova_history_talk.pptx` was open in
  PowerPoint and was not edited. Slide 13 (closing α 2.25, 6.7 %), slide 16
  (8.43, 0.9997, 0.552) and slide 17 (8.43, 5.63–28.7) need the repair_v8
  values.
- **WP13** threshold freezing (D2) is a WP13 task.
- **Pre-existing, not #19:** `wp6_external_crosschecks.py` still reads the
  repair_v6 normalization (never moved to v7). Its table is a declared input
  but feeds no macro. It was not rerun, to avoid mixing two changes.
- WP5 report tables (`wp5_imf_norm_repair_v6.csv`, gate record v6) remain the
  accepted repair_v6 gate record, as at repair_v7.
