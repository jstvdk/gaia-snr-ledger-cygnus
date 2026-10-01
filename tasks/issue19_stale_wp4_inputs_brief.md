# Issue #19 — stale pre-repair WP4 inputs: the age envelope and the frozen anchor masses

*Opened 2026-10-01. Register entry: [PROJECT_TRACE.md §9, issue 19](../PROJECT_TRACE.md).
Scope diagnostic (read-only, reproducible):
[issue19_scope_diagnostic.py](../scripts/issue19_scope_diagnostic.py) →
[issue19_envelope_by_version.csv](../tables/issue19_envelope_by_version.csv) ·
[issue19_anchor_mass_rederivation.csv](../tables/issue19_anchor_mass_rederivation.csv).*

```sh
PYTHONPATH=scripts python scripts/issue19_scope_diagnostic.py   # conda env cygob2-gaia
```

---

## 1. What is wrong, in two sentences

Two products from the **unversioned, pre-repair WP4 run (2026-07-23)** survived every
later repair:
1. the quoted **2.25–5.67 Myr "two-indicator" age envelope**, which is narrative only;
2. the **150 spectroscopic-anchor masses**, which are numerical and enter the WP5 2–8 M☉ window and the WP6 census above 8 M☉.

Neither was regenerated when WP3/WP4 were repaired (v1 → v3 → v5), so the
`repair_v7` chain mixes repaired photometric masses with anchor masses read off
isochrones at the wrong ages.

## 2. Findings, with evidence

### F1 — the age envelope is from the pre-repair run (narrative only)

| `wp4_age_posteriors` version | retained rows | retained PMS rows | MAP span |
|---|---:|---:|---|
| unversioned (pre-repair, 07-23) | 104 | 38 | **2.25–5.67 Myr** |
| repair_v1 | 66 | 0 | 1.78–4.50 |
| repair_v3 | 66 | 0 | 2.00–4.50 |
| **repair_v5 (the chain's)** | 66 | **0** | **2.00–4.01 Myr** |

- In the chain, **no pre-main-sequence (PMS) age indicator survives**: n = 15, 4 and 3 stars, all grid-railed. Every age is upper-MS only.
- Any statement built on the PMS indicator is therefore unsupported on repair_v5:
  - "both retained indicators";
  - "A's upper MS is ~1.5 Myr older than its PMS → extended star formation";
  - "B's PMS wants B older" (issue #9), cited as supporting evidence in [wp4_wp5_age_reconciliation.md §5](../reports/wp4_wp5_age_reconciliation.md).
- The ages the ledger consumes (WP5 counts-based, with repair_v5 as prior) are **not** affected. F1 changes text and figures only.

### F2 — the frozen anchor masses use pre-repair ages (numerical)

- [wp4_anchors_hrd.py:77](../scripts/wp4_anchors_hrd.py) reads the unversioned `wp4_age_posteriors.parquet`. It matches each anchor's spectroscopic T_eff and M_G0 to the isochrone at the **pre-repair** upper-MS ages.
- [wp4_mass_posteriors_repair.py:55](../scripts/wp4_mass_posteriors_repair.py) then copies those masses into every repaired mass catalogue as `spectroscopic_hrd_frozen`, **identical on every R_V branch**.

The ages used:

| | A PARSEC | B PARSEC | C PARSEC | A MIST | B MIST | C MIST |
|---|---:|---:|---:|---:|---:|---:|
| frozen (pre-repair) | 4.47 | 3.98 | 3.98 | 3.57 | 3.57 | 3.57 |
| repair_v5 | 3.98 | 3.55 | 2.51 | 4.01 | 4.01 | 3.18 |

The anchor magnitudes themselves are current: at most 0.08 mag from repair_v5. Re-deriving the masses with the same procedure at the repair_v5 ages:
- **PARSEC, subgroup A: 8 anchors (Σ membership probability = 8.0) move from the 2–8 M☉ window to above 8 M☉.**
- No crossings anywhere else, and none of the 2 M☉ boundary.
- Masses of C's most massive anchors change by up to ×2.7. This does not move any count, but it changes the stars' positions relative to m_TO.

### F3 — two smaller defects in the same procedure (decide whether to fix together)

- **(a) Present-day mass instead of initial mass.** `iso_hrd_points` returns the isochrone's present-day `Mass`, not `Mini`, but the IMF is in initial mass. The difference is negligible below ~20 M☉, but up to ~10–25 % for 40–60 M☉ stars at 4 Myr. It affects no 2 or 8 M☉ crossing, only positions relative to m_TO.
- **(b) Anchor masses are not R_V-dependent.** Photometric masses are; anchor M_G0 moves ~0.29 mag between R_V = 3.1 and 3.5.

### Not affected (verified or by construction)

- Membership and subgroups.
- Extinction: anchor M_G0 differs by at most 0.08 mag.
- The repair_v5 upper-MS ages and the 1,242 photometric masses.
- WP5 injections, the runaway traceback, and the ledger engine itself.

## 3. Expected size of the effect (estimates; the rerun decides)

| quantity | repair_v7 | expected after fix | basis |
|---|---|---|---|
| k_A, PARSEC branches | — | −1 to −3 % | 8 of ~401 window stars leave |
| A's N_death (baseline) | 4.17 | ≈ 4.05–4.15 | ∝ k_A |
| association N_death (baseline) | 8.43 | ≈ 8.3–8.4 | |
| A's census above 8 M☉ (PARSEC) | 57 | ≈ 65 | +8 anchors |
| A's closure ratio (PARSEC cells) | — | +10 to +20 % | census up, prediction down |
| A's closing slope (median of 6 cells) | 2.34 | lower, ≈ 2.2–2.3 | PARSEC cells move by ≈ −0.1; MIST cells do not |
| MIST branches, B, C | — | unchanged within ~1 % | no crossings |

**The headline numbers survive. The closure narrative may not:** "A and B want Salpeter, C wants shallower" and A = 2.34 must be recomputed before they are quoted again.

## 4. Decisions needed from the user before step 2

- **D1 — fix F3 now or separately?** Recommended: do F3(b) (per-R_V masses) now, because it is the same code path. Do F3(a) (Mini) as a **separately pre-registered, separately scored** change, so the two effects are attributable.
- **D2 — order relative to WP13.** Recommended: freeze WP13's materiality thresholds now (they don't depend on these numbers), but read M0 and M1 only on the repaired chain.
- **D3 — version name.** Recommended: `repair_v8` for everything regenerated. Nothing at v7 or earlier is overwritten (rule 3).

## 5. Work plan (each step leaves a versioned artifact plus provenance JSON)

> **Status 2026-10-01: DONE except step 9's fig04-anchor choice, the slides, and F3(a).** D1–D3 were taken at the recommended options. repair_v8 was adopted (I1–I4 PASS); P1 and P4 FAILED at R_V = 3.5 only; P2, P3, P5, P6 PASS. See [reports/issue19_completion_report.md](../reports/issue19_completion_report.md).

- [x] **1. Pre-register** (`scripts/issue19_prereg.py` → `provenance/issue19_repair_v8_prereg.json`), with SHA-256 of every input. Write the predictions before any rerun:
  - **P1:** PARSEC A loses 8.0 ± 0.5 window stars on every PARSEC branch.
  - **P2:** MIST branches and subgroups B and C: k and N_death within 1 %.
  - **P3:** baseline association N_death ∈ [8.20, 8.43].
  - **P4:** A's PARSEC closure ratios rise by 8–20 %.
  - **P5:** the association closing slope stays within [2.10, 2.30].
  - **P6:** the repair_v5 retained envelope is 2.00–4.01 Myr, and no PMS row is retained.

  Failed predictions stay recorded as failed.
- [x] **2. WP4 anchors.** Give `wp4_anchors_hrd.py` an `--age-version` argument (default must refuse the unversioned file). Use each branch's own age and R_V-specific M_G0, so there are 6 mass columns. Output `data/processed/wp4_anchor_hrd_repair_v8.parquet`. Apply D1.
- [x] **3. WP4 masses.** `wp4_mass_posteriors_repair.py` reads the versioned anchor file → `wp4_mass_posteriors_repair_v8.parquet`. Also fix the stale `mass_method` label `photometric_posterior_repair_v1` (cosmetic).
- [x] **4. WP5.** Rerun the joint age–k fit on repair_v8 masses (`wp5_fit_imf_joint.py --upstream-version …`) → normalization, posterior draws and association mass at `repair_v8`. Confirm the injection response is reusable (it does not use anchor masses); if not, rerun it too. Re-score the WP5 gate (40/54 fits and 15/36 headline branches may change).
- [x] **5. WP6.** Closure test with the 4.0 M☉ floor, massive census, orphan anchors, `wp12_closure_slopes.py`.
- [x] **6. WP7.** `wp7_ledger.py` (2,000,000 iterations), age scan, BH scan, alpha-headline branch sets, binary bound, convergence scan.
- [x] **7. WP9 / WP12.**
  - Verdict, gate map, combination gate, branch-gate table, mixed-slope ledger, scenario score.
  - `wp12_common.frozen()` checks inputs against `wp12_revision_prereg.json`. Write a **new** hash record for repair_v8; do not edit the old one.
- [x] **8. Manuscript plumbing.**
  - `wp10_inputs.py`: bump versions, and add the unversioned WP4 posteriors, the unversioned `wp4_anchor_hrd.parquet` and the two pre-repair WP4 markdown tables to `FORBIDDEN`.
  - `wp10_numbers.py`: add an age-envelope macro computed from the posteriors.
  - `wp10_validate.py`: remove `"2.25"` and `"5.67"` from the bare-number allowlist (line 55).
- [~] **9. Figures.**
  - In `wp12_figures.py`, compute the envelope instead of hard-coding it (lines 433–434, 461, 813), and relabel it "retained upper-MS envelope".
  - Make fig04 draw the 61 unlabelled anchors, or adopt the redesign draft ([draft_fig04_cmd_redesign.py](../scripts/draft_fig04_cmd_redesign.py)).
- [x] **10. Text.**
  - `wp4_report.py` is hardwired to the unversioned posteriors (the same defect class as issue #2): add `--version` and regenerate the tables in `wp4_ages.md`.
  - Correct the WP4 headline in `PROJECT_TRACE.md` (§ WP4, "2.25–5.67", "both retained indicators").
  - Correct `reports/wp3_obligations_discharge.md` (O1) and `reports/wp4_wp5_age_reconciliation.md` §5 point 2.
  - In `manuscript/main.tex`, fix lines ~594, ~610 and ~798 and any PMS or "extended star formation" claim.
  - Update `tasks/manuscript_reframing_revision_brief.md:836`, the group-meeting brief and the slides.
- [x] **11. Score and close.** Score P1–P6 into `provenance/issue19_repair_v8_outcome.json`, write `reports/issue19_completion_report.md`, update the headline numbers in `CLAUDE.md` if they moved, and close #19 in §9 with before/after values.

## 6. Acceptance criteria

1. No script in the chain reads an unversioned WP4 product. Check with `grep -n 'wp4_age_posteriors.parquet\|wp4_anchor_hrd.parquet' scripts/`: only historical scripts should remain, and they must be refused by `wp10_inputs`.
2. Anchor masses differ by R_V branch, and their ages equal the repair_v5 posteriors.
3. `wp10_validate.py` passes with no 2.25/5.67 allowlist; the envelope is a macro.
4. Every changed headline number has a before/after row in the completion report; nothing at repair_v7 is overwritten.
5. P1–P6 are scored, and failures are stated as failures.

## 7. Rules that apply (from CLAUDE.md)

- No bare numbers in `main.tex`.
- Inputs only through `wp10_inputs.py`.
- Nothing is overwritten and nothing is retuned.
- A cross-check that disagrees becomes an issue, never a reason to move a number.
- Predictions are written before the rerun.
