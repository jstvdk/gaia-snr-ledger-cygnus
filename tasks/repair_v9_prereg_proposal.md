# repair_v9 chain: pre-registration proposal (for owner sign-off)

Date: 2026-10-07.

Basis:
- brief §4 ([stage_repair_v9_wp4_age_redesign_and_wp13.md](stage_repair_v9_wp4_age_redesign_and_wp13.md));
- the Stage 1 outcome ([reports/wp4v9_age_redesign.md](../reports/wp4v9_age_redesign.md));
- owner decisions 6 and 7 (`provenance/decisions_2026_10_07.json`).

After sign-off this becomes `scripts/repair_v9_prereg.py` → `provenance/repair_v9_prereg.json`, written once and committed before any chain step.

## 1. Ages consumed

The new script `scripts/repair_v9_headline_ages.py` writes `data/processed/wp4_age_posteriors_repair_v9_headline.parquet`. The Stage 1 table is not edited.

| subgroup | source | PARSEC / MIST, R_V 3.1 (68 %) |
|---|---|---|
| A | Stage 1 `ums` row = Phase A′ test c | 3.16 (2.80–3.50) / 3.57 (3.31–3.70) |
| C | Stage 1 `ums` row = Phase A′ test c | 3.55 (3.32–3.97) / 4.01 (3.76–4.27) |
| B | **A + C combined:** the product of the two test-c likelihood curves, same family and R_V. Flagged "not measured for B". | 3.55 (3.16–3.70) / 3.57 (3.52–3.88) |

Each row is the same curve for every f_bin and dmu key. Test c has no variants for these.

## 2. Headline chain (brief §4.2, unchanged in substance)

**Code changes before the run:**
- add `repair_v9` to `chain.py`;
- give `wp4_mass_posteriors_repair.py` an `--age-version` argument, since it currently reads `wp4_age_posteriors_{WP_REPAIR_VERSION}`;
- fix the stale `mass_method` label.

**Steps:**
1. Anchor masses at the v9 ages, per R_V; present-day `Mass` as in v8.
2. Photometric mass posteriors.
3. WP5 truth-age node injections: the repair_v6/v7 rule, 9 unsnapped nodes per posterior with the isochrone interpolated. This is 162 nodes, about 2.5 h. Then the joint age–k fit.
4. WP6: mass-extension injections at the new ages, closure, census, slopes.
5. WP7: ledger at 2,000,000 iterations, with the age, BH and convergence scans and the binary bound. α = 2.3 is the headline.
6. WP8–WP12: a new `wp12_revision_prereg_repair_v9.json`; WP11 at 500,000 iterations.
7. Manuscript plumbing: `wp10_inputs`, `numbers.tex`, tables, figures, `wp10_validate`.

## 3. Age scan (decision 7) — new

**Grid:** 11 native ages per family, the same for every subgroup:
- PARSEC 2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01, 5.62, 6.31;
- MIST 2.00, 2.25, 2.52, 2.83, 3.18, 3.57, 4.01, 4.50, 5.05, 5.67, 6.37.

**At each scan age t, all three subgroups are set to that exact age:**
- anchor masses at t;
- photometric masses at t;
- one WP5 node injection at t (native age, no interpolation), then the WP5 fit, so **k is refitted at t**;
- the WP7 ledger per subgroup over all 54 branches, at 200,000 iterations.

This is 198 node injections (about 3 h) plus 11 mini-chains (about 20–30 min each, 4 in parallel). The scan does not run WP6 or WP8–12.

**Why three one-dimensional scans are enough:** subgroup populations are independent Poisson draws. So the deaths for any (t_A, t_B, t_C) are the sum of independent per-subgroup draws, built from the stored per-subgroup draws without rerunning.

**Outputs:**
- `tables/repair_v9_age_scan.csv`: per subgroup, branch and t: the mean, median, 16–84 % range and P(≥1) of N_death; k; P(last < 100 kyr);
- `data/processed/repair_v9_age_scan_draws.npz`;
- a figure of N_death against age for A, B and C.

The scan is reported as a sensitivity, never as the headline.

## 4. Integrity checks (all gate adoption)

- **I1 replay:** the repair_v5 ages through the v9 code reproduce the repair_v8 WP5 k and the WP7 baseline (8.36) within Monte Carlo noise.
- **I2:** the anchor masses are per branch, at the v9 ages.
- **I3:** every input hash matches this pre-registration.
- **I4:** no v9 script reads a v8 or unversioned age product (grep plus `wp10_inputs.FORBIDDEN`).
- **I5 (new):** for each subgroup, weighting the scan's N_death(t) by the headline posterior reproduces the headline per-subgroup mean within max(0.5, 15 %). This checks that the scan and the chain agree; the tolerance allows for the discrete grid.

## 5. Predictions (scored, never gating)

The figures in V1 are estimates made before the run, from repair_v8's k at fixed MAP ages (PARSEC, R_V 3.1, α 2.3).

- **V1, direction per subgroup:**
  - A falls (4.10 → about 0.9 at 3.16 Myr);
  - B falls (4.26 → about 2.2);
  - C rises (0.00 → about 2.6).

  **So the baseline N_death falls**, from 8.36 to about 5.7.

  Averaging over the posterior widths will raise these numbers, because N_death is convex in age. Scored as follows:
  - **V1a:** the PARSEC baseline < 8.36;
  - **V1b:** it lies in [4, 9];
  - **V1c:** A and B fall and C rises.

  This replaces the brief's V1, which held A and B fixed.
- **V2:** the brief's wording is kept: closure ratios move in the direction implied by each subgroup's turnoff change. The turnoff mass rises for A and B (younger) and falls for C (older). The expected sign of each closure ratio will be derived from `wp6_closure_test` before the pre-registration is written.
- **V3:** C's baseline N_death is above 0.5, so C is **not** below the first-death boundary. This is the novelty claim tested directly.
- **V4 (scan):** N_death(t) is non-decreasing in t for every subgroup and branch, and is 0 for t ≤ 2.5 Myr.
- **V5 (scan):** the subgroups' N_death(t) curves agree within 25 % at every t ≥ 3.2 Myr. They differ only through k, so equal ages give similar counts.

## 6. Adoption rule

repair_v9 is adopted (`ADOPTED = "repair_v9"`) if I1–I5 pass. Predictions do not gate adoption.

Then:
- `reports/repair_v9_completion_report.md`, with a before/after table for every headline number;
- close #21 and update #20;
- update the CLAUDE.md headline numbers and the withdrawn table.

## 7. After adoption

- Sensitivity items S1–S3 (brief §4.5).
- WP13, with its own pre-registration (brief §5).

## Runtime

Roughly 9–10 h unattended, run in `screen` + `caffeinate` with resumable steps:
- injections about 5.5 h;
- WP6 injections about 1–2 h (an estimate);
- scan fits about 1.5 h;
- headline chain about 1.5 h.
