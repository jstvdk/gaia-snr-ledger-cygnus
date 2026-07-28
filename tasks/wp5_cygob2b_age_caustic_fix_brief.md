# Fix brief — CygOB2-B residual is a mass-localized bump, not a slope: gated repair plan

> **Status update 2026-07-27 (steps 1–2 executed):** **G1 PASSED** — labels
> exonerated (bin-2 shift 0.009 < 0.5, no mass-dependent instability; issue #6
> closed). **G2 PASSED** — the bump residual is strictly monotone in the
> injection truth age (5.49 → 0.74 over 2.24 → 5.01 Myr) and two
> posterior-supported ages (3.162, 3.548 Myr) pass the full gate. **But the
> plain posterior-marginalized mixture FAILS** (bin-2 3.40, χ²p 4e-4; with
> 1 Myr SF spread 2.89, χ²p 2e-3) because B's WP4 posterior is bottom-heavy
> (σ_lo ≈ 0.02 vs σ_hi ≈ 0.43 Myr). The §4-authorized **joint age–k fit
> (WP4 posterior as prior) passes under both prior variants** (max|res|
> 2.48 / 2.06, χ²p 0.014 / 0.063) — step 3 must implement the joint fit, not
> plain marginalization. Full results:
> [WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](../reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md);
> gate records `provenance/wp5_stable_label_refit_execution.json`,
> `provenance/wp5_age_scan_execution.json`,
> `provenance/wp5_age_joint_fit_diagnostic_execution.json`.

*Diagnosis date: 2026-07-27. Working version `repair_v3`; WP5 BLOCKED (issue
#1c). This brief supersedes the hypothesis ordering of
`tasks/wp5_cygob2b_massfunction_brief.md` §3 — everything that brief excludes
stays excluded, and the run instructions in its §4 still apply. Evidence here:
`provenance/wp5_bump_shape_diagnostic_execution.json` (sha256 of every input),
produced by `scripts/wp5_bump_shape_diagnostic.py` (deterministic, read-only).*

---

## 1. New findings (2026-07-27 shape diagnostic)

**F1 — No single IMF slope explains B.** Rebuilding the forward rates from the
stored repair_v3 injection response and scanning α continuously over 1.8–4.0:

| family | R_V | best α | bin-2 residual at best α | min bin-2 residual over all α |
|---|---:|---:|---:|---:|
| PARSEC | 3.0 | 3.3 | 0.77 | 0.60 |
| PARSEC | 3.1 | 3.3 | **2.25** | **2.05** |
| PARSEC | 3.5 | 3.1 | 2.74 | 2.36 |
| MIST | 3.0 | 3.3 | 0.39 | 0.28 |
| MIST | 3.1 | 3.1 | 0.85 | 0.66 |
| MIST | 3.5 | 3.4 | 3.52 | 3.43 |

Two readings, both fatal for "B simply has a steeper IMF": (i) the α that best
fits B is ≈3.1–3.4 — far outside any credible cluster IMF and inconsistent
with A and C; (ii) even granting it, the 3.17–4.0 M☉ excess survives at
R_V ≥ 3.1 in PARSEC and at 3.5 in both families. B's anomaly decomposes into a
**global mass-scale displacement** (apparent steepening) **plus an R_V-coupled
localized bump**. *(Caveat: this reconstruction omits the Dirichlet response
draws and k posterior, so its absolute residuals are milder than the official
gate values — 2.87 vs 3.62 at the baseline. The α-shape conclusion is
unaffected: observed counts do not depend on α.)*

**F2 — The bump bin sits on an isochrone caustic.** G0(Mini) is non-monotonic
(the PMS/Henyey fold: one luminosity ↔ several masses) at the fitted MAP ages,
measured on the shared WP3/WP4 isochrones:

| subgroup | MAP age (baseline) | PARSEC fold (M☉) | MIST fold (M☉) |
|---|---:|---|---|
| CygOB2-A | 3.98 Myr | 2.40–2.80 | 2.39–2.77 |
| CygOB2-B | 2.82 Myr | 2.77–3.22 | 2.73–3.21 |
| CygOB2-C | 2.51 Myr | 2.90–3.36 | 2.86–3.36 |

The fold straddles the bin-1/bin-2 boundary (3.175 M☉) for B. Honesty note:
**C's fold overlaps bin 2 too, and C passes** — so fold location alone is
necessary, not sufficient. The B-specific ingredients are (a) WP4's documented
upper-MS vs PMS indicator disagreement, which exists **only for B**, and (b)
issue #1d: B's extinction is now set by broadband photometry alone (no anchors
near B), so an extinction-scale error moves B's stars across the fold
coherently. B is also the only subgroup whose bump is strongly R_V-monotone.

**F3 — The injection truth model cannot see an age error.**
`scripts/wp5_injections_repair.py` generates injected truth photometry from the
isochrone at a **single age — the upper-MS MAP** (`indicator=="ums"`, one
`age_map`; no posterior width, no SF-duration spread). The recovery side
correctly marginalizes the age posterior, but the truth side is a δ-function.
If B's real 2–4 M☉ population sits at a different age than B's UMS MAP, the
model caustic is displaced from the real one — producing exactly a sharp,
subgroup-specific, mass-localized residual that survives every instrumental
re-run at the same assumed age. This is the one class of error the end-to-end
injection design is structurally blind to.

**F4 — Contamination signatures are absent.** B's 86 bin-2 stars vs the rest
of B's window: membership P 0.852 vs 0.856, A_V 6.31 vs 6.37, parallax 0.6171
vs 0.6202 mas, RUWE 1.019 vs 1.032, median distance to B centroid 0.29° vs
0.34°. No P, parallax, astrometric-quality or spatial offset. Mislabeled A
stars would flood bins 0–1 (53% of A's window mass sits below 3.17 M☉), not
bin 2. Contamination is demoted, not excluded.

---

## 2. Considered and rejected: cutting the window at 3.3 M☉ (or any edge move)

Raising B's lower edge to exclude the bump region is rejected on five grounds:

1. **Already measured, runs the wrong way.** The lower-edge scan
   (`provenance/wp5_lower_edge_scan_execution.json`, issue #4b) shows B gets
   monotonically worse at edges 2.5/2.75 (3.62 → 5.11 → 5.23); B passes only
   at 5.0 M☉ and CygOB2-C then fails instead. "Do not retry blind."
2. **A cut on measured mass cannot excise a caustic.** The fold makes the mass
   map multi-valued: true 3.5 M☉ stars are measured below 3.3 and true 2.8 M☉
   stars above it. A cut re-samples the mis-modeled stars; it does not remove
   them.
3. **The anomaly spans the window.** The bump rides on a tilt (bins 3–5:
   −0.36/−1.30/−0.92). If the cause is a mass-scale error, the 4–8 M☉ counts
   and hence k are biased too; the cut hides the symptom the gate can see and
   keeps the disease WP6's N_SN inherits.
4. **It is the forbidden move.** The window's lower edge is Class C — measured
   by completeness (`CUTS_AND_THRESHOLDS.md` §7.1) — and §6.4 forbids tuning a
   threshold to pass a gate. 3.3 M☉ has no completeness justification; its only
   justification would be that it passes.
5. **Statistics.** Above 4 M☉ B keeps ~122 of 269 weighted sources; the k
   uncertainty grows ~50% for exactly the subgroup normalization WP5 exists to
   deliver.

---

## 3. Fix plan — every step carries a quality gate and a provenance record

No stored `repair_v3` artifact is overwritten at any step; scans use the
`wp5_lower_edge_scan.py` monkeypatch pattern, adoption runs get a new
`WP_REPAIR_VERSION`. No gate threshold moves.

### Step 0 — shape diagnostic (DONE, this brief)

- **Gate D0 (passed):** bin-2 residual ≥ 2 at *every* α on the baseline branch
  → "steeper single-slope IMF" excluded as a complete explanation.
- **Provenance:** `provenance/wp5_bump_shape_diagnostic_execution.json` ·
  `tables/wp5_bump_shape_diagnostic.csv` ·
  `scripts/wp5_bump_shape_diagnostic.py`.

### Step 1 — label-stability audit (cheap, closes issue #6 regardless of outcome)

Build the per-star A/B/C label confusion matrix over the 50 GMM seeds
(`provenance/wp2_gmm_seed_stability.csv`, `scripts/wp2_derive_subgroups.py`),
then refit the WP5 baseline restricted to stars with ≥90% seed-stable labels
(monkeypatch, no stored artifacts).

- **Gate G1:** if the stable-label refit moves B's bin-2 residual by **< 0.5**
  and no mass-dependent instability appears in B → labels exonerated, proceed
  to Step 2. Otherwise contamination becomes the primary path (Step 4a) and
  the confusion matrix is the evidence.
- **Provenance:** `provenance/wp2_label_stability_execution.json` +
  per-star CSV; refit record
  `provenance/wp5_stable_label_refit_execution.json`.

### Step 2 — age-conditional refit of B (the discriminating experiment)

Re-run B-only baseline injections (PARSEC, R_V 3.1) with the truth age
monkeypatched to ~6 nodes spanning B's WP4 age posterior plus B's
PMS-indicator age, refit each; then one run marginalizing the truth age over
the posterior nodes (posterior-weighted mixture), optionally with a 1 Myr
SF-duration spread (the plan §1.4 branch WP7 already mandates).

- **Gate G2:** B passes all three gate statistics (χ²p ≥ 0.01, trend p ≥ 0.05,
  max|res| ≤ 3.0) at an age **with non-negligible WP4 posterior support**, or
  under the posterior-marginalized mixture → mechanism confirmed, go to
  Step 3. If no age in B's central-95% posterior flattens both bump and tilt →
  age excluded, go to Step 4.
- **Provenance:** `provenance/wp5_age_scan_execution.json` recording every
  age node, per-node residual vectors, and sha256 of the response snapshots.

### Step 3 — adoption run `repair_v4` (only if G2 passes)

Implement truth-side age marginalization over the WP4 posterior (plus
SF-duration spread if G2 required it) **for all three subgroups uniformly** —
never a per-subgroup age choice — and re-run the full WP3→WP5 chain as
`repair_v4` per the §4 instructions of the previous brief.

- **Gate G3 (the unchanged WP5 gate):** baseline branch passes for all three
  subgroups; A and C must not regress; the 54-branch grid is reported with the
  issue-#5 retention policy for any remaining failures; association mass
  within factor 2 of literature.
- **Provenance:** full-chain `provenance/*_repair_v4_execution.json` set ·
  report `reports/WP5_AGE_MARGINALIZATION_repair_v4.md` stating cause and
  evidence · new gate record JSON.

### Step 4 — contingency if G2 fails

- **4a Contamination, now bump-constrained:** any contaminant must peak at
  3–4 M☉ *apparent* mass — reddened background red-clump giants do
  (M_G ≈ +0.45 at μ = 11.05). Extend the `wp4_clump.parquet` classification
  (31 rows today) to B's full window sample; flag, never delete.
  **Gate G4a:** flagged-contaminant refit flattens B *and* the classifier is
  justified by CMD/parallax evidence independent of the residual gate.
- **4b Per-subgroup extinction law:** continuous R_V reconstruction scan for B
  (F1 shows B nearly passes at R_V = 3.0 in both families).
  **Gate G4b:** adopting a B-specific R_V requires independent extinction
  evidence on B's sightlines (e.g. 2MASS colour-colour law fit), not gate
  score — issue #1d documents why anchors cannot arbitrate.
- **Provenance:** one execution JSON per sub-path, same conventions.

### Step 5 — closure bookkeeping (only after G3 or G4 passes)

Superseding WP5 completion report; `downstream_wp6_authorized` flipped in a
new gate record; `PROJECT_TRACE.md` §7/§8 + issue register (#1c closed with
stated cause, #2 regenerate `wp5_imf_norm.csv` after adding a version flag to
`wp5_report.py`, #5 policy recorded); `AUDIT.txt` refreshed via `audit.py`.

- **Gate G5:** `audit.py` clean; every number in the report traceable to a
  provenance JSON.

---

## 4. Anti-tuning rule for Step 2/3

Because B's age would be selected partly by its effect on B's residuals, the
adopted fix must be **marginalization over the existing WP4 posterior** (or a
joint age–k fit with that posterior as prior), applied identically to all
subgroups. Picking the single age that scores best on the gate is the same sin
as moving a threshold, one level up.

---

## 5. Artifacts of this diagnostic (2026-07-27)

| artifact | role |
|---|---|
| `scripts/wp5_bump_shape_diagnostic.py` | deterministic, read-only reproduction of F1/F2/F4 |
| `provenance/wp5_bump_shape_diagnostic_execution.json` | all numbers above + sha256 of the 8 inputs |
| `tables/wp5_bump_shape_diagnostic.csv` | full α-scan (6 branches × 23 α, per-bin residuals) |
| this brief | findings, rejected remedy, gated plan |

Input artifacts (hashes in the JSON): repair_v3 injection response,
mass-function bins, IMF normalization, WP4 mass and age posteriors, WP3
extinction catalogue, PARSEC and MIST isochrone tables.
