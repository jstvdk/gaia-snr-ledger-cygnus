# Repair brief — CygOB2-B's mass function is steeper than the model and blocks the WP5 gate

> **Update 2026-07-27 (later the same day):** the shape diagnostic in
> [wp5_cygob2b_age_caustic_fix_brief.md](wp5_cygob2b_age_caustic_fix_brief.md)
> supersedes the hypothesis ordering in §3 below — B's anomaly is a
> mass-localized bump on an isochrone caustic plus a tilt, no single IMF slope
> explains it, and hypothesis §3.3 (age / single-age injection truth) is now
> the front-runner, with §3.1 demoted (contamination signatures measured and
> absent). Everything excluded in §2 stays excluded; §4 run instructions still
> apply. Follow the gated plan in the new brief.

*Diagnosis date: 2026-07-27. Pipeline status: WP5 BLOCKED, WP6 not authorized.
Current working version `repair_v3`; branch grid 26/54 passing. This document
states what is already excluded, so you do not repeat it. Hand this to the
executing agent together with `paper1_execution_plan.md`, `PROJECT_TRACE.md`
and `CUTS_AND_THRESHOLDS.md`.*

---

## 1. Problem statement

The WP5 Poisson IMF fit rejects CygOB2-B on the mass-function residual gate.
Baseline branch (PARSEC, R_V = 3.1, α = 2.3), `repair_v3`:

| subgroup | χ² p | trend p | max abs residual | pass |
|---|---:|---:|---:|:--:|
| CygOB2-A | 0.350 | 0.468 | 1.37 | yes |
| **CygOB2-B** | **0.0002** | 0.111 | **3.62** | **no** |
| CygOB2-C | 0.341 | 0.957 | 1.52 | yes |

Gate (`scripts/wp5_fit_imf.py`, `fit_one`): a branch passes only if
`chi_p >= 0.01` **and** `trend_p >= 0.05` **and**
`max_abs_pearson_residual <= 3.0`, for all three subgroups.

CygOB2-B's per-bin Pearson residuals are a monotone tilt, not one bad bin:

| bin (M☉) | 2.0–2.52 | 2.52–3.17 | 3.17–4.0 | 4.0–5.04 | 5.04–6.35 | 6.35–8.0 |
|---|---:|---:|---:|---:|---:|---:|
| repair_v2 | +1.04 | +1.11 | +3.02 | −0.41 | −0.51 | +0.25 |
| **repair_v3** | +1.45 | +2.41 | **+3.62** | −0.36 | −1.30 | −0.92 |

**Statement of the anomaly: CygOB2-B has more 2.5–4 M☉ stars, relative to its
4–8 M☉ stars, than a single-slope IMF folded through the measured completeness
response predicts.** A and C do not show this.

This is not a threshold artifact. Across the 54 fits × 6 bins there are 324
residuals; chance alone predicts ~0.9 above |3|. Fourteen are observed, and
nine of those fall in the single cell (CygOB2-B, bin 2). Concentration, not
magnitude, is the evidence.

---

## 2. What is already excluded — do not redo these

Three instrumental explanations were each traced, fixed or tested end to end,
and eliminated. Two of them were **real bugs that were correctly fixed and
still did not explain B.**

### 2.1 Parent-range truncation at 8 M☉ — real bug, fixed, not the cause

The forward response integrated the parent mass distribution only to 8 M☉ while
the WP4 mass estimator has 22–25% scatter, so the ~300 living members above
8 M☉ that scatter down into the top bin had no term in the model. Fixed in
`repair_v2` by extending `MASS_GRID` to 0.5–18 M☉ (ceiling measured: the >8 M☉
contribution converges to 100.0% at 18 for every α). Top-bin residuals fell in
all three subgroups (A −1.16, B −0.82, C −1.37) and CygOB2-C cleared entirely.
**B's bin 2 was untouched.**
→ `reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md`,
`reports/WP5_PARENT_RANGE_FIX_repair_v2.md`

### 2.2 Anchor extinction prior over-constraining B — real bug, fixed, made B worse

The repair_v1 prior applied a width calibrated at the anchor density (0.071°) to
B's members sitting 0.377° away, where the fitted variogram says the required
width is 2.32× larger. It collapsed B's differential extinction 9.5×
(1.860 → 0.196 mag) versus 2.0× for A. Fixed in `repair_v3` with a per-star
variogram-based width; B's spread recovered to 0.603 mag.
**Result: B's residual went 3.02 → 3.62 and the grid 33 → 26 of 54.** The
over-tight prior had been *masking* the tilt.
→ `reports/WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md`,
`provenance/wp3_anchor_prior_diagnostic_execution.json`

### 2.3 Completeness ramp / calibration lower edge — ruled out, effect runs backwards

B recovers only 6.2% of injected 2.0 M☉ stars, so a low-completeness artifact
was the natural candidate. Raising B's lower edge makes B **monotonically
worse**: 3.62 → 5.11 → 5.23 at edges 2.0 / 2.5 / 2.75. A ramp artifact would be
relieved by removing those bins. Only at a 5.0 M☉ edge does B pass, leaving 82
sources, and CygOB2-C then fails instead. Tested with both global edges
(repair_v2) and per-subgroup completeness-driven edges (repair_v3).
→ `scripts/wp5_lower_edge_scan.py`,
`provenance/wp5_lower_edge_scan_execution.json`

### 2.4 Also settled, so you don't reopen them

- **No absolute 95% completeness edge exists anywhere on this field** (bright
  plateau 0.79–0.85, limited by the WP2 quality filter, not the magnitude
  limit). `CUTS_AND_THRESHOLDS.md` §7.1 is formally superseded; the
  `corrected_no_absolute95_edge` fallback is in force.
- **The window edge and the parent integration range are separate numbers.**
  See `CUTS_AND_THRESHOLDS.md` §7.2. Do not re-tie them.
- **The residual gate must not be relaxed.** 3.62 is a symptom; the R_V ordering
  and the nine-branch concentration are the evidence and survive any reasonable
  threshold.

---

## 3. Candidate hypotheses, in the order to test them

### 3.1 Subgroup-label contamination — cheapest, no refit needed

`membership_probability` measures **cluster-vs-field only**. It has never
measured A-vs-B-vs-C. A star can be a certain Cyg OB2 member (P = 0.98) and
still be ambiguous between subgroups. If B is absorbing lower-mass stars that
belong to A or the field, that produces exactly this tilt.

WP2 labels came from Gaussian mixtures over 50 deterministic seeds, accepted on
seed stability. The per-seed information already exists:

- `provenance/wp2_gmm_seed_stability.csv`
- `provenance/wp2_hdbscan_subgroup_scan.csv`
- `scripts/wp2_derive_subgroups.py` (consensus co-assignment matrix)

**First step:** build a per-star label-stability confusion matrix over the 50
seeds. Then re-run the WP5 baseline weighting each star by its label stability,
or restricting to stars whose label is stable in ≥90% of seeds. If B's tilt
flattens, this is the cause.

**Watch for:** a mass-dependent label stability. If B's *low-mass* stars are the
unstable ones, that is the smoking gun.

### 3.2 A second population overlapping B

Check B's CMD, parallax and proper-motion distributions for a foreground or
background group at similar kinematics. Cygnus sits along a spiral arm — plan
§WP2 caveats flag contamination by unrelated young stars at similar distance and
PM as the dominant risk. B is the subgroup with the highest median A_V (6.39)
and the tightest kinematics.

### 3.3 B's age

`repair_v3` puts B's upper-MS MAP at 2.82 Myr against C's 2.51 and A's value in
`wp4_age_posteriors_repair_v3.parquet`. An age error tilts the inferred mass
function in exactly the observed direction, because masses are read off the
isochrone at the fitted age. **WP4 already carries a documented indicator
disagreement for subgroup B** (upper-MS vs PMS, N = 19 tail) —
see `wp4_completion_report.md`.

**First step:** refit B's WP5 baseline across B's full WP4 age posterior rather
than only the upper-MS MAP, and check whether any age in the posterior flattens
the residuals. If one does, the age is the free parameter, not the IMF.

---

## 4. How to run the pipeline

Versioning is enforced by two environment variables, both defaulting to
repair_v1 behaviour so older versions stay exactly reproducible:

```bash
export WP_REPAIR_VERSION=repair_v4        # output suffix for the whole chain
export WP3_ANCHOR_PRIOR_MODE=variogram    # "global" = repair_v1 prior
PYTHONPATH=scripts bash scripts/run_repair_v3_chain.sh
```

The chain is WP3 extinction → WP4 ages → WP4 masses → WP5 injections → WP5 fit,
about 25 minutes end to end. To re-fit WP5 only, without new injections:

```bash
PYTHONPATH=scripts python3 scripts/wp5_fit_imf.py \
    --repair-version repair_v3 --wp5-version repair_v4
```

`--repair-version` selects the upstream WP3/WP4 products; `--wp5-version`
controls the WP5 inputs and every WP5 output. `scripts/wp5_lower_edge_scan.py`
is the template for a cheap what-if scan that monkeypatches one function and
re-fits all 54 branches without touching any stored artifact.

**Binding conventions.** Never overwrite a previous version — every run writes
`*_repair_vN` artifacts plus a provenance JSON with sha256 of every input and
output. Model branches (2 isochrone families × 3 R_V × 3 IMF slopes) are carried
and reported, never averaged (plan §1.4). Do not tune a threshold to pass a
gate; replace the diagnostic instead (`CUTS_AND_THRESHOLDS.md` §6.4).

---

## 5. Definition of done

1. All three subgroups pass the residual gate on the baseline branch, **and**
   the cause is a stated, evidenced model change — not a threshold move.
2. The 54-branch grid result is reported, with an explicit retention policy for
   any branch that still fails (open issue #5; failures currently cluster at
   R_V = 3.5 and α = 2.6).
3. A completion report supersedes `wp5_completion_report.md`, and
   `downstream_wp6_authorized` flips in a new gate record.
4. Before WP6 starts, open issue #3 is handled: bright-mass completeness is
   ~0.80, **not ~1.0**. WP6's closure test must divide by the injection
   response or it will manufacture a spurious ~20% massive-star deficit and
   bias N_SN high.
5. `PROJECT_TRACE.md` §8 and the issue register updated; `AUDIT.txt` refreshed
   via `audit.py`.

---

## 6. Reading order

| # | document | why |
|---|---|---|
| 1 | `PROJECT_TRACE.md` | index of every artifact, status board, issue register (#1c is this task) |
| 2 | this brief | what is excluded and what to try |
| 3 | `paper1_execution_plan.md` §WP5 + §1 | objective, gate, binding conventions |
| 4 | `CUTS_AND_THRESHOLDS.md` §3.3, §6.4, §7 | how every number is chosen; the two amendments |
| 5 | `reports/WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md` | current state, full evidence |
| 6 | `reports/WP5_PARENT_RANGE_FIX_repair_v2.md` | the truncation fix and its verified effect |
| 7 | `method_explained.md` §WP5 | physics of IMF normalization |

Optional background: `reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md`
and the frozen `reports/WP3_WP5_REPAIR_REPORT_repair_v1.md` (its verdict stands,
its diagnosis is superseded — read only to avoid repeating its dead ends).
