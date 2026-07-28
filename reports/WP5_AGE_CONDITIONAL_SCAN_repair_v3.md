# Steps 1–2 of the #1c gated plan: labels exonerated, age mechanism confirmed

*Executed 2026-07-27, on the plan in
[wp5_cygob2b_age_caustic_fix_brief.md](../tasks/wp5_cygob2b_age_caustic_fix_brief.md).
Working version `repair_v3`; no stored repair_v3 artifact was overwritten at
any step. WP5 remains BLOCKED pending the step 3 adoption run (`repair_v4`).*

---

## 1. Step 1 — label-stability audit: **gate G1 PASSED**, issue #6 closed

Per-star A/B/C label confusion over the 50 frozen WP2 GMM seeds
([wp2_label_stability.py](../scripts/wp2_label_stability.py) →
[wp2_label_stability_execution.json](../provenance/wp2_label_stability_execution.json)
+ per-star CSV):

- 1327 / 1331 clean members carry their consensus label in ≥ 90% of seeds
  (A: 2 unstable, B: 1, C: 1). Seed-averaged cross-assignment between any two
  subgroups is ≤ 0.16%.
- **No mass-dependent instability in B** — the smoking gun the superseded
  brief predicted for contamination is absent: window Spearman ρ = 0.011
  (p = 0.84); 341/342 of B's 2–8 M☉ stars are stable; the single unstable
  window star (bin 2) cannot produce a ~30-star excess.

Stable-label WP5 refit, full 54-branch grid, official `fit_one` unchanged
([wp5_stable_label_refit_execution.json](../provenance/wp5_stable_label_refit_execution.json)):

| G1 criterion | measured | limit | verdict |
|---|---:|---:|---|
| B baseline bin-2 residual shift | **0.009** | < 0.5 | pass |
| mass-dependent instability in B | none (ρ = 0.011, p = 0.84) | none | pass |

Grid stays 26/54; A/C baselines move ≤ 0.10. **Labels exonerated.**
Subgroup-label contamination is excluded as the cause of the B anomaly, and
issue #6 (label uncertainty never quantified) is closed with the per-star
matrix as the quantification.

## 2. Step 2 — age-conditional refit of B: **gate G2 PASSED**

B-only baseline (PARSEC, R_V = 3.1) injections re-run with the **truth** age
forced to eight native isochrone ages; recovery side untouched (repaired WP3
estimator + WP4 sampler marginalizing the production 9-node age posterior);
every fit through the unmodified official `fit_one`
([wp5_age_conditional_scan.py](../scripts/wp5_age_conditional_scan.py) →
[wp5_age_scan_execution.json](../provenance/wp5_age_scan_execution.json),
response snapshots hashed there). α = 2.3 baseline:

| truth age (Myr) | WP4 posterior weight (9-node) | bin-2 residual | χ² p | max abs res | gate |
|---:|---:|---:|---:|---:|:--:|
| 2.239 | 0 | **+5.49** | 1.5e-08 | 5.49 | no |
| 2.512 | 0 | +4.47 | 4.2e-06 | 4.47 | no |
| 2.818 (UMS MAP = production) | 6/9 | +3.77 | 0.0002 | 3.77 | no |
| 3.162 | 2/9 | **+2.85** | **0.0103** | 2.85 | **yes** |
| 3.548 | 1/9 | +2.29 | 0.038 | 2.29 | **yes** |
| 3.981 (A's age) | 0 | +1.81 | 0.206 | 1.81 | yes |
| 4.467 | 0 | +1.22 | 0.322 | 1.58 | yes |
| 5.012 | 0 | +0.74 | 0.607 | 1.27 | yes |

- The bump residual is **strictly monotone in truth age** and the entire tilt
  (bins 0–5) flattens with it — the caustic-displacement signature of finding
  F2/F3, not a bin-2 accident. The control node (2.818 = production truth,
  fresh RNG stream) reproduces the stored 3.62 at 3.77, so the Monte-Carlo
  noise floor (~0.15) is small against the node-to-node signal (~4.7).
- **G2 satisfied**: 3.162 Myr lies inside B's 68% posterior interval
  (2.814–3.246) and carries 2/9 of the recovery-side posterior weight;
  3.548 Myr carries 1/9. Both pass all three gate statistics.
- B's PMS indicator MAP could not serve as a node — it is unmeasurable (n = 2,
  grid-railed at 10 Myr); nodes 3.981–5.012 probe that direction instead, and
  the monotone improvement out to 5 Myr is consistent in sign with the
  UMS-vs-PMS disagreement WP4 documents for B only (issue #9).

**Plain truth-side marginalization is NOT sufficient.** Both
posterior-weighted mixtures fail, because B's WP4 posterior is extremely
bottom-heavy (σ_lo ≈ 0.02 Myr vs σ_hi ≈ 0.43 Myr, so 6/9 nodes collapse onto
the MAP):

| mixture (truth marginalized) | bin-2 residual | χ² p | gate |
|---|---:|---:|:--:|
| WP4 posterior, 9-node discretization | +3.40 | 0.0004 | no |
| same + 1 Myr SF-duration spread | +2.89 | 0.0018 | no |

## 3. Joint age–k fit forecast: **step 3 is viable**

The brief's anti-tuning rule (§4) authorizes exactly two adoptions:
posterior marginalization (fails above) or **a joint age–k fit with the WP4
posterior as prior, applied identically to all subgroups**. Forecast for B,
baseline branch
([wp5_age_joint_fit_diagnostic.py](../scripts/wp5_age_joint_fit_diagnostic.py)
→ [wp5_age_joint_fit_diagnostic_execution.json](../provenance/wp5_age_joint_fit_diagnostic_execution.json)):
per node, the Poisson marginal likelihood of B's official observed counts with
k integrated out analytically (Jeffreys prior); node weights ∝ WP4 prior × ML;
posterior-predictive composite fit with the unmodified `fit_one`. No gate
statistic enters the weights.

| prior | joint posterior over truth age | bin-2 res | χ² p | max abs res | gate |
|---|---|---:|---:|---:|:--:|
| WP4 9-node | 2.818: 0.02 · 3.162: 0.36 · **3.548: 0.62** | +2.48 | 0.014 | 2.48 | **yes** |
| + 1 Myr SF spread | 3.162: 0.12 · 3.548: 0.27 · **3.981: 0.54** · 4.467: 0.07 | +2.06 | 0.063 | 2.06 | **yes** |

The WP5 counts, under the WP4 prior, pull B's effective 2–4 M☉ truth age to
≈ 3.5–4.0 Myr — into the upper tail of the UMS posterior and toward the PMS
indicator's direction.

## 4. Decision and step 3 design implications

1. **Mechanism confirmed** (G2): B's residual is an age-conditional artifact
   of the single-age injection truth model sitting on the PMS/Henyey fold,
   not an IMF feature. Contamination (G1) and every instrumental channel
   (previous reports) are excluded.
2. **repair_v4 must implement the joint age–k fit**, not plain
   marginalization: truth-side age nodes from `age_posterior_nodes`, node
   weights updated by the WP5 Poisson likelihood under that prior, applied
   **identically to all three subgroups** (A and C already pass at their MAPs,
   so their joint fits are expected to stay put; that must be verified, not
   assumed). Carrying the 1 Myr SF-duration spread branch is favoured by both
   the mixture and joint-fit numbers and is already mandated as a WP7 branch.
3. Gate G3 is unchanged: full WP3→WP5 `repair_v4` chain, all three subgroups
   on the baseline, no regressions for A/C, 54-branch grid reported with the
   issue-#5 retention policy, association mass within factor 2.

## 4b. Steps 3a–3b executed: machinery implemented and validated

**3a — implementation.** [wp5_joint_age_fit.py](../scripts/wp5_joint_age_fit.py)
implements the joint fit; `inject_curve` gained an optional
`truth_age_override` (default `None` reproduces the frozen repair_v1–v3
behaviour exactly, and the override is consumed before any RNG draw, so a
branch's nodes differ *only* in truth age).

The node rule introduces **no new parameter**: the truth-side nodes are the
recovery side's own `age_posterior_nodes` discretization of the WP4 posterior
(nine equiprobable split-normal nodes), snapped to the native isochrone ages
the truth generator already snaps to, with prior weight = summed node count.
Across all 54 branches this needs **36 node responses**, because for narrow
posteriors several nodes collapse onto one native age.

**3b — validation** ([wp5_joint_fit_baseline_check_execution.json](../provenance/wp5_joint_fit_baseline_check_execution.json)),
zero new injections required:

- **Single-node equivalence is exact.** On the baseline branch A and C each
  span one native age, and `fit_joint` reproduces the official `fit_one`
  bit-for-bit — every gate field identical, max k-draw difference **0.0**.
  This is what makes "applied identically to all subgroups" verifiable: the
  new machinery *provably cannot move* a subgroup whose posterior spans a
  single node.
- **The baseline branch now passes for all three subgroups:**

| subgroup | nodes | truth-age posterior | bin-2 res | max abs res | χ² p | gate |
|---|---:|---|---:|---:|---:|:--:|
| CygOB2-A | 1 | 3.981 (fixed) | −0.66 | 1.36 | 0.356 | yes |
| **CygOB2-B** | 3 | 2.818: 0.02 · 3.162: 0.36 · **3.548: 0.62** | **+2.51** | **2.51** | **0.023** | **yes** |
| CygOB2-C | 1 | 2.512 (fixed) | 0.84 | 1.53 | 0.332 | yes |

  B's stored repair_v3 failure (bin-2 3.62, χ² p 0.0002) is resolved (2.51,
  χ² p 0.023) with k moving only −3.6% (1925 → 1855). A and C drift from the
  stored numbers by ≤ 0.012 in residual and ≤ 0.16% in k — pure Monte Carlo
  from RNG stream position (the production loop shares one generator across
  all 54 fits; the check uses a fresh one per fit), *not* a model difference,
  as the exact equivalence above proves.

## 4c. Step 3c executed: `repair_v4`

The chain ran 33 node injections (49 s each; 3 of the 36 nodes reused from the
G2 scan) and the 54-branch joint fit. **The baseline passes for all three
subgroups:**

| subgroup | max abs residual v3 → v4 | χ² p v3 → v4 | k v3 → v4 | fitted truth age | pass |
|---|---:|---:|---:|---:|:--:|
| CygOB2-A | 1.37 → 1.45 | 0.350 → 0.380 | 1724 → 1708 | 3.981 | yes |
| **CygOB2-B** | **3.62 → 2.51** | **0.0002 → 0.023** | 1925 → 1855 | 3.395 | **yes** |
| CygOB2-C | 1.52 → 1.52 | 0.341 → 0.492 | 1884 → 1873 | 2.512 | yes |

Grid **26 → 29/54**; association mass 30,155 [28,987, 31,312] M☉, within a
factor two. B's fitted truth-age posterior mean is 3.40 Myr on the baseline
and 3.16–3.96 Myr across branches — consistently above its 2.82 Myr upper-MS
MAP, in the direction of the PMS indicator (issue #9).

**Acceptance is nevertheless withheld, on a newly identified problem.** Gate
G3's "A and C must not regress" clause was resolved under the **strict
per-branch reading** (science decision), and two CygOB2-C cells (MIST,
R_V = 3.1, α = 2.0 and 2.3) flip pass → fail. Those flips are *not* caused by
the fix, and that is measured rather than argued
([wp5_trend_stability_check_execution.json](../provenance/wp5_trend_stability_check_execution.json)):

- Of the 7 flips between v3 and v4, **3 are the model change** (all fail →
  pass, driven by large max-residual improvements), and **4 are single-node
  cells whose estimator is provably identical to repair_v3** — two flipping
  each way, so the effect is unbiased.
- A **paired refit** of the identical estimator on the two independent
  Monte-Carlo realizations of the same model gives **opposite gate verdicts**
  in all four; `trend_p` is entirely insensitive to the fit RNG seed, so the
  driver is injection noise.
- **Mechanism:** with 6 bins the two-sided Spearman p-value is quantized
  (ρ = 0.771 → 0.072, 0.829 → 0.042, 0.943 → 0.005), so realization noise can
  push a *well-fitting* residual vector across p = 0.05 while χ² and max
  residual are unchanged or better.

This is issue **#11**, now the sole WP5 blocker. The fix must **replace the
diagnostic, never move the threshold** (`CUTS_AND_THRESHOLDS.md` §6.4): a
continuous weighted-least-squares residual slope, pre-declared and then
re-evaluated across all versions on equal footing. `accepted: false` and
`downstream_wp6_authorized: false` in
[wp5_repair_v4_gate.json](../provenance/wp5_repair_v4_gate.json).

## 5. Honest caveats

- The joint fit **uses the same observed counts to weight the age nodes and to
  score the gate** — inherent to any joint fit. The χ² dof is not reduced for
  the effectively fitted age (df stays 5); with one absorbed parameter the
  passing margins (χ² p 0.014–0.063) would tighten but the max-residual
  criterion (2.48/2.06 vs 3.0) is unaffected. State this in the repair_v4
  report.
- The ML node weights use point-estimate forward rates (no Dirichlet draws);
  the final gate statistics come from the full official fit. Consistent with
  the diagnostic caveat recorded in the shape-diagnostic JSON.
- `trend_p` is 0.111 for almost every fit — the 6-bin Spearman statistic is
  coarsely quantized and nearly insensitive here; the discriminating
  statistics are χ² p and max residual.
- The joint posterior concentrates at/beyond the WP4 90% upper bound
  (3.548 at 62%, or 3.981 at 54% with the SF spread). That is a tension with
  the UMS age likelihood, plausibly the same tension as issue #9 (UMS vs PMS
  indicator disagreement, B only) and issue #1d (unanchored extinction can
  displace B across the fold). repair_v4 should report the per-subgroup joint
  age posteriors next to the WP4 UMS posteriors so this tension stays visible.

## 6. Artifacts

| artifact | role |
|---|---|
| `scripts/wp2_label_stability.py` · `provenance/wp2_label_stability_execution.json` · `provenance/wp2_label_stability_per_star.csv` | step 1 stability matrix |
| `scripts/wp5_stable_label_refit.py` · `provenance/wp5_stable_label_refit_execution.json` | step 1 refit, gate G1 record |
| `scripts/wp5_age_conditional_scan.py` · `provenance/wp5_age_scan_execution.json` | step 2 scan, gate G2 record |
| `data/processed/wp5_age_scan_B_response_age*_repair_v3.parquet` (+ curves, + mixture composites) | per-node response snapshots (sha256 in the scan JSON) |
| `scripts/wp5_age_joint_fit_diagnostic.py` · `provenance/wp5_age_joint_fit_diagnostic_execution.json` | step 3 viability forecast |
| `scripts/wp5_joint_age_fit.py` | step 3a joint age–k fit machinery (library) |
| `scripts/wp5_joint_fit_baseline_check.py` · `provenance/wp5_joint_fit_baseline_check_execution.json` | step 3b validation: exact single-node equivalence + baseline joint fit |
| `scripts/wp5_injections_agenodes.py` · `scripts/wp5_fit_imf_joint.py` · `scripts/run_repair_v4_chain.sh` | step 3c adoption chain (written, not yet run) |
