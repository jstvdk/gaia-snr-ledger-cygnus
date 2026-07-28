# WP5 completion report (repair_v5)

**Verdict: BLOCKED.**

This report supersedes [wp5_completion_report.md](wp5_completion_report.md)
(frozen pre-repair run, 0/54 branches passing). Gate record:
[provenance/wp5_repair_v5_gate.json](provenance/wp5_repair_v5_gate.json).

## Cause of the blocking failure (issue #1c)

The WP5 injection truth model generated synthetic photometry at a single age (the upper-MS MAP).  At CygOB2-B's MAP age the PMS/Henyey isochrone fold sits on the 3.17-4.0 Msun bin boundary, so a truth-age error displaced the model caustic and produced a sharp, subgroup-specific, mass-localized residual that survived every instrumental re-run.  Evidence: provenance/wp5_bump_shape_diagnostic_execution.json (no single IMF slope explains it; the fold location), provenance/wp5_age_scan_execution.json (the residual is strictly monotone in truth age and vanishes toward older ages), provenance/wp5_stable_label_refit_execution.json (labels exonerated), provenance/wp5_lower_edge_scan_execution.json (completeness ramp excluded).

## Fix

Truth-side age marginalization by a joint age-k fit with the WP4 age posterior as prior, applied identically to all three subgroups.  No gate threshold was moved; no per-subgroup age was chosen; with a single truth-age node the estimator reduces to the unmodified wp5_fit_imf.fit_one bit-for-bit.

## Baseline gate (PARSEC, R_V = 3.1, α = 2.3)

| subgroup | χ² p | trend p | max abs residual | fitted truth age (Myr) | pass |
|---|---:|---:|---:|---:|:--:|
| CygOB2-A | 0.4476 | 0.468 | 1.34 | 3.981 | yes |
| CygOB2-B | 0.2929 | 0.111 | 2.27 | 4.150 | yes |
| CygOB2-C | 0.3504 | 0.957 | 1.56 | 2.512 | yes |

## Why this verdict

The fix itself works and is not in question: the baseline branch passes for all three subgroups, CygOB2-B's worst residual falls 3.62 -> 2.51 with k moving only -3.6%, and the branch grid improves 26 -> 29 of 54.  What blocks acceptance is a measurement property of the trend diagnostic, not the age marginalization.

**Blocking reason.** G3 sub-criteria not met: no_A_or_C_branch_regression.  Adopted reading of 'A and C must not regress' is the STRICT per-branch one (science decision, 2026-07-27).  The two flagged CygOB2-C cells (MIST, R_V=3.1, alpha=2.0 and 2.3) are single-truth-age-node cells whose estimator is provably identical to repair_v3, and their flips are measured Monte-Carlo noise on the six-bin rank trend statistic, not an effect of the age marginalization — see provenance/wp5_trend_stability_check_execution.json.  Under the strict reading this still blocks acceptance, and the trend diagnostic itself is the next work item (issue #11).

Branch grid: **38/54** subgroup-branch fits pass,
under the retention policy in `CUTS_AND_THRESHOLDS.md` §13.
Association mass 29185
[28116,
30275] M☉, within a factor two of the
16,500 M☉ literature scale: True.

## Downstream authorization

`downstream_wp6_authorized = false`.

**WP6 may not start until open issue #3 is handled**: bright-mass completeness
plateaus near 0.8, not 1.0, so WP6 step 2(a) must divide by the injection
response or it will manufacture a spurious ~20% massive-star deficit and bias
N_SN high.

## Carried limitations

- **Issue #1d** — CygOB2-B's A_V is set almost entirely by broadband photometry (no spectroscopic anchors near B), so the Teff/A_V degeneracy the anchors were meant to break is reintroduced.  An extinction-scale error moves B's stars across the isochrone fold coherently.
- **Issue #9** — The fitted truth-age posterior for CygOB2-B sits above its upper-MS MAP, in the direction of the documented upper-MS/PMS indicator disagreement.  The joint fit uses the same counts to weight the age nodes and to score the gate, and the chi-square dof is not reduced for the effectively fitted age.
- **Issue #4** — No absolute 95% completeness edge exists on this field; the corrected_no_absolute95_edge fallback is in force.

## Reproduction

```bash
bash scripts/run_repair_v4_chain.sh                      # nodes + joint fit
PYTHONPATH=scripts python3 scripts/wp5_report.py --wp5-version repair_v5
PYTHONPATH=scripts python3 scripts/wp5_repair_v4_finalize.py --wp5-version repair_v5
```

Full evidence chain: `reports/WP5_AGE_CONDITIONAL_SCAN_repair_v3.md` and the
provenance records listed in the gate JSON.
