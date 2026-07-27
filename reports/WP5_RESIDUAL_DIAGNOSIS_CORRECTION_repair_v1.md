# Correction to the WP3–WP5 repair report (repair_v1)

*Created 2026-07-27. Supersedes the "Remaining problem" section of
`reports/WP3_WP5_REPAIR_REPORT_repair_v1.md`. That report is left byte-identical
on purpose — its sha256 is registered in `provenance/wp3_repair_gate.json`, and
per the project's own convention a repair attempt is a new versioned artifact,
never an in-place edit. This document is the amendment; the frozen report is the
record of what was believed at the time.*

## Summary

The repair report's **verdict is correct** — WP5 fails its gate and WP6 is not
authorized. Its **diagnosis is wrong**. The report attributes the failure to
CygOB2-C as a subgroup and proposes investigating subgroup-label contamination,
intrinsic age spread, and spatially varying selection. The failure is not
specific to CygOB2-C and none of those three hypotheses is needed to explain it.

The cause is a specification error in the forward response: the parent mass
distribution is integrated only up to 8 M☉, while the mass estimator has a
22–25% scatter that carries stars from above 8 M☉ down into the top observed
bin. The top bin is therefore under-predicted by construction, in every
subgroup.

## Evidence that the failure is not CygOB2-C-specific

From `tables/wp5_baseline_residuals_repair_v1.csv`, the top mass bin
(6.35–8.0 M☉) holds the most positive Pearson residual in **all three**
subgroups:

| subgroup | top-bin residual | rank among its 6 bins |
|---|---:|---|
| CygOB2-A | +2.63 | most positive |
| CygOB2-B | +1.08 | most positive |
| CygOB2-C | +3.26 | most positive |

A subgroup-specific pathology does not produce the same sign, the same bin, and
the same rank in three independent subgroups. CygOB2-C is not anomalous; it is
the one that happened to cross the `max_abs_pearson_residual <= 3.0` threshold.
Note also that CygOB2-C **passes** the χ² criterion (p = 0.0397 against a 0.01
threshold) and the trend criterion (p = 0.623). One bin, in one corner of the
window, is the entire blocker.

## Root cause

`scripts/wp5_common.py:59-62` sets

```python
MASS_GRID = np.round(np.arange(0.50, 8.0001, 0.25), 2)
CALIBRATION_HI = 8.0
```

with a comment stating that the grid is deliberately extended *below* the 2 M☉
calibration edge so the response can absorb low-mass stars scattering **up**
into the window. `scripts/wp5_fit_imf.py:202-206` then clips the response to
`[MASS_GRID.min(), CALIBRATION_HI]`.

The downward extension is present and correct. **There is no upward
counterpart.** The parent population stops dead at the same 8.0 M☉ value that
defines the observed window, so no star above 8 M☉ can contribute to any
observed bin.

### Why the 8 M☉ ceiling was inherited by the wrong quantity

`CUTS_AND_THRESHOLDS.md` §7.1 classifies the upper edge as Class A:

> 8 M☉ is where stars start dying, so above it the surviving population no
> longer traces the birth population.

That reasoning is sound for the **observed counting window** — which stars you
tally. It was then silently applied to the **parent integration range** — which
stars the model allows to exist. Those are different objects, and the second
does not follow from the first:

- At the WP4 ages (upper-MS MAP 3.16–4.50 Myr), essentially nothing above 8 M☉
  has died yet; main-sequence lifetimes exceed the association age up to several
  tens of M☉. Stars at 9–15 M☉ are physically present in the field right now.
- WP4 itself counts them: `provenance/wp3_repair_gate.json` records
  `median_mass_n_gt8 = 302` and `posterior_expected_n_gt8 = 310.9`.
- Those ~300 stars are in the observed sample and can be *measured* below 8 M☉,
  but the model assigns them zero probability of existing.

The 8 M☉ edge is a convention about what the calibration window counts, not a
claim that the mass function terminates there.

### The estimator is broad enough for this to matter

Measured directly from `data/processed/wp5_injection_response_repair_v1.parquet`
(PARSEC, R_V = 3.1), the recovery kernel has a width of 22–25% in mass. For
stars injected at exactly 8.0 M☉:

| subgroup | P(recovered in 6.35–8.0) | P(recovered > 8) |
|---|---:|---:|
| CygOB2-A | 0.29 | 0.40 |
| CygOB2-B | 0.28 | 0.45 |
| CygOB2-C | 0.29 | 0.39 |

Roughly 40% of true-8 M☉ stars are measured as heavier than they are. By the
same kernel, stars at 9–12 M☉ are measured into the top bin at a comparable
rate. That flux is real, it is in the data, and the model has no term for it.

## Quantitative check

Re-fitting all 54 branches with the parent extended above 8 M☉, using the
empirical log-mass kernel taken from the existing response file:

| | before | after |
|---|---:|---:|
| CygOB2-C baseline max abs residual | 3.26 | **1.34** |
| CygOB2-C baseline χ² p | 0.040 | **0.75** |
| CygOB2-A baseline max abs residual | 2.63 | 1.13 |
| CygOB2-B baseline max abs residual | 2.75 | 2.16 |
| baseline all-subgroup gate | FAIL | **PASS** |
| full branch grid passing | 31 / 54 | **43 / 54** |
| median shift in k | — | **+9%** |

**This is a scoping estimate, not the repair.** It extrapolates the response
with an analytic log-normal kernel rather than simulating it. The actual fix is
to extend `MASS_GRID` and re-run the injections end to end; the delivered
numbers will differ in detail.

## Why this is not gate tuning

The repair report correctly refuses to keep adjusting the extinction prior or
mass-likelihood width, on the grounds that doing so would be tuning to the gate.
That caution does not apply here, for three reasons:

1. **It was identified from the residual pattern, not the p-value.** The
   signature — same bin, same sign, all three subgroups — is diagnostic on its
   own and would be visible with no gate defined at all.
2. **It is the symmetric counterpart of a correction the pipeline already
   accepts.** The downward extension to 0.5 M☉ exists for exactly this reason
   and is documented as such in the source.
3. **It changes the physical answer, not just the test statistic.** k moves by
   ~9%, which propagates into the WP5 normalization and therefore into N_SN at
   WP7. A tuning change would move the p-value and leave the science untouched.

## Disposition of the three hypotheses in the frozen report

| hypothesis | status |
|---|---|
| subgroup-label uncertainty / contamination | **Demoted, not dismissed.** Not needed to explain the residual. It is a genuine open question, but `membership_probability` measures cluster-vs-field, not A-vs-B-vs-C, so it has never been quantified. Cheap check: label-stability confusion matrix over the existing seed scan in `provenance/wp2_hdbscan_subgroup_scan.csv` and `provenance/wp2_gmm_seed_stability.csv`. Non-blocking. |
| intrinsic age spread vs single-age model | **Demoted.** No longer motivated by this residual. Remains relevant to WP4's documented indicator disagreement. |
| spatially varying selection | **Demoted.** The injections are already run per-subgroup footprint with the real extinction distribution. |

## Required actions

1. Extend `MASS_GRID` in `scripts/wp5_common.py` to ~20 M☉. Keep
   `CALIBRATION_HI = 8.0` — only the injected parent range changes, not the
   observed window.
2. Confirm the isochrone phase cut in `scripts/wp5_common.py:128-133`
   (`label <= 1` / `phase <= 1`) still returns main-sequence points at
   10–20 M☉ for the fitted ages.
3. Re-run `scripts/wp5_injections_repair.py`, then re-fit and re-gate all 54
   branches.
4. Regenerate `tables/wp5_imf_norm.csv` — the committed table is the
   **pre-repair** run (0/54 passing, residuals to 12.3) and is stale against
   `data/processed/wp5_imf_normalization_repair_v1.parquet` (31/54).
5. Add the upper-edge distinction to `CUTS_AND_THRESHOLDS.md` §7: the window
   edge and the parent integration range are separate numbers with separate
   justifications.
6. Decide and record a policy for branches that still fail after the fix — they
   cluster at R_V = 3.5 and α = 2.6, the grid corners. Plan §1.4 forbids
   silently dropping branches.

## Separate finding: the bright-mass completeness ceiling is a WP6 trap

Not part of the blocker, but it must not reach WP6 unflagged.

`tables/wp5_completeness_baseline_repair_v1.csv` records bright-plateau
completeness of 0.833 / 0.854 / 0.831 for A / B / C. Decomposing the selection
chain at ≥6.5 M☉ (PARSEC, R_V = 3.1):

| subgroup | query pass | quality pass | membership given CMD-ready |
|---|---:|---:|---:|
| CygOB2-A | 0.998 | 0.838 | 0.980 |
| CygOB2-B | 0.999 | 0.899 | 0.948 |
| CygOB2-C | 0.998 | 0.877 | 0.939 |

The ceiling is almost entirely the **WP2 quality filter** (RUWE / BP-RP excess),
not the magnitude limit and not the clustering. This is physically expected —
massive stars in Cygnus are binary-rich and crowded — and the response
correction handles it inside WP5.

It becomes a hazard at WP6. Plan WP6 step 2(a) assumes high-mass completeness
"should be ~1". It is ~0.83. A closure test that takes the observed massive-star
count at face value will manufacture a spurious ~17% deficit and attribute it to
runaways or hidden stars, **biasing N_SN high** — the exact failure mode in the
plan's own risk register. WP6 must divide by the same injection response.

Relatedly: no subgroup reaches an absolute 95% completeness edge anywhere, so
the procedure prescribed in `CUTS_AND_THRESHOLDS.md` §7.1 ("find the mass where
recovery ≥ 95%; that is your lower edge") is unachievable on this field. The
frozen report states the correct fallback — a response-corrected 2–8 M☉
likelihood, explicitly *not* a "95% complete" claim. That deviation should be
recorded as a formal supersession of §7.1 rather than left as a report remark.
