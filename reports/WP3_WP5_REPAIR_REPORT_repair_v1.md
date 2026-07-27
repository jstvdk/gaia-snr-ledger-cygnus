# WP3–WP5 extinction/mass repair report (repair_v1)

## Verdict

**Repair remains BLOCKED at the unchanged WP5 residual gate.** The root WP3
extinction invariant and WP4 posterior-quality criteria pass, but the converged
baseline response still leaves the following subgroup failure:

| subgroup | chi-square p | trend p | max abs residual |
|---|---:|---:|---:|
| CygOB2-C | 0.0397001 | 0.622787 | 3.25698 |

WP6 is therefore **not authorized**. The frozen WP5 blocking verdict remains
preserved as the trigger, and these versioned files are a separate repair
attempt rather than replacements.

## What was repaired

- WP3 now uses a 0.03-mag error floor in every band, an eight-anchor spatial
  prior with measured widths of 0.452/0.453/0.475 mag, and a full gridded
  extinction posterior. The hidden template-branch width is calibrated on the
  149 spectroscopic anchors from the asymmetric central-68% discrepancy.
- WP4 ages were refitted from the repaired extinction catalogue. Masses are
  posterior samples over the full extinction distribution, all available six
  bands, the binary branch, and the subgroup age posterior. The six-band model
  width is the anchor-measured 0.38 mag; the mass measure is log-uniform, not an
  IMF prior. Spectroscopic-HRD overrides are unchanged.
- WP5 retained 223,200 catalogue injections and the original residual gate.
  Synthetic photometry passes through the actual repaired WP3 and WP4
  estimators. Baseline response resolution was increased from 16 to 64 mass
  draws and shown not to remove the remaining residual.

## Acceptance checks

| check | result | pass |
|---|---:|:---:|
| max absolute binned median ΔA_V | 0.066 mag (<0.30) | True |
| binned rank-test p-value | 0.623 (>=0.05) | True |
| largest repeated baseline mass | 0.80% (<2%) | True |
| stars at 2.5–3.2 M_sun | 276 | True |
| baseline WP5 all-subgroup residual gate | — | False |
| baseline association mass | 31293 M_sun | True |
| median-mass N(>8 M_sun) | 302 | reported |
| posterior E[N(>8 M_sun)] | 310.9 | reported |

No baseline subgroup reaches an absolute 95% completeness edge. The bright
plateaus remain below 95%, so the valid statement remains a response-corrected
2–8 M_sun likelihood, **not** “95% complete.” Relative-to-plateau edges are
diagnostic only and are recorded in `tables/wp5_completeness_baseline_repair_v1.csv`.

## Remaining problem

After the extinction/mass repair, CygOB2-A and CygOB2-B pass the baseline
residual gate. CygOB2-C retains a localized excess at the high end of the
calibration interval, with max |residual| above 3 despite acceptable global
χ² and trend tests. Because this survives the 64-draw response convergence
check, changing the extinction prior or mass-likelihood width further would be
gate tuning. The next diagnostic must test the C subgroup model itself:
subgroup-label uncertainty/contamination, intrinsic age spread versus the
single-age model, and spatially varying selection. Until one of those is
modeled and injected end-to-end, WP6 must remain blocked.
