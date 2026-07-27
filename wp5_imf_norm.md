# WP5 — IMF normalization and completeness

**Status: BLOCKED AT THE VALIDATION GATE. Do not propagate these normalizations
to WP6/WP7 as accepted inputs.**

This work package was executed against the frozen WP2--WP4 products.  It
produced complete diagnostic artifacts and a forward-response IMF fit, but it
did not satisfy either the absolute-completeness requirement or the
mass-function residual requirement.

## 1. Frozen inputs and branch grid

The fit uses all three kinematic subgroups; PARSEC and MIST masses; R_V =
3.0/3.1/3.5; and high-mass IMF slopes alpha = 2.0/2.3/2.6.  Each observed star
contributes its WP2 membership probability.  The six logarithmic observed-mass
bins cover 2--8 Msun.  Unresolved binaries are injected with f_bin=0.4 and
q~U[0.1,1], matching the central WP4 binary branch.

The injection response extends in true primary mass from 0.5 to 8 Msun so that
lower-mass binaries and mass-estimation scatter into the observed 2--8 Msun
window are modeled rather than ignored.

## 2. Completeness experiment

Synthetic stars were generated with the branch-matched WP4 age and the same
frozen PARSEC/MIST isochrones.  Their extinction was drawn spatially from the
matching WP3 subgroup/R_V distribution.  Real Gaia DR3 observational states
were cloned from nearby sources of similar G and BP-RP in the one-degree field.
Every injection then passed the frozen G/parallax query, exact WP2 quality
filter, reconstructed cluster-versus-field posterior, P>0.5 handoff, and G/BP/RP
mass-readiness requirement.

The 128-point Sobol-normal membership integration agrees with the frozen WP2
10,000-draw decision for 98.95% of a stratified 760-source validation sample
(7 false negatives, 1 false positive; median absolute probability difference
0.0068).

Baseline completeness:

| subgroup | recovery_at_2_Msun | bright_6_8_median | max_monotone_recovery |
|---|---|---|---|
| CygOB2-A | 0.515 | 0.8 | 0.803 |
| CygOB2-B | 0.278 | 0.782 | 0.828 |
| CygOB2-C | 0.117 | 0.762 | 0.8 |

No branch/subgroup curve reaches 95% absolute end-to-end recovery anywhere in
2--8 Msun.  Gaia query recovery becomes essentially complete above about
2.25--3.25 Msun depending on subgroup, but the WP2 quality and membership stages
leave a bright plateau of only about 75--82%.  Therefore the planned absolute
95% lower edge does not exist.  The diagnostic 95%-of-plateau edges are retained
in the normalization table, but they are not relabelled as absolute
completeness.

Because the selection response is explicitly measured, a diagnostic corrected
fit over the nominal 2--8 Msun observed window remains mathematically possible.
It is labelled `corrected_no_absolute95_edge`, not "95% complete".

## 3. Forward Poisson IMF fit

The forward model predicts each *recovered* mass bin:

lambda_i = k integral R_i(m_true) m_true^(-alpha) dm_true,

where R_i contains query/quality/membership recovery, unresolved-binary
brightening, WP3 extinction error, photometric error, and migration through the
exact WP4 nearest-isochrone mass estimator.  A Jeffreys prior on k and a
Dirichlet posterior for the injection response yield 10,000 posterior draws per
subgroup and branch.

Baseline k posteriors:

| subgroup | k_median | k_lo68 | k_hi68 | raw_calibration_sources | membership_weighted_calibration_sources |
|---|---|---|---|---|---|
| CygOB2-A | 1.53e+03 | 1.44e+03 | 1.63e+03 | 349 | 306 |
| CygOB2-B | 1.46e+03 | 1.37e+03 | 1.56e+03 | 307 | 254 |
| CygOB2-C | 1.53e+03 | 1.42e+03 | 1.64e+03 | 265 | 220 |

The first observed-mass bin is the decisive failure:

| subgroup | membership_weighted_count | expected_count_at_k_median | pearson_residual |
|---|---|---|---|
| CygOB2-A | 169 | 103 | 6.44 |
| CygOB2-B | 141 | 73.5 | 7.86 |
| CygOB2-C | 101 | 50 | 7.17 |

All 54 subgroup x family x R_V x alpha fits fail the predeclared residual gate.
The best branch has p_chi2 = 3.44e-06;
the baseline subgroup p-values are
CygOB2-A: 1.3e-13, CygOB2-B: 1.8e-18, CygOB2-C: 1.4e-13.
The response-aware model still underpredicts the 2.0--2.52 Msun bin by
6.4--7.9 sigma while overpredicting several adjacent bins.  This is not a
smooth completeness slope that can be repaired by changing k or alpha.

The most direct evidence points to an incompatibility between the frozen WP4
point-mass catalogue and a power-law birth IMF in this lower-CMD regime.  WP4
assigns photometric masses by a nearest point on a single-star isochrone; the
2--5 Msun catalogue therefore needs per-star mass posteriors or a direct
CMD-space IMF likelihood before the IMF gate can be trusted.  This is a
diagnosis, not permission to rewrite the frozen WP4 products.

## 4. Association-mass diagnostic

For the baseline PARSEC, R_V=3.1, alpha=2.3 branch, the summed normalization
implies a primary-system mass of
20721
[19941,
21530] Msun and a
multiplicity-adjusted stellar mass of
25049
[24107,
26027] Msun.
That median is within a factor two of the Wright+15 16,500 Msun scale.  Across
all 18 association branches the medians span
24518--
33531 Msun, with
17/18 branches inside the
factor-two band.

This sanity check passes at baseline but cannot override the failed shape
residuals: a biased mass function can integrate to a plausible total mass.

## 5. Gate assessment

| Criterion | Result |
|---|---|
| Lower calibration edge >=95% complete | **FAIL** — no absolute 95% edge exists in any branch/subgroup |
| >=50 sources per subgroup in the corrected 2--8 diagnostic window | **PASS** — 248--368 raw per branch/subgroup |
| Residuals consistent with Poisson scatter; no hidden mass trend | **FAIL** — 0/54 fits pass; baseline first-bin excess 6.4--7.9 sigma |
| Association mass within factor about two of literature | **PASS baseline; one extreme branch fails** |

**Blocking conclusion:** WP5 is not accepted.  WP6/WP7 must not consume
`wp5_imf_normalization.parquet` as a validated normalization.

## 6. Required remediation

1. Reopen only the WP4 lower/intermediate-mass inference, replacing nearest-point
   masses with per-star mass posteriors or a direct CMD-space population
   likelihood that carries binaries and extinction covariance.
2. Preserve the current injections and rerun the response matrix against the
   revised recovered-mass estimator.
3. Reassess whether the P>0.5 WP3/WP4 handoff should be extended to the soft
   P>0.05 catalogue; the current injection correction is usable, but it cannot
   manufacture missing upstream photometry/mass posteriors.
4. Accept WP5 only after an untuned residual gate passes and the completeness
   deviation is either removed or explicitly approved as a corrected-selection
   design change.

## 7. Outputs

- `data/processed/wp5_completeness_curves.parquet`
- `data/processed/wp5_injection_response.parquet`
- `data/processed/wp5_imf_normalization.parquet` (diagnostic; blocked)
- `data/processed/wp5_mass_function_bins.parquet`
- `data/processed/wp5_association_mass.parquet`
- `data/processed/wp5_imf_posterior_draws.npz`
- `tables/wp5_imf_norm.csv`, `tables/wp5_imf_norm.md`
- `tables/wp5_association_mass.csv`
- `figures/wp5/wp5_completeness_curves.png`
- `figures/wp5/wp5_mass_function.png`
- `figures/wp5/wp5_association_mass.png`
- `notebooks/wp5_imf_normalization_and_completeness.ipynb`
- `provenance/wp5_manifest.json`, `provenance/wp5_provenance.md`
