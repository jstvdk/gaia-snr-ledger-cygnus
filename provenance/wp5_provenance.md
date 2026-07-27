# WP5 provenance — completeness and IMF normalization

**Generated:** 2026-07-23T14:43:27.459596+00:00  
**Status:** `WP5_BLOCKED_RESIDUAL_AND_ABSOLUTE_95_COMPLETENESS_GATES`  
**Downstream authority:** none; diagnostic products must not enter WP6/WP7.

## Pipeline

1. `scripts/wp5_common.py` reconstructs the frozen WP2 classifier, loads the
   branch-matched WP3/WP4 models, and defines the IMF/multiplicity integrals.
2. `scripts/wp5_injections.py` injects 400 stars at each 0.25-Msun point from
   0.5--8 Msun for every subgroup x family x R_V branch, clones real Gaia
   observational states, runs the frozen selection, and records recovered
   masses.
3. `scripts/wp5_fit_imf.py` fits all subgroup x family x R_V x alpha branches
   with a response-aware Poisson likelihood and 10,000 posterior draws.
4. `scripts/wp5_report.py` creates the figures, tables, and gate report.
5. `scripts/wp5_make_notebook.py` creates the explanatory notebook, which is
   executed in the `cygob2-gaia` environment.
6. `scripts/wp5_finalize.py` validates schemas/hashes and freezes the manifest.

## Frozen execution constants

| Quantity | Value | Class / role |
|---|---:|---|
| Injection seed | 20260723 | reproducibility |
| Mass grid | 0.5--8 Msun, step 0.25 | response support |
| Injections per mass/branch/subgroup | 400 | binomial precision |
| Binary fraction | 0.4 | central WP4 binary branch |
| q distribution | U[0.1,1] | WP4 convention |
| Membership integration | 128-point Sobol normal | validated approximation |
| Published handoff | P>0.5 | frozen WP3/WP4 selection |
| Absolute completeness target | 0.95 | required gate; failed |
| Observed calibration window | 2--8 Msun | corrected diagnostic fallback |
| Observed mass bins | 6 logarithmic | residual diagnostic |
| IMF slopes | [2.0, 2.3, 2.6] | mandatory branches |
| Posterior draws | 10000 | k/response uncertainty |
| Total-mass range | 0.08--120 Msun | Kroupa-like integration |

## Injection validation

The QMC approximation was compared directly with the frozen WP2 10,000-draw
probabilities: decision agreement =
0.98947
over 760 stratified
sources.  The injection is catalogue-level, not image-level: Gaia epoch images
and AGIS cannot be rerun.  That limitation is explicit in the execution log.

## Selection definition

Recovery requires all of:

1. synthetic observed G<19 and raw parallax in 0.35--1.10 mas;
2. exact WP2 quality state (RUWE<1.4, >=8 visibility periods, BP/RP excess
   present, positive finite covariance, Lindegren zero-point domain);
3. WP2 posterior-odds membership probability >0.5;
4. G/BP/RP present for a WP4 recovered mass.

The mass response includes photometric error, A_V error, unresolved binaries,
and the exact WP4 nearest-isochrone inverse mapper.  True masses down to 0.5
Msun are included so upward migration into the observed 2--8 Msun window is not
omitted.

## Statistical model

Observed bin counts are sums of membership probabilities.  The forward
intensity uses the injected true-to-recovered response.  A Jeffreys prior is
used for k; response rows receive a Jeffreys-multinomial Dirichlet posterior.
Primary-system mass uses a continuous two-part IMF (alpha=1.3 below 0.5 Msun;
branch alpha above), and the separately reported stellar mass adds companions
under the declared f_bin/q convention.

## Gate failure preserved

No absolute 95% completeness edge exists.  All 54 response-aware fits fail the
Poisson residual gate, including every mandatory slope/family/R_V branch.  No
threshold, bin edge, or branch was tuned after observing the failure.  The
diagnostic association mass is retained because it helps localize the failure,
but it has no downstream authority.
