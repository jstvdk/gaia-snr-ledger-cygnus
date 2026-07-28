# WP3 anchor-prior repair (repair_v3) — negative result

*Executed 2026-07-27. Tests and fixes the hypothesis in
`reports/WP5_PARENT_RANGE_FIX_repair_v2.md` for PROJECT_TRACE issue #1b.
repair_v1 and repair_v2 artifacts are untouched.*

## Verdict

**The diagnosis was confirmed. The fix did not work, and made the WP5 gate
worse.** WP5 remains BLOCKED; WP6 is not authorized.

Two separate conclusions, both true:

1. The repair_v1 anchor spatial prior **was** over-constraining CygOB2-B. That
   is established by three independent lines of evidence and is not in doubt.
2. Correcting it **does not** resolve the CygOB2-B mid-window excess. The
   residual grew from 3.02 to 3.62 and the branch grid fell from 33/54 to
   26/54. The extinction prior was not the cause; it was partly *masking* the
   real problem.

## Part 1 — the diagnosis, confirmed

`scripts/wp3_anchor_prior_diagnostic.py` →
`provenance/wp3_anchor_prior_diagnostic_execution.json`,
`tables/wp3_anchor_prior_diagnostic.csv`, `tables/wp3_anchor_variogram.csv`,
`figures/wp3_repair/wp3_anchor_prior_diagnostic.png`.

**(a) The prior width is applied far outside where it was calibrated.** The
repair_v1 prior gives each member the median A_V of its eight nearest
spectroscopic anchors, with one global width per R_V branch from leave-one-out
at the anchor density. The fitted anchor variogram (R_V = 3.1: nugget 0.000 mag,
sill 1.228 mag, range 0.853 deg) shows the field decorrelating well inside the
distance at which the prior is actually used:

| subgroup | median 8th-anchor separation | adopted σ | variogram-required σ | understated by |
|---|---:|---:|---:|---:|
| CygOB2-A | 0.089° | 0.453 | 0.635 | 1.40× |
| CygOB2-C | 0.130° | 0.453 | 0.744 | 1.64× |
| **CygOB2-B** | **0.377°** | 0.453 | **1.052** | **2.32×** |

The anchors' own median eighth-neighbour separation is 0.071°. CygOB2-B sits
5.3× further out, where the anchors carry almost no information, yet received
the same tight prior.

**(b) The prior collapses exactly one subgroup.** Against the frozen WP3
solution, which has no spatial prior:

| subgroup | frozen WP3 | repair_v1 | collapse |
|---|---:|---:|---:|
| CygOB2-A | 1.867 | 0.917 | 2.0× |
| CygOB2-C | 1.667 | 0.575 | 2.9× |
| **CygOB2-B** | **1.860** | **0.196** | **9.5×** |

B entered with essentially the same differential extinction as A (1.860 vs
1.867) and left 4.7× narrower.

**(c) Independent 3D dust cubes do not see B as uniform.** Vergely+22 and
Dharmawardena+22 have no knowledge of the anchors and rank B's spread at or
above A's (0.233 vs 0.195, and 0.078 vs 0.045). Nothing external supports a
uniformly-extinguished CygOB2-B.

## Part 2 — the fix

`ANCHOR_PRIOR_MODE = "variogram"` in `scripts/wp3_repair_common.py`: the prior
width becomes per-star, read off the fitted variogram at that star's own
eighth-nearest-anchor separation, floored at `MIN_PRIOR_SIGMA_MAG`. Stars inside
the anchor field keep a tight prior; stars extrapolated beyond the correlation
length correctly receive one near the sill, so their own photometry rather than
distant anchors sets their A_V.

Both knobs (`WP_REPAIR_VERSION`, `WP3_ANCHOR_PRIOR_MODE`) are environment
variables defaulting to the repair_v1 behaviour, so the whole chain re-runs at a
new version with no script edits and repair_v1 stays exactly reproducible.
Verified: in `global` mode `prior_sigma_at` returns the repair_v1 scalar for
every star.

Chain: `scripts/run_repair_v3_chain.sh` (WP3 → WP4 ages → WP4 masses → WP5
injections → WP5 fit).

**The extinction outcome is exactly as intended:**

| subgroup | frozen | repair_v1 | **repair_v3** |
|---|---:|---:|---:|
| CygOB2-A | 1.867 | 0.917 | 0.852 |
| CygOB2-B | 1.860 | **0.196** | **0.603** |
| CygOB2-C | 1.667 | 0.575 | 0.598 |

B's differential extinction is restored 3.1× while A and C move by −7% and +4%.
That asymmetry is the signature of a targeted correction rather than a global
loosening.

## Part 3 — but the WP5 gate got worse

| | repair_v2 | repair_v3 |
|---|---:|---:|
| branch grid passing | 33/54 | **26/54** |
| CygOB2-B baseline max abs residual | 3.02 | **3.62** |
| CygOB2-B baseline χ² p | 0.036 | **0.0002** |
| CygOB2-A / CygOB2-C baseline | pass / pass | pass / pass |

CygOB2-B's per-bin residuals sharpened into a clean monotone tilt:

| bin | 2.0–2.52 | 2.52–3.17 | 3.17–4.0 | 4.0–5.04 | 5.04–6.35 | 6.35–8.0 |
|---|---:|---:|---:|---:|---:|---:|
| repair_v2 | +1.04 | +1.11 | +3.02 | −0.41 | −0.51 | +0.25 |
| repair_v3 | +1.45 | +2.41 | **+3.62** | −0.36 | −1.30 | −0.92 |

B's observed mass function is steeper than the model across the whole window,
not high in one bin. The over-tight repair_v1 prior had been narrowing B's mass
distribution and partly hiding this.

## Part 4 — the completeness-ramp explanation is also ruled out

CygOB2-B and -C recover only 6.2% and 1.2% of injected 2.0 M☉ stars, so the
bottom of the window is reconstructed from almost nothing. That made a
per-subgroup lower edge — the procedure plan WP5 step 2 and CUTS §7.1 actually
prescribe — the obvious next candidate. It fails, and fails informatively.
`scripts/wp5_lower_edge_scan.py` →
`provenance/wp5_lower_edge_scan_execution.json`:

| rule | grid | A | B | C |
|---|---:|---|---|---|
| global 2.0 (current) | 26/54 | 2.00, r=1.37 ✓ | 2.00, **r=3.62** ✗ | 2.00, r=1.52 ✓ |
| absolute recovery ≥ 0.30 | 25/54 | 2.00, r=1.37 ✓ | 2.50, **r=5.11** ✗ | 2.75, r=1.42 ✓ |
| absolute recovery ≥ 0.50 | 28/54 | 2.25, r=1.70 ✓ | 2.75, **r=5.23** ✗ | 2.75, r=1.42 ✓ |
| 95% of bright plateau | 29/54 | 3.25, r=2.56 ✓ | 5.00, r=1.02 ✓ | 4.50, r=1.51 ✗ |

Raising B's lower edge makes B **monotonically worse** (3.62 → 5.11 → 5.23).
A completeness-ramp artifact would do the opposite — removing the unreliable
low-completeness bins would relieve it. Only at an edge of 5.0 M☉ does B pass,
and there it retains 82 sources and CygOB2-C fails instead.

## What is now excluded

| hypothesis | status |
|---|---|
| parent-range truncation at 8 M☉ | fixed in repair_v2; confirmed, but never affected B's bin 2 |
| anchor spatial prior over-constraining B | **confirmed real and fixed here — but not the cause** |
| completeness ramp / lower-edge choice | ruled out; the effect runs the wrong way |
| B's mass function is a single power law + this completeness model | **contradicted by the data** |

The excess is robust to every instrumental explanation tried. What survives is
that CygOB2-B genuinely has more 2.5–4 M☉ stars, relative to its 4–8 M☉ stars,
than a single-slope IMF predicts. Remaining candidates, none yet tested:

1. **Subgroup-label contamination.** `membership_probability` measures
   cluster-vs-field, never A-vs-B-vs-C. B may be absorbing lower-mass stars that
   belong to A or to the field. Cheap first check: label-stability confusion
   matrix over `provenance/wp2_gmm_seed_stability.csv`.
2. **A second population overlapping B in the CMD** — a foreground or background
   group at similar kinematics.
3. **B's age.** repair_v3 puts B's upper-MS MAP at 2.82 Myr against C's 2.51. An
   age error tilts the inferred mass function exactly this way, and WP4 already
   carries a documented indicator disagreement for subgroup B.

## Recommendation on which version to keep

**Keep repair_v3 as the working extinction, despite the worse gate score.** The
variogram argument stands entirely on its own — a width calibrated at 0.071° is
not valid at 0.377°, whatever it does downstream. Choosing repair_v2 because it
scores better on the residual gate would be selecting the extinction model by
its effect on the residuals, which is precisely the gate tuning the frozen
repair_v1 report warned against.

One honest caveat: with a ~1.05 mag prior width, CygOB2-B's A_V is now set
almost entirely by broadband photometry, which reintroduces the Teff/A_V
degeneracy the anchors were meant to break. repair_v3 may be trading one bias
for another for B specifically. Resolving that needs spectroscopic anchors
closer to B's footprint, which DR3 does not provide — a genuine limitation to
carry, and a concrete DR4/follow-up item.

This is a judgement call about the science, not about the pipeline, and it
should be made explicitly rather than inherited by default.

## Artifacts

| kind | path |
|---|---|
| diagnostic script | `scripts/wp3_anchor_prior_diagnostic.py` |
| diagnostic tables | `tables/wp3_anchor_prior_diagnostic.csv`, `tables/wp3_anchor_variogram.csv` |
| diagnostic figure | `figures/wp3_repair/wp3_anchor_prior_diagnostic.png` |
| lower-edge scan | `scripts/wp5_lower_edge_scan.py` |
| chain runner | `scripts/run_repair_v3_chain.sh` |
| repair_v3 data | `data/processed/{wp3_extinction,wp3_extinction_posterior,wp4_age_posteriors,wp4_mass_posteriors,wp4_mass_posterior_samples,wp5_completeness_curves,wp5_injection_response,wp5_imf_normalization,wp5_mass_function_bins,wp5_association_mass,wp5_imf_posterior_draws}_repair_v3.*` |
| provenance | `provenance/{wp3_anchor_prior_diagnostic,wp5_lower_edge_scan}_execution.json`, `provenance/*_repair_v3.json` |

Baseline association mass 30,696 M☉ multiplicity-adjusted, within a factor 2 of
the 16,500 M☉ literature scale — essentially unchanged from repair_v2's 30,453.
