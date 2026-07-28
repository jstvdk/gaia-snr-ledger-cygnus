# repair_v5 — the kriged anchor prior, and the resolution of CygOB2-B

*Adopted 2026-07-28. Decision record:
[wp3_kriging_adoption.json](../provenance/wp3_kriging_adoption.json).
This report supersedes the diagnosis in
[WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md](WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md)
(whose fix was correct but incomplete) and completes the story begun in
[WP5_AGE_CONDITIONAL_SCAN_repair_v3.md](WP5_AGE_CONDITIONAL_SCAN_repair_v3.md).*

**Adopting this method is not the same as passing the WP5 gate.** WP5 remains
unaccepted and `downstream_wp6_authorized` stays `false`; see §6.

---

## 1. The defect

The WP3 extinction fit combines each star's photometry with a spatial prior
built from the 149 spectroscopic anchors. Two pieces of that prior contradicted
each other:

- `prior_sigma_at` evaluated the fitted variogram at each star's own
  eighth-nearest-anchor separation, so a star far from anchors correctly
  received a width near the sill — a formal statement that those anchors carry
  almost no information about it.
- `evaluate` nevertheless returned the plain **median of the same eight
  anchors**, giving them weights summing to one no matter how far away they
  were.

So the prior simultaneously declared the anchors uninformative and centred
itself on them at full strength. Its first and second moments disagreed.

This matters unevenly across the association, because anchor coverage is
unequal:

| subgroup | 8th-nearest anchor | anchors matched inside the subgroup | prior width |
|---|---:|---:|---:|
| CygOB2-A | 0.089° | 59 | 0.637 mag |
| **CygOB2-B** | **0.377°** | **4** | 1.050 mag |
| CygOB2-C | 0.130° | 42 | 0.764 mag |

CygOB2-B's prior was extrapolated from anchors nearly beyond the fitted
correlation range (0.853°), and it sat **+1.024 mag above the global anchor
median**, pulling B's adopted A_V **+0.485 mag above** what B's own photometry
preferred.

## 2. Why the prior, not the photometry, was the unreliable side

The disagreement alone proves nothing — the prior exists precisely because
broadband photometry is degenerate. Four independent lines settle it:

1. **The misfit is a pure A_V offset.** Per-band residuals divided by each
   band's extinction coefficient are constant to <0.1 mag across all six bands
   (B: −0.43, −0.49, −0.47, −0.51, −0.54, −0.45). The model has the colours
   right and only the extinction wrong.
2. **No optical/near-IR conflict.** A_V from G/BP/RP alone and from J/H/Ks
   alone agree to 0.04–0.07 mag in every subgroup, and B is not the outlier
   (B +0.044 vs A +0.074). *(This retracted an earlier hypothesis of mine that
   the optical bands were outvoting the near-IR ones.)*
3. **The error grows with anchor distance — only in B.** Correlating
   (prior − photometry) against anchor separation within each subgroup gives
   ρ = **+0.42 (p = 1×10⁻¹⁸)** for B, against −0.13 for A and −0.08 (n.s.) for
   C. Anchor separation is a property of sky position and carries no
   photometric information, so a correlation implicates the prior.
4. **Independent dust maps contradict the prior's ranking.** Both 3D cubes
   place B *lowest*: Vergely 3.378 (B) vs 3.685 (A) and 3.763 (C);
   Dharmawardena 1.200 vs 1.371 and 1.278. The uncorrected prior ranked B
   highest.

## 3. The fix

Replace the neighbour median with a **simple-kriging** estimate using the
variogram **already fitted in repair_v3** (nugget ≈ 0, sill 1.228 mag, range
0.853°). Weights fall off with separation and the deficit is taken up by the
field mean, so the estimate degrades gracefully to "no local information"
instead of asserting a distant measurement. The prior **width is deliberately
left unchanged** — this corrects where the prior sits, not how confident it is.

**No new parameter is introduced**, and the correction's size follows from
anchor geometry rather than from any choice:

| subgroup | Σ kriging weights | prior-mean shift |
|---|---:|---:|
| A | 1.000 | +0.010 |
| **B** | **0.772** | **−0.512** |
| C | 0.992 | −0.110 |

`ANCHOR_PRIOR_MODE=kriging`; `variogram` and `global` remain reachable and
repair_v3's behaviour is bit-preserved.

## 4. Pre-registered predictions, all confirmed

Written to [wp3_kriging_prior_prereg.json](../provenance/wp3_kriging_prior_prereg.json)
**before** the chain was run; outcomes in
[wp3_kriging_prior_outcome.json](../provenance/wp3_kriging_prior_outcome.json).

| | prediction | measured | verdict |
|---|---|---|---|
| P1 | B's A_V falls 0.2–0.5 mag | 6.391 → 6.031 (**−0.359**) | ✓ |
| P2 | B's age rises above 2.818 Myr | **2.818 → 3.548**, 68% [3.389, 3.991] | ✓ |
| P3 | A and C essentially unchanged | ages moved **0.000 Myr**; A_V −0.051, −0.059 | ✓ |
| P4 | B's tilt flattens, \|T\| < 2 | \|T\| 2.40 → **1.31**, trend p 0.017 → 0.186 | ✓ |

P2 deserves emphasis: B's new age interval reaches **3.981 Myr**, the age at
which WP5's mass function had *independently* passed its gate in
[wp5_tilt_vs_age_diagnostic_execution.json](../provenance/wp5_tilt_vs_age_diagnostic_execution.json).
The WP4/WP5 age conflict closed without either being tuned toward the other.

## 5. Result

CygOB2-B's baseline residuals across the three versions:

| version | bin 0 | bin 1 | bin 2 | bin 3 | bin 4 | bin 5 | T | gate |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| repair_v3 | +1.45 | +2.41 | **+3.62** | −0.36 | −1.30 | −0.92 | −3.22 | fail |
| repair_v4 | +1.20 | +2.08 | **+2.51** | −0.24 | −0.79 | −0.54 | −2.40 | fail |
| **repair_v5** | +0.87 | +2.27 | **+0.36** | −0.24 | −0.03 | +0.17 | **−1.31** | **pass** |

This closes finding **F1** of the original fix brief, which had decomposed B's
anomaly into *a localized bump plus a global mass-scale displacement*. The bump
was the PMS/Henyey isochrone fold under a single-age injection truth model
(fixed in repair_v4); the displacement was this extinction error.

**repair_v5 is the first version whose baseline passes for all three subgroups
under both the incumbent and the better-calibrated replacement trend
statistic.** Branch grid **26 → 38 of 54** (37 under the replacement); of the 13
cells that gained, 11 are CygOB2-B.

## 6. What this does *not* settle

Under the strict per-branch reading of gate G3 adopted on 2026-07-27, three
A/C non-baseline cells regress against repair_v4:

- CygOB2-A PARSEC R_V=3.1 α=2.0 — trend p **0.040**
- CygOB2-A MIST R_V=3.5 α=2.6 — trend p **0.048**
- CygOB2-C MIST R_V=3.5 α=2.0 — max\|r\| **3.25** (a genuine failure)

The first two sit inside the indeterminate band around the 0.05 threshold that
**issue #11** identified and that criterion R3 has not yet resolved. WP5 is
therefore still not accepted and WP6 is not authorized. What blocks WP5 is now
a threshold-methodology question, not an astrophysical one.

## 7. Downstream consequences and obligations

Measured in [wp3_kriging_downstream_impact.json](../provenance/wp3_kriging_downstream_impact.json).

**Magnitude effect — small.** k moves +2.2% (A), **−12.2% (B)**, +1.0% (C);
association mass 30,155 → 29,185 M☉, **−3.2%**. That is far inside the ±18%
spread already carried across the 18 model branches, so N_SN shifts ~3%.

**Structural effect — this is the consequential one.** The age spread is
unchanged (1.469 Myr, since A and C set the endpoints), but B relocates from
near C to near A. The association now reads as **two older subgroups (A 3.98,
B 3.55) plus one younger (C 2.51)** rather than one older plus two younger.
That is a different star-formation history and WP7's supernova timeline depends
on it.

**Selection bias — tested and absent.** Kriging shrinks poorly-anchored stars
toward the anchor median, so an anchor sample biased toward low-extinction
sightlines would have been imported into exactly the stars this fix touches
most. Comparing each anchor against the median photometry-only A_V of its 15
nearest members (median separation 0.043°) gives **+0.500 mag**
(Wilcoxon p = 4.3×10⁻¹⁶) — anchors read *higher*, the opposite sign to the
failure mode.

**Obligations created by adoption** (tracked in the decision record):

- **O1** — report the revised star-formation history as a result, not a
  footnote.
- **O2** — carry the anchor absolute-scale systematic: broadband photometry
  sits ~0.5 mag below spectroscopically calibrated anchors at matched position,
  so the absolute mass scale and N_SN inherit the anchor calibration. Predates
  repair_v5.
- **O3** — carry B's calibration asymmetry: 4 anchors against A's 59 and C's
  42. Gaia XP spectra are the natural remedy (issue #1d).
- **O4** — resolve issue #11's R3 criterion before WP5 can be accepted.

## 8. Reproduction

```bash
bash scripts/run_repair_v5_chain.sh
```

Every repair_v1–v4 artifact is preserved byte-identical; the older prior modes
remain reachable via `WP3_ANCHOR_PRIOR_MODE=variogram|global`.
