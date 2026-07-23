# WP3 provenance log — per-star extinction & de-reddened CMDs

*Running log of every cut, threshold, artifact, and query for WP3. Companion to
`CUTS_AND_THRESHOLDS.md` (number classes) and the master manifest
`provenance/wp3_manifest.json` (sha256 of every input/output). Machine-readable
step logs are the `provenance/wp3_*_execution.json` files referenced below.*

**Date:** 2026-07-23 · **Data release:** Gaia DR3 · **Distance posterior:**
1.6245 ± 0.045 kpc (single population; `wp2_membership_manifest.json`,
`one_component_mean_kpc`) · μ = 11.05 mag.

Runs on the **full WP2 member list** (P>0.5, N=1392); no subgroups are used
(the two-distance-population claim was not confirmed in WP2).

---

## Step 1 — Isochrone grids (shared with WP4)

- **Scripts:** `wp3_fetch_parsec.py`, `wp3_fetch_mist.py`.
  **Logs:** `wp3_parsec_fetch_execution.json`, `wp3_mist_fetch_execution.json`.
  **README:** `wp3_isochrones_README.md`.
- PARSEC v1.2S (CMD 3.9) + MIST v1.2, solar metallicity, 1–10 Myr log grid
  (21 ages), Gaia EDR3 + 2MASS JHK, intrinsic (Av=0). Outputs
  `wp3_isochrones_parsec.parquet` (6020 rows), `wp3_isochrones_mist.parquet`
  (11213 rows). Generated once to prevent a WP3/WP4 inconsistency.

## Step 2 — Extinction law (R_V branches)

- **Script:** `wp3_extinction_law.py`. Cardelli+89 with O'Donnell 1994 optical
  update, R_V-dependent (Class B+E, cut #17). Band coefficients k_X=A_X/A_V at
  the band effective wavelengths for **R_V = 3.0 / 3.1 (baseline) / 3.5**, all
  three run from the start (retrofitting branches is painful).
- **Cross-check:** matches the CMD 3.9 service (Cardelli+O'Donnell, R_V=3.1,
  G2V) to 0.0% for G/BP/RP: G 0.8363, BP 1.0834, RP 0.6344.
- k_X (R_V=3.1): G 0.836, BP 1.083, RP 0.634, J 0.288, H 0.178, Ks 0.117.

## Step 3 — Member photometry assembly

- **Script:** `wp3_assemble_photometry.py`. **Log:** `wp3_assemble_photometry_execution.json`.
- 1392 members joined to Gaia DR3 G/BP/RP + 2MASS JHK. Gaia mag errors from
  flux errors via DR3 VEGAMAG zeropoints, 0.003 mag floor. 2MASS: `*_msigcom`.
- Coverage: **1336** with a 2MASS PSC match, **1261** complete JHK, **1239**
  ph_qual AAA, **56 members without any 2MASS match are flagged (not dropped)**,
  9 without Gaia photometry (bright saturated anchors carried from WP1).

## Step 4 — Anchor extinction from intrinsic colours (NOT broadband)

- **Script:** `wp3_anchor_extinction.py`. **Log:** `wp3_anchor_extinction_execution.json`.
- Rationale (plan WP3 caveat; CUTS §5): for hot OB stars broadband colours
  saturate and the Teff/extinction degeneracy is real. For the anchors the
  intrinsic SED is **fixed** from the spectroscopic Teff (or spectral type) and
  only the reddening is solved:  m_X = M_X0(Teff) + μ + k_X·A_V (two linear
  parameters μ, A_V by WLS over available bands). A wrong luminosity class is
  absorbed by μ, so A_V is set by the logg-robust colour pattern.
- M_X0(Teff) from the WP3 isochrone MS locus (PARSEC+MIST averaged; family
  colour spread recorded). Teff units normalised (anchor column mixes K and kK).
  Spectral-type→Teff fallback for anchors without a tabulated Teff.
- **WR stars (3) flagged and excluded** from the colour method (emission-line,
  non-photospheric colours); fall back to literature A_V.
- Result: 190 anchors P>0.5 → **149 via intrinsic colour**, 3 WR, 41 literature
  fallback (no usable Teff/photometry). A_V(R_V=3.1) median 5.65 mag, [4.3,7.5]
  5–95%, **0 negative**.
- **Validation vs literature (Berlanas) A_V:** N=145, median residual **+0.12
  mag**, rms **0.33 mag** — the intrinsic-colour method reproduces the
  spectroscopic literature extinctions.

## Step 5 — Broadband multiband extinction (all other members)

- **Script:** `wp3_broadband_extinction.py`. **Log:** `wp3_broadband_extinction_execution.json`.
- m_X = M_X0(template) + μ + k_X·A_V fit against the union of PARSEC+MIST
  isochrone templates (17233 points, ages 1–10 Myr). Per template A_V is the
  closed-form WLS optimum; the estimate is the χ²-posterior-weighted mean over
  templates, weighted-std as statistical error. Run per R_V branch. IR bands
  break the Teff/A_V degeneracy.
- **56 members without 2MASS:** errors inflated (×1.8 on the bands + 0.5 mag
  added in quadrature on the reported A_V), not dropped (plan WP3 step 2).
- Result (non-anchor, R_V=3.1): median A_V 5.88, [2.1,8.6] 5–95%, max 10.6;
  **2 marginally negative, 0 beyond 2σ**.

## Step 6 — Catalogue, absolute magnitudes, de-reddened CMDs

- **Script:** `wp3_build_catalogue.py`. **Log:** `wp3_build_catalogue_execution.json`.
- Per-member A_V chosen: anchors → intrinsic colour (else literature); all
  others → broadband. Per R_V branch: A_V, A_V_err, A_G=k_G·A_V, and
  de-reddened M_G0, (BP−RP)0, (G−Ks)0.
- **Distance:** group posterior μ=11.05 for all, **overridden by individual
  parallax for the 33 brightest well-measured members** (G<11 and parallax
  S/N>10) — defensible where the parallax is good.
- Outputs: `wp3_extinction.parquet` (all branches) and per-branch ECSV
  catalogues `tables/wp3_extinction_rv{30,31,35}.cat` (**the wp3_extinction.cat
  deliverable, one per R_V branch**).
- A_V medians: R_V 3.0 → 5.77, 3.1 → 5.85, 3.5 → 6.25 (all in the 4–8 core band).

## Step 7 — Cube comparison + spatial coherence

- **Scripts:** `wp3_cube_comparison.py`, `wp3_figures_and_coherence.py`.
  **Logs:** `wp3_cube_comparison_execution.json`, `wp3_figures_coherence_execution.json`.
- **Vergely+22** (EXPLORE, A0/pc density, 25 pc): line-of-sight integral to the
  group distance. Significant **positive** spatial correlation with per-star
  A_V (Spearman ρ=0.29, p=3e-28); per-star A_V ~1.6× the cube (star 5.85 vs cube
  3.65), as expected since the 25 pc reconstruction under-resolves the clumpy
  dust local to Cyg OB2.
- **Dharmawardena+22** (l,b,d cumulative, 0.14°): sits far below (median 1.31,
  ratio 4.5) with weak negative ρ over the compact member footprint — its
  smoothed extinction peak is offset onto the molecular clouds (l≈78.7, b≈−0.75)
  rather than the relatively dust-cleared stellar core; its dynamic range over
  the ~1° footprint is too small to track per-star structure. Reported honestly;
  cubes are a cross-check, not a substitute (plan WP3 caveat).
- **Direct nearest-neighbour coherence** (the plan's actual requirement): each
  star's A_V vs the mean A_V of its 8 nearest sky neighbours → Spearman
  **ρ=0.66 (p=2e-167)**, vs a label-permutation baseline of 0.00±0.036 →
  **18.3σ**. Spatial coherence decisively demonstrated.

## Step 8 — WP5 truncation-bias note (anticipated here, not in WP5)

- The member sample is **extinction-biased by the G<19 limit**. A star of mass m
  drops below G=19 at 1.62 kpc above A_V_edge = (19 − M_G0(m) − μ)/k_G:
  **2 M☉ lost above A_V≈6.5, 3 M☉ above ≈8.7, 5 M☉ above ≈10.1 mag** (4 Myr
  PARSEC). The recovered A_V distribution is therefore **truncated at the high
  end for low-mass stars**: the high-A_V members (up to 10.6) are necessarily
  the more luminous stars, visible in the CMD as the high-A_V points piling up at
  the luminous end. WP5 completeness must be computed as a function of *both*
  mass and position (extinction is patchy), per CUTS §4.3 / §7.1. Figure:
  `figures/wp3/wp3_wp5_truncation.png`.

---

## Cut/threshold register (WP3 additions to CUTS_AND_THRESHOLDS.md)

| Cut | Value | Class | Justification |
|---|---|---|---|
| Extinction law | CCM89+O'Donnell94 | B | matches CMD service; R_V-dependent |
| R_V branches | 3.0 / 3.1 / 3.5 | B+E | carried in parallel from the start (#17) |
| Isochrone families | PARSEC v1.2S, MIST v1.2 | E | both carried, never averaged (#18) |
| Isochrone age grid | 1–10 Myr, Δlog=0.05 | A | fine grid (#19) |
| Anchor A_V method | fixed-Teff intrinsic colour | A | avoids Teff/A_V degeneracy for hot stars |
| Broadband A_V method | χ²-posterior over PARSEC+MIST templates, μ fixed | A/C | IR breaks degeneracy |
| No-2MASS error inflation | ×1.8 bands + 0.5 mag on A_V | D | optical-only penalty, flagged not dropped |
| Bright-star distance | G<11 & parallax S/N>10 → individual μ | A | defensible individual parallax |
| Gaia mag error floor | 0.003 mag | B | calibration systematic |

## Gate assessment — PASSED

- A_V distribution in the literature 4–20 mag range: **yes** (median 5.85,
  span ~2–10.6, core 4–8 populated).
- No unphysical negative extinctions beyond noise: **yes** (2 marginally
  negative, 0 beyond 2σ).
- Spatial coherence demonstrated: **yes** (neighbour ρ=0.66 at 18σ; Vergely
  ρ=0.29).
- Comparison vs Vergely+22 and Dharmawardena+22 at member distances: **done**
  (documented, including the honest Dharmawardena offset).
