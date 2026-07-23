# WP3 completion report — per-star extinction & de-reddened CMDs

Completed: 2026-07-23 · Gaia DR3 · distance posterior 1.6245 ± 0.045 kpc (μ=11.05)

- **WP4 input gate: READY** · **WP3 validation gate: PASS**
- Runs on the full WP2 member list (P>0.5, **N=1392**); no subgroups (two-distance
  population unconfirmed in WP2).

## Deliverables

- `tables/wp3_extinction_rv{30,31,35}.cat` — per-star extinction catalogue, one
  ECSV per R_V branch (the `wp3_extinction.cat` deliverable).
- `data/processed/wp3_extinction.parquet` — all branches + de-reddened absolute
  mags/colours in one table.
- `data/processed/wp3_isochrones_{parsec,mist}.parquet` — **shared with WP4**
  (Gaia G/BP/RP + 2MASS JHK, solar Z, 1–10 Myr; see `provenance/wp3_isochrones_README.md`).
- Figures `figures/wp3/`: de-reddened CMDs, extinction sky map, A_V distribution
  + cube comparison, WP5 truncation edge.
- Provenance: `provenance/wp3_provenance.md` (running log),
  `provenance/wp3_manifest.json` (sha256 of every input/output), and one
  `wp3_*_execution.json` per step.

## Method (per plan WP3 + CUTS §4.3/§5/§7)

- **Isochrones first** (PARSEC v1.2S + MIST v1.2), reddened in-code with a
  **R_V-dependent CCM89+O'Donnell law**; R_V = 3.0 / 3.1 / 3.5 branches from the
  start. Law verified to 0.0% against the CMD 3.9 service.
- **Extinction split by star type:** 190 spectroscopic anchors get A_V from
  **intrinsic colours at fixed spectroscopic Teff** (149 of them; 3 WR flagged,
  rest literature) — avoiding the Teff/extinction degeneracy for hot stars.
  All other members: **broadband multiband fit** against the isochrone templates.
- **56 members without 2MASS flagged with inflated errors, not dropped.**
- **Absolute magnitudes** from the 1.62 kpc posterior; the 33 brightest
  well-measured members use their individual parallax.

## Gate results

| Criterion | Result |
|---|---|
| A_V in literature 4–20 mag range | **PASS** — median 5.85 (R_V=3.1), span ~2–10.6, core 4–8 populated |
| No unphysical negatives beyond noise | **PASS** — 2 marginally <0, 0 beyond 2σ |
| Spatial coherence demonstrated | **PASS** — nearest-neighbour A_V Spearman ρ=0.66 (18.3σ over permutation) |
| Comparison vs Vergely+22 / Dharmawardena+22 | **DONE** — Vergely ρ=0.29 (p=3e-28), ~1.6× above cube (clumpy dust under-resolved); Dharmawardena offset documented honestly |
| Anchor A_V vs literature (Berlanas) | median residual +0.12 mag, rms 0.33 mag |

## WP5 warning carried forward (do not rediscover as a completeness bug)

The member sample is **extinction-biased by G<19**. At 1.62 kpc a star drops
below the limit above A_V ≈ **6.5 (2 M☉)**, **8.7 (3 M☉)**, **10.1 (5 M☉)**
(4 Myr PARSEC). The recovered A_V distribution is therefore truncated at the
high end for low-mass stars, and completeness must be modelled as a function of
**both mass and position** (extinction is patchy). See
`figures/wp3/wp3_wp5_truncation.png` and `provenance/wp3_provenance.md` §8.

Every file and checksum is enumerated in `provenance/wp3_manifest.json`.
