# WP4 completion report — subgroup ages & per-star masses

Completed: 2026-07-23 · Gaia DR3 · distance posterior 1.6245 ± 0.045 kpc (μ=11.054)

- **WP5/WP7 input gate: READY WITH CONSTRAINTS** · **WP4 validation gate:
  SATISFIED** (one association-wide PASS with a subgroup-B limitation; one
  documented indicator disagreement)
- Runs on the WP2 members (P>0.5, **N=1392**) joined to the A/B/C subgroup labels
  (1,331 labelled; 61 anchor-exempt unlabelled). Isochrones **identical** to WP3.

## Headline result

The three kinematic subgroups remain **consistent with a coeval upper main
sequence**: retained upper-MS branch MAPs span **3.16–4.50 Myr**. The honest
envelope across both retained indicators is **2.25–5.67 Myr**, including the
younger subgroup-A PMS branches and one high-R_V subgroup-B PMS tail (N=19).
Thus ≈3.5–4.5 Myr is a central upper-MS summary, not the full WP4 envelope.
The A/B/C ordering is still unresolved, and the star-formation-duration branch
(0/1/2 Myr) remains mandatory for WP7.

<!-- BEGIN GENERATED:CLOSURE_SUMMARY -->
- Age branches: 104/132 measurable; 28 excluded; 9 grid-railed and excluded.
- Mass branches: all 6 retained; 55 common null source_ids; `mass_baseline` is the reporting-only PARSEC R_V=3.1 branch.
- HRD anchors A/B/C/unassigned: 60/5/43/48.
- WP6 proper-motion hand-off: 7 candidates, none removed from membership.
<!-- END GENERATED:CLOSURE_SUMMARY -->

## Deliverables

- `wp4_ages.md` — full results: per subgroup × family × R_V age posteriors, both
  indicators, credible intervals, distance sensitivity, gate, literature cross-check.
- `data/processed/wp4_masses.parquet` — per-star masses on 6 branches
  (`mass_{PARSEC,MIST}_rv{3.0,3.1,3.5}`) + `mass_method` provenance flag; also
  `tables/wp4_masses.cat`.
- `data/processed/wp4_age_posteriors.parquet` (132 branch rows) +
  `wp4_posterior_curves.npz`; `tables/wp4_ages_summary.csv`.
- `data/processed/wp4_anchor_hrd.parquet`, `data/processed/wp4_clump.parquet`.
- Figures `figures/wp4/`: CMD per subgroup + isochrones + binary loci + clump;
  spectroscopic HRD; age posteriors (both indicators); age money-plot.
- Provenance: `provenance/wp4_provenance.md`, `provenance/wp4_manifest.json`
  (SHA-256 of every input/script/log/output), one `wp4_*_execution.json` per step.

## Method (per plan WP4 + method_explained §4 + CUTS §1/§9)

- **Forward-model CMD likelihood.** Per trial age, synthesise the population's CMD
  density (IMF-weighted single stars + explicit unresolved-binary component,
  q~U[0.1,1], fluxes added), evaluate per-star marginal likelihood in (colour,
  M_G0) with propagated photometric + A_V-residual errors, sum weighted by
  membership P. Posterior with uniform prior in log-age.
- **Two independent indicators:** upper main sequence (M_G0 ≤ +1.5) and PMS
  turn-on (M_G0 ≥ +2.5), reported separately.
- **Branches carried, never averaged:** PARSEC/MIST × R_V 3.0/3.1/3.5 × f_bin
  0.3/0.4/0.5 × indicator. Distance fixed at 1.6245 kpc, ±0.045 kpc propagated at
  μ ± σ_μ.
- **Spectroscopic override:** 156 anchors placed on the (logTe, M_G0) HRD using
  spectroscopic Teff; masses from the isochrone; per-star provenance recorded.
- **Binaries handled in-model**, not σ-clipped; the WP3 high-A_V clump classified
  explicitly and shown not to drive the age.
- **Measurability is explicit:** N≥15, finite posterior, and no MAP within one
  native grid step of 1 or 10 Myr. Exclusion reasons are stored row by row.
- **Family treatment is symmetric:** each subgroup uses each family's own
  fitted upper-MS age. MIST is constant across anchors because all three fits
  selected the same 3.571-Myr grid age, not because it was fixed by hand.

## Gate results

| Criterion | Result |
|---|---|
| Age pattern internally consistent & literature-comparable; coeval stated if identical | **PASS** — upper-MS 3.16–4.50 Myr and A/B/C coeval within systematics; full two-indicator 2.25–5.67 Myr envelope stated |
| Anchor spectroscopic HRD consistent with isochrones within PARSEC–MIST differences | **PASS association-wide, LIMITED for B** — 131/150 (87%) within χ≤2.5, but subgroup B has only **N=5 anchors** |
| Both indicators agree where measurable, or disagreement documented | **DOCUMENTED DISAGREEMENT** — B agrees at baseline but reaches 5.67 Myr in a retained high-R_V PMS branch (N=19); A's offset documented; C PMS excluded |

## Cross-check vs literature (branch-matched)

- **Wright+15:** upper-MS (PARSEC 4.0–4.5, MIST 3.6) matches the **rotating** end
  (4–5 Myr); PMS/MIST (2.5–3.6) matches the **non-rotating** end (2–3 Myr). No tension.
- **Berlanas+20 bursts:** do **not** map onto A/B/C as spatially distinct ages
  (subgroups coeval within systematics); the burst structure instead echoes as a
  *temporal* spread (A's upper-MS/PMS offset). Carried as the SF-duration branch.

## Numbers for WP5/WP7

- Adopted upper-MS ages: **A ≈ 3.6–4.75, B ≈ 3.5–4.0, C ≈ 3.1–4.0 Myr** (family +
  branch envelope; PARSEC older, MIST younger).
- Per-star masses: **150 spectroscopic-HRD, 1,242 photometric-isochrone**; PARSEC
  0.6–40 M☉, MIST 0.4–51 M☉.
- `mass_baseline = mass_PARSEC_rv3.1` is reporting-only. WP5 must fit all six
  mass branches independently. The same 55 stars are null in every branch.
- **924 members in the 2–8 M☉ calibration window** (783 membership-weighted) — the
  WP5 IMF-normalisation sample (its lower edge to be *measured* by WP5, not adopted).

## Contingency flagged (plan WP4 caveat / WP9)

The retained two-indicator envelope reaches 2.25–5.67 Myr and the subgroup-B
high-R_V tail crosses the plan's ±1.5-Myr soft-verdict warning scale. WP4 does
not resolve that tension; WP9 must test verdict stability to it.

## WP5 hand-off constraints

1. Fit the IMF separately on all six mass branches; never fit only
   `mass_baseline`.
2. Run full-chain injection/recovery separately for A/B/C. Watch C's weak PMS
   constraint and B's five-anchor HRD coverage.
3. Measure the lower edge of the nominal 2–8 M☉ calibration window; do not
   assume it. Read the WP3 truncation warning and
   `figures/wp3/wp3_wp5_truncation.png` before implementation.

Every file and checksum is enumerated in `provenance/wp4_manifest.json`.
