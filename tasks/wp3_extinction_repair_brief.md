# Repair brief — WP3 extinction degeneracy corrupts WP4 masses and blocks the WP5 IMF gate

*Diagnosis date: 2026-07-24. Status of pipeline: WP5 BLOCKED at validation gate
(all 54 branch fits fail the mass-function residual gate; first observed-mass
bin off by 6.4–7.9σ). This document states the traced root cause, the evidence,
and the ordered fix. Hand this to the executing agent together with
`paper1_execution_plan.md`.*

---

## 1. Problem statement

The WP4 per-star masses do not form a believable IMF and cannot be repaired by
the WP5 completeness correction:

| mass bin (M☉, baseline branch, photometric stars) | observed N | α=2.3 expectation (same total in 2–8) |
|---|---|---|
| 2.0–2.52 | **454** | ~277 |
| 2.52–3.17 | 130 | ~204 |
| 3.17–5.04 | **170** | ~265 |
| 5.04–8.0 | 135 | ~144 |

An **excess** at 2–2.5 M☉ coexists with a **deficit** at 2.5–5 M☉. Membership
contamination cannot produce a deficit; **mass migration can and does**. The
massive end is also affected: part of the 8–12 M☉ census (118 stars) is inflow
from 5–8 M☉, so the current count of 334 stars > 8 M☉ is biased and WP6's
closure test must not run before the fix.

## 2. Root cause (traced, verified on the frozen data)

**`scripts/wp3_broadband_extinction.py` produces per-star A_V estimates that
are ill-conditioned along the reddening/PMS degeneracy, and the error it
reports (~0.05–0.07 mag) hides a ±2 mag bimodal ambiguity.** Chain:

1. **The fit is effectively optical-only.** WLS weights are `1/err²` with raw
   catalogue errors and **no calibration/model error floor** (line ~89). Gaia
   G/BP/RP errors are milli-mag; 2MASS J/H/Ks errors are ~0.02–0.03 mag, so the
   IR bands — the intended degeneracy breakers — carry ~100× less weight.
   (100% of the 2–8 M☉ window stars *have* 2MASS; it just never gets a vote.)
2. **The optical-only solution is degenerate.** At the fitted ~4 Myr age the
   isochrone PMS arc (2.0–2.5 M☉) has CMD slope ~3.3 mag/mag, nearly parallel
   to the reddening vector (~2.05). A true 3–5 M☉ ZAMS star at A_V≈6.5 refits
   almost equally well as a ~2.2 M☉ PMS star at A_V≈4.6. The χ²-posterior over
   templates (with underestimated errors) collapses to a single branch chosen
   by ~0.05-mag photometric details.
3. **WP4 nearest-point masses inherit the branch flip** (`scripts/wp4_masses.py`,
   nearest isochrone point in de-reddened CMD): A_V and mass slide together —
   down (3–5 → 2–2.5 pileup) and up (5–8 → 8–12 inflow).
4. **WP5 injections cannot see this error.** `scripts/wp5_injections.py`
   (line ~348) models the extinction error as `av_true + N(0, av_err)` with
   σ≈0.05 mag — a nearly diagonal response matrix — while the real estimator
   error is a bimodal ±1.5–2 mag branch flip. Hence the unfixable 6–8σ
   first-bin residual in every branch. The WP5 gate fired correctly.

## 3. Key evidence (reproduce before and after the fix)

- **Sightline test (the smoking gun).** Fitted A_V minus the median A_V of the
  8 nearest spectroscopic anchors (149 stars, `av_method != broadband_multiband`
  in `data/processed/wp3_extinction.parquet`), per assigned-mass bin:
  **−1.46 mag** at 2–2.52 M☉; ~0 at 3.2–5 M☉; **+1.5…+1.8 mag** at 6.4–12 M☉.
  Extinction is a sightline property; any correlation with assigned mass is
  unphysical and measures the migration.
- The 2–2.52 M☉ bin's A_V distribution truncates at ~6.0 while anchors on the
  same sightlines reach 7.2. A real 2.3 M☉ PMS star at A_V 6.5 would have
  G≈17.5 (well above the G<19 limit) — such stars are absent because the
  "PMS" population is largely the low-A_V misfit branch of true MS stars.
- Membership is NOT the driver: pileup stars have median P=0.97, RUWE 1.04,
  parallax 0.616 mas (on-cluster).
- Grid snapping (secondary defect): 60 stars sit at exactly mass 2.547, 35 at
  8.001, etc. — nearest-point assignment to discrete isochrone rows.
- **Remediation preview passes.** Re-deredden with neighborhood-anchor A_V and
  re-match masses: the first-bin fraction of the 2–8 window goes from 0.511 to
  **0.318 vs 0.311 expected (α=2.3)**, and the 3.2–5 range refills. The
  residual 2.5–3.2 dip is the grid-snapping artifact, removed by step F2 below.

## 4. Fix — ordered steps

**F1. Reopen WP3 broadband extinction (root fix).**
   - Add an error floor (~0.02–0.03 mag) in quadrature to *all* band errors so
     2MASS genuinely constrains the template choice.
   - Add a spatial extinction prior built from the 149 spectroscopic anchors
     (e.g. per-star A_V ~ N(local anchor-map value, σ_patchiness), σ estimated
     from the anchor variogram; keep it wide enough to respect the real
     Cygnus differential extinction — this is a prior, not a replacement).
   - Output a **full per-star A_V posterior** (explicitly multimodal where the
     branches persist), not a point value with a 0.05-mag error.

**F2. Reopen WP4 lower/intermediate-mass inference.**
   - Replace nearest-point masses with **per-star mass posteriors**
     marginalizing over: the F1 A_V posterior, photometric errors, the
     unresolved-binary branch, and the subgroup age posterior.
   - This also removes grid snapping and the artificial 2.5–3.2 M☉ hole.
   - Keep the spectroscopic-HRD override for anchors unchanged.

**F3. Rerun WP5 with an honest response matrix.**
   - Preserve the existing injections; change one thing: pass injected
     photometry through the **actual F1 estimator** (fit A_V on the injected
     magnitudes) instead of `av_true + N(0, av_err)`, then through the F2 mass
     estimator. The response matrix then contains the true migration kernel.
   - Do NOT relax the WP5 residual gate; it was correct.

**F4. Add a permanent WP3 gate (cheap invariant that would have caught this).**
   - Residual of fitted A_V against the local anchor map must be independent
     of downstream assigned mass: |median Δ A_V per mass bin| < ~0.3 mag and
     rank-correlation(ΔA_V, mass) consistent with 0.

**F5. Then, and only then:** re-run WP6 closure (the >8 M☉ census will change)
   and proceed to WP7.

## 5. Provenance discipline

This reopens frozen WP3/WP4 products. Follow the WP4-schema-repair precedent:
version the repair (do not overwrite artifacts in place), keep the old files
and the WP5 blocking verdict in the record as the trigger, log every new
threshold (error floor, prior width) with its class in
`CUTS_AND_THRESHOLDS.md`, and write `wp3_repair_execution.json` /
updated manifests. The WP5 gate history must remain visible: the paper's
auditability claim depends on showing the gate fired and what it caught.

## 6. Acceptance criteria

1. WP3 gate F4 passes (mass-independent A_V residuals vs anchor map).
2. WP4 mass posteriors: no single mass value holds >2% of the catalogue; the
   2.5–3.2 M☉ region is populated.
3. WP5: at least the baseline branch passes the untuned Poisson residual gate;
   completeness statement re-derived (the 95% absolute edge question from the
   current WP5 report must be re-evaluated, not assumed fixed).
4. Association-mass sanity check still within factor ~2 of the Wright+15
   1.65×10⁴ M☉ scale.
5. N(>8 M☉) re-reported with the migration removed, before WP6 consumes it.
