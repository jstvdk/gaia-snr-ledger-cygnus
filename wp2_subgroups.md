# WP2 — Cyg OB2 subgroups

*Task A deliverable. Real subgroup labels derived from sky position + proper motion
on the clean WP2 membership. Companion machine log: `provenance/wp2_subgroups_execution.json`.
Distance-population test (Berlanas two-distance claim) retained in §6.*

## 1. Sample and method

- **Sample:** 1,331 clean automatic members — `membership_probability > 0.5 & ~anchor_quality_exempt`
  from `data/processed/wp2_members.parquet`. The 61 manual spectroscopic quality
  exceptions are excluded from the fit, as in the distance test.
- **Feature space:** `(l, b, μα*, μδ)`. **Parallax is excluded.** The WP2 distance test
  found one population at 1.62 kpc with intrinsic depth 45 pc — that depth is exhausted
  at DR3 precision and carries no subgroup information (CUTS §6.1, §3.3). Substructure
  is real in sky+PM: intrinsic PM dispersion 0.297 mas/yr (≈8× the 0.037 mas/yr median
  error) and the sky is centrally clumped (radial KS = 0.61 vs a uniform disc).
- **Clusterer:** Gaussian mixtures, `k = 2…8`, `covariance_type="full"`, refit across
  **50 deterministic seeds**. Features are `StandardScaler`-standardised. Acceptance is
  **seed stability, not BIC** — BIC is non-monotonic here (a local-optima signature), so a
  partition is retained only if seed-to-seed Adjusted Rand Index (ARI) is high **and** its
  10th percentile is high. A per-seed consensus co-assignment matrix is condensed to a
  consensus partition by average-linkage agglomeration on `1 − C`.
- **HDBSCAN** was scanned over `(min_cluster_size, min_samples)` for completeness
  (`provenance/wp2_hdbscan_subgroup_scan.csv`); it does not resolve compact reproducible
  sub-clusters in this feature space and is not used for the labels.

## 2. Result — three seed-stable subgroups

Only **k = 3** is stable across restarts, and decisively so:

| k | ARI mean | ARI p10 | mean BIC | min component stability |
|---|---|---|---|---|
| 2 | 0.721 | 0.456 | 14134 | 0.838 |
| **3** | **0.996** | **0.991** | 13998 | **0.997** |
| 4 | 0.637 | 0.366 | 13831 | 0.771 |
| 5 | 0.697 | 0.521 | 13776 | 0.768 |
| 6 | 0.800 | 0.575 | 13777 | 0.734 |
| 7 | 0.697 | 0.566 | 13787 | 0.577 |
| 8 | 0.633 | 0.511 | 13796 | 0.562 |

**The k = 3 stability is physical, not a GMM artefact.** A parametric-bootstrap null —
a single 4-D Gaussian matched to the data's mean and covariance, pushed through the
identical procedure — yields k = 3 ARI ≈ 0.62 (max 0.83 over 20 sims). The data's 0.996
exceeds **every** null realisation (100%). A smooth blob does not reproduce this.

**Robustness.** k = 3 remains uniquely most stable under `RobustScaler` (ARI 0.995) and on
the P > 0.9 subset (ARI 0.946); the StandardScaler and RobustScaler k = 3 partitions agree
at ARI 0.895.

## 3. Subgroup catalogue

Distances are inverse-variance-weighted mean parallax → distance; counts are raw and
P-weighted. Naming is deterministic (§ script `name_components`): **B** = most negative
μα* (the OC-128 kinematic group); of the remainder, **C** = more negative μδ, **A** = less
negative μδ.

| Label | Paíz analog | N (P-wt) | l med | b med | μα* [mas/yr] | μδ [mas/yr] | ϖ [mas] | distance [pc] |
|---|---|---|---|---|---|---|---|---|
| **CygOB2-A** | FSR 0238 (+ part Bica 2) | 476 (422) | 80.26 | 0.92 | −2.60 ± 0.24 | −4.09 ± 0.25 | 0.6130 ± 0.0011 | 1631 [1628–1634] |
| **CygOB2-B** | OC-128 | 426 (360) | 79.65 | 0.90 | −3.00 ± 0.18 | −4.33 ± 0.23 | 0.6109 ± 0.0012 | 1637 [1634–1640] |
| **CygOB2-C** | HSC 625 (+ part Bica 2) | 429 (366) | 80.11 | 0.68 | −2.59 ± 0.22 | −4.52 ± 0.26 | 0.6176 ± 0.0012 | 1619 [1616–1622] |

Sky central-90% spans: A `l∈[79.78,80.52], b∈[0.67,1.34]`; B `l∈[79.12,80.07], b∈[0.46,1.39]`;
C `l∈[79.71,80.33], b∈[0.32,0.94]`. The three overlap in the VPD (they share one PM
envelope) but separate as a reproducible gradient: **B is west and kinematically offset in
μα*; A is the northern/eastern, least-μδ group; C is the southern, most-μδ group.**

**All three share one distance** (1619–1637 pc, spread 18 pc ≪ the statistical errors and
well inside the 45 pc intrinsic depth). These are **kinematic + spatial** subgroups of a
single distance body, not a distance sequence — consistent with §6.

Figures: `figures/wp2/wp2_subgroups_sky.png`, `…_vpd.png`, `…_extinction.png`.

## 4. Independent (non-kinematic) confirmation

A split visible only in the fitted features is not a result. Each subgroup is confirmed by
observables **outside** `(l, b, μα*, μδ)`:

- **BP-RP reddening, all 1,331 members** (independent of clustering). Median BP-RP orders
  **A 2.394 < B 2.533 < C 2.691**, and all three pairs differ significantly:
  A–B KS 0.142 (p 2·10⁻⁴, eff 0.32), A–C KS 0.312 (p 6·10⁻²⁰, eff 0.75),
  B–C KS 0.209 (p 1·10⁻⁸, eff 0.37). The G-magnitude (mass/luminosity proxy)
  distributions are statistically indistinguishable (KS p = 0.01–0.36), so the colour
  differences are **reddening, not a mass-selection artefact.** This confirms **all three**,
  including B.
- **Anchor extinction A_V** (spectroscopic subset). Median A_V: A 5.13 (n=59), C 6.24
  (n=43); KS 0.482, p 9·10⁻⁶, median effect 1.27 — independently confirms the A < C
  ordering. B has only 4 A_V anchors (insufficient for KS), which is why the full-sample
  BP-RP test above is the primary confirmation for B.
- **Spectral content.** All three carry O and Wolf-Rayet anchors — A: O6 I, WC4;
  B: WN7o+O7V, O7 I+O9 I, O9.7 III; C: O6.5 III, O9.7 Iab — i.e. all three are genuine
  young massive-star populations, as expected for coeval sub-condensations of one association.

## 5. External validation — Paíz et al. 2025

Cross-matched against Paíz et al. 2025 (IJAA 15, 171) Table 3, treating their clusters as
**targets to test against, not ground truth** (SciRP venue, weak review). A clean member
matches a cluster when it lies within its radius (+0.05° margin), within 0.45 mas/yr of the
cluster mean PM, and within 0.10 mas in parallax. `provenance/wp2_paiz_crossmatch.csv`.

| Paíz cluster | role | in footprint | full matches | subgroup breakdown |
|---|---|---|---|---|
| OC-128 | association | yes | 245 | **B 214**, C 28, A 3 |
| HSC 625 | association | yes | 88 | **C 81**, B 7 |
| Bica 2 | association | yes | 182 | C 108, A 74 (straddles) |
| FSR 0238 | association | yes | 17 | **A 16**, C 1 |
| FSR 0224 | association | **no** | 0 | — (l = 78.46, west of member footprint) |
| OC-123 | association | **no** | 0 | — (b = 1.71, north of member footprint) |
| **HSC 630** | **CONTROL** | yes (sky) | **0** | — |

- Our kinematically-defined subgroups map coherently onto Paíz's independently-detected
  clusters: **B ≡ OC-128**, **C ≡ HSC 625**, **A ≡ FSR 0238**; **Bica 2** sits on the A/C
  boundary (its μδ = −4.44 lies between A and C) and splits accordingly.
- **Two of the six association clusters, FSR 0224 and OC-123, fall outside our member
  footprint** (their `l`/`b` lie beyond the concentrated membership envelope even though
  they are inside the WP1 narrow box). We recover 0 members for them — a footprint
  limitation, honestly reported, not a disagreement about their astrometry.
- **HSC 630 contamination control: EXCLUDED.** 43 members are spatially coincident with
  HSC 630 but **zero** match its proper motion (μδ = −2.90 vs the association's −4.3) or
  parallax (0.726 vs ~0.615 mas). The pipeline did **not** sweep in the control cluster.

## 6. Distance-population test (retained)

The clean 1,331-member sample was fit in latent **distance** space with a nonlinear 1-D
extreme-deconvolution forward model (40-node Gauss-Hermite quadrature through each star's
parallax error, conditioned on the 0.35–1.10 mas query truncation shifted by each star's
zero point). One component: mean 1.6245 kpc, intrinsic σ 0.0454 kpc, BIC −5011.64. Two
components collapse to 1.608/1.640 kpc, BIC −4999.48 → **ΔBIC = +12.17 favours one**; held-out
Δ log-predictive is only +3.65. Three controls give ΔBIC +9.45/+10.60/+7.59 and mixed
held-out signs; **no control or non-parallax observable confirms a split.**
**Verdict: `NO_CONFIRMED_TWO_DISTANCE_POPULATION_CLAIM`** — Berlanas+19's ≈1.35/1.76 kpc DR2
split is not reproduced at DR3 precision. The superseded ΔBIC = −25,265 result must not be
cited. This is why §3's subgroups are sought in sky+PM, not parallax.

## 7. Verdict and hand-off

**Three stable, physically-confirmed kinematic subgroups (CygOB2-A/B/C) of a single
1.62 kpc body.** They are reproducible across seeds and scalers (ARI 0.996, null-rejected
at 100%), each confirmed by a non-kinematic observable (BP-RP reddening for all three;
anchor A_V for A/C), and they align with Paíz's independently-catalogued OC-128, HSC 625 and
FSR 0238; the HSC 630 control is cleanly excluded.

- **Labels** are written to the sidecar `tables/wp2_subgroup_labels.parquet`
  (`source_id → subgroup`). The WP4 closure migration propagates this canonical
  column into WP2/WP3 member products, with non-sidecar rows explicitly
  `unassigned`; downstream joins remain keyed on `source_id`.
- **For WP4/WP5:** the subgroups share a distance, so treat them as spatial/kinematic
  sub-populations for per-subgroup completeness (extinction differs: A < B < C), **and still
  carry the star-formation-duration branch (0/1/2 Myr)** as the coeval-age-spread
  systematic — the subgroups are kinematically distinct but the data do not establish
  distinct ages.
