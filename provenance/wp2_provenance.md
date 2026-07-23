# WP2 provenance log

## Failed baseline run — 2026-07-22

- Notebook: `notebooks/wp2_membership_and_substructure.ipynb`.
- Frozen input: `data/processed/wp1_gaia_narrow.parquet`, 245,843 rows.
- Quality-filtered analysis sample: 199,041 rows.
- Baseline clustering: five-dimensional `(l,b,parallax,pmra,pmdec)` RobustScaler space; DBSCAN `eps=0.42`, `min_samples=15`.
- Baseline Monte Carlo: 100 draws, seed 20260722, diagonal errors only.
- Result: 187,235 initial candidates; 159,450 with P > 0.05; 156,414 with P > 0.5. The result was one percolated structure spanning the selection box, not a physical association catalogue.
- The reported one-versus-two parallax mixture result, `Delta BIC=-25,264.98`, is invalid for physical interpretation because it used the contaminated 160k-star sample, ignored measurement-error deconvolution and the 0.35–1.10 mas truncation, and had no control-field comparison.
- Verdict: **FAILED and superseded**. `data/processed/wp2_members.parquet` and `provenance/wp2_membership_manifest.json` from this run must not be consumed downstream.

## Required-method review — 2026-07-22

- Read `CUTS_AND_THRESHOLDS.md` Sections 6 and 8 before changes.
- `eps` is treated as measured/free: derive a k-distance knee for `k=min_samples`, then scan a range bracketing the knee beginning near 0.03.
- Diagnostics are largest-cluster fraction, spatial extent, Berlanas recall, and control-field yield; cluster count alone is not an acceptance diagnostic.
- The circular three-robust-scale acceptance locus is prohibited. Membership must be a cluster-versus-field posterior odds ratio.
- Publication membership uncertainty target is `sigma_P<=0.005`, requiring 10,000 Monte Carlo draws at P=0.5.
- Monte Carlo must use Gaia's full astrometric covariance, not independent diagonal perturbations.
- Distance-population testing must use the clean member sample, deconvolve per-star errors, include the query truncation in the likelihood, operate in distance space, run identically on controls, and require independent non-parallax confirmation.

## Gate amendment — landed before rerun

- `paper1_execution_plan.md` WP2 gate changed from recall-only to a conjunction of recall, control-field precision, population sanity, spatial compactness, no-percolation, and published-structure comparison.
- Numeric acceptance rules: Berlanas recovery ≥80% at P>0.5; mean area-normalized control yield ≤10% of association yield (target ≤5%); P>0.5 members between 10² and 10⁴; central-90% l and b spans each <80% of their selection-box widths; convex-hull area <50% of selection-box area; largest cluster ≤10% of the quality-filtered analysis sample.
- Counterfactual check against the failed run: member-count and largest-cluster criteria both fail by more than an order of magnitude, and the full-box extent criterion fails. Therefore the amended gate would have blocked the old result.

## Density diagnostics and rejected branches — 2026-07-22

- Global three-dimensional `(parallax, pmRA, pmDec)` DBSCAN scan: `eps=0.02–0.20`, `min_samples={10,15,20}`. The nominal full-sample k-distance knee was `eps=0.1916`, already deep in the field-percolated regime. At `eps=0.03`, `min_samples=10`, the largest component was already 27,646 stars (13.9% of the analysis sample), spanned the complete box, and had a 40.9% control/target ratio. Frozen in `wp2_density_scan_execution.json` and `wp2_dbscan_scan.csv`; rejected.
- HDBSCAN EOM branch: largest components contained 34.9–86.0% of the analysis sample and no component met the fixed association-overdensity rule. Frozen in `wp2_hdbscan_scan_execution.json` and `wp2_hdbscan_scan.csv`; rejected for percolation.
- HDBSCAN leaf branch: largest components were small (0.27–0.47%), but Berlanas seed recall was at most 5.95% and selected components still occupied almost the full box (convex hull about 31 deg2). Frozen in `wp2_hdbscan_leaf_scan_execution.json` and `wp2_hdbscan_leaf_scan.csv`; rejected for fragmentation and spatial contamination.
- A preliminary 1-deg2 Wright-core scan separated the dense-regime knee (`eps=0.0607`) from the field-tail knee (`eps=0.5406`). It established that a target-local seed is necessary, but it was not the published Berlanas footprint used by the gate. Frozen in `wp2_target_density_scan_execution.json` and `wp2_target_dbscan_scan.csv`; diagnostic only.

## Final membership rerun — 2026-07-22

- Adopted the external Berlanas+19 Section 2.1 footprint: radius 1 degree around `(l,b)=(79.8,0.8)` degrees. Three non-overlapping radius-1-degree controls are centred at `(78,3)`, `(82,3)`, and `(82,0.8)` degrees.
- Repeated the DBSCAN scan at `eps={0.03,0.04,0.05,0.06,0.07,0.08,0.10}` and `min_samples={10,15,20}`. The target dense-regime knee is `eps=0.0690`; the global `eps=0.4534` knee is again a field tail. Selected `eps=0.05`, `min_samples=15`, inside the knee bracket, as a training seed only: 1,588 rows, 0.798% of the 199,041-star analysis sample.
- Published membership is not the DBSCAN label. A five-component cluster GMM (selected by BIC in relative sky position plus corrected parallax and proper motions) is compared with a 24-component control-field GMM (also BIC-selected). Membership is posterior cluster-versus-field odds.
- Scanned the cluster prior from 0.01 to 0.08 without consulting Berlanas recovery. The largest prior satisfying the analytic control-leakage ceiling was 0.04. The frozen scan is `wp2_mixture_prior_control_scan.csv`.
- Propagated the full Gaia `(parallax,pmRA,pmDec)` covariance, including all three correlation terms, through 10,000 realizations per plausible star. Membership probability is the fraction of draws with cluster posterior odds greater than field odds. The worst-case binomial standard error is 0.005; no covariance matrices required PSD repair.
- A rejected development attempt subtracted measurement covariance from the empirical membership GMM and then classified repeated noisy realizations against those intrinsic components. It passed precision/spatial checks but produced fewer than 100 `P>0.5` objects and failed recall. This double-use of an intrinsic-population model was rejected before replacing canonical artifacts; see `wp2_membership_attempt_deconvolved_failure.json`.
- Final gate: 1,392 unique `P>0.5` members including manual quality exceptions; Berlanas recovery 189/229 = 82.53% (128 automatic quality-sample recoveries plus 61 explicitly flagged manual quality exceptions); control yields 80, 45, 21 versus 1,331 automatic target members, mean ratio 3.656%; member central-90% spans 1.070 degrees in l and 0.889 degrees in b; convex hull 2.082 deg2; density seed 0.798% of the analysis sample. Every amended criterion passes.
- Published-map check: the automatic members have mean `(l,b)=(79.987,0.842)` degrees and lie within the Wright+15 core / Berlanas+19 footprint. The compact sky structure agrees qualitatively with those maps; the distance split does not, and that disagreement is carried explicitly into the subgroup verdict.
- `wp2_berlanas_recovery_audit.csv` contains all 229 published benchmark rows, including a numerical posterior or explicit per-row quality/soft-floor disposition for every one of the 40 misses. `wp2_anchor_assignments.parquet` is canonicalized to 252 rows and 252 unique Gaia source IDs; the 540-row literature evidence table is not used as a census table.
- The failed 159,450-member baseline remains frozen as `wp2_members_failed_20260722.parquet` with `wp2_membership_manifest_failed_20260722.json`. It is not a downstream input.

## Distance-population test — 2026-07-22

- Fit only the 1,331 clean automatic `P>0.5` members; the 61 manual spectroscopic quality exceptions were deliberately excluded from the distribution fit.
- Implemented a nonlinear one-dimensional extreme-deconvolution forward model in latent distance space. Forty-node Gauss-Hermite quadrature integrates each latent Gaussian through each star's parallax error. The likelihood is conditioned on the raw-query 0.35–1.10 mas truncation after shifting the bounds by each star's zero point.
- One latent component: mean 1.6245 kpc, intrinsic sigma 0.0454 kpc, BIC -5011.64. Two components collapse to nearby means 1.6082/1.6396 kpc with BIC -4999.48. Thus Delta BIC(two minus one) = +12.17 favours one component; three-fold held-out prediction gives only a small +3.65 log-density advantage to the more flexible model.
- Identical control results: Delta BIC = +9.45, +10.60, +7.59 and held-out Delta log predictive density = -2.41, -0.23, +0.64. No control meets both preference criteria.
- The forced two-component responsibilities show no independent confirmation: standardized differences are 0.095 in l, 0.168 in b, 0.001 in pmRA, 0.015 in pmDec, and 0.086 in extinction. The predeclared confirmation rule (effect >=0.30 and KS p<0.01) is not met.
- Verdict: **no confirmed two-distance-population claim in the clean DR3 membership sample**. The old Delta BIC = -25,265 result is superseded and must not be cited.

## Subgroup derivation (Task A) — 2026-07-23

- Script: `scripts/wp2_derive_subgroups.py`. Command: `python scripts/wp2_derive_subgroups.py`.
  Machine log: `provenance/wp2_subgroups_execution.json`; deliverable: `wp2_subgroups.md`.
- Sample: the 1,331 clean automatic members (`membership_probability>0.5 & ~anchor_quality_exempt`).
  Feature space `(l,b,pmra,pmdec)` — **parallax deliberately excluded** (distance depth 45 pc is
  exhausted at DR3; the distance test found one population). StandardScaler.
- Substructure diagnostics reproduced: intrinsic PM dispersion 0.297 mas/yr (≈8× the 0.037 mas/yr
  median error); sky radial KS = 0.61 vs a uniform disc.
- Clusterer: full-covariance GMM, k=2..8, 50 deterministic seeds (1000–1049). Acceptance is
  **seed stability, not BIC** (BIC non-monotonic → local optima). A partition is kept only if
  seed-to-seed ARI mean ≥ 0.90 and its 10th percentile ≥ 0.90. Consensus co-assignment matrix →
  average-linkage agglomeration on 1−C. Frozen: `provenance/wp2_gmm_seed_stability.csv`.
- Only **k=3** qualifies (ARI mean 0.996, p10 0.991); k=2 borderline (0.721/0.456), k≥4 unstable.
- **Parametric-bootstrap null:** a single 4-D Gaussian matched to the data's mean/covariance,
  run through the identical procedure, gives k=3 ARI ≈ 0.62 (max 0.83 over 20 sims). Data 0.996
  exceeds 100% of nulls → the k=3 stability is physical, not GMM regularity on a smooth blob.
- **Robustness:** k=3 uniquely most stable under RobustScaler (0.995) and the P>0.9 subset (0.946);
  StandardScaler vs RobustScaler k=3 partitions agree at ARI 0.895.
- HDBSCAN `(min_cluster_size, min_samples)` scan resolves no compact reproducible sub-clusters;
  reported only, not used for labels. Frozen: `provenance/wp2_hdbscan_subgroup_scan.csv`.
- Deterministic physical naming: B = most negative median μα* (OC-128 group); of the rest, C =
  more negative median μδ, A = less negative μδ. Counts A/B/C = 476/426/429 (P-weighted 422/360/366).
  All three at one distance: A 1631, B 1637, C 1619 pc (spread 18 pc ≪ errors, inside the 45 pc depth).
- **Independent (non-kinematic) confirmation.** (a) BP-RP reddening on all 1,331 members (not a
  clustering feature): medians A 2.394 < B 2.533 < C 2.691; A–B KS 0.142 p 2e-4 eff 0.32,
  A–C KS 0.312 p 6e-20 eff 0.75, B–C KS 0.209 p 1e-8 eff 0.37; G-band (mass proxy) KS p 0.01–0.36
  (indistinguishable) → colour differences are reddening, not mass selection → confirms all three.
  (b) Anchor A_V: A 5.13 (n=59) vs C 6.24 (n=43), KS 0.482 p 9e-6 eff 1.27; B has only 4 A_V anchors.
  (c) All three carry O and WR anchors.
- **External validation — Paíz et al. 2025 (IJAA 15, 171) Table 3** (local paper:
  `papers/Paiz_2025.pdf`; treated as targets, not
  ground truth). Match = within radius+0.05°, within 0.45 mas/yr of cluster mean PM, within 0.10 mas
  parallax. Frozen: `provenance/wp2_paiz_crossmatch.csv`. OC-128 → B (214/245), HSC 625 → C (81/88),
  FSR 0238 → A (16/17), Bica 2 straddles A/C (74/108). FSR 0224 and OC-123 fall **outside** the member
  footprint (l=78.46 / b=1.71) → 0 recovered; reported as a footprint limitation, not a disagreement.
- **HSC 630 contamination control: EXCLUDED.** 43 members are spatially coincident but **0** match its
  PM (μδ −2.90 vs −4.3) or parallax (0.726 vs ~0.615). The pipeline did not sweep in the control.
- **Verdict: three seed-stable, physically-confirmed kinematic subgroups of one 1.62 kpc body.**
  Labels are authoritative in the **sidecar** `tables/wp2_subgroup_labels.parquet`
  (`source_id → subgroup`). The 2026-07-23 WP4 closure migration
  (`scripts/wp4_schema_repair.py`) propagates this one canonical column into the
  WP2/WP3 member products and replaces every missing sidecar match with the
  explicit category `unassigned`. For the P>0.5 gate sample, those missing
  matches are exactly the 61 spectroscopic `anchor_quality_exempt` stars; no
  subgroup is imputed for them. Figures:
  `figures/wp2/wp2_subgroups_{sky,vpd,extinction}.png`.
- Hand-off to WP4/WP5: subgroups share a distance → use them as spatial/kinematic sub-populations
  with per-subgroup extinction (A<B<C) for completeness, **and still carry the SF-duration branch
  (0/1/2 Myr)** as the coeval age-spread systematic (subgroups are kinematically, not age, distinct).

## Execution environment and commands

- Host: macOS 15.6.1 arm64; repository HEAD before these uncommitted changes: `d95f814f295bbc512e283b13133a4c1752d0b4eb`; Conda environment: `cygob2-gaia`; Python 3.11.15.
- Packages used: NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1, scikit-learn 1.9.0, Astropy 8.0.1, Matplotlib 3.11.0, PyArrow 25.0.0, plus `hdbscan` and the Gaia DR3 zero-point package installed in the same environment.
- Diagnostic commands: `python scripts/wp2_membership_pipeline.py --stage scan`, `--stage hdbscan-scan`, `--stage hdbscan-leaf-scan`, and `--stage target-scan`, all via `conda run -n cygob2-gaia --no-capture-output`.
- Canonical commands: `python scripts/wp2_finalize_membership.py`, `python scripts/wp2_distance_population_test.py`, and `python scripts/wp2_finalize_audits.py`, all via the same Conda invocation.
- Validation commands: `python scripts/wp1_validate.py`; Python byte-compilation of all changed WP1/WP2 scripts; explicit Parquet uniqueness, gate, hash, and row-count audits.
