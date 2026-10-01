# Issue #20 — is Cyg OB2-C really young? Resolve it, and study C if it is

*Opened 2026-10-01. Register entry: [PROJECT_TRACE.md §9, issue 20](../PROJECT_TRACE.md).
Scope diagnostic (read-only, reproducible):
[issue20_c_age_diagnostic.py](../scripts/issue20_c_age_diagnostic.py) →
[issue20_c_age_by_version.csv](../tables/issue20_c_age_by_version.csv) ·
[issue20_bright_members.csv](../tables/issue20_bright_members.csv) ·
[issue20_hrd_age_scan.csv](../tables/issue20_hrd_age_scan.csv).
Related: [issue #19](issue19_stale_wp4_inputs_brief.md) (the WP4 anchor-HRD gate was scored at pre-repair ages) and
[WP13](wp13_pooled_vs_resolved_ablation_brief.md) (its premise depends on C's age).*

```sh
PYTHONPATH=scripts python scripts/issue20_c_age_diagnostic.py   # conda env cygob2-gaia
```

---

## 1. Why this matters

The project's surviving novelty claim is that resolving Cyg OB2 into subgroups changes the death history *because the young C subgroup sits below the first-death boundary*. On the baseline, C is at 2.52 Myr and contributes **0** deaths, against 4.17 for A and 4.26 for B.

**If C is not young, both the zero and the claim go.** If C *is* young, or has a young component, then C is a scientific target in its own right.

## 2. What the scope diagnostic found

### 2.1 C's young age appeared with the extinction repair

Upper-MS photometric age (MAP), baseline settings (R_V = 3.1, f_bin = 0.4, no distance offset):

| | pre-repair | repair_v1 | repair_v3 | repair_v5 |
|---|---:|---:|---:|---:|
| C, PARSEC | 3.98 | **2.51** | 2.51 | 2.51 |
| C, MIST | 3.57 | 3.18 | 3.18 | 3.18 |
| A, PARSEC | 4.47 | 3.98 | 3.98 | 3.98 |
| B, PARSEC | 3.98 | 2.82 | 2.82 | 3.55 |

### 2.2 C's brightest members are mostly evolved supergiants, with very young O stars mixed in

| M_G,0 | spectral type | T_eff (K) | reading |
|---:|---|---:|---|
| −7.51 | B0 Ia | 30,900 | evolved, ≳3.5 Myr if single |
| −7.29 | O5 I + O3.5 III | 38,900 | binary; the O3.5 III companion is young |
| −7.13 | B1 Ib | 24,500 | evolved |
| −6.98 | B0 Ia | 30,900 | evolved |
| −5.92 | O9.7 Iab | 27,700 | evolved |
| −5.86 | **O3 If** | 42,700 | among the youngest massive stars |
| −5.82 | O6.5 III | 37,200 | giant |

A colour–magnitude fit cannot tell a hot supergiant from a hot main-sequence star, because optical colours saturate above ~20,000 K.
- PARSEC can put main-sequence stars at M_G ≈ −8.3 at 2.5 Myr, so it reads C's brightest stars as young main-sequence stars, which gives 2.51 Myr.
- MIST's very massive main sequence never gets brighter than about −6.8, so it reads them as evolved stars, which gives 3.18 Myr.

### 2.3 The spectroscopic stars prefer ~3.5–4 Myr for C, the same as A

Each anchor's spectroscopic T_eff and M_G,0 is matched to isochrones with the project's own nearest-point metric (`wp4_anchors_hrd.match_hrd`). The table gives total χ²; lower is better.

| | PARSEC 2.51 | PARSEC 3.55 | PARSEC 3.98 | MIST 2.52 | MIST 3.57 | MIST 4.01 |
|---|---:|---:|---:|---:|---:|---:|
| C, 43 stars | 176.8 | **144.0** | 151.5 | 292.7 | 141.8 | **134.8** |
| C, 9 brightest | 29.5 | **16.2** | 30.1 | 171.8 | 33.6 | **30.3** |
| A, 59 stars (control) | 211.1 | **195.1** | 206.7 | 231.2 | **176.1** | 181.8 |

**These numbers are diagnostic only:**
- the matching is nearest-point, with no IMF weighting and no binary model;
- the reduced χ² is ~3–4, so the Δχ² values overstate significance;
- the anchor sample is biased toward luminous stars;
- converting spectral type to T_eff is uncertain for supergiants.

The original WP4 HRD consistency gate passed **at the pre-repair ages** (C at 3.98 Myr). It was never re-scored at 2.51 Myr (issue #19, F2).

## 3. Three hypotheses (pre-register these before any test)

| | hypothesis | what it predicts |
|---|---|---|
| **H1 — C is young** (~2.5 Myr) | The supergiants are interlopers (A/B members seen in projection, or field), or binary/merger products. | The supergiants' proper motions and positions fit A/B better than C. Per-star ages of the remaining C stars are young. Removing the supergiants keeps a young photometric fit. |
| **H2 — C is coeval with A/B** (~3.5–4 Myr) | The young age is an artefact: evolved supergiants read as main-sequence stars, possibly helped by the repair_v1 extinction change. | The supergiants are kinematic C members. The O3 stars are the anomaly: rejuvenated merger products or "blue stragglers". A hybrid CMD+HRD fit gives ≥3.3 Myr in both families. |
| **H3 — C has two components** | A younger episode (O3 stars) superposed on an older one (B supergiants): an extended or multi-episode history, as Berlanas et al. (2019, 2020) report for the region. | Per-star spectroscopic ages are bimodal. A two-age mixture beats a single age (ΔBIC > 10). The two components may differ spatially or kinematically. |

## 4. Work plan

### Phase A — decide which hypothesis holds (read-only, ~2–3 days, do this first)

- [ ] **A1. Pre-register** (`scripts/issue20_prereg.py` → `provenance/issue20_prereg.json`, with input SHA-256s). Fix the hypotheses, every test statistic and the decision rule (§5) before running A3–A6. Failed predictions stay recorded as failed.
- [ ] **A2. Star-by-star audit** of C's 43 spectroscopic stars and its 15 brightest members (`tables/issue20_c_star_audit.csv`). For each star record:
  - spectral type and luminosity class;
  - binarity (RUWE, plus SB/eclipsing flags from the literature);
  - membership probability;
  - proper-motion Mahalanobis distance to the A, B and C centroids;
  - position relative to each centroid;
  - spectroscopic versus photometric A_V;
  - literature identifier (Cyg OB2 #, MT91, Schulte numbers).
- [ ] **A3. A proper spectroscopic-HRD age likelihood** per subgroup, both families and all three R_V values:
  - per-star likelihood over the isochrone, IMF-weighted;
  - an unresolved-binary component, as in WP4;
  - T_eff calibration error of σ_logTe 0.03–0.05, using a supergiant T_eff scale where the luminosity class is I/II;
  - output: the posterior per subgroup, plus a per-star single-isochrone age where defined.

  Then test single-age against two-age mixtures for C (H3).
- [ ] **A4. Photometric sensitivity refits of C** (sensitivity only, never adopted from this step):
  - (a) a hybrid likelihood in which spectroscopic stars use their spectroscopic T_eff;
  - (b) the upper-MS fit with luminosity class I–II stars excluded.
- [ ] **A5. Attribute the 3.98 → 2.51 jump.** Compare A_V and M_G,0 of C's brightest stars between the pre-repair and repair_v1 extinction. Did the repair brighten them?
- [ ] **A6. Closure cross-check.** Compute C's closure ratio and closing slope at C ages of 2.51, 3.16, 3.55 and 3.98 Myr.
  - Expected direction: an older C lowers the turnoff, so fewer massive stars are predicted alive, so C's existing excess (closing slope 2.06) grows.
  - If so, an older C makes C's massive-star excess *harder* to explain. Record this as evidence either way.
- [ ] **A7. Literature check.** Published ages for C's region and members: Wright+2015, Berlanas+2019/2020, and spectroscopic studies of the O3 If stars. Treat it as a cross-check, never a calibration (rule 3).

### Phase B — act on the verdict

- **H2 (coeval).**
  - Re-derive C's age prior from the hybrid CMD+HRD likelihood.
  - Rerun WP5 → WP6 → WP7 as a separately attributable step after #19's repair_v8 (call it `repair_v9`, or a labelled stage inside v8 with its own predictions).
  - C will then contribute deaths and the baseline total rises.
  - Re-scope WP13: the pooled-versus-resolved difference is expected to shrink, so the "equivalent" outcome becomes likely.
  - Methodological lesson worth stating: photometric upper-MS ages are biased young when evolved hot supergiants are present and spectroscopic T_eff is not used.
- **H1 (young).** The ledger stands, WP13's premise is strengthened, and Phase C opens.
- **H3 (two components).** Treat C as non-coeval:
  - Check whether the existing formation-duration branch (δ up to 2 Myr) captures it, or whether a two-component C is needed in the ledger.
  - The first deaths of the older component are dated by its own age.
  - Phase C opens with the age spread as a central result.

### Phase C — a focus study of C (only if H1 or H3)

Science questions:
1. **Why is C young?**
   - Sequential or triggered formation by A/B feedback? Test the geometry: is C on the side of the A/B cavity edge? Is C's relative motion consistent with gas pushed outward?
   - Or a separate cloud that collapsed later?
2. **Its feedback state.**
   - On the PARSEC relation, C's first death is about τ(120 M☉) ≈ 3.0 Myr after birth, i.e. ~0.5 Myr in the future at 2.52 Myr.
   - C is then a purely wind-driven site, a clean laboratory for wind-only cosmic-ray acceleration in the Cygnus cocoon. Connects to the wind-luminosity program in `method_explained.md`.
3. **The top of its IMF.** C's closing slope of 2.06 is the shallowest of the three. Does a genuinely young C carry a top-heavy IMF, or is it mass segregation or an age effect?
4. **Its dynamical state.** An expansion or traceback age from Gaia proper motions; 3-D kinematics with radial velocities; mass segregation of the O stars.
5. **Star-formation history of Cygnus X.** Do A, B and C form a sequence in time and space?

Data needs:
- spectroscopic follow-up of C's photometric-only bright members;
- Gaia DR4 epoch astrometry and XP spectra (better T_eff for hot stars);
- radial velocities for 3-D kinematics;
- gas maps (CO, dust) for the triggering geometry.

## 5. Decision rule (to be frozen in A1; a starting proposal)

Using A3's posterior, with binaries modelled, in **both** families:
- **H1:** C's 95% upper bound lies below A's 95% lower bound.
- **H2:** C's and A's 68% intervals overlap.
- **H3:** a two-age mixture improves the C likelihood by ΔBIC > 10 over the best single age, *and* the two ages differ by more than 1 Myr.
- **Otherwise inconclusive.** C is then carried as a branch (young / coeval) through the ledger, and no novelty claim rests on C's age.

## 6. Ordering with #19 and WP13

- Phase A is cheap, read-only and decides the paper's direction: **run it first.**
- The #19 fix (repair_v8) can proceed in parallel.
- WP13: freeze its thresholds now, but read M0/M1 only after #20's verdict, because M1's C age depends on it.

## 7. Rules that apply (from CLAUDE.md)

- Nothing is overwritten or retuned.
- A cross-check that disagrees becomes an issue, not a reason to move a number.
- Predictions are written before tests.
- Inputs go through `wp10_inputs.py`. Exception: this diagnostic reads the unversioned `wp4_anchor_hrd.parquet` for spectral types and T_eff only. Phase A must register a versioned spectroscopic table instead.
