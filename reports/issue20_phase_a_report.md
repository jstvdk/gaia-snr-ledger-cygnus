# Issue #20, Phase A: is Cyg OB2-C really young?

*2026-10-02. Brief: [tasks/issue20_subgroup_c_age_brief.md](../tasks/issue20_subgroup_c_age_brief.md).
Pre-registration, committed in `de62dd8` before any test ran:
[provenance/issue20_prereg.json](../provenance/issue20_prereg.json) (with the
method module [scripts/issue20_common.py](../scripts/issue20_common.py), whose hash it records).
Outcome: [provenance/issue20_phase_a_outcome.json](../provenance/issue20_phase_a_outcome.json).
Phase A is read-only. No chain product changed, and repair_v8 is still the quoted chain.*

## 1. Verdict

**By the frozen rule the verdict is H3 ("C has two components"). It carries the
pre-registered flag "mixture preference not specific to C", so Phase B needs a
user decision.**

- **H3 test.** For C, the best two-age mixture beats the best single age by
  ΔBIC = 22.7 (PARSEC) and 15.4 (MIST). The two ages are more than 1 Myr apart,
  and each component holds at least 3 stars. This holds on all 18 robustness
  cells (R_V × f_bin × family) and under the secondary (joint) statistic.
- **The same test on A also passes**, in both families at the decision cell
  (ΔBIC 22.6 PARSEC, 10.006 MIST) and on 14 of 18 cells. The qualifier written
  for this case therefore applies.

What the numbers actually show, in plain terms:

1. **The spectroscopy does not support a young C (H1).** C's single-age
   spectroscopic-HRD age is **3.98 Myr (PARSEC; 68 % 3.54–4.23)** and **4.01 Myr
   (MIST; 3.94–4.49)**. A's is 3.16 (2.97–3.58) and 3.57 (3.39–3.76). C is *at
   least as old as A* in every cell. H1 fails in all 30 conditional runs.
   - Taken alone, the single-age rule gives H2 in PARSEC, because the 68 %
     intervals overlap.
   - In MIST it gives neither hypothesis: C is *older* than A, the 68 % intervals
     are disjoint, and the 95 % intervals overlap.
   - The single-age verdict is therefore INCONCLUSIVE, in the direction of C
     being older, not younger.
2. **The "two components" are not the H3 the brief describes.** The brief's H3
   is young O3 stars plus older B supergiants. What the test found is different:
   - The older component sits at the **oldest age on the grid (10.0 / 10.12 Myr)
     in both A and C, in every cell**.
   - It soaks up a few stars that are too luminous for their T_eff: Schulte 64
     (B2 Ve, which SIMBAD lists as a Be star), Schulte 2 (B1 I at only
     M_G0 = −4.64), and MT91 575 in A.
   - C's "young" component is just C's own age (3.55 / 4.01 Myr).

   This is outlier absorption, a sign of a mis-specified model. It is not a
   second star-formation episode. The frozen H3 criteria had no guard against a
   component sitting on the grid edge. That is a design gap in the rule, and it
   is recorded here rather than repaired after the fact.
3. **Photometry and spectroscopy genuinely disagree, and the cause is not the
   one H2 proposed** (§3). C's photometric young age:
   - survives removing its supergiants;
   - is carried entirely by the extinction of its *photometric-only*
     upper-main-sequence members;
   - is not affected at all by the spectroscopic stars.

## 2. A3: spectroscopic-HRD age likelihood

**Method** (as registered):
- per star, the conditional likelihood p(log T_eff | M_G0, t), using an IMF-weighted particle model on every native isochrone age;
- all evolutionary phases included;
- unresolved binaries with q ~ U(0.1, 1);
- σ_logTe of 0.03 for classes III–V and 0.05 for classes I–II;
- σ_MG0 with a 0.25 mag floor;
- type-derived T_eff of class I/II stars replaced by the supergiant scale (Martins+2005 Table 6; Crowther+2006 Table 4).

The secondary statistic is the joint density within M_G0 ≤ 0.

| decision cell (R_V 3.1, f_bin 0.4) | A: MAP [68 %] {95 %} | C: MAP [68 %] {95 %} | per-family label |
|---|---|---|---|
| PARSEC | 3.16 [2.97, 3.58] {2.72, 4.04} | 3.98 [3.54, 4.23] {3.28, 4.53} | H2 (single-age) · H3 (mixture) |
| MIST | 3.57 [3.39, 3.76] {3.20, 4.08} | 4.01 [3.94, 4.49] {3.77, 4.72} | none (single-age) · H3 (mixture) |

Source: [tables/issue20_hrd_posteriors.csv](../tables/issue20_hrd_posteriors.csv).
C's MAP never falls below 3.55 Myr in any of the 30 conditional runs, including:
- every R_V and f_bin;
- σ_logTe = 0.05 for all stars;
- σ_MG0 floor 0.40;
- distance ±0.06 mag;
- f_bin = 0.7;
- table T_eff without the supergiant scale.

Full grid: [tables/issue20_hrd_mixture.csv](../tables/issue20_hrd_mixture.csv),
robustness cells in the outcome JSON.

**Per-star ages** (report only; [tables/issue20_hrd_per_star.csv](../tables/issue20_hrd_per_star.csv)), PARSEC / MIST MAP in Myr:

| star | type | M_G0 | age (PARSEC / MIST) |
|---|---|---:|---|
| Schulte 7 | O3 If | −5.86 | 1.8 / 2.0 |
| Schulte 9 | O5 I + O3.5 III | −7.29 | 1.6 / 3.2 |
| Schulte 8C | O5 III | −5.76 | 2.2 / 2.8 |
| [CPR2002] A27 | B0 Ia | −6.98 | 4.5 / 5.7 |
| Schulte 18 | B1 Ib | −7.13 | 5.6 / 7.2 |
| A29 | O9.7 Iab | −5.92 | 6.3 / 7.2 |

The hottest O stars look young and the B supergiants look old. That spread
matches Berlanas+2020's bursts at ~3 and ~5–6 Myr, with a possible ~1.5 Myr
group made of the hottest stars (§6). It cannot tell H3 (several episodes) from
H2 (young-looking O stars that are merger products or blue stragglers in an
older population). Berlanas+2020 leave the same question open.

## 3. A4 and A5: photometric sensitivity and the 3.98 → 2.51 Myr jump

Source: [tables/issue20_photometric_refits.csv](../tables/issue20_photometric_refits.csv).
All values are at R_V 3.1.

| fit of C | PARSEC MAP | MIST MAP |
|---|---:|---:|
| repair_v5 upper-MS, unchanged (I1 reproduces the stored value exactly) | 2.51 | 3.18 |
| A4(b): the 7 class I/II anchors removed | 2.51 | 3.18 |
| A4(a): hybrid, photometric CMD + spectroscopic HRD | 2.82 | 3.57 |
| A4(a) control: A, hybrid | 4.47 | 4.01 |

A5 attribution ([tables/issue20_extinction_attribution.csv](../tables/issue20_extinction_attribution.csv)):
- **The repair did not brighten C's brightest stars.** All 43 spectroscopic stars have identical A_V and M_G0 before and after repair_v1, because their A_V comes from intrinsic colours.
- The one bright non-anchor (A12) became *fainter*: A_V 8.39 → 7.02, M_G0 −7.43 → −6.20.
- Counterfactual fits of C (PARSEC):
  - repair_v1 with the anchors' rows taken from the pre-repair run: still **2.51**;
  - repair_v1 with the non-anchors' rows taken from the pre-repair run: **3.98**.
- The whole jump is therefore carried by the extinction of C's photometric-only members. The same holds for MIST: 3.18 vs 3.57.

So the mechanism H2 proposed (evolved supergiants read as main-sequence stars)
is **not** what makes C photometrically young. The photometric young age comes
from C's 369 photometric-only upper-main-sequence stars (R_V 3.1), through
their repair_v1 extinction. How many of them moved, and by how much, was not
tabulated in Phase A. The 43 spectroscopic stars say about 4 Myr. In the hybrid fit the
photometric stars outweigh the spectroscopic ones in PARSEC (2.82), but less so
in MIST (3.57).

This tension sits between the repaired photometric extinction of C's fainter
upper-MS stars and C's spectroscopy. It is the same class of question that
issue #1d settled for B with a near-IR cross-check.

## 4. A6: closure at fixed C ages

[tables/issue20_closure_by_c_age.csv](../tables/issue20_closure_by_c_age.csv).
Only the integral's upper limit varies. The responses, k and the census are held
at their stored values, so this does not capture how k depends on age.

I4 reproduces the stored repair_v8 rows exactly.

| C age (Myr) | 2.51 | 3.16 | 3.55 | 3.98 |
|---|---:|---:|---:|---:|
| closure ratio, PARSEC R_V 3.1, α 2.3 | 1.360 | 1.377 | 1.401 | 1.417 |
| median closing slope over 6 cells | 2.048 | 2.038 | 2.034 | 2.021 |

The direction is as expected: an older C makes C's excess slightly *larger*.
The effect is small, because the predicted living stars above 8 M☉ are dominated
by 8–30 M☉.

In MIST, raising the cap above about 70 M☉ changes nothing. C's injected
responses do not extend above its node turnoffs.

## 5. A2: star-by-star audit

[tables/issue20_c_star_audit.csv](../tables/issue20_c_star_audit.csv): 43 spectroscopic stars plus A12, the one bright non-anchor.

- **Kinematics.** Of C's 7 class I/II stars, 5 are nearest to C in proper motion.
  - The two exceptions are both B supergiants. Schulte 18 (B1 Ib) has Mahalanobis distance 0.73 to A against 1.51 to C. [CPR2002] A27 (B0 Ia) has 0.36 to A against 1.52 to C.
  - Two of C's four B supergiants therefore move like A members. This is partial support for H1's "interloper" reading, for those two stars.
- **Binarity.**
  - SIMBAD flags 11 of the 44 stars as SB / EB / double, among them Schulte 9 (SB) and Schulte 8B/8C (double).
  - No star has RUWE > 1.4.
  - Schulte 64 is a Be star, the main contributor to the 10 Myr mixture component.

## 6. A7: literature (cross-check only)

- **Wright+2015** ([arXiv:1502.05718](https://arxiv.org/abs/1502.05718)): star formation "more or less continuously between 1 and 7 Myr ago".
- **Berlanas+2020** ([arXiv:2008.09917](https://arxiv.org/abs/2008.09917)):
  - O stars show "at least two star-forming bursts at ∼3 and ∼5 Myr", with the older group on the ~6 Myr (Geneva) isochrone;
  - A29 (in our C) is one of the older group's supergiants;
  - a ~1.5 Myr group "containing the hottest stars" (Cyg OB2 #7 = Schulte 7, #22A) is possible but "depends on the possible evolutionary paths";
  - merger products would instead place them at 3–4 Myr;
  - the age groups cannot be separated spatially.
- **Negueruela+2008** (as summarised by Berlanas+2020): 2.5 Myr for the association, with an older population; they judged a separate younger population for the most massive stars "not the most probable case".

The literature is consistent with our per-star spread. It does not single out C
as a young subgroup.

## 7. Secondary predictions

| | prediction | measured | outcome |
|---|---|---|---|
| S1 | C's class I/II stars are kinematic C members (H2/H3) or not (H1) | 5/7 nearest C (0.71) | consistent with H2/H3; weak, because the labels are PM-based |
| S2 | H1: removing class I/II keeps C young | 2.51 / 3.18 vs A − 0.5 = 3.48 / 3.51 | **young persists** (consistent with H1) |
| S3 | H2: the hybrid fit gives ≥ 3.3 Myr in both families | 2.82 / 3.57 | **FAIL** |
| S4 | an older C increases C's closure excess | monotone on 6/6 cells; slope 2.048 → 2.021 | PASS |

S2 and S3 point the other way from A3. The photometric data keep C young; the
spectroscopic data do not. Read together:
- the disagreement is between two data sets;
- it is not about how the supergiants are treated.

## 8. Integrity, deviations and expectations

**Integrity.** I1–I5 PASS, all with exact reproduction:
- I1: the repair_v5 C MAPs;
- I2: the pre-repair 3.98 and repair_v1 2.51;
- I3: the diagnostic's χ² scan, from the versioned anchor table;
- I4: the stored closure;
- I5: the IMF bookkeeping.

For I6, every input was read through the hash-checked `frozen()`.

**Deviations, recorded:**
1. The method module was not changed after the pre-registration. Before it was frozen, two pre-test edits were made:
   - the MIST gap rule changed from a mass jump to an observable jump, because the mass jump wrongly cut the contiguous join from the pre-main sequence to the main sequence;
   - star-chunking was added for memory, with no numerical effect.
2. The analysis drivers (`issue20_hrd_likelihood.py`, `issue20_photometric.py`, `issue20_closure_c_age.py`, `issue20_star_audit.py`, `issue20_score.py`) were written *after* the pre-registration commit and are not hashed in it. They implement its text. The scorer applies the decision rule verbatim.
3. Robustness was evaluated over the full R_V × f_bin cross (8 cells per family), a superset of the registered "R_V 3.0/3.5; f_bin 0.3/0.5".
4. SIMBAD was queried live on 2026-10-02. The response is cached at `data/raw/issue20_simbad_otypes.csv` and hashed in the A2 execution record.
5. **Design gap in the frozen rule.** The H3 criteria do not reject a mixture component at the grid edge. The verdict is still reported as H3 and is not re-scored. The registered specificity qualifier is what flags the problem.

**Expectations recorded before the test.**
- The previous agent expected H2 or H3 to be plausible. The present agent judged H2 or INCONCLUSIVE more likely than H1, with H3 possible.
- The single-age outcome (INCONCLUSIVE, with C older) and the formal H3 both fall within those expectations.
- H2's own predictions failed: S3 failed, and S2 came out the way H1 predicts.

## 9. What this means for N_death and for WP13

**Registered consequence of H3:** the expected direction for N_death is **up**,
by less than under H2.

Because the H3 preference is not specific to C, Phase B is a user decision.
Options:

| option | what it means | effect on the chain |
|---|---|---|
| **(a) accept H3 as registered** | two-component C in a separately pre-registered repair_v9 after WP13 | N_death up; the size is set by the old component's weight and age |
| **(b) treat the non-specific H3 as mis-specification and use the single-age verdict, INCONCLUSIVE** (the case the qualifier was written for) | carry C as a young / coeval branch; no novelty claim rests on C's age | baseline 8.36 unchanged; a coeval-C branch widens the headline range upward |
| **(c) pre-register a Phase A′** | an HRD likelihood with an explicit outlier component and a guard against grid-edge components, plus a near-IR extinction cross-check of C's photometric-only upper-MS stars (as issue #1d did for B) | decides between "photometry right" and "spectroscopy right"; read-only, about 1–2 days |

Whichever option is chosen, **"C is young" is not established**:
- C's own spectroscopy places it at or beyond A's age;
- the surviving novelty claim's premise (C below the first-death boundary) rests on the photometric age alone, which now has an identified, testable dependence on the repair_v1 extinction of C's fainter stars.

WP13 still runs on repair_v8 as the stage brief requires. Its framing must not
assume a young C.
