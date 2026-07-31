# B1 — the two association masses, reconciled

*Written 2026-07-30 to discharge item **B1** of
[pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
Definitional work only: no fit was re-run, no stored number moved, and the
script asserts that it reproduces the stored products before it reports
anything. Evidence:
[wp5_association_mass_reconciliation_execution.json](../provenance/wp5_association_mass_reconciliation_execution.json)
· [wp5_association_mass_reconciliation.csv](../tables/wp5_association_mass_reconciliation.csv)
· [wp5_association_mass_reconciliation.py](../scripts/wp5_association_mass_reconciliation.py).*

---

## 1. The complaint

Two numbers circulate under the same name:

| where | number | claim attached |
|---|---:|---|
| WP5 baseline ([wp5_imf_norm_repair_v6.md](../wp5_imf_norm_repair_v6.md) §4) | **29,122 M☉** | "within a factor two of the 16,500 M☉ literature scale" |
| WP6 cross-check table ([PROJECT_TRACE §10b](../PROJECT_TRACE.md), [harer_2025_supernova_rate.md](../cross_checks/harer_2025_supernova_rate.md) §3) | **1.74 × 10⁴ M☉** | "agree, 5%" with Wright+15's 1.65 × 10⁴ M☉ |

Same phrase, same upstream fit, opposite margins — one comfortable inside a
factor-two gate, the other apparently a 5% bullseye. A referee who notices this
before we do has found a real inconsistency in the bookkeeping, even though
neither arithmetic is wrong.

## 2. What each number is

Both are integrals of the **same** fitted normalization `k` from the same
accepted WP5 products. They differ by exactly two terms:

| symbol | definition | 0.08–0.5 M☉? | companions? |
|---|---|:--:|:--:|
| **M1** | `k · ∫₀.₅¹²⁰ M^(1−α) dM` — single power law, primaries only | no | no |
| **M2** | M1 + the 0.08–0.5 M☉ segment of the Kroupa-like broken power law (slope 1.3 below the 0.5 M☉ break, continuous there) | yes | no |
| **M3** | M2 + unresolved companion mass, `f_bin = 0.40`, `q ~ U(0.1, 1)`, companion counted only where `q·m ≥ 0.08 M☉` | yes | yes |

M1 is what [wp6_external_crosschecks.py](../scripts/wp6_external_crosschecks.py)
computes (`IMF_MASS_LO = 0.5`, the frozen `MASS_GRID` lower edge — an artifact of
the grid, not a statement about the IMF). M2 is stored as
`primary_system_mass_median_Msun`. M3 is stored as
`multiplicity_adjusted_mass_median_Msun` and is what `wp5_report.py` prints as
*the* association mass.

**At the accepted baseline (repair_v7, PARSEC, R_V = 3.1, α = 2.3):**

| step | M☉ | running total |
|---|---:|---:|
| M1 — primaries, 0.5–120 M☉ | 17,481 | **17,481** |
| + primaries, 0.08–0.5 M☉ (Kroupa segment) | +6,711 (+38%) | **24,192** = M2 |
| + unresolved companions | +5,053 (+21% of M2) | **29,246** = M3 |

M3/M1 = **1.673**, and it is 1.44–1.98 across all 18 association branches. The
recomputation reproduces the stored `primary_system_mass_median_Msun` to
**7 × 10⁻¹² M☉** over all 36 cells (two versions × 18 branches), so this is the
decomposition of the stored numbers, not a re-derivation of them.

Version drift is negligible: repair_v6 → repair_v7 moves M3 from 29,122 to
**29,246 M☉** (+0.4%), the `k` shift already recorded for repair_v7.

Across the baseline family:

| α | k | M1 | M2 | M3 |
|---|---:|---:|---:|---:|
| 2.0 | 3,668 | 20,105 | 23,892 | 28,996 |
| **2.3** | **5,280** | **17,481** | **24,192** | **29,246** |
| 2.6 | 7,399 | 17,994 | 29,574 | 35,616 |

## 3. Which one belongs next to Wright+2015?

Wright, Drew & Mohr-Smith (2015) §6.3 obtain 16,500 (+3,800/−2,800) M☉ by Monte
Carlo: draw individual stellar masses from a Maschberger (2013) universal IMF
(α = 2.30, β = 1.40) **over 0.01–150 M☉**, and count the total mass produced
until the number of drawn stars at 20–40 M☉ matches their observed 36 (+1/−4).

Two features of that construction decide the pairing:

1. **It integrates the whole IMF, not 0.5–120 M☉.** So M1 is *not* its
   counterpart. M1's near-coincidence with 16,500 M☉ is a cancellation: cutting
   the integral at 0.5 M☉ removes about as much mass as our normalization sits
   above theirs.
2. **The census counts primary OB stars**, and the Monte Carlo matches drawn
   individual stars to that primary count, so unresolved companions are never
   added on top. That makes it a **primary-star** total.

Both point to **M2** as the like-for-like partner:

| comparison | ratio to Wright+15 | offset in their +3,800 M☉ error | verdict |
|---|---:|---:|---|
| M1 / 16,500 | 1.06 | +0.3σ | **not a valid pairing** — different mass range |
| **M2 / 16,500** | **1.47** | +2.0σ | the like-for-like number; a real ~50% offset, comfortably inside the factor-2 gate |
| M3 / 16,500 | 1.77 | +3.4σ | valid only if Wright's figure is read as including companions; still inside the factor-2 gate |

The like-for-like comparison is therefore an **agreement at the factor level and
a ~50% offset in the mean** — not a 5% match, and not a failure either. Wright's
own §6.3 notes that shifting the high-mass exponent by ±0.1 moves their total
mass by ±5,600/−3,900 M☉, i.e. ±30%, which is the resolution at which any two
IMF-extrapolated cluster masses can be compared at all. Our own α = 2.0 → 2.6
branch spread on M2 (23,892 → 29,574 M☉) is the same size and in the same
direction.

Residual definitional ambiguity, stated rather than hidden: Wright's lower
integration limit is 0.01 M☉ against our grid's 0.08 M☉. Continuing our
slope-1.3 segment down to 0.01 M☉ would add ~1,970 M☉ (M2 → ~26,200, ratio 1.59),
but a real Kroupa/Maschberger IMF flattens again in the brown-dwarf regime, so
that is an overestimate. It is reported for scale in the CSV
(`brown_dwarf_extension_0p01_to_0p08_Msun`) and adopted nowhere.

## 4. Consequences

1. **The "agree, 5%" claim is withdrawn as stated.** It compares M1 against a
   full-range literature mass. The honest statement is: *our normalization
   implies a primary-star total mass 1.47× Wright+15's — an agreement at the
   factor level, at the resolution IMF-extrapolated cluster masses admit.* The
   [§10b](../PROJECT_TRACE.md) row and
   [harer_2025_supernova_rate.md](../cross_checks/harer_2025_supernova_rate.md) §3
   are corrected accordingly.
2. **Nothing downstream moves.** The mass comparison is a sanity check that
   enters no likelihood; N_SN depends on `k` and the turnoff, not on the mass
   integral's lower limit. WP6's *rate* comparison — the load-bearing half of
   the Härer cross-check — is untouched, because the rate samples the IMF only
   at the turnoff and never integrates below 8 M☉.
3. **The strongest form of the Härer cross-check survives, and is arguably
   sharpened.** Its whole point (§5 of that document) is that mass and rate
   depend on `k` in *different* ways, so a broken extrapolation would break them
   inconsistently. That argument needs both to agree at the factor level, which
   they do at every definition; it never needed 5%.
4. **The paper must define the quantity where it quotes it.** WP10 quotes M3 as
   the association's stellar mass, states the range and the multiplicity
   convention in the same sentence, and gives M2 as the number to compare with
   census-based literature values.

## 5. What was actually wrong, and why it survived this long

Neither number, and no line of code. What failed is that three nested integrals
of one normalization were all called "the association mass" in prose, in two
documents written a day apart for different purposes, and the factor-1.67
between them was never written down anywhere. The lesson is the same one issue
#17 taught at the other end of the same integral: **an integration limit is part
of a quantity's definition, and a quantity whose limits are not quoted with it
is not yet a number.**
