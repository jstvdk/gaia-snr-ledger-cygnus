# S4 — the WP3 adoption obligations, discharged

*Item **S4** of [pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
The `repair_v5` extinction adoption
([wp3_kriging_adoption.json](../provenance/wp3_kriging_adoption.json)) was
accepted on 2026-07-28 subject to four obligations. **O4 was discharged the same
day** ([WP5_TREND_STATISTIC_R3_RESOLUTION.md](WP5_TREND_STATISTIC_R3_RESOLUTION.md)).
This document discharges **O1, O2 and O3** by fixing exactly what the manuscript
must say, and where. It creates no new number.*

---

## O1 — report the revised star-formation history as a result

**The obligation.** *"The age spread is unchanged (1.469 Myr) but CygOB2-B moves
from near C to near A, so the association reads as two older subgroups plus one
younger rather than one older plus two younger. WP7's supernova timeline depends
on this."*

**Discharged as follows.** The star-formation history is a **result** of this
paper, stated in the results section and not buried in the methods, in the form:

> Cyg OB2 resolves into three kinematic subgroups whose counts-based ages are
> **A 4.00, B 4.09 and C 2.52 Myr**: **two older subgroups and one younger**,
> with a total spread of 1.47 Myr. The structure is not coeval, and the
> supernova timeline follows from it — A and B supply the entire supernova
> budget, while C's turnoff has not yet crossed the IMF ceiling and it has lost
> no stars at all.

Three points must accompany it, and all three are consequences of item **B2**
([wp4_wp5_age_reconciliation.md](wp4_wp5_age_reconciliation.md)):

1. **The quoted ages are the WP5 joint age–`k` truth-age posterior means** — the
   counts-based ages, which are what the ledger consumes. The upper-MS CMD MAPs
   (3.98 / 3.55 / 2.51 Myr) appear once, in the methods, labelled as the prior.
   Only B differs materially, by +0.54 Myr, and that difference is worth a
   factor **1.92** in B's supernova count.
2. **B's counts-based age rails against the top of its prior support** — 53% of
   the posterior weight on the topmost of nine nodes, against 0.17 and 0.12 for
   A and C. B's age is therefore a **one-sided lower bound**, and so is its
   contribution to N_SN (+15% if the true age were the top node). This must be
   stated where the ages are quoted, not only in an appendix.
3. **The ordering was not chosen; it changed under a pre-registered fix.** The
   pre-repair chain read as one older subgroup plus two younger. The kriged
   anchor prior moved B by +0.73 Myr with A and C moving by 0.000 Myr, all four
   predictions confirmed in advance. That provenance is what makes the revised
   history a result rather than a re-labelling, and the paper says so.

**The honest envelope stays.** The 2.25–5.67 Myr span across both retained WP4
indicators remains the age uncertainty the paper reports; the three central
values above are not a claim of 0.01 Myr precision.

## O2 — carry the anchor absolute-scale systematic

**The obligation.** *"Broadband photometry sits about 0.5 mag below
spectroscopically calibrated anchors at matched sky position, so the absolute
mass scale and hence N_SN inherit the anchor calibration. Predates repair_v5 and
is not introduced by it."*

**Discharged as follows.** The systematics section states:

> Extinctions for the ~150 spectroscopic anchors are set by fixed-T_eff
> intrinsic colours, and all other members inherit an anchor-informed spatial
> prior. At matched sky position, prior-free broadband photometry prefers an
> A_V about **0.5 mag lower** than the anchors (Wilcoxon p = 4.3 × 10⁻¹⁶). The
> absolute extinction scale — and therefore the absolute mass scale, the IMF
> normalization and N_SN — rests on the anchor calibration rather than on the
> photometry alone. **This is a systematic on the scale, not on the relative
> structure**: it moves all three subgroups in the same direction and does not
> produce the A/B/C ordering.

Two things must be said with it:

- **The sign was checked and is the reassuring one.** Anchors read *higher*, not
  lower, than co-located members' photometry-only A_V — the opposite of the
  failure mode in which a selection-biased anchor set drags the scale down.
  Anchor selection bias was tested for and is absent.
- **It predates the repair.** It is a property of the WP3 method as designed,
  not a cost of the kriged prior, and it is why the paper reports N_SN with a
  branch spread rather than a formal error bar.

## O3 — carry CygOB2-B's calibration asymmetry

**The obligation.** *"B still has only 4 spectroscopic anchors against A's 59
and C's 42. Gaia XP spectra are the natural way to remove it."*

**Discharged as follows.** Stated where B's age and mass function are discussed,
and again in the DR4 outlook:

> CygOB2-B contains **4 spectroscopic anchors inside the member sample, against
> 59 for A and 42 for C**, and its eighth-nearest anchor lies **0.374°** away
> against 0.089° and 0.139°. B's extinction therefore rests on weaker local
> calibration than its siblings'. This is an observational limitation of the
> DR3-era anchor set, not a modelling choice: the kriging weights sum to 1.000
> for A and 0.992 for C but only **0.772** for B, so the correction the kriged
> prior applies is intrinsically B-specific and no per-subgroup decision was
> made.

The paper connects this explicitly to three things it already reports, because
they are the same fact seen from different sides:

| appears as | where |
|---|---|
| B's extinction needed the largest correction (−0.359 mag) | methods, WP3 |
| B's mass-function residual was the last gate blocker (issue #1c/#1d) | methods, WP5 |
| B's counts-based age rails against its prior support (item B2) | results, ages |

**Remedy named, not promised.** Gaia XP low-resolution spectra give
spectrophotometric temperatures for exactly the stars B lacks anchors for, and
DR4 extends them. This is stated as the concrete DR4 improvement for B — in
contrast to the α hinge, which DR4 will *not* settle.

## Where each obligation lands in the manuscript

| obligation | section | form |
|---|---|---|
| **O1** | Results — "The star-formation history" | a headline result with the counts-based ages, the two-older-plus-one-younger structure, B's one-sided age bound, and the pre-registered provenance of the change |
| **O2** | Discussion — systematics budget | a scale systematic on N_SN, with its sign checked and its independence from the subgroup ordering stated |
| **O3** | Results (B's age), Discussion (systematics), and DR4 outlook | the 4-vs-59-vs-42 anchor asymmetry, its three downstream appearances, and Gaia XP as the named remedy |

All three are load-bearing enough that a reader who skips the appendix still
meets them.
