# `repair_v7` — costed recommendation, not a decision

*2026-07-28. Written after issue #15 closed. **Not executed** — this changes an
accepted WP5 result, so it is the principal investigator's call.*

> ## OUTCOME 2026-07-29 — the discriminator overturned this recommendation
>
> The test in §4 was run and pre-registered
> ([wp5_fbin_discriminator_prereg.json](../provenance/wp5_fbin_discriminator_prereg.json)
> → [wp5_fbin_discriminator_execution.json](../provenance/wp5_fbin_discriminator_execution.json)).
>
> | measurement | result | threshold |
> |---|---:|---:|
> | **D1** — recovery over the 2–8 M☉ calibration window | **+0.16% ± 0.09%** | 2% |
> | **D2** — up-scatter R(est > 8 \| M) over 4–8 M☉ | **+9.87% ± 1.58%** | 2% |
>
> **D2 is five times the threshold. `repair_v7` is JUSTIFIED.** G1 and G2 both
> passed; the control arm reproduced the accepted `repair_v6` node byte for byte.
>
> **The recommendation below to defer was wrong, and its own proposed test is
> what caught it.** The reasoning in §5 — that the sub-8 effect should be
> *smaller* than the 0.4% measured above 8 M☉ because the f_bin change is
> smaller — failed to account for the mechanism being a **threshold** effect.
> Inside the calibration window most stars are recovered either way, so the
> recovery fraction is insensitive (D1 ≈ 0). At the 8 M☉ boundary a small
> brightness shift converts directly into a large probability of crossing, which
> is exactly what G2 predicted and what the 25× ratio between D2 and D1 shows.
>
> Everything below is retained unedited as the pre-test reasoning.

---

**Superseded recommendation: do not run the full chain yet. Run a ~15-minute
discriminator first (§4), and let its answer decide.**

---

## 1. What would trigger it

Issue #15 measured mass-dependent multiplicity **above 8 M☉** and found it
absorbs **3.7%** of the WP6 closure excess. Under the pre-registered decision
rule, that is the "M3 fails" branch and **no `repair_v7` is triggered**.

The open gap is different, and was exposed by issue #17 rather than by #15:

> The injection truth model applies **f_bin = 0.40 at every mass below 8 M☉**,
> and that is where the WP5 calibration window lives.

This matters for two distinct reasons, which should not be conflated:

1. **The normalization `k`.** WP5 fits `k` from 2–8 M☉ counts through the
   injection response. If real 4–8 M☉ stars are more multiple than 0.40, the
   response is wrong in the window that sets `k` — and `k` scales *everything*
   downstream, including N_SN.
2. **The up-scatter channel.** Issue #17 showed that stars with true mass 4–8 M☉
   contribute roughly **24%** of the predicted observable count above 8 M☉
   (R = 0.23 at 7 M☉, 0.09 at 6, 0.03 at 5). Multiplicity acts *hardest* exactly
   here, because this is the population that can be pushed *across* the
   threshold by an unresolved companion. A star already above 8 M☉ cannot be.

**This is why the issue #15 null does not settle the question.** The test was
scoped to the range where the mechanism has least room to act.

## 2. What the literature supports below 8 M☉

Weaker than at the top, and that weakness is the honest reason for caution:

| mass range | multiplicity | source |
|---|---|---|
| O type (≳16 M☉) | >90%, f_b ≈ 0.70 | Sana+2012, Duchêne & Kraus 2013 |
| early B (~8–16 M☉) | ~60–70% | Duchêne & Kraus 2013 |
| late B / A (~2–5 M☉) | ~50% | Duchêne & Kraus 2013 |
| solar | ~50% | Duchêne & Kraus 2013 |

So 0.40 is probably a little low at 4–8 M☉ — but by perhaps 0.10–0.20, not the
0.30 gap that applies at the O-star end. **The anchors are softer here**, and a
`repair_v7` built on them would carry a weaker justification than `repair_v6`
did.

## 3. Cost

| stage | estimate |
|---|---|
| WP3 extinction repair | ~1 h |
| WP4 ages + masses | ~1 h |
| WP5 injections, 162 nodes × full `MASS_GRID` | **~3 h** |
| WP5 IMF fit | ~30 min |
| verdict stability (R3a/R3b/R3c) | ~1 h (OOM-killed once; needs the cache fix) |
| WP6 re-run end to end | ~1.5 h |
| **total** | **~8 h**, plus re-acceptance through G3 |

**The real cost is not compute.** It is that `repair_v6` is an *accepted*
artifact behind a passed gate, and `repair_v7` would have to clear G3, verdict
stability, and the V1 byte-identity check afresh. That machinery exists and
works, but it is the reason this is not a casual re-run.

## 4. The discriminator to run first — ~15 minutes

**Do not spend 8 hours to find out whether the effect is 1% or 15%. Measure that
first, on one branch.**

Take the baseline branch (PARSEC, R_V = 3.1, CygOB2-A) and a single truth-age
node. Inject the existing `MASS_GRID` twice, paired and bit-identical exactly as
issue #15 did, with the *only* change being that `truth_binary_fraction` is
extended below 8 M☉:

```
f_bin(M) = 0.40 at 2 M☉  →  0.55 at 8 M☉  →  0.70 at 16 M☉
```

Then measure two things, which answer two different questions:

| quantity | what it decides |
|---|---|
| shift in R(recovered \| M) across **2–8 M☉** | whether `k` moves — the WP5 question |
| shift in R(estimated > 8 \| M) for **4–8 M☉** | whether up-scatter moves — the WP6 question |

**Pre-declare the threshold before running**, in the project's usual style. A
defensible one: *if either response shifts by more than 2% (comparable to the
+2.4% Orellana distance systematic already carried), `repair_v7` is justified;
below that, the effect is recorded as a systematic and the chain stands.*

Cost: one node, two arms, ~2 minutes of compute plus the pre-registration.

## 5. Why I recommend deferring the full run

1. **The measured analogue was small.** Above 8 M☉ a 0.40 → 0.70 change in f_bin
   moved the response by ~0.4%. Below 8 the f_bin change is smaller
   (0.40 → ~0.55), so a first-order expectation is a *smaller* effect, though it
   acts through two channels instead of one.
2. **It does not address the real open question.** CygOB2-C's residual of 1.405
   against A's 0.894 is a *direction* disagreement between subgroups.
   A uniform change to the truth model moves all subgroups the same way and
   cannot produce that pattern — as issue #15's M2 result already demonstrated.
3. **The anchors are weaker below 8 M☉** (§2), so the scientific justification
   is thinner than for the change already tested.
4. **Higher-value work is queued**: the parallax-blind membership test (issue #8)
   and the distance-contamination test for CygOB2-C, neither of which requires
   touching WP5.

## 6. What would change my recommendation

- the §4 discriminator exceeding its pre-declared threshold;
- a decision to publish the IMF slope as a *primary* result rather than as a
  branch discriminator — that raises the bar on the normalization;
- referee pressure on the constant-f_bin assumption, which is a fair thing to
  ask about and is now documented in [CUTS §18](../CUTS_AND_THRESHOLDS.md) and
  master-table row 20.

## 7. If it is run anyway

Non-negotiables, from what this session learned:

- **pre-register before running** — predictions, thresholds, and the decision
  rule, as `repair_v6` and issues #15/#17 all did;
- **keep `repair_v6` bit-preserved** behind default arguments, and re-run V1;
- **fix the verdict-stability cache first** — it was OOM-killed once by holding
  all 18 branches' responses simultaneously;
- **carry both f_bin models as a Class E branch** rather than replacing 0.40,
  unless the discriminator shows the old value is simply wrong. Master-table
  row 20 already carries f_bin as Class E; it has never actually been carried on
  the response side.


---

## 8. Post-discriminator scope (added 2026-07-29)

**Estimated effect on WP6.** The 4–8 M☉ segment contributes **23.5%** of the
predicted observable count above 8 M☉ (from issue #17's convergence scan:
predicted rises by a factor 1.444/1.105 when the floor drops from 8 to 4). A
+9.87% shift in that segment raises the total predicted count by **2.3%**, moving
the grid-median closure ratio **1.099 → ~1.074**. Modest, and in the direction of
better closure.

**Do not take the shortcut D1 appears to offer.** D1 ≈ 0 suggests the WP5
normalization `k` is stable and that only WP6 needs re-running. That inference is
**not safe**:

- D1 measures the **recovery fraction**, not **mass migration within the
  window**. A star can be recovered in both arms while being assigned a
  different mass.
- D2's size is direct evidence that mass estimates *are* shifting by enough to
  move stars across a threshold. There is no reason that stops at 8 M☉ — the
  same shift redistributes stars among the 2–8 M☉ bins that the WP5 likelihood
  actually fits.
- The WP5 fit consumes the full response matrix `R(observed bin | true mass)`,
  which D1 deliberately compresses to a scalar.

So the honest scope is the **full chain**, as costed in §3, not a WP6-only
re-run. What D1 does buy is a *prior expectation* that `k` will move very little
— which makes the re-run a confirmation rather than a gamble, and means the
headline N_SN is unlikely to shift much.

**Sequencing.** This does not block WP7. The plan's WP7 is *"pure computation on
frozen inputs"* and consumes `k`, the WP4 ages and the runaway correction — none
of which D1 suggests will move materially. Running WP7 on the accepted chain and
`repair_v7` in parallel is defensible, provided WP7's inputs are re-checked
against `repair_v7` before anything is published.
