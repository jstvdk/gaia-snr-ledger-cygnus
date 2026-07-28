# Cross-check — Berlanas et al. 2019: the two-distance population

**Source.** S. R. Berlanas, N. J. Wright, A. Herrero, J. E. Drew, D. J. Lennon,
MNRAS **484**, 1838 (2019). Analysed the Gaia DR2 parallax distribution of 229
known OB members of Cyg OB2.

**Their claim.** Line-of-sight substructure within the association: a **main
group at ~1760 pc** and a **foreground group at ~1350 pc** containing **19 stars,
seven of them O-type**. The foreground group is more spatially dispersed than the
main one. (Description as summarised by Orellana et al. 2021 §5, who reproduce
the numbers.)

**Why it matters.** WP2 treats Cyg OB2 as a single population at one distance and
splits it into subgroups **A/B/C in sky position and proper motion, not in
parallax**. If a substantial fraction of our members actually sat 270 pc closer,
their absolute magnitudes — and therefore masses, ages and the WP5 normalization
— would be systematically wrong.

---

## 1. What we tested

`provenance/wp2_distance_population_execution.json`, produced by
`scripts/wp2_distance_population_test.py`.

The 1,331-member clean sample was fitted in **latent distance space** with a
one-dimensional nonlinear extreme-deconvolution forward model: 40-node
Gauss–Hermite quadrature through each star's own parallax error, conditioned on
the 0.35–1.10 mas query truncation shifted by that star's zero point. One
component versus two, compared by BIC, with three control fields and held-out
predictive likelihood.

## 2. Result

| model | outcome |
|---|---|
| one component | mean **1.6245 kpc**, intrinsic depth **σ = 45.4 pc**, BIC −5011.64 |
| two components | collapse to **1.608 / 1.640 kpc**, BIC −4999.48 |
| **ΔBIC** | **+12.17, favours one component** |
| held-out Δ log-predictive | +3.65 (weak, opposite sign to BIC) |
| three control fields | ΔBIC +9.45 / +10.60 / +7.59, mixed held-out signs |
| non-parallax observables | no confirmation — see below |

**The two-component fit does not find 1.35/1.76 kpc.** It finds two components
**32 pc apart**, which is the fit resolving the association's own line-of-sight
depth, not a second population.

### Independent diagnostics on the two fitted components

If the split were physical, the two components should differ in something other
than parallax. They do not, materially:

| observable | standardized mean difference | KS p |
|---|---:|---:|
| l | 0.095 | 0.213 |
| b | 0.168 | 0.0004 |
| μ_α* | 0.0006 | 0.772 |
| μ_δ | ~0.00 | — |

Proper motions are **identical** between the two components. Only b shows a
marginal difference. A real foreground population 270 pc closer would show a
proper-motion offset.

### No foreground tail in the members at all

Direct check on the member parallaxes:

| subgroup | N | median parallax | implied distance | fraction closer than 1430 pc | fraction closer than 1350 pc |
|---|---:|---:|---:|---:|---:|
| CygOB2-A | 476 | 0.6159 | 1624 pc | 2.7% | **0.0%** |
| CygOB2-B | 426 | 0.6171 | 1620 pc | 3.5% | **0.0%** |
| CygOB2-C | 429 | 0.6235 | 1604 pc | 3.3% | **0.0%** |

**Not one member sits at Berlanas's 1350 pc.** The three subgroups agree in
distance to within 20 pc — less than half the fitted 45 pc depth.

**Recorded verdict:** `NO_CONFIRMED_TWO_DISTANCE_POPULATION_CLAIM`.

---

## 3. The caveat that matters — this test is partly circular

**WP2's membership classifier uses `FEATURES = ['parallax_corrected', 'pmra',
'pmdec']`. Parallax is one of the three clustering features.**

A foreground group at 1350 pc would therefore have been assigned low membership
probability and removed **before** the distance test ran. Testing the survivors of
a parallax-based selection for a parallax split will tend to return "one
population" almost regardless of the truth.

So the honest scope of the verdict is:

> **Within the WP2 member sample there is no evidence of two distances.**

It is **not** evidence against Berlanas's foreground group. The register's
one-line summary ("Berlanas two-distance population not confirmed") reads
stronger than the test supports and has been tightened.

### What would settle it

A test on a **parallax-blind** selection — members chosen on sky position and
proper motion only, then examined in parallax. That is a real experiment, not
done here, and it is the only way to break the circularity.

---

## 4. Corroboration from a third source

Orellana et al. 2021 independently find foreground structure in the same field:
two proper-motion groups of 179 and 188 stars at **~1280 pc**, plus the open
cluster UCB585 at ~1460 pc, all with proper motions distinct from Cyg OB2.

This supports a consistent picture across all three analyses:

- **there is foreground structure in the Cyg OB2 field** — Berlanas and Orellana
  both find it, and they broadly agree on where (1280–1350 pc);
- **it is separable by proper motion** — which is how Orellana separate it and
  how our WP2 clustering would;
- **it is not inside our member sample** — zero members closer than 1350 pc.

The disagreement with Berlanas is therefore probably **not** about whether
foreground stars exist, but about whether they should be counted as Cyg OB2
members. Berlanas's 19 stars were drawn from a *known-OB* list; ours from an
astrometric clustering that excludes them by construction.

Note also that Berlanas's main group sits at **1760 pc** against our 1620 pc and
Orellana's 1669 pc — a spread of ~8%, part of the broader distance disagreement
documented in [orellana_2021_astrometry.md](orellana_2021_astrometry.md).

---

## 5. Why this is not currently a threat to WP6 — and where it could become one

**Direction of the bias.** A star truly at 1350 pc but assumed to be at 1620 pc
is assigned a distance modulus **0.40 mag too large**, so it is inferred to be
more luminous and therefore **more massive** than it is. Unrecognised foreground
contamination therefore **inflates the count above 8 M☉** and **raises the
closure ratio**.

That is the direction of CygOB2-C's residual (ratio 1.405, the WP6 outlier). So
the mechanism is a live candidate in principle.

**But our data do not support it here:** C's foreground fraction (3.3% closer
than 1430 pc) is indistinguishable from A's (2.7%) and B's (3.5%), and no
subgroup contains a single star at 1350 pc.

**Not covered by the existing attribution.** WP6's alternative A4
("contamination") tested *membership-probability weighting versus a hard cut* —
a different mechanism entirely. **Distance contamination is not covered by any
WP6 alternative**, and is registered as an open question for the CygOB2-C
residual.

---

## 6. What a failure would have looked like

| failure mode | what we would have seen | actual |
|---|---|---|
| foreground group inside our members | bimodal member parallaxes, components ~250 pc apart | components 32 pc apart |
| real split masked by noise | two components differing in proper motion or position | PM identical, KS p = 0.77 |
| control fields showing the same split | ΔBIC favouring two in controls | all three favour one |

---

## 7. Verdict

**NOT CONFIRMED — with the scope limit stated.** No two-distance structure exists
within the WP2 member sample, by BIC, by control fields, and by direct inspection
of member parallaxes. But because parallax is a membership clustering feature,
this test **cannot** rule out Berlanas's foreground population; it can only show
that it is not in our members.

Two things carried forward:

1. A **parallax-blind membership test** is the experiment that would settle it.
2. **Distance contamination** is an untested candidate for the CygOB2-C closure
   residual, not covered by WP6's alternatives A1–A6.
