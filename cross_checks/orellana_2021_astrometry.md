# Cross-check — Orellana et al. 2021: systemic proper motion and distance

**Source.** R. B. Orellana, M. S. De Biasi, L. G. Páiz, *"New members of Cygnus
OB2 from Gaia DR2"*, MNRAS **502**, 6080–6093 (2021).
Local copy: `/Users/vdk/Downloads/Orellana_etal_2021.pdf`.

**What they did.** Gaia **DR2**, circular field of radius 1° centred on
(l, b) = (79.8°, +0.8°), G ≤ 17.5. Four overdensities identified in the vector
point diagram; a parametric proper-motion mixture model assigns membership; mean
parallax of the most probable members gives the distance.

**Why this matters to us.** It independently constrains two numbers the chain
leans on hard: the **systemic proper motion**, which WP6's runaway traceback
subtracts (issue #16), and the **distance**, which sets every absolute magnitude
in WP3/WP4 and therefore every mass in WP5.

---

## 1. Systemic proper motion — the strong result

| | Orellana 2021 (DR2) | this work (DR3) | difference |
|---|---|---|---:|
| μ_α* | −2.71 ± 0.02 mas/yr | **−2.7067** | **0.003** |
| μ_δ | −4.24 ± 0.02 mas/yr | **−4.3168** | 0.077 |

**μ_α* agrees to 0.003 mas/yr** with an entirely independent analysis on a
different data release using a different membership method.

This directly validates the fix for [issue #16](../PROJECT_TRACE.md). The first
WP6 runaway traceback used *absolute* proper motions and recovered the canonical
ejected star BD+43 3654 with probability **0.000**, because Cyg OB2's systemic
motion is larger than a typical ejection signature. The corrected traceback
subtracts the systemic motion — and the vector it subtracts is confirmed here
against an external source.

**The μ_δ difference is 0.077 mas/yr = 0.59 km/s at 1.62 kpc.** Formally larger
than their quoted statistical error, but negligible against the 10–100 km/s
ejection window the runaway search operates in. Plausible causes: DR2→DR3
systematics, or different member samples (they use a 1° circle to G ≤ 17.5; we
use a subgroup-labelled sample from a wider box). Carried as a small systematic,
not an issue.

---

## 2. Distance — a real 3% offset

| source | distance |
|---|---:|
| Orellana 2021, astrometric members (2767) | 1683 ± 5 pc |
| Orellana 2021, astro-photometric members (300 in common) | **1669 ± 6 pc** |
| Páiz & Orellana (priv. comm., G ≤ 13.5) | 1670 ± 7 pc |
| **this work (WP2 adopted)** | **1620 pc** |
| this work, median member parallax | 1615 pc |
| this work, inverse-variance weighted | 1629 pc |
| this work, latent-distance mixture fit | 1624.5 pc |

Our four internal estimates cluster at 1615–1629 pc. Orellana's cluster at
1669–1683 pc. **The offset is ~3%, or ~50 pc.**

### Almost certainly the parallax zero point

Orellana apply the **DR2 global offset of −0.029 mas** (Lindegren et al. 2018;
Lu ri et al. 2018). We apply the **DR3 per-star Lindegren et al. 2021**
correction. Their mean member parallax is 0.599 mas; ours is 0.614 mas. The
0.015 mas difference is exactly the scale of a DR2→DR3 zero-point revision.

This is a methodological difference, not an error on either side. Ours is the
more modern treatment; theirs was correct for DR2.

### Context — the literature spans 1330 to 2100 pc

Orellana's own Table 4:

| Morgan+54 | Reddish+66 | Humphreys & McElroy 84 | Massey & Thompson 91 | Comerón & Pasquali 12 | Kiminki+15 | Berlanas+19 | Lim+19 | Orellana+21 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1500 | 2100 | 1820 | 1700 | 1445 | 1330 ± 60 | 1760 | 1600 ± 100 | 1669 ± 6 |

**Our 1620 pc sits mid-range**, closest to Lim et al. 2019 (1600 ± 100). The 3%
offset from Orellana is small compared with the spread of published values.

### Propagated effect on WP6 — about 2.4%

Adopting 1669 pc instead of 1620 pc:

1. distance modulus rises by 5·log₁₀(1669/1616) = **0.070 mag**
2. at fixed apparent magnitude, stars become intrinsically brighter → more massive
3. on the upper main sequence (L ∝ M^3.5), M_bol = −8.75·log₁₀M, so this is
   0.0080 dex = **+1.86% in mass**
4. scaling every mass by (1 + f) is equivalent to lowering the 8 M☉ threshold by
   the same factor, so with dN/dM ∝ M^(−2.3) the count above it rises by
   (α − 1)·f = **+2.4%**

So every closure ratio would rise ~2.4%. That **helps CygOB2-A** (currently 0.894,
below unity) and **hurts CygOB2-C** (currently 1.405). Same order as other
carried systematics and well below the ~10% residual, but it belongs in the
systematics budget.

**Not adopted.** Changing the distance to match a cross-check would be tuning.
It is recorded as a systematic instead.

---

## 3. The other three overdensities — consistent, not contradictory

Orellana identify four VPD overdensities in the field:

| structure | members | distance |
|---|---:|---:|
| **Cygnus OB2** | 2767 | 1683 ± 5 pc |
| UCB585 (open cluster) | 8 additional | ~1460 pc |
| "right" overdensity | 179 | ~1280 pc |
| "left" overdensity | 188 | ~1280 pc |

These are **separate structures with distinct proper motions**, which is exactly
how Orellana separates them. Our WP2 membership clusters in
(`parallax_corrected`, `pmra`, `pmdec`) and would likewise assign them low
membership probability.

Their existence is therefore **not** evidence against our single-population
result — it is evidence that the field contains foreground structure that both
analyses remove. It is, however, directly relevant to the Berlanas question:
see [berlanas_2019_two_distance.md](berlanas_2019_two_distance.md).

---

## 4. Membership counts

| | Orellana 2021 | this work |
|---|---:|---:|
| members | 2767 | 2112 |
| magnitude limit | G ≤ 17.5 | G < 19 + quality cuts |
| footprint | 1° circle at (79.8, +0.8) | wider box, subgroup-labelled |

Same order of magnitude, but the selections differ enough that the counts are not
directly comparable and no conclusion is drawn from them.

Worth noting for WP2's own gate: Orellana find that **33 of 333** previously
published photometric/spectroscopic members are *not* astrometric members, 16 of
them O–B stars. Independent confirmation that published member lists for Cyg OB2
carry real contamination — the same reason our census weights by membership
probability rather than hard-cutting.

---

## 5. What a failure would have looked like

| failure mode | what we would have seen | actual |
|---|---|---|
| systemic PM wrong (issue #16 class) | μ disagreeing by ≫0.1 mas/yr | 0.003 in μ_α* |
| distance badly wrong | offset ≫10%, outside the literature spread | 3%, mid-range |
| membership sweeping in a foreground group | our distance pulled *toward* 1280 pc | not seen |

---

## 6. Verdict

**AGREE on proper motion**, at a precision that independently validates the
issue #16 fix.

**3% distance offset**, attributable to the DR2 versus DR3 parallax zero point,
propagating to ~2.4% on the WP6 closure ratios. Recorded as a systematic; nothing
adjusted.
