# repair_v6 — truth-age isochrone interpolation, and the acceptance of WP5

*2026-07-28. Pre-registration:
[wp5_node_interpolation_prereg.json](../provenance/wp5_node_interpolation_prereg.json),
written before any injection existed. Outcome:
[wp5_node_interpolation_outcome.json](../provenance/wp5_node_interpolation_outcome.json).
Gate record: [wp5_repair_v6_gate.json](../provenance/wp5_repair_v6_gate.json) ·
[wp5_completion_report_repair_v6.md](../wp5_completion_report_repair_v6.md).*

**WP5 is accepted. `downstream_wp6_authorized: true`.** This is the first
version to clear gate G3, and it does so under the strict per-branch reading —
the refined reading of §14.7 was not needed.

---

## 1. The defect

`repair_v4` introduced truth-side age marginalization: instead of generating
injected photometry at a single age, WP5 spreads the truth age over nine
equal-probability quantiles of the WP4 posterior. Those nine quantiles were then
**snapped to the nearest tabulated isochrone age**.

The isochrone grid is coarse — 0.05 dex, about 12% in age — so the snap
collapses nine distinct ages onto two or three, and it does so
**discontinuously**. An arbitrarily small change in the WP4 posterior can push a
quantile across a snap boundary and delete a whole node.

This was not hypothetical. It is what produced the one genuine
`repair_v4 → repair_v5` regression:

> CygOB2-C MIST R_V = 3.5 α = 2.0 kept an **identical** fitted age
> (2.002 → 2.002 Myr) and moved only −0.059 mag in A_V, but its 68% upper edge
> fell 2.084 → 2.062 Myr — **0.022 Myr** — which deleted the 2.248 Myr node
> holding **41.6%** of the weight. Bin 0's residual jumped +2.27 → **+3.25**.

The truth and recovery sides were also **inconsistent with each other**:
`wp4_repair_common._interpolate_age_sequence` already interpolated between
tabulated ages on the recovery side while the truth side snapped.

## 2. The discontinuity, measured

Translating each branch's WP4 posterior by a continuous offset and measuring the
1-Wasserstein distance between the resulting node distributions
([wp5_node_rule_continuity_execution.json](../provenance/wp5_node_rule_continuity_execution.json)):

| rule | max W₁ move per 0.0005 Myr step | amplification | node count within a branch |
|---|---:|---:|:--:|
| snapped (v4/v5) | **0.0546 Myr** | **109×** | 1 – 4 |
| unsnapped (v6) | 0.0005 Myr | **1.0×** | 9 – 9 |

The unsnapped rule is **Lipschitz-1**: it moves by exactly the shift that caused
it, which is the tightest continuity any rule can have.

## 3. The fix

The nine quantiles stay where the posterior puts them, each with weight 1/9, and
`wp5_common.load_isochrone_between_ages` builds the truth isochrone there using
**the same bracketing and linear age blend the recovery side already used**. The
two sides now discretize the same posterior the same way.

**No new parameter.** The node count is still `N_AGE_NODES = 9`, inherited from
the recovery side. At a node that happens to land on the grid the interpolated
loader reproduces the tabulated isochrone to **0.00e+00 mag**, verified at all
35 native ages in both families.

**Cost:** 37 → **162** node injections, 2.71 h. No reuse was possible — every
node now sits at an age no previous run had injected at.

**Backward compatibility is proven, not asserted**
([wp5_v6_backward_compatibility.json](../provenance/wp5_v6_backward_compatibility.json)):
the default node rule still reproduces the exact repair_v4 and repair_v5 node
sets found on disk; the interpolating loader equals the native table at every
native age; the vectorized Jeffreys `k` is bitwise identical to the scalar loop;
and the interpolation flag is off for every pre-v6 version.

## 4. Pre-registered predictions, all four confirmed

| | prediction | measured | verdict |
|---|---|---|---|
| P1 | the node rule is continuous | Lipschitz-1, ratio 1.0 | ✓ *(established before the run)* |
| P2 | the blocking cell clears, max\|r\| < 3.0 | **3.25 → 2.94** | ✓ *(does not gate adoption)* |
| P3 | baseline passes all three subgroups under both statistics, B's \|T\| < 2 | ✓ / ✓, \|T\| = **1.43** | ✓ |
| P4 | the change is local | age ≤ **0.118** Myr, mass **−0.22%**, grid +3 | ✓ |

**Adoption rested on P3 and P4**, which bound the *size* of the change without
prescribing its direction, and therefore cannot be satisfied by moving the
result toward passing. **P2 — the cell that motivated the investigation — was
deliberately excluded from the adoption rule** and is reported as an outcome
rather than a criterion. It confirmed anyway, and by the predicted mechanism:
restoring the deleted node pulled bin 0 from +3.25 to +2.94.

P4's age-shift bound deserves a note: the limit was 0.15 Myr and the worst
branch came in at **0.118 Myr**. That passes, but not by a wide margin, and it is
reported rather than buried.

## 5. Result

**The baseline barely moved**, which is the reassuring part — a nine-node
mixture reproduced what a two-node mixture had found, so `repair_v5`'s physics
was not an artifact of the coarse approximation:

| subgroup | residuals (v5) | residuals (v6) |
|---|---|---|
| CygOB2-A | +1.23 +1.34 −0.36 −0.33 +0.11 +1.09 | +1.10 +1.41 −0.37 −0.29 +0.17 +1.16 |
| CygOB2-B | +0.87 +2.27 +0.36 −0.24 −0.03 +0.17 | +0.87 +2.43 +0.42 −0.27 −0.09 +0.13 |
| CygOB2-C | +1.50 +0.55 +0.74 −0.07 +0.21 +1.56 | +1.34 +0.50 +0.79 −0.05 +0.25 +1.57 |

| quantity | repair_v5 | repair_v6 |
|---|---:|---:|
| truth-age nodes, baseline A / B / C | 2 / 4 / 1 | **9 / 9 / 9** |
| k, CygOB2-A | 1746 | 1726 (−1.1%) |
| k, CygOB2-B | 1629 | 1639 (+0.6%) |
| k, CygOB2-C | 1891 | 1884 (−0.3%) |
| association mass | 29,185 M☉ | **29,122 M☉** (−0.22%) |
| branch grid | 38/54 | **40/54** |

## 6. An unanticipated benefit: the grid stopped being underpowered

Issue #11's analysis found that **41–44% of the 54-cell grid had an
indeterminate verdict** in every version — a cell whose pass/fail is a coin flip
under the injection experiment's own uncertainty. That triggered §14.7(5) as a
finding about the grid.

`repair_v6` collapses it:

| version | grid | determinate pass | determinate fail | indeterminate |
|---|---:|---:|---:|---:|
| repair_v3 | 25/54 | 16 | 15 | 23 (43%) |
| repair_v4 | 27/54 | 19 | 11 | 24 (44%) |
| repair_v5 | 37/54 | 27 | 5 | 22 (41%) |
| **repair_v6** | **40/54** | **34** | **11** | **9 (17%)** |

This was **not predicted** and is worth stating as a mechanism rather than a
happy accident. A nine-node mixture averages the forward response over nine
independent injection realizations instead of two, so the response's Dirichlet
uncertainty falls by roughly √(9/2). The fix bought a factor of about 2 in
effective injection statistics for free, and 17% is below the 25% limit — **the
underpowered-grid finding does not apply to the accepted version**, though it
stands for v3–v5 and is recorded as such.

## 7. Gate G3

| clause | result |
|---|---|
| baseline, all three subgroups, residual gate | **pass** |
| no CygOB2-A or CygOB2-C branch regression | **pass — `A_or_C_regressions: []`** |
| baseline mass within factor two of literature | **pass** |
| ≥ 50 sources per subgroup branch | **pass** |
| **G3** | **PASS** |

Both readings agree and neither discounts anything
(`readings_agree: true`, `discounted_as_indeterminate_in_both: 0`). The strict
per-branch reading chosen on 2026-07-27 stays binding and was sufficient on its
own; §14.7's refinement changed nothing here, and the decision it raises does
not have to be taken under pressure.

**`accepted: true`, `blocking_reason: null`, `downstream_wp6_authorized: true`.**

## 8. What WP5 hands to WP6

- **Normalization** `k` per subgroup and branch, with the 54-branch spread as
  the systematic.
- **Association stellar mass 29,122 M☉** (baseline), within a factor two of the
  literature 16,500 M☉.
- **The bright-mass response**, which WP6 must apply per subgroup rather than as
  a scalar — see [CUTS §16](../CUTS_AND_THRESHOLDS.md). Effective completeness
  above 8 M☉ is **0.872** (grid median), and CygOB2-B's **0.962** against A's
  0.834 and C's 0.846 is a real per-subgroup systematic, not noise.
- **Carried caveats:** the anchor absolute-scale systematic (O2), CygOB2-B's
  4-anchor calibration asymmetry (O3), and the 14 branches that still fail
  under the §13 retention policy.

## 9. Reproduction

```bash
bash scripts/run_repair_v6_chain.sh
```

Every repair_v1–v5 artifact is preserved byte-identical.
