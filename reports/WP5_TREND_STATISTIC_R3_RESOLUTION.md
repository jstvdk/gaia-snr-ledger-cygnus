# Issue #11 resolved: the trend-gate criterion was measuring the wrong thing

*2026-07-28. Pre-declaration: [CUTS_AND_THRESHOLDS.md §14.6–14.7](../CUTS_AND_THRESHOLDS.md).
Evidence: [wp5_verdict_stability_execution.json](../provenance/wp5_verdict_stability_execution.json),
[wp5_verdict_stability.py](../scripts/wp5_verdict_stability.py).
Supersedes the R3 row of §14.4 and the "R3 stability FAIL" line in issue #11.*

---

## 1. The criterion was unsatisfiable, not merely strict

R3 as written required that on the four identical-model cell pairs of issue #11,
**both** Monte-Carlo realizations return the same gate verdict. The replacement
statistic failed it on one cell, and I recorded that as a defect in the
statistic. That reading was wrong.

For a cell whose model gives it probability π of passing, two independent
injection realizations disagree with probability **2π(1 − π)**. A test that
never disagrees would need π ∈ {0, 1} for every cell — its verdict would have to
be a deterministic function of the model, with no sampling variation at all. No
test of finite data has that property. **R3 was unsatisfiable by any candidate
statistic**, and the cell it failed on was simply the one whose true p sits
nearest 0.05, where 2π(1 − π) is largest.

R3 is therefore withdrawn and replaced by three criteria that test the defect
that is actually present. R1, R2 and R4 are untouched and their recorded
outcomes stand.

## 2. R3a — the real defect: a forbidden region at the threshold

With n = 6 the Spearman p-value lives on a lattice. Simulating null data from
every baseline cell's own fitted λ and collecting the achievable p-values:

| statistic | widest gap in [0.01, 0.20] | gap straddling 0.05 |
|---|---:|---:|
| incumbent rank test | 0.0454 | **0.0308** (0.0416 → 0.0724) |
| replacement slope test | 0.0003 | **0.00005** |

**No p-value exists between 0.042 and 0.072.** A cell landing there has no
verdict the data can support — it is not near a threshold, it is on a cliff with
nothing underneath. The replacement is continuous in the residuals and its
achievable values are dense. **R3a PASS** (limit 0.005).

## 3. R3b — the flip rate, and an independent conviction of the incumbent

Each cell's pass-probability π is computed from the fit's **own** model of
injection noise: the Dirichlet posterior over the response matrix's category
counts, already propagated into every published k. 400 response replicates are
drawn, the normalization is refitted on each by the production Jeffreys rule
against the **unchanged** observed counts, and the full three-way gate is
re-evaluated. No new parameter, no new noise model.

The implied flip probability 2π(1 − π) then predicts how often the four
identical-model pairs should disagree:

| statistic | observed flips | expected | 95% interval | verdict |
|---|---:|---:|---|---|
| incumbent | **4** | 1.43 | [0, 3] | **INCONSISTENT** |
| replacement | **1** | 0.74 | [0, 2] | consistent |

This is the sharpest single result of the investigation, and it was not
guaranteed — the π model could have predicted four flips and exonerated the
incumbent. Instead it convicts it: the incumbent flipped **more often than its
own injection uncertainty can explain**. The excess is exactly the lattice of
§2, which converts small residual changes into large p-value jumps. **R3b PASS.**

## 4. R3c — 41–44% of the branch grid has no determinate verdict

Cells with 0.05 < π < 0.95 are labelled **indeterminate**: their verdict is a
coin flip on the injection realization.

| version | grid | determinate pass | determinate fail | **indeterminate** |
|---|---:|---:|---:|---:|
| repair_v3 | 25/54 | 16 | 15 | **23 (43%)** |
| repair_v4 | 27/54 | 19 | 11 | **24 (44%)** |
| repair_v5 | 37/54 | 27 | 5 | **22 (41%)** |

R3c's requirement — compute π, label, report — is met. But §14.7(5)'s reporting
obligation is **triggered**, and it must be stated plainly rather than absorbed:

> **The 54-cell branch grid is underpowered at `N_INJECT_PER_MASS = 400`.
> Roughly two cells in five have a verdict that is noise. Per-cell verdicts must
> not be read as independent evidence, and the grid count must be reported with
> the indeterminate set broken out.**

This is a finding about the *grid*, not about the statistic, and it is the
deeper form of issue #11: the original observation ("gate flips in
provably-identical cells") was a symptom of injection statistics, and replacing
the statistic removes the lattice but not the underlying noise.

Note the trend across versions is nevertheless real and in the right direction:
determinate **failures** fall 15 → 11 → **5** while determinate passes rise
16 → 19 → **27**. The improvement from `repair_v4` to `repair_v5` is not an
artifact of indeterminacy.

## 5. Consequence for gate G3

Under §14.7 a cell indeterminate in **both** versions counts as neither a
regression nor an improvement, symmetrically and for every version.

| comparison | strict reading | refined (§14.7) |
|---|---|---|
| repair_v3 → repair_v4 | 1 regression — **BLOCK** | 0 — **PASS** (1 discounted) |
| repair_v4 → repair_v5 | 3 regressions — **BLOCK** | **1 — BLOCK** (2 discounted) |

The two CygOB2-A cells that regress at trend p = 0.040 and 0.048 are
indeterminate in both versions (π = 0.46 → 0.30 and 0.79 → 0.78) and are
correctly discounted. **One regression survives:**

**CygOB2-C MIST R_V = 3.5 α = 2.0** — determinate pass in `repair_v4`
(π = 1.000) and indeterminate in `repair_v5` (π = 0.900, max|r| = 3.25). It is
*not* discounted, because it is not indeterminate in both.

That cell is the node-snapping discontinuity of **issue #13**: its fitted age did
not move at all, but a 0.022 Myr shift of its 68% upper edge deleted a truth-age
node carrying 41.6% of the weight. Resolving R3 was never going to clear it, and
the pre-registered `repair_v6` fix
([wp5_node_interpolation_prereg.json](../provenance/wp5_node_interpolation_prereg.json))
is what decides it.

**WP5 remains unaccepted and WP6 remains unauthorized.** What blocks it is now a
single cell with a single, diagnosed, already-implemented cause.

## 6. A decision that belongs to the project owner

§14.7 refines a reading of gate G3 that was chosen explicitly on 2026-07-27
("strict reading — stay blocked"). The refinement is pre-declared, symmetric,
retroactive to every version, and it **does not change the current verdict** —
`repair_v4 → repair_v5` blocks under both readings. But it would matter if a
future comparison differed between them, so every G3 evaluation from here reports
both, and WP6 authorization is not flipped on the refined reading alone without
that difference being stated.

## 7. Reproduction

```bash
WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
PYTHONPATH=scripts python3 scripts/wp5_verdict_stability.py \
    --versions repair_v3 repair_v4 repair_v5 --replicates 400
```
