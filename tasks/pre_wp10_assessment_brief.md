# Pre-WP10 assessment and next steps

> ## STATUS — EXECUTED 2026-07-30. All items discharged; WP10 drafted.
>
> Every checkbox below is closed. Reports are linked at each item. Two things
> in this brief turned out to be **wrong and are corrected where they appear**:
>
> 1. **D1's "median ≈ 9" is wrong.** 8.79 is the median of the *full* 54-branch
>    set; the retained 36-branch set's ensemble median is **13.29**. Dropping
>    α = 2.6 therefore *raises* the centre by ~50% rather than only trimming a
>    low tail — which is why the manuscript leads with the baseline branch
>    value and never an ensemble median.
> 2. **Item 3's "entire SN budget above ~52 M☉" is over-general.** True on the
>    baseline and coeval branches; on the 1–2 Myr formation-window branches the
>    lowest progenitor is 33.9 M☉. The defensible claim is **above 30 M☉ on
>    every branch** — which still lies above the whole Sukhbold island
>    structure, so the argument is unchanged. This was caught by pre-registered
>    prediction D1-P1, which **failed** and is recorded as failed.
>
> Two further things were found while executing: **B1's "agree, 5%" cross-check
> compared mismatched definitions** and is withdrawn as stated, and **CygOB2-B's
> counts-based age rails against the top of its prior grid**, making its
> supernova contribution a one-sided lower bound.
>
> Summary and links: [PROJECT_TRACE.md §1](../PROJECT_TRACE.md).


*Written 2026-07-30, after an independent read of the full WP0–WP9 chain
([PROJECT_TRACE.md](../PROJECT_TRACE.md), [paper1_execution_plan.md](../paper1_execution_plan.md),
[wp7_ledger.md](../reports/wp7_ledger.md), [wp8_crosschecks.md](../reports/wp8_crosschecks.md),
[wp9_verdict.md](../reports/wp9_verdict.md), [wp5_alpha_plausibility.md](../reports/wp5_alpha_plausibility.md),
[wp0_requirements_table.md](../wp0_requirements_table.md)). This brief is the
handoff for the agent executing the pre-WP10 fixes and then WP10 itself. It is
an assessment document — it changes no number and adopts no branch.*

---

## 1. Overall verdict

**The project is publishable and close to submission-ready.** The chain is
internally consistent where it matters — WP7, WP8 and WP9 all ran on the same
`repair_v7` chain — the framing decision (regular A&A article) was made by a
pre-registered mechanical rule, and the pre-registration/honesty discipline
(failed predictions recorded as failed: L1, L6, M2, F3; withdrawn results kept
on the books) is well above field standard.

What remains is a short list of reconciliations, one strategic branch decision,
one bounded piece of new work, and then manuscript assembly. Nothing below
requires new science beyond item T3.

## 2. The publishable results, ranked

1. **The pulsar result (WP8) — the paper's sharpest finding, bill it
   prominently.** PSR J2032+4127's companion MT91 213 is this project's *own
   census star* (orphan anchor `gaia_dr3:2067835682818358400`, B0V 17 M☉,
   positionally in CygOB2-A). The pulsar's existence excludes the branch on
   which N_SN = 0 everywhere (islands explodability: P(≥1 SN) = 0.000000 vs
   0.9997 all-explode). Existence, age (P = 0.960 within the widened 151–401 kyr
   window) and kinematics (27.6 km/s peculiar, as a bound binary requires) all
   agree with a ledger built without looking at the pulsar. The (a)/(b)/(c)
   degeneracy — high-mass explodability / older population / binary stripping —
   is honestly stated and must stay.

2. **The ledger (WP7).** N_SN = 8.43 baseline (PARSEC, R_V 3.1, α 2.3, coeval,
   all-explode); branch spread reported as the real uncertainty;
   P(last SN < 100 kyr) = 0.552; first explosion 1.30 Myr ago; rate ~8 Myr⁻¹.
   Knödlseder+02 / Martin+10 comparisons done via input differences, not
   hand-waving.

3. **The structural finding: the entire SN budget lives above ~52 M☉.** For any
   BH threshold ≤ 40 M☉ the ledger returns exactly zero on every branch. Cyg
   OB2's SN history is entirely conditional on very-massive-star explodability;
   the Sukhbold 15–25 M☉ island structure is irrelevant at this age. Clean,
   quotable, generalizes to any ~4 Myr association, and is what makes the
   pulsar test sharp.

4. **The Härer verdict, correctly conditional (WP9).** P_verdict = 0.323–0.736
   (median 0.533) over 36 headline branches; supported on every branch under
   Härer's own permissive "few hundred kyr" window (0.727–0.854); unresolved at
   their preferred ~50 kyr; hinges almost entirely on α = 2.0 vs 2.3 (axis
   spread 0.264 vs ≤0.065 for everything else). C3 = 1.000 — every ledger SN
   came from a >52 M☉ stripped-envelope progenitor, exactly Härer's required
   Ib/c channel. Energy is a stated conditional, never a probability.

5. **Secondary results worth their own subsections:** the census closure test
   as out-of-sample validation of the IMF extrapolation (closes at α ≈ 2.25,
   6.7% at Salpeter — this is the prepared answer to the strongest referee
   objection); the four-fold independent "CygOB2-C is different" signal
   (closure ratio 1.405, M2 failure, closing α 2.056, E1+E2 slope preference);
   single-population membership at 1.62 kpc with Berlanas two-distance not
   confirmed *within the member sample* (keep the circularity caveat — parallax
   was a WP2 clustering feature, issue #8).

## 3. Blocking reconciliations (fix before WP10 assembly)

- [x] **DONE** → [reports/wp5_association_mass_reconciliation.md](../reports/wp5_association_mass_reconciliation.md) · **B1 — Two association masses circulate and are not reconciled.**
  WP5 baseline reports **29,122 M☉** ("within factor 2 of the 16,500 M☉
  literature scale"); the WP6 cross-check table reports **1.74 × 10⁴ M☉**
  ("agree, 5%" with Wright/Härer's 1.65 × 10⁴). Presumably different
  integration ranges or footprints, but the same phrase "association mass"
  covers both with opposite margins. Find the definitions in
  [wp5_imf_norm_repair_v6.md](../wp5_imf_norm_repair_v6.md) /
  [cross_checks/harer_2025_supernova_rate.md](../cross_checks/harer_2025_supernova_rate.md),
  write one reconciliation paragraph, and make the paper define each quantity
  explicitly. **This is the one item I treat as blocking.**

- [x] **DONE** → [reports/wp4_wp5_age_reconciliation.md](../reports/wp4_wp5_age_reconciliation.md) · **B2 — CygOB2-B's age appears as both 3.548 Myr and 4.07 Myr.**
  WP4-adopted (repair_v5 upper-MS MAP): A 3.981, B **3.548**, C 2.512.
  Obligation O1 in PROJECT_TRACE §10 item 6: "A 4.00, B **4.07**, C 2.52".
  If 4.07 is the WP5 joint truth-age posterior and 3.548 the upper-MS MAP, say
  so where O1 is discharged; otherwise one is stale. The paper must quote one
  number per subgroup with a stated definition.

- [x] **DONE** → [PROJECT_TRACE.md](../PROJECT_TRACE.md) · **B3 — Stale status blocks in PROJECT_TRACE.** §7 header still reads
  "WP5 … **BLOCKED**" and the §1 narrative still says "WP5 is nevertheless
  still not accepted" — both superseded by the status board's PASS
  (repair_v6, re-passed under repair_v7). Mark superseded blocks as historical
  rather than current status.

## 4. Strategic decision for the PI (pre-register, then adopt)

- [x] **DONE** → [reports/wp7_alpha_headline_adoption.md](../reports/wp7_alpha_headline_adoption.md) · **D1 — Drop α = 2.6 from the headline range; keep α = 2.3-only OFF the
  table.** Per [wp5_alpha_plausibility.md](../reports/wp5_alpha_plausibility.md):
  dropping 2.6 is free (E1 alone supports it — 1/18 cells won, worst median χ²;
  E2 agrees without being spent) and cuts the reported N_SN spread from factor
  14.9 to 5.1 (**5.63–28.74, median ≈ 9**), removing the "maybe only 2 SNe"
  tail. Restricting further to α = 2.3 would spend E2 — the project's only
  out-of-sample IMF validation — and is not worth it (standing decision,
  PROJECT_TRACE §10 item 2). The plausibility report requires the adoption to
  have **its own pre-registration**; α = 2.6 stays in the sensitivity table
  exactly as WP9 §7 already reports it (P = 0.142–0.250, reported not deleted).

## 5. Highest-value remaining science (≈ 1 day)

- [x] **DONE** → [reports/wp7_binary_bound.md](../reports/wp7_binary_bound.md) · **T3 — Convert binary mass transfer from "acknowledged" to "bounded".**
  This is the paper's largest genuine vulnerability. Three independent places
  point at the same unmodeled physics: WP6 measured f_bin ≈ 0.7 above 8 M☉;
  Härer's own comparison is BPASS-based and Ic-dominated (a binary channel);
  pulsar degeneracy reading (c) is binary stripping. A referee will ask for at
  least a magnitude estimate. No BPASS integration needed — a
  literature-scaled bracket (e.g., the fraction of stripped-envelope SNe at
  3–5 Myr that BPASS attributes to binary channels, applied as a crude ± on
  N_SN) defuses the most likely major-revision request. Do it as a labelled
  literature-based branch or a discussion-section bound, not a pipeline change.

## 6. Small loose ends

- [x] **DONE** → [reports/wp0_dedup_resweep_2026-07-30.md](../reports/wp0_dedup_resweep_2026-07-30.md) · **S1 — Re-run the WP0 dedup sweep immediately before submission** (WP0's
  own instruction; the 2026-07-20 sweep is stale and the Härer+25 citing
  literature has grown).
- [x] **DONE** → [reports/wp7_ledger.md](../reports/wp7_ledger.md) · **S2 — WP7 §5 arithmetic:** 322.4 retained + 54.9 runaways vs the 380.6
  ledger — the reconciliation exists in §6 but the sum shown isn't exact; one
  sentence closes it.
- [x] **DONE** → [scripts/wp10_inputs.py](../scripts/wp10_inputs.py) · **S3 — Ensure WP10 pulls versioned tables only.** The unversioned
  `tables/wp5_imf_norm.csv` / `wp5_imf_norm.md` are the frozen 0/54 pre-repair
  record (issue #2) and must not reach the manuscript.
- [x] **DONE** → [reports/wp3_obligations_discharge.md](../reports/wp3_obligations_discharge.md) · **S4 — Discharge remaining WP3 adoption obligations in the text:** O1
  revised SF history (two older subgroups + one younger — after B2 fixes the
  numbers), O2 anchor absolute-scale systematic (~0.5 mag), O3 B's 4-anchor
  calibration asymmetry (Gaia XP as the DR4-era remedy).

## 7. WP10 framing guidance

- **Regular A&A article** — the mechanically-derived outcome; do not re-litigate.
- **Abstract skeleton** (from WP8 §10, sharpened): *Cyg OB2 has produced of
  order 8 supernovae (5.6–29 across carried branches), the most recent probably
  within the last ~100 kyr; an independent compact remnant — whose companion is
  a star in our own census — confirms at least one successful explosion; the
  Härer et al. PeV-bubble scenario is supported under its permissive age window
  and unresolved at its preferred ~50 kyr, hinging on the high-mass IMF slope —
  which Gaia DR4 will not settle, but a second association analysed identically
  would.*
- **Do not let P ≈ 0.53 read as a null result.** WP9 devil's-advocate
  objection 5 has the rebuttal: the analysis returns 0.055 / 0.53 / 0.84 for
  7 kyr / 100 kyr / permissive windows, and 0.55 vs an ignorance baseline of
  0.077 is a **7× measurement**. Put that comparison near the abstract, not
  buried in a rebuttal.
- **Energy condition:** WP0 requirement D7 asked for a joint statement
  including energy; WP9 carries it as a stated conditional. The paper must say
  explicitly that energy is conditional and why (progenitor masses are
  measured, explosion energies are not) — and that the energy and progenitor
  requirements point the same way.
- **The pre-registration discipline is itself a contribution.** One discussion
  paragraph: every marker comparison was pre-registered, markers were frozen at
  WP1 before WP5–WP7 existed, so "did you tune to the markers?" has a provable
  answer.
- **Keep the DR4 asymmetry** (WP9 §9): what DR4 fixes *and what it doesn't*
  (the α hinge). It reads as credibility.
- **Name the open items as Paper 2 / DR4 work, don't hold Paper 1 for them:**
  parallax-blind membership test, CygOB2-C distance-contamination test
  (PROJECT_TRACE §10c), C's 0-vs-7 SNe family ambiguity.

## 8. Suggested execution order — as executed

Executed in this order on 2026-07-30, with WP10 assembled last.

1. B1, B2, B3 (reconciliations — half a day)
2. D1 (pre-register + adopt the α = 2.6 drop — an hour plus the prereg record)
3. T3 (binary bound — up to a day)
4. S2, S3, S4 (as WP10 setup)
5. WP10 assembly per [paper1_execution_plan.md](../paper1_execution_plan.md)
   §WP10, using the regular-article structure
6. S1 (dedup re-run) + WP0 §5 disclosure protocol (Vink first, Brian second)
   immediately before submission

## 9. What remains before submission

Nothing in this brief. What is left is outside it and is listed in
[manuscript/README.md](../manuscript/README.md):

- **compile the manuscript** — no LaTeX toolchain exists in the working
  environment and `aa.cls` is not vendored, so `main.tex` has never been run
  through TeX. `scripts/wp10_validate.py` stands in for the checks a compile
  would make (undefined macros, unresolved refs, missing citations and
  figures), and it passes, but it is not a compile;
- **author list and affiliations**, currently placeholders;
- **Appendix A's branch/gate table**, stubbed;
- **re-run S1's dedup sweep** immediately before arXiv posting — the 2026-07-30
  sweep cannot cover the window between now and submission;
- **the disclosure protocol**: Vink first with the draft, Brian second with the
  verdict and draft, before arXiv.
