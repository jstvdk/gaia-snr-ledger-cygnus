# Task B — WP4 closure, schema repair and provenance commit

**Project:** `/Users/vdk/science/gaia_snr_history_cygnus`
**Status going in:** WP0–WP4 executed; WP4 self-reports `WP4_COMPLETE_GATE_PASSED`.
**Purpose of this task:** close WP4 properly before WP5 opens. This is a *repair and
verification* task, not a re-analysis. Do not re-derive membership, extinction or ages
unless a step below explicitly requires it.

Read first: `paper1_execution_plan.md` §1 (global conventions) and WP4/WP5;
`CUTS_AND_THRESHOLDS.md` §1 (number classes), §7.1, §9; `wp4_ages.md`;
`provenance/wp4_provenance.md`; `wp3_completion_report.md` (WP5 warning).

---

## 0. Findings already verified — do not re-litigate these

An external audit (`audit.py`, run 2026-07-23T15:33 under Python 3.11.15) confirmed the
following. Treat them as established; if your work contradicts any of them, stop and
report the contradiction rather than silently changing a number.

Row counts reconcile end to end:

```
wp2_members.parquet             2,112 rows (soft floor P>0.05)
  membership_probability >0.5 → 1,392    = gate table member count
  ..._astrometric        >0.5 → 1,331    = clean automatic members
tables/wp2_subgroup_labels      1,331 = A 476 + B 426 + C 429
  + 61 anchor_quality_exempt  → 1,392
wp3_* (all four products)       1,392
wp4_masses                      1,392 (476/426/429 + 61 unassigned)
wp1_wright15_gaia_crossmatch      167   (167/167 gate)
wp2_berlanas_recovery_audit       229   (189/229 = 82.5% denominator)
```

The WP4 branch grid is complete and symmetric: 44 rows per subgroup = 36 science branches
(2 families × 3 R_V × 3 f_bin × 2 indicators) + 8 distance-sensitivity rows at baseline
R_V. This predicts mean R_V = (12×3.0 + 20×3.1 + 12×3.5)/44 = 3.1818; the file reports
3.182. `dmu = ±0.0602` matches σ_μ = 0.060.

Distance chain is internally consistent: DM 11.054 ↔ 1.6245 kpc (5·log₁₀(162.45) = 11.0537).

`clump_class = "unresolved_binary"` (17 rows) is a genuine physical class, not a
placeholder. The audit's placeholder regex false-positived on the substring "unresolved".

---

## 1. Resolve the `subgroup_label` / `subgroup` collision  — BLOCKING

Two columns currently encode subgroup membership with opposite truth values:

| Column | Files | Content |
|---|---|---|
| `subgroup_label` | `data/processed/wp2_members.parquet` (2,112) | `CygOB2_distance_structure_unresolved`, **all rows** |
| | `data/processed/wp2_anchor_assignments.parquet` (190) | placeholder, all rows |
| | `data/processed/wp3_extinction.parquet` (1,392) | placeholder, all rows |
| | `data/processed/wp3_member_photometry.parquet` (1,392) | placeholder, all rows |
| `subgroup` | `data/processed/wp4_masses.parquet`, `wp4_anchor_hrd.parquet`, `wp4_age_posteriors.parquet` | real: `CygOB2-A/B/C`, `unassigned` |
| `subgroup_label` | `tables/wp2_subgroup_labels.parquet` (1,331) | real: A/B/C — the authoritative sidecar |

Nothing is numerically wrong today: WP4 joined the sidecar and wrote a new column name.
But WP5–WP7 will join across exactly these files, and a join that picks up
`subgroup_label` from a WP2/WP3 product will silently return the placeholder for every
star rather than failing.

**Do:** pick one name (recommend `subgroup`) and one source of truth
(`tables/wp2_subgroup_labels.parquet`). Then either populate it correctly in all WP2/WP3
products, or drop the stale column entirely from them. Do not leave both names in the
pipeline. Record the choice in `GLOSSARY.md §6 "Columns you added"` and in
`provenance/wp2_provenance.md`.

The 61 `anchor_quality_exempt` stars have no subgroup by construction (labels were derived
on the 1,331 clean automatic members only). Keep them as an explicit `unassigned`
category — do not impute a subgroup for them.

---

## 2. Rename `mass_best` and protect the branch structure — BLOCKING

`mass_best` in `wp4_masses.parquet` is byte-identical in summary statistics to
`mass_PARSEC_rv3.1` (mean 6.268, std 6.634, min 0.562, max 40.43). It is not a
combination of branches; it is one branch.

That is a defensible baseline, but the name invites WP5 to fit the IMF on a single column,
which would discard the Class-E branch structure that `paper1_execution_plan.md` §1.4 and
`CUTS_AND_THRESHOLDS.md` §1 require be carried to the end and never averaged.

The stakes are visible in the same file: `mass_parsec_mist_spread` reaches **23.25 M☉**,
and PARSEC tops out at 40.4 M☉ where MIST reaches 50.8 M☉. The branch disagreement is
concentrated at the top of the IMF — precisely the stars that set N_SN in WP7.

**Do:**
1. Rename `mass_best` → `mass_baseline` throughout (`scripts/wp4_masses.py`,
   `wp4_masses.parquet`, `tables/wp4_masses.cat`, `wp4_ages.md §8`, any WP5 stub).
2. Add a one-line note in `provenance/wp4_provenance.md §3` stating that
   `mass_baseline = PARSEC, R_V = 3.1` and that it is a reporting convenience with no
   downstream authority.
3. Confirm all six `mass_{PARSEC,MIST}_rv{3.0,3.1,3.5}` columns are retained and
   non-null-consistent (currently 55 NaN in each — verify it is the *same* 55 source_ids
   in all six, and record why they have no mass).

---

## 3. Re-derive the age envelope with the `measurable` filter applied — BLOCKING

Across the 132 rows of `wp4_age_posteriors.parquet`, `age_map` spans **1.123 to 10.12 Myr**,
against a headline of 3.5–4.5 Myr. Two specific defects:

- **A posterior railed the grid.** The age grid is 1–10 Myr log-spaced
  (`CUTS_AND_THRESHOLDS.md` #19). A MAP of 10.12 sits at or beyond the ceiling. A
  boundary-pinned posterior is not a measurement.
- **`n_stars` reaches a minimum of 1.** Some branch fitted an age from a single star —
  almost certainly subgroup C's PMS turn-on indicator, which `wp4_ages.md §5` already
  states is unmeasurable (~3 members in the faint window).

**Do:**
1. Verify `tables/wp4_ages_envelope.md` and the quoted 3.5–4.5 Myr headline are computed
   over rows with `measurable == True` only.
2. Report the count of excluded rows and the reason for each exclusion (not measurable /
   grid-railed / n_stars below a stated floor). Set and document that floor as a Class B
   number in `CUTS_AND_THRESHOLDS.md`.
3. For any row where `age_map` is within one grid step of 1 or 10 Myr, either extend the
   grid and refit or exclude it explicitly. Do not report a railed MAP as an age.
4. State the honest envelope in `wp4_ages.md §3` — the span over retained branches, not
   the baseline branch's credible interval.

---

## 4. Document the two branch asymmetries — BLOCKING for the gate text

**4a. Anchor coverage is severely lopsided across subgroups.** In `wp4_anchor_hrd.parquet`:

```
CygOB2-A     60
unassigned   48
CygOB2-C     43
CygOB2-B      5
```

The subgroups are near-equal in size (476/426/429), so subgroup B's spectroscopic HRD
indicator rests on **five stars** — yet `wp4_ages.md §10` records "**PASS** — B agrees" for
the criterion "both indicators agree per subgroup where measurable." Five stars is thin
evidence carrying a gate criterion.

Additionally, 48 of 156 anchors are `unassigned`, because subgroup labels were derived on
the 1,331 clean members and the anchor-exempt stars were excluded by construction. The
most informative stars are systematically outside the subgroup scheme. This propagates
into WP7, where per-subgroup ages become supernova clocks.

**Do:** restate the §10 gate row for B with the N made explicit, and add the anchor-count
imbalance to `provenance/wp4_provenance.md §6` (known limitations). If B's agreement does
not survive being stated as "5 anchors," downgrade the criterion to *documented
disagreement* rather than PASS. A downgrade is an acceptable outcome — see §6 below.

**4b. The two isochrone families were not treated symmetrically.** In the same file,
`age_used_PARSEC` varies 3.981–4.467 across anchors, while `age_used_MIST` has
std = 1.3e-15 — a single fixed value of 3.571 Myr for all 156 anchors. One family received
a per-subgroup age, the other a constant.

**Do:** either re-run the MIST anchor check per subgroup so the families are treated
identically, or document in `wp4_ages.md §1` and `provenance/wp4_provenance.md §3` exactly
why MIST uses a fixed age and what that does to the χ comparison. Do not leave it
undeclared — the anchor gate compares χ_PARSEC against χ_MIST, and the comparison is not
like-for-like as it stands.

---

## 5. Explain the nine astrometry-less rows

`wp2_members.parquet` has 2,112 rows but only 2,103 non-null in `ra`, `dec`, `l_deg`,
`b_deg`, `parallax_raw`, `parallax_corrected`, `parallax_error`, `pmra`, `pmdec`, `ruwe`,
`zero_point_boundary_flag`. Nine members carry a `membership_probability` with no
astrometry at all. They propagate into WP3 (1,392 vs 1,383 non-null).

**Do:** identify the nine `source_id`s, establish how a probability was assigned without
astrometry, and either document the mechanism in `provenance/wp2_provenance.md` or remove
them and re-run the affected counts. If removed, every count in §0 above changes and the
WP2 gate table must be regenerated.

Related, lower priority: `ruwe` reaches 24.03 in the member list. This is expected —
anchor-exempt stars bypass `quality_pass` — and is arguably correct given
`CUTS_AND_THRESHOLDS.md` §5.1 on RUWE's binary bias with ~70% O-star multiplicity.
Confirm the high-RUWE members are all `anchor_quality_exempt == True` and note it.

---

## 6. High proper-motion outliers — hand to WP6, do not silently discard

In `wp2_members.parquet`, `pmdec` reaches +4.895 against a mean of −4.31 with σ = 0.42
(≈22σ); `pmra` reaches +0.439 against a mean of −2.70 with σ = 0.34.

These may be contamination, but `paper1_execution_plan.md` WP6 has a runaway search as a
deliverable, and runaways are direct kinematic evidence of past supernovae.

**Do:** extract the outlier source_ids to `provenance/wp4_pm_outliers.csv` with their
membership probabilities, subgroup, anchor status and RUWE. Add a pointer in
`provenance/wp2_provenance.md` flagging them as WP6 runaway candidates. Do **not** cut
them from the member list in this task.

---

## 7. Regenerate stale reports

All three WP4 documents predate the products they describe:

```
15:05:04  wp4_ages.md
15:07:19  provenance/wp4_provenance.md
15:07:55  wp4_completion_report.md
15:12:35–15:12:47   every WP4 parquet, npz, figure, table and manifest
```

**Do:** after steps 1–6, re-run `scripts/wp4_report.py` and diff the regenerated documents
against the current ones. If any headline number moves, the WP4 gate must be re-assessed
against `paper1_execution_plan.md` WP4, not patched. Record the diff outcome in
`provenance/wp4_provenance.md §5`.

---

## 8. Commit and repair reference paths

`git HEAD` is `c0cf2af "WP2 - WP3"` (2026-07-23 11:39). Everything WP4 is untracked: seven
scripts, six provenance JSONs, four tables, four figures, both reports, `figures/wp4/`.
`paper1_execution_plan.md` §1.6–1.7 require per-artifact provenance and blocking gates; an
uncommitted WP4 has no anchor.

Separately, the paper library was renamed wholesale and is uncommitted: 15 PDFs show as
deleted, ~18 new names untracked. **`papers/ijaa_4501403.pdf` (Paíz et al. 2025) is
deleted and now exists as `papers/Paiz_2025.pdf`.** `provenance/wp2_paiz_crossmatch.csv`
was produced before the rename.

**Do:**
1. `grep -rn "ijaa_4501403\|0101509v1\|2201.05124\|2211.11625\|2308.01295\|2311.09089\|2402.07784\|2403.16650\|2508.21644\|2512.05854\|2603.27741\|aa38384-20\|aa55531-25" scripts/ provenance/ *.md`
   and repoint every hit to the new filename.
2. `papers/Berlanas_2020.pdf` and `papers/Berlanas_2020_v2.pdf` are both exactly
   4,453,133 bytes — same file twice. Keep one, delete the other, fix references.
3. Commit in three separate commits, in this order: (a) paper rename + reference repair,
   (b) schema repair from §1–§2, (c) WP4 closure — scripts, products, provenance, reports.
4. Add `audit.py` to the repo. Note in its docstring that it must be run under the project
   environment's Python; a bare Python 3.14 has no pandas and silently skips the entire
   data-products section.

---

## Acceptance criteria

WP4 is closed when all of the following hold:

1. Exactly one subgroup column name exists across all WP2/WP3/WP4 products, sourced from
   `tables/wp2_subgroup_labels.parquet`, with the 61 exempt stars explicitly `unassigned`.
   A fresh `audit.py` run reports **zero** placeholder hits in any `subgroup*` column.
2. `mass_best` no longer exists by that name; all six branch mass columns are retained;
   the 55 massless stars are the same set in every branch and their absence is explained.
3. The quoted age envelope is computed over `measurable == True` rows with no
   grid-railed MAP values, the exclusion count and reasons are stated, and the `n_stars`
   floor is registered as a Class B number.
4. Subgroup B's gate row states N = 5 explicitly, and the MIST fixed-age choice is either
   removed or documented in both `wp4_ages.md §1` and `provenance/wp4_provenance.md §3`.
5. The nine astrometry-less members are explained or removed; if removed, every downstream
   count and the WP2 gate table are regenerated.
6. PM outliers are extracted to `provenance/wp4_pm_outliers.csv` and flagged for WP6, with
   no stars cut from the member list.
7. All three WP4 documents postdate every WP4 product, and the regeneration diff is
   recorded.
8. `git status` is clean; no reference in `scripts/` or `provenance/` points at a deleted
   paper filename.

---

## Explicitly out of scope — do not do these

- Do not re-derive membership, subgroups, extinction or age posteriors from scratch. If a
  step above requires a refit, refit only that step and say so.
- Do not change the coeval verdict. `wp4_ages.md` reports A/B/C consistent at ≈3.5–4.5 Myr
  and invokes the plan's "if identical, state the coeval result" clause. That is compliant.
  If step 3 or 4 widens the envelope enough to threaten it, **report the tension; do not
  resolve it** — it belongs to WP9's verdict framing.
- Do not delete or overwrite `data/processed/wp2_members_failed_20260722.parquet`,
  `provenance/wp2_membership_manifest_failed_20260722.json`, or either
  `*_attempt_*_failure.json`. Frozen failures are deliberate.
- Do not revisit the Berlanas two-distance hypothesis. It was tested and not confirmed.
- Do not tune any threshold to make a gate pass. If a gate fails after repair, it fails,
  and that is the deliverable.

**Note on fragility:** Berlanas recall is 189/229 = 82.5% against a gate of ≥80% — six
stars of margin. If step 5 removes members, this is the criterion that breaks first. Check
it explicitly before declaring WP2 unaffected.

---

## Hand-off to WP5

When WP4 closes, WP5 (IMF normalization and completeness) opens with two constraints that
come out of this audit and should be written into its brief:

1. **Fit the IMF once per mass branch — all six columns — not on the baseline alone.**
   The PARSEC/MIST spread reaches 23 M☉ at the top of the IMF, which is where N_SN is set.
2. **Run injection-recovery separately per subgroup.** A/B/C differ in extinction and sky
   density, so a single association-wide selection function will bias the sparse subgroup.
   Watch C (weak age constraint, ~3 stars in the PMS window) and B (5 spectroscopic
   anchors, so any completeness correction there is nearly unvalidated).

The calibration window (2–8 M☉) is Class C per `CUTS_AND_THRESHOLDS.md` §7.1 and **must be
measured, not assumed** — §7.1 already warns it likely needs raising, and the injection
must run through the full chain (query cuts, quality cuts, extinction, membership
probability), not just the magnitude limit. `wp3_completion_report.md` carries a
truncation-bias warning and `figures/wp3/wp3_wp5_truncation.png`: read both before
starting, and do not rediscover it as a completeness bug.
