# Cygnus OB2 — Gaia DR3 massive-star and supernova history

A Gaia DR3 census of Cygnus OB2 (d = 1.62 kpc), resolved into three kinematic
subgroups, normalised through the IMF, and converted into a branch-resolved
ledger of how many massive stars have already died and when. Work packages
WP0–WP12; every step has a pre-registered gate written before its results
existed, and failed gates stay recorded as failed.

## ⚠️ Never read AUDIT.txt

It is **16 MB / ~4.3 million tokens** — a generated checksum inventory, not
prose. Reading it, or grepping `*.txt` at the repo root, destroys the context
window. Regenerate it with `audit.py`; never open it.

## Read these first

| file | ~tok | what it gives you |
|---|---|---|
| `PROJECT_TRACE.md` | 31k | The navigation spine. §1's status board carries all 13 work packages with gate verdict and headline numbers inline. |
| `manuscript/numbers.tex` | 3k | Densest file here: all 216 quoted quantities as macros, each with its source artifact in a trailing comment. |
| `manuscript/README.md` | 2k | Build state, what the WP12 revision changed, figure list, limitations stated plainly. |
| `reports/novelty_prior_art_and_incremental_value_audit_2026-08-12.md` | 7k | The most recent work and the honest account of where this stands. **Read before making any novelty claim.** |
| `cross_checks/README.md` | 0.3k | The external-validation rules. |

Then by task: `CUTS_AND_THRESHOLDS.md` (every cut and why, 15k) ·
`GLOSSARY.md` (11k) · `method_explained.md` (why the method is shaped this way,
13k) · `reports/wp7_ledger.md` + `tables/wp12_appendix_a_full.csv` (the ledger
and the 54-branch grid) · `wp2_subgroups.md` (membership) ·
`wp6_closure.md` (closure and runaways).

## Rules that are not enforced by code

1. **Never type a number into `manuscript/main.tex`.** Add it to
   `scripts/wp10_numbers.py`, where it is read from a versioned artifact, and
   use the macro. Check V7 in `wp10_validate.py` fails the build on any bare
   number in running text; V2 fails on any macro that stops being used.
2. **Resolve inputs through `scripts/wp10_inputs.py`**, never by opening files
   directly. It declares the authorized version-suffixed artifacts and refuses
   11 superseded ones by name. The WP12 analyses go further —
   `wp12_common.frozen()` verifies each input against the SHA-256 recorded in
   `provenance/wp12_revision_prereg.json` and raises if the file moved.
3. **Nothing is overwritten and nothing is retuned.** Superseded products stay
   on disk; a cross-check that disagrees becomes an issue in `PROJECT_TRACE.md`
   §9, never a reason to move a number. Cross-checks are validations, not
   calibrations.
4. **The chain runs on `repair_v9`** (since 2026-10-07, issue #21; selected in
   `scripts/chain.py`). WP4 ages are the Phase A′ test-c spectroscopic-HRD
   posteriors (A, C; B borrows A + C, flagged not measured) in
   `wp4_age_posteriors_repair_v9_headline.parquet`. repair_v8 and earlier
   products are preserved but must not be quoted; WP6–WP12 tables carry the
   `_repair_v9` suffix, and `wp10_inputs` refuses the superseded ones.

## Withdrawn — never quote these

From the pre-`repair_v7` chain, the fixes to issues #16 and #17, issue #19
(repair_v7 → repair_v8: anchor masses read at pre-repair ages) and issue #21
(repair_v8 → repair_v9: photometric WP4 ages replaced by spectroscopic ones):

| withdrawn | correct (repair_v9) |
|---|---|
| N_death 8.36 (A 4.10, B 4.26, **C 0.00**); earlier 8.43 | **6.66** (A 2.09, B 2.85, C 1.72) |
| "C is young / below the first-death boundary"; C 2.51–2.52 Myr | C 3.55 Myr (PARSEC) / 4.01 (MIST), coeval with A; 1.72 deaths |
| headline 5.57–28.7 (36 branches); 54-branch 1.94–28.74 | 6.29–34.5; 2.02–34.53 |
| P(last < 100 kyr) 0.547 | 0.697 |
| WP4 ages A 3.98 / B 3.55 / C 2.51 (photometric upper MS); WP5 4.01 / 4.09 / 2.52 | WP4 3.16 / 3.55 / 3.55 (spectroscopic; B borrowed); WP5 3.48 / 3.70 / 3.38 |
| 2.00–4.01 Myr upper-MS envelope | 3.16–4.01 Myr adopted WP4 age envelope (MAP) |
| WP5 40/54 cells, 25/36 headline cells | 28/54 (C 6/18), 18/36 |
| closure A/B/C 1.006/1.068/1.360, grid median 1.073; slopes A 2.29, B 2.25, C 2.05 | 0.924/1.036/1.498, 1.049; A 2.33, B 2.27, C 2.00 (association 2.19) |
| living ledger 388.6; association 29,014 M☉; progenitors above ~34 M☉ | 386.5; 28,949 M☉; above ~33 M☉ |
| 45% closure excess, grid median 1.444; later 1.067 | grid median 1.049 (repair_v9) |
| closing slope α = 2.070; later A 2.34 | α ≈ 2.19 association-wide (A 2.33, B 2.27, C 2.00) |
| 260 raw / 109 corrected runaways | 119 raw / 54.9 corrected |
| living ledger 471.9; later 380.6 | 386.5 above 8 M☉ |
| association mass "agrees with Wright+15 to 5%"; later 1.47× | 1.45× like-for-like (28,949 M☉ total) |
| "retained 36-branch median ≈ 9"; later 13.29, 13.16 | 15.95 (9.57 is the *full* 54-branch median) |
| 2.25–5.67 Myr "two-indicator" age envelope; any PMS age | no PMS row survives; see the repair_v9 rows above |
| A baseline closure 0.865 ("A falls short") | 0.924 (repair_v9) |

## Headline numbers

Chain `repair_v9` (adopted 2026-10-07). Members 1,392 (1,331 subgroup-labelled:
A 476, B 426, C 429) · WP4 ages (spectroscopic HRD, PARSEC) A 3.16, B 3.55
(borrowed from A + C, not measured), C 3.55 Myr; WP5 counts-based ages A 3.48,
B 3.70, C 3.38 · association stellar mass 28,949 M☉ · baseline deaths
**6.66** (A 2.09, B 2.85, C 1.72) · headline range 6.29–34.5 over 36 branches,
2.02–34.53 over all 54 · P(≥1 event) 0.997 · P(last < 100 kyr) 0.697 · every
progenitor above ~33 M☉ · WP5 **28/54** cells (C 6/18: C's low-mass counts
want ≤ 3.2–3.6 Myr, issue #22), gate G3 **fails** its no-regression clause.

Always write `N_death`, or `N_SN | all explode` with the conditioning explicit.
The all-explode ledger is **not** an observed supernova history.

## Where it actually stands

WP0–WP9 complete, WP10 drafted and revised, WP11 Part B complete, WP12
complete; repair_v9 adopted 2026-10-07
([reports/repair_v9_completion_report.md](reports/repair_v9_completion_report.md)).
The 12 Aug novelty audit found that much of the emphasised result is already
established (Wright+2015; Menchiari+2024 got 7 ± 2.5 deaths at 3 Myr). Its one
surviving claim — that resolving the population changes the death history
because a young C sits below the first-death boundary — **does not survive
repair_v9**: on spectroscopic ages C is coeval with A (3.55 vs 3.16 Myr PARSEC)
and contributes 1.72 deaths; at equal ages the three subgroups' death curves
agree within 7 %. WP13 (pooled vs resolved, M0 vs M1, thresholds frozen before
M0 is read) still decides the paper's framing, now expecting "equivalent".

Open, issue #22: the WP5 low-mass counts disagree with the spectroscopic ages
(gate 28/54 cells; C passes only at ≤ 3.2–3.6 Myr; A and B rail at the top of
their age prior). The ages behind the ledger use Gaia G plus literature
spectral types, not Gaia photometry alone (Stage 1: the Gaia CMD does not
measure these ages).

Everything downstream still hinges on the IMF slope α: it moves the count by a
factor of five; the verdict score is above 0.5 on 18/18 branches at α = 2.0
and 7/18 at α = 2.3 (0/18 on repair_v8).

## Environment and entry points

Conda env `cygob2-gaia` (`provenance/environment_cygob2-gaia.yml`).

```sh
bash scripts/run_manuscript_chain.sh          # the whole chain

PYTHONPATH=scripts python scripts/wp10_inputs.py    # authorize + audit inputs
PYTHONPATH=scripts python scripts/wp12_tables.py    # -> tables_generated.tex
PYTHONPATH=scripts python scripts/wp10_numbers.py   # -> numbers.tex
PYTHONPATH=scripts python scripts/wp12_figures.py   # -> figures/paper/fig01..11
PYTHONPATH=scripts python scripts/wp10_validate.py  # source validation
```

The ~8.4 GB data tree under `data/` is gitignored and lives only on this
machine; it has no permanent archive.

## Known documentation inconsistency

`PROJECT_TRACE.md` §1 still says the manuscript is "STILL NOT COMPILED IN THIS
ENVIRONMENT". That is stale — `provenance/manuscript_build_execution.json`
records a clean Tectonic 0.17.0 build on 2026-08-04, 21 pages, no undefined
references, all pages visually inspected. Trust the provenance JSON.
