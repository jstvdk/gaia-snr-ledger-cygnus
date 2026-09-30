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
| `manuscript/numbers.tex` | 3k | Densest file here: all 212 quoted quantities as macros, each with its source artifact in a trailing comment. |
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
4. **The chain runs on `repair_v7`.** Products from earlier repair rounds are
   preserved but must not be quoted.

## Withdrawn — never quote these

From the pre-`repair_v7` chain and the fixes to issues #16 and #17:

| withdrawn | correct |
|---|---|
| 45% closure excess, grid median 1.444 | 6.7% excess, grid median 1.067 |
| closing slope α = 2.070 | α ≈ 2.20 association-wide (A 2.34, B 2.25, C 2.06) |
| 260 raw / 109 corrected runaways | 119 raw / 54.9 corrected |
| living ledger 471.9 | 380.6 above 8 M☉ |
| association mass "agrees with Wright+15 to 5%" | 1.47× like-for-like (29,122 vs 16,500 M☉) |
| "retained 36-branch median ≈ 9" | 13.29 (8.79 is the *full* 54-branch median) |

## Headline numbers

Members 1,392 (1,331 subgroup-labelled: A 476, B 426, C 429) · ages A 4.00,
B 4.09, C 2.52 Myr · association stellar mass 29,122 M☉ · baseline deaths
**8.43** (A 4.17, B 4.26, C 0.00) · headline range 5.63–28.7 over 36 branches,
1.93–28.74 over all 54 · P(≥1 event) 0.9997 · P(last < 100 kyr) 0.552 · every
progenitor above ~34 M☉.

Always write `N_death`, or `N_SN | all explode` with the conditioning explicit.
The all-explode ledger is **not** an observed supernova history.

## Where it actually stands

WP0–WP9 complete, WP10 drafted and revised, WP11 Part B complete, WP12
complete. The open question is not a method one: the 12 Aug novelty audit found
that much of the emphasised result is already established (Wright+2015;
Menchiari+2024 got 7 ± 2.5 deaths at 3 Myr). The surviving claim — that
resolving the population into subgroups materially changes the death history,
because the young C subgroup sits below the first-death boundary — is
**a hypothesis under test, not a result**. The decisive experiment is the
pooled-vs-resolved (M0 vs M1) ablation with materiality thresholds frozen
before the M0 result is read; §12 of the audit is the work plan.

Everything downstream hinges on the IMF slope α: it moves the count by a factor
of five and splits the verdict 18/18 at α = 2.0 against 0/18 at α = 2.3.

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
