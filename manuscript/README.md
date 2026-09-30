# manuscript — build, provenance and known limitations

Regular A&A article. The framing was decided mechanically at WP9 (the
pre-registered rule returned INCONCLUSIVE), not chosen here.

**Revised 2026-08-03** against
[`tasks/manuscript_reframing_revision_brief.md`](../tasks/manuscript_reframing_revision_brief.md),
with the supporting analysis in
[`reports/wp12_manuscript_revision.md`](../reports/wp12_manuscript_revision.md).
The paper is a *measurement* of the massive-star and supernova history of
Cygnus OB2 first; the cocoon evaluation and the COSI forecast are applications
of that ledger, confronted in the Discussion, not the premise of the work.

## What the revision changed

| area | before | after |
|---|---|---|
| membership description | "clustering in (l, b, μα\*, μδ), parallax excluded" | two distinct operations: association selection **with** parallax, subgroup labelling **without** it |
| WP2 recall | "82.5 % recall" | automatic **76.2 %** and manual exceptions reported separately; control fields flagged as prior-calibrating, not fully external |
| WP5 gate | "40/54 pass, failures in the R_V = 3.5 / α = 2.6 corners" | the real `repair_v7` failure map (A 4, B 7, C 3; PARSEC 4, MIST 10), plus **15/36** headline branches built entirely on passing cells |
| closure | "The census closes at Salpeter" | "Association-wide closure is near Salpeter, **but the subgroups differ**"; closing slopes A 2.34, B 2.25, C 2.06; mixed-slope sensitivity added |
| WP7 convergence | not stated | the relative gate **failed as registered**; absolute material-cell drift reported in addition, not in substitution |
| explodability branch | "islands of implosion" | **hard direct-collapse cutoff**; the threshold scan is the primary result |
| pulsar | "settles that branch" | "**independently excludes an entirely explosion-free history**" |
| C3 | "= 1, all were type Ib/c" | a **deterministic model mapping**; bracketed to 0.68–0.96 at a 60 M☉ threshold, with primary sources |
| P_verdict | an unconditional probability | a **conditional scenario-availability score** shown *against* C4 |
| C4 thresholds | "below 0.60 moves any branch" | **0.720** moves the first, **0.581** moves all |
| WP8.5 neighbours | implied | attempted, its pre-registered gate **failed**, and the claim is explicitly scoped out |
| binary bound | "bounded" | a **conservative sensitivity bracket**, not a formal statistical bound |
| WP11 provenance | "committed, hashed, before execution" | described accurately: a timestamped local specification, no input hashes |
| Appendix A | empty heading | generated 54-branch table + electronic full version |
| figures | 6 | **11**, including the mandatory CMD/age figure that was missing |

## Files

| file | what it is |
|---|---|
| `main.tex` | the manuscript. **Contains no hand-typed numbers.** |
| `numbers.tex` | **generated** — every quoted quantity as a LaTeX macro |
| `tables_generated.tex` | **generated** — Tables 1–3 and Appendix A |
| `references.bib` | bibliography (29 entries, all cited) |
| `aa.cls`, `aa.bst` | A&A class and bibliography style, vendored |

## Build

The whole chain is one script:

```sh
bash scripts/run_manuscript_chain.sh
```

which runs, in order: input authorization and audit → the revision analyses →
the generated tables → `numbers.tex` → the figures → source validation →
the file inventory → PDF compilation with `latexmk` or, when it is unavailable,
Tectonic.

To run only the manuscript-facing steps:

```sh
PYTHONPATH=scripts python scripts/wp10_inputs.py     # authorize + audit inputs
PYTHONPATH=scripts python scripts/wp12_tables.py     # -> tables_generated.tex
PYTHONPATH=scripts python scripts/wp10_numbers.py    # -> numbers.tex
PYTHONPATH=scripts python scripts/wp12_figures.py    # -> figures/paper/fig01..11
PYTHONPATH=scripts python scripts/wp10_validate.py   # source validation
cd manuscript && tectonic --keep-logs --keep-intermediates main.tex
```

Environment: `provenance/environment_cygob2-gaia.yml` (full, pinned) or
`provenance/environment_cygob2-gaia_from-history.yml` (the declared
dependencies only).

## Current build status

- `wp10_validate.py` — **PASS**, all seven checks, 212 macros defined and 212
  used, 29 bibliography entries all cited, and 11 figures all present.
- Tectonic 0.17.0 — **PASS**, 21-page PDF produced with no undefined
  references/citations and no overfull boxes.
- Poppler 26.08.0 — all 21 pages rendered and visually inspected; **PASS**.
- Final PDF SHA-256:
  `c80d2928a95291920954d0516a726cf38ec045803e3b9f7d536e81979feed474`.
- The exact command, tool versions, checksums and inspection result are recorded
  in `provenance/manuscript_build_execution.json`.

## The rule that makes this reproducible

**Never type a number into `main.tex`.** Add it to `scripts/wp10_numbers.py`,
where it is read from a versioned artifact resolved through
`scripts/wp10_inputs.py`, and use the macro. `wp10_validate.py` check V7 fails
the build on any bare number in running text that is not a whitelisted
definition or literature value, and check V2 fails on any macro that stops being
used — so the text cannot silently drift away from the pipeline in either
direction. Checks V1–V5 see `main.tex` *plus* `tables_generated.tex`, because
that is what the reader sees; V7 deliberately sees `main.tex` alone, since a
generated table is numbers by construction.

Inputs are resolved through `wp10_inputs.py` rather than opened directly. It
declares the authorized (version-suffixed) artifacts and a forbidden list of
superseded ones — chiefly the pre-repair WP5 products in which 0 of 54 branches
passed the mass-function gate, which would contradict every accepted number in
the chain. Asking for one raises.

The WP12 analyses go one step further: `wp12_common.frozen()` verifies each
input against the SHA-256 recorded in `provenance/wp12_revision_prereg.json`,
which was written **before** any WP12 result existed, and raises if the file
moved.

## Figures

| figure | content | built by |
|---|---|---|
| 1 | analysis flow and evidence hierarchy | `wp12_figures.py` |
| 2 | membership and kinematic subgroups | `wp12_figures.py` |
| 3 | control fields and membership calibration | `wp12_figures.py` |
| 4 | de-reddened CMD and the age evidence | `wp12_figures.py` |
| 5 | extinction map and calibration support | `wp12_figures.py` |
| 6 | mass-function fits and the 54-cell residual-gate map | `wp12_figures.py` |
| 7 | massive-star closure by subgroup and slope | `wp12_figures.py` |
| 8 | branch-resolved supernova history | `wp12_figures.py` |
| 9 | age and explodability sensitivity | `wp12_figures.py` |
| 10 | conditional scenario score against C₄ | `wp12_figures.py` |
| 11 | isotope forecast and yield uncertainty | `wp12_figures.py` |

The earlier six-figure set (`fig1`–`fig6`, built by `wp10_figures.py` and
`make_paper_figures_wp1_wp2.py`) is **preserved on disk and no longer
referenced** by the manuscript, per the project's nothing-is-overwritten rule.

Colour rules are fixed across the whole set and not varied per panel: each
subgroup keeps one hue everywhere, each headline IMF slope keeps one hue,
sequential quantities use one perceptually uniform ramp, and nothing is encoded
by colour alone. The categorical hues are the Okabe–Ito colourblind-safe set,
checked for adjacent-pair separation under deuteranopia, protanopia and
tritanopia.

## Still to do before submission

1. **Author list and affiliations** — placeholders in `main.tex`.
2. **Re-run the prior-art deduplication sweep** immediately before submission
   (`reports/wp0_dedup_resweep_2026-07-30.md` is the current one).
3. **Archive the data tree** — the ignored ~8.4 GB working tree needs a DOI'd
   archival location with hashes, or the recorded queries need to be verified as
   sufficient to regenerate it.
4. **Disclosure protocol** — Vink first (draft), Brian second (verdict + draft),
   before arXiv; arXiv posting simultaneous with journal submission.
5. **Decide on notebooks as a release requirement.** The project does *not*
   currently have one executed, commented notebook per work package, and must
   not claim to. If notebooks are required, WP4 and WP6–WP12 need dedicated ones
   written.
