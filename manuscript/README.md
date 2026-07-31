# WP10 — manuscript

Regular A&A article. The framing was decided mechanically at WP9 (the
pre-registered rule returned INCONCLUSIVE), not chosen here.

**Reframed 2026-07-30 (PI decision).** The paper is a *measurement* of the
massive-star / supernova history of Cygnus OB2 first; the Härer et al.
PeV-bubble evaluation is one application, confronted in the Discussion
(`sec:verdict`), not the premise of the work. The Introduction lists all
candidate high-energy sources (winds, supernova-into-cavity, Cyg X-3
microquasar, γ Cygni SNR, PSR J2032+4127) and motivates the ledger as the
population property they all depend on. This implements, rather than
contradicts, WP9's rule: "lead with the ledger; the verdict is one result
among several." No number changed.

## Files

| file | what it is |
|---|---|
| `main.tex` | the manuscript. **Contains no hand-typed numbers.** |
| `numbers.tex` | **generated** — every quoted quantity as a LaTeX macro |
| `references.bib` | bibliography |

## Build

```sh
PYTHONPATH=scripts python3 scripts/wp10_inputs.py     # authorize + audit inputs
PYTHONPATH=scripts python3 scripts/wp10_numbers.py    # -> manuscript/numbers.tex
PYTHONPATH=scripts python3 scripts/wp10_figures.py    # -> figures/paper/fig3..6
PYTHONPATH=scripts python3 scripts/wp10_validate.py   # the WP10 gate
cd manuscript && latexmk -pdf main.tex
```

**`aa.cls` is not vendored here.** Fetch it from
<https://www.aanda.org/for-authors> and place it beside `main.tex` before
compiling.

**No LaTeX toolchain is installed in this environment, so the manuscript has
not been compiled.** `wp10_validate.py` stands in for the failures a compile
would surface — undefined macros, unresolved cross-references, missing
citations, missing figures — but a real compile is still required before
submission.

## The rule that makes this reproducible

**Never type a number into `main.tex`.** Add it to `scripts/wp10_numbers.py`,
where it is read from a versioned artifact resolved through
`scripts/wp10_inputs.py`, and use the macro. `wp10_validate.py` check V7 fails
the build on any bare number in running text that is not a whitelisted
definition or literature value, and check V2 fails on any macro that stops
being used — so the text cannot silently drift away from the pipeline in
either direction.

Inputs are resolved through `wp10_inputs.py` rather than opened directly. It
declares the authorized (version-suffixed) artifacts and a forbidden list of
superseded ones — chiefly the pre-repair WP5 products in which 0 of 54 branches
passed the mass-function gate (issue #2), which would contradict every accepted
number in the chain. Asking for one raises.

## Figures

| figure | source | built by |
|---|---|---|
| 1 membership and substructure | WP2 | `make_paper_figures_wp1_wp2.py` |
| 2 control fields | WP2 | `make_paper_figures_wp1_wp2.py` |
| 3 mass function + IMF fit | WP5 `repair_v7` | `wp10_figures.py` |
| 4 R_SN(t) explosion history | WP7 | `wp10_figures.py` |
| 5 age-sensitivity honesty plot | WP7 | `wp10_figures.py` |
| 6 verdict against branches | WP9 | `wp10_figures.py` |

## Still to do before submission

1. **Author list and affiliations** — placeholders in `main.tex`.
2. **Appendix A** — the branch/gate table is stubbed; generate it from
   `tables/wp5_imf_norm_repair_v6.csv`, `tables/wp6_closure_repair_v7.csv`,
   `tables/wp7_ledger.csv` and `tables/wp9_verdict.csv`.
3. **Compile** with `aa.cls` and fix whatever only a real TeX run finds.
4. **Re-run the WP0 dedup sweep** immediately before submission
   (`reports/wp0_dedup_resweep_2026-07-30.md` is the current one).
5. **Disclosure protocol** — Vink first (draft), Brian second (verdict +
   draft), before arXiv; arXiv posting simultaneous with journal submission.
6. **Data and code availability statement** — the pipeline is intended to be
   public on release.
