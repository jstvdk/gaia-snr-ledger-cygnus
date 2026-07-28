> **SUPERSEDED 2026-07-27.** This report describes the frozen pre-repair WP5
> run (0/54 branches passing) and is retained as a historical record. The
> current WP5 verdict is in
> [wp5_completion_report_repair_v6.md](wp5_completion_report_repair_v6.md);
> the repair history is in [PROJECT_TRACE.md](PROJECT_TRACE.md) §8.

# WP5 completion report

**Verdict: BLOCKED — validation gate failed.**

WP5 was executed end to end and its diagnostic products are reproducible, but
it is not scientifically accepted.  No branch reaches the required 95%
end-to-end completeness, and all 54 response-aware Poisson IMF fits fail the
mass-function residual gate.  The baseline association mass is plausible
(about 25,000 Msun including the declared multiplicity convention), but that
integral sanity check does not repair the rejected mass-function shape.

The next authorized step is not WP6 or WP7.  It is a scoped revision of the WP4
2--5 Msun mass inference (posterior or direct CMD-space population model),
followed by a rerun of the preserved WP5 injection response.

See `wp5_imf_norm.md` for the scientific diagnosis and
`provenance/wp5_provenance.md` for exact reproducibility details.
