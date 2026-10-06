# Issue #20, Phase A′: three independent age tests for A, B and C

*2026-10-06. Owner decision: [provenance/decisions_2026_10_02.json](../provenance/decisions_2026_10_02.json), item 4.
Pre-registration, committed in `99bbe90` before any test touched real stars:
[provenance/issue20b_prereg.json](../provenance/issue20b_prereg.json).
Results:
[tables/issue20b_age_tests.csv](../tables/issue20b_age_tests.csv) ·
[tables/issue20b_adopted_ages.csv](../tables/issue20b_adopted_ages.csv) ·
[provenance/issue20b_outcome.json](../provenance/issue20b_outcome.json).
Post-hoc diagnostic, not pre-registered:
[tables/issue21_error_model_diagnostic.csv](../tables/issue21_error_model_diagnostic.csv).
Read-only: no chain product changed.*

## 1. Outcome by the frozen rule

**Only one cell of the 18 is resolved: C in MIST at R_V 3.0 and 3.1, adopted at
4.45 Myr (68 % 3.89–4.51) from tests a and c.**

Everything else is NOT RESOLVED. The registered comparison of C with A is
therefore **undetermined**, because A is not adopted in either family.

| R_V 3.1 | test a: near-IR free-A_V CMD | test b: ESP-HS HRD | test c: spectroscopic HRD (outlier-robust) | status |
|---|---|---|---|---|
| A, PARSEC | 5.01 [4.97, 5.02] | (invalid) | 3.16 [2.80, 3.50] | NOT RESOLVED |
| A, MIST | 7.15 [7.02, 7.28] | (invalid) | 3.57 [3.31, 3.70] | NOT RESOLVED |
| B, PARSEC | 7.94 [7.91, 8.34] | (invalid) | (invalid: 4 stars) | NOT RESOLVED |
| B, MIST | 8.03 [7.60, 8.12] | (invalid) | (invalid: 4 stars) | NOT RESOLVED |
| C, PARSEC | 4.47 [4.46, 4.62] | (invalid) | 3.55 [3.32, 3.97] | NOT RESOLVED |
| C, MIST | 4.50 [4.20, 4.57] | (invalid) | 4.01 [3.76, 4.27] | **ADOPTED 4.45** |

Values are MAP [68 %] in Myr.

- **Test b failed its registered validity check** in every cell. On the 30
  calibration anchors, ESP-HS T_eff scatters by 0.087 dex against a 0.08 limit.
  It is recorded as failed; the threshold was not relaxed. For reference, its
  C posteriors were broad, at 1.7–6.6 Myr.
- **Test a passed its registered validity check.** The near-IR scale against the
  anchors is s = 0.93 with 0.30 mag scatter. However, its ages for A and B
  (5–8 Myr) are physically impossible (§2).
- **No two-component C.** With the outlier term and the grid-edge guard, the
  two-age mixture gives ΔBIC −2.6 (PARSEC) and −5.3 (MIST) for C. A's
  edge-of-grid component is now rejected by the guard. Phase A's formal H3 does
  not survive a model that allows outlier stars.
- **Integrity checks pass:**
  - J1: test a recovers synthetic 2.51 and 3.98 Myr populations exactly, in both families.
  - J2: test c with ε = 0 reproduces Phase A exactly.
  - J3: every input matched its pre-registered hash.

## 2. Why test a disagrees with everything (issue #21, post-hoc)

Test a's *control* uses the same machinery but WP3's extinction instead of the
near-IR extinction. It should have reproduced WP4. It did not: PARSEC A 7.08,
C 5.01, against WP4's 3.98 and 2.51.

The cause is isolated in [tables/issue21_error_model_diagnostic.csv](../tables/issue21_error_model_diagnostic.csv). That diagnostic:
- uses WP4's own stars, extinction and model particles;
- changes only how each star's A_V uncertainty enters the likelihood;
- reproduces every stored WP4 age exactly when WP4's error model is used.

| R_V 3.1, f_bin 0.4 | WP4 model (A_V error independent in G and colour) | A_V error correlated along the reddening vector | bright window M_G0 ≤ −1, WP4 model | bright window, correlated |
|---|---:|---:|---:|---:|
| A, PARSEC / MIST | 3.98 / 4.01 | 5.62 / 6.37 | 2.82 / 4.01 (broad) | 2.82 / 4.01 (broad) |
| B, PARSEC / MIST | 3.55 / 4.01 | 5.62 / 8.03 | 3.98 / 2.83 | 2.82 / 4.01 |
| C, PARSEC / MIST | **2.51** / 3.18 | 5.01 / 4.50 | **3.98** / 3.18 | **3.98** / 3.18 |

What this shows:

1. **The photometric ages depend on an approximation.**
   `wp4_common.branch_photometry` adds a star's A_V uncertainty to σ(M_G0) and
   σ(colour) as independent errors. A real A_V error moves the star *along* the
   reddening vector, so the two errors are fully correlated.

   The photometric stars' A_V errors are large (median 0.6–0.9 mag), so the
   choice matters. Propagating them correctly moves every subgroup to 4.5–8 Myr.
2. **Those correctly propagated full-window ages cannot be right.** A contains
   O5–O7 supergiants (O5 If + B0 V, O6 I + O5.5 III, O7 I + O6 I + O9 V). These
   cannot exist at 6–7 Myr, when the turnoff is 26–29 M☉ (`wp6_mass_extension_decision.turnoff_mass`); O5–O6 stars need about 40–60 M☉.

   So the full window (M_G0 ≤ 1.5) is mis-modelled in its *faint* part, whatever
   the error model. The cause is open: candidates are the pre-main-sequence
   (PMS) to main-sequence transition, binaries, contamination, or A_V outliers.
3. **The bright upper main sequence is stable.** For M_G0 ≤ −1, where the
   turnoff information is, the error model barely matters. In PARSEC, C comes
   out at 3.98 Myr, matching its spectroscopy (3.55–3.98). In MIST it comes out
   at 3.18. A's bright window is broad (68 % 1.1–2.8 Myr in PARSEC).
4. **C's 2.51 Myr is produced by the faint part of the window under the
   independent-error approximation.** It is not produced by C's bright stars
   (Phase A §3 already showed it is not produced by the supergiants either).

Test a is affected by the same faint-end problem: its likelihood is correctly
correlated by construction. Its tiny intervals (±0.05 Myr) are the usual sign
of a mis-specified model.

## 3. What can and cannot be said now

- **Not established: "C is young."** None of the measurements that are stable
  to the error model puts C below about 3.2 Myr. Those are the spectroscopic HRD
  (Phase A and test c) and the bright-window photometry.
- **Not established: the WP4 photometric ages of any subgroup.** A 3.98, B 3.55
  and C 2.51 Myr hold only under an error approximation that is wrong in
  principle. Every downstream product inherits them: WP5 truth-age priors, WP6
  turnoffs, the WP7 ledger and N_death 8.36.
- **Most robust measurement in the project:** the spectroscopic HRD for A
  (3.16 / 3.57) and C (3.55 / 4.01). B has only 4 spectroscopic stars.

## 4. Decision needed (issue #21)

The frozen rule leaves 17 of 18 cells unresolved. It registered a two-branch
fallback: keep the repair_v5 photometric posterior, plus test c's posterior.
That fallback would carry a photometric age now known to rest on an incorrect
error approximation.

| option | what it means | cost |
|---|---|---|
| **(a) Follow the registered fallback literally** | repair_v9 carries two age branches per subgroup: photometric repair_v5 and spectroscopic. B's spectroscopic branch rests on 4 stars | quick; but it keeps a known-defective branch alive |
| **(b) Pre-register a WP4 age redesign first** (repair_v9 WP4 step) | correlated A_V error propagation; a bright upper-MS indicator window (e.g. M_G0 ≤ −1), with the faint part's mis-modelling studied separately; an acceptance check that A's and C's bright-window ages agree with their spectroscopic HRD within 68 %. Then propagate WP5 → WP12 | about 2–4 days plus about 1–2 h of injections. This is the honest route to "use one age throughout" |
| **(c) Adopt the spectroscopic HRD ages for A and C now; B from its bright window** | ages from 43–59 stars; B on photometry alone | quick; a methods change that needs its own pre-registration, and a weaker basis for B |

Recommendation: **(b)**. The defect sits in WP4 itself, not just in C. A fix
pre-registered at the WP4 level, checked against the spectroscopy, is the only
route that gives one defensible age per subgroup for the whole chain.

Whichever option is chosen, nothing quoted today changes until a pre-registered
chain version replaces repair_v8.
