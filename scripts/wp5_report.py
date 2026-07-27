#!/usr/bin/env python3
"""WP5 figures, tables, blocked-gate report, and provenance narrative."""
from __future__ import annotations

import json
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wp5_common as w


COLORS = {
    "CygOB2-A": "#4477AA",
    "CygOB2-B": "#EE6677",
    "CygOB2-C": "#228833",
}


def completeness_figure(curves: pd.DataFrame) -> None:
    base = curves[curves["family"].eq("PARSEC") & curves["R_V"].eq(3.1)]
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
    for axis, subgroup in zip(axes, w.SUBGROUPS, strict=True):
        data = base[base["subgroup"].eq(subgroup)].sort_values("primary_mass")
        axis.plot(
            data["primary_mass"],
            data["query_fraction"],
            color="#999999",
            lw=1.8,
            label="Gaia query",
        )
        axis.plot(
            data["primary_mass"],
            data["quality_fraction"],
            color="#CCBB44",
            lw=1.8,
            label="+ WP2 quality",
        )
        axis.plot(
            data["primary_mass"],
            data["recovery_fraction"],
            "o",
            ms=3,
            color=COLORS[subgroup],
            alpha=0.65,
            label="full recovery (raw)",
        )
        axis.plot(
            data["primary_mass"],
            data["recovery_isotonic"],
            color=COLORS[subgroup],
            lw=2.2,
            label="full recovery (monotone)",
        )
        axis.axhline(0.95, color="black", ls="--", lw=1, label="95% target")
        axis.axvline(2.0, color="black", ls=":", lw=1)
        axis.set(
            title=subgroup,
            xlabel=r"true primary mass [$M_\odot$]",
            xlim=(0.5, 8.0),
            ylim=(0.0, 1.03),
        )
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("recovery fraction")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=5,
        frameon=False,
    )
    figure.suptitle(
        "WP5 catalogue-level completeness — PARSEC, $R_V=3.1$",
        y=0.98,
    )
    figure.tight_layout(rect=(0.0, 0.10, 1.0, 0.93))
    path = w.FIGURES / "wp5_completeness_curves.png"
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def mass_function_figure(
    normalization: pd.DataFrame, bins: pd.DataFrame
) -> None:
    base_n = normalization[
        normalization["family"].eq("PARSEC")
        & normalization["R_V"].eq(3.1)
        & normalization["alpha"].eq(2.3)
    ]
    base_b = bins[
        bins["family"].eq("PARSEC")
        & bins["R_V"].eq(3.1)
        & bins["alpha"].eq(2.3)
    ]
    figure, axes = plt.subplots(2, 3, figsize=(14, 7.2), sharex="col")
    for column, subgroup in enumerate(w.SUBGROUPS):
        data = base_b[base_b["subgroup"].eq(subgroup)].sort_values("bin_index")
        summary = base_n[base_n["subgroup"].eq(subgroup)].iloc[0]
        center = data["mass_geometric_center"].to_numpy(float)
        observed = data["membership_weighted_count"].to_numpy(float)
        expected = data["expected_count_at_k_median"].to_numpy(float)
        axes[0, column].errorbar(
            center,
            observed,
            yerr=np.sqrt(np.maximum(observed, 1e-12)),
            fmt="o",
            color=COLORS[subgroup],
            label="membership-weighted data",
        )
        axes[0, column].step(
            center,
            expected,
            where="mid",
            color="black",
            lw=2,
            label="forward IMF model",
        )
        axes[0, column].set_title(
            f"{subgroup}\n"
            + r"$p_{\chi^2}=$"
            + f"{summary['poisson_chi_square_p']:.1e}"
        )
        axes[0, column].set_yscale("log")
        axes[0, column].grid(alpha=0.2)
        axes[1, column].axhline(0.0, color="black", lw=1)
        axes[1, column].axhspan(-3.0, 3.0, color="#DDDDDD", alpha=0.4)
        axes[1, column].plot(
            center,
            data["pearson_residual"],
            "o-",
            color=COLORS[subgroup],
        )
        axes[1, column].set_xscale("log")
        axes[1, column].set_ylim(-5.0, 9.0)
        axes[1, column].grid(alpha=0.2)
        axes[1, column].set_xlabel(r"recovered mass [$M_\odot$]")
    axes[0, 0].set_ylabel("count per observed-mass bin")
    axes[1, 0].set_ylabel("Pearson residual")
    axes[0, 0].legend(frameon=False, fontsize=8)
    figure.suptitle(
        r"WP5 IMF residual gate — PARSEC, $R_V=3.1$, $\alpha=2.3$",
        y=0.99,
    )
    figure.tight_layout()
    path = w.FIGURES / "wp5_mass_function.png"
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def association_mass_figure(association: pd.DataFrame) -> None:
    figure, axis = plt.subplots(figsize=(10, 5.5))
    families = {"PARSEC": "o", "MIST": "s"}
    rv_offset = {3.0: -0.12, 3.1: 0.0, 3.5: 0.12}
    for family, marker in families.items():
        for rv in w.R_V_BRANCHES:
            data = association[
                association["family"].eq(family)
                & association["R_V"].eq(rv)
            ].sort_values("alpha")
            x = data["alpha"].to_numpy(float) + rv_offset[rv]
            y = data["multiplicity_adjusted_mass_median_Msun"].to_numpy(float)
            lo = data["multiplicity_adjusted_mass_lo68_Msun"].to_numpy(float)
            hi = data["multiplicity_adjusted_mass_hi68_Msun"].to_numpy(float)
            axis.errorbar(
                x,
                y,
                yerr=np.vstack([y - lo, hi - y]),
                fmt=marker,
                ms=6,
                capsize=2,
                label=f"{family}, $R_V={rv}$",
            )
    axis.axhline(
        w.LITERATURE_MASS_MSUN,
        color="black",
        lw=1.5,
        label="Wright+15 scale",
    )
    axis.axhspan(
        w.LITERATURE_MASS_MSUN / 2.0,
        w.LITERATURE_MASS_MSUN * 2.0,
        color="#BBBBBB",
        alpha=0.25,
        label="factor-two gate",
    )
    axis.set(
        xlabel=r"high-mass IMF slope $\alpha$",
        ylabel=r"multiplicity-adjusted association mass [$M_\odot$]",
        xticks=w.IMF_SLOPES,
        title="WP5 association-mass sensitivity",
    )
    axis.grid(alpha=0.2)
    axis.legend(ncol=2, fontsize=8, frameon=False)
    figure.tight_layout()
    path = w.FIGURES / "wp5_association_mass.png"
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame[columns].copy()
    def render(value):
        if isinstance(value, (float, np.floating)):
            return f"{value:.3g}" if np.isfinite(value) else "—"
        return str(value)

    header = "| " + " | ".join(columns) + " |"
    separator = "|" + "|".join("---" for _ in columns) + "|"
    rows = [
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in view.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def write_reports(
    curves: pd.DataFrame,
    normalization: pd.DataFrame,
    bins: pd.DataFrame,
    association: pd.DataFrame,
) -> None:
    base_curve = curves[
        curves["family"].eq("PARSEC") & curves["R_V"].eq(3.1)
    ]
    completeness_rows = []
    for subgroup in w.SUBGROUPS:
        data = base_curve[base_curve["subgroup"].eq(subgroup)]
        at_two = data.loc[data["primary_mass"].eq(2.0)].iloc[0]
        completeness_rows.append(
            {
                "subgroup": subgroup,
                "recovery_at_2_Msun": at_two["recovery_fraction"],
                "bright_6_8_median": data.loc[
                    data["primary_mass"].between(6.0, 8.0),
                    "recovery_fraction",
                ].median(),
                "max_monotone_recovery": data.loc[
                    data["primary_mass"].ge(2.0), "recovery_isotonic"
                ].max(),
            }
        )
    completeness = pd.DataFrame(completeness_rows)
    base_norm = normalization[
        normalization["family"].eq("PARSEC")
        & normalization["R_V"].eq(3.1)
        & normalization["alpha"].eq(2.3)
    ]
    base_assoc = association[
        association["family"].eq("PARSEC")
        & association["R_V"].eq(3.1)
        & association["alpha"].eq(2.3)
    ].iloc[0]
    first_bins = bins[
        bins["family"].eq("PARSEC")
        & bins["R_V"].eq(3.1)
        & bins["alpha"].eq(2.3)
        & bins["bin_index"].eq(0)
    ][
        [
            "subgroup",
            "membership_weighted_count",
            "expected_count_at_k_median",
            "pearson_residual",
        ]
    ]

    norm_table = normalization.copy()
    norm_table.to_csv(w.TABLES / "wp5_imf_norm.csv", index=False)
    (w.TABLES / "wp5_imf_norm.md").write_text(
        "# WP5 IMF-normalization branch table\n\n"
        + markdown_table(
            norm_table,
            [
                "subgroup",
                "family",
                "R_V",
                "alpha",
                "k_median",
                "k_lo68",
                "k_hi68",
                "poisson_chi_square_p",
                "residual_gate_pass",
            ],
        )
        + "\n",
        encoding="utf-8",
    )
    association.to_csv(w.TABLES / "wp5_association_mass.csv", index=False)

    report = f"""# WP5 — IMF normalization and completeness

**Status: BLOCKED AT THE VALIDATION GATE. Do not propagate these normalizations
to WP6/WP7 as accepted inputs.**

This work package was executed against the frozen WP2--WP4 products.  It
produced complete diagnostic artifacts and a forward-response IMF fit, but it
did not satisfy either the absolute-completeness requirement or the
mass-function residual requirement.

## 1. Frozen inputs and branch grid

The fit uses all three kinematic subgroups; PARSEC and MIST masses; R_V =
3.0/3.1/3.5; and high-mass IMF slopes alpha = 2.0/2.3/2.6.  Each observed star
contributes its WP2 membership probability.  The six logarithmic observed-mass
bins cover 2--8 Msun.  Unresolved binaries are injected with f_bin=0.4 and
q~U[0.1,1], matching the central WP4 binary branch.

The injection response extends in true primary mass from 0.5 to 8 Msun so that
lower-mass binaries and mass-estimation scatter into the observed 2--8 Msun
window are modeled rather than ignored.

## 2. Completeness experiment

Synthetic stars were generated with the branch-matched WP4 age and the same
frozen PARSEC/MIST isochrones.  Their extinction was drawn spatially from the
matching WP3 subgroup/R_V distribution.  Real Gaia DR3 observational states
were cloned from nearby sources of similar G and BP-RP in the one-degree field.
Every injection then passed the frozen G/parallax query, exact WP2 quality
filter, reconstructed cluster-versus-field posterior, P>0.5 handoff, and G/BP/RP
mass-readiness requirement.

The 128-point Sobol-normal membership integration agrees with the frozen WP2
10,000-draw decision for 98.95% of a stratified 760-source validation sample
(7 false negatives, 1 false positive; median absolute probability difference
0.0068).

Baseline completeness:

{markdown_table(completeness, list(completeness.columns))}

No branch/subgroup curve reaches 95% absolute end-to-end recovery anywhere in
2--8 Msun.  Gaia query recovery becomes essentially complete above about
2.25--3.25 Msun depending on subgroup, but the WP2 quality and membership stages
leave a bright plateau of only about 75--82%.  Therefore the planned absolute
95% lower edge does not exist.  The diagnostic 95%-of-plateau edges are retained
in the normalization table, but they are not relabelled as absolute
completeness.

Because the selection response is explicitly measured, a diagnostic corrected
fit over the nominal 2--8 Msun observed window remains mathematically possible.
It is labelled `corrected_no_absolute95_edge`, not "95% complete".

## 3. Forward Poisson IMF fit

The forward model predicts each *recovered* mass bin:

lambda_i = k integral R_i(m_true) m_true^(-alpha) dm_true,

where R_i contains query/quality/membership recovery, unresolved-binary
brightening, WP3 extinction error, photometric error, and migration through the
exact WP4 nearest-isochrone mass estimator.  A Jeffreys prior on k and a
Dirichlet posterior for the injection response yield 10,000 posterior draws per
subgroup and branch.

Baseline k posteriors:

{markdown_table(base_norm, ["subgroup", "k_median", "k_lo68", "k_hi68", "raw_calibration_sources", "membership_weighted_calibration_sources"])}

The first observed-mass bin is the decisive failure:

{markdown_table(first_bins, list(first_bins.columns))}

All 54 subgroup x family x R_V x alpha fits fail the predeclared residual gate.
The best branch has p_chi2 = {normalization['poisson_chi_square_p'].max():.3g};
the baseline subgroup p-values are
{", ".join(f"{row.subgroup}: {row.poisson_chi_square_p:.2g}" for row in base_norm.itertuples())}.
The response-aware model still underpredicts the 2.0--2.52 Msun bin by
6.4--7.9 sigma while overpredicting several adjacent bins.  This is not a
smooth completeness slope that can be repaired by changing k or alpha.

The most direct evidence points to an incompatibility between the frozen WP4
point-mass catalogue and a power-law birth IMF in this lower-CMD regime.  WP4
assigns photometric masses by a nearest point on a single-star isochrone; the
2--5 Msun catalogue therefore needs per-star mass posteriors or a direct
CMD-space IMF likelihood before the IMF gate can be trusted.  This is a
diagnosis, not permission to rewrite the frozen WP4 products.

## 4. Association-mass diagnostic

For the baseline PARSEC, R_V=3.1, alpha=2.3 branch, the summed normalization
implies a primary-system mass of
{base_assoc['primary_system_mass_median_Msun']:.0f}
[{base_assoc['primary_system_mass_lo68_Msun']:.0f},
{base_assoc['primary_system_mass_hi68_Msun']:.0f}] Msun and a
multiplicity-adjusted stellar mass of
{base_assoc['multiplicity_adjusted_mass_median_Msun']:.0f}
[{base_assoc['multiplicity_adjusted_mass_lo68_Msun']:.0f},
{base_assoc['multiplicity_adjusted_mass_hi68_Msun']:.0f}] Msun.
That median is within a factor two of the Wright+15 16,500 Msun scale.  Across
all 18 association branches the medians span
{association['multiplicity_adjusted_mass_median_Msun'].min():.0f}--
{association['multiplicity_adjusted_mass_median_Msun'].max():.0f} Msun, with
{int(association['within_factor_two_literature'].sum())}/18 branches inside the
factor-two band.

This sanity check passes at baseline but cannot override the failed shape
residuals: a biased mass function can integrate to a plausible total mass.

## 5. Gate assessment

| Criterion | Result |
|---|---|
| Lower calibration edge >=95% complete | **FAIL** — no absolute 95% edge exists in any branch/subgroup |
| >=50 sources per subgroup in the corrected 2--8 diagnostic window | **PASS** — 248--368 raw per branch/subgroup |
| Residuals consistent with Poisson scatter; no hidden mass trend | **FAIL** — 0/54 fits pass; baseline first-bin excess 6.4--7.9 sigma |
| Association mass within factor about two of literature | **PASS baseline; one extreme branch fails** |

**Blocking conclusion:** WP5 is not accepted.  WP6/WP7 must not consume
`wp5_imf_normalization.parquet` as a validated normalization.

## 6. Required remediation

1. Reopen only the WP4 lower/intermediate-mass inference, replacing nearest-point
   masses with per-star mass posteriors or a direct CMD-space population
   likelihood that carries binaries and extinction covariance.
2. Preserve the current injections and rerun the response matrix against the
   revised recovered-mass estimator.
3. Reassess whether the P>0.5 WP3/WP4 handoff should be extended to the soft
   P>0.05 catalogue; the current injection correction is usable, but it cannot
   manufacture missing upstream photometry/mass posteriors.
4. Accept WP5 only after an untuned residual gate passes and the completeness
   deviation is either removed or explicitly approved as a corrected-selection
   design change.

## 7. Outputs

- `data/processed/wp5_completeness_curves.parquet`
- `data/processed/wp5_injection_response.parquet`
- `data/processed/wp5_imf_normalization.parquet` (diagnostic; blocked)
- `data/processed/wp5_mass_function_bins.parquet`
- `data/processed/wp5_association_mass.parquet`
- `data/processed/wp5_imf_posterior_draws.npz`
- `tables/wp5_imf_norm.csv`, `tables/wp5_imf_norm.md`
- `tables/wp5_association_mass.csv`
- `figures/wp5/wp5_completeness_curves.png`
- `figures/wp5/wp5_mass_function.png`
- `figures/wp5/wp5_association_mass.png`
- `notebooks/wp5_imf_normalization_and_completeness.ipynb`
- `provenance/wp5_manifest.json`, `provenance/wp5_provenance.md`
"""
    (w.ROOT / "wp5_imf_norm.md").write_text(report, encoding="utf-8")

    completion = """# WP5 completion report

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
"""
    (w.ROOT / "wp5_completion_report.md").write_text(
        completion, encoding="utf-8"
    )

    injection_log = json.loads(
        (w.PROVENANCE / "wp5_injections_execution.json").read_text(
            encoding="utf-8"
        )
    )
    provenance = f"""# WP5 provenance — completeness and IMF normalization

**Generated:** {datetime.now(timezone.utc).isoformat()}  
**Status:** `WP5_BLOCKED_RESIDUAL_AND_ABSOLUTE_95_COMPLETENESS_GATES`  
**Downstream authority:** none; diagnostic products must not enter WP6/WP7.

## Pipeline

1. `scripts/wp5_common.py` reconstructs the frozen WP2 classifier, loads the
   branch-matched WP3/WP4 models, and defines the IMF/multiplicity integrals.
2. `scripts/wp5_injections.py` injects 400 stars at each 0.25-Msun point from
   0.5--8 Msun for every subgroup x family x R_V branch, clones real Gaia
   observational states, runs the frozen selection, and records recovered
   masses.
3. `scripts/wp5_fit_imf.py` fits all subgroup x family x R_V x alpha branches
   with a response-aware Poisson likelihood and 10,000 posterior draws.
4. `scripts/wp5_report.py` creates the figures, tables, and gate report.
5. `scripts/wp5_make_notebook.py` creates the explanatory notebook, which is
   executed in the `cygob2-gaia` environment.
6. `scripts/wp5_finalize.py` validates schemas/hashes and freezes the manifest.

## Frozen execution constants

| Quantity | Value | Class / role |
|---|---:|---|
| Injection seed | {w.SEED} | reproducibility |
| Mass grid | 0.5--8 Msun, step 0.25 | response support |
| Injections per mass/branch/subgroup | {w.N_INJECT_PER_MASS} | binomial precision |
| Binary fraction | {w.F_BINARY} | central WP4 binary branch |
| q distribution | U[{w.Q_MIN},1] | WP4 convention |
| Membership integration | {w.MEMBERSHIP_QMC_POINTS}-point Sobol normal | validated approximation |
| Published handoff | P>{w.MEMBERSHIP_THRESHOLD} | frozen WP3/WP4 selection |
| Absolute completeness target | {w.COMPLETENESS_TARGET} | required gate; failed |
| Observed calibration window | 2--8 Msun | corrected diagnostic fallback |
| Observed mass bins | {w.N_IMF_BINS} logarithmic | residual diagnostic |
| IMF slopes | {w.IMF_SLOPES} | mandatory branches |
| Posterior draws | {w.N_POSTERIOR_DRAWS} | k/response uncertainty |
| Total-mass range | 0.08--120 Msun | Kroupa-like integration |

## Injection validation

The QMC approximation was compared directly with the frozen WP2 10,000-draw
probabilities: decision agreement =
{injection_log['qmc_validation_against_frozen_wp2_10k']['decision_agreement']:.5f}
over {injection_log['qmc_validation_against_frozen_wp2_10k']['rows']} stratified
sources.  The injection is catalogue-level, not image-level: Gaia epoch images
and AGIS cannot be rerun.  That limitation is explicit in the execution log.

## Selection definition

Recovery requires all of:

1. synthetic observed G<19 and raw parallax in 0.35--1.10 mas;
2. exact WP2 quality state (RUWE<1.4, >=8 visibility periods, BP/RP excess
   present, positive finite covariance, Lindegren zero-point domain);
3. WP2 posterior-odds membership probability >0.5;
4. G/BP/RP present for a WP4 recovered mass.

The mass response includes photometric error, A_V error, unresolved binaries,
and the exact WP4 nearest-isochrone inverse mapper.  True masses down to 0.5
Msun are included so upward migration into the observed 2--8 Msun window is not
omitted.

## Statistical model

Observed bin counts are sums of membership probabilities.  The forward
intensity uses the injected true-to-recovered response.  A Jeffreys prior is
used for k; response rows receive a Jeffreys-multinomial Dirichlet posterior.
Primary-system mass uses a continuous two-part IMF (alpha=1.3 below 0.5 Msun;
branch alpha above), and the separately reported stellar mass adds companions
under the declared f_bin/q convention.

## Gate failure preserved

No absolute 95% completeness edge exists.  All 54 response-aware fits fail the
Poisson residual gate, including every mandatory slope/family/R_V branch.  No
threshold, bin edge, or branch was tuned after observing the failure.  The
diagnostic association mass is retained because it helps localize the failure,
but it has no downstream authority.
"""
    (w.PROVENANCE / "wp5_provenance.md").write_text(
        provenance, encoding="utf-8"
    )


def main() -> None:
    w.FIGURES.mkdir(parents=True, exist_ok=True)
    w.TABLES.mkdir(parents=True, exist_ok=True)
    curves = pd.read_parquet(w.PROC / "wp5_completeness_curves.parquet")
    normalization = pd.read_parquet(w.PROC / "wp5_imf_normalization.parquet")
    bins = pd.read_parquet(w.PROC / "wp5_mass_function_bins.parquet")
    association = pd.read_parquet(w.PROC / "wp5_association_mass.parquet")
    completeness_figure(curves)
    mass_function_figure(normalization, bins)
    association_mass_figure(association)
    write_reports(curves, normalization, bins, association)
    outputs = [
        w.FIGURES / "wp5_completeness_curves.png",
        w.FIGURES / "wp5_mass_function.png",
        w.FIGURES / "wp5_association_mass.png",
        w.TABLES / "wp5_imf_norm.csv",
        w.TABLES / "wp5_imf_norm.md",
        w.TABLES / "wp5_association_mass.csv",
        w.ROOT / "wp5_imf_norm.md",
        w.ROOT / "wp5_completion_report.md",
        w.PROVENANCE / "wp5_provenance.md",
    ]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_report.py",
        "status": "SUCCESS_REPORTING_BLOCKED_GATE",
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / "wp5_completeness_curves.parquet",
                w.PROC / "wp5_injection_response.parquet",
                w.PROC / "wp5_imf_normalization.parquet",
                w.PROC / "wp5_mass_function_bins.parquet",
                w.PROC / "wp5_association_mass.parquet",
                w.PROVENANCE / "wp5_injections_execution.json",
                w.PROVENANCE / "wp5_imf_fit_execution.json",
            ]
        },
        "outputs": {
            str(path.relative_to(w.ROOT)): {
                "sha256": w.sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in outputs
        },
    }
    w.write_json(w.PROVENANCE / "wp5_report_execution.json", record)


if __name__ == "__main__":
    main()
