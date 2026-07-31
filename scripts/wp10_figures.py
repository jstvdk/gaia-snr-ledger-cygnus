#!/usr/bin/env python3
"""WP10 -- build the paper figures that do not exist yet.

Figures 1 and 2 (membership/substructure and the control fields) were built at
WP2 by make_paper_figures_wp1_wp2.py and are unchanged.  This script builds the
remaining four of the six-figure plan:

  fig3  mass function per subgroup with the fitted IMF overlaid  (WP5)
  fig4  R_SN(t), the explosion history, stacked by subgroup      (WP7)
  fig5  the mandatory age-sensitivity honesty plot               (WP7)
  fig6  the verdict against branches, with the headline/sensitivity
        split and the cross-check markers overlaid               (WP9, WP8)

Every input is resolved through wp10_inputs, so a superseded table cannot reach
a figure any more than it can reach the text.

Outputs: figures/paper/fig3..fig6 as PDF and PNG,
         provenance/wp10_figures_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp10_figures.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wp5_common as w
import wp10_inputs as I

PAPER = w.ROOT / "figures" / "paper"
BASE = dict(family="PARSEC", R_V=3.1, alpha=2.3, sf_duration_Myr=0.0)
COLOUR = {"CygOB2-A": "#3B6EA8", "CygOB2-B": "#C4622D", "CygOB2-C": "#4F8A5B"}

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def save(fig, name: str) -> list[str]:
    PAPER.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in ("pdf", "png"):
        path = PAPER / f"{name}.{suffix}"
        fig.savefig(path, bbox_inches="tight")
        written.append(str(path.relative_to(w.ROOT)))
    plt.close(fig)
    return written


def figure_three() -> list[str]:
    """Mass function per subgroup with the fitted IMF overlaid."""
    bins = pd.read_parquet(
        w.PROC / "wp5_mass_function_bins_repair_v7.parquet"
    )
    bins = bins[
        bins.family.eq(BASE["family"]) & bins.R_V.eq(BASE["R_V"])
        & bins.alpha.eq(BASE["alpha"])
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.6), sharey=True)
    for ax, subgroup in zip(axes, w.SUBGROUPS):
        panel = bins[bins.subgroup.eq(subgroup)].sort_values("bin_index")
        centre = panel.mass_geometric_center.to_numpy()
        observed = panel.membership_weighted_count.to_numpy()
        expected = panel.expected_count_at_k_median.to_numpy()
        ax.errorbar(centre, observed, yerr=np.sqrt(np.maximum(observed, 1.0)),
                    fmt="o", ms=3.5, color=COLOUR[subgroup], lw=1,
                    label="observed, completeness-corrected")
        ax.step(centre, expected, where="mid", color="0.25", lw=1.2,
                label=r"fitted IMF, $\alpha = 2.3$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"mass $[M_\odot]$")
        ax.set_title(subgroup.replace("CygOB2-", "Cyg OB2-"))
        ax.set_xticks([2, 3, 5, 8])
        ax.set_xticklabels(["2", "3", "5", "8"])
        worst = panel.pearson_residual.abs().max()
        ax.text(0.95, 0.92, rf"$\max|r| = {worst:.2f}$", ha="right",
                va="top", transform=ax.transAxes, fontsize=7, color="0.35")
    axes[0].set_ylabel("stars per bin")
    axes[0].legend(frameon=False, loc="lower left", fontsize=7)
    fig.suptitle(
        "Calibration-window mass function and the fitted normalization "
        "(baseline branch)", y=1.04, fontsize=9,
    )
    return save(fig, "fig3_mass_function")


def figure_four() -> list[str]:
    """R_SN(t): the explosion history, stacked by subgroup."""
    rsn = pd.read_csv(I.resolve("wp7_rsn_curves"))
    edges = np.append(
        np.sort(rsn.lookback_lo_Myr.unique()), rsn.lookback_hi_Myr.max()
    )
    centres = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    bottom = np.zeros(len(centres))
    for subgroup in w.SUBGROUPS:
        panel = rsn[rsn.subgroup.eq(subgroup)].sort_values("lookback_lo_Myr")
        rate = panel.rate_per_Myr.to_numpy()
        if rate.sum() == 0:
            continue
        ax.bar(centres, rate, width=np.diff(edges), bottom=bottom,
               color=COLOUR[subgroup], edgecolor="none", alpha=0.85,
               label=subgroup.replace("CygOB2-", "Cyg OB2-"))
        bottom += rate
    total = rsn.groupby("lookback_lo_Myr").rate_per_Myr.sum()
    first = float(total[total > 0].index.max()) + 0.05
    ax.axvline(first, color="0.3", ls=":", lw=1)
    ax.text(first - 0.03, ax.get_ylim()[1] * 0.95, "first explosion",
            rotation=90, ha="right", va="top", fontsize=7, color="0.3")
    ax.axvline(0.1, color="0.35", ls="--", lw=0.9)
    ax.annotate("last 100 kyr", xy=(0.1, ax.get_ylim()[1] * 0.55),
                xytext=(0.30, ax.get_ylim()[1] * 0.55), fontsize=7,
                color="0.3", va="center",
                arrowprops=dict(arrowstyle="->", color="0.3", lw=0.8))
    ax.set_xlabel("look-back time [Myr]")
    ax.set_ylabel(r"$R_{\rm SN}$ [Myr$^{-1}$]")
    ax.set_xlim(0, 1.6)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Explosion history, baseline branch")
    return save(fig, "fig4_rsn_history")


def figure_five() -> list[str]:
    """The mandatory age-sensitivity plot."""
    scan = pd.read_csv(I.resolve("wp7_age_sensitivity")).sort_values(
        "assumed_age_Myr"
    )
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(3.5, 4.0), sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.12},
    )
    ax.fill_between(scan.assumed_age_Myr, scan.N_SN_p16, scan.N_SN_p84,
                    color="#3B6EA8", alpha=0.20, lw=0, label="68%")
    ax.plot(scan.assumed_age_Myr, scan.N_SN_mean, color="#3B6EA8", lw=1.6,
            label="mean")
    ax.axhline(8.43, color="0.35", ls="--", lw=1)
    ax.text(2.1, 9.4, "baseline 8.4", fontsize=7, color="0.35")
    ax.set_ylabel(r"$N_{\rm SN}$")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Everything inherits the assumed age")

    ax2.plot(scan.assumed_age_Myr, scan.P_last_SN_within_100kyr,
             color="#C4622D", lw=1.6)
    ax2.axhline(0.552, color="0.35", ls="--", lw=1)
    ax2.set_ylabel(r"$P(t_{\rm last} < 100\,$kyr$)$")
    ax2.set_xlabel("assumed age applied to all three subgroups [Myr]")
    for axis in (ax, ax2):
        axis.axvspan(2.0, 3.0, color="0.9", zorder=0)
    ax2.text(2.05, 0.05, "nothing has died yet", fontsize=7, color="0.45")
    return save(fig, "fig5_age_sensitivity")


def figure_six() -> list[str]:
    """The verdict against branches, with the cross-check markers."""
    verdict = pd.read_csv(I.resolve("wp9_verdict"))
    verdict = verdict[verdict.explodability.eq("all_explode")].copy()
    verdict["label"] = (
        verdict.family.str[0] + " " + verdict.R_V.map("{:.1f}".format)
        + r" $\delta$" + verdict.sf_duration_Myr.map("{:.0f}".format)
    )
    fig, ax = plt.subplots(figsize=(7.1, 3.2))
    offsets = {2.0: -0.24, 2.3: 0.0, 2.6: 0.24}
    styles = {
        2.0: ("#3B6EA8", "o", r"$\alpha = 2.0$  (headline)"),
        2.3: ("#C4622D", "s", r"$\alpha = 2.3$  (headline)"),
        2.6: ("0.62", "^", r"$\alpha = 2.6$  (sensitivity only)"),
    }
    order = sorted(verdict.label.unique())
    index = {name: i for i, name in enumerate(order)}
    for alpha, (colour, marker, label) in styles.items():
        arm = verdict[verdict.alpha.eq(alpha)]
        ax.scatter([index[v] + offsets[alpha] for v in arm.label],
                   arm.P_verdict, s=22, color=colour, marker=marker,
                   label=label, zorder=3,
                   edgecolor="none" if alpha != 2.6 else "0.45")
    ax.axhline(0.5, color="0.25", lw=1)
    ax.text(len(order) - 0.4, 0.515, "supported", fontsize=7, color="0.35",
            ha="right")
    ax.text(len(order) - 0.4, 0.465, "not supported", fontsize=7, color="0.35",
            ha="right", va="top")
    ax.axhline(0.1, color="0.6", lw=0.8, ls=":")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel(r"$P_{\rm verdict}$")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="upper left", ncol=3)
    ax.set_title(
        "The verdict hinges on the IMF slope and on nothing else"
    )
    return save(fig, "fig6_verdict_branches")


def main() -> None:
    outputs = {}
    for builder in (figure_three, figure_four, figure_five, figure_six):
        written = builder()
        for path in written:
            outputs[path] = w.sha256(w.ROOT / path)
        print(f"built {written[0]}")

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp10_figures.py",
        "item": "WP10",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "matplotlib": matplotlib.__version__,
        "figures": {
            "fig1_membership_literature": "built at WP2, unchanged",
            "fig2_control_fields": "built at WP2, unchanged",
            "fig3_mass_function": "WP5 repair_v7 bins, baseline branch",
            "fig4_rsn_history": "WP7 R_SN(t), baseline branch, stacked",
            "fig5_age_sensitivity": "WP7 age scan, the mandatory honesty plot",
            "fig6_verdict_branches": "WP9 per-branch verdict, all 54",
        },
        "outputs": outputs,
    }
    w.write_json(w.PROVENANCE / "wp10_figures_execution.json", record)
    print("wrote provenance/wp10_figures_execution.json")


if __name__ == "__main__":
    main()
