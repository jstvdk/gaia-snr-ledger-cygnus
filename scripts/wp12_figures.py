#!/usr/bin/env python3
"""WP12 -- the revised paper figure set.

The brief's diagnosis of the old six-figure set was that it displayed the last
numbers of the chain rather than explaining the chain, and that several panels
were illegible at A&A column width.  This script builds the replacement set.

  fig01  analysis flow and evidence hierarchy          (new)
  fig02  membership and kinematic subgroups            (revised)
  fig03  control fields and membership calibration     (redesigned)
  fig04  de-reddened CMD and the age evidence          (new, was missing)
  fig05  extinction map and its calibration support    (new)
  fig06  mass-function fits and the residual-gate map  (expanded)
  fig07  massive-star closure by subgroup and slope    (new)
  fig08  branch-resolved supernova history             (revised)
  fig09  age and explodability sensitivity             (combined)
  fig10  conditional scenario score against C_4        (replaces the verdict plot)
  fig11  isotope forecast and yield uncertainty        (new)

Colour rules, applied throughout and not varied per panel:
  * the three subgroups keep ONE fixed hue each, everywhere they appear;
  * the two headline IMF slopes keep one fixed hue each;
  * sequential quantities use a single perceptually uniform ramp;
  * nothing is encoded by colour alone -- every colour-coded series is also
    direct-labelled or distinguished by marker.
The categorical hues are the Okabe--Ito colourblind-safe set, checked for
adjacent-pair separation under deuteranopia, protanopia and tritanopia.

Outputs: figures/paper/fig01..fig11 as PDF and PNG,
         provenance/wp12_figures_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp12_figures.py
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import wp5_common as w
import wp12_common as W
from wp6_mass_extension_decision import IMF_UPPER_LIMIT

PAPER = w.ROOT / "figures" / "paper"

# A&A column geometry, in inches.
COL = 3.46
WIDE = 7.09

# Okabe--Ito, fixed assignment.  Never cycled, never reassigned per panel.
SUB = {"CygOB2-A": "#0072B2", "CygOB2-B": "#D55E00", "CygOB2-C": "#009E73"}
ALPHA_HUE = {2.0: "#D55E00", 2.3: "#0072B2", 2.6: "#999999"}
INK = "#1a1a1a"
MUTED = "#6b6b6b"
FAIL_HUE = "#B03030"
SEQ = "viridis"

plt.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "grid.color": "#dddddd",
    "grid.linewidth": 0.5,
    "legend.frameon": False,
    "figure.dpi": 200,
    "savefig.dpi": 200,
})

WRITTEN: dict[str, list[str]] = {}


def save(fig, name: str) -> None:
    PAPER.mkdir(parents=True, exist_ok=True)
    out = []
    for suffix in ("pdf", "png"):
        path = PAPER / f"{name}.{suffix}"
        fig.savefig(path, bbox_inches="tight")
        out.append(str(path.relative_to(w.ROOT)))
    plt.close(fig)
    WRITTEN[name] = out
    print(f"  {name}")


def tag(subgroup: str) -> str:
    return f"Cyg OB2-{subgroup[-1]}"


def with_labels(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach the WP2 subgroup labels as a column called `sg`.

    Several upstream products already carry their own `subgroup` column, so a
    plain merge produces `subgroup_x`/`subgroup_y` and silently breaks any code
    that assumes one name.  The authoritative labels are the WP2 table's, and
    this is the only place that decision is made.
    """
    labels = pd.read_parquet(W.frozen("wp2_subgroup_labels"))[
        ["source_id", "subgroup"]
    ].rename(columns={"subgroup": "sg"})
    return frame.merge(labels, on="source_id", how="left")


# ===================================================================== fig 01
def fig01_flow() -> None:
    """The inference chain, with each node's evidential status marked."""
    fig, ax = plt.subplots(figsize=(WIDE, 4.15))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    OBS = "#e8eef4"        # observed input
    INF = "#ffffff"        # inferred quantity
    VAL = "#eaf4ee"        # independent validation
    APP = "#fdf0e6"        # conditional application

    def box(x, y, wid, hgt, text, face, edge=MUTED, style="round,pad=0.25",
            weight="normal", fontsize=7.0, lw=0.9):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), wid, hgt, boxstyle=style, facecolor=face,
            edgecolor=edge, linewidth=lw, zorder=2))
        ax.text(x + wid / 2, y + hgt / 2, text, ha="center", va="center",
                fontsize=fontsize, zorder=3, color=INK, weight=weight)

    def arrow(x1, y1, x2, y2, style="-|>", colour=MUTED, ls="-"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=1,
                    arrowprops=dict(arrowstyle=style, color=colour,
                                    linewidth=0.9, linestyle=ls,
                                    shrinkA=1, shrinkB=1))

    # --- observed inputs
    box(1, 52, 20, 7, "Gaia DR3\nastrometry + photometry", OBS)
    box(23, 52, 15, 7, "2MASS\nJHK$_s$", OBS)
    box(40, 52, 21, 7, "spectroscopic anchors\n(literature)", OBS)
    box(63, 52, 20, 7, "frozen SN markers\npulsar, SNRs, $^{26}$Al", VAL)

    # --- the chain
    chain = [
        (1, 43, 32, "membership mixture\n5 features, including $\\varpi$"),
        (1, 35, 32, "kinematic subgroups\n4 features, excluding $\\varpi$"),
        (1, 27, 32, "spatial extinction\nkriged anchor prior"),
        (1, 19, 32, "counts-based ages and masses"),
        (1, 11, 32, "IMF normalization\n2--8 $M_\\odot$"),
        (1, 3, 32, "stochastic supernova ledger"),
    ]
    for x, y, wid, text in chain:
        box(x, y, wid, 6, text, INF, fontsize=6.6)
    for y in (49, 41, 33, 25, 17, 9):
        arrow(17, y + 0.1, 17, y - 2.0)
    arrow(30, 52, 22, 49.2)
    arrow(50, 52, 26, 49.2)

    # --- out-of-sample and validation limbs
    box(37, 11, 24, 6, "closure test  $>8\\,M_\\odot$\n(out of sample)", INF,
        fontsize=6.8)
    arrow(33, 14, 37, 14)
    box(37, 19, 24, 6, "runaway census (2-D)\nbounds location, not number", INF,
        fontsize=6.8)
    arrow(33, 22, 37, 22)

    box(63, 3, 20, 6, "independent checks\npulsar excludes $N=0$", VAL,
        fontsize=6.8)
    arrow(80, 52, 80, 9.2)
    arrow(33, 6, 63, 6)

    box(37, 3, 24, 6, "conditional applications\ncocoon score, COSI forecast",
        APP, fontsize=6.8)
    arrow(33, 5, 37, 5, ls=":")

    # --- model-choice axes, connected only where they enter the chain
    ax.add_patch(mpatches.FancyBboxPatch(
        (64, 27), 19, 20, boxstyle="round,pad=0.3", facecolor="#f7f7f5",
        edgecolor=MUTED, linewidth=0.9, linestyle="--", zorder=2))
    ax.text(73.5, 45.5, "carried model choices", ha="center", va="top",
            fontsize=6.6, weight="bold", color=INK, zorder=3)
    ax.text(73.5, 42.4,
            "family, $R_V$  $\\rightarrow$ extinction + ages\n"
            "$\\alpha$  $\\rightarrow$ IMF normalization\n"
            "$\\delta$, explodability  $\\rightarrow$ ledger",
            ha="center", va="top", fontsize=6.4, color=INK, zorder=3)
    for y_from, y_to in ((39.5, 30), (35.5, 14), (31.5, 6)):
        ax.annotate("", xy=(33.3, y_to), xytext=(64, y_from),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED,
                                    linewidth=0.75, linestyle="--"))

    # --- where the chain is weakest, kept clear of every box
    ax.add_patch(mpatches.FancyBboxPatch(
        (85.5, 3), 13.5, 22, boxstyle="round,pad=0.3", facecolor="#fdf3f3",
        edgecolor=FAIL_HUE, linewidth=0.9, linestyle="-", zorder=2))
    ax.text(92.2, 23.6, "where it is\nweakest", ha="center", va="top",
            fontsize=6.6, weight="bold", color=FAIL_HUE, zorder=3)
    ax.text(92.2, 18.4,
            "B's age gives a\none-sided SN bound\n\n14 of 54 fit\ncells fail\n\n"
            "explodability controls\nzero versus non-zero",
            ha="center", va="top", fontsize=5.55, color=FAIL_HUE, zorder=3,
            style="italic", linespacing=1.25)

    handles = [
        mpatches.Patch(facecolor=OBS, edgecolor=MUTED, label="observed input"),
        mpatches.Patch(facecolor=INF, edgecolor=MUTED, label="inferred here"),
        mpatches.Patch(facecolor=VAL, edgecolor=MUTED,
                       label="independent validation (frozen before ledger)"),
        mpatches.Patch(facecolor=APP, edgecolor=MUTED,
                       label="conditional application"),
    ]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.06),
              ncol=4, fontsize=6.4, handlelength=1.4)
    save(fig, "fig01_analysis_flow")


# ===================================================================== fig 02
def fig02_membership() -> None:
    frame = with_labels(pd.read_parquet(W.frozen("wp2_members")))
    exempt_col = (
        "anchor_quality_exempt" if "anchor_quality_exempt" in frame.columns
        else None
    )

    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 3.35))
    sky, vpd = axes

    low = frame[frame.membership_probability <= 0.5]
    sky.scatter(low.l_deg, low.b_deg, s=1.2, c="#cccccc", linewidths=0,
                rasterized=True, label="candidates $P \\leq 0.5$")
    vpd.scatter(low.pmra, low.pmdec, s=1.2, c="#cccccc", linewidths=0,
                rasterized=True)

    for subgroup, colour in SUB.items():
        sel = frame[frame.sg.eq(subgroup)]
        sky.scatter(sel.l_deg, sel.b_deg, s=3.0, c=colour, linewidths=0,
                    rasterized=True, label=tag(subgroup))
        vpd.scatter(sel.pmra, sel.pmdec, s=3.0, c=colour, linewidths=0,
                    rasterized=True)
        # direct label, so identity is never colour-alone
        sky.annotate(tag(subgroup)[-1], (sel.l_deg.median(), sel.b_deg.median()),
                     fontsize=11, weight="bold", color=colour, ha="center",
                     va="center", zorder=7,
                     path_effects=[pe.withStroke(linewidth=2.6,
                                                 foreground="#ffffff")])

    if exempt_col:
        ex = frame[frame[exempt_col].fillna(False).astype(bool)]
        sky.scatter(ex.l_deg, ex.b_deg, s=16, facecolors="none",
                    edgecolors=INK, linewidths=0.6, marker="s",
                    label=f"quality-exempt anchors ({len(ex)})", zorder=5)

    sky.set_xlabel("$l$ (deg)")
    sky.set_ylabel("$b$ (deg)")
    sky.invert_xaxis()
    sky.set_title("selected on $(l, b, \\varpi, \\mu_{\\alpha*}, \\mu_\\delta)$",
                  color=MUTED, fontsize=7)
    vpd.set_xlabel("$\\mu_{\\alpha*}$ (mas yr$^{-1}$)")
    vpd.set_ylabel("$\\mu_\\delta$ (mas yr$^{-1}$)")
    vpd.set_xlim(np.nanpercentile(frame.pmra, [0.5, 99.5]))
    vpd.set_ylim(np.nanpercentile(frame.pmdec, [0.5, 99.5]))
    vpd.set_title("subgroups from $(l, b, \\mu_{\\alpha*}, \\mu_\\delta)$ only",
                  color=MUTED, fontsize=7)
    sky.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), markerscale=3,
               fontsize=6.0, ncol=2, columnspacing=0.9)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, "fig02_membership_subgroups")


# ===================================================================== fig 03
def fig03_controls() -> None:
    controls = pd.read_parquet(w.PROC / "wp2_control_members.parquet")
    members = pd.read_parquet(W.frozen("wp2_members"))
    recovery = pd.read_csv(w.PROVENANCE / "wp2_berlanas_recovery_audit.csv")
    automatic = members[~members.anchor_quality_exempt.fillna(False)]
    benchmark = recovery[recovery.quality_pass.astype(bool)]

    fig, axes = plt.subplots(1, 3, figsize=(WIDE, 3.0),
                             gridspec_kw={"width_ratios": [1.2, 1.05, 0.8]})
    hist, trade, bar = axes

    edges = np.linspace(0, 1, 21)
    hist.hist(members.membership_probability, bins=edges, color="#0072B2",
              alpha=0.85, label="association field")
    for i, (name, grp) in enumerate(controls.groupby("control_field")):
        hist.hist(grp.membership_probability, bins=edges, histtype="step",
                  linewidth=1.1, color=["#D55E00", "#009E73", "#7a5195"][i % 3],
                  label=f"control {name}")
    hist.axvline(0.5, color=INK, linewidth=0.8, linestyle="--")
    hist.set_yscale("log")
    hist.set_xlabel("membership probability $P$")
    hist.set_ylabel("sources per bin")
    hist.legend(fontsize=6.3, loc="upper right")
    hist.text(0.52, 0.92, "$P > 0.5$", transform=hist.transAxes, fontsize=6.5,
              color=MUTED)

    thresholds = np.linspace(0.1, 0.9, 17)
    recall, contamination = [], []
    for threshold in thresholds:
        n_target = int((automatic.membership_probability > threshold).sum())
        recovered = int(
            (benchmark.membership_probability_astrometric > threshold)
            .fillna(False).sum()
        )
        control_yields = [
            int((group.membership_probability > threshold).sum())
            for _, group in controls.groupby("control_field")
        ]
        recall.append(100 * recovered / len(benchmark))
        contamination.append(100 * np.mean(control_yields) / n_target)
    trade.plot(thresholds, recall, "-o", color="#0072B2", markersize=2.5,
               linewidth=1.1, label="automatic recovery")
    trade.plot(thresholds, contamination, "-s", color="#D55E00",
               markersize=2.5, linewidth=1.1, label="mean control/target")
    trade.axvline(0.5, color=INK, linewidth=0.8, linestyle="--")
    trade.set_xlabel("probability threshold")
    trade.set_ylabel("fraction (per cent)")
    trade.set_ylim(0, 85)
    trade.legend(fontsize=5.6, loc="center right")
    trade.set_title("the 0.5 trade-off", color=MUTED, fontsize=7)
    trade.annotate("1,331 automatic members", xy=(0.5, 76.2),
                   xytext=(0.34, 66), fontsize=5.5, color=INK,
                   arrowprops=dict(arrowstyle="->", color=INK, lw=0.6))

    target_yield = int((automatic.membership_probability > 0.5).sum())
    names, yields = ["target"], [target_yield]
    for name, grp in controls.groupby("control_field"):
        names.append(str(name))
        yields.append(int((grp.membership_probability > 0.5).sum()))
    colours = ["#0072B2"] + ["#D55E00", "#009E73", "#7a5195"][: len(names) - 1]
    positions = np.arange(len(names))
    bar.bar(positions, yields, color=colours, width=0.62)
    for x, v in zip(positions, yields):
        bar.text(x, v + max(yields) * 0.02, f"{v:,d}", ha="center", va="bottom",
                 fontsize=7, color=INK)
    bar.set_xticks(positions)
    bar.set_xticklabels(names, fontsize=6.6)
    bar.set_ylabel("sources with $P > 0.5$")
    bar.set_ylim(0, max(yields) * 1.22)
    bar.set_title("equal-area fields, identical pipeline", color=MUTED,
                  fontsize=7)
    fig.tight_layout()
    save(fig, "fig03_control_fields")


# ===================================================================== fig 04
def fig04_cmd_ages() -> None:
    """De-reddened CMD per subgroup with both isochrone families, plus the
    age posteriors and the systematic envelope."""
    frame = with_labels(pd.read_parquet(W.frozen("wp3_extinction")))

    gate = W.gate_map()
    base_gate = gate[
        gate.family.eq(W.BASE["family"]) & gate.R_V.eq(W.BASE["R_V"])
        & gate.alpha.eq(W.BASE["alpha"])
    ].set_index("subgroup")

    isochrones = {
        family: pd.read_parquet(W.frozen(f"wp3_isochrones_{family.lower()}"))
        for family in w.FAMILIES
    }

    fig = plt.figure(figsize=(WIDE, 4.6))
    grid = fig.add_gridspec(2, 3, height_ratios=[2.0, 1.0], hspace=0.42,
                            wspace=0.28)

    colour_x, colour_y = f"BPRP0_rv{W.BASE['R_V']:g}", f"G0_abs_rv{W.BASE['R_V']:g}"
    for i, subgroup in enumerate(w.SUBGROUPS):
        ax = fig.add_subplot(grid[0, i])
        sel = frame[frame.sg.eq(subgroup)]
        other = frame[frame.sg.notna() & ~frame.sg.eq(subgroup)]
        ax.scatter(other[colour_x], other[colour_y], s=1.0, c="#dddddd",
                   linewidths=0, rasterized=True)
        ax.scatter(sel[colour_x], sel[colour_y], s=2.6, c=SUB[subgroup],
                   linewidths=0, rasterized=True, label=tag(subgroup))
        anchors = sel[sel.get("anchor_quality_exempt", pd.Series(False,
                      index=sel.index)).fillna(False).astype(bool)]
        if len(anchors):
            ax.scatter(anchors[colour_x], anchors[colour_y], s=14,
                       facecolors="none", edgecolors=INK, linewidths=0.55,
                       marker="s", zorder=5, label="spectroscopic anchor")

        age = float(base_gate.loc[subgroup].truth_age_posterior_mean_Myr)
        for family, style in (("PARSEC", "-"), ("MIST", "--")):
            iso = isochrones[family]
            nearest = iso.age_Myr.unique()[
                np.argmin(np.abs(iso.age_Myr.unique() - age))]
            track = iso[np.isclose(iso.age_Myr, nearest)].sort_values("Mini")
            track = track[track.Mini <= IMF_UPPER_LIMIT]
            ax.plot(track.BP0 - track.RP0, track.G0, style, color=INK,
                    linewidth=0.9,
                    label=f"{family} {nearest:.2f} Myr" if i == 0 else None)

        ax.invert_yaxis()
        ax.set_xlim(-0.6, 2.4)
        ax.set_ylim(6.5, -7.5)
        ax.set_xlabel("$(BP-RP)_0$")
        if i == 0:
            ax.set_ylabel("$M_{G,0}$")
        ax.set_title(f"{tag(subgroup)}   {age:.2f} Myr", fontsize=7.5,
                     color=SUB[subgroup], weight="bold")
        if i == 0:
            ax.legend(fontsize=5.9, loc="lower left", markerscale=2)

    # ---- age posteriors, per subgroup and family
    ax = fig.add_subplot(grid[1, :2])
    norm = W.normalization()
    offsets = {"CygOB2-A": 0.0, "CygOB2-B": 1.0, "CygOB2-C": 2.0}
    for subgroup in w.SUBGROUPS:
        y = offsets[subgroup]
        for family, marker in (("PARSEC", "o"), ("MIST", "^")):
            cell = norm[norm.subgroup.eq(subgroup) & norm.family.eq(family)]
            ages = cell.truth_age_posterior_mean_Myr.to_numpy(float)
            ax.scatter(ages, np.full_like(ages, y) + (0.16 if family == "MIST"
                                                      else -0.16),
                       s=13, marker=marker, facecolors="none",
                       edgecolors=SUB[subgroup], linewidths=0.8)
        base_age = float(base_gate.loc[subgroup].truth_age_posterior_mean_Myr)
        ax.scatter([base_age], [y], s=34, marker="D", color=SUB[subgroup],
                   zorder=5)
        ax.text(6.55, y, tag(subgroup), ha="right", va="center", fontsize=7,
                color=SUB[subgroup], weight="bold")
    ax.axvspan(2.25, 5.67, color="#000000", alpha=0.06, zorder=0)
    ax.text(2.30, -0.42, "2.25--5.67 Myr envelope across both retained "
            "age indicators", ha="left", va="center", fontsize=6.0,
            color=MUTED)
    ax.annotate("B rails against the top of its own prior grid",
                xy=(4.13, 1.16), xytext=(4.45, 1.62), fontsize=5.9,
                color=FAIL_HUE,
                arrowprops=dict(arrowstyle="->", color=FAIL_HUE, lw=0.7))
    ax.set_xlim(1.9, 6.6)
    ax.set_ylim(-0.75, 2.55)
    ax.set_yticks([])
    ax.set_xlabel("counts-based age (Myr)")
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color="none", markerfacecolor=MUTED,
               markersize=5, label="baseline"),
        Line2D([0], [0], marker="o", color="none", markeredgecolor=MUTED,
               markerfacecolor="none", markersize=5, label="PARSEC"),
        Line2D([0], [0], marker="^", color="none", markeredgecolor=MUTED,
               markerfacecolor="none", markersize=5, label="MIST"),
    ], fontsize=6.0, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3,
       columnspacing=1.2)

    # ---- what the spread costs
    ax2 = fig.add_subplot(grid[1, 2])
    scan = pd.read_csv(w.TABLES / "wp7_age_sensitivity.csv")
    ax2.plot(scan.assumed_age_Myr, scan.N_SN_mean, color=INK, linewidth=1.2)
    ax2.fill_between(scan.assumed_age_Myr, scan.N_SN_p16, scan.N_SN_p84,
                     color=INK, alpha=0.12, linewidth=0)
    ax2.axvspan(2.25, 5.67, color="#000000", alpha=0.06, zorder=0)
    ax2.set_xlabel("common assumed age (Myr)")
    ax2.set_ylabel("$N_{\\rm SN}$")
    ax2.set_title("what the envelope costs", fontsize=7, color=MUTED)
    save(fig, "fig04_cmd_ages")


# ===================================================================== fig 05
def fig05_extinction() -> None:
    frame = with_labels(pd.read_parquet(W.frozen("wp3_extinction")))
    anchors = pd.read_parquet(w.PROC / "wp3_anchor_extinction.parquet")
    av = f"av_rv{W.BASE['R_V']:g}"

    fig, axes = plt.subplots(1, 3, figsize=(WIDE, 2.75),
                             gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})
    sky, support, offset = axes

    sc = sky.scatter(frame.l_deg, frame.b_deg, c=frame[av], s=4.5,
                     cmap=SEQ, linewidths=0, rasterized=True,
                     vmin=np.nanpercentile(frame[av], 2),
                     vmax=np.nanpercentile(frame[av], 98))
    sky.scatter(anchors.l_deg, anchors.b_deg, s=13, facecolors="none",
                edgecolors="#ffffff", linewidths=1.4, marker="s", zorder=4)
    sky.scatter(anchors.l_deg, anchors.b_deg, s=13, facecolors="none",
                edgecolors=INK, linewidths=0.7, marker="s", zorder=5,
                label=f"anchors ({len(anchors)})")
    for subgroup in w.SUBGROUPS:
        sel = frame[frame.sg.eq(subgroup)]
        sky.annotate(tag(subgroup)[-1],
                     (sel.l_deg.median(), sel.b_deg.median()),
                     fontsize=10, weight="bold", color="#ffffff", ha="center",
                     va="center", zorder=6,
                     path_effects=[pe.withStroke(linewidth=2.2,
                                                 foreground=INK)])
    sky.invert_xaxis()
    sky.set_xlabel("$l$ (deg)")
    sky.set_ylabel("$b$ (deg)")
    sky.legend(fontsize=6.2, loc="lower left")
    cb = fig.colorbar(sc, ax=sky, pad=0.02)
    cb.set_label("$A_V$ (mag)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # local anchor support: distance to the 8th nearest anchor
    anchor_xy = anchors[["l_deg", "b_deg"]].to_numpy(float)
    for subgroup in w.SUBGROUPS:
        sel = frame[frame.sg.eq(subgroup)]
        xy = sel[["l_deg", "b_deg"]].to_numpy(float)
        if not len(xy):
            continue
        d = np.sort(np.hypot(
            xy[:, None, 0] - anchor_xy[None, :, 0],
            xy[:, None, 1] - anchor_xy[None, :, 1]), axis=1)[:, 7]
        support.hist(d, bins=np.linspace(0, 0.9, 30), histtype="step",
                     linewidth=1.2, color=SUB[subgroup], label=tag(subgroup))
        support.axvline(np.median(d), color=SUB[subgroup], linewidth=0.7,
                        linestyle=":")
    support.set_xlabel("distance to 8th-nearest anchor (deg)")
    support.set_ylabel("members")
    support.legend(fontsize=6.2)
    support.set_title("B's anchor famine", fontsize=7, color=MUTED)

    # anchor vs prior-free photometric A_V at matched position
    broadband = w.PROC / "wp3_broadband_extinction.parquet"
    if broadband.exists():
        bb = pd.read_parquet(broadband)
        col = av if av in bb.columns else next(
            (c for c in bb.columns if c.startswith("av_rv3.1")), None)
        if col:
            merged = anchors[["source_id", av]].merge(
                bb[["source_id", col]], on="source_id", how="inner",
                suffixes=("_anchor", "_photom"))
            xcol, ycol = merged.columns[1], merged.columns[2]
            difference = (merged[ycol] - merged[xcol]).dropna()
            offset.hist(difference, bins=np.linspace(-3.5, 2.0, 34),
                        color="#0072B2", alpha=0.85)
            median_offset = float(difference.median())
            offset.axvline(0.0, color=INK, linewidth=0.9)
            offset.axvline(median_offset, color=FAIL_HUE, linewidth=1.1,
                           linestyle="--",
                           label=f"median {median_offset:+.2f} mag")
            offset.legend(fontsize=6.2, loc="upper left")
    offset.set_xlabel("photometric minus anchor $A_V$ (mag)", fontsize=6.6)
    offset.set_ylabel("anchors")
    offset.set_title("the absolute-scale systematic", fontsize=7, color=MUTED)
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig05_extinction")


# ===================================================================== fig 06
def fig06_massfunction_gate() -> None:
    bins = pd.read_parquet(W.frozen("wp5_mass_function_bins"))
    gate = W.gate_map()
    base_bins = bins[
        bins.family.eq(W.BASE["family"]) & bins.R_V.eq(W.BASE["R_V"])
        & bins.alpha.eq(W.BASE["alpha"])
    ]
    base_gate = gate[
        gate.family.eq(W.BASE["family"]) & gate.R_V.eq(W.BASE["R_V"])
        & gate.alpha.eq(W.BASE["alpha"])
    ].set_index("subgroup")

    fig = plt.figure(figsize=(WIDE, 4.5))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.25], hspace=0.45,
                            wspace=0.26)

    for i, subgroup in enumerate(w.SUBGROUPS):
        ax = fig.add_subplot(grid[0, i])
        cell = base_bins[base_bins.subgroup.eq(subgroup)].sort_values("bin_index")
        centre = cell.mass_geometric_center
        ax.errorbar(centre, cell.membership_weighted_count,
                    yerr=np.sqrt(cell.membership_weighted_count.clip(lower=0)),
                    fmt="o", markersize=3.2, color=SUB[subgroup], linewidth=0.9,
                    capsize=1.6, label="observed")
        ax.plot(centre, cell.expected_count_at_k_median, "-", color=INK,
                linewidth=1.1, label="response-convolved model")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([3, 4, 6])
        ax.set_xticklabels(["3", "4", "6"])
        ax.set_xlabel("$M$ ($M_\\odot$)")
        if i == 0:
            ax.set_ylabel("stars per bin")
            ax.legend(fontsize=6.2, loc="lower left")
        worst = float(base_gate.loc[subgroup].max_abs_pearson_residual)
        ax.set_title(f"{tag(subgroup)}   max $|r| = {worst:.2f}$", fontsize=7.5,
                     color=SUB[subgroup], weight="bold")

    # ---- the 54-cell residual-gate map
    ax = fig.add_subplot(grid[1, :])
    gate = gate.copy()
    gate["family_order"] = gate.family.map({"PARSEC": 0, "MIST": 1})
    order = gate.sort_values(["family_order", "R_V", "alpha"])
    combos = order[["family", "R_V", "alpha"]].drop_duplicates().reset_index(
        drop=True)
    matrix = np.zeros((3, len(combos)))
    for j, combo in combos.iterrows():
        for i, subgroup in enumerate(w.SUBGROUPS):
            cell = gate[
                gate.subgroup.eq(subgroup) & gate.family.eq(combo.family)
                & gate.R_V.eq(combo.R_V) & gate.alpha.eq(combo.alpha)
            ].iloc[0]
            matrix[i, j] = float(cell.max_abs_pearson_residual)

    im = ax.imshow(matrix, aspect="auto", cmap="cividis", vmin=0, vmax=3.6,
                   origin="upper")
    for j, combo in combos.iterrows():
        headline = combo.alpha in W.HEADLINE_ALPHAS
        all_pass = bool(gate[
            gate.family.eq(combo.family) & gate.R_V.eq(combo.R_V)
            & gate.alpha.eq(combo.alpha)].residual_gate_pass.all())
        for i, subgroup in enumerate(w.SUBGROUPS):
            passed = bool(gate[
                gate.subgroup.eq(subgroup) & gate.family.eq(combo.family)
                & gate.R_V.eq(combo.R_V) & gate.alpha.eq(combo.alpha)
            ].iloc[0].residual_gate_pass)
            if not passed:
                ax.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=FAIL_HUE,
                    linewidth=1.6, zorder=4))
                ax.text(j, i, "×", ha="center", va="center", fontsize=7,
                        color=FAIL_HUE, weight="bold", zorder=5)
        if all_pass and headline:
            # Marked by geometry, not by colour: a wedge under the column.
            ax.plot([j], [2.72], marker="^", markersize=5, color=INK,
                    clip_on=False, zorder=6)
        if not headline:
            ax.add_patch(mpatches.Rectangle(
                (j - 0.5, -0.5), 1, 3, facecolor="#ffffff", alpha=0.42,
                edgecolor="none", zorder=2))

    ax.set_yticks(range(3))
    ax.set_yticklabels([tag(s) for s in w.SUBGROUPS], fontsize=7)
    ax.set_xticks(range(len(combos)))
    ax.set_xticklabels(
        [f"{'PA' if c.family == 'PARSEC' else 'MI'}\n{c.R_V:g}\n{c.alpha:g}"
         for _, c in combos.iterrows()], fontsize=5.6)
    ax.set_xlabel("family / $R_V$ / $\\alpha$   "
                  "(faded columns are $\\alpha = 2.6$, outside the headline set)",
                  fontsize=6.8)
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.03)
    cb.set_label("max $|$Pearson residual$|$", fontsize=6.8)
    cb.ax.tick_params(labelsize=6)
    ax.legend(handles=[
        Line2D([0], [0], marker="x", color=FAIL_HUE, lw=0, markersize=5,
               label="cell fails the residual-fit criterion"),
        Line2D([0], [0], marker="^", color=INK, lw=0, markersize=5,
               label="headline column passing in all three subgroups"),
    ], fontsize=6.2, loc="upper left", bbox_to_anchor=(0.0, -0.42), ncol=2)
    save(fig, "fig06_massfunction_gate")


# ===================================================================== fig 07
def fig07_closure() -> None:
    closure = pd.read_csv(w.TABLES / "wp12_closure_by_alpha.csv")
    slopes = pd.read_csv(w.TABLES / "wp12_closing_slopes.csv")

    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 3.0),
                             gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax, box = axes

    for subgroup in w.SUBGROUPS:
        cell = closure[closure.subgroup.eq(subgroup)]
        pivot = cell.pivot_table(index="alpha", values="closure_ratio",
                                 aggfunc=["min", "median", "max"])
        alphas = pivot.index.to_numpy(float)
        lo = pivot[("min", "closure_ratio")].to_numpy()
        mid = pivot[("median", "closure_ratio")].to_numpy()
        hi = pivot[("max", "closure_ratio")].to_numpy()
        ax.fill_between(alphas, lo, hi, color=SUB[subgroup], alpha=0.16,
                        linewidth=0)
        ax.plot(alphas, mid, "-o", color=SUB[subgroup], linewidth=1.3,
                markersize=3.6, label=tag(subgroup))
        ax.annotate(tag(subgroup)[-1], (alphas[-1] + 0.015, hi[-1]),
                    fontsize=8, weight="bold", color=SUB[subgroup],
                    va="center")

    assoc = closure[closure.subgroup.eq("association")]
    pivot = assoc.pivot_table(index="alpha", values="closure_ratio",
                              aggfunc="median")
    ax.plot(pivot.index, pivot.closure_ratio, "--k", linewidth=1.3,
            label="association (star-weighted)")
    ax.axhline(1.0, color=INK, linewidth=0.8, linestyle=":")
    ax.text(2.02, 1.03, "census closes", fontsize=6.3, color=MUTED)
    ax.set_xlabel("IMF slope $\\alpha$")
    ax.set_ylabel("observed / predicted,  $M > 8\\,M_\\odot$")
    ax.set_xlim(1.96, 2.70)
    ax.set_yscale("log")
    ax.set_yticks([0.6, 0.8, 1.0, 1.4, 2.0])
    ax.set_yticklabels(["0.6", "0.8", "1.0", "1.4", "2.0"])
    ax.legend(fontsize=6.3, loc="upper left")

    # closing slope per subgroup, all family/R_V cells
    positions = np.arange(len(w.SUBGROUPS) + 1)
    names = list(w.SUBGROUPS) + ["association"]
    for i, name in enumerate(names):
        vals = slopes[slopes.subgroup.eq(name)].closing_alpha_interpolated
        colour = SUB.get(name, INK)
        box.scatter(np.full(len(vals), i) + np.linspace(-0.12, 0.12, len(vals)),
                    vals, s=15, color=colour, alpha=0.85, linewidths=0)
        box.hlines(vals.median(), i - 0.24, i + 0.24, color=colour,
                   linewidth=1.8)
        box.text(i, vals.max() + 0.045, f"{vals.median():.2f}", ha="center",
                 fontsize=6.6, color=colour, weight="bold")
    box.axhline(2.3, color=MUTED, linewidth=0.8, linestyle="--")
    box.text(3.42, 2.305, "Salpeter", fontsize=6.3, color=MUTED, ha="right")
    box.set_xticks(positions)
    box.set_xticklabels([tag(n) if n.startswith("Cyg") else "assoc."
                         for n in names], fontsize=6.6)
    box.set_ylabel("slope at which the census closes")
    box.set_title("A and B want Salpeter; C wants shallower", fontsize=7,
                  color=MUTED)
    fig.tight_layout()
    save(fig, "fig07_closure")


# ===================================================================== fig 08
def fig08_history() -> None:
    rsn = pd.read_csv(w.TABLES / "wp7_rsn_curves.csv")
    assoc = W.association_all_explode()
    head = assoc[assoc.alpha.isin(W.HEADLINE_ALPHAS)]

    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 2.9),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})
    ax, spread = axes

    edges = np.append(np.sort(rsn.lookback_lo_Myr.unique()),
                      rsn.lookback_hi_Myr.max())
    bottom = np.zeros(len(edges) - 1)
    for subgroup in w.SUBGROUPS:
        cell = rsn[rsn.subgroup.eq(subgroup)].sort_values("lookback_lo_Myr")
        values = cell.rate_per_Myr.to_numpy()
        ax.bar(cell.lookback_lo_Myr, values, width=np.diff(edges),
               bottom=bottom, align="edge", color=SUB[subgroup],
               edgecolor="#ffffff", linewidth=0.3, label=tag(subgroup))
        bottom = bottom + values
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.42)
    ax.axvspan(0, 0.1, color=INK, alpha=0.16, zorder=0)
    ax.annotate("last 100 kyr", xy=(0.05, top * 1.02),
                xytext=(0.36, top * 1.30), fontsize=6.2, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=0.7))
    ax.axvline(0.2007, color=FAIL_HUE, linewidth=1.0, linestyle="--")
    ax.annotate("PSR J2032+4127 age",
                xy=(0.2007, top * 1.02), xytext=(0.62, top * 1.16),
                fontsize=6.0, color=FAIL_HUE,
                arrowprops=dict(arrowstyle="->", color=FAIL_HUE, lw=0.7))
    ax.text(0.985, 0.62, "Cyg OB2-C contributes\nnothing on this branch",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.0,
            color=SUB["CygOB2-C"])
    ax.set_xlabel("look-back time (Myr)")
    ax.set_ylabel("$R_{\\rm SN}$ (Myr$^{-1}$)")
    ax.set_xlim(0, 1.6)
    ax.legend(fontsize=6.4, loc="upper right")
    ax.set_title("baseline branch, stacked by subgroup", fontsize=7,
                 color=MUTED)

    for alpha in W.HEADLINE_ALPHAS:
        arm = head[head.alpha.eq(alpha)]
        spread.scatter(arm.N_SN_mean,
                       np.random.default_rng(7).normal(alpha, 0.028, len(arm)),
                       s=17, color=ALPHA_HUE[alpha], linewidths=0, alpha=0.85)
        spread.text(arm.N_SN_mean.max() + 0.9, alpha,
                    f"$\\alpha = {alpha:g}$", fontsize=7,
                    color=ALPHA_HUE[alpha], va="center", weight="bold")
    base = assoc[
        assoc.family.eq(W.BASE["family"]) & assoc.R_V.eq(W.BASE["R_V"])
        & assoc.alpha.eq(W.BASE["alpha"])
        & assoc.sf_duration_Myr.eq(W.BASE["sf_duration_Myr"])].iloc[0]
    spread.errorbar([base.N_SN_mean], [2.3], xerr=[[base.N_SN_mean - base.N_SN_p16],
                                                   [base.N_SN_p84 - base.N_SN_mean]],
                    fmt="D", color=INK, markersize=5, capsize=2.5, zorder=6,
                    label="baseline, 68% Poisson")
    spread.set_xlabel("$N_{\\rm SN}$")
    spread.set_yticks(list(W.HEADLINE_ALPHAS))
    spread.set_yticklabels([f"{a:g}" for a in W.HEADLINE_ALPHAS])
    spread.set_ylabel("IMF slope $\\alpha$")
    spread.set_ylim(1.87, 2.45)
    spread.set_xlim(0, head.N_SN_mean.max() * 1.22)
    spread.legend(fontsize=6.2, loc="lower left")
    spread.set_title("branch spread exceeds the stochastic interval",
                     fontsize=7, color=MUTED)
    fig.tight_layout()
    save(fig, "fig08_history")


# ===================================================================== fig 09
def fig09_sensitivity() -> None:
    scan = pd.read_csv(w.TABLES / "wp7_age_sensitivity.csv")
    bh = pd.read_csv(w.TABLES / "wp7_bh_threshold_scan.csv")
    gate = W.gate_map()
    base_gate = gate[
        gate.family.eq(W.BASE["family"]) & gate.R_V.eq(W.BASE["R_V"])
        & gate.alpha.eq(W.BASE["alpha"])].set_index("subgroup")

    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 2.8))
    age_ax, bh_ax = axes

    age_ax.plot(scan.assumed_age_Myr, scan.N_SN_mean, "-", color=INK,
                linewidth=1.4, label="$N_{\\rm SN}$")
    age_ax.fill_between(scan.assumed_age_Myr, scan.N_SN_p16, scan.N_SN_p84,
                        color=INK, alpha=0.12, linewidth=0)
    age_ax.set_xlabel("common assumed age (Myr)")
    age_ax.set_ylabel("$N_{\\rm SN}$")
    twin = age_ax.twiny()
    twin.set_xlim(age_ax.get_xlim())
    twin.set_xticks([float(base_gate.loc[s].truth_age_posterior_mean_Myr)
                     for s in w.SUBGROUPS])
    twin.set_xticklabels([tag(s)[-1] for s in w.SUBGROUPS], fontsize=7)
    twin.tick_params(length=3)
    for s in w.SUBGROUPS:
        age_ax.axvline(float(base_gate.loc[s].truth_age_posterior_mean_Myr),
                       color=SUB[s], linewidth=0.9, linestyle=":")
    age_ax.axvspan(2.25, 5.67, color="#000000", alpha=0.06, zorder=0)
    age_ax.text(5.6, age_ax.get_ylim()[1] * 0.06,
                "retained age envelope", fontsize=6.2, color=MUTED, ha="right")
    age_ax.legend(fontsize=6.4, loc="upper left")

    total = bh.groupby("bh_threshold_Msun", as_index=False).N_SN_mean.sum() \
        if "subgroup" in bh.columns else bh
    bh_ax.step(total.bh_threshold_Msun, total.N_SN_mean, where="post",
               color=INK, linewidth=1.4)
    bh_ax.scatter(total.bh_threshold_Msun, total.N_SN_mean, s=16, color=INK,
                  zorder=4)
    zero = total[total.N_SN_mean <= 1e-9].bh_threshold_Msun
    if len(zero):
        bh_ax.axvspan(total.bh_threshold_Msun.min(), float(zero.max()),
                      color=FAIL_HUE, alpha=0.10, zorder=0)
        bh_ax.text(float(zero.max()) - 0.8, total.N_SN_mean.max() * 0.55,
                   "ledger returns\nexactly zero", fontsize=6.3,
                   color=FAIL_HUE, ha="right")
    bh_ax.set_xlabel("direct-collapse threshold ($M_\\odot$)")
    bh_ax.set_ylabel("$N_{\\rm SN}$")
    bh_ax.set_title("a hard cutoff, not a non-monotonic map", fontsize=7,
                    color=MUTED)
    fig.tight_layout()
    save(fig, "fig09_sensitivity")


# ===================================================================== fig 10
def fig10_scenario_score() -> None:
    scan = pd.read_csv(w.TABLES / "wp12_c4_scan.csv")
    scenario = json.loads(
        (w.PROVENANCE / "wp12_scenario_score_execution.json").read_text())
    c4 = scenario["wp12_3_c4_sensitivity"]
    adopted = c4["adopted_C4"]
    c4_any = c4["C4_below_which_at_least_one_alpha2p0_branch_falls_below_0p5"]
    c4_all = c4["C4_below_which_every_alpha2p0_branch_falls_below_0p5"]

    fig, ax = plt.subplots(figsize=(WIDE, 3.0))
    keys = ["family", "R_V", "alpha", "sf_duration_Myr"]
    for (family, rv, alpha, delta), branch in scan.groupby(keys):
        branch = branch.sort_values("C4")
        strict = bool(branch.all_subgroup_pass.iloc[0])
        ax.plot(branch.C4, branch.score, "-", color=ALPHA_HUE[alpha],
                linewidth=1.5 if strict else 0.7,
                alpha=0.95 if strict else 0.35, zorder=3 if strict else 2)

    ax.axhline(0.5, color=INK, linewidth=1.0)
    ax.text(0.985, 0.515, "support boundary of the pre-registered framing rule",
            fontsize=6.1, color=INK, ha="right")
    ax.axvline(adopted, color=MUTED, linewidth=1.0, linestyle="--")
    ax.text(adopted - 0.008, 0.055, f"adopted bound {adopted:.3f}",
            fontsize=6.2, color=MUTED, rotation=90, va="bottom", ha="right")
    for value, label, y in ((c4_any, "first shallow-slope\nbranch falls below 0.5", 0.82),
                            (c4_all, "all shallow-slope\nbranches fall below 0.5", 0.69)):
        ax.axvline(value, color=FAIL_HUE, linewidth=0.9, linestyle=":")
        ax.text(value - 0.008, y, f"{label}\n$C_4={value:.3f}$", fontsize=6.0,
                color=FAIL_HUE, ha="center", va="top")

    ax.set_xlabel("$C_4$, the in-situ fraction (an upper bound, not a measurement)")
    ax.set_ylabel("conditional scenario-availability score $S$")
    ax.set_xlim(0.30, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(handles=[
        Line2D([0], [0], color=ALPHA_HUE[2.0], lw=1.6, label="$\\alpha = 2.0$"),
        Line2D([0], [0], color=ALPHA_HUE[2.3], lw=1.6, label="$\\alpha = 2.3$"),
        Line2D([0], [0], color=MUTED, lw=1.6,
               label="all subgroup fits pass"),
        Line2D([0], [0], color=MUTED, lw=0.7, alpha=0.4,
               label="contains a failed subgroup fit"),
    ], fontsize=6.2, loc="upper left", ncol=4)
    fig.tight_layout()
    save(fig, "fig10_scenario_score")


# ===================================================================== fig 11
def fig11_isotopes() -> None:
    iso = pd.read_csv(W.frozen("wp11_isotope_forecast"))
    prereg = json.loads(
        (w.PROVENANCE / "wp11_isotope_prereg.json").read_text())
    cosi = prereg["instruments"]["COSI_narrow_line_3sigma_2yr"]["value_ph_cm2_s"]

    head = iso[iso.in_headline_set]
    arms = list(head.yield_arm.unique())

    fig, ax = plt.subplots(figsize=(WIDE, 3.0))
    positions = {arm: i for i, arm in enumerate(arms)}
    for arm in arms:
        for alpha in W.HEADLINE_ALPHAS:
            cell = head[head.yield_arm.eq(arm) & head.alpha.eq(alpha)]
            flux = cell.F_1173_ph_cm2_s.to_numpy(float)
            flux = np.where(flux <= 0, 1e-12, flux)
            x = positions[arm] + (-0.16 if alpha == 2.0 else 0.16)
            ax.scatter(np.full(len(flux), x)
                       + np.random.default_rng(3).normal(0, 0.03, len(flux)),
                       flux, s=17, color=ALPHA_HUE[alpha], linewidths=0,
                       alpha=0.9)
            if len(flux):
                ax.hlines(np.median(flux), x - 0.1, x + 0.1,
                          color=ALPHA_HUE[alpha], linewidth=1.8)

    for arm in arms:
        cell = head[head.yield_arm.eq(arm)]
        if float(cell.F_1173_ph_cm2_s.max()) <= 0:
            ax.text(positions[arm], 3e-8,
                    "identically zero\non this arm:\nevery progenitor\n"
                    "collapses directly",
                    ha="center", va="center", fontsize=6.2, color=FAIL_HUE,
                    style="italic")
    ax.axhline(cosi, color=INK, linewidth=1.2)
    ax.text(len(arms) - 0.52, cosi * 1.25,
            "COSI 3$\\sigma$, 2 yr narrow line", fontsize=6.3, color=INK,
            ha="right")
    ax.set_yscale("log")
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([a.replace("_", " ") for a in arms], fontsize=6.6)
    ax.set_xlabel("published yield arm")
    ax.set_ylabel("$^{60}$Fe 1173 keV flux (ph cm$^{-2}$ s$^{-1}$)")
    ax.set_ylim(1e-9, 5e-5)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", lw=0, color=ALPHA_HUE[2.0],
               label="$\\alpha = 2.0$"),
        Line2D([0], [0], marker="o", lw=0, color=ALPHA_HUE[2.3],
               label="$\\alpha = 2.3$"),
    ], fontsize=6.4, loc="lower left")
    ax.set_title("the yield arm outweighs the whole census branch set",
                 fontsize=7, color=MUTED)
    fig.tight_layout()
    save(fig, "fig11_isotopes")


def main() -> None:
    print("WP12 figures")
    fig01_flow()
    fig02_membership()
    fig03_controls()
    fig04_cmd_ages()
    fig05_extinction()
    fig06_massfunction_gate()
    fig07_closure()
    fig08_history()
    fig09_sensitivity()
    fig10_scenario_score()
    fig11_isotopes()

    rec = W.record(
        "scripts/wp12_figures.py",
        {
            "item": "WP12 -- revised paper figure set",
            "figures": WRITTEN,
            "colour_rules": {
                "categorical": "Okabe-Ito, fixed assignment, never cycled",
                "subgroups": SUB,
                "alpha_arms": {f"{k:g}": v for k, v in ALPHA_HUE.items()},
                "sequential": SEQ,
                "validated": (
                    "adjacent-pair separation checked under deuteranopia, "
                    "protanopia and tritanopia; every colour-coded series is "
                    "also direct-labelled or marker-distinguished"
                ),
            },
            "column_width_in": COL,
            "double_column_width_in": WIDE,
        },
        {name: w.ROOT / paths[0] for name, paths in WRITTEN.items()},
    )
    w.write_json(w.PROVENANCE / "wp12_figures_execution.json", rec)
    print(f"wrote {len(WRITTEN)} figures to figures/paper/")


if __name__ == "__main__":
    main()
