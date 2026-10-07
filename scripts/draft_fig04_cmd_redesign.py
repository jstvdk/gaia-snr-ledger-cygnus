#!/usr/bin/env python3
"""DRAFT -- a redesign of the fig04 CMD panels.  Not part of the manuscript chain.

The published fig04 draws one black isochrone over single-colour dots.  This
draft puts on the same dots what the analysis actually takes from the CMD:

  * each star's baseline mass class (PARSEC, R_V = 3.1) -- below the
    calibration window, inside it (sets k), or above 8 Msun (the living census
    the closure test compares against);
  * the isochrone as a mass ruler: tick marks at fixed initial masses along the
    main sequence, the end of the main sequence, and the turnoff mass above
    which every star is dead;
  * the upper-MS magnitude window that the retained age fit uses;
  * the MIST isochrone at MIST's own fitted age, main sequence only, because
    the family disagreement lives at the bright end;
  * spectroscopic stars ringed, and a fourth panel for the 61 quality-exempt
    anchors that carry no subgroup label and are absent from the published
    panels.

Inputs are read through wp12_common.frozen() and wp10_inputs.resolve(); nothing
is refitted and no stored number moves.

Outputs: figures/drafts/fig04_cmd_redesign_draft.{pdf,png}
         figures/drafts/fig04_cmd_only_draft.{pdf,png}   (--cmd-only: CMD panels and legend only)

Run:
  PYTHONPATH=scripts python3 scripts/draft_fig04_cmd_redesign.py [--cmd-only]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import wp5_common as w
import wp12_common as W
from wp10_inputs import resolve
from wp12_figures import INK, MUTED, SUB, WIDE, tag, with_labels
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

OUT = w.ROOT / "figures" / "drafts"

# Mass classes.  Okabe--Ito hues not used for the subgroups, so a mass class is
# never mistaken for a subgroup identity.
CLASSES = (
    ("$M < 2\\,M_\\odot$: below the window", 0.0, 2.0, "#c8c8c8"),
    ("$2$–$8\\,M_\\odot$: calibration window, sets $k$", 2.0, 8.0, "#56B4E9"),
    ("$M \\geq 8\\,M_\\odot$: living massive census", 8.0, np.inf, "#E69F00"),
)
TICK_MASSES = (2, 3, 5, 8, 15, 25, 40)
UMS_LIMIT = 1.5          # WP4 upper-MS indicator window, M_G0 <= +1.5
YLIM = (6.2, -9.9)
XLIM = (-0.95, 2.3)


def nearest_track(iso: pd.DataFrame, age: float) -> tuple[float, pd.DataFrame]:
    ages = iso.age_Myr.unique()
    native = float(ages[np.argmin(np.abs(ages - age))])
    track = iso[np.isclose(iso.age_Myr, native) & (iso.Mini <= IMF_UPPER_LIMIT)]
    track = track.sort_values("Mini").reset_index(drop=True)
    return native, track.assign(colour=track.BP0 - track.RP0)


def split_main_sequence(track: pd.DataFrame) -> int:
    """Index of the terminal-age main sequence: the 'hook', i.e. the first local
    minimum of log Teff reached after the upper MS has cooled by more than
    0.05 dex from its hottest point, or the last row if there is none.

    Along the upper MS log Teff first RISES with mass and then falls as the
    most massive stars evolve; the hook is where it turns up again.
    """
    logte = track.logTe.to_numpy()
    hottest, lowest, where = -np.inf, np.inf, None
    for i in np.flatnonzero(track.Mini.to_numpy() > 10.0):
        hottest = max(hottest, logte[i])
        if logte[i] < hottest - 0.05:
            if logte[i] < lowest:
                lowest, where = logte[i], i
            elif logte[i] - lowest > 0.02:
                return where
    return len(track) - 1


def draw_isochrone(ax, track: pd.DataFrame, family: str, age: float,
                   ticks: bool) -> None:
    tams = split_main_sequence(track)
    ms, post = track.iloc[: tams + 1], track.iloc[tams:]
    if family == "PARSEC":
        ax.plot(ms.colour, ms.G0, "-", color=INK, linewidth=1.0, zorder=3)
        ax.plot(post.colour, post.G0, ":", color=MUTED, linewidth=0.7, zorder=2)
    else:
        ax.plot(ms.colour, ms.G0, "--", color=MUTED, linewidth=0.9, zorder=2)
    if not ticks:
        return
    for mass in TICK_MASSES:
        if mass > ms.Mini.max():
            continue
        x = np.interp(mass, ms.Mini, ms.colour)
        y = np.interp(mass, ms.Mini, ms.G0)
        ax.plot([x], [y], marker="_", markersize=7, color=INK,
                markeredgewidth=1.2, zorder=4)
        ax.text(x + 0.09, y, f"{mass:g}", fontsize=5.8, color=INK,
                va="center", zorder=4)
    end = ms.iloc[-1]
    mto = turnoff_mass(family, age)
    note = (f"MS ends at {end.Mini:.0f} $M_\\odot$\n"
            f"born above {min(mto, IMF_UPPER_LIMIT):.0f} $M_\\odot$: dead"
            if mto < IMF_UPPER_LIMIT else
            f"MS reaches the {IMF_UPPER_LIMIT:.0f} $M_\\odot$ cap\n"
            f"$m_{{\\rm TO}}$ = {mto:.0f} $M_\\odot$: nothing has died")
    ax.annotate(note, xy=(end.colour, end.G0), xytext=(0.30, -9.1),
                fontsize=5.8, color=INK, va="center",
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.6))


def mass_class_scatter(ax, stars: pd.DataFrame, size: float) -> None:
    missing = stars.mass_baseline.isna()
    ax.scatter(stars.loc[missing, "x"], stars.loc[missing, "y"], s=size,
               marker="x", c=MUTED, linewidths=0.5, zorder=3,
               rasterized=True)
    for _, lo, hi, colour in CLASSES:
        sel = stars[stars.mass_baseline.ge(lo) & stars.mass_baseline.lt(hi)]
        ax.scatter(sel.x, sel.y, s=size, c=colour, linewidths=0, zorder=3,
                   rasterized=True)
    spec = stars[stars.is_spectroscopic_anchor.fillna(False).astype(bool)]
    ax.scatter(spec.x, spec.y, s=size * 4.2, facecolors="none",
               edgecolors=INK, linewidths=0.5, zorder=4)


def class_counts(stars: pd.DataFrame) -> str:
    parts = []
    for label, lo, hi, _ in CLASSES[1:]:
        n = int((stars.mass_baseline.ge(lo) & stars.mass_baseline.lt(hi)).sum())
        parts.append(("2–8" if lo == 2.0 else "≥8") + f": {n}")
    return "   ".join(parts)


def style(ax, first_col: bool, bottom_row: bool) -> None:
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.axhline(UMS_LIMIT, color=MUTED, linewidth=0.6, linestyle="--", zorder=1)
    ax.text(XLIM[1] - 0.05, UMS_LIMIT - 0.25, "age fit uses stars above",
            fontsize=5.6, color=MUTED, ha="right")
    if bottom_row:
        ax.set_xlabel("$(BP-RP)_0$")
    if first_col:
        ax.set_ylabel("$M_{G,0}$")


def retained_envelope() -> tuple[float, float]:
    """Span of the retained (measurable, not grid-railed) WP4 age MAPs in the
    chain's own posterior file -- computed, never typed (issue #19)."""
    post = pd.read_parquet(resolve("wp4_age_posteriors"))
    kept = post[post.measurable.astype(bool) & ~post.grid_railed.astype(bool)]
    return float(kept.age_map.min()), float(kept.age_map.max())


def age_strip(ax, norm: pd.DataFrame, lo: float, hi: float) -> None:
    """Counts-based age per branch.  The y axis only separates the subgroups."""
    rows = {"CygOB2-A": 0.0, "CygOB2-B": 1.0, "CygOB2-C": 2.0}
    for subgroup, y in rows.items():
        for family, marker, dy in (("PARSEC", "o", -0.16), ("MIST", "^", 0.16)):
            cell = norm[norm.subgroup.eq(subgroup) & norm.family.eq(family)]
            ages = cell.truth_age_posterior_mean_Myr.to_numpy(float)
            ax.scatter(ages, np.full_like(ages, y + dy), s=13, marker=marker,
                       facecolors="none", edgecolors=SUB[subgroup],
                       linewidths=0.8)
        base = norm[norm.subgroup.eq(subgroup) & norm.family.eq("PARSEC")
                    & norm.R_V.eq(W.BASE["R_V"]) & norm.alpha.eq(W.BASE["alpha"])]
        ax.scatter(base.truth_age_posterior_mean_Myr, [y], s=34, marker="D",
                   color=SUB[subgroup], zorder=5)
        ax.text(5.05, y, tag(subgroup), ha="right", va="center", fontsize=7,
                color=SUB[subgroup], weight="bold")
    ax.axvspan(lo, hi, color="#000000", alpha=0.06, zorder=0)
    ax.text(lo + 0.03, -0.55, f"{lo:.2f}–{hi:.2f} Myr: retained upper-MS "
            "envelope (repair_v5)", fontsize=6.0, color=MUTED, va="center")
    ax.set_xlim(1.8, 5.1)
    ax.set_ylim(-0.8, 2.5)
    ax.set_yticks([])
    ax.set_xlabel("counts-based age (Myr)")
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color="none", markerfacecolor=MUTED,
               markersize=5, label="baseline (PARSEC, $R_V$ = 3.1, $\\alpha$ = 2.3)"),
        Line2D([0], [0], marker="o", color="none", markeredgecolor=MUTED,
               markerfacecolor="none", markersize=5, label="PARSEC branches"),
        Line2D([0], [0], marker="^", color="none", markeredgecolor=MUTED,
               markerfacecolor="none", markersize=5, label="MIST branches"),
    ], fontsize=5.8, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3,
       frameon=False, columnspacing=1.0)


def age_cost(ax, lo: float, hi: float) -> None:
    scan = pd.read_csv(resolve("wp7_age_sensitivity"))
    ax.plot(scan.assumed_age_Myr, scan.N_SN_mean, color=INK, linewidth=1.2)
    ax.fill_between(scan.assumed_age_Myr, scan.N_SN_p16, scan.N_SN_p84,
                    color=INK, alpha=0.12, linewidth=0)
    ax.axvspan(lo, hi, color="#000000", alpha=0.06, zorder=0)
    ax.set_xlabel("common assumed age (Myr)")
    ax.set_ylabel("$N_{\\rm death}$ (all explode)")
    ax.set_title("what the age range costs", fontsize=7, color=MUTED)


def main(cmd_only: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rv = W.BASE["R_V"]
    photo = with_labels(pd.read_parquet(W.frozen("wp3_extinction")))
    masses = pd.read_parquet(resolve("wp4_masses"))[
        ["source_id", "mass_baseline", "is_spectroscopic_anchor"]]
    stars = photo.merge(masses, on="source_id", how="left").assign(
        x=lambda f: f[f"BPRP0_rv{rv:g}"], y=lambda f: f[f"G0_abs_rv{rv:g}"])

    norm = W.normalization()
    base_age = norm[norm.family.eq("PARSEC") & norm.R_V.eq(rv)
                    & norm.alpha.eq(W.BASE["alpha"])].set_index("subgroup")
    mist_age = norm[norm.family.eq("MIST") & norm.R_V.eq(rv)
                    & norm.alpha.eq(W.BASE["alpha"])].set_index("subgroup")
    isochrones = {family: pd.read_parquet(W.frozen(f"wp3_isochrones_{family.lower()}"))
                  for family in w.FAMILIES}

    if cmd_only:
        fig = plt.figure(figsize=(WIDE, 8.3))
        outer = fig.add_gridspec(2, 1, height_ratios=[4.7, 0.42], hspace=0.14)
    else:
        fig = plt.figure(figsize=(WIDE, 10.4))
        outer = fig.add_gridspec(3, 1, height_ratios=[4.7, 0.42, 1.35], hspace=0.16)
    top = outer[0].subgridspec(2, 2, hspace=0.16, wspace=0.08)
    first = fig.add_subplot(top[0, 0])
    axes = [first] + [fig.add_subplot(top[i // 2, i % 2], sharex=first,
                                      sharey=first) for i in (1, 2, 3)]
    for ax in axes[1:]:
        if ax is axes[1] or ax is axes[3]:
            ax.tick_params(labelleft=False)
    for ax in axes[:2]:
        ax.tick_params(labelbottom=False)

    for i, subgroup in enumerate(w.SUBGROUPS):
        ax = axes[i]
        sel = stars[stars.sg.eq(subgroup)]
        age = float(base_age.loc[subgroup].truth_age_posterior_mean_Myr)
        age_m = float(mist_age.loc[subgroup].truth_age_posterior_mean_Myr)
        native, track = nearest_track(isochrones["PARSEC"], age)
        native_m, track_m = nearest_track(isochrones["MIST"], age_m)
        draw_isochrone(ax, track_m, "MIST", native_m, ticks=False)
        draw_isochrone(ax, track, "PARSEC", native, ticks=True)
        mass_class_scatter(ax, sel, 3.2)
        hidden = int((sel.y < YLIM[1]).sum())
        ax.set_title(f"{tag(subgroup)}   PARSEC {age:.2f} Myr · MIST {age_m:.2f} Myr",
                     fontsize=7.2, color=SUB[subgroup], weight="bold")
        ax.text(XLIM[0] + 0.05, 5.8, class_counts(sel)
                + (f"   ({hidden} above the frame)" if hidden else ""),
                fontsize=5.8, color=INK)
        style(ax, i % 2 == 0, i >= 2)

    # The 61 quality-exempt anchors: members, spectroscopic, no subgroup label.
    ax = axes[3]
    rest = stars[stars.sg.notna()]
    ax.scatter(rest.x, rest.y, s=1.0, c="#e3e3e3", linewidths=0, zorder=1,
               rasterized=True)
    exempt = stars[stars.anchor_quality_exempt.fillna(False).astype(bool)]
    shown = exempt[exempt.x.notna() & exempt.y.notna()]
    mass_class_scatter(ax, shown, 7.0)
    native, track = nearest_track(isochrones["PARSEC"], 4.0)
    draw_isochrone(ax, track, "PARSEC", native, ticks=False)
    off = int((shown.x < XLIM[0]).sum())
    ax.set_title(f"unlabelled anchors ({len(shown)} of {len(exempt)} plotted)",
                 fontsize=7.2, color=INK, weight="bold")
    ax.text(XLIM[0] + 0.05, 5.8,
            "not in the published panels; grey = labelled members"
            + (f"; {off} left of the frame" if off else ""),
            fontsize=5.6, color=INK)
    style(ax, False, True)

    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                      markersize=4.5, label=label) for label, _, _, c in CLASSES]
    handles += [
        Line2D([0], [0], marker="x", color="none", markeredgecolor=MUTED,
               markersize=4, label="no mass (no de-reddened colour)"),
        Line2D([0], [0], marker="o", color="none", markeredgecolor=INK,
               markerfacecolor="none", markersize=6,
               label="spectroscopic mass"),
        Line2D([0], [0], color=INK, linewidth=1.0,
               label="PARSEC main sequence, ticks = initial mass ($M_\\odot$)"),
        Line2D([0], [0], color=MUTED, linewidth=0.7, linestyle=":",
               label="PARSEC post-main-sequence"),
        Line2D([0], [0], color=MUTED, linewidth=0.9, linestyle="--",
               label="MIST main sequence at its own fitted age"),
    ]
    legend_ax = fig.add_subplot(outer[1])
    legend_ax.axis("off")
    legend_ax.legend(handles=handles, loc="center", ncol=3, fontsize=6.0,
                     frameon=False)

    # ---- bottom row: the age evidence and what the age range costs ---------
    if not cmd_only:
        lo, hi = retained_envelope()
        bottom = outer[2].subgridspec(1, 3, wspace=0.32)
        age_strip(fig.add_subplot(bottom[0, :2]), norm, lo, hi)
        age_cost(fig.add_subplot(bottom[0, 2]), lo, hi)

    stem = "fig04_cmd_only_draft" if cmd_only else "fig04_cmd_redesign_draft"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight")
    print(f"wrote {OUT / (stem + '.png')}")


if __name__ == "__main__":
    import sys
    main(cmd_only="--cmd-only" in sys.argv[1:])
