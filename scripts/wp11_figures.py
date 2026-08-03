#!/usr/bin/env python3
"""WP11 Part B figure -- the 60Fe forecast against COSI.

The one WP11 result with no figure in the paper.  Built for the conference talk
and reusing the paper's style and colour mapping exactly: alpha = 2.0 is
#3B6EA8 and alpha = 2.3 is #C4622D in Fig. 6, so they stay that way here.
Colour follows the entity, not the chart.

UNITS -- the one thing that is easy to get wrong here:
  * COSI's 3-sigma narrow-line sensitivity is quoted PER LINE (Tomsick+2023).
  * Martin+2009's SPI Cygnus upper limit is for the COMBINED 1173+1332 keV
    emission.
The plot is in per-line 1173 keV flux, so the SPI limit is drawn at half its
published value and the label says so.  Mixing the two conventions on one axis
would overstate the SPI margin by a factor 2.

Output: figures/wp11/wp11_cosi_forecast.{pdf,png}

Run:
  PYTHONPATH=scripts python3 scripts/wp11_figures.py
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wp5_common as w

OUT = w.ROOT / "figures" / "wp11"
PRIMARY_ARM = "LC06_NL"
COLOUR = {2.0: "#3B6EA8", 2.3: "#C4622D"}
MARKER = {2.0: "o", 2.3: "s"}

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 200,
})


def main() -> None:
    table = pd.read_csv(w.TABLES / "wp11_isotope_forecast.csv")
    prereg = json.loads(
        (w.PROVENANCE / "wp11_isotope_prereg.json").read_text()
    )
    cosi = prereg["instruments"]["COSI_narrow_line_3sigma_2yr"]["value_ph_cm2_s"]
    spi_combined = prereg[
        "instruments"]["SPI_Cygnus_60Fe_upper_limit"]["value_ph_cm2_s"]
    spi_per_line = spi_combined / 2.0

    head = table[table.in_headline_set & table.yield_arm.eq(PRIMARY_ARM)]
    langer = table[table.in_headline_set & table.yield_arm.eq("LC06_Langer")]

    fig, ax = plt.subplots(figsize=(7.6, 4.4))

    # Two yield arms shown directly rather than described: the gap between the
    # filled and hollow clusters IS the result (I4 -- the yield model moves the
    # forecast further than the whole branch set does).
    for i, alpha in enumerate((2.0, 2.3)):
        for frame, arm_name in ((True, PRIMARY_ARM), (False, "LC06_Langer")):
            source = head if frame else langer
            arm = source[source.alpha.eq(alpha)].sort_values("F_1173_ph_cm2_s")
            x = i + np.linspace(-0.26, 0.26, len(arm))
            ax.scatter(
                x, arm.F_1173_ph_cm2_s, s=34, marker=MARKER[alpha], zorder=3,
                color=COLOUR[alpha] if frame else "none",
                edgecolor="none" if frame else COLOUR[alpha],
                linewidth=0 if frame else 1.1,
            )

    ax.axhline(cosi, color="0.20", lw=1.6, zorder=2)
    ax.text(-1.02, cosi * 1.13, "COSI  $3\\sigma$ narrow line, 2-yr survey",
            fontsize=8, color="0.20")
    ax.axhline(spi_per_line, color="0.45", lw=1.1, ls="--", zorder=2)
    ax.text(-1.02, spi_per_line * 1.13,
            "INTEGRAL/SPI upper limit, Cygnus (per line)",
            fontsize=8, color="0.45")

    ax.set_yscale("log")
    ax.set_ylim(3.5e-8, 3.2e-5)
    ax.set_xlim(-1.08, 1.62)

    # The null arm cannot be drawn on a log axis: it is exactly zero.
    ax.text(1.58, 1.15e-5,
            "LC18 recommended yields: everything\n"
            "$>25\\,M_\\odot$ collapses to a black hole\n"
            r"$\Rightarrow$ predicted flux is identically ZERO",
            fontsize=7.5, color="0.35", ha="right", va="bottom",
            linespacing=1.5)

    handles = [
        plt.Line2D([], [], marker="o", ls="none", color="0.35",
                   label="LC06 — primary yield arm"),
        plt.Line2D([], [], marker="o", ls="none", mfc="none", mec="0.35",
                   color="0.35", label="LC06 — Langer WR mass loss"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False,
              bbox_to_anchor=(0.0, 0.02))

    n_above = int((head[head.alpha.eq(2.0)].F_1173_ph_cm2_s >= cosi).sum())
    n_total = int(head.alpha.eq(2.0).sum())
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"$\alpha = 2.0$" + f"\n{n_total} branches",
                        r"$\alpha = 2.3$" + f"\n{n_total} branches"],
                       fontsize=10)
    ax.set_ylabel(r"$^{60}$Fe 1173 keV flux  (ph cm$^{-2}$ s$^{-1}$)")
    ax.set_title(
        f"COSI separates what Gaia cannot: {n_above}/{n_total} branches "
        f"detectable at $\\alpha=2.0$, 0/{n_total} at $\\alpha=2.3$\n"
        "— but the yield model moves the prediction further than the "
        "branch set does",
        fontsize=9.5, loc="left", pad=10,
    )
    ax.grid(axis="y", color="0.92", lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        path = OUT / f"wp11_cosi_forecast.{suffix}"
        fig.savefig(path, bbox_inches="tight")
        print(f"wrote {path.relative_to(w.ROOT)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
