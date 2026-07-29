#!/usr/bin/env python3
"""WP7 figures: the age-sensitivity honesty plot and its companions.

Four panels, each answering a question the ledger must not bury:

  a  N_SN against ASSUMED age -- the mandatory honesty plot.  Everything
     inherits the WP4 ages, and this shows exactly how much.
  b  N_SN against the black-hole threshold -- the explodability dependence,
     scanned rather than hidden behind a branch label.
  c  R_SN(t), explosions per Myr against look-back time.
  d  the branch spread on N_SN, which is the honest error bar.

Outputs: figures/wp7/wp7_ledger_panels.png

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp7_figures.py
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wp5_common as w

OUT = w.ROOT / "figures" / "wp7"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    age = pd.read_csv(w.TABLES / "wp7_age_sensitivity.csv")
    bh = pd.read_csv(w.TABLES / "wp7_bh_threshold_scan.csv")
    rsn = pd.read_csv(w.TABLES / "wp7_rsn_curves.csv")
    ledger = pd.read_csv(w.TABLES / "wp7_ledger.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    fig.suptitle(
        "WP7 — the supernova ledger of Cyg OB2 "
        "(repair_v7, 2,000,000 iterations)",
        fontsize=13, fontweight="bold",
    )

    # ---- a: the honesty plot ----------------------------------------------
    ax = axes[0, 0]
    ax.fill_between(
        age.assumed_age_Myr, age.N_SN_p16, age.N_SN_p84,
        alpha=0.25, color="#c0392b", label="68% credible",
    )
    ax.plot(age.assumed_age_Myr, age.N_SN_mean, color="#c0392b", lw=2,
            label="mean $N_{\\rm SN}$")
    for value, name, colour in (
        (4.00, "A", "#2c3e50"), (4.09, "B", "#16a085"), (2.52, "C", "#8e44ad"),
    ):
        ax.axvline(value, ls="--", lw=1.1, color=colour, alpha=0.8)
        ax.text(value, ax.get_ylim()[1] * 0.55, name, ha="center", fontsize=10,
                color=colour, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=colour,
                          alpha=0.95))
    ax.set_xlabel("assumed age applied to all three subgroups  [Myr]")
    ax.set_ylabel("$N_{\\rm SN}$, association total")
    ax.set_title("(a) everything inherits the WP4 ages", fontsize=11)
    ax.legend(fontsize=9, loc="center right")
    ax.grid(alpha=0.25)

    # ---- b: explodability --------------------------------------------------
    ax = axes[0, 1]
    total = bh.groupby("bh_threshold_Msun").N_SN_mean.sum()
    ax.plot(total.index, total.values, "o-", color="#2980b9", lw=2)
    ax.axvspan(0, 40, color="#7f8c8d", alpha=0.18)
    ax.text(21, total.values.max() * 0.55,
            "no supernovae:\nnothing below\n40 M$_\\odot$ has died",
            fontsize=9, ha="center", color="#2c3e50")
    ax.set_xlabel("black-hole threshold  [M$_\\odot$]")
    ax.set_ylabel("$N_{\\rm SN}$, association total")
    ax.set_title("(b) the whole budget lives above 40 M$_\\odot$", fontsize=11)
    ax.grid(alpha=0.25)

    # ---- c: R_SN(t) --------------------------------------------------------
    ax = axes[1, 0]
    grouped = rsn.groupby(["lookback_lo_Myr", "lookback_hi_Myr"]).rate_per_Myr
    curve = grouped.sum().reset_index()
    mid = 0.5 * (curve.lookback_lo_Myr + curve.lookback_hi_Myr)
    ax.step(mid, curve.rate_per_Myr, where="mid", color="#27ae60", lw=1.8)
    ax.fill_between(mid, 0, curve.rate_per_Myr, step="mid", alpha=0.25,
                    color="#27ae60")
    ax.set_xlim(0, 1.6)
    ax.invert_xaxis()
    ax.set_xlabel("look-back time  [Myr]   (present day at right)")
    ax.set_ylabel("$R_{\\rm SN}$  [Myr$^{-1}$]")
    ax.set_title("(c) explosion history, baseline branch", fontsize=11)
    ax.grid(alpha=0.25)

    # ---- d: branch spread --------------------------------------------------
    ax = axes[1, 1]
    assoc = ledger[
        ledger.scope.eq("association") & ledger.explodability.eq("all_explode")
    ]
    positions, labels, colours = [], [], []
    data = []
    for i, alpha in enumerate(sorted(assoc.alpha.unique())):
        block = assoc[assoc.alpha.eq(alpha)]
        data.append(block.N_SN_mean.to_numpy())
        positions.append(i)
        labels.append(f"$\\alpha$ = {alpha}")
        colours.append("#e67e22" if alpha == 2.3 else "#95a5a6")
    parts = ax.violinplot(data, positions=positions, showmedians=True,
                          widths=0.7)
    for body, colour in zip(parts["bodies"], colours):
        body.set_facecolor(colour)
        body.set_alpha(0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("$N_{\\rm SN}$, association total")
    ax.set_title("(d) the IMF slope is the dominant lever", fontsize=11)
    ax.grid(alpha=0.25, axis="y")
    ax.text(
        0.5, 0.96,
        "WP6's census disfavours $\\alpha$ = 2.6\nand closes at $\\alpha$ = 2.25",
        transform=ax.transAxes, ha="center", va="top", fontsize=9,
        bbox=dict(boxstyle="round", fc="#fdf2e9", ec="#e67e22", alpha=0.9),
    )

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = OUT / "wp7_ledger_panels.png"
    fig.savefig(path, dpi=160)
    print(f"wrote {path.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
