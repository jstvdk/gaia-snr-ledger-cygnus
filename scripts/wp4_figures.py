#!/usr/bin/env python3
"""WP4 figures.

  wp4_cmd_subgroups.png   - de-reddened CMD per subgroup, PARSEC+MIST best-age
                            isochrones (single + binary locus), clump highlighted
  wp4_hrd_anchors.png     - spectroscopic HRD (logTe, M_G0): anchors + isochrones
  wp4_age_posteriors.png  - upper-MS age posteriors per subgroup, both families
  wp4_age_summary.png     - age vs subgroup x branch (the money plot), both
                            indicators, with credible intervals
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wp4_common as w
from pathlib import Path

FIGD = w.ROOT / "figures" / "wp4"
FIGD.mkdir(parents=True, exist_ok=True)

C_PAR = "#1f77b4"    # PARSEC
C_MIST = "#d62728"   # MIST
C_A, C_B, C_C = "#4c72b0", "#55a868", "#c44e52"
SUBCOL = {"CygOB2-A": C_A, "CygOB2-B": C_B, "CygOB2-C": C_C}


def _iso_at(iso, family, age, single_only=True):
    a = iso[family]
    uages = np.unique(a["age_Myr"].values)
    anear = uages[np.argmin(np.abs(uages - age))]
    d = a[np.isclose(a["age_Myr"], anear)].sort_values("Mini")
    if single_only:
        d = d[d["label"] <= 1] if family == "PARSEC" else d[d["phase"] <= 0]
        # keep MS dwarfs (log g >= 3.5) and cap at M_ini <= 20 Msun: above that
        # the young turnoff hook swings wildly in colour and the upper MS is
        # colour-degenerate anyway (that regime is shown on the HRD figure with
        # spectroscopic Teff).  This yields a clean ZAMS->turnoff line.
        d = d[(d["logg"] >= 3.5) & (d["Mini"] <= 20.0)]
    col = (d["BP0"] - d["RP0"]).values
    mg = d["G0"].values
    return col, mg, float(anear)


def fig_cmd_subgroups(m, iso, post, clump_ids):
    fit = w.branch_photometry(m, 3.1)
    ubase = post[(post.indicator == "ums") & (post.R_V == 3.1)
                 & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharex=True, sharey=True)
    for ax, sub in zip(axes, w.SUBGROUPS):
        s = fit[fit.subgroup == sub]
        s = s[np.isfinite(s.colour) & np.isfinite(s.MG0)]
        inclump = s["source_id"].isin(clump_ids)
        ax.scatter(s.loc[~inclump, "colour"], s.loc[~inclump, "MG0"], s=9,
                   c="0.55", alpha=0.6, lw=0, label="members")
        ax.scatter(s.loc[inclump, "colour"], s.loc[inclump, "MG0"], s=22,
                   facecolors="none", edgecolors="orange", lw=1.1,
                   label="high-$A_V$ clump")
        for fam, col in [("PARSEC", C_PAR), ("MIST", C_MIST)]:
            age = float(ubase[(ubase.subgroup == sub) & (ubase.family == fam)].age_map.iloc[0])
            c, g, anear = _iso_at(iso, fam, age)
            o = np.argsort(g)
            ax.plot(c[o], g[o], color=col, lw=1.8,
                    label=f"{fam} {anear:.1f} Myr")
            # equal-mass binary locus (0.752 mag brighter)
            ax.plot(c[o], g[o] - 0.752, color=col, lw=1.0, ls="--", alpha=0.6)
        ax.set_title(f"{sub}", fontsize=12)
        ax.set_xlabel(r"$(BP-RP)_0$")
        ax.set_xlim(-1.4, 2.2)
        ax.set_ylim(6.5, -8.5)
        ax.grid(alpha=0.15)
        if sub == "CygOB2-A":
            ax.set_ylabel(r"$M_{G,0}$")
            ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("WP4 de-reddened CMD per subgroup (R_V=3.1); dashed = equal-mass binary locus",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGD / "wp4_cmd_subgroups.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_hrd_anchors(iso, post):
    anc = pd.read_parquet(w.PROC / "wp4_anchor_hrd.parquet")
    ubase = post[(post.indicator == "ums") & (post.R_V == 3.1)
                 & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    # isochrones at the ensemble ages
    for fam, col in [("PARSEC", C_PAR), ("MIST", C_MIST)]:
        age = float(ubase[ubase.family == fam].age_map.median())
        a = iso[fam]
        uages = np.unique(a["age_Myr"].values)
        anear = uages[np.argmin(np.abs(uages - age))]
        d = a[np.isclose(a["age_Myr"], anear)]
        d = d[d["label"] <= 2] if fam == "PARSEC" else d[d["phase"] <= 2]
        d = d.sort_values("Mini")
        ax.plot(d["logTe"], d["G0"], color=col, lw=1.7, label=f"{fam} {anear:.1f} Myr")
    good = anc[anc.MG0_obs.notna() & ~anc.extreme_hot]
    cons = good[good.hrd_consistent]
    inc = good[~good.hrd_consistent]
    ax.scatter(cons.logTe_spec, cons.MG0_obs, s=26, c="k", alpha=0.75,
               label=f"anchors consistent ({len(cons)})")
    ax.scatter(inc.logTe_spec, inc.MG0_obs, s=34, marker="^", c="orange",
               edgecolors="k", lw=0.4,
               label=f"anchors over-luminous / binary ({len(inc)})")
    ext = anc[anc.extreme_hot]
    if len(ext):
        ax.scatter(ext.logTe_spec, ext.MG0_obs, s=60, marker="*", c="magenta",
                   edgecolors="k", label=f"extreme-hot / WR ({len(ext)})")
    ax.invert_xaxis(); ax.invert_yaxis()
    ax.set_xlabel(r"$\log_{10}(T_{\rm eff}/{\rm K})$ (spectroscopic)")
    ax.set_ylabel(r"$M_{G,0}$")
    ax.set_title("WP4 spectroscopic HRD - anchors vs fitted-age isochrones")
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGD / "wp4_hrd_anchors.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_age_posteriors():
    npz = np.load(w.PROC / "wp4_posterior_curves.npz")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
    for ax, sub in zip(axes, w.SUBGROUPS):
        for fam, col in [("PARSEC", C_PAR), ("MIST", C_MIST)]:
            key = f"{sub}|{fam}|ums"
            if key + "|age" in npz.files:
                ax.plot(npz[key + "|age"], npz[key + "|post"], color=col, lw=2,
                        label=f"{fam} (upper-MS)")
            kp = f"{sub}|{fam}|pms"
            if kp + "|age" in npz.files:
                ax.plot(npz[kp + "|age"], npz[kp + "|post"], color=col, lw=1.3,
                        ls=":", alpha=0.8, label=f"{fam} (PMS turn-on)")
        ax.set_title(sub); ax.set_xlabel("age (Myr)"); ax.set_xlim(1, 10)
        ax.grid(alpha=0.15)
        if sub == "CygOB2-A":
            ax.set_ylabel("posterior density"); ax.legend(fontsize=8)
    fig.suptitle("WP4 age posteriors (baseline R_V=3.1, f_bin=0.4); solid=upper-MS, dotted=PMS turn-on")
    fig.tight_layout()
    fig.savefig(FIGD / "wp4_age_posteriors.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_age_summary(post):
    """Money plot: MAP age with 68% CI per subgroup, both families, both
    indicators, plus the full branch envelope (R_V x f_bin) as a light band."""
    fig, ax = plt.subplots(figsize=(10, 5.8))
    xpos = {"CygOB2-A": 0, "CygOB2-B": 1, "CygOB2-C": 2}
    offs = {("PARSEC", "ums"): -0.18, ("MIST", "ums"): -0.06,
            ("PARSEC", "pms"): +0.06, ("MIST", "pms"): +0.18}
    styles = {("PARSEC", "ums"): (C_PAR, "o", "PARSEC upper-MS"),
              ("MIST", "ums"): (C_MIST, "o", "MIST upper-MS"),
              ("PARSEC", "pms"): (C_PAR, "s", "PARSEC PMS turn-on"),
              ("MIST", "pms"): (C_MIST, "s", "MIST PMS turn-on")}
    seen = set()
    for sub in w.SUBGROUPS:
        for (fam, ind), (col, mk, lab) in styles.items():
            base = post[(post.subgroup == sub) & (post.family == fam)
                        & (post.indicator == ind) & (post.R_V == 3.1)
                        & (post.f_bin == 0.4) & (post.dmu == 0.0)]
            if not len(base):
                continue
            r = base.iloc[0]
            x = xpos[sub] + offs[(fam, ind)]
            meas = bool(r.measurable)
            lo, hi = r.age_lo68, r.age_hi68
            mfc = col if meas else "white"
            ax.errorbar(x, r.age_map, yerr=[[r.age_map - lo], [hi - r.age_map]],
                        fmt=mk, color=col, mfc=mfc, ms=8, capsize=3, lw=1.5,
                        label=lab if lab not in seen else None,
                        alpha=1.0 if meas else 0.45)
            seen.add(lab)
            # full R_V x f_bin envelope for the upper-MS as a faint bar
            if ind == "ums":
                env = post[(post.subgroup == sub) & (post.family == fam)
                           & (post.indicator == "ums") & (post.dmu == 0.0)]
                ax.plot([x, x], [env.age_map.min(), env.age_map.max()],
                        color=col, lw=6, alpha=0.15, solid_capstyle="round")
    ax.axhspan(2, 3, color="0.8", alpha=0.3, label="Wright+15 non-rotating (2-3)")
    ax.axhspan(4, 5, color="0.6", alpha=0.2, label="Wright+15 rotating (4-5)")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(w.SUBGROUPS)
    ax.set_ylabel("age (Myr)"); ax.set_ylim(1, 7)
    ax.set_title("WP4 subgroup ages - filled=measurable, open=excluded; "
                 "thick faint bar = R_V x f_bin envelope")
    ax.grid(alpha=0.15, axis="y")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGD / "wp4_age_summary.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def run():
    m = w.load_members()
    iso = w.load_isochrones()
    post = pd.read_parquet(w.PROC / "wp4_age_posteriors.parquet")
    clump = pd.read_parquet(w.PROC / "wp4_clump.parquet")
    clump_ids = set(clump["source_id"])
    fig_cmd_subgroups(m, iso, post, clump_ids)
    fig_hrd_anchors(iso, post)
    fig_age_posteriors()
    fig_age_summary(post)
    print("wrote figures to", FIGD)
    for p in sorted(FIGD.glob("*.png")):
        print("  ", p.name, p.stat().st_size, "bytes")


if __name__ == "__main__":
    run()
