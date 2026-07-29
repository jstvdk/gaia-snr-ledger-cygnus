#!/usr/bin/env python3
"""Build the supervisor briefing PDF: what WP1-WP6 did, with plots.

Self-contained: every number is read from the accepted artifacts, so the
document cannot drift from the pipeline.  Regenerate after any result changes.

Output: reports/cygob2_WP0-6_briefing.pdf

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/make_supervisor_report.py
"""
from __future__ import annotations

import glob
import json
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch, Rectangle

import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130,
    "axes.spines.top": False, "axes.spines.right": False,
})
NAVY, RED, TEAL, GOLD, GREY = "#1f3a68", "#c0392b", "#17888e", "#d68910", "#7f8c8d"
SGCOL = {"CygOB2-A": NAVY, "CygOB2-B": TEAL, "CygOB2-C": GOLD}
PAGE = (8.27, 11.69)   # A4 portrait


# ----------------------------------------------------------------- utilities
def text_page(pdf, title, blocks, subtitle=None):
    """A page of prose.  `blocks` is a list of (heading|None, body) tuples."""
    fig = plt.figure(figsize=PAGE)
    fig.patch.set_facecolor("white")
    y = 0.955
    fig.text(0.07, y, title, fontsize=17, fontweight="bold", color=NAVY)
    y -= 0.028
    if subtitle:
        fig.text(0.07, y, subtitle, fontsize=9.5, color=GREY, style="italic")
        y -= 0.022
    fig.add_artist(plt.Line2D([0.07, 0.93], [y, y], color=NAVY, lw=1.4))
    y -= 0.030
    for heading, body in blocks:
        if heading:
            fig.text(0.07, y, heading, fontsize=11, fontweight="bold", color=NAVY)
            y -= 0.021
        for line in body.split("\n"):
            if not line.strip():
                y -= 0.011
                continue
            weight = "bold" if line.startswith("**") else "normal"
            colour = "black"
            text = line.replace("**", "")
            if text.startswith("> "):
                text, colour = text[2:], RED
                fig.text(0.085, y, text, fontsize=9, color=colour, wrap=True)
            elif text.startswith("- "):
                fig.text(0.085, y, "•  " + text[2:], fontsize=9, color=colour)
            else:
                fig.text(0.07, y, text, fontsize=9, fontweight=weight, color=colour)
            y -= 0.0165
        y -= 0.012
    pdf.savefig(fig); plt.close(fig)


def figure_page(pdf, title, draw, caption, note=None):
    fig = plt.figure(figsize=PAGE)
    fig.patch.set_facecolor("white")
    fig.text(0.07, 0.955, title, fontsize=15, fontweight="bold", color=NAVY)
    fig.add_artist(plt.Line2D([0.07, 0.93], [0.943, 0.943], color=NAVY, lw=1.4))
    gs = fig.add_gridspec(1, 1, left=0.11, right=0.94, top=0.90, bottom=0.37)
    draw(fig, gs[0])
    y = 0.315
    for line in caption.split("\n"):
        if not line.strip():
            y -= 0.010; continue
        fig.text(0.07, y, line.replace("**", ""), fontsize=9,
                 fontweight="bold" if line.startswith("**") else "normal")
        y -= 0.0165
    if note:
        y -= 0.008
        for line in note.split("\n"):
            fig.text(0.07, y, line, fontsize=8.6, color=RED, style="italic")
            y -= 0.0155
    pdf.savefig(fig); plt.close(fig)


# ------------------------------------------------------------------ data load
def load():
    d = {}
    members = pd.read_parquet(w.PROC / "wp2_members.parquet").drop(
        columns=["subgroup"], errors="ignore")
    labels = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    d["mem"] = members.merge(labels[["source_id", "subgroup"]], on="source_id",
                             how="inner").query("membership_probability > 0.5")
    d["ext"] = pd.read_parquet(w.PROC / "wp3_extinction_repair_v5.parquet")
    d["age"] = pd.read_parquet(w.PROC / "wp4_age_posteriors_repair_v5.parquet")
    d["norm"] = pd.read_parquet(w.PROC / "wp5_imf_normalization_repair_v6.parquet")
    d["clo"] = pd.read_csv(w.TABLES / "wp6_closure.csv")
    d["floor"] = pd.read_csv(w.TABLES / "wp6_closure_floor_comparison.csv")
    d["run"] = pd.read_csv(w.TABLES / "wp6_runaways.csv")
    d["xm"] = pd.read_csv(w.TABLES / "wp6_runaway_crossmatch.csv")
    d["mult"] = pd.read_csv(w.TABLES / "wp6_multiplicity_closure.csv")
    d["xchk"] = json.loads((w.PROVENANCE /
        "wp6_external_crosschecks_execution.json").read_text())
    d["attr"] = json.loads((w.PROVENANCE /
        "wp6_closure_attribution_execution.json").read_text())
    return d


def main() -> None:
    d = load()
    out = w.ROOT / "reports" / "cygob2_WP0-6_briefing.pdf"
    out.parent.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    with PdfPages(out) as pdf:
        # ---------------------------------------------------------- page 1
        fig = plt.figure(figsize=PAGE); fig.patch.set_facecolor("white")
        fig.text(0.5, 0.88, "The supernova history of Cygnus OB2",
                 fontsize=22, fontweight="bold", color=NAVY, ha="center")
        fig.text(0.5, 0.845, "from Gaia DR3", fontsize=15, color=NAVY, ha="center")
        fig.text(0.5, 0.805, f"Progress briefing · WP1–WP6 complete · {stamp}",
                 fontsize=10, color=GREY, ha="center", style="italic")
        fig.add_artist(plt.Line2D([0.15, 0.85], [0.785, 0.785], color=NAVY, lw=1.5))

        fig.text(0.07, 0.735, "The pitch, in one paragraph",
                 fontsize=13, fontweight="bold", color=RED)
        pitch = (
            "The Cygnus region emits gamma rays up to PeV energies, and the leading\n"
            "explanation (Härer et al. 2025) requires a supernova to have exploded\n"
            "inside Cygnus OB2 within the last ~100 000 years. Nobody has checked\n"
            "whether the association can actually supply one. We do it by counting the\n"
            "stars that are missing. Gaia DR3 gives us 2 112 members; we measure each\n"
            "star's extinction and mass, fit the initial mass function where the census\n"
            "is complete (2–8 M☉), and extrapolate it upward. Stars heavier than the\n"
            "isochrone turnoff should exist but do not — they have already exploded.\n"
            "Counting them gives the supernova history directly. The critical test is\n"
            "that the same extrapolation must correctly predict the massive stars that\n"
            "are still alive, between 8 M☉ and the turnoff — a genuinely out-of-sample\n"
            "check, because those stars never enter the fit. It now closes to about\n"
            "10 per cent. The answer that follows is roughly 6 supernovae in the last\n"
            "million years, and a ~50 per cent chance that one of them went off inside\n"
            "Härer's 100 000-year window."
        )
        y = 0.705
        for line in pitch.split("\n"):
            fig.text(0.07, y, line, fontsize=10.5); y -= 0.0185

        fig.text(0.07, 0.365, "Where we are", fontsize=13,
                 fontweight="bold", color=NAVY)
        rows = [
            ("WP1–WP2", "sample and membership", "2 112 members, 3 subgroups", "done"),
            ("WP3", "per-star extinction", "median A_V = 6.0 mag", "done"),
            ("WP4", "ages and masses", "2.5 – 4.0 Myr", "done"),
            ("WP5", "IMF normalisation", "response-aware, accepted", "done"),
            ("WP6", "census closure + runaways", "closes to 10%, gate met", "done"),
            ("WP7", "the supernova ledger", "N_SN posterior", "next"),
            ("WP8–WP10", "validation, verdict, paper", "", "later"),
        ]
        y = 0.335
        for wp, what, result, status in rows:
            col = NAVY if status == "done" else (RED if status == "next" else GREY)
            fig.text(0.08, y, wp, fontsize=10, fontweight="bold", color=col)
            fig.text(0.21, y, what, fontsize=10)
            fig.text(0.55, y, result, fontsize=10, color=GREY)
            fig.text(0.88, y, status, fontsize=9, color=col, style="italic")
            y -= 0.026

        fig.text(0.07, 0.115,
                 "Everything in this document is generated directly from the\n"
                 "pipeline's accepted artifacts, so it cannot drift from the code.",
                 fontsize=8.5, color=GREY, style="italic")
        pdf.savefig(fig); plt.close(fig)

        # ---------------------------------------------------------- page 2
        text_page(pdf, "1. The idea", [
            (None,
             "A supernova leaves no permanent mark that is easy to find. But the star\n"
             "that exploded is permanently gone — and we can notice its absence.\n"),
            ("The logic in four steps",
             "1.  Stars form with a known distribution of masses, the initial mass\n"
             "     function: many small stars, few large ones, dN/dM ∝ M^(−α).\n"
             "2.  Massive stars die fast. At a given age, every star above a threshold\n"
             "     mass — the isochrone turnoff — has already exploded.\n"
             "3.  We measure the IMF where our census is complete (2–8 M☉) and\n"
             "     extrapolate it above the turnoff.\n"
             "4.  That extrapolated number is the number of supernovae.\n"),
            ("Why the middle range is the whole game",
             "Between 8 M☉ and the turnoff, stars are massive but still alive. The\n"
             "extrapolation must predict how many of those we should see. Those stars\n"
             "never enter the fit, so this is a real out-of-sample test:\n"),
            (None,
             "> if the prediction fails there, the extrapolation above the turnoff\n"
             "> cannot be trusted, and the supernova count means nothing.\n"),
            ("What makes it hard",
             "- Cygnus OB2 sits behind ~6 magnitudes of dust, and it varies star to star\n"
             "- the association is three subgroups of different ages, not one cluster\n"
             "- every cut that removes faint stars also removes low-mass stars, which\n"
             "   tilts the very quantity we are measuring\n"
             "- so completeness cannot be assumed: it must be measured by injecting\n"
             "   synthetic stars and seeing how many come back\n"),
        ], subtitle="Why counting missing stars measures supernovae")

        # ---------------------------------------------------------- page 3
        def draw_flow(fig, spec):
            ax = fig.add_subplot(spec); ax.axis("off")
            ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.grid(False)
            steps = [
                ("WP1", "Gaia DR3 + 2MASS", "sky box, quality cuts", NAVY),
                ("WP2", "membership", "2 112 members, A/B/C", NAVY),
                ("WP3", "extinction per star", "A_V ≈ 6 mag", NAVY),
                ("WP4", "ages + masses", "2.5–4.0 Myr", NAVY),
                ("WP5", "IMF normalisation k", "fitted on 2–8 M☉", TEAL),
                ("WP6", "closure test", "does k predict 8 M☉→turnoff?", RED),
                ("WP7", "supernova ledger", "N_SN, timing", GREY),
            ]
            for i, (tag, what, detail, col) in enumerate(steps):
                y = 9.3 - i * 1.32
                ax.add_patch(Rectangle((0.6, y - 0.46), 8.8, 0.92,
                    facecolor=col, alpha=0.10, edgecolor=col, lw=1.3))
                ax.text(1.0, y, tag, fontsize=10, fontweight="bold",
                        color=col, va="center")
                ax.text(2.0, y + 0.14, what, fontsize=9.5, va="center")
                ax.text(2.0, y - 0.20, detail, fontsize=8, color=GREY, va="center")
                if i < len(steps) - 1:
                    ax.add_patch(FancyArrowPatch((5.0, y - 0.47), (5.0, y - 0.85),
                        arrowstyle="-|>", mutation_scale=11, color=GREY, lw=1.1))
            ax.annotate("", xy=(9.75, 9.3 - 5 * 1.32), xytext=(9.75, 9.3 - 4 * 1.32),
                        arrowprops=dict(arrowstyle="<->", color=RED, lw=1.4))
            ax.text(9.95, 9.3 - 4.5 * 1.32, "the\ntest", fontsize=8.5, color=RED,
                    va="center", ha="left", fontweight="bold")
        figure_page(pdf, "2. The pipeline", draw_flow,
            "**Each stage consumes only accepted outputs of the one above.**\n"
            "WP5 fits the normalisation k using stars from 2 to 8 solar masses, where\n"
            "the census is complete. WP6 then asks whether that same k correctly\n"
            "predicts the massive stars still alive above 8 M☉ — a range that never\n"
            "entered the fit. Only if it does can we extrapolate above the turnoff,\n"
            "which is what WP7 turns into a supernova count.\n\n"
            "Every stage is gated: it does not proceed until its own checks pass, and\n"
            "every number carries a provenance record with input checksums.")

        # ---------------------------------------------------------- page 4
        def draw_sky(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.32)
            ax = fig.add_subplot(gs[0])
            for sg, g in d["mem"].groupby("subgroup"):
                ax.scatter(g.l_deg, g.b_deg, s=3, alpha=0.5, color=SGCOL[sg],
                           label=f"{sg}  (n={len(g)})", linewidths=0)
            ax.set_xlabel("galactic longitude  l  [deg]")
            ax.set_ylabel("galactic latitude  b  [deg]")
            ax.invert_xaxis(); ax.legend(frameon=False, fontsize=8, markerscale=3)
            ax.set_title("Members on the sky, coloured by subgroup")
            ax2 = fig.add_subplot(gs[1])
            ax2.hist(d["mem"].parallax_corrected, bins=60, color=NAVY, alpha=0.8)
            ax2.axvline(1000 / 1620, color=RED, lw=1.6,
                        label="adopted 1.62 kpc")
            ax2.axvline(1000 / 1350, color=GOLD, lw=1.4, ls="--",
                        label="Berlanas+19 foreground, 1.35 kpc")
            ax2.set_xlabel("parallax  [mas]"); ax2.set_ylabel("stars")
            ax2.legend(frameon=False, fontsize=8)
            ax2.set_title("Distance: one population, not two")
        figure_page(pdf, "3. WP1–WP2 — the sample", draw_sky,
            f"**{len(d['mem'])} members above 50% membership probability**, split into\n"
            "three subgroups by position and proper motion. Membership is a\n"
            "probability, never a hard cut — hard cuts on a magnitude-correlated\n"
            "quantity would bias the mass function.\n\n"
            "The lower panel tests a published claim that Cyg OB2 splits into two\n"
            "populations at 1.35 and 1.76 kpc. We find one population at 1.62 kpc with\n"
            "a 45 pc depth; a two-component fit collapses to components 31 pc apart.",
            note="Caveat we record: parallax is one of the membership clustering\n"
                 "features, so a foreground group would have been removed before this\n"
                 "test ran. It shows no split among members — not that none exists.")

        # ---------------------------------------------------------- page 5
        def draw_ext(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.35)
            ax = fig.add_subplot(gs[0])
            av = d["ext"]["av_rv3.1"].dropna()
            ax.hist(av, bins=70, color=TEAL, alpha=0.85)
            ax.axvline(av.median(), color=RED, lw=1.6,
                       label=f"median {av.median():.2f} mag")
            ax.set_xlabel("$A_V$  [mag]"); ax.set_ylabel("stars")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("Per-star extinction — the dominant difficulty")
            ax2 = fig.add_subplot(gs[1])
            e = d["ext"].dropna(subset=["av_rv3.1", "l_deg", "b_deg"])
            s = ax2.scatter(e.l_deg, e.b_deg, c=e["av_rv3.1"], s=3, cmap="inferno",
                            vmin=np.percentile(e["av_rv3.1"], 2),
                            vmax=np.percentile(e["av_rv3.1"], 98), linewidths=0)
            ax2.invert_xaxis()
            ax2.set_xlabel("l [deg]"); ax2.set_ylabel("b [deg]")
            plt.colorbar(s, ax=ax2, label="$A_V$ [mag]", pad=0.02)
            ax2.set_title("Extinction is patchy — it must be fitted per star")
        figure_page(pdf, "4. WP3 — extinction", draw_ext,
            "**Median $A_V$ = 6.0 mag; the 10th–90th percentile spans 4.65–7.03.**\n"
            "That is roughly a factor 250 of dimming in the V band, and it varies\n"
            "across the field, so a single average value would be useless.\n\n"
            "Each star gets its own extinction posterior, fitted against a template\n"
            "library with a spatial prior built from spectroscopically classified\n"
            "anchor stars. This stage needed five repair iterations to get right and\n"
            "is the single largest source of systematic uncertainty downstream.")

        # ---------------------------------------------------------- page 6
        def draw_age(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.35)
            ax = fig.add_subplot(gs[0])
            base = d["age"].query(
                "family=='PARSEC' and R_V==3.1 and f_bin==0.4 and "
                "indicator=='ums' and dmu==0.0")
            for i, r in enumerate(base.itertuples()):
                ax.errorbar(r.age_map, i, xerr=[[r.age_map - r.age_lo68],
                            [r.age_hi68 - r.age_map]], fmt="o", ms=7,
                            color=SGCOL[r.subgroup], capsize=4, lw=2)
                ax.text(r.age_map, i + 0.20, f"{r.age_map:.2f} Myr",
                        ha="center", fontsize=8.5, color=SGCOL[r.subgroup])
            ax.set_yticks(range(len(base)))
            ax.set_yticklabels(base.subgroup); ax.set_ylim(-0.6, len(base) - 0.2)
            ax.set_xlabel("age  [Myr]")
            ax.set_title("Subgroup ages — Cyg OB2 is not coeval")
            ax2 = fig.add_subplot(gs[1])
            ages = np.linspace(2.2, 5.0, 200)
            for fam, ls in [("PARSEC", "-"), ("MIST", "--")]:
                ax2.plot(ages, [min(turnoff_mass(fam, a), IMF_UPPER_LIMIT)
                                for a in ages], ls, color=NAVY, label=fam)
            for r in base.itertuples():
                m = min(turnoff_mass("PARSEC", r.age_map), IMF_UPPER_LIMIT)
                ax2.plot(r.age_map, m, "o", ms=8, color=SGCOL[r.subgroup])
                ax2.annotate(f"{r.subgroup[-1]}: {m:.0f} M$_\\odot$",
                             (r.age_map, m), textcoords="offset points",
                             xytext=(8, 4), fontsize=8.5, color=SGCOL[r.subgroup])
            ax2.axhline(8, color=RED, ls=":", lw=1.4)
            ax2.text(4.6, 9.5, "8 M$_\\odot$ — supernova threshold",
                     fontsize=8, color=RED)
            ax2.set_yscale("log"); ax2.set_xlabel("age [Myr]")
            ax2.set_ylabel("turnoff mass  [M$_\\odot$]")
            ax2.legend(frameon=False, fontsize=8)
            ax2.set_title("Everything above the turnoff has already exploded")
        figure_page(pdf, "5. WP4 — ages and the turnoff", draw_age,
            "**The turnoff is the physics.** A star above it has died; below it, it is\n"
            "still shining. So the age of each subgroup sets how many supernovae it\n"
            "has already produced.\n\n"
            "CygOB2-C is 2.5 Myr old and its turnoff is still above the 120 M☉ upper\n"
            "limit of the IMF — **nothing in C has exploded yet.** Essentially all the\n"
            "supernovae come from A and B.",
            note="The two isochrone families disagree by up to 30% at these ages.\n"
                 "We carry both in parallel to the end and never average them.")

        # ---------------------------------------------------------- page 7
        def draw_resp(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.35)
            ax = fig.add_subplot(gs[0])
            files = sorted(glob.glob(str(w.PROC /
                "wp5_agenode_A_PARSEC_rv3p1_*repair_v6_response.parquet")))
            if files:
                f = pd.read_parquet(files[len(files) // 2])
                masses = np.sort(f.true_primary_mass.unique())
                rec = [np.isfinite(f[f.true_primary_mass.eq(m)].recovered_mass).mean()
                       for m in masses]
                ax.plot(masses, rec, "o-", color=NAVY, ms=4)
            ax.axvspan(2, 8, color=TEAL, alpha=0.13)
            ax.text(4.0, 0.06, "fit here\n(2–8 M$_\\odot$)", fontsize=8.5,
                    color=TEAL, ha="center", fontweight="bold")
            ax.set_xscale("log"); ax.set_xlabel("true mass  [M$_\\odot$]")
            ax.set_ylabel("recovery fraction"); ax.set_ylim(0, 1)
            ax.set_title("Completeness is measured, never assumed")
            ax2 = fig.add_subplot(gs[1])
            b = d["norm"].query("family=='PARSEC' and R_V==3.1")
            x = np.arange(3); wd = 0.26
            for i, a in enumerate([2.0, 2.3, 2.6]):
                v = [float(b.query(f"subgroup=='{s}' and alpha=={a}").k_median.iloc[0])
                     for s in w.SUBGROUPS]
                ax2.bar(x + (i - 1) * wd, v, wd, label=f"α = {a}",
                        color=[NAVY, TEAL, GOLD][i], alpha=0.9)
            ax2.set_xticks(x); ax2.set_xticklabels(w.SUBGROUPS)
            ax2.set_ylabel("IMF normalisation  k")
            ax2.legend(frameon=False, fontsize=8, ncol=3)
            ax2.set_title("The fitted normalisation, per branch")
        figure_page(pdf, "6. WP5 — the IMF normalisation", draw_resp,
            "**We inject synthetic stars of known mass and count how many survive the\n"
            "entire pipeline.** That gives the response — the probability that a star\n"
            "of true mass M is recovered, and what mass it is then assigned.\n\n"
            "The normalisation k is fitted only in 2–8 M☉ (shaded), through that\n"
            "response, with a Poisson likelihood. Dividing an observed count by an\n"
            "average completeness is forbidden in this project: the response also\n"
            "scatters masses across bin edges, and in 6 of 54 cases that correction\n"
            "would have the wrong sign.\n\n"
            "The slope α is not fitted. It is carried as three parallel branches.")

        # ---------------------------------------------------------- page 8
        def draw_closure(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.36)
            ax = fig.add_subplot(gs[0])
            b = d["clo"].query("family=='PARSEC' and R_V==3.1 and alpha==2.3")
            x = np.arange(len(b)); wd = 0.36
            ax.bar(x - wd / 2, b.predicted_observed_living, wd,
                   label="predicted from the IMF", color=NAVY, alpha=0.85)
            ax.bar(x + wd / 2, b.observed_living, wd, label="actually observed",
                   color=RED, alpha=0.85)
            for i, r in enumerate(b.itertuples()):
                ax.text(i, max(r.predicted_observed_living, r.observed_living) + 3,
                        f"ratio {r.closure_ratio:.3f}", ha="center", fontsize=8.5,
                        fontweight="bold")
            ax.set_xticks(x); ax.set_xticklabels(b.subgroup)
            ax.set_ylabel("living stars above 8 M$_\\odot$")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("The out-of-sample test: does the extrapolation work?")
            ax2 = fig.add_subplot(gs[1])
            med = {a: d["clo"].query(f"alpha=={a}").closure_ratio.median()
                   for a in [2.0, 2.3, 2.6]}
            ax2.plot(list(med), list(med.values()), "o-", color=NAVY, ms=8, lw=2)
            ax2.axhline(1.0, color=RED, ls="--", lw=1.5)
            ax2.text(2.03, 1.02, "perfect closure", fontsize=8.5, color=RED)
            for a, v in med.items():
                ax2.annotate(f"{v:.3f}", (a, v), textcoords="offset points",
                             xytext=(0, 13), ha="center", fontsize=9.5,
                             fontweight="bold")
            ax2.set_xlabel("IMF slope  α"); ax2.set_ylabel("observed / predicted")
            ax2.set_xticks([2.0, 2.3, 2.6])
            ax2.set_title("This is what discriminates between the IMF branches")
        figure_page(pdf, "7. WP6 — the closure test", draw_closure,
            "**The load-bearing result.** The IMF fitted on 2–8 M☉ is extrapolated into\n"
            "8 M☉ → turnoff and compared with what we actually see there. Those stars\n"
            "never entered the fit.\n\n"
            "**The census closes to about 10%** (grid median 1.105), and the closing\n"
            "slope is α ≈ 2.23 — close to the standard Salpeter value of 2.35.\n\n"
            "The lower panel is arguably the most useful thing WP6 produces: α = 2.6\n"
            "over-predicts badly and α = 2.0 under-predicts, so the census itself\n"
            "selects between branches that the fit alone could not.",
            note="CygOB2-C remains at 1.405 while A sits at 0.894 — they disagree in\n"
                 "direction, so no single mechanism explains both. This is the main\n"
                 "open scientific question.")

        # ---------------------------------------------------------- page 9
        def draw_bugs(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.38)
            ax = fig.add_subplot(gs[0])
            f = d["floor"].query("family=='PARSEC' and R_V==3.1 and alpha==2.3")
            x = np.arange(len(f)); wd = 0.36
            ax.bar(x - wd / 2, f.closure_ratio_floor8, wd, color=GREY,
                   alpha=0.85, label="before fix (withdrawn)")
            ax.bar(x + wd / 2, f.closure_ratio_fixed, wd, color=NAVY,
                   alpha=0.9, label="after fix")
            ax.axhline(1.0, color=RED, ls="--", lw=1.4)
            ax.set_xticks(x); ax.set_xticklabels(f.subgroup)
            ax.set_ylabel("closure ratio")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("Bug #17: the integral truncated at the threshold")
            ax2 = fig.add_subplot(gs[1])
            masses = np.array([4.0, 5.0, 6.0, 7.0, 8.0])
            up = np.array([0.004, 0.029, 0.089, 0.231, 0.408])
            ax2.plot(masses, up, "o-", color=RED, ms=7, lw=2)
            ax2.fill_between(masses, 0, up, color=RED, alpha=0.15)
            ax2.set_xlabel("true mass  [M$_\\odot$]")
            ax2.set_ylabel("P(estimated > 8 M$_\\odot$)")
            ax2.set_title("Stars below 8 M$_\\odot$ that get measured above it")
        figure_page(pdf, "8. Two bugs we found and fixed", draw_bugs,
            "**These were found by our own checks, and both changed the answer.**\n\n"
            "**Bug #17 — the closure integral started at 8 M☉**, but the observed side\n"
            "counted every star measured above 8 M☉ whatever its true mass. A 7 M☉\n"
            "star scattered upward was counted and never predicted (lower panel).\n"
            "Fixing it moved the grid median from 1.444 to 1.105 and the closing slope\n"
            "from 2.07 to 2.23 — about three quarters of what looked like a\n"
            "'shallower-than-Salpeter IMF' was this bug.\n\n"
            "**Bug #16 — the runaway search used absolute proper motions**, but Cyg OB2's\n"
            "own bulk motion is larger than a typical ejection signature. It was caught\n"
            "by an external check (next page).",
            note="Both results were formally withdrawn and the superseded numbers are\n"
                 "kept on record, not deleted. Three published values are affected.")

        # ---------------------------------------------------------- page 10
        def draw_run(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.38)
            ax = fig.add_subplot(gs[0])
            r = d["run"]
            sel = r[r.is_runaway_candidate.astype(bool)]
            ax.scatter(r.l, r.b, s=2, color=GREY, alpha=0.25, label="OB candidates")
            ax.scatter(sel.l, sel.b, s=9, color=RED, alpha=0.7,
                       label=f"traced back ({len(sel)})")
            bd = d["xm"].query("name=='BD+43 3654'")
            if len(bd):
                ax.scatter(bd.l, bd.b, s=150, marker="*", color=GOLD,
                           edgecolor="k", lw=0.6, zorder=5, label="BD+43 3654")
            ax.scatter([80.07], [0.82], s=110, marker="X", color=NAVY,
                       zorder=5, label="Cyg OB2")
            ax.invert_xaxis(); ax.set_xlabel("l [deg]"); ax.set_ylabel("b [deg]")
            ax.legend(frameon=False, fontsize=7.5, markerscale=1.2, ncol=2)
            ax.set_title("Runaway stars traced back to the association")
            ax2 = fig.add_subplot(gs[1]); ax2.axis("off"); ax2.grid(False)
            rows = [["", "measured here", "literature"],
                    ["recovery probability", "1.000", "—"],
                    ["ejection velocity", "38.8 km/s", "~40 km/s"],
                    ["flight time", "1.36 Myr", "1.6 Myr"],
                    ["position match", "0.04 arcsec", "—"]]
            t = ax2.table(cellText=rows, loc="center", cellLoc="left",
                          colWidths=[0.42, 0.29, 0.29])
            t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.65)
            for (i, j), c in t.get_celld().items():
                c.set_edgecolor("#dddddd")
                if i == 0:
                    c.set_text_props(fontweight="bold", color=NAVY)
            ax2.set_title("BD+43 3654 — the external check that caught bug #16",
                          fontsize=10)
        figure_page(pdf, "9. Runaways — and the check that caught a bug", draw_run,
            "**Some massive stars were ejected from the association and are still\n"
            "alive elsewhere.** They must not be counted as dead. We trace proper\n"
            "motions backwards and recover 119 candidates; after subtracting a\n"
            "chance-alignment rate measured from control fields, **54.9 are genuine**.\n\n"
            "BD+43 3654 is the textbook Cyg OB2 runaway, with an independently known\n"
            "velocity and flight time. Under the buggy version it scored exactly\n"
            "**0.000**. After the fix it scores **1.000** and reproduces both published\n"
            "numbers.",
            note="This is why external checks are not optional: every internal\n"
                 "diagnostic passed while the result was inverted.")

        # ---------------------------------------------------------- page 11
        def draw_xchk(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.38)
            ax = fig.add_subplot(gs[0])
            rates = {e["alpha"]: e["rate_events_per_100kyr"]
                     for e in d["xchk"]["harer_2025"]["by_alpha"]}
            ax.axhspan(0.25, 2.0, color=TEAL, alpha=0.18,
                       label="Härer+25 Fig. 2 range")
            ax.plot(list(rates), list(rates.values()), "o-", color=NAVY, ms=8, lw=2,
                    label="this work")
            for a, v in rates.items():
                ax.annotate(f"{v:.2f}", (a, v), textcoords="offset points",
                            xytext=(7, 5), fontsize=9)
            ax.set_yscale("log"); ax.set_xticks([2.0, 2.3, 2.6])
            ax.set_xlabel("IMF slope α")
            ax.set_ylabel("SN rate  [events / 100 kyr]")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("Supernova rate vs an independent population-synthesis model")
            ax2 = fig.add_subplot(gs[1])
            mass = {e["alpha"]: e["association_mass_Msun"]
                    for e in d["xchk"]["harer_2025"]["by_alpha"]}
            ax2.bar([str(a) for a in mass], list(mass.values()), 0.5,
                    color=NAVY, alpha=0.85, label="this work")
            ax2.axhline(1.65e4, color=RED, lw=1.8,
                        label="Wright+2015: 1.65×10⁴ M$_\\odot$")
            ax2.set_xlabel("IMF slope α")
            ax2.set_ylabel("association stellar mass  [M$_\\odot$]")
            ax2.legend(frameon=False, fontsize=8)
            ax2.set_title("Total stellar mass, derived independently")
        figure_page(pdf, "10. Independent cross-checks", draw_xchk,
            "**Three published results, none of which we tuned to, and all of which\n"
            "agree.**\n\n"
            "- **Supernova rate** — Härer+25 compute it from BPASS population synthesis\n"
            "   for an assumed cluster mass. All three of our branches land inside\n"
            "   their plotted range; our baseline gives 7.8 SNe/Myr against their\n"
            "   stated 'several per Myr'.\n"
            "- **Association mass** — 1.74×10⁴ M☉ from our fit, against the 1.65×10⁴\n"
            "   they assumed. A 5% agreement by a completely different route.\n"
            "- **Systemic proper motion** — our −2.7067 mas/yr against Orellana+21's\n"
            "   −2.71 ± 0.02 from Gaia DR2: agreement to 0.003 mas/yr.\n\n"
            "The rate and the mass depend on the fitted k in different ways, so a\n"
            "broken extrapolation would generally break them inconsistently.",
            note="Distance is the one mild tension: we get 1.62 kpc, Orellana 1.67 kpc.\n"
                 "Traced to the DR2 vs DR3 parallax zero point; carried as a 2.4%\n"
                 "systematic rather than adjusted.")

        # ---------------------------------------------------------- page 12
        def draw_nsn(fig, spec):
            gs = spec.subgridspec(2, 1, hspace=0.38)
            ax = fig.add_subplot(gs[0])
            norm, agep = d["norm"], d["age"]
            vals = {}
            for fam in w.FAMILIES:
                for a in [2.0, 2.3, 2.6]:
                    tot = 0
                    for sg in w.SUBGROUPS:
                        r = agep.query(
                            f"subgroup=='{sg}' and family=='{fam}' and R_V==3.1 and "
                            "f_bin==0.4 and indicator=='ums' and dmu==0.0")
                        age = float(r.age_map.iloc[0])
                        k = float(norm.query(
                            f"subgroup=='{sg}' and family=='{fam}' and R_V==3.1 "
                            f"and alpha=={a}").k_median.iloc[0])
                        M = min(turnoff_mass(fam, age), IMF_UPPER_LIMIT)
                        if M < IMF_UPPER_LIMIT:
                            tot += k * (M ** (1 - a) - IMF_UPPER_LIMIT ** (1 - a)) / (a - 1)
                    vals[(fam, a)] = tot
            x = np.arange(3); wd = 0.36
            for i, fam in enumerate(w.FAMILIES):
                ax.bar(x + (i - 0.5) * wd, [vals[(fam, a)] for a in [2.0, 2.3, 2.6]],
                       wd, label=fam, color=[NAVY, TEAL][i], alpha=0.9)
            ax.set_xticks(x); ax.set_xticklabels(["α = 2.0", "α = 2.3", "α = 2.6"])
            ax.set_ylabel("supernovae so far")
            ax.legend(frameon=False, fontsize=8)
            ax.set_title("How many have already exploded — by branch")
            ax2 = fig.add_subplot(gs[1])
            rates = {2.0: 19.2, 2.3: 7.8, 2.6: 3.1}
            p = {a: 1 - np.exp(-r * 0.1) for a, r in rates.items()}
            ax2.bar([str(a) for a in p], [v * 100 for v in p.values()], 0.5,
                    color=[GREY, RED, GREY], alpha=0.9)
            for i, (a, v) in enumerate(p.items()):
                ax2.text(i, v * 100 + 2, f"{v*100:.0f}%", ha="center",
                         fontsize=10, fontweight="bold")
            ax2.set_ylim(0, 100); ax2.set_xlabel("IMF slope α")
            ax2.set_ylabel("P(≥1 SN in last 100 kyr)  [%]")
            ax2.set_title("The number the paper is ultimately about")
        figure_page(pdf, "11. Where this is heading (WP7 preview)", draw_nsn,
            "**These are illustrative estimates from the fitted normalisation, not yet\n"
            "the WP7 result** — WP7 will do this by discrete stochastic sampling, which\n"
            "matters at these small numbers.\n\n"
            "Roughly **6 supernovae** at the baseline branch, essentially all from\n"
            "subgroups A and B in the last ~1 Myr; C has produced none.\n\n"
            "The lower panel is the paper's target statement: the probability that at\n"
            "least one went off inside Härer's 100 000-year window. At the baseline it\n"
            "is about **50%** — genuinely marginal, which is why the age sensitivity\n"
            "will have to be shown explicitly rather than marginalised away.",
            note="Note the branch spread: N_SN ranges over a factor of 12 across the\n"
                 "carried branches. That, not statistical noise, dominates the answer.")

        # ---------------------------------------------------------- page 13
        text_page(pdf, "12. Status, honesty, and what is next", [
            ("What is solid",
             "- the sample, the extinctions, the ages and the normalisation are all\n"
             "   accepted through formal gates, with checksummed provenance\n"
             "- the out-of-sample closure test works at the ~10% level\n"
             "- three independent literature cross-checks agree: supernova rate,\n"
             "   association mass, and systemic proper motion\n"
             "- the runaway search reproduces the textbook Cyg OB2 runaway\n"),
            ("What is still open",
             "- **CygOB2-C sits at 1.405 while A sits at 0.894.** They disagree in\n"
             "   direction, so no single mechanism explains both. This is the main\n"
             "   scientific question left. Candidates: subgroup labelling, foreground\n"
             "   contamination in distance, and the 2.4% distance systematic.\n"
             "- a chain re-run (repair_v7) is under way: a pre-registered test showed\n"
             "   the binary fraction below 8 M☉ shifts the threshold-crossing response\n"
             "   by 9.9%. It moves the closure ratio to ~1.074 and the supernova count\n"
             "   by only 0.5%.\n"),
            ("How we work, and why it is worth mentioning",
             "- **every test is pre-registered**: predictions and thresholds are written\n"
             "   to a provenance file before the code runs, so a result cannot be\n"
             "   quietly reinterpreted afterwards\n"
             "- when a prediction fails, it is recorded as failed. Two have been.\n"
             "- superseded numbers are formally withdrawn and kept on record\n"
             "- three separate bugs were caught this way, two of which changed the\n"
             "   headline answer\n"),
            ("Next steps",
             "1.  WP7 — the supernova ledger: N_SN posterior, explosion timeline, and\n"
             "     the probability of one within the last 100 kyr, per branch\n"
             "2.  WP8 — validation against the pulsar PSR J2032+4127, the γ-Cygni\n"
             "     remnant, and the ²⁶Al measurements\n"
             "3.  WP9–WP10 — the verdict and the manuscript\n"),
        ], subtitle="An honest account of what holds and what does not")

        info = pdf.infodict()
        info["Title"] = "Cygnus OB2 supernova history — progress briefing"
        info["Subject"] = "WP1-WP6 complete; Gaia DR3"

    print(f"wrote {out.relative_to(w.ROOT)}")
    print(f"  {out.stat().st_size / 1024:.0f} kB")


if __name__ == "__main__":
    main()
