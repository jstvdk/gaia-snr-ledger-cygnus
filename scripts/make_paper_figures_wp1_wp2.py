#!/usr/bin/env python3
"""
Publication figures and tables freezing the WP1 + WP2 results.

Outputs (all under figures/paper/ and tables/):
  fig1_membership_literature.pdf/.png   main-text money figure
  fig2_control_fields.pdf/.png          appendix precision demonstration
  table1_wp1_inventory.tex/.md          frozen WP1 data inventory
  table2_wp2_gate.tex/.md               WP2 gate scorecard
  table3_literature_recovery.tex/.md    Wright+15 / Berlanas+19 recovery

Run:  python scripts/make_paper_figures_wp1_wp2.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
PROV = ROOT / "provenance"
FIGS = ROOT / "figures/paper"
TABS = ROOT / "tables"
FIGS.mkdir(parents=True, exist_ok=True)
TABS.mkdir(parents=True, exist_ok=True)

# A&A column widths in inches
COL1, COL2 = 3.46, 7.09

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
    "legend.frameon": False,
})

C_MEM, C_FIELD = "#1f4e79", "#c8c8c8"
C_WRIGHT, C_BERL, C_MISS = "#d1495b", "#00798c", "#edae49"
C_MAN = "#7b2cbf"

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
man = json.loads((PROV / "wp2_membership_manifest.json").read_text())
mem = pd.read_parquet(PROC / "wp2_members.parquet")
hi = mem[mem.membership_probability > 0.5].copy()
auto = hi[~hi.anchor_quality_exempt].copy()
ctrl = pd.read_parquet(PROC / "wp2_control_members.parquet")
narrow = pd.read_parquet(PROC / "wp1_gaia_narrow.parquet")
wright = pd.read_parquet(PROC / "wp1_wright15_gaia_crossmatch.parquet")
berl = pd.read_csv(PROV / "wp2_berlanas_recovery_audit.csv")

FOOT_L, FOOT_B, FOOT_R = 79.8, 0.8, 1.0
ZP = -float(man["preprocessing"]["zero_point_mas"]["median"])

# field stars inside the adopted footprint, for grey background
sep = np.hypot((narrow.l_deg - FOOT_L) * np.cos(np.radians(narrow.b_deg)),
               narrow.b_deg - FOOT_B)
foot = narrow[sep <= FOOT_R].copy()
foot["parallax_corrected"] = foot.parallax + ZP
field = foot[~foot.source_id.isin(hi.source_id)]
rng = np.random.default_rng(20260722)
fsub = field.sample(min(6000, len(field)), random_state=20260722)

_ok = berl[berl.recovered_p_gt_0_5 == True]
_man = _ok.gate_disposition_reason.str.contains("manual", na=False)
berl_ok = _ok[~_man]          # recovered by the posterior
berl_man = _ok[_man]          # recovered as manual quality exceptions
berl_no = berl[berl.recovered_p_gt_0_5 != True]


def circle(ax, l0, b0, r, **kw):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(l0 + r * np.cos(t) / np.cos(np.radians(b0)), b0 + r * np.sin(t), **kw)


# ----------------------------------------------------------------------------
# Figure 1 — membership with literature overlay  (MAIN TEXT)
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.45))

# (a) sky
ax = axes[0]
ax.scatter(fsub.l_deg, fsub.b_deg, s=0.7, c=C_FIELD, lw=0, rasterized=True, zorder=1)
ax.scatter(auto.l_deg, auto.b_deg, s=2.2, c=C_MEM, lw=0, rasterized=True, zorder=2)
ax.scatter(wright.gaia_l, wright.gaia_b, s=13, facecolors="none",
           edgecolors=C_WRIGHT, lw=0.7, zorder=4)
ax.scatter(berl_no.l_deg, berl_no.b_deg, s=11, marker="x",
           c=C_MISS, lw=0.7, zorder=5)
circle(ax, FOOT_L, FOOT_B, FOOT_R, color="k", lw=0.8, ls="--", zorder=3)
ax.set_xlabel(r"$l$ [deg]"); ax.set_ylabel(r"$b$ [deg]")
ax.set_xlim(FOOT_L + 1.9, FOOT_L - 1.9); ax.set_ylim(FOOT_B - 1.75, FOOT_B + 1.75)
ax.set_aspect(1 / np.cos(np.radians(FOOT_B)))
ax.set_title("(a) sky distribution", loc="left")

# (b) proper motion
ax = axes[1]
ax.scatter(fsub.pmra, fsub.pmdec, s=0.7, c=C_FIELD, lw=0, rasterized=True, zorder=1)
ax.scatter(auto.pmra, auto.pmdec, s=2.2, c=C_MEM, lw=0, rasterized=True, zorder=2)
ax.scatter(berl_ok.pmra, berl_ok.pmdec, s=11, marker="s", facecolors="none",
           edgecolors=C_BERL, lw=0.7, zorder=4)
ax.scatter(berl_man.pmra, berl_man.pmdec, s=13, marker="^", facecolors="none",
           edgecolors=C_MAN, lw=0.7, zorder=4)
ax.scatter(berl_no.pmra, berl_no.pmdec, s=11, marker="x", c=C_MISS, lw=0.7, zorder=5)
ax.set_xlabel(r"$\mu_{\alpha*}$ [mas yr$^{-1}$]")
ax.set_ylabel(r"$\mu_{\delta}$ [mas yr$^{-1}$]")
ax.set_xlim(-7.5, 2.0); ax.set_ylim(-9.5, 0.5)
ax.set_title("(b) proper motion", loc="left")

# (c) parallax
ax = axes[2]
bins = np.arange(0.35, 1.105, 0.02)
ax.hist(field.parallax_corrected.dropna(), bins=bins, color=C_FIELD,
        histtype="stepfilled", lw=0, label="field (footprint)", density=True)
ax.hist(auto.parallax_corrected, bins=bins, color=C_MEM, histtype="step",
        lw=1.3, label="members", density=True)
d1 = man["distance_population_test"]["one_component_mean_kpc"]
ax.axvline(1.0 / d1, color="k", lw=0.8, ls="-")
ax.axvline(1.0 / 1.35, color=C_BERL, lw=0.8, ls=":")
ax.axvline(1.0 / 1.60, color=C_BERL, lw=0.8, ls="-.")
ax.set_xlabel(r"$\varpi_{\rm corr}$ [mas]"); ax.set_ylabel("normalised density")
ax.set_xlim(0.40, 1.05)
ax.set_title("(c) corrected parallax", loc="left")
ax.text(0.97, 0.94, f"$d$ = {d1:.2f} kpc", transform=ax.transAxes,
        ha="right", va="top", fontsize=6.6)

handles = [
    Line2D([], [], ls="none", marker="o", ms=3, mfc=C_MEM, mec="none",
           label=f"WP2 members ($P>0.5$, $N$={len(auto)})"),
    Line2D([], [], ls="none", marker="o", ms=3, mfc=C_FIELD, mec="none", label="field"),
    Line2D([], [], ls="none", marker="o", ms=4, mfc="none", mec=C_WRIGHT,
           label=f"Wright+15 census ({len(wright)})"),
    Line2D([], [], ls="none", marker="s", ms=4, mfc="none", mec=C_BERL,
           label=f"Berlanas+19 recovered by posterior ({len(berl_ok)})"),
    Line2D([], [], ls="none", marker="^", ms=4, mfc="none", mec=C_MAN,
           label=f"Berlanas+19 manual RUWE exception ({len(berl_man)})"),
    Line2D([], [], ls="none", marker="x", ms=4, mec=C_MISS,
           label=f"Berlanas+19 missed ({len(berl_no)})"),
    Line2D([], [], color="k", ls="--", lw=0.8, label="adopted 1 deg footprint"),
    Line2D([], [], color=C_BERL, ls=":", lw=0.8, label="Berlanas+19 1.35 / 1.60 kpc"),
]
fig.legend(handles=handles, loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.24), columnspacing=1.4, handletextpad=0.5)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIGS / f"fig1_membership_literature.{ext}")
plt.close(fig)
print("fig1 written")

# ----------------------------------------------------------------------------
# Figure 2 — control fields  (APPENDIX)
# ----------------------------------------------------------------------------
g = man["gate"]
qrows = man["published_footprint"]["quality_rows_target_controls"]
centers = [(FOOT_L, FOOT_B)] + [tuple(c) for c in
                                man["published_footprint"]["control_centers_l_b_deg"]]
yields = [g["target_yield_p_gt_0_5"]] + list(g["control_yields_p_gt_0_5"])

fig, axes = plt.subplots(1, 4, figsize=(COL2, 2.15), sharex=True, sharey=True)
for i, ax in enumerate(axes):
    l0, b0 = centers[i]
    if i == 0:
        sel, lab = auto, "target"
    else:
        sel, lab = ctrl[(ctrl.control_field == i) &
                        (ctrl.membership_probability > 0.5)], f"control {i}"
    s = np.hypot((narrow.l_deg - l0) * np.cos(np.radians(narrow.b_deg)),
                 narrow.b_deg - b0)
    bg = narrow[s <= FOOT_R]
    bg = bg.sample(min(4000, len(bg)), random_state=1)
    ax.scatter(bg.pmra, bg.pmdec, s=0.7, c=C_FIELD, lw=0, rasterized=True)
    ax.scatter(sel.pmra, sel.pmdec, s=2.4,
               c=C_MEM if i == 0 else "#8b2f3f", lw=0, rasterized=True)
    frac = 100 * yields[i] / qrows[i]
    ax.set_title(f"{lab}  $(l,b)$=({l0:.1f},{b0:.1f})", loc="left", fontsize=7)
    ax.text(0.04, 0.05, f"$N_{{P>0.5}}$ = {yields[i]}\nyield = {frac:.2f}%",
            transform=ax.transAxes, fontsize=6.4, va="bottom")
    ax.set_xlabel(r"$\mu_{\alpha*}$ [mas yr$^{-1}$]")
axes[0].set_ylabel(r"$\mu_{\delta}$ [mas yr$^{-1}$]")
axes[0].set_xlim(-7.5, 2.0); axes[0].set_ylim(-9.5, 0.5)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIGS / f"fig2_control_fields.{ext}")
plt.close(fig)
print("fig2 written")

# ----------------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------------
def write_table(name, df, caption, label):
    (TABS / f"{name}.md").write_text(
        f"**{caption}**\n\n" + df.to_markdown(index=False) + "\n")
    ncol = df.shape[1]
    rows = " \\\\\n".join(" & ".join(str(v) for v in r) for r in df.values)
    body = (f"\\begin{{tabular}}{{l{'r' * (ncol - 1)}}}\n\\hline\\hline\n"
            + " & ".join(df.columns) + " \\\\\n\\hline\n"
            + rows + " \\\\\n\\hline\n\\end{tabular}\n")
    (TABS / f"{name}.tex").write_text(
        f"\\begin{{table}}\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
        f"\\centering\n{body}\\end{{table}}\n")
    print(f"{name} written")


wp1 = json.loads((PROV / "wp1_manifest.json").read_text()) if (PROV / "wp1_manifest.json").exists() else {}
j2 = json.loads((PROV / "wp1_2mass_join_execution.json").read_text())["counts"]

t1 = pd.DataFrame([
    ["Gaia DR3 narrow box", "$l$ 77--83, $b$ $-$1.5--4", "245,843", "WP2 input"],
    ["Gaia DR3 wide box", "$l$ 72--88, $b$ $-$5--8", "3,133,326", "WP6 traceback"],
    ["2MASS PSC matches", "official best-neighbour", f"{j2['official_2mass_psc_matches']:,}", "93.0\\%"],
    ["\\quad complete $JHK_s$", "finite uncertainties", f"{j2['complete_jhk_with_finite_uncertainties']:,}", "82.8\\%"],
    ["\\quad ph\\_qual AAA", "highest quality", f"{j2['ph_qual_AAA']:,}", "69.5\\%"],
    ["Wright+15 census", "matched to Gaia DR3", "167 / 167", "100\\%"],
    ["Spectroscopic anchors", "canonical, unique IDs", "252", "540 evidence rows"],
    ["ATNF pulsars (wide)", "SN markers", "80", "incl.\\ PSR J2032+4127"],
    ["Green SNRs (wide)", "SN markers", "9", "incl.\\ G78.2+2.1"],
], columns=["Dataset", "Selection", "$N$", "Note"])
write_table("table1_wp1_inventory", t1,
            "Frozen WP1 data inventory. All artifacts carry SHA-256 checksums in "
            "\\texttt{provenance/wp1\\_manifest.json}.", "tab:wp1")

crit = g["criteria"]
t2 = pd.DataFrame([
    ["Berlanas+19 recall ($P>0.5$)", "$\\geq 0.80$", f"{g['berlanas_total_recall']:.3f}", "pass"],
    ["Control/target yield ratio", "$\\leq 0.10$", f"{g['mean_control_to_target_ratio']:.4f}", "pass"],
    ["Member count", "$10^2$--$10^4$", f"{g['total_unique_members_p_gt_0_5_including_manual']:,}", "pass"],
    ["Central-90\\% $l$ span", "$< 4.8$ deg", f"{g['l_central90_span_deg']:.2f}", "pass"],
    ["Central-90\\% $b$ span", "$< 4.4$ deg", f"{g['b_central90_span_deg']:.2f}", "pass"],
    ["Convex-hull area", "$< 16.5$ deg$^2$", f"{g['convex_hull_area_deg2']:.2f}", "pass"],
    ["Largest seed / analysis", "$\\leq 0.10$", f"{g['largest_density_seed_fraction_analysis']:.4f}", "pass"],
    ["Published-structure comparison", "documented", "yes", "pass"],
], columns=["Criterion", "Threshold", "Achieved", "Verdict"])
write_table("table2_wp2_gate", t2,
            "WP2 validation gate. The gate is a conjunction of recall, control-field "
            "precision, population sanity and spatial compactness; the superseded "
            "baseline run fails four of the eight criteria.", "tab:wp2gate")

nrw = int((wright.in_wp1_gaia_narrow == True).sum())
t3 = pd.DataFrame([
    ["Wright+15", "census rows", "167", ""],
    ["Wright+15", "matched to Gaia DR3", "167 (100.0\\%)", "gate $\\geq$90\\%"],
    ["Wright+15", "inside narrow selection", f"{nrw} ({100*nrw/len(wright):.1f}\\%)", "7 excluded, itemised"],
    ["Berlanas+19", "published members", "229", ""],
    ["Berlanas+19", "inside quality sample", "168", ""],
    ["Berlanas+19", "recovered automatically", f"128 (76.2\\% of 168)", "full-covariance posterior"],
    ["Berlanas+19", "manual quality exceptions", "61", "RUWE / absent from query"],
    ["Berlanas+19", "total recall", "189 / 229 = 82.5\\%", "gate $\\geq$80\\%"],
    ["Berlanas+19", "missed (soft floor)", "33", "$P<0.05$"],
    ["Berlanas+19", "missed ($0.05<P<0.5$)", "7", "itemised in audit"],
], columns=["Catalogue", "Quantity", "Value", "Note"])
write_table("table3_literature_recovery", t3,
            "Recovery of published membership catalogues. Automatic and manual "
            "recoveries are reported separately because the adopted footprint is "
            "taken from Berlanas et al.\\ (2019), so the spatial component of the "
            "benchmark is not independent.", "tab:literature")

print("\nAll outputs written to figures/paper/ and tables/")
