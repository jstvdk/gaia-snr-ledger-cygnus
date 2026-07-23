#!/usr/bin/env python3
"""WP3 figures + spatial-coherence test + WP5 truncation-bias analysis.

Produces:
  figures/wp3/wp3_dereddened_cmd.png     de-reddened CMDs (both colours) + isochrones
  figures/wp3/wp3_extinction_map.png     sky map coloured by A_V (spatial coherence)
  figures/wp3/wp3_av_distribution.png    A_V histogram, anchor vs broadband, cube comparison
  figures/wp3/wp3_wp5_truncation.png     the G<19 extinction-completeness edge

Also computes the nearest-neighbour A_V coherence statistic (the plan's actual
'extinction correlates across neighbouring stars' requirement) and the A_V
threshold above which stars of a given mass drop below G=19 at 1.62 kpc.

Provenance: provenance/wp3_figures_coherence_execution.json
"""
import hashlib, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from wp3_common import DIST_MODULUS, D_KPC
from wp3_extinction_law import band_coefficients

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

PAR = "data/processed/wp3_isochrones_parsec.parquet"
MIST = "data/processed/wp3_isochrones_mist.parquet"

def neighbour_coherence(df, k=8):
    """Spearman correlation between each star's A_V and the mean A_V of its k
    nearest sky neighbours; a random-label permutation baseline for context."""
    m = df["av_rv3.1"].notna() & df["l_deg"].notna() & df["b_deg"].notna()
    d = df[m]
    xy = np.column_stack([d.l_deg.values * np.cos(np.radians(d.b_deg.values)), d.b_deg.values])
    av = d["av_rv3.1"].values
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=k + 1)
    neigh_mean = np.array([np.mean(av[idx[i, 1:]]) for i in range(len(av))])
    rho, p = spearmanr(av, neigh_mean)
    # permutation baseline
    rng = np.random.default_rng(3)
    perm_rhos = []
    for _ in range(200):
        avp = rng.permutation(av)
        nm = np.array([np.mean(avp[idx[i, 1:]]) for i in range(len(avp))])
        perm_rhos.append(spearmanr(avp, nm)[0])
    return {"k": k, "n": int(len(av)), "spearman_rho": float(rho), "spearman_p": float(p),
            "perm_rho_mean": float(np.mean(perm_rhos)), "perm_rho_std": float(np.std(perm_rhos)),
            "z_over_permutation": float((rho - np.mean(perm_rhos)) / np.std(perm_rhos))}

def truncation_edge(masses=(2.0, 3.0, 5.0), age_myr=4.0, Glim=19.0):
    """A_V at which a star of given mass reaches G=Glim at the group distance."""
    par = pd.read_parquet(PAR)
    iso = par[np.abs(par.age_Myr - age_myr) < 0.6]
    kG = band_coefficients(3.1)["G"]
    edges = {}
    for mtarget in masses:
        row = iso.iloc[(iso.Mini - mtarget).abs().argsort().iloc[0]]
        MG0 = row.G0
        # Glim = MG0 + mu + kG*A_V  ->  A_V = (Glim - MG0 - mu)/kG
        av_edge = (Glim - MG0 - DIST_MODULUS) / kG
        edges[f"{mtarget:.1f}"] = {"M_G0": float(MG0), "A_V_edge": float(av_edge)}
    return edges, kG

def main():
    import os
    os.makedirs("figures/wp3", exist_ok=True)
    df = pd.read_parquet("data/processed/wp3_extinction.parquet")
    cube = pd.read_parquet("data/processed/wp3_cube_comparison.parquet")
    par = pd.read_parquet(PAR); mist = pd.read_parquet(MIST)

    # ---------- Figure 1: de-reddened CMDs ----------
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, (cx, mx, xlab) in zip(axes, [("BPRP0_rv3.1", "G0_abs_rv3.1", r"$(BP-RP)_0$"),
                                          ("GKs0_rv3.1", "G0_abs_rv3.1", r"$(G-K_s)_0$")]):
        sc = ax.scatter(df[cx], df[mx], c=df["av_rv3.1"], s=8, cmap="inferno_r",
                        vmin=3, vmax=9, alpha=0.7, rasterized=True)
        for iso, col, lab in [(par, "tab:cyan", "PARSEC 4 Myr"), (mist, "tab:green", "MIST 4 Myr")]:
            g = iso[np.abs(iso.age_Myr - 4.0) < 0.3].copy()
            if "label" in g:
                g = g[g.label <= 3]        # PMS+MS+early post-MS
            g = g.sort_values("Mini")
            xcol = (g.BP0 - g.RP0) if "BPRP" in cx else (g.G0 - g.Ks0)
            ax.plot(xcol, g.G0, color=col, lw=1.5, label=lab)
        ax.set_xlabel(xlab); ax.set_ylabel(r"$M_{G,0}$"); ax.invert_yaxis()
        ax.legend(loc="upper right", fontsize=8); ax.grid(alpha=0.2)
    cb = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.02); cb.set_label(r"$A_V$ (mag), $R_V=3.1$")
    fig.suptitle(f"WP3 de-reddened CMDs — Cyg OB2 members (d={D_KPC} kpc, N={len(df)})")
    fig.savefig("figures/wp3/wp3_dereddened_cmd.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ---------- Figure 2: extinction sky map ----------
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sc = ax.scatter(df.l_deg, df.b_deg, c=df["av_rv3.1"], s=22, cmap="inferno_r",
                    vmin=3, vmax=9, edgecolor="k", linewidth=0.15)
    ax.set_xlabel("Galactic $l$ (deg)"); ax.set_ylabel("Galactic $b$ (deg)")
    ax.invert_xaxis(); ax.set_title("WP3 per-star extinction map (spatial coherence)")
    cb = fig.colorbar(sc); cb.set_label(r"$A_V$ (mag), $R_V=3.1$")
    fig.savefig("figures/wp3/wp3_extinction_map.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ---------- Figure 3: A_V distribution ----------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    anc = df[df.av_source.str.startswith("anchor")]; bb = df[df.av_source == "broadband"]
    ax.hist(bb["av_rv3.1"].dropna(), bins=40, range=(0, 14), alpha=0.6, label=f"broadband (N={len(bb)})", color="tab:blue")
    ax.hist(anc["av_rv3.1"].dropna(), bins=40, range=(0, 14), alpha=0.7, label=f"anchor intrinsic-colour (N={len(anc)})", color="tab:orange")
    ax.axvspan(4, 8, color="green", alpha=0.12, label="Cyg OB2 core 4–8 mag")
    ax.set_xlabel(r"$A_V$ (mag)"); ax.set_ylabel("N"); ax.legend(fontsize=8); ax.set_title("A_V distribution")
    ax = axes[1]
    for rv, c in [("3.0", "tab:red"), ("3.1", "k"), ("3.5", "tab:purple")]:
        ax.hist(df[f"av_rv{rv}"].dropna(), bins=40, range=(0, 14), histtype="step", lw=1.6, label=f"$R_V$={rv}", color=c)
    ax.set_xlabel(r"$A_V$ (mag)"); ax.set_ylabel("N"); ax.legend(); ax.set_title("R_V branches")
    ax = axes[2]
    m = cube.vergely_A0.notna()
    ax.scatter(cube.vergely_A0[m], cube["av_rv3.1"][m], s=6, alpha=0.4)
    lim = [0, 12]; ax.plot(lim, lim, "k--", lw=1, label="1:1")
    ax.set_xlabel(r"Vergely+22 $A_0$ (mag)"); ax.set_ylabel(r"per-star $A_V$ (mag)")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.legend(); ax.set_title("vs Vergely+22 cube")
    fig.savefig("figures/wp3/wp3_av_distribution.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ---------- coherence + truncation ----------
    coh = neighbour_coherence(df)
    edges, kG = truncation_edge()

    # ---------- Figure 4: WP5 truncation edge ----------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.hist(df["av_rv3.1"].dropna(), bins=40, range=(0, 14), color="0.7", label="member $A_V$")
    for m, col in zip(edges, ["tab:red", "tab:orange", "tab:green"]):
        ax.axvline(edges[m]["A_V_edge"], color=col, lw=2,
                   label=f"{m} $M_\\odot$ lost above $A_V$={edges[m]['A_V_edge']:.1f}")
    ax.set_xlabel(r"$A_V$ (mag)"); ax.set_ylabel("N members")
    ax.set_title("WP5 preview: G<19 extinction-completeness edge at 1.62 kpc")
    ax.legend(fontsize=9)
    fig.savefig("figures/wp3/wp3_wp5_truncation.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    prov = {
        "script": "scripts/wp3_figures_and_coherence.py",
        "inputs": {"wp3_extinction": sha("data/processed/wp3_extinction.parquet"),
                   "wp3_cube_comparison": sha("data/processed/wp3_cube_comparison.parquet")},
        "neighbour_coherence": coh,
        "wp5_truncation_edges_Glt19_at_1.62kpc": edges,
        "kG_over_AV": kG,
        "figures": [f"figures/wp3/{f}" for f in
                    ["wp3_dereddened_cmd.png", "wp3_extinction_map.png",
                     "wp3_av_distribution.png", "wp3_wp5_truncation.png"]],
    }
    json.dump(prov, open("provenance/wp3_figures_coherence_execution.json", "w"), indent=2)
    print("Nearest-neighbour A_V coherence: rho=%.3f (p=%.1e), permutation baseline %.3f±%.3f -> z=%.1f"
          % (coh["spearman_rho"], coh["spearman_p"], coh["perm_rho_mean"], coh["perm_rho_std"], coh["z_over_permutation"]))
    print("WP5 truncation edges (A_V where G=19 at 1.62 kpc):")
    for m, e in edges.items():
        print(f"  {m} Msun: M_G0={e['M_G0']:+.2f} -> A_V_edge={e['A_V_edge']:.1f} mag")
    print("Figures written to figures/wp3/")

if __name__ == "__main__":
    main()
