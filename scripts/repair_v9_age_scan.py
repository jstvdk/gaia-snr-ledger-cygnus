#!/usr/bin/env python3
"""repair_v9 -- the age scan (owner decision 7 of 2026-10-07; prereg §3).

Every subgroup is scanned, independently, over the eleven native isochrone
ages from 2.0 to 6.3 Myr of each family:
  PARSEC 2.00 2.24 2.51 2.82 3.16 3.55 3.98 4.47 5.01 5.62 6.31
  MIST   2.00 2.25 2.52 2.83 3.18 3.57 4.01 4.50 5.05 5.67 6.37
At scan point i all three subgroups are set to the i-th native age of their
family.  Each point is a mini-chain registered in chain.py as
repair_v9_scan<ii>:
  tables   the age table: every subgroup a point mass at t (map = 68/90 % edges = t)
  anchors  wp4_anchors_hrd.py at t                     (run by the runner)
  masses   wp4_mass_posteriors_repair.py at t           (run by the runner)
  inject   repair_v9_injections.py, ONE native node per subgroup x family x R_V
           (snapped rule, no interpolation, never a reused G2 snapshot)
  fit      wp5_fit_imf_joint.py -- k is refitted at t
  ledger   the WP7 engine (wp7_ledger.run_population, unmodified) per subgroup and
           branch, 200,000 iterations, all-explode; the FULL distribution of
           N_death is stored (as a probability mass function), so the deaths of
           any (t_A, t_B, t_C) follow by convolving independent subgroups -- the
           subgroup populations are independent Poisson draws
  aggregate tables/repair_v9_age_scan.csv, data/processed/repair_v9_age_scan_draws.npz,
           figures/repair_v9/repair_v9_age_scan.png

The scan is a sensitivity, never the headline.  It runs no WP6 and no WP8-12.

Run:  PYTHONPATH=scripts python3 scripts/repair_v9_age_scan.py <stage> [--index i]
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import chain as C
import wp5_common as w
import wp5_joint_age_fit as J
from wp7_ledger import TurnoffRelation, draw_key, run_population
from wp7_ledger_prereg import RECENT_WINDOW_MYR, SF_DURATIONS_MYR

SCAN_ITERATIONS = 200_000
HEADLINE_AGES = w.PROC / "wp4_age_posteriors_repair_v9_headline.parquet"
OUT_CSV = w.TABLES / "repair_v9_age_scan.csv"
OUT_NPZ = w.PROC / "repair_v9_age_scan_draws.npz"
OUT_FIG = w.ROOT / "figures" / "repair_v9" / "repair_v9_age_scan.png"
OUT_JSON = w.PROVENANCE / "repair_v9_age_scan_execution.json"
PARTS = w.PROC / "repair_v9_age_scan_parts"


def scan_version(i: int) -> str:
    return f"repair_v9_scan{i:02d}"


def scan_ages() -> dict[str, np.ndarray]:
    out = {}
    for fam in w.FAMILIES:
        a = J.native_isochrone_ages(fam)
        a = a[(a >= 1.99) & (a <= 6.4)]
        if len(a) != C.SCAN_POINTS:
            raise RuntimeError(f"{fam}: {len(a)} native ages in 2.0-6.4 Myr, expected {C.SCAN_POINTS}")
        out[fam] = a
    return out


def stage_tables() -> None:
    head = pd.read_parquet(HEADLINE_AGES)
    ages = scan_ages()
    for i in range(C.SCAN_POINTS):
        path = w.PROC / f"wp4_age_posteriors_{scan_version(i)}.parquet"
        if path.exists():
            print(f"present: {path.name}")
            continue
        t = head.copy()
        t_age = t.family.map({f: float(ages[f][i]) for f in w.FAMILIES})
        for col in ("age_map", "age_lo68", "age_hi68", "age_lo90", "age_hi90", "age_mean"):
            t[col] = t_age
        t["window"] = "age_scan_point"
        t["error_model"] = "age_scan"
        t["fallback"] = ""
        t["grid_railed"] = False
        t["measurable"] = True
        t["exclusion_reason"] = ""
        t.to_parquet(path, index=False)
        print(f"wrote {path.name}: PARSEC {ages['PARSEC'][i]:.3f} / MIST {ages['MIST'][i]:.3f} Myr")


def stage_check_nodes(i: int) -> None:
    """The WP5 fit silently falls back to a gate-G2 scan snapshot when a snapped
    node file is missing; refuse to fit unless every scan node exists."""
    version = scan_version(i)
    ages = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{version}.parquet")
    native = {f: J.native_isochrone_ages(f) for f in w.FAMILIES}
    missing = []
    for fam in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for sg in w.SUBGROUPS:
                nodes = J.truth_age_nodes(ages, sg, fam, rv, native[fam], snap=True)
                if len(nodes) != 1:
                    raise RuntimeError(f"{version} {sg}/{fam}/{rv}: {len(nodes)} nodes, expected 1")
                for age in nodes:
                    for p in (J.node_response_path(sg, fam, rv, age, version),
                              J.node_curve_path(sg, fam, rv, age, version)):
                        if not p.exists():
                            missing.append(p.name)
    if missing:
        raise SystemExit(f"{version}: {len(missing)} node files missing, e.g. {missing[:3]}")
    print(f"{version}: all 18 single-node responses present")


def stage_ledger(i: int) -> None:
    version = scan_version(i)
    PARTS.mkdir(exist_ok=True)
    part = PARTS / f"{version}.npz"
    if part.exists():
        print(f"present: {part.name}")
        return
    draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{version}.npz")
    norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{version}.parquet")
    relations = {f: TurnoffRelation(f) for f in w.FAMILIES}
    rows, pmfs = [], {}
    for fi, fam in enumerate(w.FAMILIES):
        for ri, rv in enumerate(w.R_V_BRANCHES):
            for ai, alpha in enumerate(w.IMF_SLOPES):
                for di, delta in enumerate(SF_DURATIONS_MYR):
                    for si, sg in enumerate(w.SUBGROUPS):
                        key = draw_key(sg, fam, rv, alpha)
                        k_all, age_all = draws[f"k__{key}"], draws[f"truth_age_draws__{key}"]
                        rng = np.random.default_rng([w.SEED, i, fi, ri, ai, di, si])
                        pick = rng.integers(0, k_all.size, SCAN_ITERATIONS)
                        res = run_population(rng, k_all[pick], age_all[pick], alpha, delta,
                                             relations[fam])
                        n = res["n_sn"]["all_explode"]
                        last = res["t_last"]["all_explode"]
                        cell = norm[norm.subgroup.eq(sg) & norm.family.eq(fam)
                                    & np.isclose(norm.R_V, rv) & np.isclose(norm.alpha, alpha)]
                        pmf = np.bincount(n) / n.size
                        pmfs[f"{fam}|{rv}|{alpha}|{delta}|{sg}"] = pmf
                        rows.append({
                            "scan_index": i, "subgroup": sg, "family": fam, "R_V": rv,
                            "alpha": alpha, "sf_duration_Myr": delta,
                            "age_Myr": float(np.median(age_all)),
                            "k_median": float(np.median(k_all)),
                            "residual_gate_pass": (bool(cell.residual_gate_pass.iloc[0])
                                                   if len(cell) == 1 else None),
                            "N_death_mean": float(n.mean()),
                            "N_death_median": float(np.median(n)),
                            "N_death_p16": float(np.percentile(n, 16)),
                            "N_death_p84": float(np.percentile(n, 84)),
                            "P_at_least_one": float((n >= 1).mean()),
                            "P_last_within_100kyr": float((last < RECENT_WINDOW_MYR).mean()),
                            "iterations": SCAN_ITERATIONS,
                        })
    keys = sorted(pmfs)
    tmp = part.with_name(part.name + ".tmp.npz")
    np.savez_compressed(tmp, keys=np.array(keys), rows=json.dumps(rows),
                        **{f"pmf_{j}": pmfs[k] for j, k in enumerate(keys)})
    tmp.replace(part)
    print(f"wrote {part.name}")


def stage_aggregate() -> None:
    if OUT_CSV.exists():
        raise SystemExit(f"{OUT_CSV.relative_to(w.ROOT)} exists; nothing is overwritten")
    rows, store = [], {}
    for i in range(C.SCAN_POINTS):
        z = np.load(PARTS / f"{scan_version(i)}.npz")
        rows += json.loads(str(z["rows"]))
        for j, k in enumerate(z["keys"]):
            store[f"{k}|{i}"] = z[f"pmf_{j}"]
    table = pd.DataFrame(rows).sort_values(
        ["family", "R_V", "alpha", "sf_duration_Myr", "subgroup", "scan_index"])
    table.to_csv(OUT_CSV, index=False)
    keys = sorted(store)
    np.savez_compressed(OUT_NPZ, keys=np.array(keys),
                        note=np.array("pmf of N_death (all-explode) per family|R_V|alpha|"
                                      "sf_duration|subgroup|scan_index; 200,000 iterations"),
                        **{f"pmf_{j}": store[k] for j, k in enumerate(keys)})
    figure(table)
    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_age_scan.py",
        "preregistration": "provenance/repair_v9_prereg.json",
        "scan_ages_Myr": {f: [float(a) for a in v] for f, v in scan_ages().items()},
        "iterations": SCAN_ITERATIONS,
        "status": "sensitivity -- never the headline",
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in sorted(PARTS.glob("*.npz"))},
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (OUT_CSV, OUT_NPZ, OUT_FIG)},
    })
    print(f"wrote {OUT_CSV.relative_to(w.ROOT)} ({len(table)} rows), {OUT_NPZ.name}, {OUT_FIG.name}")


# reference categorical palette, slots 1-3 (fixed order: A, B, C)
COLOURS = {"CygOB2-A": "#2a78d6", "CygOB2-B": "#eb6834", "CygOB2-C": "#1baf7a"}


def figure(table: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    head = pd.read_parquet(HEADLINE_AGES)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True, constrained_layout=True)
    for ax, fam in zip(axes, w.FAMILIES):
        sel = table[table.family.eq(fam) & np.isclose(table.R_V, 3.1)
                    & np.isclose(table.alpha, 2.3) & np.isclose(table.sf_duration_Myr, 0.0)]
        for sg in w.SUBGROUPS:
            s = sel[sel.subgroup.eq(sg)].sort_values("age_Myr")
            c = COLOURS[sg]
            ax.fill_between(s.age_Myr, s.N_death_p16, s.N_death_p84, color=c, alpha=0.15, lw=0)
            ax.plot(s.age_Myr, s.N_death_mean, color=c, lw=2, marker="o", ms=4,
                    label=f"{sg[-1]}")
            ax.annotate(sg[-1], (s.age_Myr.iloc[-1], s.N_death_mean.iloc[-1]),
                        xytext=(4, 0), textcoords="offset points", va="center",
                        color="#333333", fontsize=9)
            h = head[head.subgroup.eq(sg) & head.family.eq(fam) & np.isclose(head.R_V, 3.1)
                     & np.isclose(head.f_bin, 0.4) & np.isclose(head.dmu, 0.0)]
            ax.axvspan(float(h.age_lo68.iloc[0]), float(h.age_hi68.iloc[0]), ymin=0,
                       ymax=0.025 + 0.025 * w.SUBGROUPS.index(sg), color=c, alpha=0.5, lw=0)
        ax.set_title(f"{fam}, R_V 3.1, α 2.3, coeval", fontsize=10, color="#333333")
        ax.set_xlabel("assumed age of the subgroup (Myr)")
        ax.grid(alpha=0.25, lw=0.6)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    axes[0].set_ylabel("N_death | all explode (mean, 16–84 %)")
    axes[0].legend(title="subgroup", frameon=False, loc="upper left")
    fig.suptitle("repair_v9 age scan — sensitivity, not the headline (k refitted at each age; "
                 "bars at the foot: headline 68 % age intervals)", fontsize=9, color="#555555")
    fig.savefig(OUT_FIG, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["tables", "check-nodes", "ledger", "aggregate"])
    ap.add_argument("--index", type=int)
    args = ap.parse_args()
    if args.stage == "tables":
        stage_tables()
    elif args.stage == "aggregate":
        stage_aggregate()
    else:
        if args.index is None:
            raise SystemExit("--index is required")
        {"check-nodes": stage_check_nodes, "ledger": stage_ledger}[args.stage](args.index)


if __name__ == "__main__":
    main()
