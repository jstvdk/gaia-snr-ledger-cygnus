#!/usr/bin/env python3
"""repair_v9 sensitivity S2 -- subgroup-assignment uncertainty (brief §4.5;
pre-registered in provenance/repair_v9_s1_s2_prereg.json).

The A/B/C labels are hard, but the WP2 Gaussian mixture behind them is soft.
S2 replays the frozen WP2 procedure exactly (wp2_derive_subgroups: clean
members, StandardScaler on l, b, pmra, pmdec, GaussianMixture k = 3 full
covariance, the 50 frozen seeds), names each seed's components with the frozen
physical rule (name_components), and averages each star's predict_proba over
the 50 seeds -> responsibilities r_A, r_B, r_C.

Propagation (weights, not resampling): every labelled star appears once per
subgroup with weight membership_probability x r_subgroup.  That is exactly how
the WP5 fit already counts stars (sum of membership_probability x soft bin
occupancy, wp5_joint_age_fit._observed_counts), so the fit runs unchanged on a
"soft-label" mass table.  Each star keeps its own repair_v9 mass posterior:
re-deriving it at another subgroup's age would move stars by at most the A-C
headline age difference (0.4 Myr), which matters near the turnoff and not in
WP5's calibration window (< 8 Msun).  Stars outside the clean WP2 set keep
their label; unassigned stars stay unassigned.

Stages:
  responsibilities  tables/repair_v9_s2_responsibilities.csv
  masses            data/processed/wp4_mass_posteriors_repair_v9_s2soft.parquet (+ samples)
  fit               wp5_fit_imf_joint.py on it, repair_v9 ages and responses
  ledger            WP7 engine (unmodified) on the S2 draws, 54 branches x 3
                    subgroups, 500,000 iterations -> tables/repair_v9_s2_label_uncertainty.csv

Run:  PYTHONPATH=scripts python3 scripts/repair_v9_s2_label_uncertainty.py <stage>
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import wp5_common as w
import wp7_ledger as L
from wp2_derive_subgroups import FEATURES, SEEDS, fit_labels, load_clean, name_components
from wp7_ledger_prereg import RECENT_WINDOW_MYR, SF_DURATIONS_MYR

VERSION = "repair_v9_s2soft"
LABELS = ["CygOB2-A", "CygOB2-B", "CygOB2-C"]
ITERATIONS = 500_000
RESP = w.TABLES / "repair_v9_s2_responsibilities.csv"
OUT = w.TABLES / "repair_v9_s2_label_uncertainty.csv"
OUT_JSON = w.PROVENANCE / "repair_v9_s2_execution.json"
PREREG = w.PROVENANCE / "repair_v9_s1_s2_prereg.json"


def stage_responsibilities() -> None:
    if RESP.exists():
        raise SystemExit(f"{RESP.relative_to(w.ROOT)} exists")
    _, clean = load_clean()
    side = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    hard = clean[["source_id"]].merge(side[["source_id", "subgroup"]], on="source_id",
                                      how="left", validate="one_to_one")
    X = StandardScaler().fit_transform(clean[FEATURES].values)
    acc = np.zeros((len(clean), 3))
    for seed in SEEDS:
        gm = fit_labels(X, 3, seed)
        naming = name_components(clean, gm.predict(X))
        proba = gm.predict_proba(X)
        for comp in range(3):
            acc[:, LABELS.index(naming[comp])] += proba[:, comp]
    acc /= len(SEEDS)
    out = pd.DataFrame({"source_id": clean["source_id"].to_numpy("int64"),
                        "subgroup_hard": hard["subgroup"].to_numpy(),
                        **{f"r_{lab[-1]}": acc[:, i] for i, lab in enumerate(LABELS)}})
    out["r_max"] = acc.max(axis=1)
    out["argmax_matches_hard"] = [LABELS[i] == h for i, h in zip(acc.argmax(axis=1), out.subgroup_hard)]
    out.to_csv(RESP, index=False)
    print(f"wrote {RESP.name}: {len(out)} stars, median r_max {out.r_max.median():.3f}, "
          f"r_max>0.9 {np.mean(out.r_max > 0.9):.3f}, argmax=hard {out.argmax_matches_hard.mean():.3f}")


def stage_masses() -> None:
    target = w.PROC / f"wp4_mass_posteriors_{VERSION}.parquet"
    if target.exists():
        raise SystemExit(f"{target.name} exists")
    masses = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v9.parquet")
    store = np.load(w.PROC / "wp4_mass_posterior_samples_repair_v9.npz")
    resp = pd.read_csv(RESP).set_index("source_id")
    soft = masses.source_id.isin(resp.index) & masses.subgroup.isin(LABELS)
    parts, idx = [masses[~soft]], [np.flatnonzero(~soft)]
    for lab in LABELS:
        block = masses[soft].copy()
        r = resp.loc[block.source_id, f"r_{lab[-1]}"].to_numpy(float)
        block["subgroup"] = lab
        block["membership_probability"] = block["membership_probability"].to_numpy(float) * r
        parts.append(block)
        idx.append(np.flatnonzero(soft))
    out = pd.concat(parts, ignore_index=True)
    order = np.concatenate(idx)
    out.to_parquet(target, index=False)
    np.savez_compressed(w.PROC / f"wp4_mass_posterior_samples_{VERSION}.npz",
                        source_id=out["source_id"].to_numpy("int64"), family=store["family"],
                        rv=store["rv"], samples=store["samples"][order])
    for lab in LABELS:
        hard = masses.loc[masses.subgroup.eq(lab), "membership_probability"].sum()
        softw = out.loc[out.subgroup.eq(lab), "membership_probability"].sum()
        print(f"{lab}: membership weight hard {hard:.1f} -> soft {softw:.1f}")


def stage_fit() -> None:
    subprocess.run([sys.executable, "scripts/wp5_fit_imf_joint.py", "--upstream-version", "repair_v5",
                    "--mass-version", VERSION, "--age-version", "repair_v9_headline",
                    "--response-version", "repair_v9", "--wp5-version", VERSION,
                    "--compare-version", "repair_v9"], cwd=w.ROOT, check=True)


def stage_ledger() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists")
    rows = []
    for version in ("repair_v9", VERSION):   # same engine, same seeds: a paired comparison
        draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{version}.npz")
        relations = {f: L.TurnoffRelation(f) for f in w.FAMILIES}
        for fi, fam in enumerate(w.FAMILIES):
            for ri, rv in enumerate(w.R_V_BRANCHES):
                for ai, alpha in enumerate(w.IMF_SLOPES):
                    for di, delta in enumerate(SF_DURATIONS_MYR):
                        total = np.zeros(ITERATIONS, dtype=int)
                        last = np.full(ITERATIONS, np.inf)
                        for si, sg in enumerate(LABELS):
                            key = L.draw_key(sg, fam, rv, alpha)
                            k_all, age_all = draws[f"k__{key}"], draws[f"truth_age_draws__{key}"]
                            rng = np.random.default_rng([w.SEED, 2, fi, ri, ai, di, si])
                            pick = rng.integers(0, k_all.size, ITERATIONS)
                            res = L.run_population(rng, k_all[pick], age_all[pick], alpha, delta,
                                                   relations[fam])
                            n, t = res["n_sn"]["all_explode"], res["t_last"]["all_explode"]
                            total += n
                            last = np.minimum(last, t)
                            rows.append({"labels": "hard" if version == "repair_v9" else "soft",
                                         "family": fam, "R_V": rv, "alpha": alpha,
                                         "sf_duration_Myr": delta, "subgroup": sg,
                                         "k_median": float(np.median(k_all)),
                                         "N_death_mean": float(n.mean()),
                                         "P_last_within_100kyr": float((t < RECENT_WINDOW_MYR).mean())})
                        rows.append({"labels": "hard" if version == "repair_v9" else "soft",
                                     "family": fam, "R_V": rv, "alpha": alpha,
                                     "sf_duration_Myr": delta, "subgroup": "ALL", "k_median": np.nan,
                                     "N_death_mean": float(total.mean()),
                                     "P_last_within_100kyr": float((last < RECENT_WINDOW_MYR).mean())})
        print(f"{version} done", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    resp = pd.read_csv(RESP)
    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_s2_label_uncertainty.py",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "responsibilities": {"stars": int(len(resp)), "median_r_max": float(resp.r_max.median()),
                             "fraction_r_max_gt_0p9": float(np.mean(resp.r_max > 0.9)),
                             "fraction_argmax_equals_hard_label": float(resp.argmax_matches_hard.mean())},
        "iterations": ITERATIONS,
        "inputs": {p: w.sha256(w.ROOT / p) for p in (
            "data/processed/wp5_imf_posterior_draws_repair_v9.npz",
            f"data/processed/wp5_imf_posterior_draws_{VERSION}.npz",
            "data/processed/wp4_mass_posteriors_repair_v9.parquet",
            "tables/wp2_subgroup_labels.parquet", "data/processed/wp2_members.parquet")},
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (RESP, OUT)},
    })
    base = table[table.family.eq("PARSEC") & np.isclose(table.R_V, 3.1) & np.isclose(table.alpha, 2.3)
                 & np.isclose(table.sf_duration_Myr, 0.0)]
    print(base.pivot(index="subgroup", columns="labels", values="N_death_mean").round(3))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["responsibilities", "masses", "fit", "ledger"])
    stage = ap.parse_args().stage
    {"responsibilities": stage_responsibilities, "masses": stage_masses,
     "fit": stage_fit, "ledger": stage_ledger}[stage]()


if __name__ == "__main__":
    main()
