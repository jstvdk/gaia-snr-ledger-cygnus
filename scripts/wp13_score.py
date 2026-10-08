#!/usr/bin/env python3
"""WP13 -- integrity checks, T1, the paired M0/M1 ledger and the verdict.

Every threshold and rule is imported from wp13_rules (hashed in
provenance/wp13_prereg.json); this file only feeds them.

Stages (each refuses to overwrite its outputs):
  I3        every consumed input and frozen-code hash matches the prereg
  I1        identity replay (pooling off) vs repair_v9: anchors, masses,
            27 PARSEC R_V 3.1 nodes, WP5 k_median, within repair_v9_integrity.TOL
  t1        ln BF(M1:M0) per family x R_V from the stored test-c curves (A1)
            -> tables/wp13_t1.csv
  ledger    M1 replayed with the stored WP7 seeding (W13-I2) and M0 / M0-lite on
            their own seeds, paired by iteration index, 54 branches, 2,000,000
            iterations -> tables/wp13_ablation.csv, data/processed/wp13_baseline_curves.npz
  score     T1-T7 and the verdict -> provenance/wp13_outcome.json
            (needs scripts/wp13_t7.py to have run)

Run:  PYTHONPATH=scripts python3 scripts/wp13_score.py <stage>
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import repair_v9_integrity as RI
import wp13_rules as R
import wp4v9_common as V
import wp5_common as w
import wp5_joint_age_fit as J
import wp7_ledger as L
from wp7_ledger_prereg import SF_DURATIONS_MYR

PREREG = w.PROVENANCE / "wp13_prereg.json"
INTEGRITY = w.PROVENANCE / "wp13_integrity.json"
OUTCOME = w.PROVENANCE / "wp13_outcome.json"
T1_OUT = w.TABLES / "wp13_t1.csv"
ABLATION = w.TABLES / "wp13_ablation.csv"
CURVES = w.PROC / "wp13_baseline_curves.npz"
T7_RECORD = w.PROVENANCE / "wp13_t7_execution.json"
M1_DRAWS = w.PROC / "wp5_imf_posterior_draws_repair_v9.npz"
M0_DRAWS = w.PROC / "wp5_imf_posterior_draws_wp13_m0.npz"
M1_LEDGER = w.TABLES / "wp7_ledger_repair_v9.csv"
M1_NORM = w.PROC / "wp5_imf_normalization_repair_v9.parquet"
LABELLED = ("CygOB2-A", "CygOB2-B", "CygOB2-C")
SN = L.SN_THRESHOLD_MSUN
EPOCH_GRID = np.linspace(0.0, 10.0, 1001)


def prereg() -> dict:
    if not PREREG.exists():
        raise SystemExit("provenance/wp13_prereg.json missing -- pre-register and commit first")
    return json.loads(PREREG.read_text())


def save_check(name: str, result: dict) -> None:
    out = json.loads(INTEGRITY.read_text()) if INTEGRITY.exists() else {"checks": {}}
    result["run_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out["checks"][name] = result
    out.update({"script": "scripts/wp13_score.py", "preregistration": str(PREREG.relative_to(w.ROOT)),
                "preregistration_sha256": w.sha256(PREREG)})
    w.write_json(INTEGRITY, out)
    print(f"{name}: {'PASS' if result['pass'] else 'FAIL'}")


# ----------------------------------------------------------------- I3
def stage_i3() -> None:
    rec = prereg()
    bad_in = {n: v["path"] for n, v in rec["consumed_inputs"].items()
              if w.sha256(w.ROOT / v["path"]) != v["sha256"]}
    bad_code = [p for p, h in rec["frozen_code"].items() if w.sha256(w.ROOT / p) != h]
    save_check("W13-I3", {"pass": not bad_in and not bad_code, "inputs_changed": bad_in,
                          "code_changed": bad_code})


# ----------------------------------------------------------------- I1
def stage_i1() -> None:
    prereg()
    TOL = RI.TOL
    parts = {}
    # anchors
    d = RI.numeric_diff(pd.read_parquet(w.PROC / "wp4_anchor_hrd_wp13_identity.parquet"),
                        pd.read_parquet(w.PROC / "wp4_anchor_hrd_repair_v9.parquet"), "source_id")
    cols = d.get("columns", {})
    worst = max((c["max_abs"] for c in cols.values()), default=0.0)
    nan_ok = all(c["nan_pattern_same"] for c in cols.values())
    parts["anchors"] = {"pass": bool(d["rows_match"] and nan_ok and worst <= TOL["I1a_anchor_abs"]),
                        "max_abs_any_column": worst}
    # masses
    a = pd.read_parquet(w.PROC / "wp4_mass_posteriors_wp13_identity.parquet")
    b = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v9.parquet")
    d = RI.numeric_diff(a, b, "source_id")
    cols = d.get("columns", {})
    mass_cols = [c for c in cols if c.startswith("mass_")]
    frac = max((cols[c]["frac_rel_gt_tol"] for c in mass_cols), default=0.0)
    gt8 = abs(float(a["mass_baseline_p_gt8"].sum()) - float(b["mass_baseline_p_gt8"].sum()))
    parts["masses"] = {"pass": bool(d["rows_match"] and all(cols[c]["nan_pattern_same"] for c in mass_cols)
                                    and frac <= TOL["I1b_mass_frac_beyond"]
                                    and gt8 <= TOL["I1b_expected_gt8_abs"]),
                       "max_frac_entries_beyond_rel_tol": frac, "expected_gt8_abs_diff": gt8}
    # nodes
    ages = pd.read_parquet(w.PROC / "wp4_age_posteriors_repair_v9_headline.parquet")
    native = J.native_isochrone_ages("PARSEC")
    node_ok, n_nodes = True, 0
    for sg in LABELLED:
        for age in J.truth_age_nodes(ages, sg, "PARSEC", 3.1, native, snap=False):
            ra = pd.read_parquet(J.node_response_path(sg, "PARSEC", 3.1, age, "wp13_identity"))
            rb = pd.read_parquet(J.node_response_path(sg, "PARSEC", 3.1, age, "repair_v9"))
            ca = pd.read_parquet(J.node_curve_path(sg, "PARSEC", 3.1, age, "wp13_identity"))
            cb = pd.read_parquet(J.node_curve_path(sg, "PARSEC", 3.1, age, "repair_v9"))
            draws = [c for c in rb.columns if c.startswith("recovered_mass_draw_")]
            x, y = ra[draws].to_numpy(float), rb[draws].to_numpy(float)
            ok = x.shape == y.shape
            if ok:
                both = np.isfinite(x) & np.isfinite(y)
                eq = (np.isnan(x) & np.isnan(y)) | (both & (np.abs(x - y) <= TOL["I1c_draw_rel"] * np.abs(y)))
                cdiff = float(np.nanmax(np.abs(
                    ca.sort_values("primary_mass").recovery_isotonic.to_numpy(float)
                    - cb.sort_values("primary_mass").recovery_isotonic.to_numpy(float))))
                ok = float(eq.mean()) >= TOL["I1c_draw_frac_equal"] and cdiff <= TOL["I1c_curve_abs"]
            node_ok &= ok
            n_nodes += 1
    parts["nodes"] = {"pass": bool(node_ok and n_nodes == 27), "nodes": n_nodes}
    # WP5 fit
    a = pd.read_parquet(w.PROC / "wp5_imf_normalization_wp13_identity.parquet")
    b = pd.read_parquet(M1_NORM)
    key = ["subgroup", "family", "R_V", "alpha"]
    m = a.merge(b, on=key, suffixes=("_a", "_b"), validate="one_to_one")
    rel = float(((m.k_median_a - m.k_median_b).abs() / m.k_median_b).max())
    flips = int((m.residual_gate_pass_a.astype(bool) != m.residual_gate_pass_b.astype(bool)).sum())
    parts["wp5"] = {"pass": bool(len(m) == 54 and rel <= TOL["I1d_k_rel"]
                                 and flips <= TOL["I1d_gate_flips_max"]),
                    "k_median_max_rel_diff": rel, "residual_gate_flips": flips}
    save_check("W13-I1", {"pass": all(p["pass"] for p in parts.values()), "parts": parts,
                          "tolerances": "repair_v9_integrity.TOL"})


# ----------------------------------------------------------------- T1
def stage_t1() -> None:
    prereg()
    if T1_OUT.exists():
        raise SystemExit(f"{T1_OUT.relative_to(w.ROOT)} exists")
    curves = json.loads(V.frozen("spectroscopic_curves").read_text())["loglike_curves"]
    native = {f: np.array(v) for f, v in json.loads(
        (w.PROVENANCE / "wp4v9_fit_execution.json").read_text())["native_ages"].items()}
    rows = []
    for fam in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            ca = np.array(curves[f"c|{R.EPS}|CygOB2-A|{fam}|{float(rv)}"], float)
            cc = np.array(curves[f"c|{R.EPS}|CygOB2-C|{fam}|{float(rv)}"], float)
            rows.append({"family": fam, "R_V": float(rv),
                         "lnZ_A": R.log_evidence(native[fam], ca),
                         "lnZ_C": R.log_evidence(native[fam], cc),
                         "lnZ_AC": R.log_evidence(native[fam], ca + cc),
                         "ln_BF_M1_M0": R.t1_ln_bf(native[fam], ca, cc)})
    table = pd.DataFrame(rows)
    table.to_csv(T1_OUT, index=False)
    print(table.round(3).to_string(index=False))


# ----------------------------------------------------------------- ledger
def mu_of(k, age, alpha, delta, relation):
    return k * L.imf_integral(alpha, relation.turnoff(age + delta / 2.0), L.IMF_UPPER_LIMIT)


def population(rng, k, age, alpha, delta, relation) -> dict:
    res = L.run_population(rng, k, age, alpha, delta, relation)
    mu = mu_of(k, age, alpha, delta, relation)
    if not np.isclose(mu.mean(), res["mean_expected"], rtol=1e-12):
        raise RuntimeError("expected-count recomputation disagrees with run_population")
    keep = res["dead_masses"] >= SN
    return {"n": res["n_sn"]["all_explode"], "n_islands": res["n_sn"]["islands"],
            "last": res["t_last"]["all_explode"], "mu": mu,
            "epochs": res["epochs"][keep], "iter": res["dead_iteration"][keep]}


def combine(parts: list[dict], n_iter: int) -> dict:
    out = {"n": sum(p["n"] for p in parts), "n_islands": sum(p["n_islands"] for p in parts),
           "mu": sum(p["mu"] for p in parts),
           "last": np.minimum.reduce([p["last"] for p in parts]),
           "epochs": np.concatenate([p["epochs"] for p in parts]),
           "iter": np.concatenate([p["iter"] for p in parts])}
    return out


def cdf(epochs: np.ndarray) -> np.ndarray:
    e = np.sort(epochs)
    return np.searchsorted(e, EPOCH_GRID, side="right") / max(e.size, 1)


def stage_ledger(n_iter: int) -> None:
    rec = prereg()
    if ABLATION.exists():
        raise SystemExit(f"{ABLATION.relative_to(w.ROOT)} exists")
    m1, m0 = np.load(M1_DRAWS), np.load(M0_DRAWS)
    stored = pd.read_csv(M1_LEDGER)
    stored = stored[stored.explodability.eq("all_explode")]
    relations = {f: L.TurnoffRelation(f) for f in w.FAMILIES}
    rng_master = np.random.default_rng(w.SEED)          # the stored WP7 seeding (W13-I2)
    rows, replay, baseline_extra = [], [], {}
    for fi, fam in enumerate(w.FAMILIES):
        rel = relations[fam]
        for ri, rv in enumerate(w.R_V_BRANCHES):
            for ai, alpha in enumerate(w.IMF_SLOPES):
                for di, delta in enumerate(SF_DURATIONS_MYR):
                    subs = {}
                    for sg in LABELLED:
                        key = L.draw_key(sg, fam, rv, alpha)
                        k_all, age_all = m1[f"k__{key}"], m1[f"truth_age_draws__{key}"]
                        rng = np.random.default_rng(rng_master.integers(0, 2 ** 63 - 1))
                        pick = rng.integers(0, k_all.size, n_iter)
                        subs[sg] = population(rng, k_all[pick], age_all[pick], alpha, delta, rel)
                    one = combine(list(subs.values()), n_iter)
                    # M0 and M0-lite: own seeds, same engine
                    key0 = L.draw_key(R.POOLED_LABEL, fam, rv, alpha)
                    k0, a0 = m0[f"k__{key0}"], m0[f"truth_age_draws__{key0}"]
                    rng = np.random.default_rng([w.SEED, 13, fi, ri, ai, di])
                    pick = rng.integers(0, k0.size, n_iter)
                    zero = population(rng, k0[pick], a0[pick], alpha, delta, rel)
                    ksum = sum(m1[f"k__{L.draw_key(sg, fam, rv, alpha)}"] for sg in LABELLED)
                    n_lite = min(ksum.size, a0.size)
                    rng = np.random.default_rng([w.SEED, 13, 99, fi, ri, ai, di])
                    pick = rng.integers(0, n_lite, n_iter)
                    lite = population(rng, ksum[pick], a0[pick], alpha, delta, rel)

                    t2 = R.t2(zero["mu"], one["mu"])
                    t3 = R.t3(zero["last"] < R.RECENT_WINDOW_MYR, one["last"] < R.RECENT_WINDOW_MYR)
                    first1 = R.first_death_epoch(one["epochs"], one["iter"], n_iter)
                    first0 = R.first_death_epoch(zero["epochs"], zero["iter"], n_iter)
                    D = R.sup_distance(zero["epochs"], one["epochs"])
                    row = {"family": fam, "R_V": float(rv), "alpha": float(alpha),
                           "sf_duration_Myr": float(delta),
                           "scored_T6": bool(np.isclose(alpha, R.SCORED_ALPHA)),
                           "reported_36": bool(any(np.isclose(alpha, a) for a in R.REPORTED_ALPHAS)),
                           "N_M1": float(one["n"].mean()), "N_M0": float(zero["n"].mean()),
                           "N_M0_lite": float(lite["n"].mean()),
                           "N_M1_islands": float(one["n_islands"].mean()),
                           "N_M0_islands": float(zero["n_islands"].mean()),
                           "mu_M1": float(one["mu"].mean()), "mu_M0": float(zero["mu"].mean()),
                           "P_M1": float((one["last"] < R.RECENT_WINDOW_MYR).mean()),
                           "P_M0": float((zero["last"] < R.RECENT_WINDOW_MYR).mean()),
                           "P_M0_lite": float((lite["last"] < R.RECENT_WINDOW_MYR).mean()),
                           **{k: v for k, v in t2.items() if k != "pass"}, "T2": t2["pass"],
                           **{k: v for k, v in t3.items() if k != "pass"}, "T3": t3["pass"],
                           "D": D, "first_death_M1_Myr": first1, "first_death_M0_Myr": first0}
                    rows.append(row)
                    # W13-I2: the replay against the stored ledger
                    for sg, part in [*subs.items(), ("ALL", one)]:
                        s = stored[stored.subgroup.eq(sg) & stored.family.eq(fam)
                                   & np.isclose(stored.R_V, rv) & np.isclose(stored.alpha, alpha)
                                   & np.isclose(stored.sf_duration_Myr, delta)]
                        if len(s) != 1:
                            raise RuntimeError(f"stored ledger row missing for {sg} {fam} {rv} {alpha} {delta}")
                        replay.append({"subgroup": sg, "family": fam, "R_V": rv, "alpha": alpha,
                                       "sf_duration_Myr": delta,
                                       "dN": abs(float(part["n"].mean()) - float(s.N_SN_mean.iloc[0])),
                                       "dP": abs(float((part["last"] < R.RECENT_WINDOW_MYR).mean())
                                                 - float(s.P_last_SN_within_100kyr.iloc[0]))})
                    if (fam, float(rv), float(alpha), float(delta)) == tuple(R.BASELINE.values()):
                        baseline_extra = {
                            "N_M1_sub": {sg: float(p["n"].mean()) for sg, p in subs.items()},
                            "cdf_M1": cdf(one["epochs"]), "cdf_M0": cdf(zero["epochs"]),
                            "dmu_sample": (zero["mu"] - one["mu"])[:200_000],
                        }
                    print(f"  {fam} R_V {rv} a {alpha} d {delta}: N_M1 {row['N_M1']:.3f} "
                          f"N_M0 {row['N_M0']:.3f} P {row['P_M1']:.3f}/{row['P_M0']:.3f}", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(ABLATION, index=False)
    np.savez_compressed(CURVES, epoch_grid=EPOCH_GRID, cdf_M1=baseline_extra["cdf_M1"],
                        cdf_M0=baseline_extra["cdf_M0"], dmu_sample=baseline_extra["dmu_sample"])
    rp = pd.DataFrame(replay)
    tol_n = RI.TOL["I1e_ndeath_abs"]
    i2 = {"pass": bool(rp.dN.max() <= tol_n and rp.dP.max() <= 0.005),
          "max_abs_dN": float(rp.dN.max()), "max_abs_dP": float(rp.dP.max()),
          "rows": int(len(rp)), "tolerance": {"N": tol_n, "P": 0.005}}
    save_check("W13-I2", i2)
    w.write_json(w.PROVENANCE / "wp13_ledger_execution.json", {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_score.py ledger", "iterations": n_iter,
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "seeds": {"M1": "stored WP7 seeding (default_rng(SEED), one child per branch x subgroup)",
                  "M0": "default_rng([SEED, 13, fi, ri, ai, di])",
                  "M0_lite": "default_rng([SEED, 13, 99, fi, ri, ai, di])"},
        "baseline_N_M1_by_subgroup": baseline_extra["N_M1_sub"],
        "held_fixed_W13_I5": {
            "labelled_members_pooled": 1331, "unlabelled_excluded_from_both": 61,
            "extinction": "repair_v5", "imf_upper_limit_Msun": L.IMF_UPPER_LIMIT,
            "imf_slopes": list(w.IMF_SLOPES), "sf_durations_Myr": list(SF_DURATIONS_MYR),
            "sn_threshold_Msun": SN, "turnoff_relation": "wp7_ledger.TurnoffRelation",
            "engine": "wp7_ledger.run_population, unmodified"},
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (M1_DRAWS, M0_DRAWS, M1_LEDGER)},
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (ABLATION, CURVES)},
        "prereg_iterations": rec["thresholds"]["LEDGER_ITERATIONS"],
    })


# ----------------------------------------------------------------- score
def stage_score() -> None:
    prereg()
    if OUTCOME.exists():
        raise SystemExit(f"{OUTCOME.relative_to(w.ROOT)} exists")
    integ = json.loads(INTEGRITY.read_text())["checks"]
    t1tab = pd.read_csv(T1_OUT)
    abl = pd.read_csv(ABLATION)
    t7 = json.loads(T7_RECORD.read_text())
    led = json.loads((w.PROVENANCE / "wp13_ledger_execution.json").read_text())

    t1 = R.t1_pass({(r.family, float(r.R_V)): float(r.ln_BF_M1_M0) for r in t1tab.itertuples()})
    b = R.BASELINE
    base = abl[abl.family.eq(b["family"]) & np.isclose(abl.R_V, b["R_V"]) & np.isclose(abl.alpha, b["alpha"])
               & np.isclose(abl.sf_duration_Myr, b["sf_duration_Myr"])].iloc[0]
    t4 = R.t4(float(base.D), float(t7["T7a"]["D_null95"]), float(base.first_death_M0_Myr),
              float(base.first_death_M1_Myr))
    norm = pd.read_parquet(M1_NORM)
    k_sub = {sg: float(norm[norm.subgroup.eq(sg) & norm.family.eq(b["family"]) & np.isclose(norm.R_V, b["R_V"])
                            & np.isclose(norm.alpha, b["alpha"])].k_median.iloc[0]) for sg in LABELLED}
    t5 = R.t5(led["baseline_N_M1_by_subgroup"], float(base.N_M0), k_sub)
    base_flags = {"dN": float(base.dN), "dP": float(base.dP), "t2": bool(base.T2), "t3": bool(base.T3)}

    def t6_on(sel):
        return R.t6([{"family": r.family, "dN": r.dN, "dP": r.dP, "t2": bool(r.T2), "t3": bool(r.T3)}
                     for r in sel.itertuples()], base_flags)
    t6 = t6_on(abl[abl.scored_T6])
    t6_36 = t6_on(abl[abl.reported_36])
    flags = {"T1": t1["pass"], "T2": bool(base.T2), "T3": bool(base.T3), "T4": t4["pass"],
             "T5": t5["pass"], "T6": t6["pass"], "T7a": t7["T7a"]["pass"], "T7b": t7["T7b"]["pass"]}
    outcome = R.verdict(flags)
    w.write_json(OUTCOME, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_score.py score",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "integrity": {k: v["pass"] for k, v in integ.items()},
        "verdict": outcome, "flags": flags,
        "T1": t1, "T2": {k: float(base[k]) for k in ("dN", "dN_lo95", "dN_hi95", "N_M1")},
        "T3": {k: float(base[k]) for k in ("dP", "dP_lo95", "dP_hi95")}, "T4": t4, "T5": t5,
        "T6_scored_18": t6, "T6_reported_36": t6_36, "T7": t7,
        "baseline": {k: (v.item() if hasattr(v, "item") else v) for k, v in base.to_dict().items()},
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p)
                   for p in (T1_OUT, ABLATION, T7_RECORD, INTEGRITY)},
    })
    print(f"verdict: {outcome}")
    print(json.dumps(flags))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["I3", "I1", "t1", "ledger", "score"])
    ap.add_argument("--iterations", type=int, default=R.LEDGER_ITERATIONS)
    args = ap.parse_args()
    if args.stage == "ledger" and args.iterations != R.LEDGER_ITERATIONS:
        sys.exit("the ledger runs at the pre-registered iteration count only")
    {"I3": stage_i3, "I1": stage_i1, "t1": stage_t1, "score": stage_score,
     "ledger": lambda: stage_ledger(args.iterations)}[args.stage]()


if __name__ == "__main__":
    main()
