#!/usr/bin/env python3
"""WP13 -- T7 ledger-level injections and the D6 (age gap x young fraction)
grid (prereg A3, A4).  Runs after the M0 fit: sigma_0 is M0's measured width.

Baseline cell only (PARSEC, R_V 3.1, alpha 2.3, delta 0), all-explode, the
unmodified WP7 engine.  Per realisation every model gets a fitted posterior
centred on a noisy draw of its truth with the measured width:
    t_hat ~ N(t_true, sigma),  ages ~ N(t_hat, sigma)            (clipped 1-10 Myr)
    ln k_hat ~ N(ln k_true, s_k),  ln k ~ N(ln k_hat, s_k)
and is run through wp7_ledger.run_population (T7_ITERATIONS draws).  The truth
is mu_true = sum_s k_s int_{M_to(t_s)}^{120} M^-alpha, and P_true comes from the
engine at the point-mass truth (200,000 iterations).

  T7a  coeval: all subgroups at 3.516 Myr.  False pass = wp13_rules.t2 or t3
       passes; the same realisations give the null distribution of T4's D.
  T7b  multi-age: true ages = the M1 WP5 posterior means; M0's fitted truth is
       the k-weighted mean true age.  Bias of the posterior-mean mu and of P,
       and M1's 68 % coverage of mu_true.
  D6   descriptive map, see wp13_rules.D6_*.

Outputs: tables/wp13_t7.csv (per realisation), tables/wp13_d6_grid.csv,
provenance/wp13_t7_execution.json.

Run:  PYTHONPATH=scripts python3 scripts/wp13_t7.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp13_rules as R
import wp5_common as w
import wp7_ledger as L
from wp13_score import LABELLED, M0_DRAWS, M1_DRAWS, PREREG, population, prereg

OUT = w.TABLES / "wp13_t7.csv"
OUT_D6 = w.TABLES / "wp13_d6_grid.csv"
OUT_JSON = w.PROVENANCE / "wp13_t7_execution.json"
T_COEVAL = 3.516                 # the disclosed k-weighted mean age (brief §2, re-based)
TRUTH_ITERATIONS = 200_000
FAM, RV, ALPHA, DELTA = (R.BASELINE[k] for k in ("family", "R_V", "alpha", "sf_duration_Myr"))


def widths(draws, label: str) -> dict:
    key = L.draw_key(label, FAM, RV, ALPHA)
    k, age = draws[f"k__{key}"], draws[f"truth_age_draws__{key}"]
    return {"k_median": float(np.median(k)), "s_k": float(np.std(np.log(k))),
            "age_mean": float(np.mean(age)), "sigma": float(np.std(age))}


def fitted(rng, t_true, sigma, k_true, s_k, n):
    t_hat = t_true + sigma * rng.standard_normal()
    lk_hat = np.log(k_true) + s_k * rng.standard_normal() if k_true > 0 else -np.inf
    ages = np.clip(t_hat + sigma * rng.standard_normal(n), 1.0, 10.0)
    ks = np.exp(lk_hat + s_k * rng.standard_normal(n)) if k_true > 0 else np.zeros(n)
    return ks, ages


def run_model(rng, comps, relation, n):
    """comps: list of (t_true, sigma, k_true, s_k); returns the summed population."""
    parts = [population(rng, *fitted(rng, t, s, k, sk, n), ALPHA, DELTA, relation)
             for t, s, k, sk in comps]
    return {"mu": sum(p["mu"] for p in parts),
            "last": np.minimum.reduce([p["last"] for p in parts]),
            "epochs": np.concatenate([p["epochs"] for p in parts])}


def truth(rng, comps, relation):
    mu = sum(k * float(L.imf_integral(ALPHA, relation.turnoff(np.array([t + DELTA / 2])), L.IMF_UPPER_LIMIT)[0])
             for t, k in comps)
    last = np.full(TRUTH_ITERATIONS, np.inf)
    for t, k in comps:
        res = L.run_population(rng, np.full(TRUTH_ITERATIONS, k), np.full(TRUTH_ITERATIONS, t),
                               ALPHA, DELTA, relation)
        last = np.minimum(last, res["t_last"]["all_explode"])
    return mu, float((last < R.RECENT_WINDOW_MYR).mean())


def main() -> None:
    prereg()
    for p in (OUT, OUT_D6, OUT_JSON):
        if p.exists():
            raise SystemExit(f"{p.relative_to(w.ROOT)} exists")
    relation = L.TurnoffRelation(FAM)
    m1, m0 = np.load(M1_DRAWS), np.load(M0_DRAWS)
    wid = {sg: widths(m1, sg) for sg in LABELLED}
    w0 = widths(m0, R.POOLED_LABEL)
    k_true = {sg: wid[sg]["k_median"] for sg in LABELLED}
    k_tot = sum(k_true.values())
    n = R.T7_ITERATIONS
    rows = []

    designs = {"T7a": {sg: T_COEVAL for sg in LABELLED},
               "T7b": {sg: wid[sg]["age_mean"] for sg in LABELLED}}
    summary = {}
    for di, (name, t_true) in enumerate(designs.items()):
        t0 = sum(k_true[sg] * t_true[sg] for sg in LABELLED) / k_tot
        mu_true, p_true = truth(np.random.default_rng([w.SEED, 13, 7, di, 10 ** 6]),
                                [(t_true[sg], k_true[sg]) for sg in LABELLED], relation)
        for r in range(R.T7_REALISATIONS):
            rng = np.random.default_rng([w.SEED, 13, 7, di, r])
            one = run_model(rng, [(t_true[sg], wid[sg]["sigma"], k_true[sg], wid[sg]["s_k"])
                                  for sg in LABELLED], relation, n)
            zero = run_model(rng, [(t0, w0["sigma"], k_tot, w0["s_k"])], relation, n)
            i1, i0 = one["last"] < R.RECENT_WINDOW_MYR, zero["last"] < R.RECENT_WINDOW_MYR
            q16, q84 = np.percentile(one["mu"], [16, 84])
            rows.append({"design": name, "realisation": r, "mu_true": mu_true, "P_true": p_true,
                         "N_hat_M1": float(one["mu"].mean()), "N_hat_M0": float(zero["mu"].mean()),
                         "P_hat_M1": float(i1.mean()), "P_hat_M0": float(i0.mean()),
                         "M1_covers_68": bool(q16 <= mu_true <= q84),
                         "T2_pass": R.t2(zero["mu"], one["mu"])["pass"],
                         "T3_pass": R.t3(i0, i1)["pass"],
                         "D": R.sup_distance(zero["epochs"], one["epochs"])})
        d = pd.DataFrame([x for x in rows if x["design"] == name])
        summary[name] = {
            "true_ages_Myr": t_true, "M0_truth_age_Myr": t0, "mu_true": mu_true, "P_true": p_true,
            "false_pass_rate": float((d.T2_pass | d.T3_pass).mean()),
            "D_null95": float(np.percentile(d.D, R.T4_NULL_QUANTILE)),
            "bias_N_M1": float((d.N_hat_M1 - mu_true).mean()), "bias_N_M0": float((d.N_hat_M0 - mu_true).mean()),
            "bias_P_M1": float((d.P_hat_M1 - p_true).mean()), "bias_P_M0": float((d.P_hat_M0 - p_true).mean()),
            "M1_coverage_68": float(d.M1_covers_68.mean())}
        print(name, {k: v for k, v in summary[name].items() if k != "true_ages_Myr"}, flush=True)

    a, b = summary["T7a"], summary["T7b"]
    summary["T7a"]["pass"] = bool(a["false_pass_rate"] < R.T7A_MAX_FALSE_PASS)
    summary["T7b"]["pass"] = bool(abs(b["bias_N_M1"]) < abs(b["bias_N_M0"])
                                  and abs(b["bias_P_M1"]) < abs(b["bias_P_M0"])
                                  and R.T7B_COVERAGE[0] <= b["M1_coverage_68"] <= R.T7B_COVERAGE[1])

    # ---- D6: descriptive grid
    sig_med = float(np.median([wid[sg]["sigma"] for sg in LABELLED]))
    sk_med = float(np.median([wid[sg]["s_k"] for sg in LABELLED]))
    d6 = []
    for gi, gap in enumerate(R.D6_AGE_GAP_MYR):
        for fi, f in enumerate(R.D6_YOUNG_FRACTION):
            comps = [(T_COEVAL, (1 - f) * k_tot), (T_COEVAL - gap, f * k_tot)]
            t0 = T_COEVAL - f * gap
            mu_true, p_true = truth(np.random.default_rng([w.SEED, 13, 6, gi, fi, 10 ** 6]), comps, relation)
            est = []
            for r in range(R.T7_REALISATIONS):
                rng = np.random.default_rng([w.SEED, 13, 6, gi, fi, r])
                one = run_model(rng, [(t, sig_med, k, sk_med) for t, k in comps], relation, n)
                zero = run_model(rng, [(t0, w0["sigma"], k_tot, w0["s_k"])], relation, n)
                est.append((one["mu"].mean(), zero["mu"].mean(),
                            (one["last"] < R.RECENT_WINDOW_MYR).mean(), (zero["last"] < R.RECENT_WINDOW_MYR).mean()))
            e = np.array(est)
            bias_n0 = float(e[:, 1].mean() - mu_true)
            bias_p0 = float(e[:, 3].mean() - p_true)
            d6.append({"age_gap_Myr": gap, "young_fraction": f, "mu_true": mu_true, "P_true": p_true,
                       "bias_N_M1_rel": float(e[:, 0].mean() - mu_true) / mu_true,
                       "bias_N_M0_rel": bias_n0 / mu_true,
                       "bias_P_M1": float(e[:, 2].mean() - p_true), "bias_P_M0": bias_p0,
                       "one_age_adequate": bool(abs(bias_n0) / mu_true <= R.D6_ADEQUATE_N_REL
                                                and abs(bias_p0) <= R.D6_ADEQUATE_P_ABS)})
        print(f"D6 gap {gap} done", flush=True)

    pd.DataFrame(rows).to_csv(OUT, index=False)
    pd.DataFrame(d6).to_csv(OUT_D6, index=False)
    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_t7.py",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "widths": {"M1": wid, "M0": w0, "D6_sigma_median": sig_med, "D6_s_k_median": sk_med},
        "iterations": {"per_model_per_realisation": n, "truth": TRUTH_ITERATIONS,
                       "realisations": R.T7_REALISATIONS},
        "T7a": summary["T7a"], "T7b": summary["T7b"],
        "D6_adequate_cells": int(sum(x["one_age_adequate"] for x in d6)),
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (M1_DRAWS, M0_DRAWS)},
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (OUT, OUT_D6)},
    })
    print(f"T7a {'PASS' if summary['T7a']['pass'] else 'FAIL'}  T7b {'PASS' if summary['T7b']['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
