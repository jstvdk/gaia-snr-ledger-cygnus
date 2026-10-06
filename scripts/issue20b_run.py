#!/usr/bin/env python3
"""Issue #20 Phase A' -- run tests a, b, c for A, B and C (provenance/issue20b_prereg.json).

Per family x R_V x subgroup: the near-IR free-extinction CMD fit (test a, plus
its WP3-extinction control), the ESP-HS hot-star HRD (test b) and the
outlier-robust spectroscopic HRD (test c), each with its registered validity
check; the calibrations of tests a and b; the guarded two-age mixture of test
c (report only); integrity checks J1 (synthetic recovery) and J2 (test c with
eps = 0 reproduces Phase A).  Per-test log-likelihood curves are written so
that the scorer can form the adopted mixture posteriors.

Outputs:
  tables/issue20b_age_tests.csv
  tables/issue20b_calibrations.csv
  provenance/issue20b_run_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20b_run.py
"""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I20
import issue20_hrd_likelihood as A3
import issue20b_common as B
import wp5_common as w
from wp3_extinction_law import band_coefficients

OUT_TESTS = w.TABLES / "issue20b_age_tests.csv"
OUT_CAL = w.TABLES / "issue20b_calibrations.csv"
OUT_EXEC = w.PROVENANCE / "issue20b_run_execution.json"
EPS_VALUES = (B.EPS_OUTLIER, 0.02, 0.10)


def guarded_mixture(joint, marg, weights, ages, eps):
    """Two-age mixture on the outlier-robust conditional likelihood, with the
    registered grid-edge guard (no component within one native step of an edge)."""
    lo, hi = B.OUTLIER_LOGTE
    log_u = np.log(eps) - np.log(hi - lo)
    single = B.with_outliers(joint - marg, eps)
    j_single = int(np.argmax(weights @ single))
    lnl_single = float((weights @ single)[j_single])
    best = {"lnL": -np.inf}
    n = len(ages)
    for j1 in range(n):
        for j2 in range(j1 + 1, n):
            for phi in I20.PHI_GRID:
                a, b = np.log(phi), np.log(1 - phi)
                num = np.logaddexp(a + joint[:, j1], b + joint[:, j2])
                den = np.logaddexp(a + marg[:, j1], b + marg[:, j2])
                ll = np.logaddexp(np.log1p(-eps) + num - den, log_u)
                tot = float(weights @ ll)
                if tot <= best["lnL"]:
                    continue
                resp = np.exp(a + joint[:, j1] - num)
                n1 = float(weights @ resp)
                n2 = float(weights.sum() - n1)
                if min(n1, n2) < I20.MIN_COMPONENT_STARS:
                    continue
                best = {"lnL": tot, "j1": j1, "j2": j2, "age_young": float(ages[j1]),
                        "age_old": float(ages[j2]), "phi_young": float(phi), "n_young": n1, "n_old": n2}
    if not np.isfinite(best["lnL"]):
        return {"detected": False, "reason": "no admissible mixture"}
    dbic = I20.delta_bic(best["lnL"], lnl_single, weights.sum())
    edge = any(j in (0, 1, n - 2, n - 1) for j in (best["j1"], best["j2"]))
    detected = bool(dbic > 10 and best["age_old"] - best["age_young"] > 1.0 and not edge)
    return {**best, "single_best_age": float(ages[j_single]), "lnL_single": lnl_single,
            "delta_BIC": float(dbic), "component_at_grid_edge": bool(edge), "detected": detected}


def j1_synthetic(iso, rng) -> list[dict]:
    rows = []
    for fam in B.FAMILIES:
        ages = np.unique(iso[fam].age_Myr.to_numpy(float))
        parts = {a: B.cmd_particles(I20.isochrone_at(iso[fam], a), fam) for a in ages}
        k = band_coefficients(3.1)
        for true in (2.51, 3.98):
            p = B.cmd_particles(I20.isochrone_at(iso[fam], true), fam)
            inside = np.where(p["G"] <= B.UMS_EDGE)[0]
            j = rng.choice(inside, 400, p=p["w"][inside] / p["w"][inside].sum())
            av = rng.uniform(4, 8, 400)
            df = pd.DataFrame({"g": p["G"][j] + k["G"] * av + rng.normal(0, .07, 400),
                               "c": p["c"][j] + (k["BP"] - k["RP"]) * av + rng.normal(0, .04, 400),
                               "x": p["x"][j] + (k["J"] - k["Ks"]) * av + rng.normal(0, .05, 400),
                               "sig_g": .07, "sig_c": .04, "sig_x": .05})
            df.attrs.update(kG=k["G"], kc=k["BP"] - k["RP"], kx=k["J"] - k["Ks"])
            ll = np.array([B.cmd_loglike(df, parts[a], "nir", 3.1).sum() for a in ages])
            r = I20.posterior_summary(ages, ll)
            native = float(ages[np.argmin(np.abs(ages - true))])
            step = int(abs(np.argmin(np.abs(ages - r["map"])) - np.argmin(np.abs(ages - native))))
            rows.append({"family": fam, "true": native, "map": r["map"], "grid_steps_off": step,
                         "pass": step <= 1})
    return rows


def main() -> None:
    t0 = time.time()
    method = B.check_method_unchanged()
    iso = {"PARSEC": pd.read_parquet(B.frozen("isochrones_parsec")),
           "MIST": pd.read_parquet(B.frozen("isochrones_mist"))}
    anchors = pd.read_parquet(B.frozen("anchor_hrd"))
    ext = pd.read_parquet(B.frozen("extinction_v5")).drop(columns=["subgroup"], errors="ignore")
    labels = pd.read_parquet(B.frozen("subgroup_labels"))
    esphs = pd.read_csv(B.frozen("esphs"))
    members = ext.merge(labels, on="source_id", how="inner")
    members = members[members.subgroup.isin(B.SUBGROUPS)]
    phase_a = json.loads(B.frozen("issue20_outcome").read_text())

    cal_rows, rows, curves = [], [], {}
    for rv in B.R_V_BRANCHES:
        nir = B.nir_scale(anchors, ext, rv)
        esc = B.esphs_calibration(anchors, esphs, ext, rv)
        cal_rows += [{"test": "a", **nir}, {"test": "b", **esc}]
        for fam in B.FAMILIES:
            ages = np.unique(iso[fam].age_Myr.to_numpy(float))
            cmd_parts = [B.cmd_particles(I20.isochrone_at(iso[fam], a), fam) for a in ages]
            hrd_parts = [I20.hrd_particles(I20.isochrone_at(iso[fam], a), fam, B.F_BIN) for a in ages]
            spec = A3.star_inputs(A3.load_sample(), {
                "R_V": rv, "f_bin": B.F_BIN, "sigma_mode": "primary", "cal": I20.SIG_CAL_MAG,
                "dmu": 0.0, "teff_mode": "supergiant_scale"})
            for sg in B.SUBGROUPS:
                sub = members[members.subgroup.eq(sg)]
                # ---- test a and its control
                stars = B.nir_inputs(sub, rv, nir["scale"])
                wts = stars.membership_probability.to_numpy(float)
                for mode, test in (("nir", "a"), ("wp3", "a_control")):
                    ll = np.column_stack([B.cmd_loglike(stars, p, mode, rv) for p in cmd_parts])
                    s = B.subgroup_posterior(ages, ll, wts)
                    valid = bool(nir["valid"] and len(stars) >= B.MIN_STARS and not s["railed"])
                    rows.append({"test": test, "eps": np.nan, "subgroup": sg, "family": fam,
                                 "R_V": rv, "n_stars": len(stars), "n_eff": float(wts.sum()),
                                 **s, "valid": valid if test == "a" else np.nan})
                    curves[(test, np.nan, sg, fam, rv)] = (wts @ ll).tolist()
                # ---- test b
                bst = B.esphs_inputs(sub, esphs, esc)
                bw = bst.membership_probability.to_numpy(float)
                if len(bst):
                    jm = [I20.star_terms(bst.logTe.values, bst.MG0.values, bst.sig_logTe.values,
                                         bst.sig_MG0.values, p) for p in hrd_parts]
                    joint = np.column_stack([x[0] for x in jm])
                    marg = np.column_stack([x[1] for x in jm])
                for eps in EPS_VALUES:
                    if len(bst) == 0:
                        rows.append({"test": "b", "eps": eps, "subgroup": sg, "family": fam,
                                     "R_V": rv, "n_stars": 0, "valid": False})
                        continue
                    ll = B.with_outliers(joint - marg, eps)
                    s = B.subgroup_posterior(ages, ll, bw)
                    valid = bool(esc["valid"] and len(bst) >= B.MIN_STARS and not s["railed"])
                    rows.append({"test": "b", "eps": eps, "subgroup": sg, "family": fam, "R_V": rv,
                                 "n_stars": len(bst), "n_eff": float(bw.sum()), **s, "valid": valid})
                    curves[("b", eps, sg, fam, rv)] = (bw @ ll).tolist()
                # ---- test c
                cst = spec[spec.subgroup.eq(sg)]
                cw = cst.membership_probability.to_numpy(float)
                jm = [I20.star_terms(cst.logTe_used.values, cst.MG0.values, cst.sigma_logTe.values,
                                     cst.sigma_MG0.values, p) for p in hrd_parts]
                joint = np.column_stack([x[0] for x in jm])
                marg = np.column_stack([x[1] for x in jm])
                for eps in (0.0, *EPS_VALUES):
                    ll = (joint - marg) if eps == 0.0 else B.with_outliers(joint - marg, eps)
                    s = B.subgroup_posterior(ages, ll, cw)
                    valid = bool(len(cst) >= B.MIN_STARS and not s["railed"])
                    rows.append({"test": "c", "eps": eps, "subgroup": sg, "family": fam, "R_V": rv,
                                 "n_stars": len(cst), "n_eff": float(cw.sum()), **s,
                                 "valid": valid if eps > 0 else np.nan})
                    curves[("c", eps, sg, fam, rv)] = (cw @ ll).tolist()
                if sg != "CygOB2-B" and rv == 3.1:
                    mix = guarded_mixture(joint, marg, cw, ages, B.EPS_OUTLIER)
                    rows.append({"test": "c_mixture", "eps": B.EPS_OUTLIER, "subgroup": sg,
                                 "family": fam, "R_V": rv, "n_stars": len(cst),
                                 **{k: v for k, v in mix.items() if k not in ("j1", "j2")}})
                print(f"R_V {rv} {fam} {sg}: a {len(stars)}  b {len(bst)}  c {len(cst)}", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(OUT_TESTS, index=False)
    pd.DataFrame(cal_rows).to_csv(OUT_CAL, index=False)

    # ---- integrity
    j1 = j1_synthetic(iso, np.random.default_rng(20261006))
    j2 = []
    for fam in B.FAMILIES:
        for sg, key in (("CygOB2-A", "A"), ("CygOB2-C", "C")):
            r = table[table.test.eq("c") & table.eps.eq(0.0) & table.subgroup.eq(sg)
                      & table.family.eq(fam) & table.R_V.eq(3.1)].iloc[0]
            ref = phase_a["decision_cells"][fam][key]["map"]
            j2.append({"family": fam, "subgroup": sg, "c_eps0_map": float(r["map"]),
                       "phase_a_map": ref, "abs_diff": abs(float(r["map"]) - ref)})
    w.write_json(OUT_EXEC, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20b_run.py",
        "preregistration": {"path": "provenance/issue20b_prereg.json",
                            "sha256": w.sha256(B.PREREG_PATH)},
        "method_check": method,
        "inputs": B.prereg()["consumed_inputs"],
        "calibrations": cal_rows,
        "integrity": {"J1": {"pass": all(r["pass"] for r in j1), "rows": j1},
                      "J2": {"pass": max(r["abs_diff"] for r in j2) < 1e-6, "rows": j2},
                      "J3": "inputs via issue20b_common.frozen; method hash checked"},
        "loglike_curves": {"|".join(map(str, k)): v for k, v in curves.items()},
        "native_ages": {f: np.unique(iso[f].age_Myr.to_numpy(float)).tolist() for f in B.FAMILIES},
        "runtime_s": round(time.time() - t0, 1),
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (OUT_TESTS, OUT_CAL)},
    })
    print(f"J1 {'PASS' if all(r['pass'] for r in j1) else 'FAIL'}; "
          f"J2 {'PASS' if max(r['abs_diff'] for r in j2) < 1e-6 else 'FAIL'}; "
          f"{round(time.time() - t0, 1)} s")


if __name__ == "__main__":
    main()
