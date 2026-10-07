#!/usr/bin/env python3
"""repair_v9 Stage 1 -- score GA1-GA3 and W1-W4, apply the adoption / fallback
rule, and write the repair_v9 WP4 age table (provenance/wp4v9_age_prereg.json).

Outputs:
  data/processed/wp4_age_posteriors_repair_v9.parquet
  provenance/wp4v9_age_outcome.json

Run:
  PYTHONPATH=scripts python3 scripts/wp4v9_score.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

import issue20_common as I20
import wp4v9_common as V
import wp5_common as w
from wp6_mass_extension_decision import turnoff_mass

OUT = w.PROC / "wp4_age_posteriors_repair_v9.parquet"
OUT_JSON = w.PROVENANCE / "wp4v9_age_outcome.json"
DEC_RV, DEC_FB = V.DECISION["R_V"], V.DECISION["f_bin"]
SCHEMA = ["subgroup", "family", "R_V", "f_bin", "indicator", "dmu", "n_stars", "measurable",
          "grid_railed", "exclusion_reason", "age_map", "age_lo68", "age_hi68", "age_lo90",
          "age_hi90", "age_mean"]


def fine_posterior(ages, total):
    la = np.log10(np.asarray(ages, float))
    laf = np.linspace(la.min(), la.max(), 4001)
    f = PchipInterpolator(la, np.asarray(total) - np.max(total))(laf)
    p = np.exp(f)
    p /= np.trapezoid(p, laf)
    return laf, p


def ga3_fraction(fam, ages, total, m_lo):
    laf, p = fine_posterior(ages, total)
    bad = np.array([turnoff_mass(fam, 10 ** x) < m_lo for x in laf])
    return float(np.trapezoid(p * bad, laf))


def m_lo_by_subgroup(fam, rv):
    anchors = pd.read_parquet(V.frozen("anchor_hrd_v8"))
    iso = pd.read_parquet(V.frozen(f"isochrones_{fam.lower()}"))
    iso = iso.dropna(subset=["logTe", "G0", "Mini"])
    out = {}
    for sg in V.SUBGROUPS:
        a = anchors[anchors.subgroup.eq(sg) & ~anchors.extreme_hot
                    & anchors[f"MG0_obs_rv{rv:.1f}"].notna()]
        rows = []
        for _, r in a.iterrows():
            p = I20.parse_spectral_type(r.spectral_type)
            if p["letter"] not in ("O", "B") or p["composite"] or p["lum_class"] in ("I", "II"):
                continue
            m = V.minimum_initial_mass(r.logTe_spec, r[f"MG0_obs_rv{rv:.1f}"], iso)
            rows.append((m, r.spectral_type, int(r.source_id)))
        rows = [x for x in rows if np.isfinite(x[0])]
        best = max(rows) if rows else (np.nan, None, None)
        out[sg] = {"m_lo": best[0], "star": best[1], "source_id": best[2], "n_candidates": len(rows)}
    return out


def summary(ages, total, n):
    return V.posterior(np.asarray(ages), np.asarray(total), n)


def main() -> None:
    V.check_method_unchanged()
    fx = json.loads((w.PROVENANCE / "wp4v9_fit_execution.json").read_text())
    curves, nat = fx["loglike_curves"], {f: np.array(v) for f, v in fx["native_ages"].items()}
    fits = pd.read_csv(w.TABLES / "wp4v9_real_fits.csv")
    faint = pd.read_csv(w.TABLES / "wp4v9_faint_end_study.csv")
    inj = pd.read_csv(w.TABLES / "wp4v9_injection_validation.csv")
    spec_tests = pd.read_csv(V.frozen("spectroscopic_tests"))
    spec_curves = json.loads(V.frozen("spectroscopic_curves").read_text())["loglike_curves"]

    def curve(fit, window, model, sg, fam, rv, fb, dmu=0.0):
        return np.array(curves[f"{fit}|{window}|{model}|{sg}|{fam}|{rv}|{fb}|{dmu}"])

    def row(fit, window, model, sg, fam, rv, fb, dmu=0.0):
        r = fits[fits.fit.eq(fit) & fits.window.eq(window) & fits.error_model.eq(model)
                 & fits.subgroup.eq(sg) & fits.family.eq(fam) & np.isclose(fits.R_V, rv)
                 & np.isclose(fits.f_bin, fb) & np.isclose(fits.dmu, dmu)]
        return r.iloc[0]

    # ---- GA1
    ga1 = {}
    for (fam, sg), blk in inj[inj.window.eq("bright") & inj.model.eq("M1")].groupby(["family", "subgroup"]):
        bias = blk.groupby("true_age").bias_median.mean()
        cov = float(blk.covered68.mean())
        ga1[f"{sg}|{fam}"] = {"bias_by_age": {f"{k:.2f}": float(v) for k, v in bias.items()},
                              "coverage68": cov, "n": int(len(blk)),
                              "pass": bool((bias.abs() < 0.2).all() and 0.55 <= cov <= 0.80)}
    old_bias = {}
    for (fam, window), blk in inj[inj.model.eq("old")].groupby(["family", "window"]):
        old_bias[f"{fam}|{window}"] = float(blk.groupby(["subgroup", "true_age"]).bias_median.mean().mean())
    m1_other = {f"{fam}|{window}": float(b.groupby(["subgroup", "true_age"]).bias_median.mean().mean())
                for (fam, window), b in inj[inj.model.eq("M1")].groupby(["family", "window"])}

    # ---- GA2, GA3, adoption per subgroup x family (decision cell), qualifiers on the rest
    m_lo = {(fam, rv): m_lo_by_subgroup(fam, rv) for fam in V.FAMILIES for rv in V.R_V_BRANCHES}

    def spec_row(sg, fam, rv):
        s = spec_tests[spec_tests.test.eq("c") & np.isclose(spec_tests.eps, 0.05) & spec_tests.subgroup.eq(sg)
                       & spec_tests.family.eq(fam) & np.isclose(spec_tests.R_V, rv)]
        return s.iloc[0]

    def checks(sg, fam, rv, fb):
        r = row("M1", "bright", "M1", sg, fam, rv, fb)
        out = {"M1": {k: float(r[k]) for k in ("age_map", "age_median", "age_lo68", "age_hi68", "age_lo95", "age_hi95")},
               "measurable": bool(r.measurable)}
        if sg != "CygOB2-B":
            s = spec_row(sg, fam, rv)
            out["spec68"] = [float(s.lo68), float(s.hi68)]
            out["GA2"] = bool(not (r.age_hi68 < s.lo68 or s.hi68 < r.age_lo68))
        bound = m_lo[(fam, rv)][sg]["m_lo"]
        out["GA3_fraction"] = ga3_fraction(fam, nat[fam], curve("M1", "bright", "M1", sg, fam, rv, fb), bound)
        out["GA3"] = bool(out["GA3_fraction"] < 0.05)
        return out

    decisions, qual, rows_out = {}, {}, []
    for fam in V.FAMILIES:
        for sg in V.SUBGROUPS:
            d = checks(sg, fam, DEC_RV, DEC_FB)
            ga1_ok = ga1[f"{sg}|{fam}"]["pass"]
            if not ga1_ok or not d["GA3"] or not d["measurable"]:
                mode = "not_adopted"
            elif sg == "CygOB2-B" or d["GA2"]:
                mode = "M1"
            else:
                mode = "joint"
            d.update(GA1=ga1_ok, mode=mode)
            if mode == "joint":
                jc = (curve("M1_nonanchor", "bright", "M1", sg, fam, DEC_RV, DEC_FB)
                      + np.array(spec_curves[f"c|0.05|{sg}|{fam}|{DEC_RV}"]))
                d["joint"] = {k: float(v) for k, v in summary(nat[fam], jc, 0).items()
                              if k.startswith("age_")}
                d["joint_bimodal"] = V.bimodal(nat[fam], jc)
                d["joint_GA3_fraction"] = ga3_fraction(fam, nat[fam], jc, m_lo[(fam, DEC_RV)][sg]["m_lo"])
                if d["joint_bimodal"] or d["joint_GA3_fraction"] >= 0.05:
                    d["mode"] = mode = "branches"
            decisions[f"{sg}|{fam}"] = d
            # qualifiers
            other = {}
            for rv in V.R_V_BRANCHES:
                for fb in V.F_BINS:
                    if rv == DEC_RV and fb == DEC_FB:
                        continue
                    c = checks(sg, fam, rv, fb)
                    other[f"{rv}|{fb}"] = {k: c.get(k) for k in ("GA2", "GA3")}
            same = all(v.get("GA2") == d.get("GA2") and v["GA3"] == d["GA3"] for v in other.values())
            qual[f"{sg}|{fam}"] = {"robust": same, "cells": other}

            # ---- rows of the repair_v9 table
            for rv in V.R_V_BRANCHES:
                for fb in V.F_BINS:
                    dmus = (0.0, 0.06, -0.06) if (rv == DEC_RV and fb == DEC_FB) else (0.0,)
                    for dmu in dmus:
                        base = row("M1", "bright", "M1", sg, fam, rv, fb, dmu)
                        if mode == "not_adopted" and sg != "CygOB2-B":
                            # prereg 'if_GA1_or_GA3_fails': repair_v9 carries the test-c spectroscopic
                            # posterior (eps 0.05, same family and R_V) as the only age; test c has no
                            # distance or f_bin variants, so every dmu / f_bin key gets the same curve
                            sr = spec_row(sg, fam, rv)
                            s = summary(nat[fam], np.array(spec_curves[f"c|0.05|{sg}|{fam}|{rv}"]),
                                        int(sr.n_stars))
                            rows_out.append({"subgroup": sg, "family": fam, "R_V": rv, "f_bin": fb,
                                             "indicator": "ums", "dmu": dmu,
                                             **{k: s[k] for k in SCHEMA if k in s},
                                             "n_stars": int(sr.n_stars),
                                             "window": "spectroscopic_hrd", "error_model": "test_c",
                                             "fallback": "GA1_or_GA3_failed:spectroscopic"
                                                         + ("" if dmu == 0.0 else ":no_dmu_refit"),
                                             "adopted": True, "branch": "", "photometry_only": False})
                            rows_out.append({"subgroup": sg, "family": fam, "R_V": rv, "f_bin": fb,
                                             "indicator": "ums_not_adopted", "dmu": dmu,
                                             **{k: base[k] for k in SCHEMA if k in base.index
                                                and k not in ("subgroup", "family", "R_V", "f_bin",
                                                              "indicator", "dmu")},
                                             "window": "bright", "error_model": "M1_marginalised",
                                             "fallback": "GA1_or_GA3_failed", "adopted": False,
                                             "branch": "", "photometry_only": False})
                            continue
                        if mode in ("M1", "branches", "not_adopted"):
                            # not_adopted here is B only: the M1 posterior, flagged
                            src = {"window": "bright", "error_model": "M1_marginalised",
                                   "fallback": {"M1": "", "branches": "branches",
                                                "not_adopted": "failed_GA1_or_GA3"}[mode]}
                            post = {k: base[k] for k in SCHEMA if k in base.index}
                        else:
                            jc = (curve("M1_nonanchor", "bright", "M1", sg, fam, rv, fb, dmu)
                                  + np.array(spec_curves[f"c|0.05|{sg}|{fam}|{rv}"]))
                            s = summary(nat[fam], jc, int(base.n_stars))
                            post = {**{k: s[k] for k in SCHEMA if k in s}, "n_stars": int(base.n_stars)}
                            src = {"window": "bright", "error_model": "M1_x_spectroscopic_joint",
                                   "fallback": "joint"}
                        rows_out.append({**{k: post[k] for k in post if k in SCHEMA},
                                         "subgroup": sg, "family": fam, "R_V": rv, "f_bin": fb,
                                         "indicator": "ums", "dmu": dmu,
                                         **src, "adopted": True,
                                         "branch": "photometric" if mode == "branches" else "",
                                         "photometry_only": sg == "CygOB2-B"})
                        if mode == "branches" and dmu == 0.0:
                            sc = np.array(spec_curves[f"c|0.05|{sg}|{fam}|{rv}"])
                            s = summary(nat[fam], sc, int(spec_row(sg, fam, rv).n_stars))
                            rows_out.append({"subgroup": sg, "family": fam, "R_V": rv, "f_bin": fb,
                                             "indicator": "ums_branch_spectroscopic", "dmu": dmu,
                                             **{k: s[k] for k in SCHEMA if k in s},
                                             "window": "spectroscopic_hrd", "error_model": "test_c",
                                             "fallback": "branches", "adopted": True,
                                             "branch": "spectroscopic", "photometry_only": False})
    # every other fit as sensitivity rows
    for _, r in fits.iterrows():
        if r.fit == "M1" and r.window == "bright":
            continue
        rows_out.append({**{k: r[k] for k in SCHEMA if k in r.index and k != "indicator"},
                         "indicator": "ums_sensitivity", "window": r.window,
                         "error_model": f"{r.fit}:{r.error_model}", "fallback": "", "adopted": False,
                         "branch": "", "photometry_only": r.subgroup == "CygOB2-B"})
    table = pd.DataFrame(rows_out)
    key = ["subgroup", "family", "R_V", "f_bin", "dmu"]
    dup = table[table.indicator.eq("ums")].duplicated(key).any()
    if dup:
        raise RuntimeError("more than one adopted 'ums' row per key")
    table.to_parquet(OUT, index=False)

    # ---- predictions
    def adopted_map(sg, fam):
        # the age repair_v9 carries at the decision cell (the 'ums' row)
        r = table[table.indicator.eq("ums") & table.subgroup.eq(sg) & table.family.eq(fam)
                  & np.isclose(table.R_V, DEC_RV) & np.isclose(table.f_bin, DEC_FB) & np.isclose(table.dmu, 0.0)]
        return float(r.age_map.iloc[0])
    w1 = all(adopted_map("CygOB2-C", f) >= 3.2 for f in V.FAMILIES)
    w2 = all(abs(adopted_map("CygOB2-A", f) - adopted_map("CygOB2-C", f)) < 1.0 for f in V.FAMILIES)
    w3 = all(old_bias[f"{f}|faint"] < -0.2 for f in V.FAMILIES)
    bright_dec = {(sg, f): row("M1", "bright", "M1", sg, f, DEC_RV, DEC_FB) for sg in V.SUBGROUPS for f in V.FAMILIES}

    def overlaps(r, b):
        return not (r.age_hi68 < b.age_lo68 or b.age_hi68 < r.age_lo68)
    base_faint = faint[faint.fit.eq("M4_base")]
    moot = all(overlaps(r, bright_dec[(r.subgroup, r.family)]) for _, r in base_faint.iterrows() if r.measurable)
    causes = {}
    for cause, blk in faint[~faint.fit.eq("M4_base")].groupby("fit"):
        ok = [overlaps(r, bright_dec[(r.subgroup, r.family)]) for _, r in blk.iterrows() if r.measurable]
        causes[cause] = {"reconciles": bool(ok and all(ok)), "n_cells": len(ok),
                         "cells_overlapping": int(sum(ok))}
    w4 = "moot" if moot else ("PASS" if not any(c["reconciles"] for c in causes.values()) else "FAIL")

    widened = {}
    for sg in V.SUBGROUPS:
        for fam in V.FAMILIES:
            for rv in V.R_V_BRANCHES:
                p = row("M1", "bright", "M1", sg, fam, rv, 0.4)
                q = row("M1_widened_x2", "bright", "M1", sg, fam, rv, 0.4)
                half = 0.5 * (p.age_hi68 - p.age_lo68)
                widened[f"{sg}|{fam}|{rv}"] = {
                    "median_shift": float(q.age_median - p.age_median),
                    "width68_primary": float(p.age_hi68 - p.age_lo68),
                    "width68_widened": float(q.age_hi68 - q.age_lo68),
                    "double_counting_sensitive": bool(abs(q.age_median - p.age_median) > half)}

    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp4v9_score.py",
        "preregistration": {"path": "provenance/wp4v9_age_prereg.json", "sha256": w.sha256(V.PREREG_PATH)},
        "GA1": ga1, "old_model_injection_bias": old_bias, "M1_injection_bias_other_windows": m1_other,
        "m_lo": {f"{k[0]}|{k[1]}": v for k, v in m_lo.items()},
        "decisions": decisions, "qualifiers": qual,
        "adopted_age_map_decision_cell": {k: adopted_map(*k.split("|")) for k in decisions},
        "predictions": {"W1": w1, "W2": w2, "W3": w3, "W4": w4, "W4_causes": causes},
        "double_counting_sensitivity": widened,
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    print("GA1", {k: v["pass"] for k, v in ga1.items()})
    for k, d in decisions.items():
        print(k, d["mode"], "M1", round(d["M1"]["age_map"], 2), [round(d["M1"]["age_lo68"], 2), round(d["M1"]["age_hi68"], 2)],
              "GA2", d.get("GA2"), "spec", d.get("spec68"), "GA3", round(d["GA3_fraction"], 3),
              "joint", {k2: round(v, 2) for k2, v in d.get("joint", {}).items() if k2 in ("age_map", "age_lo68", "age_hi68")})
    print("W1", w1, "W2", w2, "W3", w3, "W4", w4)


if __name__ == "__main__":
    main()
