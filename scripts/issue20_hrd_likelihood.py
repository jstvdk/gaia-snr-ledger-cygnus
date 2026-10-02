#!/usr/bin/env python3
"""Issue #20 Phase A, step A3 -- spectroscopic-HRD age likelihood per subgroup.

Implements provenance/issue20_prereg.json "a3_method" with the method module
scripts/issue20_common.py (hash-checked).  For every run in the registered
settings grid it evaluates, per star and per native isochrone age, the joint
and marginal particle densities, and from them

  * the subgroup posteriors (conditional = primary, windowed joint = secondary),
  * the two-age mixture for A and C (H3, and A as the specificity control),
  * per-star ages on the decision cells.

Integrity checks I3 (the scope diagnostic's chi^2 scan, recomputed from the
versioned anchor table) and I5 (particle IMF bookkeeping) are run here.

Outputs:
  tables/issue20_hrd_posteriors.csv
  tables/issue20_hrd_mixture.csv
  tables/issue20_hrd_per_star.csv
  provenance/issue20_hrd_likelihood_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20_hrd_likelihood.py
"""
from __future__ import annotations

import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I
import wp4_anchors_hrd as H
import wp5_common as w
from wp3_extinction_law import band_coefficients

OUT_POST = w.TABLES / "issue20_hrd_posteriors.csv"
OUT_MIX = w.TABLES / "issue20_hrd_mixture.csv"
OUT_STAR = w.TABLES / "issue20_hrd_per_star.csv"
OUT_EXEC = w.PROVENANCE / "issue20_hrd_likelihood_execution.json"

JOINT_EDGE = 0.0
DIAG_GRID = {
    "PARSEC": (2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01),
    "MIST": (2.00, 2.25, 2.52, 2.83, 3.18, 3.57, 4.01, 4.50, 5.04),
}


def base_settings() -> dict:
    return {"f_bin": 0.4, "sigma_mode": "primary", "cal": I.SIG_CAL_MAG,
            "dmu": 0.0, "teff_mode": "supergiant_scale"}


def run_list() -> list[dict]:
    runs = []
    for family in I.FAMILIES:
        for rv in I.R_V_BRANCHES:
            for f_bin in (0.3, 0.4, 0.5):
                tag = "decision" if (rv == 3.1 and f_bin == 0.4) else "robustness"
                runs.append({**base_settings(), "family": family, "R_V": rv,
                             "f_bin": f_bin, "role": tag})
        sens = [
            {"sigma_mode": "all_0.05"},
            {"cal": I.SIG_CAL_MAG_WIDE},
            {"dmu": +0.060},
            {"dmu": -0.060},
            {"f_bin": 0.7},
            {"teff_mode": "table"},
        ]
        for change in sens:
            runs.append({**base_settings(), "family": family, "R_V": 3.1,
                         **change, "role": "sensitivity"})
    for k, r in enumerate(runs):
        r["run_id"] = k
    return runs


def load_sample() -> pd.DataFrame:
    anchors = pd.read_parquet(I.frozen("anchor_hrd"))
    ext = pd.read_parquet(I.frozen("extinction_v5"))[
        ["source_id", "G_err", "av_err_rv3.0", "av_err_rv3.1", "av_err_rv3.5"]]
    lit = pd.read_parquet(I.frozen("spectroscopic_anchors"))
    lit["source_id"] = pd.to_numeric(lit["source_id"], errors="coerce")
    lit = lit.dropna(subset=["source_id"]).astype({"source_id": "int64"})
    s = anchors[anchors.subgroup.isin(I.SUBGROUPS)].merge(ext, on="source_id", how="left")
    s = s.merge(lit[["source_id", "object_name"]], on="source_id", how="left")
    parsed = s.spectral_type.map(I.parse_spectral_type)
    s["letter"] = [p["letter"] for p in parsed]
    s["lum_class"] = [p["lum_class"] for p in parsed]
    s["excluded"] = np.where(s.extreme_hot, "extreme_hot",
                             np.where(s.letter.isna(), "not_OB_primary", ""))
    return s.reset_index(drop=True)


def star_inputs(sample: pd.DataFrame, run: dict) -> pd.DataFrame:
    rv = run["R_V"]
    keep = sample[(sample.excluded == "") & sample[f"MG0_obs_rv{rv:.1f}"].notna()].copy()
    rows = []
    for _, r in keep.iterrows():
        a = I.assign_teff(r.spectral_type, r.teff_spec, run["sigma_mode"])
        if run["teff_mode"] == "table":
            a["teff_used"], a["logTe_used"] = float(r.teff_spec), float(np.log10(r.teff_spec))
            a["teff_source"] = "table_all"
        rows.append(a)
    t = pd.DataFrame(rows, index=keep.index)
    keep["logTe_used"] = t["logTe_used"]
    keep["teff_used"] = t["teff_used"]
    keep["teff_source"] = t["teff_source"]
    keep["sigma_logTe"] = t["sigma_logTe"]
    keep["MG0"] = keep[f"MG0_obs_rv{rv:.1f}"] - run["dmu"]
    keep["sigma_MG0"] = I.sigma_mg0(keep["G_err"], keep[f"av_err_rv{rv:.1f}"],
                                    band_coefficients(rv)["G"], cal=run["cal"])
    return keep


def i3_diagnostic_reproduction(sample: pd.DataFrame, iso: dict) -> dict:
    diag = pd.read_csv(I.frozen("diag_hrd_age_scan"))
    hrd = sample[~sample.extreme_hot & sample["MG0_obs_rv3.1"].notna()]
    worst = 0.0
    rows = []
    for subgroup in ("CygOB2-A", "CygOB2-C"):
        stars = hrd[hrd.subgroup.eq(subgroup)]
        for family, ages in DIAG_GRID.items():
            table = iso[family]
            native = np.unique(table.age_Myr.values)
            for age in ages:
                near = float(native[np.argmin(np.abs(native - age))])
                pts = H.iso_hrd_points(table[np.isclose(table.age_Myr, near)], family)
                chi = np.array([H.match_hrd(r.logTe_spec, r["MG0_obs_rv3.1"], pts)[2]
                                for _, r in stars.iterrows()])
                ours = round(float(np.sum(chi ** 2)), 1)
                ref = diag[diag.subgroup.eq(subgroup) & diag.family.eq(family)
                           & np.isclose(diag.age_Myr, round(near, 2))]
                ref_val = float(ref.chi2_all.iloc[0])
                worst = max(worst, abs(ours - ref_val))
                rows.append({"subgroup": subgroup, "family": family,
                             "age": round(near, 2), "recomputed": ours,
                             "diagnostic": ref_val, "n": int(len(stars)),
                             "n_diag": int(ref.n.iloc[0])})
    return {"pass": bool(worst <= 0.1), "worst_abs_diff": worst, "rows": rows}


def main() -> None:
    t0 = time.time()
    method = I.check_method_unchanged()
    iso = {"PARSEC": pd.read_parquet(I.frozen("isochrones_parsec")),
           "MIST": pd.read_parquet(I.frozen("isochrones_mist"))}
    ages = {f: np.unique(iso[f].age_Myr.to_numpy(float)) for f in I.FAMILIES}
    sample = load_sample()

    # ---- particle cache + I5
    particles, i5 = {}, []
    for family in I.FAMILIES:
        for f_bin in (0.3, 0.4, 0.5, 0.7):
            for age in ages[family]:
                p = I.hrd_particles(I.isochrone_at(iso[family], age), family, f_bin)
                particles[(family, f_bin, float(age))] = p
                if f_bin == 0.4:
                    i5.append({"family": family, "age": float(age),
                               "imf_check": p["imf_check"],
                               "dropped_gap_fraction": p["dropped_gap_fraction"],
                               "alive_fraction": p["alive_fraction"],
                               "n_particles": int(len(p["w"]))})
    i5_pass = bool(max(r["imf_check"] for r in i5) < 1e-9)

    post_rows, mix_rows, star_rows = [], [], []
    for run in run_list():
        family = run["family"]
        stars = star_inputs(sample, run)
        grid = ages[family]
        joint = np.empty((len(stars), len(grid)))
        marg = np.empty_like(joint)
        lnorm = np.empty(len(grid))
        for j, age in enumerate(grid):
            p = particles[(family, run["f_bin"], float(age))]
            joint[:, j], marg[:, j] = I.star_terms(
                stars.logTe_used.values, stars.MG0.values,
                stars.sigma_logTe.values, stars.sigma_MG0.values, p)
            lnorm[j] = I.window_log_norm(p, JOINT_EDGE)
        stars = stars.reset_index(drop=True)
        keys = {k: run[k] for k in ("run_id", "role", "family", "R_V", "f_bin",
                                    "sigma_mode", "cal", "dmu", "teff_mode")}
        for statistic in ("conditional", "joint"):
            ll = joint - marg if statistic == "conditional" else joint - lnorm[None, :]
            for subgroup in I.SUBGROUPS:
                sel = stars.subgroup.eq(subgroup).values
                if sel.sum() == 0:
                    continue
                weights = stars.membership_probability.values[sel]
                total = weights @ ll[sel]
                summary = I.posterior_summary(grid, total)
                post_rows.append({**keys, "statistic": statistic, "subgroup": subgroup,
                                  "n_stars": int(sel.sum()),
                                  "n_eff": float(weights.sum()), **summary})
                if subgroup == "CygOB2-B":
                    continue
                j_single, lnl_single = I.single_best(ll[sel], weights)
                mix = I.mixture_best(joint[sel], marg[sel], weights, grid,
                                     statistic=statistic, log_norm=lnorm)
                dbic = (I.delta_bic(mix["lnL"], lnl_single, weights.sum())
                        if np.isfinite(mix["lnL"]) else np.nan)
                h3 = bool(np.isfinite(dbic) and dbic > 10.0
                          and mix["age_old"] - mix["age_young"] > 1.0)
                mix_rows.append({**keys, "statistic": statistic, "subgroup": subgroup,
                                 "n_eff": float(weights.sum()),
                                 "single_best_age": float(grid[j_single]),
                                 "lnL_single": lnl_single,
                                 **{k: mix.get(k, np.nan) for k in
                                    ("age_young", "age_old", "phi_young",
                                     "n_young", "n_old")},
                                 "lnL_mix": mix["lnL"], "delta_BIC": dbic,
                                 "h3_criteria_met": h3})
            if run["role"] == "decision" and statistic == "conditional":
                for i, r in stars.iterrows():
                    s = I.posterior_summary(grid, ll[i])
                    spread = float(ll[i].max() - ll[i].min())
                    star_rows.append({
                        "family": family, "R_V": run["R_V"], "source_id": int(r.source_id),
                        "object_name": r.object_name, "subgroup": r.subgroup,
                        "membership_probability": r.membership_probability,
                        "spectral_type": r.spectral_type, "lum_class": r.lum_class,
                        "teff_table": r.teff_spec, "teff_used": r.teff_used,
                        "teff_source": r.teff_source, "sigma_logTe": r.sigma_logTe,
                        "MG0": r.MG0, "sigma_MG0": r.sigma_MG0,
                        "age_map": s["map"], "age_lo68": s["lo68"], "age_hi68": s["hi68"],
                        "lnL_spread": spread, "railed": s["railed"],
                        "defined": bool(spread >= 2.0 and not s["railed"]),
                        **{f"cond_t{a:.2f}": float(ll[i, j]) for j, a in enumerate(grid)}})
        print(f"run {run['run_id']:2d} {run['role']:11s} {family:6s} R_V={run['R_V']} "
              f"f_bin={run['f_bin']} {run['sigma_mode']} cal={run['cal']} "
              f"dmu={run['dmu']:+.3f} {run['teff_mode']}", flush=True)

    i3 = i3_diagnostic_reproduction(sample, iso)
    pd.DataFrame(post_rows).to_csv(OUT_POST, index=False)
    pd.DataFrame(mix_rows).to_csv(OUT_MIX, index=False)
    pd.DataFrame(star_rows).to_csv(OUT_STAR, index=False)

    excluded = sample[sample.excluded != ""][["source_id", "subgroup", "spectral_type", "excluded"]]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_hrd_likelihood.py",
        "preregistration": "provenance/issue20_prereg.json",
        "method_check": method,
        "inputs": {n: I.prereg()["consumed_inputs"][n] for n in
                   ("anchor_hrd", "extinction_v5", "spectroscopic_anchors",
                    "isochrones_parsec", "isochrones_mist", "diag_hrd_age_scan")},
        "n_runs": len(run_list()),
        "sample_counts": sample[sample.excluded == ""].subgroup.value_counts().to_dict(),
        "excluded": excluded.to_dict(orient="records"),
        "integrity": {
            "I3": i3,
            "I5": {"pass": i5_pass, "per_age_fbin0.4": i5},
            "I6": "every input opened through issue20_common.frozen (hash-verified)",
        },
        "runtime_s": round(time.time() - t0, 1),
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p)
                    for p in (OUT_POST, OUT_MIX, OUT_STAR)},
    }
    w.write_json(OUT_EXEC, record)
    print(f"I3 {'PASS' if i3['pass'] else 'FAIL'} (worst |diff| {i3['worst_abs_diff']:.2f}); "
          f"I5 {'PASS' if i5_pass else 'FAIL'}; {record['runtime_s']} s")


if __name__ == "__main__":
    main()
