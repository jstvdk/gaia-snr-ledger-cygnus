#!/usr/bin/env python3
"""Issue #20 Phase A, steps A4 and A5 -- photometric sensitivity refits of C
and attribution of the 3.98 -> 2.51 Myr jump.  Sensitivity only: nothing here
is ever adopted as C's age.

A4(a) hybrid: C's upper-MS likelihood (wp4_common, unchanged) for its
      photometric-only members, plus, for its spectroscopic members, the joint
      (logTe, M_G0) density in the same window (M_G0 <= 1.5) from
      issue20_common.hrd_particles with the A3 T_eff and sigmas.  A is run as
      a control (report only).
A4(b) the WP4 upper-MS fit of C with its luminosity class I/II anchors removed.
A5    per-star Delta A_V and Delta M_G0 between the pre-repair and repair_v1
      extinction for C's anchors and 15 brightest members, and counterfactual
      C fits on repair_v1 with the anchors' (or the non-anchors') rows taken
      from the pre-repair run.

Integrity: I1 (repair_v5 C fits reproduce the stored MAPs) and I2 (pre-repair
and repair_v1 C fits reproduce theirs).

Outputs:
  tables/issue20_photometric_refits.csv
  tables/issue20_extinction_attribution.csv
  provenance/issue20_photometric_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20_photometric.py
"""
from __future__ import annotations

import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I
import wp4_common as W4
import wp5_common as w
from wp3_extinction_law import band_coefficients

OUT_REFIT = w.TABLES / "issue20_photometric_refits.csv"
OUT_ATTR = w.TABLES / "issue20_extinction_attribution.csv"
OUT_EXEC = w.PROVENANCE / "issue20_photometric_execution.json"

F_BIN = 0.4
HYBRID_EDGE = W4.UMS_FAINT_EDGE   # 1.5


def members(extinction_name: str) -> pd.DataFrame:
    ext = pd.read_parquet(I.frozen(extinction_name)).drop(
        columns=["subgroup", "subgroup_label"], errors="ignore")
    labels = pd.read_parquet(I.frozen("subgroup_labels"))
    m = ext.merge(labels[["source_id", "subgroup"]], on="source_id", how="left",
                  validate="one_to_one")
    m["subgroup"] = m["subgroup"].fillna("unassigned")
    return m


def anchors() -> pd.DataFrame:
    a = pd.read_parquet(I.frozen("anchor_hrd"))
    a = a[a.subgroup.isin(I.SUBGROUPS)].copy()
    parsed = a.spectral_type.map(I.parse_spectral_type)
    a["letter"] = [p["letter"] for p in parsed]
    a["lum_class"] = [p["lum_class"] for p in parsed]
    a["usable"] = ~a.extreme_hot & a.letter.notna()
    return a


def ums_loglike(fit_sub: pd.DataFrame, iso_family: pd.DataFrame):
    """Per-star upper-MS log-likelihood on the native grid (wp4_common)."""
    grid = np.sort(iso_family["age_Myr"].unique())
    mask = W4._window_mask_stars(fit_sub, "ums")
    stars = fit_sub[mask & np.isfinite(fit_sub["colour"]) & np.isfinite(fit_sub["MG0"])]
    out = np.empty((len(stars), len(grid)))
    for j, a in enumerate(grid):
        parts = W4.build_model_particles(iso_family[np.isclose(iso_family["age_Myr"], a)], F_BIN)
        out[:, j] = W4.star_loglike(stars, parts, "ums") if parts is not None else -50.0
    return grid, stars, out


def fit(fit_sub, iso_family):
    grid, stars, ll = ums_loglike(fit_sub, iso_family)
    total = stars["P"].values @ ll
    post = W4.posterior_from_loglike(grid, total)
    return post, int(len(stars))


def stored_map(name: str, subgroup: str, family: str, rv: float) -> float:
    post = pd.read_parquet(I.frozen(name))
    row = post[post.subgroup.eq(subgroup) & post.family.eq(family) & post.R_V.eq(rv)
               & post.f_bin.eq(F_BIN) & post.indicator.eq("ums") & post.dmu.eq(0.0)]
    return float(row.age_map.iloc[0])


def hybrid(m: pd.DataFrame, anc: pd.DataFrame, subgroup: str, family: str,
           rv: float, iso_family: pd.DataFrame) -> tuple[dict, int, int]:
    fit_all = W4.branch_photometry(m, rv)
    fit_sub = fit_all[fit_all.subgroup.eq(subgroup)]
    spec = anc[anc.subgroup.eq(subgroup) & anc.usable
               & anc[f"MG0_obs_rv{rv:.1f}"].notna()].copy()
    phot = fit_sub[~fit_sub.source_id.isin(spec.source_id)]
    grid, phot_stars, ll_phot = ums_loglike(phot, iso_family)

    ext = m.set_index("source_id")
    teff = [I.assign_teff(r.spectral_type, r.teff_spec) for _, r in spec.iterrows()]
    logte = np.array([t["logTe_used"] for t in teff])
    sig_t = np.array([t["sigma_logTe"] for t in teff])
    mg = spec[f"MG0_obs_rv{rv:.1f}"].to_numpy(float)
    sig_m = I.sigma_mg0(ext.loc[spec.source_id, "G_err"].values,
                        ext.loc[spec.source_id, f"av_err_rv{rv:.1f}"].values,
                        band_coefficients(rv)["G"])
    in_window = mg <= HYBRID_EDGE
    ll_spec = np.empty((int(in_window.sum()), len(grid)))
    for j, a in enumerate(grid):
        p = I.hrd_particles(I.isochrone_at(iso_family, a), family, F_BIN)
        joint, _ = I.star_terms(logte[in_window], mg[in_window], sig_t[in_window],
                                sig_m[in_window], p)
        ll_spec[:, j] = joint - I.window_log_norm(p, HYBRID_EDGE)
    total = (phot_stars["P"].values @ ll_phot
             + spec.membership_probability.values[in_window] @ ll_spec)
    return W4.posterior_from_loglike(grid, total), len(phot_stars), int(in_window.sum())


def main() -> None:
    t0 = time.time()
    method = I.check_method_unchanged()
    iso = {"PARSEC": pd.read_parquet(I.frozen("isochrones_parsec")),
           "MIST": pd.read_parquet(I.frozen("isochrones_mist"))}
    anc = anchors()
    m5 = members("extinction_v5")
    rows, i1 = [], []

    for family in I.FAMILIES:
        for rv in I.R_V_BRANCHES:
            fit_all = W4.branch_photometry(m5, rv)
            # I1 -- reproduction of the stored repair_v5 MAP
            post, n = fit(fit_all[fit_all.subgroup.eq("CygOB2-C")], iso[family])
            ref = stored_map("age_posteriors_v5", "CygOB2-C", family, rv)
            i1.append({"family": family, "R_V": rv, "recomputed": post["map"],
                       "stored": ref, "abs_diff": abs(post["map"] - ref)})
            rows.append({"step": "A4_reference", "subgroup": "CygOB2-C", "family": family,
                         "R_V": rv, "variant": "repair_v5 upper-MS, unchanged",
                         "n_photometric": n, "n_spectroscopic": 0,
                         **{f"age_{k}": post[k] for k in ("map", "lo68", "hi68", "lo90", "hi90")}})
            # A4(b) -- class I/II anchors of C removed
            drop = anc[anc.subgroup.eq("CygOB2-C") & anc.lum_class.isin(["I", "II"])].source_id
            sub = fit_all[fit_all.subgroup.eq("CygOB2-C") & ~fit_all.source_id.isin(drop)]
            post, n = fit(sub, iso[family])
            rows.append({"step": "A4b", "subgroup": "CygOB2-C", "family": family, "R_V": rv,
                         "variant": f"class I/II anchors removed ({len(drop)})",
                         "n_photometric": n, "n_spectroscopic": 0,
                         **{f"age_{k}": post[k] for k in ("map", "lo68", "hi68", "lo90", "hi90")}})
            # A4(a) -- hybrid, C and the A control
            for subgroup in ("CygOB2-C", "CygOB2-A"):
                post, n_phot, n_spec = hybrid(m5, anc, subgroup, family, rv, iso[family])
                rows.append({"step": "A4a" if subgroup == "CygOB2-C" else "A4a_control",
                             "subgroup": subgroup, "family": family, "R_V": rv,
                             "variant": "hybrid CMD + spectroscopic HRD",
                             "n_photometric": n_phot, "n_spectroscopic": n_spec,
                             **{f"age_{k}": post[k] for k in ("map", "lo68", "hi68", "lo90", "hi90")}})
            # reference: A's stored repair_v5 MAP (S2 threshold)
            rows.append({"step": "A_reference", "subgroup": "CygOB2-A", "family": family,
                         "R_V": rv, "variant": "stored repair_v5 upper-MS MAP",
                         "age_map": stored_map("age_posteriors_v5", "CygOB2-A", family, rv)})
            print(f"A4 {family} R_V={rv} done", flush=True)

    # ---- A5
    m_pre = members("extinction_pre_repair")
    m_v1 = members("extinction_v1")
    c_anchor_ids = set(anc[anc.subgroup.eq("CygOB2-C")].source_id)
    c5 = m5[m5.subgroup.eq("CygOB2-C")]
    brightest = set(c5.nsmallest(15, "G0_abs_rv3.1").source_id)
    ids = sorted(c_anchor_ids | brightest)
    cols = ["source_id", "av_rv3.1", "G0_abs_rv3.1", "BPRP0_rv3.1", "av_method"]
    a = (m_pre[cols].set_index("source_id").add_suffix("_pre")
         .join(m_v1[cols].set_index("source_id").add_suffix("_v1"))
         .join(m5[cols].set_index("source_id").add_suffix("_v5")).loc[ids])
    a["dAV_v1_minus_pre"] = a["av_rv3.1_v1"] - a["av_rv3.1_pre"]
    a["dMG0_v1_minus_pre"] = a["G0_abs_rv3.1_v1"] - a["G0_abs_rv3.1_pre"]
    a["is_anchor"] = a.index.isin(c_anchor_ids)
    a["in_15_brightest_v5"] = a.index.isin(brightest)
    a = a.join(anc.set_index("source_id")[["spectral_type", "lum_class"]], how="left")
    a.reset_index().sort_values("G0_abs_rv3.1_v5").to_csv(OUT_ATTR, index=False)

    i2 = []
    rv = 3.1
    for family in I.FAMILIES:
        for label, frame, ref_name in (("pre_repair", m_pre, "age_posteriors_pre_repair"),
                                       ("repair_v1", m_v1, "age_posteriors_v1")):
            fit_all = W4.branch_photometry(frame, rv)
            post, n = fit(fit_all[fit_all.subgroup.eq("CygOB2-C")], iso[family])
            ref = stored_map(ref_name, "CygOB2-C", family, rv)
            i2.append({"family": family, "version": label, "recomputed": post["map"],
                       "stored": ref, "abs_diff": abs(post["map"] - ref)})
            rows.append({"step": "A5", "subgroup": "CygOB2-C", "family": family, "R_V": rv,
                         "variant": f"{label} extinction, unchanged", "n_photometric": n,
                         **{f"age_{k}": post[k] for k in ("map", "lo68", "hi68", "lo90", "hi90")}})
        for label, swap in (("repair_v1, anchors' rows from pre-repair", True),
                            ("repair_v1, non-anchors' rows from pre-repair", False)):
            pre = m_pre.set_index("source_id")
            hyb = m_v1.set_index("source_id").copy()
            chosen = [s for s in hyb.index if (s in c_anchor_ids) == swap]
            shared = [c for c in hyb.columns if c in pre.columns and c != "subgroup"]
            hyb.loc[chosen, shared] = pre.loc[chosen, shared].values
            fit_all = W4.branch_photometry(hyb.reset_index(), rv)
            post, n = fit(fit_all[fit_all.subgroup.eq("CygOB2-C")], iso[family])
            rows.append({"step": "A5", "subgroup": "CygOB2-C", "family": family, "R_V": rv,
                         "variant": label, "n_photometric": n,
                         **{f"age_{k}": post[k] for k in ("map", "lo68", "hi68", "lo90", "hi90")}})
        print(f"A5 {family} done", flush=True)

    pd.DataFrame(rows).to_csv(OUT_REFIT, index=False)
    i1_pass = bool(max(r["abs_diff"] for r in i1) <= 0.01)
    i2_pass = bool(max(r["abs_diff"] for r in i2) <= 0.01)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_photometric.py",
        "preregistration": "provenance/issue20_prereg.json",
        "method_check": method,
        "inputs": {n: I.prereg()["consumed_inputs"][n] for n in
                   ("extinction_v5", "extinction_pre_repair", "extinction_v1",
                    "subgroup_labels", "anchor_hrd", "age_posteriors_v5",
                    "age_posteriors_pre_repair", "age_posteriors_v1",
                    "isochrones_parsec", "isochrones_mist")},
        "f_bin": F_BIN,
        "integrity": {"I1": {"pass": i1_pass, "rows": i1},
                      "I2": {"pass": i2_pass, "rows": i2}},
        "runtime_s": round(time.time() - t0, 1),
        "outputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in (OUT_REFIT, OUT_ATTR)},
    }
    w.write_json(OUT_EXEC, record)
    print(f"I1 {'PASS' if i1_pass else 'FAIL'}; I2 {'PASS' if i2_pass else 'FAIL'}; "
          f"{record['runtime_s']} s")


if __name__ == "__main__":
    main()
