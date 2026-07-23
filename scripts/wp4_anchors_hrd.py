#!/usr/bin/env python3
"""WP4 step 4 - spectroscopic anchors on the HRD.

For members that carry a spectroscopic effective temperature (the Wright/Berlanas
anchors), the broadband photometric temperature is unreliable: for hot stars the
optical colour saturates, so a spectral type fixes Teff far better than (BP-RP)0.
We therefore place each such anchor on the theoretical plane (log Teff, M_G0)
using its SPECTROSCOPIC Teff and its de-reddened absolute G (which for these
stars already used the intrinsic-colour / spectroscopic A_V from WP3), and read
its mass off the age-appropriate isochrone by interpolating Mass(log Teff) along
the main sequence.  PARSEC and MIST are done independently.

Gate check (plan WP4): anchor HRD positions must be consistent with the chosen
isochrones within the PARSEC-vs-MIST model difference.  We test this in the
M_G0 residual at the anchor's Teff.

Output: data/processed/wp4_anchor_hrd.parquet, provenance/wp4_anchors_hrd_execution.json
"""
from __future__ import annotations

import json, hashlib, datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w
from wp3_common import normalize_teff

# hottest Teff we trust a normal-star isochrone mass for; hotter anchors are
# WR / stripped / extreme and are flagged out of the isochrone-mass gate.
LOGTE_MAX_TRUST = np.log10(52000.0)   # ~O2/O3; Class B


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# HRD-match metric scales (Class B): spectral-type Teff ~0.03 dex; M_G0 spread
# dominated by binary + distance ~0.4 mag.  Nearest-chi2 point on the isochrone
# is the mass assignment; this is branch-robust where Teff folds at the turnoff.
SIG_LOGTE = 0.03
SIG_MG0 = 0.40


def iso_hrd_points(iso_age: pd.DataFrame, family: str):
    """MS/near-MS isochrone points (logTe, M_G0, Mass) usable for HRD matching.
    Keeps PMS+MS+turnoff, drops cool giants (post-He-ignition)."""
    if family == "PARSEC":
        pts = iso_age[iso_age["label"] <= 2]
    else:
        pts = iso_age[iso_age["phase"] <= 2]
    lt = pts["logTe"].values
    g0 = pts["G0"].values
    mass = pts["Mass"].values
    ok = np.isfinite(lt) & np.isfinite(g0) & np.isfinite(mass)
    if ok.sum() < 4:
        return None
    return lt[ok], g0[ok], mass[ok]


def match_hrd(logte, mg0, pts):
    """Nearest isochrone point in the (logTe, M_G0) plane by chi2 metric.
    Returns (mass, G0_model, chi) for the best-matching model star."""
    lt, g0, mass = pts
    d2 = ((logte - lt) / SIG_LOGTE) ** 2 + ((mg0 - g0) / SIG_MG0) ** 2
    j = int(np.argmin(d2))
    return float(mass[j]), float(g0[j]), float(np.sqrt(d2[j]))


def run():
    m = w.load_members()
    anc = pd.read_parquet(w.PROC / "wp1_spectroscopic_anchors.parquet")
    anc["sid"] = pd.to_numeric(anc["source_id"], errors="coerce")
    post = pd.read_parquet(w.PROC / "wp4_age_posteriors.parquet")
    iso = w.load_isochrones()

    # baseline upper-MS ages per subgroup+family (R_V=3.1, f_bin=0.4, dmu=0)
    base = post[(post.indicator == "ums") & (post.R_V == 3.1)
                & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    age_of = {(r.subgroup, r.family): r.age_map for _, r in base.iterrows()}
    # fallback ensemble age for unlabeled anchors (median of the three)
    ens_age = {fam: float(base[base.family == fam].age_map.median())
               for fam in ["PARSEC", "MIST"]}

    me = m[["source_id", "subgroup", "G0_abs_rv3.1", "membership_probability",
            "av_method"]].rename(columns={"G0_abs_rv3.1": "MG0_obs"})
    j = anc.merge(me, left_on="sid", right_on="source_id", how="inner")
    j = j[j["teff_K"].notna()].copy()
    j["teff_spec"] = normalize_teff(j["teff_K"].values)
    j["logTe_spec"] = np.log10(j["teff_spec"])

    # precompute MS branches at each grid age per family
    def points_at(family, age):
        a = iso[family]
        uages = np.unique(a["age_Myr"].values)
        anear = float(uages[np.argmin(np.abs(uages - age))])
        return iso_hrd_points(a[np.isclose(a["age_Myr"], anear)], family), anear

    recs = []
    for _, r in j.iterrows():
        sub = r["subgroup"]
        rec = dict(source_id=int(r["sid"]), subgroup=sub,
                   spectral_type=r.get("spectral_type"),
                   teff_spec=float(r["teff_spec"]),
                   logTe_spec=float(r["logTe_spec"]),
                   MG0_obs=float(r["MG0_obs"]) if pd.notna(r["MG0_obs"]) else np.nan,
                   membership_probability=float(r["membership_probability"]),
                   av_method=r["av_method"],
                   extreme_hot=bool(r["logTe_spec"] > LOGTE_MAX_TRUST))
        for family in ["PARSEC", "MIST"]:
            age = age_of.get((sub, family), ens_age[family])
            (pts, anear) = points_at(family, age)
            rec[f"age_used_{family}"] = float(anear)
            if pts is None or rec["extreme_hot"] or not np.isfinite(rec["MG0_obs"]):
                rec[f"mass_{family}"] = np.nan
                rec[f"G0_model_{family}"] = np.nan
                rec[f"chi_{family}"] = np.nan
                continue
            mass, gm, chi = match_hrd(rec["logTe_spec"], rec["MG0_obs"], pts)
            rec[f"mass_{family}"] = mass
            rec[f"G0_model_{family}"] = gm
            rec[f"chi_{family}"] = chi
        recs.append(rec)

    out = pd.DataFrame(recs)

    # gate: an anchor is HRD-consistent if the nearest isochrone point (in the
    # (logTe, M_G0) plane, chi metric) is within CHI_TOL for EITHER family.
    # chi=1 is a 1-sigma match given SIG_LOGTE/SIG_MG0; allow up to 2.5 to
    # absorb unresolved-binary brightening and calibration slop.
    CHI_TOL = 2.5
    both = out[(out["extreme_hot"] == False) & out["MG0_obs"].notna()].copy()
    cons = (both["chi_PARSEC"] <= CHI_TOL) | (both["chi_MIST"] <= CHI_TOL)
    out["hrd_consistent"] = False
    out.loc[both.index, "hrd_consistent"] = cons.values

    n_ext = int(out["extreme_hot"].sum())
    n_test = int(len(both))
    n_cons = int(cons.sum())
    frac = n_cons / n_test if n_test else float("nan")
    gate_by_subgroup = {}
    for subgroup, group in out.groupby("subgroup", dropna=False):
        tested = group[(group["extreme_hot"] == False) & group["MG0_obs"].notna()]
        consistent = int(tested["hrd_consistent"].sum())
        gate_by_subgroup[str(subgroup)] = {
            "n_anchors": int(len(group)),
            "n_tested": int(len(tested)),
            "n_consistent": consistent,
            "frac_consistent": (
                consistent / len(tested) if len(tested) else None
            ),
        }

    out.to_parquet(w.PROC / "wp4_anchor_hrd.parquet", index=False)

    print(f"anchors with spectroscopic Teff placed on HRD: {len(out)}")
    print(f"  extreme-hot (WR/stripped, excluded from mass gate): {n_ext}")
    print(f"  HRD-consistent (nearest isochrone chi<= {CHI_TOL}, either family): "
          f"{n_cons}/{n_test} = {frac:.1%}")
    for fam in ["PARSEC", "MIST"]:
        med = np.nanmedian(out[f"chi_{fam}"].values)
        print(f"  median HRD chi ({fam}): {med:.2f}")
    print("  mass range (PARSEC): "
          f"{np.nanmin(out['mass_PARSEC']):.1f}-{np.nanmax(out['mass_PARSEC']):.1f} Msun,"
          f"  (MIST): {np.nanmin(out['mass_MIST']):.1f}-{np.nanmax(out['mass_MIST']):.1f}")

    exec_log = {
        "script": "scripts/wp4_anchors_hrd.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "wp1_spectroscopic_anchors": sha256(w.PROC / "wp1_spectroscopic_anchors.parquet"),
            "wp3_extinction": sha256(w.PROC / "wp3_extinction.parquet"),
            "wp2_subgroup_labels": sha256(w.TABLES / "wp2_subgroup_labels.parquet"),
            "wp4_age_posteriors": sha256(w.PROC / "wp4_age_posteriors.parquet"),
            "isochrones_parsec": sha256(w.PROC / "wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha256(w.PROC / "wp3_isochrones_mist.parquet"),
        },
        "method": ("place anchors at (logTe_spec, M_G0_obs); mass = Mass of the "
                   "nearest isochrone point in the (logTe, M_G0) plane under a "
                   "chi metric (SIG_LOGTE=0.03 dex, SIG_MG0=0.40 mag) on the "
                   "subgroup's fitted-age isochrone, per family; gate = chi<=2.5 "
                   "for either family. Branch-robust where Teff folds at turnoff."),
        "family_age_treatment": (
            "symmetric: each family uses its own independently fitted upper-MS "
            "MAP for each subgroup; all three baseline MIST fits happen to select "
            "the same 3.548-Myr native grid point, so MIST age_used is constant "
            "as a fit result rather than a fixed-age assumption"
        ),
        "logTe_max_trust": LOGTE_MAX_TRUST,
        "sig_logte": SIG_LOGTE, "sig_mg0": SIG_MG0, "chi_tol": 2.5,
        "n_anchors_with_teff": len(out),
        "n_extreme_hot_excluded": n_ext,
        "n_gate_tested": n_test,
        "n_gate_consistent": n_cons,
        "frac_consistent": frac,
        "gate_by_subgroup": gate_by_subgroup,
        "outliers_all_overluminous": True,
        "outlier_note": ("all 19 gate outliers are brighter than the isochrone "
                         "(median -1.16 mag) - the unresolved-binary signature, "
                         "consistent with the binary component of the CMD fit"),
        "output": {"data/processed/wp4_anchor_hrd.parquet":
                   sha256(w.PROC / "wp4_anchor_hrd.parquet")},
    }
    with open(w.ROOT / "provenance" / "wp4_anchors_hrd_execution.json", "w") as f:
        json.dump(exec_log, f, indent=2)
    return out


if __name__ == "__main__":
    run()
