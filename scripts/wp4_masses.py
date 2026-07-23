#!/usr/bin/env python3
"""WP4 step 5 - per-star masses.

For every member we read a mass off the subgroup's fitted-age isochrone, per
branch (family x R_V), by nearest-point matching in the de-reddened CMD plane
(colour, M_G0).  For the spectroscopic anchors we OVERRIDE the photometric mass
with the spectroscopic-HRD mass (wp4_anchors_hrd), which breaks the hot-star
colour degeneracy.  Every star carries a provenance flag recording which method
set its mass and which subgroup age was used.

Mass columns: mass_{PARSEC,MIST}_rv{3.0,3.1,3.5}  (6 per star), plus a reporting
mass_baseline (PARSEC R_V=3.1) and the PARSEC-MIST spread at baseline.

Output: data/processed/wp4_masses.parquet, provenance/wp4_masses_execution.json
"""
from __future__ import annotations

import json, hashlib, datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w

SIG_COL = 0.10    # CMD-match colour scale (Class B)
SIG_MG = 0.30     # CMD-match magnitude scale (Class B)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iso_cmd_points(iso, family, age):
    """(colour, M_G0, Mass) for MS+near-MS points at the nearest grid age."""
    a = iso[family]
    uages = np.unique(a["age_Myr"].values)
    anear = float(uages[np.argmin(np.abs(uages - age))])
    d = a[np.isclose(a["age_Myr"], anear)]
    d = d[d["label"] <= 2] if family == "PARSEC" else d[d["phase"] <= 2]
    col = (d["BP0"] - d["RP0"]).values
    mg = d["G0"].values
    mass = d["Mass"].values
    ok = np.isfinite(col) & np.isfinite(mg) & np.isfinite(mass)
    return col[ok], mg[ok], mass[ok], anear


def match_mass(col, mg, pts):
    c, g, mass = pts
    d2 = ((col - c) / SIG_COL) ** 2 + ((mg - g) / SIG_MG) ** 2
    j = int(np.argmin(d2))
    return float(mass[j]), float(np.sqrt(d2[j]))


def run():
    m = w.load_members()
    iso = w.load_isochrones()
    post = pd.read_parquet(w.PROC / "wp4_age_posteriors.parquet")
    anchor = pd.read_parquet(w.PROC / "wp4_anchor_hrd.parquet")

    # subgroup fitted upper-MS ages per (family, R_V) at f_bin=0.4, dmu=0
    ubase = post[(post.indicator == "ums") & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    def age_lookup(sub, fam, rv):
        r = ubase[(ubase.subgroup == sub) & (ubase.family == fam) & (ubase.R_V == rv)]
        if len(r):
            return float(r.age_map.iloc[0])
        return float(ubase[(ubase.family == fam) & (ubase.R_V == rv)].age_map.median())

    branches = [(fam, rv) for fam in ["PARSEC", "MIST"] for rv in w.R_V_BRANCHES]

    # anchor spectroscopic masses keyed by source_id (family-specific, RV-indep)
    anc_mass = {int(r.source_id): {"PARSEC": r.mass_PARSEC, "MIST": r.mass_MIST}
                for _, r in anchor.iterrows()}
    anc_extreme = {int(r.source_id): bool(r.extreme_hot) for _, r in anchor.iterrows()}

    rows = []
    fits = {rv: w.branch_photometry(m, rv) for rv in w.R_V_BRANCHES}
    # precompute isochrone points per (family, rv, subgroup-age)
    for i in range(len(m)):
        sid = int(m["source_id"].iloc[i])
        sub = m["subgroup"].iloc[i]
        rec = {"source_id": sid, "subgroup": sub,
               "membership_probability": float(m["membership_probability"].iloc[i]),
               "av_method": m["av_method"].iloc[i]}
        is_anchor = sid in anc_mass
        spec_ok = is_anchor and not anc_extreme.get(sid, False) \
            and np.isfinite(anc_mass[sid]["PARSEC"])
        rec["mass_method"] = "spectroscopic_hrd" if spec_ok else "photometric_isochrone"
        rec["is_spectroscopic_anchor"] = bool(is_anchor)
        for (fam, rv) in branches:
            colname = f"mass_{fam}_rv{rv:.1f}"
            if spec_ok:
                rec[colname] = float(anc_mass[sid][fam])
                rec[f"chi_{fam}_rv{rv:.1f}"] = np.nan
                continue
            frow = fits[rv].iloc[i]
            col, mg = frow["colour"], frow["MG0"]
            if not (np.isfinite(col) and np.isfinite(mg)):
                rec[colname] = np.nan
                rec[f"chi_{fam}_rv{rv:.1f}"] = np.nan
                continue
            age = age_lookup(sub, fam, rv)
            c, g, mass, anear = iso_cmd_points(iso, fam, age)
            mm, chi = match_mass(col, mg, (c, g, mass))
            rec[colname] = mm
            rec[f"chi_{fam}_rv{rv:.1f}"] = chi
        rows.append(rec)

    df = pd.DataFrame(rows)
    # headline mass + PARSEC-MIST model spread at baseline R_V=3.1
    df["mass_baseline"] = df["mass_PARSEC_rv3.1"]
    df["mass_parsec_mist_spread"] = (df["mass_PARSEC_rv3.1"] - df["mass_MIST_rv3.1"]).abs()

    df.to_parquet(w.PROC / "wp4_masses.parquet", index=False)

    n_spec = int((df.mass_method == "spectroscopic_hrd").sum())
    n_phot = int((df.mass_method == "photometric_isochrone").sum())
    print(f"per-star masses written for {len(df)} members")
    print(f"  spectroscopic-HRD masses: {n_spec}   photometric-isochrone: {n_phot}")
    for (fam, rv) in [("PARSEC", 3.1), ("MIST", 3.1)]:
        c = f"mass_{fam}_rv{rv:.1f}"
        v = df[c].dropna()
        print(f"  {c}: median={v.median():.2f}  range [{v.min():.2f},{v.max():.2f}] Msun")
    print(f"  median |PARSEC-MIST| mass spread (rv3.1): "
          f"{df['mass_parsec_mist_spread'].median():.2f} Msun")
    # calibration-window (2-8 Msun) count preview for WP5
    inwin = df["mass_baseline"].between(2, 8)
    print(f"  members with 2-8 Msun (baseline): {int(inwin.sum())} "
          f"(weighted {float(df.loc[inwin,'membership_probability'].sum()):.0f})")

    exec_log = {
        "script": "scripts/wp4_masses.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "wp3_extinction": sha256(w.PROC / "wp3_extinction.parquet"),
            "wp2_subgroup_labels": sha256(w.TABLES / "wp2_subgroup_labels.parquet"),
            "wp4_age_posteriors": sha256(w.PROC / "wp4_age_posteriors.parquet"),
            "wp4_anchor_hrd": sha256(w.PROC / "wp4_anchor_hrd.parquet"),
            "isochrones_parsec": sha256(w.PROC / "wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha256(w.PROC / "wp3_isochrones_mist.parquet"),
        },
        "method": ("nearest-point mass in (colour, M_G0) on subgroup fitted-age "
                   "isochrone per (family, R_V); spectroscopic-HRD override for "
                   "anchors; 6 branch mass columns per star + provenance flag"),
        "sig_col": SIG_COL, "sig_mg": SIG_MG,
        "n_members": len(df), "n_spectroscopic": n_spec, "n_photometric": n_phot,
        "median_parsec_mist_mass_spread_rv31": float(df["mass_parsec_mist_spread"].median()),
        "output": {"data/processed/wp4_masses.parquet":
                   sha256(w.PROC / "wp4_masses.parquet")},
    }
    with open(w.ROOT / "provenance" / "wp4_masses_execution.json", "w") as f:
        json.dump(exec_log, f, indent=2)
    return df


if __name__ == "__main__":
    run()
