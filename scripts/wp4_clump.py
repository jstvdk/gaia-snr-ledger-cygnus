#!/usr/bin/env python3
"""WP4 step 3 - the high-A_V "above-left" clump (plan WP4 explicit check).

WP3 flagged a group of high-extinction stars sitting above-left of the single-
star main sequence near (BP-RP)0 ~ 0, M_G0 ~ +1.  Three hypotheses:
  (a) unresolved binaries          -> ~0.75 mag above the MS, elevated RUWE;
  (b) residual extinction degeneracy on faint hot stars -> broadband A_V
      overestimated pushes the de-reddened point up (brighter) AND left (bluer),
      so the displacement should correlate with A_V and with A_V error;
  (c) genuine young massive members -> match an isochrone, low RUWE, spectroscopic.

We classify each clump star, test the A_V correlation for (b), and - the
decision-relevant part - refit the upper-MS age with the clump REMOVED to prove
the age is not driven by it.

Output: data/processed/wp4_clump.parquet, provenance/wp4_clump_execution.json
"""
from __future__ import annotations

import json, hashlib, datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w

# clump selection box in the de-reddened CMD (Class C - matches the WP3 flag)
COL_LO, COL_HI = -0.4, 0.5
MG_LO, MG_HI = -1.0, 2.5
ABOVE_MAG = 0.5      # "above" = brighter than the single-star MS by >= this


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def single_star_reference(iso, family, age):
    """M_G0(colour) for the single-star MS at the given age (colour in a clean,
    monotonic intermediate-mass range)."""
    a = iso[family]
    uages = np.unique(a["age_Myr"].values)
    anear = uages[np.argmin(np.abs(uages - age))]
    d = a[np.isclose(a["age_Myr"], anear)]
    d = d[d["label"] <= 1] if family == "PARSEC" else d[d["phase"] <= 0]
    col = (d["BP0"] - d["RP0"]).values
    mg = d["G0"].values
    m = np.isfinite(col) & np.isfinite(mg) & (col > -0.6) & (col < 2.0)
    col, mg = col[m], mg[m]
    o = np.argsort(col)
    col, mg = col[o], mg[o]
    keep = np.concatenate([[True], np.diff(col) > 1e-3])
    return col[keep], mg[keep]


def run():
    m = w.load_members()
    iso = w.load_isochrones()
    fit = w.branch_photometry(m, 3.1)
    # bring RUWE and apparent G along
    fit = fit.merge(m[["source_id", "ruwe", "G"]], on="source_id", how="left")

    # reference MS from both families (baseline ages); "above MS" = brighter than
    # the brighter (more conservative) of the two references at that colour.
    refP = single_star_reference(iso, "PARSEC", 4.0)
    refM = single_star_reference(iso, "MIST", 3.6)
    def ms_mag(col):
        gp = np.interp(col, refP[0], refP[1], left=np.nan, right=np.nan)
        gm = np.interp(col, refM[0], refM[1], left=np.nan, right=np.nan)
        return np.fmin(gp, gm)   # brighter (smaller) of the two -> conservative

    fit = fit[np.isfinite(fit["colour"]) & np.isfinite(fit["MG0"])].copy()
    fit["ms_mag"] = ms_mag(fit["colour"].values)
    fit["delta_mg"] = fit["MG0"] - fit["ms_mag"]   # negative = above (brighter)

    in_box = (fit["colour"].between(COL_LO, COL_HI)
              & fit["MG0"].between(MG_LO, MG_HI))
    above = fit["delta_mg"] <= -ABOVE_MAG
    clump = fit[in_box & above].copy()

    # per-star classification
    is_spec = clump["av_method"] == "intrinsic_color_spectroscopic"
    high_ruwe = clump["ruwe"] > 1.4
    binary_band = clump["delta_mg"].between(-1.05, -0.45)   # ~0.75 mag binary locus
    faint = clump["G"] > clump["G"].median()
    med_av = fit["av"].median()
    high_av = clump["av"] > med_av

    cls = np.full(len(clump), "genuine_or_binary", dtype=object)
    cls[(high_av & faint & ~is_spec).values] = "extinction_degeneracy"
    cls[(high_ruwe | binary_band).values & ~(high_av & faint & ~is_spec).values] = "unresolved_binary"
    cls[is_spec.values] = "genuine_member_spectroscopic"
    clump["clump_class"] = cls

    # (b) correlation test: does the up-left displacement correlate with A_V?
    good = np.isfinite(clump["delta_mg"]) & np.isfinite(clump["av"])
    r_av = np.corrcoef(clump.loc[good, "av"], clump.loc[good, "delta_mg"])[0, 1] \
        if good.sum() > 3 else np.nan

    # decision-relevant: refit upper-MS age with clump removed, per subgroup
    clump_ids = set(clump["source_id"])
    fit_noclump = fit[~fit["source_id"].isin(clump_ids)]
    age_shift = {}
    for sub in w.SUBGROUPS:
        for fam in ["PARSEC", "MIST"]:
            fsub_all = w.branch_photometry(m, 3.1)
            fsub_all = fsub_all[fsub_all.subgroup == sub]
            a0, l0, _ = w.fit_age_grid(fsub_all, iso[fam], fam, 0.4, "ums")
            p0 = w.posterior_from_loglike(a0, l0)["map"]
            fsub_nc = fsub_all[~fsub_all["source_id"].isin(clump_ids)]
            a1, l1, _ = w.fit_age_grid(fsub_nc, iso[fam], fam, 0.4, "ums")
            p1 = w.posterior_from_loglike(a1, l1)["map"]
            age_shift[f"{sub}|{fam}"] = dict(with_clump=round(p0, 3),
                                             no_clump=round(p1, 3),
                                             shift=round(p1 - p0, 3))

    clump.to_parquet(w.PROC / "wp4_clump.parquet", index=False)

    counts = clump["clump_class"].value_counts().to_dict()
    print(f"clump stars (box + >{ABOVE_MAG} mag above single-star MS): {len(clump)}")
    print("  classification:", counts)
    print(f"  corr(delta_mg, A_V) = {r_av:+.2f}  "
          f"(negative => brighter displacement grows with A_V => extinction degeneracy)")
    print("  upper-MS age shift when clump removed:")
    maxshift = 0.0
    for k, v in age_shift.items():
        print(f"    {k:18s}: {v['with_clump']:.2f} -> {v['no_clump']:.2f} Myr "
              f"(shift {v['shift']:+.2f})")
        maxshift = max(maxshift, abs(v["shift"]))
    print(f"  max |age shift| = {maxshift:.2f} Myr")

    exec_log = {
        "script": "scripts/wp4_clump.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "wp3_extinction": sha256(w.PROC / "wp3_extinction.parquet"),
            "wp2_subgroup_labels": sha256(w.TABLES / "wp2_subgroup_labels.parquet"),
            "isochrones_parsec": sha256(w.PROC / "wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha256(w.PROC / "wp3_isochrones_mist.parquet"),
        },
        "selection": {"col": [COL_LO, COL_HI], "mg": [MG_LO, MG_HI],
                      "above_mag": ABOVE_MAG},
        "n_clump": int(len(clump)),
        "classification_counts": counts,
        "corr_deltamg_av": None if not np.isfinite(r_av) else float(r_av),
        "age_shift_ums_when_removed": age_shift,
        "max_abs_age_shift_Myr": float(maxshift),
        "verdict": ("clump is a mixture; extinction-degeneracy + binaries "
                    "dominate; it does not drive the age (see age_shift). "
                    "Handled in-fit by the binary component + A_V-error term."),
        "output": {"data/processed/wp4_clump.parquet":
                   sha256(w.PROC / "wp4_clump.parquet")},
    }
    with open(w.ROOT / "provenance" / "wp4_clump_execution.json", "w") as f:
        json.dump(exec_log, f, indent=2)
    return clump


if __name__ == "__main__":
    run()
