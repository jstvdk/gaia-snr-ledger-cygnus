#!/usr/bin/env python3
"""Issue #20 -- is CygOB2-C really young?  Read-only scope diagnostic.

Three measurements, nothing refitted:

  D1  C's upper-MS photometric age (PARSEC, R_V = 3.1, f_bin = 0.4, dmu = 0)
      in every WP4 version, to show when the young age appeared;
  D2  the brightest members of each subgroup with their spectral types: C's
      are luminosity-class I supergiants, i.e. evolved stars that a colour--
      magnitude fit cannot tell from main-sequence stars;
  D3  an age scan of the spectroscopic HRD positions (spectroscopic Teff,
      repair_v5 M_G0) against PARSEC and MIST isochrones, with the project's
      own nearest-point metric from wp4_anchors_hrd.py (SIG_LOGTE = 0.03,
      SIG_MG0 = 0.40).  Sum of chi^2 per age, all anchors and the bright
      (M_G0 < -5) ones.

D3 is a diagnostic, not an age measurement: no IMF weighting, no explicit
binary model, an anchor sample biased toward luminous stars, and reduced
chi^2 of ~3-4, so Delta chi^2 overstates significance.  The spectroscopic
Teff are read from the unversioned wp4_anchor_hrd.parquet; they depend only
on spectral types, not on any repaired product.

Outputs:
  tables/issue20_c_age_by_version.csv
  tables/issue20_bright_members.csv
  tables/issue20_hrd_age_scan.csv

Run:
  PYTHONPATH=scripts python3 scripts/issue20_c_age_diagnostic.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import wp4_anchors_hrd as H
import wp5_common as w
import wp12_common as W
from wp10_inputs import resolve
from wp12_figures import with_labels

VERSIONS = ("", "_repair_v1", "_repair_v3", "_repair_v5")
GRID = {
    "PARSEC": (2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01),
    "MIST": (2.00, 2.25, 2.52, 2.83, 3.18, 3.57, 4.01, 4.50, 5.04),
}
BRIGHT = -5.0


def age_by_version() -> pd.DataFrame:
    rows = []
    for suffix in VERSIONS:
        post = pd.read_parquet(w.PROC / f"wp4_age_posteriors{suffix}.parquet")
        base = post[(post.indicator == "ums") & (post.R_V == 3.1)
                    & (post.f_bin == 0.4) & (post.dmu == 0.0)]
        for _, r in base.iterrows():
            rows.append({"version": suffix.lstrip("_") or "unversioned_pre_repair",
                         "subgroup": r.subgroup, "family": r.family,
                         "ums_age_map_Myr": round(float(r.age_map), 2)})
    return pd.DataFrame(rows)


def anchors() -> pd.DataFrame:
    hrd = pd.read_parquet(w.PROC / "wp4_anchor_hrd.parquet")
    mg0 = pd.read_parquet(resolve("wp3_extinction"))[["source_id", "G0_abs_rv3.1"]]
    return hrd.merge(mg0, on="source_id")


def bright_members() -> pd.DataFrame:
    photo = with_labels(pd.read_parquet(W.frozen("wp3_extinction")))
    spec = pd.read_parquet(w.PROC / "wp4_anchor_hrd.parquet")[
        ["source_id", "spectral_type", "teff_spec"]]
    m = photo[photo.sg.notna()].merge(spec, on="source_id", how="left")
    cols = ["source_id", "sg", "membership_probability", "G0_abs_rv3.1",
            "BPRP0_rv3.1", "av_rv3.1", "ruwe", "spectral_type", "teff_spec"]
    return (m.sort_values("G0_abs_rv3.1").groupby("sg").head(8)[cols]
            .sort_values(["sg", "G0_abs_rv3.1"]))


def hrd_age_scan() -> pd.DataFrame:
    hrd = anchors()
    hrd = hrd[~hrd.extreme_hot & hrd["G0_abs_rv3.1"].notna()]
    iso = H.w.load_isochrones()
    rows = []
    for subgroup in w.SUBGROUPS:
        stars = hrd[hrd.subgroup.eq(subgroup)]
        bright = stars["G0_abs_rv3.1"] < BRIGHT
        for family, ages in GRID.items():
            table = iso[family]
            native = np.unique(table.age_Myr.values)
            for age in ages:
                near = float(native[np.argmin(np.abs(native - age))])
                pts = H.iso_hrd_points(table[np.isclose(table.age_Myr, near)], family)
                chi = np.array([H.match_hrd(r.logTe_spec, r["G0_abs_rv3.1"], pts)[2]
                                for _, r in stars.iterrows()])
                rows.append({"subgroup": subgroup, "family": family,
                             "age_Myr": round(near, 2), "n": len(stars),
                             "n_bright": int(bright.sum()),
                             "chi2_all": round(float(np.sum(chi ** 2)), 1),
                             "chi2_bright": round(float(np.sum(chi[bright.values] ** 2)), 1)})
    return pd.DataFrame(rows)


def main() -> None:
    by_version = age_by_version()
    by_version.to_csv(w.TABLES / "issue20_c_age_by_version.csv", index=False)
    print(by_version.pivot_table(index=["subgroup", "family"], columns="version",
                                 values="ums_age_map_Myr").to_string())

    bright = bright_members()
    bright.to_csv(w.TABLES / "issue20_bright_members.csv", index=False)
    print(bright[bright.sg.eq("CygOB2-C")].round(2).to_string(index=False))

    scan = hrd_age_scan()
    scan.to_csv(w.TABLES / "issue20_hrd_age_scan.csv", index=False)
    best = scan.loc[scan.groupby(["subgroup", "family"]).chi2_all.idxmin()]
    print(best.to_string(index=False))


if __name__ == "__main__":
    main()
