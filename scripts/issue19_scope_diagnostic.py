#!/usr/bin/env python3
"""Issue #19 -- scope diagnostic.  Read-only: nothing is refitted or written
except the two outputs below.

Measures the two findings behind issue #19:

  F1  the quoted 2.25-5.67 Myr "two-indicator" age envelope reproduces only
      from the unversioned pre-repair wp4_age_posteriors.parquet; every
      repaired version retains zero PMS rows and repair_v5 spans 2.00-4.01;
  F2  the 150 spectroscopic-HRD anchor masses carried as
      ``spectroscopic_hrd_frozen`` into wp4_mass_posteriors_repair_v5 were read
      off isochrones at the PRE-REPAIR ages (wp4_anchors_hrd.py reads the
      unversioned posterior).  They are re-derived here with the same
      procedure at the repair_v5 ages, like for like, to count how many cross
      the 2 and 8 Msun boundaries that the WP5 window and WP6 census use.

F2's re-derivation is a scope estimate, not a replacement product: it uses the
R_V = 3.1 baseline ages for every branch, exactly as the frozen procedure did.

Outputs:
  tables/issue19_envelope_by_version.csv
  tables/issue19_anchor_mass_rederivation.csv

Run:
  PYTHONPATH=scripts python3 scripts/issue19_scope_diagnostic.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import wp4_anchors_hrd as H
import wp5_common as w
from wp10_inputs import resolve

VERSIONS = ("", "_repair_v1", "_repair_v3", "_repair_v5")


def envelope_by_version() -> pd.DataFrame:
    rows = []
    for suffix in VERSIONS:
        post = pd.read_parquet(w.PROC / f"wp4_age_posteriors{suffix}.parquet")
        kept = post[post.measurable.astype(bool) & ~post.grid_railed.astype(bool)]
        rows.append({
            "version": suffix.lstrip("_") or "unversioned_pre_repair",
            "retained_rows": len(kept),
            "retained_pms_rows": int(kept.indicator.eq("pms").sum()),
            "map_min_Myr": round(float(kept.age_map.min()), 2),
            "map_max_Myr": round(float(kept.age_map.max()), 2),
        })
    return pd.DataFrame(rows)


def baseline_ages(post: pd.DataFrame) -> dict:
    base = post[(post.indicator == "ums") & (post.R_V == 3.1)
                & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    return {(r.subgroup, r.family): float(r.age_map) for _, r in base.iterrows()}


def anchor_rederivation() -> pd.DataFrame:
    ages = baseline_ages(pd.read_parquet(resolve("wp4_age_posteriors")))
    fallback = {fam: float(np.median([a for (_, f), a in ages.items() if f == fam]))
                for fam in ("PARSEC", "MIST")}
    hrd = pd.read_parquet(w.PROC / "wp4_anchor_hrd.parquet")
    mg0 = pd.read_parquet(resolve("wp3_extinction"))[["source_id", "G0_abs_rv3.1"]]
    hrd = hrd.merge(mg0, on="source_id", how="left")
    iso = H.w.load_isochrones()
    rows = []
    for _, r in hrd.iterrows():
        mag = r["G0_abs_rv3.1"] if np.isfinite(r["G0_abs_rv3.1"]) else r.MG0_obs
        if r.extreme_hot or not np.isfinite(mag):
            continue
        for fam in ("PARSEC", "MIST"):
            age = ages.get((r.subgroup, fam), fallback[fam])
            table = iso[fam]
            native = np.unique(table.age_Myr.values)
            near = float(native[np.argmin(np.abs(native - age))])
            pts = H.iso_hrd_points(table[np.isclose(table.age_Myr, near)], fam)
            new = H.match_hrd(r.logTe_spec, mag, pts)[0] if pts is not None else np.nan
            rows.append({
                "source_id": int(r.source_id),
                "subgroup": r.subgroup if isinstance(r.subgroup, str) else "unassigned",
                "family": fam,
                "membership_probability": float(r.membership_probability),
                "age_frozen_Myr": float(r[f"age_used_{fam}"]),
                "age_repair_v5_Myr": near,
                "mass_frozen": float(r[f"mass_{fam}"]),
                "mass_rederived": float(new),
            })
    out = pd.DataFrame(rows).dropna()
    out["crosses_8"] = (out.mass_frozen >= 8) != (out.mass_rederived >= 8)
    out["crosses_2"] = (out.mass_frozen >= 2) != (out.mass_rederived >= 2)
    return out


def main() -> None:
    env = envelope_by_version()
    env.to_csv(w.TABLES / "issue19_envelope_by_version.csv", index=False)
    print(env.to_string(index=False))

    anchors = anchor_rederivation()
    anchors.to_csv(w.TABLES / "issue19_anchor_mass_rederivation.csv", index=False)
    summary = anchors.groupby(["family", "subgroup"]).agg(
        n=("source_id", "size"),
        median_ratio=("mass_rederived", lambda s: float(np.median(
            s / anchors.loc[s.index, "mass_frozen"]))),
        cross_8=("crosses_8", "sum"),
        cross_2=("crosses_2", "sum"),
    ).round(3)
    print(summary.to_string())


if __name__ == "__main__":
    main()
