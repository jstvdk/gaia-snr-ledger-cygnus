#!/usr/bin/env python3
"""Issue #20 Phase A, step A2 -- star-by-star audit of CygOB2-C's 43
spectroscopic stars and its 15 brightest members (repair_v5 M_G0, R_V 3.1).

Per star: spectral type, luminosity class, composite flag, the A3 T_eff and its
source, RUWE, literature binarity (SIMBAD object types SB*/EB*/El*, cached),
membership probability, proper-motion and sky Mahalanobis distances to the A,
B and C centroids, spectroscopic (adopted) vs literature vs pre-repair A_V,
literature identifiers.  Scores secondary prediction S1's measured quantity.

Centroids: membership-weighted mean and covariance of (pmra, pmdec) and of
(l, b) over the 1,331 labelled members of each subgroup.  The labels come from
a GMM in (l, b, pmra, pmdec), so distances to C are small for C members by
construction (S1's registered caveat).

SIMBAD is queried once by Gaia DR3 identifier; the response is cached to
data/raw/issue20_simbad_otypes.csv and hashed.  If SIMBAD is unreachable the
column is left empty and the execution record says so.

Outputs:
  tables/issue20_c_star_audit.csv
  provenance/issue20_star_audit_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20_star_audit.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I
import wp5_common as w

OUT = w.TABLES / "issue20_c_star_audit.csv"
OUT_EXEC = w.PROVENANCE / "issue20_star_audit_execution.json"
SIMBAD_CACHE = w.ROOT / "data" / "raw" / "issue20_simbad_otypes.csv"
BINARY_OTYPES = {"SB*", "EB*", "El*", "**"}


def centroids(members: pd.DataFrame) -> dict:
    out = {}
    for sg in I.SUBGROUPS:
        g = members[members.subgroup.eq(sg)]
        wts = g.membership_probability.to_numpy(float)
        for key, cols in (("pm", ["pmra", "pmdec"]), ("sky", ["l_deg", "b_deg"])):
            x = g[cols].to_numpy(float)
            mu = np.average(x, axis=0, weights=wts)
            cov = np.cov(x.T, aweights=wts)
            out[(sg, key)] = (mu, np.linalg.inv(cov))
    return out


def mahalanobis(x, mu, icov):
    d = x - mu
    return float(np.sqrt(d @ icov @ d))


def simbad_otypes(ids: list[int]) -> tuple[pd.DataFrame, str]:
    if SIMBAD_CACHE.exists():
        return pd.read_csv(SIMBAD_CACHE), "cache"
    try:
        from astroquery.simbad import Simbad
        s = Simbad()
        s.add_votable_fields("otypes")
        r = s.query_objects([f"Gaia DR3 {i}" for i in ids]).to_pandas()
        r["source_id"] = r["user_specified_id"].str.replace("Gaia DR3 ", "").astype("int64")
        table = (r.groupby("source_id")
                 .agg(simbad_main_id=("main_id", "first"),
                      simbad_otypes=("otypes.otype",
                                     lambda v: "|".join(sorted({str(x) for x in v if pd.notna(x)}))))
                 .reset_index())
        table["queried_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        SIMBAD_CACHE.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(SIMBAD_CACHE, index=False)
        return table, "queried"
    except Exception as exc:          # network or service failure: record, do not guess
        return pd.DataFrame({"source_id": ids}), f"unavailable: {exc!r}"


def main() -> None:
    method = I.check_method_unchanged()
    anchors = pd.read_parquet(I.frozen("anchor_hrd"))
    lit = pd.read_parquet(I.frozen("spectroscopic_anchors"))
    lit["source_id"] = pd.to_numeric(lit["source_id"], errors="coerce")
    lit = lit.dropna(subset=["source_id"]).astype({"source_id": "int64"})
    ext5 = pd.read_parquet(I.frozen("extinction_v5"))
    ext_pre = pd.read_parquet(I.frozen("extinction_pre_repair"))
    labels = pd.read_parquet(I.frozen("subgroup_labels"))
    mem = pd.read_parquet(I.frozen("members")).drop(columns=["subgroup"], errors="ignore")
    mem = mem.merge(labels, on="source_id", how="inner")
    cents = centroids(mem)

    c_anchor = set(anchors[anchors.subgroup.eq("CygOB2-C")].source_id)
    c5 = ext5.drop(columns=["subgroup"], errors="ignore").merge(labels, on="source_id")
    c5 = c5[c5.subgroup.eq("CygOB2-C")]
    bright = c5.nsmallest(15, "G0_abs_rv3.1").source_id.tolist()
    ids = sorted(c_anchor | set(bright))
    otypes, simbad_status = simbad_otypes(ids)

    rows = []
    for sid in ids:
        e = ext5[ext5.source_id.eq(sid)].iloc[0]
        m = mem[mem.source_id.eq(sid)].iloc[0]
        a = anchors[anchors.source_id.eq(sid)]
        l = lit[lit.source_id.eq(sid)]
        st = a.spectral_type.iloc[0] if len(a) else (l.spectral_type.iloc[0] if len(l) else None)
        rec = {"source_id": sid, "in_spectroscopic_set": sid in c_anchor,
               "rank_brightest_v5": (bright.index(sid) + 1) if sid in bright else None,
               "object_name": l.object_name.iloc[0] if len(l) else None,
               "aliases": l.object_aliases_json.iloc[0] if len(l) else None,
               "spectral_type": st, "ra": e.ra, "dec": e.dec, "l_deg": e.l_deg, "b_deg": e.b_deg,
               "G": e.G, "MG0_rv3.1": e["G0_abs_rv3.1"], "BPRP0_rv3.1": e["BPRP0_rv3.1"],
               "ruwe": e.ruwe, "ruwe_gt_1p4": bool(e.ruwe > 1.4),
               "membership_probability": e.membership_probability,
               "pmra": m.pmra, "pmdec": m.pmdec,
               "av_adopted_rv3.1": e["av_rv3.1"], "av_method": e.av_method,
               "av_pre_repair_rv3.1": float(ext_pre[ext_pre.source_id.eq(sid)]["av_rv3.1"].iloc[0]),
               "av_literature": l.extinction_av_mag.iloc[0] if len(l) else np.nan,
               "teff_table": a.teff_spec.iloc[0] if len(a) else np.nan,
               "logg_literature": l.logg_cgs.iloc[0] if len(l) else np.nan}
        parsed = I.parse_spectral_type(st)
        rec.update(lum_class=parsed["lum_class"], composite_type=parsed["composite"])
        if len(a):
            t = I.assign_teff(st, a.teff_spec.iloc[0])
            rec.update(teff_used=t["teff_used"], teff_source=t["teff_source"])
        for sg in I.SUBGROUPS:
            tag = sg[-1]
            mu, ic = cents[(sg, "pm")]
            rec[f"pm_mahal_{tag}"] = mahalanobis(np.array([m.pmra, m.pmdec]), mu, ic)
            mu, ic = cents[(sg, "sky")]
            rec[f"sky_mahal_{tag}"] = mahalanobis(np.array([e.l_deg, e.b_deg]), mu, ic)
        rec["pm_nearest"] = min("ABC", key=lambda t: rec[f"pm_mahal_{t}"])
        rec["sky_nearest"] = min("ABC", key=lambda t: rec[f"sky_mahal_{t}"])
        o = otypes[otypes.source_id.eq(sid)]
        rec["simbad_main_id"] = o.simbad_main_id.iloc[0] if len(o) and "simbad_main_id" in o else None
        ot = o.simbad_otypes.iloc[0] if len(o) and "simbad_otypes" in o else None
        rec["simbad_otypes"] = ot
        rec["simbad_binary_flag"] = (bool(set(str(ot).split("|")) & BINARY_OTYPES)
                                     if isinstance(ot, str) else None)
        rows.append(rec)
    table = pd.DataFrame(rows).sort_values("MG0_rv3.1")
    table.to_csv(OUT, index=False)

    sg = table[table.in_spectroscopic_set & table.lum_class.isin(["I", "II"])]
    frac_c = float((sg.pm_nearest == "C").mean()) if len(sg) else np.nan
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_star_audit.py",
        "preregistration": "provenance/issue20_prereg.json",
        "method_check": method,
        "inputs": {n: I.prereg()["consumed_inputs"][n] for n in
                   ("anchor_hrd", "spectroscopic_anchors", "extinction_v5",
                    "extinction_pre_repair", "subgroup_labels", "members")},
        "simbad": {"status": simbad_status,
                   "cache": str(SIMBAD_CACHE.relative_to(w.ROOT)),
                   "cache_sha256": w.sha256(SIMBAD_CACHE) if SIMBAD_CACHE.exists() else None},
        "n_rows": int(len(table)), "n_spectroscopic": int(table.in_spectroscopic_set.sum()),
        "n_brightest_not_spectroscopic": int((~table.in_spectroscopic_set).sum()),
        "S1_measured": {"class_I_II_stars": sg[["source_id", "object_name", "spectral_type",
                                                 "pm_mahal_A", "pm_mahal_B", "pm_mahal_C",
                                                 "pm_nearest"]].to_dict(orient="records"),
                        "fraction_nearest_C": frac_c},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    }
    w.write_json(OUT_EXEC, record)
    print(f"{len(table)} stars; SIMBAD {simbad_status}; S1 fraction nearest C = {frac_c:.2f} "
          f"({len(sg)} class I/II)")


if __name__ == "__main__":
    main()
