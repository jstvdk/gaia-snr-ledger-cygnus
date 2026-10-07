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

Issue #19 (repair_v8).  The original run read the unversioned PRE-REPAIR
wp4_age_posteriors.parquet and the unversioned WP3 extinction, and assigned one
mass per family to every R_V branch.  Those inputs are now explicit and
versioned, the unversioned age posterior is refused, and every (family, R_V)
branch gets its own mass, read at that branch's own upper-MS age with that
branch's own de-reddened M_G0 (decision D1, F3(b)).  The mass is still the
isochrone's present-day ``Mass``; switching to ``Mini`` (F3(a)) is a separately
pre-registered change.  The legacy single-branch columns (mass_PARSEC, ...,
MG0_obs) are kept and equal the R_V = 3.1 baseline branch.

Output: data/processed/wp4_anchor_hrd_<version>.parquet,
        provenance/wp4_anchors_hrd_execution_<version>.json

Run (repair_v8):
  PYTHONPATH=scripts python3 scripts/wp4_anchors_hrd.py \
      --age-version repair_v5 --extinction-version repair_v5 \
      --output-version repair_v8
"""
from __future__ import annotations

import argparse
import json, hashlib, datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w
from wp3_common import normalize_teff
from wp3_extinction_law import R_V_BRANCHES

# hottest Teff we trust a normal-star isochrone mass for; hotter anchors are
# WR / stripped / extreme and are flagged out of the isochrone-mass gate.
LOGTE_MAX_TRUST = np.log10(52000.0)   # ~O2/O3; Class B

FAMILIES = ["PARSEC", "MIST"]
BASELINE_RV = 3.1


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


def branch(family: str, rv: float) -> str:
    return f"{family}_rv{rv:.1f}"


def versioned(stem: str, version: str):
    if not version:
        raise SystemExit(
            f"{stem}: an explicit version is required.  The unversioned "
            f"{stem}.parquet is the pre-repair run of 2026-07-23 (issue #19) "
            "and must not be consumed."
        )
    path = w.PROC / f"{stem}_{version}.parquet"
    if not path.exists():
        raise SystemExit(f"missing input {path.relative_to(w.ROOT)}")
    return path


def load_members(extinction_path) -> pd.DataFrame:
    """Versioned WP3 extinction joined to the canonical A/B/C labels."""
    ext = pd.read_parquet(extinction_path)
    lab = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    ext = ext.drop(columns=["subgroup", "subgroup_label"], errors="ignore")
    m = ext.merge(lab[["source_id", "subgroup"]], on="source_id", how="left")
    m["subgroup"] = m["subgroup"].fillna("unassigned")
    return m


def run(age_version: str, extinction_version: str, output_version: str,
        preregistration: str = "provenance/issue19_repair_v8_prereg.json"):
    age_path = versioned("wp4_age_posteriors", age_version)
    ext_path = versioned("wp3_extinction", extinction_version)
    if not output_version:
        raise SystemExit("--output-version is required")
    out_path = w.PROC / f"wp4_anchor_hrd_{output_version}.parquet"
    if out_path.exists():
        raise SystemExit(f"{out_path.relative_to(w.ROOT)} exists; nothing is overwritten")

    m = load_members(ext_path)
    anc = pd.read_parquet(w.PROC / "wp1_spectroscopic_anchors.parquet")
    anc["sid"] = pd.to_numeric(anc["source_id"], errors="coerce")
    post = pd.read_parquet(age_path)
    iso = w.load_isochrones()

    # upper-MS MAP per subgroup, family AND R_V branch (f_bin=0.4, dmu=0);
    # fallback for unlabelled anchors: median of the three subgroup MAPs of
    # that family and R_V.
    base = post[(post.indicator == "ums") & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    age_of = {(r.subgroup, r.family, float(r.R_V)): float(r.age_map)
              for _, r in base.iterrows()}
    ens_age = {
        (fam, float(rv)): float(base[(base.family == fam) & (base.R_V == rv)].age_map.median())
        for fam in FAMILIES for rv in R_V_BRANCHES
    }

    mag_cols = {rv: f"G0_abs_rv{rv:.1f}" for rv in R_V_BRANCHES}
    me = m[["source_id", "subgroup", "membership_probability", "av_method",
            *mag_cols.values()]]
    j = anc.merge(me, left_on="sid", right_on="source_id", how="inner")
    j = j[j["teff_K"].notna()].copy()
    j["teff_spec"] = normalize_teff(j["teff_K"].values)
    j["logTe_spec"] = np.log10(j["teff_spec"])

    cache: dict = {}

    def points_at(family, age):
        a = iso[family]
        uages = np.unique(a["age_Myr"].values)
        anear = float(uages[np.argmin(np.abs(uages - age))])
        key = (family, anear)
        if key not in cache:
            cache[key] = iso_hrd_points(a[np.isclose(a["age_Myr"], anear)], family)
        return cache[key], anear

    recs = []
    for _, r in j.iterrows():
        sub = r["subgroup"]
        rec = dict(source_id=int(r["sid"]), subgroup=sub,
                   spectral_type=r.get("spectral_type"),
                   teff_spec=float(r["teff_spec"]),
                   logTe_spec=float(r["logTe_spec"]),
                   membership_probability=float(r["membership_probability"]),
                   av_method=r["av_method"],
                   extreme_hot=bool(r["logTe_spec"] > LOGTE_MAX_TRUST))
        for rv in R_V_BRANCHES:
            mag = r[mag_cols[rv]]
            rec[f"MG0_obs_rv{rv:.1f}"] = float(mag) if pd.notna(mag) else np.nan
        for family in FAMILIES:
            for rv in R_V_BRANCHES:
                tag = branch(family, rv)
                age = age_of.get((sub, family, float(rv)), ens_age[(family, float(rv))])
                (pts, anear) = points_at(family, age)
                mag = rec[f"MG0_obs_rv{rv:.1f}"]
                rec[f"age_used_{tag}"] = float(anear)
                if pts is None or rec["extreme_hot"] or not np.isfinite(mag):
                    rec[f"mass_{tag}"] = np.nan
                    rec[f"G0_model_{tag}"] = np.nan
                    rec[f"chi_{tag}"] = np.nan
                    continue
                mass, gm, chi = match_hrd(rec["logTe_spec"], mag, pts)
                rec[f"mass_{tag}"] = mass
                rec[f"G0_model_{tag}"] = gm
                rec[f"chi_{tag}"] = chi
        # legacy single-branch columns = the R_V = 3.1 baseline branch
        rec["MG0_obs"] = rec[f"MG0_obs_rv{BASELINE_RV:.1f}"]
        for family in FAMILIES:
            tag = branch(family, BASELINE_RV)
            for stem in ("age_used", "mass", "G0_model", "chi"):
                rec[f"{stem}_{family}"] = rec[f"{stem}_{tag}"]
        recs.append(rec)

    out = pd.DataFrame(recs)

    # gate: an anchor is HRD-consistent if the nearest isochrone point (in the
    # (logTe, M_G0) plane, chi metric) is within CHI_TOL for EITHER family.
    # chi=1 is a 1-sigma match given SIG_LOGTE/SIG_MG0; allow up to 2.5 to
    # absorb unresolved-binary brightening and calibration slop.  Scored on the
    # R_V = 3.1 baseline branch as before; every branch is reported.
    CHI_TOL = 2.5
    both = out[(out["extreme_hot"] == False) & out["MG0_obs"].notna()].copy()
    cons = (both["chi_PARSEC"] <= CHI_TOL) | (both["chi_MIST"] <= CHI_TOL)
    out["hrd_consistent"] = False
    out.loc[both.index, "hrd_consistent"] = cons.values

    n_ext = int(out["extreme_hot"].sum())
    n_test = int(len(both))
    n_cons = int(cons.sum())
    frac = n_cons / n_test if n_test else float("nan")
    gate_by_rv = {}
    for rv in R_V_BRANCHES:
        tested = out[(out["extreme_hot"] == False) & out[f"MG0_obs_rv{rv:.1f}"].notna()]
        ok = ((tested[f"chi_{branch('PARSEC', rv)}"] <= CHI_TOL)
              | (tested[f"chi_{branch('MIST', rv)}"] <= CHI_TOL))
        gate_by_rv[f"rv{rv:.1f}"] = {"n_tested": int(len(tested)),
                                     "n_consistent": int(ok.sum())}
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
    ages_used = {
        f"{sub}|{branch(fam, rv)}": round(age_of.get((sub, fam, float(rv)),
                                                     ens_age[(fam, float(rv))]), 4)
        for sub in [*w.SUBGROUPS, "unassigned"] for fam in FAMILIES for rv in R_V_BRANCHES
    }

    out.to_parquet(out_path, index=False)

    print(f"anchors with spectroscopic Teff placed on HRD: {len(out)}")
    print(f"  extreme-hot (WR/stripped, excluded from mass gate): {n_ext}")
    print(f"  HRD-consistent at R_V=3.1 (chi<= {CHI_TOL}, either family): "
          f"{n_cons}/{n_test} = {frac:.1%}")
    for fam in FAMILIES:
        for rv in R_V_BRANCHES:
            tag = branch(fam, rv)
            print(f"  {tag}: median chi {np.nanmedian(out[f'chi_{tag}']):.2f}, mass "
                  f"{np.nanmin(out[f'mass_{tag}']):.1f}-{np.nanmax(out[f'mass_{tag}']):.1f} Msun")

    exec_log = {
        "script": "scripts/wp4_anchors_hrd.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output_version": output_version,
        "age_version": age_version,
        "extinction_version": extinction_version,
        "issue": "#19 -- replaces the anchor masses read at pre-repair ages",
        "preregistration": preregistration,
        "inputs": {
            "wp1_spectroscopic_anchors": sha256(w.PROC / "wp1_spectroscopic_anchors.parquet"),
            str(ext_path.relative_to(w.ROOT)): sha256(ext_path),
            "wp2_subgroup_labels": sha256(w.TABLES / "wp2_subgroup_labels.parquet"),
            str(age_path.relative_to(w.ROOT)): sha256(age_path),
            "isochrones_parsec": sha256(w.PROC / "wp3_isochrones_parsec.parquet"),
            "isochrones_mist": sha256(w.PROC / "wp3_isochrones_mist.parquet"),
        },
        "method": ("place anchors at (logTe_spec, M_G0_obs[R_V]); mass = present-day "
                   "Mass of the nearest isochrone point in the (logTe, M_G0) plane "
                   "under a chi metric (SIG_LOGTE=0.03 dex, SIG_MG0=0.40 mag) on the "
                   "isochrone at the subgroup's upper-MS MAP for that family AND R_V "
                   "branch, snapped to the nearest native grid age; gate = chi<=2.5 "
                   "for either family on the R_V=3.1 branch."),
        "mass_quantity": ("present-day isochrone Mass, unchanged from the frozen "
                          "procedure; Mini (F3(a)) is deferred to a separately "
                          "pre-registered change"),
        "branch_columns": [f"mass_{branch(f, rv)}" for f in FAMILIES for rv in R_V_BRANCHES],
        "legacy_columns": "mass_<family>, age_used_<family>, chi_<family>, MG0_obs = R_V 3.1 branch",
        "ages_used_Myr_before_snapping": ages_used,
        "logTe_max_trust": LOGTE_MAX_TRUST,
        "sig_logte": SIG_LOGTE, "sig_mg0": SIG_MG0, "chi_tol": CHI_TOL,
        "n_anchors_with_teff": len(out),
        "n_extreme_hot_excluded": n_ext,
        "n_gate_tested": n_test,
        "n_gate_consistent": n_cons,
        "frac_consistent": frac,
        "gate_by_rv": gate_by_rv,
        "gate_by_subgroup": gate_by_subgroup,
        "output": {str(out_path.relative_to(w.ROOT)): sha256(out_path)},
    }
    log_path = w.ROOT / "provenance" / f"wp4_anchors_hrd_execution_{output_version}.json"
    with open(log_path, "w") as f:
        json.dump(exec_log, f, indent=2)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--age-version", required=True,
                        help="WP4 age-posterior version, e.g. repair_v5 (required)")
    parser.add_argument("--extinction-version", required=True,
                        help="WP3 extinction version, e.g. repair_v5 (required)")
    parser.add_argument("--output-version", required=True,
                        help="version suffix of the output, e.g. repair_v8")
    parser.add_argument("--preregistration",
                        default="provenance/issue19_repair_v8_prereg.json",
                        help="pre-registration this run executes (provenance only)")
    args = parser.parse_args()
    run(args.age_version, args.extinction_version, args.output_version,
        args.preregistration)


if __name__ == "__main__":
    main()
