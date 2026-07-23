#!/usr/bin/env python3
"""WP3 step: assemble the per-member photometry table.

Full WP2 member list (P>0.5, 1392 stars) joined to Gaia DR3 G/BP/RP (+ errors,
BP/RP excess, RUWE) and to the 2MASS PSC cross-match (J/H/Ks + errors + quality
flags). Flags the members without a 2MASS match. No stars are dropped.

Output: data/processed/wp3_member_photometry.parquet
Provenance: provenance/wp3_assemble_photometry_execution.json
"""
import hashlib, json
import numpy as np, pandas as pd
from wp3_common import gaia_mag_error

def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

def main():
    members = pd.read_parquet("data/processed/wp2_members.parquet")
    members = members[members.membership_probability > 0.5].copy()
    n0 = len(members)

    narrow = pd.read_parquet("data/processed/wp1_gaia_narrow.parquet", columns=[
        "source_id", "phot_g_mean_mag", "phot_g_mean_flux_error",
        "phot_bp_mean_mag", "phot_bp_mean_flux_error",
        "phot_rp_mean_mag", "phot_rp_mean_flux_error",
        "phot_bp_rp_excess_factor", "ruwe", "parallax_error"])
    tmass = pd.read_parquet("data/processed/wp1_2mass_join.parquet", columns=[
        "source_id", "j_m", "j_msigcom", "h_m", "h_msigcom", "ks_m", "ks_msigcom",
        "ph_qual", "ph_qual_aaa", "has_2mass_psc_match", "has_j", "has_h", "has_ks",
        "has_complete_jhk", "tmass_designation", "tmass_angular_distance"])

    df = members.merge(narrow, on="source_id", how="left", suffixes=("", "_narrow"))
    df = df.merge(tmass, on="source_id", how="left")

    # Gaia magnitudes + errors
    df["G"] = df["phot_g_mean_mag"]
    df["BP"] = df["phot_bp_mean_mag"]
    df["RP"] = df["phot_rp_mean_mag"]
    df["G_err"] = gaia_mag_error(df["G"], df["phot_g_mean_flux_error"], "G")
    df["BP_err"] = gaia_mag_error(df["BP"], df["phot_bp_mean_flux_error"], "BP")
    df["RP_err"] = gaia_mag_error(df["RP"], df["phot_rp_mean_flux_error"], "RP")

    # 2MASS magnitudes + errors
    df["J"] = df["j_m"]; df["H"] = df["h_m"]; df["Ks"] = df["ks_m"]
    df["J_err"] = df["j_msigcom"]; df["H_err"] = df["h_msigcom"]; df["Ks_err"] = df["ks_msigcom"]

    # Flags
    df["has_2mass"] = df["has_2mass_psc_match"].fillna(False).astype(bool)
    df["has_complete_jhk"] = df["has_complete_jhk"].fillna(False).astype(bool)
    df["ph_qual_aaa"] = df["ph_qual_aaa"].fillna(False).astype(bool)
    df["in_narrow_catalogue"] = df["phot_g_mean_mag"].notna()

    # per-band validity
    for b in ["G", "BP", "RP"]:
        df[f"has_{b}"] = df[b].notna()
    df["has_J"] = df["has_j"].fillna(False).astype(bool) & df["J"].notna()
    df["has_H"] = df["has_h"].fillna(False).astype(bool) & df["H"].notna()
    df["has_Ks"] = df["has_ks"].fillna(False).astype(bool) & df["Ks"].notna()
    df["n_bands"] = df[[f"has_{b}" for b in ["G", "BP", "RP", "J", "H", "Ks"]]].sum(axis=1)

    keep = ["source_id", "ra", "dec", "l_deg", "b_deg",
            "parallax_corrected", "parallax_error", "ruwe",
            "membership_probability", "anchor_quality_exempt", "subgroup_label",
            "G", "G_err", "BP", "BP_err", "RP", "RP_err",
            "J", "J_err", "H", "H_err", "Ks", "Ks_err",
            "phot_bp_rp_excess_factor", "ph_qual", "ph_qual_aaa",
            "has_G", "has_BP", "has_RP", "has_J", "has_H", "has_Ks",
            "has_2mass", "has_complete_jhk", "in_narrow_catalogue", "n_bands",
            "tmass_angular_distance"]
    keep = [c for c in keep if c in df.columns]
    out = df[keep].copy()

    n_no2mass = int((~out.has_2mass).sum())
    n_complete = int(out.has_complete_jhk.sum())
    n_aaa = int(out.ph_qual_aaa.sum())
    n_no_gaia = int((~out.in_narrow_catalogue).sum())

    p = "data/processed/wp3_member_photometry.parquet"
    out.to_parquet(p, index=False)

    prov = {
        "script": "scripts/wp3_assemble_photometry.py",
        "inputs": {
            "data/processed/wp2_members.parquet": sha("data/processed/wp2_members.parquet"),
            "data/processed/wp1_gaia_narrow.parquet": sha("data/processed/wp1_gaia_narrow.parquet"),
            "data/processed/wp1_2mass_join.parquet": sha("data/processed/wp1_2mass_join.parquet"),
        },
        "member_selection": "membership_probability > 0.5",
        "n_members": n0,
        "n_with_2mass_psc": int(out.has_2mass.sum()),
        "n_without_2mass_flagged": n_no2mass,
        "n_complete_jhk": n_complete,
        "n_ph_qual_aaa": n_aaa,
        "n_without_gaia_photometry": n_no_gaia,
        "gaia_mag_error_model": "sigma_mag=1.0857*flux_err/flux via DR3 VEGAMAG ZP; 0.003 mag floor",
        "output": {p: sha(p), "rows": int(len(out))},
    }
    json.dump(prov, open("provenance/wp3_assemble_photometry_execution.json", "w"), indent=2)
    print(f"members P>0.5: {n0}")
    print(f"  with 2MASS PSC match : {int(out.has_2mass.sum())}")
    print(f"  without 2MASS (flag) : {n_no2mass}")
    print(f"  complete JHK         : {n_complete}")
    print(f"  ph_qual AAA          : {n_aaa}")
    print(f"  no Gaia phot (anchors): {n_no_gaia}")
    print(f"WROTE {p}  {out.shape}")

if __name__ == "__main__":
    main()
