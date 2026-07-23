#!/usr/bin/env python3
"""WP4 reporting & provenance.

Emits:
  tables/wp4_ages_summary.csv        - machine-readable per-branch age table
  tables/wp4_ages_table.md           - the full per-branch markdown table (embedded in wp4_ages.md)
  tables/wp4_masses.cat              - whitespace catalogue of per-star masses (WP3 .cat pattern)
  provenance/wp4_manifest.json       - SHA-256 per input/output, WP1-WP3 pattern
  provenance/wp4_provenance.md       - the reproducibility narrative
"""
from __future__ import annotations

import json, hashlib, datetime as dt
import numpy as np
import pandas as pd

import wp4_common as w


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_tables():
    post = pd.read_parquet(w.PROC / "wp4_age_posteriors.parquet")
    # machine-readable summary (all branch rows)
    cols = ["subgroup", "family", "R_V", "f_bin", "indicator", "dmu", "n_stars",
            "measurable", "age_map", "age_lo68", "age_hi68", "age_lo90", "age_hi90"]
    post[cols].to_csv(w.TABLES / "wp4_ages_summary.csv", index=False)

    # compact markdown table: baseline (R_V=3.1, f_bin=0.4, dmu=0) + envelope
    base = post[(post.R_V == 3.1) & (post.f_bin == 0.4) & (post.dmu == 0.0)]
    lines = ["| Subgroup | Indicator | Family | Age MAP (Myr) | 68% CI | n | Measurable |",
             "|---|---|---|---|---|---|---|"]
    for sub in w.SUBGROUPS:
        for ind in ["ums", "pms"]:
            for fam in ["PARSEC", "MIST"]:
                r = base[(base.subgroup == sub) & (base.indicator == ind)
                         & (base.family == fam)]
                if not len(r):
                    continue
                r = r.iloc[0]
                meas = "yes" if r.measurable else "**no (starved)**"
                lines.append(
                    f"| {sub} | {ind} | {fam} | {r.age_map:.2f} | "
                    f"[{r.age_lo68:.2f}, {r.age_hi68:.2f}] | {int(r.n_stars)} | {meas} |")
    (w.TABLES / "wp4_ages_table.md").write_text("\n".join(lines) + "\n")

    # full-envelope table across R_V x f_bin x distance per subgroup/indicator/family
    env_lines = ["| Subgroup | Indicator | Family | MAP range (Myr) | 68% CI union |",
                 "|---|---|---|---|---|"]
    meas = post[post.measurable]
    for sub in w.SUBGROUPS:
        for ind in ["ums", "pms"]:
            for fam in ["PARSEC", "MIST"]:
                g = meas[(meas.subgroup == sub) & (meas.indicator == ind)
                         & (meas.family == fam)]
                if not len(g):
                    env_lines.append(f"| {sub} | {ind} | {fam} | - | not measurable |")
                    continue
                env_lines.append(
                    f"| {sub} | {ind} | {fam} | [{g.age_map.min():.2f}, {g.age_map.max():.2f}] "
                    f"| [{g.age_lo68.min():.2f}, {g.age_hi68.max():.2f}] |")
    (w.TABLES / "wp4_ages_envelope.md").write_text("\n".join(env_lines) + "\n")
    return post


def build_masses_cat():
    df = pd.read_parquet(w.PROC / "wp4_masses.parquet")
    cols = ["source_id", "subgroup", "membership_probability", "mass_method",
            "is_spectroscopic_anchor", "mass_PARSEC_rv3.1", "mass_MIST_rv3.1",
            "mass_PARSEC_rv3.0", "mass_PARSEC_rv3.5", "mass_MIST_rv3.0",
            "mass_MIST_rv3.5", "mass_baseline", "mass_parsec_mist_spread"]
    df[cols].to_csv(w.TABLES / "wp4_masses.cat", sep=" ", index=False,
                    float_format="%.4f", na_rep="nan")


def build_manifest():
    created = dt.datetime.now(dt.timezone.utc).isoformat()
    inputs = {
        "wp3_extinction": w.PROC / "wp3_extinction.parquet",
        "wp3_isochrones_parsec": w.PROC / "wp3_isochrones_parsec.parquet",
        "wp3_isochrones_mist": w.PROC / "wp3_isochrones_mist.parquet",
        "wp2_subgroup_labels": w.TABLES / "wp2_subgroup_labels.parquet",
        "wp2_anchor_assignments": w.PROC / "wp2_anchor_assignments.parquet",
        "wp1_spectroscopic_anchors": w.PROC / "wp1_spectroscopic_anchors.parquet",
    }
    scripts = ["wp4_common.py", "wp4_fit_ages.py", "wp4_anchors_hrd.py",
               "wp4_clump.py", "wp4_masses.py", "wp4_figures.py", "wp4_report.py"]
    exec_logs = ["wp4_fit_ages_execution.json", "wp4_anchors_hrd_execution.json",
                 "wp4_clump_execution.json", "wp4_masses_execution.json"]
    outputs = {
        "data/processed/wp4_age_posteriors.parquet": w.PROC / "wp4_age_posteriors.parquet",
        "data/processed/wp4_posterior_curves.npz": w.PROC / "wp4_posterior_curves.npz",
        "data/processed/wp4_anchor_hrd.parquet": w.PROC / "wp4_anchor_hrd.parquet",
        "data/processed/wp4_clump.parquet": w.PROC / "wp4_clump.parquet",
        "data/processed/wp4_masses.parquet": w.PROC / "wp4_masses.parquet",
        "tables/wp4_ages_summary.csv": w.TABLES / "wp4_ages_summary.csv",
        "tables/wp4_masses.cat": w.TABLES / "wp4_masses.cat",
        "figures/wp4/wp4_cmd_subgroups.png": w.ROOT / "figures/wp4/wp4_cmd_subgroups.png",
        "figures/wp4/wp4_hrd_anchors.png": w.ROOT / "figures/wp4/wp4_hrd_anchors.png",
        "figures/wp4/wp4_age_posteriors.png": w.ROOT / "figures/wp4/wp4_age_posteriors.png",
        "figures/wp4/wp4_age_summary.png": w.ROOT / "figures/wp4/wp4_age_summary.png",
    }
    manifest = {
        "work_package": "WP4 - subgroup ages & per-star masses",
        "created_utc": created,
        "status": "WP4_COMPLETE_GATE_PASSED",
        "distance_posterior_kpc": [w.D_KPC, w.D_KPC_ERR],
        "dist_modulus": w.DIST_MODULUS, "dist_modulus_err": w.DIST_MODULUS_ERR,
        "branches_carried": {
            "isochrone_family": ["PARSEC", "MIST"],
            "R_V": w.R_V_BRANCHES,
            "f_bin": w.F_BIN_BRANCHES,
            "age_indicator": ["upper_main_sequence", "pms_turn_on"],
            "distance_shift_mag": [0.0, round(w.DIST_MODULUS_ERR, 4), round(-w.DIST_MODULUS_ERR, 4)],
            "star_formation_duration_Myr_carried_to_WP7": [0, 1, 2],
        },
        "upstream_inputs": {k: sha256(v) for k, v in inputs.items()},
        "scripts": {f"scripts/{s}": sha256(w.ROOT / "scripts" / s) for s in scripts},
        "execution_logs": {f"provenance/{e}": sha256(w.ROOT / "provenance" / e)
                           for e in exec_logs},
        "outputs": {k: {"sha256": sha256(v), "bytes": v.stat().st_size}
                    for k, v in outputs.items()},
        "shared_isochrones_with_WP3": True,
        "reports": ["wp4_ages.md", "wp4_completion_report.md",
                    "provenance/wp4_provenance.md"],
    }
    with open(w.ROOT / "provenance" / "wp4_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


if __name__ == "__main__":
    build_tables()
    build_masses_cat()
    mani = build_manifest()
    print("wrote tables, masses.cat, manifest")
    print("outputs:", len(mani["outputs"]), "inputs:", len(mani["upstream_inputs"]))
