#!/usr/bin/env python3
"""Gate G1: WP5 baseline refit restricted to seed-stable subgroup labels.

Step 1 of the gated plan in tasks/wp5_cygob2b_age_caustic_fix_brief.md.
Consumes provenance/wp2_label_stability_per_star.csv (written by
scripts/wp2_label_stability.py) and refits the frozen repair_v3 WP5
normalization keeping only stars whose A/B/C label is stable in >=90% of the
50 WP2 GMM seeds.  Monkeypatch pattern of wp5_lower_edge_scan.py: the stored
repair_v3 artifacts are read, never written.

Gate G1 (brief section 3, step 1): labels are exonerated if the stable-label
refit moves CygOB2-B's baseline bin-2 residual by < 0.5 and the stability
audit shows no mass-dependent instability in B.

Output (new file only): provenance/wp5_stable_label_refit_execution.json

Run:  PYTHONPATH=scripts python3 scripts/wp5_stable_label_refit.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp5_fit_imf as F

VERSION = "repair_v3"
STABLE_FRACTION = 0.90
G1_RESIDUAL_SHIFT_LIMIT = 0.5
BUMP_BIN = 2

masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{VERSION}.parquet")
curves = pd.read_parquet(w.PROC / f"wp5_completeness_curves_{VERSION}.parquet")
responses = pd.read_parquet(w.PROC / f"wp5_injection_response_{VERSION}.parquet")
store = np.load(w.PROC / f"wp4_mass_posterior_samples_{VERSION}.npz")
samples = store["samples"]
stability = pd.read_csv(w.PROVENANCE / "wp2_label_stability_per_star.csv")
stored_bins = pd.read_parquet(w.PROC / f"wp5_mass_function_bins_{VERSION}.parquet")
stored_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{VERSION}.parquet")
label_audit = json.loads(
    (w.PROVENANCE / "wp2_label_stability_execution.json").read_text(encoding="utf-8")
)

stable_ids = set(
    stability.loc[
        stability["consensus_label_stability"].ge(STABLE_FRACTION), "source_id"
    ].astype("int64")
)
labelled = masses["subgroup"].isin(w.SUBGROUPS).to_numpy()
keep = ~labelled | masses["source_id"].astype("int64").isin(stable_ids).to_numpy()
removed = masses.loc[~keep, ["source_id", "subgroup"]]
masses_stable = masses.loc[keep].reset_index(drop=True)
samples_stable = samples[keep]

draw_cols = sorted(c for c in responses.columns if c.startswith("recovered_mass_draw_"))
rng = np.random.default_rng(w.SEED)
rows = []
bin_rows = []
for family in w.FAMILIES:
    for rv in w.R_V_BRANCHES:
        for alpha in w.IMF_SLOPES:
            for sg in w.SUBGROUPS:
                cu = curves[curves.family.eq(family) & curves.R_V.eq(rv) & curves.subgroup.eq(sg)]
                rp = responses[responses.family.eq(family) & responses.R_V.eq(rv) & responses.subgroup.eq(sg)]
                s, b, _ = F.fit_one(
                    masses_stable, cu, rp, sg, family, rv, alpha, rng,
                    samples_stable[:, w.FAMILIES.index(family) * len(w.R_V_BRANCHES)
                                   + w.R_V_BRANCHES.index(rv), :], draw_cols)
                rows.append(s)
                bin_rows.append(b)
refit = pd.DataFrame(rows)
refit_bins = pd.concat(bin_rows, ignore_index=True)


def baseline_residuals(bins: pd.DataFrame, sg: str) -> list[float]:
    sel = bins[
        bins.subgroup.eq(sg) & bins.family.eq("PARSEC") & bins.R_V.eq(3.1)
        & bins.alpha.eq(2.3)
    ].sort_values("bin_index")
    return [float(r) for r in sel["pearson_residual"]]


comparison = {}
for sg in w.SUBGROUPS:
    stored_res = baseline_residuals(stored_bins, sg)
    refit_res = baseline_residuals(refit_bins, sg)
    stored_row = stored_norm[
        stored_norm.subgroup.eq(sg) & stored_norm.family.eq("PARSEC")
        & stored_norm.R_V.eq(3.1) & stored_norm.alpha.eq(2.3)
    ].iloc[0]
    refit_row = refit[
        refit.subgroup.eq(sg) & refit.family.eq("PARSEC")
        & refit.R_V.eq(3.1) & refit.alpha.eq(2.3)
    ].iloc[0]
    comparison[sg] = {
        "stored_residuals": stored_res,
        "stable_label_residuals": refit_res,
        "residual_shift": [float(a - b) for a, b in zip(refit_res, stored_res)],
        "stored_gate": {
            "chi2_p": float(stored_row.poisson_chi_square_p),
            "trend_p": float(stored_row.residual_trend_p),
            "max_abs_residual": float(stored_row.max_abs_pearson_residual),
            "pass": bool(stored_row.residual_gate_pass),
        },
        "stable_label_gate": {
            "chi2_p": float(refit_row.poisson_chi_square_p),
            "trend_p": float(refit_row.residual_trend_p),
            "max_abs_residual": float(refit_row.max_abs_pearson_residual),
            "pass": bool(refit_row.residual_gate_pass),
        },
        "n_sources_removed": int((removed.subgroup.eq(sg)).sum()),
    }

b_shift = abs(comparison["CygOB2-B"]["residual_shift"][BUMP_BIN])
b_profile = label_audit["mass_profile_baseline_branch"]["CygOB2-B"]
trend = b_profile["stability_vs_mass_trend"]
no_mass_dependence = bool(
    (trend["p"] is None or trend["p"] >= 0.05)
    and all(
        cell["fraction_stable_ge_090"] is None
        or cell["fraction_stable_ge_090"] >= 0.95
        for cell in b_profile["window_bins"]
    )
)
gate_g1 = {
    "bump_bin_residual_shift_abs": float(b_shift),
    "shift_limit": G1_RESIDUAL_SHIFT_LIMIT,
    "shift_within_limit": bool(b_shift < G1_RESIDUAL_SHIFT_LIMIT),
    "no_mass_dependent_instability_in_B": no_mass_dependence,
    "labels_exonerated": bool(
        b_shift < G1_RESIDUAL_SHIFT_LIMIT and no_mass_dependence
    ),
    "verdict": (
        "PASS — labels exonerated, proceed to step 2 (age-conditional refit)"
        if b_shift < G1_RESIDUAL_SHIFT_LIMIT and no_mass_dependence
        else "FAIL — contamination becomes the primary path (step 4a)"
    ),
}

grid_pass_stored = int(stored_norm.residual_gate_pass.sum())
grid_pass_refit = int(refit.residual_gate_pass.sum())

record = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "script": "scripts/wp5_stable_label_refit.py",
    "status": "SUCCESS",
    "issue": "#1c step 1 / #6 closure — stable-label WP5 refit (gate G1)",
    "wp5_version_consumed": VERSION,
    "stored_artifacts_overwritten": False,
    "seed": w.SEED,
    "environment": {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
    "inputs": {
        str(p.relative_to(w.ROOT)): w.sha256(p) for p in [
            w.PROC / f"wp4_mass_posteriors_{VERSION}.parquet",
            w.PROC / f"wp5_completeness_curves_{VERSION}.parquet",
            w.PROC / f"wp5_injection_response_{VERSION}.parquet",
            w.PROC / f"wp4_mass_posterior_samples_{VERSION}.npz",
            w.PROC / f"wp5_mass_function_bins_{VERSION}.parquet",
            w.PROC / f"wp5_imf_normalization_{VERSION}.parquet",
            w.PROVENANCE / "wp2_label_stability_per_star.csv",
            w.PROVENANCE / "wp2_label_stability_execution.json",
        ]
    },
    "selection": {
        "stable_threshold": STABLE_FRACTION,
        "n_labelled_sources": int(labelled.sum()),
        "n_removed": int(len(removed)),
        "removed": removed.assign(
            source_id=removed.source_id.astype("int64")
        ).to_dict(orient="records"),
    },
    "baseline_comparison": comparison,
    "grid_54": {
        "stored_branches_passing": grid_pass_stored,
        "stable_label_branches_passing": grid_pass_refit,
    },
    "gate_G1": gate_g1,
}
w.write_json(w.PROVENANCE / "wp5_stable_label_refit_execution.json", record)
print(json.dumps({"baseline_comparison": comparison, "gate_G1": gate_g1,
                  "grid_54": record["grid_54"]}, indent=2))
print("wrote provenance/wp5_stable_label_refit_execution.json")
