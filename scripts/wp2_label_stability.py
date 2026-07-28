#!/usr/bin/env python3
"""Per-star A/B/C label-stability audit over the 50 frozen WP2 GMM seeds.

Step 1 of the gated plan in tasks/wp5_cygob2b_age_caustic_fix_brief.md
(issue #6: subgroup-label uncertainty was never quantified; issue #1c needs
labels exonerated or indicted before the age-conditional refit).

Read-only with respect to every stored pipeline artifact.  The frozen WP2
procedure is replayed exactly: same clean-member selection, same
StandardScaler feature space (l, b, pmra, pmdec), same GaussianMixture
hyperparameters, same 50 deterministic seeds, adopted k = 3.  Each seed's
components are named CygOB2-A/B/C with the same deterministic physical rule
used for the consensus (wp2_derive_subgroups.name_components), so per-seed
labels are comparable star by star with the frozen sidecar labels.

Outputs (new files only):
  provenance/wp2_label_stability_per_star.csv
  provenance/wp2_label_stability_execution.json

Run:  PYTHONPATH=scripts python3 scripts/wp2_label_stability.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import sklearn
from sklearn.preprocessing import StandardScaler

from wp2_derive_subgroups import (
    FEATURES,
    SEEDS,
    fit_labels,
    load_clean,
    name_components,
)
import wp5_common as w

ADOPTED_K = 3
STABLE_FRACTION = 0.90
BASELINE_MASS_COLUMN = "mass_PARSEC_rv3.1"
WINDOW_BIN_EDGES = np.geomspace(2.0, 8.0, 7)
REPAIR_VERSION = "repair_v3"
LABELS = ["CygOB2-A", "CygOB2-B", "CygOB2-C"]

INPUTS = [
    w.PROC / "wp2_members.parquet",
    w.TABLES / "wp2_subgroup_labels.parquet",
    w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
    w.PROVENANCE / "wp2_subgroups_execution.json",
]


def main() -> None:
    _, clean = load_clean()
    sidecar = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    consensus = clean[["source_id"]].merge(
        sidecar[["source_id", "subgroup"]], on="source_id", how="left", validate="one_to_one"
    )
    if consensus["subgroup"].isna().any() or len(consensus) != len(clean):
        raise RuntimeError("clean members do not align with the frozen sidecar labels")
    consensus_label = consensus["subgroup"].to_numpy()

    scaled = StandardScaler().fit_transform(clean[FEATURES].values)
    per_seed = np.empty((len(SEEDS), len(clean)), dtype=object)
    for seed_index, seed in enumerate(SEEDS):
        raw = fit_labels(scaled, ADOPTED_K, seed).predict(scaled)
        naming = name_components(clean, raw)
        per_seed[seed_index] = np.array([naming[c] for c in raw])

    fraction = {
        label: (per_seed == label).mean(axis=0) for label in LABELS
    }
    stability = np.zeros(len(clean))
    for label in LABELS:
        stability = np.where(consensus_label == label, fraction[label], stability)

    per_star = pd.DataFrame(
        {
            "source_id": clean["source_id"].to_numpy("int64"),
            "subgroup_consensus": consensus_label,
            **{f"seed_fraction_{label[-1]}": fraction[label] for label in LABELS},
            "consensus_label_stability": stability,
            "stable_ge_090": stability >= STABLE_FRACTION,
        }
    )
    csv_path = w.PROVENANCE / "wp2_label_stability_per_star.csv"
    per_star.to_csv(csv_path, index=False)

    # Seed-averaged confusion matrix: rows = consensus label, columns = the
    # label the same star receives in an individual seed run.
    confusion = {
        row: {
            col: float(
                (per_seed[:, consensus_label == row] == col).mean()
            )
            for col in LABELS
        }
        for row in LABELS
    }

    summary_by_subgroup = {
        label: {
            "n_stars": int((consensus_label == label).sum()),
            "stability_mean": float(stability[consensus_label == label].mean()),
            "stability_min": float(stability[consensus_label == label].min()),
            "n_below_090": int(
                (stability[consensus_label == label] < STABLE_FRACTION).sum()
            ),
        }
        for label in LABELS
    }

    # Mass-dependent instability check for CygOB2-B (the smoking-gun test the
    # superseded brief prescribed): stability profile across the baseline
    # observed-mass window bins.
    masses = pd.read_parquet(
        w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet"
    )[["source_id", "subgroup", BASELINE_MASS_COLUMN]]
    joined = per_star.merge(masses, on="source_id", how="inner", validate="one_to_one")
    mass_profile = {}
    for label in LABELS:
        sub = joined[joined["subgroup_consensus"].eq(label)]
        bins = []
        for lo, hi in zip(WINDOW_BIN_EDGES[:-1], WINDOW_BIN_EDGES[1:]):
            cell = sub[
                sub[BASELINE_MASS_COLUMN].ge(lo) & sub[BASELINE_MASS_COLUMN].lt(hi)
            ]
            bins.append(
                {
                    "mass_lo": float(lo),
                    "mass_hi": float(hi),
                    "n": int(len(cell)),
                    "stability_median": (
                        float(cell["consensus_label_stability"].median())
                        if len(cell)
                        else None
                    ),
                    "fraction_stable_ge_090": (
                        float(cell["stable_ge_090"].mean()) if len(cell) else None
                    ),
                }
            )
        window = sub[
            sub[BASELINE_MASS_COLUMN].ge(WINDOW_BIN_EDGES[0])
            & sub[BASELINE_MASS_COLUMN].lt(WINDOW_BIN_EDGES[-1])
        ]
        from scipy import stats

        if len(window) > 2 and window["consensus_label_stability"].nunique() > 1:
            rho = stats.spearmanr(
                window[BASELINE_MASS_COLUMN], window["consensus_label_stability"]
            )
            trend = {"spearman_rho": float(rho.statistic), "p": float(rho.pvalue)}
        else:
            trend = {"spearman_rho": None, "p": None, "note": "stability constant in window"}
        mass_profile[label] = {
            "window_bins": bins,
            "window_n": int(len(window)),
            "window_fraction_stable_ge_090": (
                float(window["stable_ge_090"].mean()) if len(window) else None
            ),
            "stability_vs_mass_trend": trend,
        }

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp2_label_stability.py",
        "status": "SUCCESS",
        "issue": "#6 subgroup-label uncertainty / step 1 of the #1c gated plan",
        "deterministic": True,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "sklearn": sklearn.__version__,
        },
        "method": {
            "replay": (
                "frozen WP2 clean-member selection (P>0.5, not anchor-exempt), "
                "StandardScaler on (l_deg, b_deg, pmra, pmdec), GaussianMixture "
                "k=3 full covariance reg_covar=1e-4 max_iter=500 n_init=1, the "
                "50 frozen seeds 1000..1049, per-seed components named with "
                "wp2_derive_subgroups.name_components"
            ),
            "stability_definition": (
                "fraction of the 50 seed runs assigning the star its frozen "
                "consensus (sidecar) label"
            ),
            "stable_threshold": STABLE_FRACTION,
            "adopted_k": ADOPTED_K,
        },
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in INPUTS},
        "n_clean_members": int(len(clean)),
        "seed_averaged_confusion_matrix_rows_consensus": confusion,
        "stability_by_subgroup": summary_by_subgroup,
        "mass_profile_baseline_branch": mass_profile,
        "outputs": {
            str(csv_path.relative_to(w.ROOT)): w.sha256(csv_path),
        },
    }
    w.write_json(w.PROVENANCE / "wp2_label_stability_execution.json", record)
    print(json.dumps({k: record[k] for k in [
        "seed_averaged_confusion_matrix_rows_consensus",
        "stability_by_subgroup",
    ]}, indent=2))
    print(json.dumps(mass_profile["CygOB2-B"], indent=2))
    print("wrote", csv_path.relative_to(w.ROOT))
    print("wrote provenance/wp2_label_stability_execution.json")


if __name__ == "__main__":
    main()
