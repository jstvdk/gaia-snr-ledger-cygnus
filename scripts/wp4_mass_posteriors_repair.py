#!/usr/bin/env python3
"""Replace frozen nearest-grid masses with versioned posterior summaries."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from wp3_extinction_law import R_V_BRANCHES
from wp3_repair_common import PROC, REPAIR_VERSION, ROOT
from wp4_repair_common import (
    F_BINARY,
    MASS_GRID_SIZE,
    MASS_POSTERIOR_DRAWS,
    N_AGE_NODES,
    N_Q_COMPONENTS,
    Q_MIN,
    RANDOM_SEED,
    MODEL_BAND_SCALE,
    age_posterior_nodes,
    infer_mass_samples,
)


FAMILIES = ["PARSEC", "MIST"]


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    extinction = pd.read_parquet(PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet")
    labels = pd.read_parquet(ROOT / "tables" / "wp2_subgroup_labels.parquet")
    extinction = extinction.drop(columns=["subgroup"], errors="ignore").merge(
        labels[["source_id", "subgroup"]],
        on="source_id",
        how="left",
        validate="one_to_one",
    )
    age_posterior = pd.read_parquet(
        PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet"
    )
    anchor_hrd = pd.read_parquet(PROC / "wp4_anchor_hrd.parquet")
    anchor_lookup = anchor_hrd.set_index("source_id")
    av_store = np.load(PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    posterior_source_id = av_store["source_id"].astype("int64")
    if not np.array_equal(
        posterior_source_id, extinction["source_id"].to_numpy("int64")
    ):
        raise RuntimeError("WP3 posterior cube and repaired extinction row order differ")
    av_probability = av_store["probability"]

    branches = [
        (family, float(rv))
        for family in FAMILIES
        for rv in R_V_BRANCHES
    ]
    mass_cube = np.full(
        (len(extinction), len(branches), MASS_POSTERIOR_DRAWS),
        np.nan,
        dtype=np.float32,
    )
    rows: list[dict] = []
    for star_index, row in extinction.iterrows():
        source_id = int(row["source_id"])
        anchor = anchor_lookup.loc[source_id] if source_id in anchor_lookup.index else None
        spectroscopic_ok = (
            anchor is not None
            and not bool(anchor["extreme_hot"])
            and np.isfinite(anchor["mass_PARSEC"])
            and np.isfinite(anchor["mass_MIST"])
        )
        record: dict[str, object] = {
            "source_id": source_id,
            "subgroup": row["subgroup"],
            "membership_probability": float(row["membership_probability"]),
            "av_method": row["av_method"],
            "mass_method": (
                "spectroscopic_hrd_frozen"
                if spectroscopic_ok
                else "photometric_posterior_repair_v1"
            ),
            "is_spectroscopic_anchor": bool(anchor is not None),
        }
        for branch_index, (family, rv) in enumerate(branches):
            column = f"mass_{family}_rv{rv:.1f}"
            if spectroscopic_ok:
                samples = np.full(
                    MASS_POSTERIOR_DRAWS, float(anchor[f"mass_{family}"])
                )
                age_nodes = np.array([float(anchor[f"age_used_{family}"])])
            else:
                nodes = age_posterior_nodes(
                    age_posterior, str(row["subgroup"]), family, rv
                )
                samples = infer_mass_samples(
                    row,
                    av_probability[star_index, list(R_V_BRANCHES).index(rv)],
                    family,
                    rv,
                    nodes,
                    seed_offset=100_000 * branch_index,
                )
                age_nodes = nodes
            mass_cube[star_index, branch_index] = samples.astype(np.float32)
            finite = samples[np.isfinite(samples)]
            if len(finite):
                q16, q50, q84 = np.quantile(finite, [0.16, 0.50, 0.84])
                record[column] = float(q50)
                record[f"{column}_q16"] = float(q16)
                record[f"{column}_q84"] = float(q84)
                record[f"{column}_sd"] = float(np.std(finite))
                record[f"{column}_p_gt8"] = float(np.mean(finite > 8.0))
            else:
                for suffix in ["", "_q16", "_q84", "_sd", "_p_gt8"]:
                    record[f"{column}{suffix}"] = np.nan
            record[f"age_median_{family}_rv{rv:.1f}"] = float(np.median(age_nodes))
        record["mass_baseline"] = record["mass_PARSEC_rv3.1"]
        record["mass_baseline_q16"] = record["mass_PARSEC_rv3.1_q16"]
        record["mass_baseline_q84"] = record["mass_PARSEC_rv3.1_q84"]
        record["mass_baseline_p_gt8"] = record["mass_PARSEC_rv3.1_p_gt8"]
        record["mass_parsec_mist_spread"] = abs(
            record["mass_PARSEC_rv3.1"] - record["mass_MIST_rv3.1"]
        )
        rows.append(record)

    output = pd.DataFrame(rows)
    output_path = PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet"
    output.to_parquet(output_path, index=False)
    samples_path = PROC / f"wp4_mass_posterior_samples_{REPAIR_VERSION}.npz"
    np.savez_compressed(
        samples_path,
        source_id=output["source_id"].to_numpy("int64"),
        family=np.array([family for family, _ in branches]),
        rv=np.array([rv for _, rv in branches]),
        samples=mass_cube,
    )

    baseline = output["mass_baseline"].dropna()
    exact_fraction = float(baseline.value_counts(normalize=True).max())
    provenance = {
        "repair_version": REPAIR_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp4_mass_posteriors_repair.py",
        "method": (
            "Monte Carlo posterior propagation of the complete WP3 A_V grid, "
            "G/BP/RP/J/H/Ks photometric errors plus floor, unresolved-binary state and q, "
            "and branch-specific subgroup age posterior; continuous parabolic "
            "projection on a dense initial-mass grid; frozen spectroscopic HRD override"
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "inputs": {
            str(path): sha256(ROOT / path)
            for path in [
                f"data/processed/wp3_extinction_{REPAIR_VERSION}.parquet",
                f"data/processed/wp3_extinction_posterior_{REPAIR_VERSION}.npz",
                f"data/processed/wp4_age_posteriors_{REPAIR_VERSION}.parquet",
                "data/processed/wp4_anchor_hrd.parquet",
                "data/processed/wp3_isochrones_parsec.parquet",
                "data/processed/wp3_isochrones_mist.parquet",
                "tables/wp2_subgroup_labels.parquet",
            ]
        },
        "configuration": {
            "posterior_draws_per_star_branch": MASS_POSTERIOR_DRAWS,
            "dense_initial_mass_grid_points": MASS_GRID_SIZE,
            "age_posterior_nodes": N_AGE_NODES,
            "binary_fraction": F_BINARY,
            "minimum_mass_ratio": Q_MIN,
            "binary_q_components": N_Q_COMPONENTS,
            "random_seed": RANDOM_SEED,
            "mass_reporting_statistic": "posterior median",
            "six_band_model_scale_mag": MODEL_BAND_SCALE,
            "six_band_model_scale_calibration": (
                "median robust scatter across G/BP/RP/J/H/Ks for 107 "
                "non-extreme spectroscopic-HRD anchors; per-band range 0.35-0.39 mag"
            ),
            "mass_prior": (
                "log-uniform scale prior on the dense mass grid; no IMF slope "
                "is imposed before WP5"
            ),
        },
        "counts": {
            "members": int(len(output)),
            "finite_baseline_mass": int(output["mass_baseline"].notna().sum()),
            "spectroscopic_hrd_overrides": int(
                output["mass_method"].eq("spectroscopic_hrd_frozen").sum()
            ),
            "photometric_posteriors": int(
                output["mass_method"].eq("photometric_posterior_repair_v1").sum()
            ),
            "baseline_between_2p5_3p2": int(
                output["mass_baseline"].between(2.5, 3.2).sum()
            ),
            "baseline_posterior_median_gt8": int((output["mass_baseline"] > 8).sum()),
            "baseline_expected_gt8_from_posterior": float(
                output["mass_baseline_p_gt8"].sum(skipna=True)
            ),
        },
        "acceptance_diagnostics": {
            "largest_exact_baseline_mass_fraction": exact_fraction,
            "largest_exact_mass_fraction_gate": 0.02,
            "largest_exact_mass_fraction_pass": bool(exact_fraction < 0.02),
        },
        "outputs": {
            str(path.relative_to(ROOT)): {
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in [output_path, samples_path]
        },
        "frozen_outputs_overwritten": False,
    }
    provenance_path = ROOT / "provenance" / "wp4_mass_repair_execution.json"
    temporary = provenance_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    temporary.replace(provenance_path)
    print(f"wrote {output_path.relative_to(ROOT)} ({len(output)} rows)")
    print(
        f"baseline finite={len(baseline)}; exact-value max={100*exact_fraction:.2f}%; "
        f"2.5-3.2 Msun={provenance['counts']['baseline_between_2p5_3p2']}; "
        f"median >8={provenance['counts']['baseline_posterior_median_gt8']}; "
        f"posterior E[N>8]={provenance['counts']['baseline_expected_gt8_from_posterior']:.1f}"
    )


if __name__ == "__main__":
    main()
