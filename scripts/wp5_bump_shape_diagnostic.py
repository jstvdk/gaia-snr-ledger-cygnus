"""Shape diagnostic for the CygOB2-B mass-function residual (issue #1c).

Three read-only measurements on the stored repair_v3 artifacts; nothing is
refit and no pipeline artifact is touched.

1. **Alpha shape scan** — rebuild the forward Poisson rates from the stored
   injection response for CygOB2-B at every family/R_V and scan the IMF slope
   alpha continuously.  Question: is there ANY single power law that flattens
   B's residuals?  (If the bin-2 excess survives all alpha, the anomaly is a
   localized bump, not a slope, and "B has a steeper IMF" is excluded.)

   This reconstruction is point-estimate only: it omits the Dirichlet response
   draws and the k posterior of the real fit, so absolute residuals are milder
   than the official ones (2.87 vs 3.62 at the baseline).  The alpha-dependence
   of the *shape* is unaffected, because the observed counts do not depend on
   alpha at all.

2. **Bump-star property audit** — compare B's bin-2 stars (posterior-median
   mass in 3.17--4.0 Msun, baseline branch) against the rest of B's window on
   membership P, A_V, parallax, RUWE and angular position.  Question: do the
   excess stars look like contaminants?

3. **Isochrone caustic location** — on each family's isochrone at each
   subgroup's fitted upper-MS MAP age, find the mass interval where G0(Mini)
   is non-monotonic (the PMS/Henyey fold, where one luminosity maps to several
   masses).  Question: does the fold coincide with B's worst bin at B's age,
   and move out of the window at A's and C's older ages?

Outputs: tables/wp5_bump_shape_diagnostic.csv (alpha scan) and
provenance/wp5_bump_shape_diagnostic_execution.json (everything, with sha256
of every input).  Run:  PYTHONPATH=scripts python3 scripts/wp5_bump_shape_diagnostic.py
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
REPAIR_VERSION = "repair_v3"
BASELINE = ("PARSEC", 3.1)
ALPHA_GRID = np.round(np.arange(1.8, 4.01, 0.1), 2)
BIN_EDGES = np.array([2.0, 2.519842, 3.174802, 4.0, 5.039684, 6.349604, 8.0])
BUMP_BIN = 2  # 3.17--4.0 Msun

INPUTS = {
    "response": PROC / f"wp5_injection_response_{REPAIR_VERSION}.parquet",
    "bins": PROC / f"wp5_mass_function_bins_{REPAIR_VERSION}.parquet",
    "norm": PROC / f"wp5_imf_normalization_{REPAIR_VERSION}.parquet",
    "masses": PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
    "ages": PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet",
    "extinction": PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet",
    "iso_parsec": PROC / "wp3_isochrones_parsec.parquet",
    "iso_mist": PROC / "wp3_isochrones_mist.parquet",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def response_matrix(resp: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """P(observed bin | true mass), including selection, from one branch."""
    draw_cols = [c for c in resp.columns if c.startswith("recovered_mass_draw")]
    masses = np.sort(resp["true_primary_mass"].unique())
    matrix = np.zeros((len(masses), len(BIN_EDGES) - 1))
    for i, mass in enumerate(masses):
        rows = resp[resp["true_primary_mass"] == mass]
        selected = rows[rows["membership_pass"]]
        if selected.empty:
            continue
        draws = selected[draw_cols].to_numpy(float).ravel()
        draws = draws[np.isfinite(draws)]
        if draws.size == 0:
            continue
        hist, _ = np.histogram(draws, bins=BIN_EDGES)
        matrix[i] = hist / draws.size * (len(selected) / len(rows))
    return masses, matrix


def alpha_scan(resp: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (family, rv), branch in resp[resp["subgroup"] == "CygOB2-B"].groupby(
        ["family", "R_V"]
    ):
        masses, matrix = response_matrix(branch)
        weight = np.empty(len(masses))
        weight[0] = 0.5 * (masses[1] - masses[0])
        weight[-1] = 0.5 * (masses[-1] - masses[-2])
        weight[1:-1] = 0.5 * (masses[2:] - masses[:-2])
        observed = (
            bins[
                (bins["subgroup"] == "CygOB2-B")
                & (bins["family"] == family)
                & (bins["R_V"] == rv)
                & (bins["alpha"] == 2.3)
            ]
            .sort_values("bin_index")["membership_weighted_count"]
            .to_numpy()
        )
        for alpha in ALPHA_GRID:
            rate = (matrix * (weight * masses ** (-alpha))[:, None]).sum(axis=0)
            expected = observed.sum() / rate.sum() * rate
            residual = (observed - expected) / np.sqrt(expected)
            rows.append(
                {
                    "family": family,
                    "R_V": rv,
                    "alpha": alpha,
                    "chi2": float((residual**2).sum()),
                    "max_abs_residual": float(np.abs(residual).max()),
                    "bin2_residual": float(residual[BUMP_BIN]),
                    **{f"res_bin{i}": float(r) for i, r in enumerate(residual)},
                }
            )
    return pd.DataFrame(rows)


def bump_star_audit(masses: pd.DataFrame, extinction: pd.DataFrame) -> dict:
    col = "mass_PARSEC_rv3.1"
    joined = masses.merge(
        extinction[
            [
                "source_id",
                "av_rv3.1",
                "parallax_corrected",
                "ruwe",
                "l_deg",
                "b_deg",
            ]
        ],
        on="source_id",
    )
    b = joined[joined["subgroup"] == "CygOB2-B"]
    window = b[(b[col] >= 2.0) & (b[col] < 8.0)]
    in_bump = (window[col] >= BIN_EDGES[BUMP_BIN]) & (
        window[col] < BIN_EDGES[BUMP_BIN + 1]
    )
    centre_l, centre_b = b["l_deg"].median(), b["b_deg"].median()

    def summary(sample: pd.DataFrame) -> dict:
        return {
            "n": int(len(sample)),
            "membership_p_median": float(sample["membership_probability"].median()),
            "av_median": float(sample["av_rv3.1"].median()),
            "parallax_median_mas": float(sample["parallax_corrected"].median()),
            "ruwe_median": float(sample["ruwe"].median()),
            "dist_to_B_centroid_median_deg": float(
                np.hypot(
                    sample["l_deg"] - centre_l, sample["b_deg"] - centre_b
                ).median()
            ),
        }

    return {
        "baseline_mass_column": col,
        "bump_bin_Msun": [float(BIN_EDGES[BUMP_BIN]), float(BIN_EDGES[BUMP_BIN + 1])],
        "bin2_stars": summary(window[in_bump]),
        "rest_of_window": summary(window[~in_bump]),
    }


def fold_interval(iso: pd.DataFrame, age_myr: float) -> dict:
    """Mass interval where G0(Mini) is non-monotonic (brightens with mass)."""
    ages = np.sort(iso["age_Myr"].unique())
    nearest = float(ages[np.argmin(np.abs(ages - age_myr))])
    seq = (
        iso[(iso["age_Myr"] == nearest) & iso["Mini"].between(1.5, 6.0)]
        .sort_values("Mini")[["Mini", "G0"]]
        .to_numpy()
    )
    rising = np.diff(seq[:, 1]) > 0  # G0 increasing with mass = inverted mapping
    amplitude = float(np.sum(np.diff(seq[:, 1])[rising]))
    if not rising.any() or amplitude < 0.05:
        return {"age_used_Myr": nearest, "fold": None, "fold_amplitude_mag": amplitude}
    masses = seq[1:, 0][rising]
    return {
        "age_used_Myr": nearest,
        "fold": [float(masses.min()), float(masses.max())],
        "fold_amplitude_mag": amplitude,
    }


def main() -> None:
    resp = pd.read_parquet(INPUTS["response"])
    bins = pd.read_parquet(INPUTS["bins"])
    norm = pd.read_parquet(INPUTS["norm"])
    masses = pd.read_parquet(INPUTS["masses"])
    ages = pd.read_parquet(INPUTS["ages"])
    extinction = pd.read_parquet(INPUTS["extinction"])

    scan = alpha_scan(resp, bins)
    table_path = ROOT / "tables" / "wp5_bump_shape_diagnostic.csv"
    scan.to_csv(table_path, index=False)

    scan_summary = []
    for (family, rv), grp in scan.groupby(["family", "R_V"]):
        best = grp.loc[grp["chi2"].idxmin()]
        official = norm[
            (norm["subgroup"] == "CygOB2-B")
            & (norm["family"] == family)
            & (norm["R_V"] == rv)
            & (norm["alpha"] == 2.3)
        ]["max_abs_pearson_residual"].iloc[0]
        scan_summary.append(
            {
                "family": family,
                "R_V": float(rv),
                "official_maxres_alpha2.3": float(official),
                "best_alpha": float(best["alpha"]),
                "chi2_at_best_alpha": float(best["chi2"]),
                "bin2_residual_at_best_alpha": float(best["bin2_residual"]),
                "min_bin2_residual_over_all_alpha": float(grp["bin2_residual"].min()),
            }
        )

    map_ages = (
        ages[
            (ages["indicator"] == "ums")
            & (ages["dmu"] == 0.0)
            & (ages["f_bin"] == 0.4)
            & (ages["family"] == BASELINE[0])
            & (ages["R_V"] == BASELINE[1])
        ]
        .set_index("subgroup")["age_map"]
        .to_dict()
    )
    isochrones = {
        "PARSEC": pd.read_parquet(INPUTS["iso_parsec"]),
        "MIST": pd.read_parquet(INPUTS["iso_mist"]),
    }
    caustic = {
        family: {
            subgroup: fold_interval(iso, age)
            for subgroup, age in sorted(map_ages.items())
        }
        for family, iso in isochrones.items()
    }

    record = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": "scripts/wp5_bump_shape_diagnostic.py",
        "repair_version": REPAIR_VERSION,
        "issue": "#1c CygOB2-B mass-function residual — shape diagnostic",
        "deterministic": True,
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in INPUTS.values()},
        "method_caveat": (
            "Point-estimate reconstruction of the forward rates; omits the "
            "Dirichlet response draws and k posterior of wp5_fit_imf.py, so "
            "absolute residuals are milder than the official gate values. "
            "Alpha-dependence of the shape is unaffected."
        ),
        "alpha_scan": {
            "grid": [float(a) for a in ALPHA_GRID],
            "per_branch": scan_summary,
            "table": str(table_path.relative_to(ROOT)),
        },
        "bump_star_audit": bump_star_audit(masses, extinction),
        "upper_ms_map_ages_Myr_baseline": {
            k: float(v) for k, v in sorted(map_ages.items())
        },
        "isochrone_caustic_G0_fold": caustic,
        "observed_bin_edges_Msun": [float(e) for e in BIN_EDGES],
    }

    out = ROOT / "provenance" / "wp5_bump_shape_diagnostic_execution.json"
    out.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["alpha_scan"]["per_branch"], indent=2))
    print(json.dumps(record["isochrone_caustic_G0_fold"], indent=2))
    print(f"wrote {table_path.relative_to(ROOT)} and {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
