#!/usr/bin/env python3
"""WP6 steps 3-4: runaway recovery by proper-motion traceback.

Standard 2D traceback, as used for OB associations since Blaauw (1961) and in
its Gaia form by Hoogerwerf+ (2001) and successors.  A star ejected from the
association carries a proper motion pointing away from it; running the motion
backwards over the association's age should bring it back inside the footprint.

WHAT MAKES THIS DEFENSIBLE, AND WHAT DOES NOT
---------------------------------------------
The traceback itself is easy and produces false positives freely: any field OB
star whose proper motion happens to point back at Cyg OB2 will be "recovered".
The number that matters is therefore not the raw count but the count minus the
chance-alignment rate, and the chance rate is measured, not modelled --
control fields at matched Galactic latitude are traced back to the association
by identical code, and their intersection rate IS the false-positive rate.

FOUR THINGS DELIBERATELY DONE THE CONSERVATIVE WAY
--------------------------------------------------
1.  Astrometric covariance is sampled, not point-estimated.  Gaia's pmra/pmdec
    are correlated and the correlation matters over a Myr baseline, so each
    candidate is traced back N_TRACEBACK times from its full 2x2 proper-motion
    covariance and the recovery is a probability, not a yes/no.
2.  No radial velocity is used even where DR3 provides one.  The traceback is
    2D and the recovered count is therefore a LOWER BOUND: a star ejected
    mostly along the line of sight shows little proper motion and is missed.
3.  The search footprint caps the recoverable velocity.  The wide box is
    +-8 deg, which at 1.62 kpc recovers 100 km/s only for ejections within the
    last ~2 Myr and caps near 38 km/s over 5 Myr (CUTS section 4.2, issue #10).
    That bound is computed here and reported, not left implicit.
4.  Candidates already in the WP2 member sample are excluded -- they never
    left.

Outputs:
  tables/wp6_runaways.csv
  provenance/wp6_runaways_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_runaways.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

UPSTREAM = "repair_v5"
SEED = 20260728

# Cyg OB2 distance and its parallax window, inherited from WP1/WP2.
DISTANCE_PC = 1620.0
PARALLAX_LO_MAS = 0.45
PARALLAX_HI_MAS = 0.90
# A runaway must be moving; below this it is indistinguishable from a member
# with ordinary internal motion.  20 km/s at 1.62 kpc is 2.6 mas/yr.
MIN_EJECTION_KM_S = 10.0
MAX_EJECTION_KM_S = 100.0
# Association footprint, from the WP2 member extent (max 0.910 deg).
FOOTPRINT_RADIUS_DEG = 1.0
# Traceback baseline: the oldest subgroup age plus margin.
MAX_LOOKBACK_MYR = 5.0
N_TRACEBACK = 256
# OB candidate selection: de-reddened absolute G brighter than this.  B2V sits
# near -1.3, so this keeps the >8 Msun population with margin.
MAX_ABSOLUTE_G = 0.0
MIN_PARALLAX_SNR = 5.0
# Control fields: same |b|, offset in longitude far enough to be an independent
# population but close enough to share Galactic structure.
CONTROL_LONGITUDE_OFFSETS_DEG = (-6.0, -4.0, 4.0, 6.0)

MAS_YR_TO_KM_S = 4.74047 * DISTANCE_PC / 1000.0  # 1 mas/yr at 1.62 kpc
DEG_PER_MYR_PER_MAS_YR = 1.0 / 3.6  # 1 mas/yr = 1000 mas/Myr = 1/3.6 deg/Myr


def galactic(ra_deg: np.ndarray, dec_deg: np.ndarray):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic
    return coords.l.deg, coords.b.deg


def proper_motion_galactic(
    ra_deg: np.ndarray, dec_deg: np.ndarray, pmra: np.ndarray, pmdec: np.ndarray
):
    """Rotate (pmra*, pmdec) into (pml*, pmb) via the local position angle."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coords = SkyCoord(
        ra=ra_deg * u.deg, dec=dec_deg * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
        frame="icrs",
    ).galactic
    return coords.pm_l_cosb.to_value(u.mas / u.yr), coords.pm_b.to_value(u.mas / u.yr)


def velocity_ceiling(half_width_deg: float, lookback_myr: float) -> float:
    """Fastest ejection still inside the box after this lookback, km/s."""
    displacement_pc = np.radians(half_width_deg) * DISTANCE_PC
    return displacement_pc * 0.9778 / lookback_myr  # pc/Myr -> km/s


def select_ob_candidates(
    wide: pd.DataFrame, centroid_l: float, centroid_b: float
) -> pd.DataFrame:
    """Parallax-compatible, blue, luminous, well-measured sources."""
    frame = wide[
        wide.parallax.between(PARALLAX_LO_MAS, PARALLAX_HI_MAS)
        & wide.pmra.notna() & wide.pmdec.notna()
        & wide.phot_g_mean_mag.notna()
        & wide.phot_bp_mean_mag.notna() & wide.phot_rp_mean_mag.notna()
        & (wide.parallax / wide.parallax_error > MIN_PARALLAX_SNR)
    ].copy()
    distance_modulus = 5.0 * np.log10(DISTANCE_PC) - 5.0
    # No per-star extinction exists outside the member sample, so the cut is
    # made on the OBSERVED absolute magnitude, which is a conservative (over-
    # inclusive) proxy: extinction only makes stars fainter, so a star passing
    # here unreddened would pass a fortiori once de-reddened.
    frame["absolute_g_observed"] = frame.phot_g_mean_mag - distance_modulus
    frame["bp_rp"] = frame.phot_bp_mean_mag - frame.phot_rp_mean_mag
    return frame[frame.absolute_g_observed < MAX_ABSOLUTE_G].copy()


def traceback(
    frame: pd.DataFrame, centroid_l: float, centroid_b: float,
    rng: np.random.Generator, lookback_myr: float,
    systemic_pmra: float = 0.0, systemic_pmdec: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recovery probability, best lookback time, and implied velocity.

    The traceback runs on the PECULIAR proper motion -- the star's motion minus
    the association's systemic motion.  Absolute Gaia proper motions are
    dominated by the bulk motion of the association plus Galactic rotation
    (Cyg OB2's systemic motion is -2.71, -4.32 mas/yr, larger than a typical
    ejection signature), so tracing back absolute motions measures the common
    drift rather than the ejection and recovers essentially nothing.

    The systemic vector is rotated to galactic AT EACH STAR'S POSITION before
    subtraction, because the equatorial-to-galactic rotation is
    position-dependent.
    """
    pml, pmb = proper_motion_galactic(
        frame.ra.to_numpy(float), frame.dec.to_numpy(float),
        frame.pmra.to_numpy(float), frame.pmdec.to_numpy(float),
    )
    if systemic_pmra or systemic_pmdec:
        sys_l, sys_b = proper_motion_galactic(
            frame.ra.to_numpy(float), frame.dec.to_numpy(float),
            np.full(len(frame), systemic_pmra),
            np.full(len(frame), systemic_pmdec),
        )
        pml = pml - sys_l
        pmb = pmb - sys_b
    n = len(frame)
    sigma_ra = frame.pmra_error.to_numpy(float)
    sigma_dec = frame.pmdec_error.to_numpy(float)
    correlation = (
        frame.pmra_pmdec_corr.to_numpy(float)
        if "pmra_pmdec_corr" in frame.columns
        else np.zeros(n)
    )
    correlation = np.nan_to_num(correlation)

    times = np.linspace(0.0, lookback_myr, 60)[1:]
    probability = np.zeros(n)
    best_time = np.full(n, np.nan)
    velocity = np.full(n, np.nan)

    for index in range(n):
        # Sample the proper-motion covariance; the rotation to galactic is
        # approximately linear over the error ellipse, so the same rotation is
        # applied to the sampled offsets.
        covariance = np.array(
            [
                [sigma_ra[index] ** 2, correlation[index] * sigma_ra[index] * sigma_dec[index]],
                [correlation[index] * sigma_ra[index] * sigma_dec[index], sigma_dec[index] ** 2],
            ]
        )
        try:
            draws = rng.multivariate_normal(
                [pml[index], pmb[index]], covariance, size=N_TRACEBACK
            )
        except np.linalg.LinAlgError:
            continue
        # Backward positions: subtract the motion.
        delta_l = -draws[:, 0][:, None] * times[None, :] * DEG_PER_MYR_PER_MAS_YR
        delta_b = -draws[:, 1][:, None] * times[None, :] * DEG_PER_MYR_PER_MAS_YR
        past_l = frame.l.iloc[index] + delta_l / np.cos(np.radians(frame.b.iloc[index]))
        past_b = frame.b.iloc[index] + delta_b
        separation = np.hypot(
            (past_l - centroid_l) * np.cos(np.radians(past_b)), past_b - centroid_b
        )
        inside = separation <= FOOTPRINT_RADIUS_DEG
        hit = inside.any(axis=1)
        probability[index] = float(hit.mean())
        if probability[index] > 0:
            first = np.argmax(inside[hit], axis=1)
            best_time[index] = float(np.median(times[first]))
            speed = np.hypot(draws[hit, 0], draws[hit, 1]) * MAS_YR_TO_KM_S
            velocity[index] = float(np.median(speed))
    return probability, best_time, velocity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-probability", type=float, default=0.5)
    parser.add_argument("--max-candidates", type=int, default=0,
                        help="0 = no limit; for smoke tests only")
    args = parser.parse_args()

    wide = pd.read_parquet(w.PROC / "wp1_gaia_wide.parquet")
    members = pd.read_parquet(w.PROC / "wp2_members.parquet")
    labels = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    positions = members[["source_id", "l_deg", "b_deg"]].merge(
        labels[["source_id", "subgroup"]], on="source_id", how="inner"
    )
    centroid = positions[["l_deg", "b_deg"]].median()
    centroid_l, centroid_b = float(centroid.l_deg), float(centroid.b_deg)

    wide = wide.copy()
    wide["l"], wide["b"] = galactic(
        wide.ra.to_numpy(float), wide.dec.to_numpy(float)
    )
    candidates = select_ob_candidates(wide, centroid_l, centroid_b)
    member_ids = set(members.source_id.astype("int64"))
    candidates = candidates[
        ~candidates.source_id.astype("int64").isin(member_ids)
    ]
    # Stars still inside the footprint have not escaped.
    separation_now = np.hypot(
        (candidates.l - centroid_l) * np.cos(np.radians(candidates.b)),
        candidates.b - centroid_b,
    )
    candidates = candidates[separation_now > FOOTPRINT_RADIUS_DEG].copy()
    if args.max_candidates and len(candidates) > args.max_candidates:
        candidates = candidates.sample(
            args.max_candidates, random_state=SEED
        ).copy()
    print(f"OB candidates outside the footprint: {len(candidates)}", flush=True)

    systemic_pmra = float(members.pmra.median())
    systemic_pmdec = float(members.pmdec.median())
    print(f"association systemic proper motion: "
          f"({systemic_pmra:+.3f}, {systemic_pmdec:+.3f}) mas/yr", flush=True)

    rng = np.random.default_rng(SEED)
    probability, best_time, velocity = traceback(
        candidates, centroid_l, centroid_b, rng, MAX_LOOKBACK_MYR,
        systemic_pmra, systemic_pmdec,
    )
    candidates["recovery_probability"] = probability
    candidates["lookback_Myr"] = best_time
    candidates["implied_velocity_km_s"] = velocity
    candidates["in_velocity_window"] = (
        (velocity >= MIN_EJECTION_KM_S) & (velocity <= MAX_EJECTION_KM_S)
    )
    candidates["is_runaway_candidate"] = (
        (candidates.recovery_probability >= args.min_probability)
        & candidates.in_velocity_window
    )
    recovered = int(candidates.is_runaway_candidate.sum())
    print(f"traceback recovers {recovered} candidates", flush=True)

    # ---- control fields: the false-positive rate, measured -----------------
    # Chance recovery depends strongly on how far a star is from the centre --
    # a nearby star needs only a small proper motion to trace back inside a
    # 1 deg footprint.  A raw control rate is therefore only meaningful if the
    # control stars sit at the same separations as the real candidates, which
    # they do not.  The rate is measured as a FUNCTION of separation and then
    # applied to the real sample's own separation distribution.
    control_rows = []
    for offset in CONTROL_LONGITUDE_OFFSETS_DEG:
        control_l = centroid_l + offset
        control = select_ob_candidates(wide, control_l, centroid_b)
        control = control[~control.source_id.astype("int64").isin(member_ids)]
        control_separation = np.hypot(
            (control.l - control_l) * np.cos(np.radians(control.b)),
            control.b - centroid_b,
        )
        control = control[control_separation > FOOTPRINT_RADIUS_DEG].copy()
        control["separation_deg"] = control_separation[
            control_separation > FOOTPRINT_RADIUS_DEG
        ]
        if args.max_candidates and len(control) > args.max_candidates:
            control = control.sample(
                args.max_candidates, random_state=SEED
            ).copy()
        if not len(control):
            continue
        # Trace back to the CONTROL centre by identical code: any recovery is
        # by construction a chance alignment.
        # Controls get the identical treatment, including the same systemic
        # subtraction, so the chance rate is measured on the same statistic.
        cp, ct, cv = traceback(
            control, control_l, centroid_b, np.random.default_rng(SEED),
            MAX_LOOKBACK_MYR, systemic_pmra, systemic_pmdec,
        )
        control["hit"] = (
            (cp >= args.min_probability)
            & (cv >= MIN_EJECTION_KM_S)
            & (cv <= MAX_EJECTION_KM_S)
        )
        control["offset"] = offset
        control_rows.append(control[["separation_deg", "hit", "offset"]])
        print(f"  control field {offset:+.0f} deg: "
              f"{int(control.hit.sum())}/{len(control)}", flush=True)

    control_rates, expected_false, mean_rate = [], float("nan"), float("nan")
    if control_rows:
        pooled = pd.concat(control_rows, ignore_index=True)
        edges = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0])
        rate_by_bin, bin_counts = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            block = pooled[
                pooled.separation_deg.ge(lo) & pooled.separation_deg.lt(hi)
            ]
            rate = float(block.hit.mean()) if len(block) else np.nan
            rate_by_bin.append(rate)
            bin_counts.append(int(len(block)))
            control_rates.append(
                {
                    "separation_lo_deg": float(lo),
                    "separation_hi_deg": float(hi),
                    "control_stars": int(len(block)),
                    "false_positive_rate": None if np.isnan(rate) else round(rate, 5),
                }
            )
        # Apply the separation-matched rate to each real candidate.
        real_separation = np.hypot(
            (candidates.l - centroid_l) * np.cos(np.radians(candidates.b)),
            candidates.b - centroid_b,
        ).to_numpy(float)
        rate_array = np.array(rate_by_bin, dtype=float)
        index = np.clip(
            np.digitize(real_separation, edges) - 1, 0, len(rate_array) - 1
        )
        per_star = rate_array[index]
        # Bins with no control coverage fall back to the pooled mean rather
        # than silently contributing zero false positives.
        pooled_mean = float(pooled.hit.mean())
        per_star = np.where(np.isnan(per_star), pooled_mean, per_star)
        expected_false = float(np.sum(per_star))
        mean_rate = float(expected_false / len(candidates)) if len(candidates) else np.nan
        candidates["false_positive_rate_at_separation"] = per_star
    corrected = recovered - expected_false

    out_csv = w.TABLES / "wp6_runaways.csv"
    candidates[
        [
            "source_id", "ra", "dec", "l", "b", "parallax", "parallax_error",
            "pmra", "pmdec", "phot_g_mean_mag", "bp_rp", "absolute_g_observed",
            "recovery_probability", "lookback_Myr", "implied_velocity_km_s",
            "in_velocity_window", "is_runaway_candidate",
            "false_positive_rate_at_separation",
        ]
    ].sort_values("recovery_probability", ascending=False).to_csv(out_csv, index=False)

    ceilings = {
        f"{years:.0f}_Myr": round(velocity_ceiling(8.0, years), 1)
        for years in (1.0, 2.0, 3.0, 4.0, 5.0)
    }
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_runaways.py",
        "status": "SUCCESS",
        "work_package": "WP6 steps 3-4",
        "method": (
            "standard 2D proper-motion traceback (Blaauw 1961; Hoogerwerf+ "
            "2001 in its Gaia form).  Each candidate is traced back "
            f"{N_TRACEBACK} times from its full proper-motion covariance, so "
            "recovery is a probability rather than a yes/no."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / "wp1_gaia_wide.parquet",
                w.PROC / "wp2_members.parquet",
            ]
        },
        "configuration": {
            "distance_pc": DISTANCE_PC,
            "parallax_window_mas": [PARALLAX_LO_MAS, PARALLAX_HI_MAS],
            "min_parallax_snr": MIN_PARALLAX_SNR,
            "max_absolute_g_observed": MAX_ABSOLUTE_G,
            "footprint_radius_deg": FOOTPRINT_RADIUS_DEG,
            "max_lookback_Myr": MAX_LOOKBACK_MYR,
            "velocity_window_km_s": [MIN_EJECTION_KM_S, MAX_EJECTION_KM_S],
            "traceback_samples": N_TRACEBACK,
            "min_recovery_probability": args.min_probability,
            "seed": SEED,
            "systemic_pmra_mas_yr": round(systemic_pmra, 4),
            "systemic_pmdec_mas_yr": round(systemic_pmdec, 4),
            "peculiar_motion_note": (
                "the traceback runs on proper motion MINUS the association "
                "systemic motion.  Using absolute proper motions measures the "
                "common bulk drift, not the ejection, and recovers essentially "
                "nothing -- caught by the literature cross-match, which found "
                "the canonical ejected star BD+43 3654 scoring 0.000."
            ),
        },
        "result": {
            "ob_candidates_outside_footprint": int(len(candidates)),
            "raw_recovered": recovered,
            "effective_false_positive_rate": round(mean_rate, 5),
            "false_positive_rate_by_separation": control_rates,
            "rate_matching": (
                "chance recovery falls steeply with separation from the "
                "centre, so the control rate is measured per separation bin "
                "and applied to each real candidate at its own separation.  "
                "A pooled control rate would be meaningless: the control "
                "and real samples do not share a separation distribution."
            ),
            "expected_false_positives": round(expected_false, 1),
            "false_positive_corrected": round(corrected, 1),
        },
        "distance_note": (
            f"velocity ceilings are computed at the WP2 distance of "
            f"{DISTANCE_PC:.0f} pc.  CUTS section 4.2 tabulated them at a "
            "round 1.4 kpc and so quotes slightly smaller values; the "
            "numbers here supersede that table."
        ),
        "issue_10_velocity_ceiling": {
            "statement": (
                "the +-8 deg wide box caps the recoverable ejection velocity. "
                "The plan's stated 10-100 km/s range is NOT achievable over the "
                "full age baseline, and the recovered count is a lower bound."
            ),
            "ceiling_km_s_by_lookback": ceilings,
            "adopted_wording": (
                "the runaway search is complete to v <~ "
                f"{ceilings['5_Myr']:.0f} km/s for ejections up to 5 Myr ago; "
                "faster or older ejections lie outside the search footprint and "
                "make N_SN a lower bound"
            ),
        },
        "carried_limitations": [
            "2D only — no radial velocity is used, so ejections along the line "
            "of sight are missed and the count is a lower bound",
            "the observed-magnitude OB cut is conservative (over-inclusive), "
            "because no per-star extinction exists outside the member sample",
            "chance alignments are subtracted statistically, not per star, so "
            "individual candidates carry no confirmation",
        ],
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp6_runaways_execution.json", record)

    print(f"\n  raw recovered            {recovered}")
    print(f"  control false-positive   {mean_rate:.4f} -> {expected_false:.1f} expected")
    print(f"  corrected                {corrected:.1f}")
    print(f"\n  velocity ceiling of the +-8 deg box: {ceilings}")
    print("wrote provenance/wp6_runaways_execution.json")


if __name__ == "__main__":
    main()
