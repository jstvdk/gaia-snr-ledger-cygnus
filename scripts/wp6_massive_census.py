#!/usr/bin/env python3
"""WP6 step 0b: the observed living massive-star census.

Two things this fixes before any closure ratio is computed.

COUNTING CONVENTION.  Thresholding point-estimate masses at 8 Msun is
Eddington-biased: the mass function is steep, so more stars scatter up across
the threshold than down.  WP4 already carries P(M > 8 Msun) per star per
branch, so the count is

    N_obs = sum over members of  p_membership * P(M > 8)

which is unbiased and carries the mass uncertainty instead of discarding it.

WHAT IS COMPARED TO WHAT.  The injection response models recovery through the
WP2 pipeline, so the quantity it predicts is the count of stars that PASSED
that pipeline -- the member sample.  The comparison is therefore

    predicted observed  (forward-modelled)   vs   member count

and the spectroscopic anchors that are absent from the member sample are kept
as a SEPARATE and independent check: the response predicts how many massive
stars the pipeline should lose, and those anchors are how many it actually
lost.  Folding them into the same number would double-correct.

62 of the 252 countable anchors are absent from the member sample -- 30 O
stars, 13 B stars, 6 Wolf-Rayets, 21 of them supergiants.  They are bright
stars that failed the WP2 quality filter, which is exactly the loss channel the
response describes.

Outputs:
  tables/wp6_massive_census.csv           per-branch observed counts
  tables/wp6_orphan_anchors.csv           anchors absent from the member sample
  provenance/wp6_massive_census_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_massive_census.py
"""
from __future__ import annotations

import platform
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

UPSTREAM = "repair_v5"
SN_THRESHOLD_MSUN = 8.0
# Association footprint: the member sample reaches 0.910 deg from its own
# subgroup centroid at most, 0.774 deg at the 99th percentile.  An anchor
# further out than the members themselves is not a Cyg OB2 star being
# missed -- the field contains Cyg OB9 and other associations -- so it is
# reported separately rather than assigned to the nearest subgroup.
FOOTPRINT_RADIUS_DEG = 1.0

# Spectral type -> initial mass, Msun.  Martins, Schaerer & Hillier (2005)
# theoretical scale for O stars; conventional main-sequence and supergiant
# masses for B.  These are used ONLY for anchors absent from the member sample,
# and only to decide whether a star is above 8 Msun and to place it in the
# closure window -- every O star and every Wolf-Rayet clears 8 Msun regardless
# of calibration detail, so the census count is far more robust than the
# individual masses.
O_DWARF = {3.0: 58, 4.0: 47, 5.0: 40, 5.5: 37, 6.0: 34, 6.5: 31, 7.0: 28,
           7.5: 26, 8.0: 24, 8.5: 22, 9.0: 20, 9.5: 19}
O_GIANT = {3.0: 62, 4.0: 55, 5.0: 49, 6.0: 43, 7.0: 37, 8.0: 32, 9.0: 27,
           9.5: 25}
O_SUPERGIANT = {3.0: 68, 4.0: 62, 5.0: 56, 6.0: 50, 7.0: 44, 8.0: 38, 9.0: 33,
                9.5: 31}
B_DWARF = {0.0: 17, 0.5: 15, 1.0: 13, 2.0: 10, 3.0: 7.0, 5.0: 5.5, 8.0: 3.8,
           9.0: 3.2}
B_SUPERGIANT = {0.0: 28, 1.0: 24, 2.0: 20, 3.0: 17, 5.0: 14, 8.0: 11, 9.0: 10}
# Wolf-Rayets are evolved massive stars; their initial masses are model
# dependent but unambiguously well above the core-collapse threshold.
WR_INITIAL_MASS = 40.0

# The OC/ON/BC/BN carbon- and nitrogen-rich subclasses carry a second
# letter before the subtype ("OC9.7Iab"); it does not change the mass scale.
TYPE_PATTERN = re.compile(r"^\s*(?P<cls>[OB])[CN]?\s*(?P<sub>\d+(?:\.\d+)?)")
WR_PATTERN = re.compile(r"^\s*(?:WR|WN|WC|WO)", re.IGNORECASE)
# Anchored at the start of the post-subtype tail, because peculiarity codes
# follow the luminosity class directly: "O7Ib(f)", "O4If+", "B1Ia", "O4V:n".
# IV must precede the bare-I alternative or "O9IV" reads as class I.
LUM_PATTERN = re.compile(r"^[\s:.]*(IV|III|II|I[ab]{0,2}|V)")


def interpolate_class(table: dict[float, float], subtype: float) -> float:
    keys = np.array(sorted(table))
    values = np.array([table[k] for k in keys])
    return float(np.interp(subtype, keys, values))


def spectral_type_mass(spectral_type: str) -> tuple[float, str]:
    """Initial mass in Msun and the rule used, or (nan, reason)."""
    if not isinstance(spectral_type, str) or not spectral_type.strip():
        return float("nan"), "no_spectral_type"
    text = spectral_type.strip()
    if WR_PATTERN.match(text):
        return WR_INITIAL_MASS, "wolf_rayet_nominal"
    match = TYPE_PATTERN.match(text)
    if not match:
        return float("nan"), "unparsed"
    cls = match.group("cls")
    subtype = float(match.group("sub"))
    # Look for the luminosity class only after the spectral subtype.
    tail = text[match.end():]
    lum = LUM_PATTERN.match(tail)
    luminosity = lum.group(1) if lum else "V"
    # Ia / Iab / Ib are all luminosity class I.
    if luminosity.startswith("I") and luminosity not in ("IV", "III", "II"):
        luminosity = "I"
    if cls == "O":
        if luminosity in ("I", "II"):
            return interpolate_class(O_SUPERGIANT, subtype), f"O{luminosity}"
        if luminosity in ("III", "IV"):
            return interpolate_class(O_GIANT, subtype), f"O{luminosity}"
        return interpolate_class(O_DWARF, subtype), "OV"
    if luminosity in ("I", "II"):
        return interpolate_class(B_SUPERGIANT, subtype), f"B{luminosity}"
    return interpolate_class(B_DWARF, subtype), "BV"


def assign_subgroup(
    frame: pd.DataFrame, centroids: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest subgroup centroid on the sky, with the separation reported."""
    labels, separations = [], []
    for row in frame.itertuples():
        d = np.hypot(
            (centroids.l_deg.to_numpy() - row.l_deg)
            * np.cos(np.radians(row.b_deg)),
            centroids.b_deg.to_numpy() - row.b_deg,
        )
        index = int(np.argmin(d))
        labels.append(centroids.subgroup.iloc[index])
        separations.append(float(d[index]))
    return np.array(labels), np.array(separations)


def main() -> None:
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet")
    anchors = pd.read_parquet(w.PROC / "wp1_spectroscopic_anchors.parquet")
    members = pd.read_parquet(w.PROC / "wp2_members.parquet")

    # ---- member side: the quantity the injection response predicts ---------
    rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            column = f"mass_{family}_rv{rv:.1f}_p_gt8"
            point = w.mass_column(family, rv)
            for subgroup in w.SUBGROUPS:
                block = masses[masses.subgroup.eq(subgroup)]
                probability = block[column].to_numpy(float)
                weight = block["membership_probability"].to_numpy(float)
                finite = np.isfinite(probability)
                expected = float(np.sum(weight[finite] * probability[finite]))
                thresholded = int(
                    np.sum(block[point].to_numpy(float) > SN_THRESHOLD_MSUN)
                )
                rows.append(
                    {
                        "subgroup": subgroup, "family": family, "R_V": float(rv),
                        "n_members": int(len(block)),
                        "observed_above_8_probabilistic": expected,
                        "observed_above_8_thresholded": thresholded,
                        "eddington_inflation": (
                            thresholded / expected if expected > 0 else np.nan
                        ),
                        "stars_with_finite_p_gt8": int(finite.sum()),
                    }
                )
    census = pd.DataFrame(rows)

    # ---- anchor side: the independent check --------------------------------
    member_ids = set(masses.source_id.astype("int64"))
    countable = anchors[anchors.countable_in_wp6.astype(bool)].copy()
    countable["in_member_sample"] = countable.source_id.astype("int64").isin(
        member_ids
    )
    orphans = countable[~countable.in_member_sample].copy()

    estimated = [spectral_type_mass(value) for value in orphans.spectral_type]
    orphans["initial_mass_Msun"] = [value for value, _ in estimated]
    orphans["mass_rule"] = [rule for _, rule in estimated]
    orphans["above_8_Msun"] = orphans.initial_mass_Msun > SN_THRESHOLD_MSUN

    # Sky positions: anchors carry ra/dec, members carry galactic.
    galactic = members[["source_id", "l_deg", "b_deg"]]
    labelled = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    positions = galactic.merge(
        labelled[["source_id", "subgroup"]], on="source_id", how="inner"
    )
    centroids = (
        positions.groupby("subgroup")[["l_deg", "b_deg"]].median().reset_index()
    )
    coords = SkyCoordLite(orphans.ra_deg.to_numpy(), orphans.dec_deg.to_numpy())
    orphans["l_deg"], orphans["b_deg"] = coords
    label, separation = assign_subgroup(orphans, centroids)
    orphans["subgroup_positional"] = label
    orphans["centroid_separation_deg"] = separation
    orphans["inside_footprint"] = separation <= FOOTPRINT_RADIUS_DEG
    orphans["classifiable"] = orphans.initial_mass_Msun.notna()
    # Only anchors inside the footprint AND with a spectral type can enter
    # the census; the rest are reported as explicit unknowns.
    orphans["counts_in_census"] = (
        orphans.inside_footprint & orphans.classifiable & orphans.above_8_Msun
    )

    orphan_summary = []
    for subgroup in w.SUBGROUPS:
        block = orphans[
            orphans.subgroup_positional.eq(subgroup) & orphans.inside_footprint
        ]
        orphan_summary.append(
            {
                "subgroup": subgroup,
                "orphan_anchors": int(len(block)),
                "orphan_anchors_above_8": int(block.counts_in_census.sum()),
                "unclassifiable": int((~block.classifiable).sum()),
                "median_initial_mass_Msun": (
                    round(float(block.initial_mass_Msun.median()), 1)
                    if len(block) else None
                ),
            }
        )

    out_census = w.TABLES / "wp6_massive_census.csv"
    out_orphans = w.TABLES / "wp6_orphan_anchors.csv"
    census.to_csv(out_census, index=False)
    orphans[
        [
            "anchor_uid", "object_name", "source_id", "spectral_type",
            "initial_mass_Msun", "mass_rule", "above_8_Msun",
            "inside_footprint", "classifiable", "counts_in_census",
            "l_deg", "b_deg", "subgroup_positional", "centroid_separation_deg",
        ]
    ].to_csv(out_orphans, index=False)

    baseline = census[census.family.eq("PARSEC") & census.R_V.eq(3.1)]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_massive_census.py",
        "status": "SUCCESS",
        "work_package": "WP6 step 0b",
        "counting_convention": (
            "N_obs = sum over members of p_membership * P(M > 8 Msun), per "
            "branch.  Thresholding point-estimate masses is Eddington-biased "
            "because the mass function is steep; the inflation factor is "
            "reported per branch so the size of the bias is on the record."
        ),
        "comparison_rule": (
            "the injection response models recovery through the WP2 pipeline, "
            "so it predicts the MEMBER count.  Anchors absent from the member "
            "sample are kept separate as an independent check on the response, "
            "never added to the same number -- that would double-correct."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet",
                w.PROC / "wp1_spectroscopic_anchors.parquet",
                w.PROC / "wp2_members.parquet",
            ]
        },
        "spectral_type_calibration": {
            "source": (
                "Martins, Schaerer & Hillier (2005) theoretical scale for O "
                "stars; conventional masses for B dwarfs and supergiants"
            ),
            "applies_to": "orphan anchors only, never to member stars",
            "robustness": (
                "every O star and every Wolf-Rayet clears 8 Msun regardless of "
                "calibration detail, so the count above threshold is far more "
                "robust than the individual masses.  Only B3 and later fall "
                "below, and the classification there is unambiguous."
            ),
            "wolf_rayet_nominal_initial_mass_Msun": WR_INITIAL_MASS,
        },
        "member_side": {
            "baseline_PARSEC_rv3.1": [
                {
                    "subgroup": row.subgroup,
                    "n_members": row.n_members,
                    "observed_above_8_probabilistic": round(
                        row.observed_above_8_probabilistic, 2
                    ),
                    "observed_above_8_thresholded": row.observed_above_8_thresholded,
                    "eddington_inflation": round(row.eddington_inflation, 3),
                }
                for row in baseline.itertuples()
            ],
            "eddington_inflation_grid_range": [
                round(float(census.eddington_inflation.min()), 3),
                round(float(census.eddington_inflation.max()), 3),
            ],
        },
        "anchor_side": {
            "countable_anchors": int(len(countable)),
            "in_member_sample": int(countable.in_member_sample.sum()),
            "orphans": int(len(orphans)),
            "orphans_above_8_in_footprint": int(orphans.counts_in_census.sum()),
            "orphans_outside_footprint": int((~orphans.inside_footprint).sum()),
            "orphans_unclassifiable": int((~orphans.classifiable).sum()),
            "footprint_radius_deg": FOOTPRINT_RADIUS_DEG,
            "unclassifiable_caveat": (
                "orphans with no spectral type in the anchor table cannot be "
                "placed above or below 8 Msun.  They are counted as unknown, "
                "not as sub-threshold, and bound the orphan census from below."
            ),
            "orphan_composition": {
                "O": int(orphans.spectral_type.astype(str).str.match("^O").sum()),
                "B": int(orphans.spectral_type.astype(str).str.match("^B").sum()),
                "WR": int(
                    orphans.spectral_type.astype(str)
                    .str.match("^W", case=False).sum()
                ),
                "unparsed": int(orphans.mass_rule.eq("unparsed").sum()),
            },
            "by_subgroup": orphan_summary,
            "subgroup_assignment": (
                "nearest WP2 subgroup centroid on the sky; positional only, "
                "flagged as such, and reported with the separation"
            ),
            "interpretation": (
                "these are spectroscopically confirmed massive stars that the "
                "Gaia pipeline lost, mostly bright supergiants failing the WP2 "
                "quality filter.  They are the empirical counterpart of the "
                "response's predicted loss and are used to test it, not to "
                "inflate the census."
            ),
        },
        "outputs": {
            str(out_census.relative_to(w.ROOT)): w.sha256(out_census),
            str(out_orphans.relative_to(w.ROOT)): w.sha256(out_orphans),
        },
    }
    w.write_json(w.PROVENANCE / "wp6_massive_census_execution.json", record)

    print("WP6 step 0b — observed massive-star census\n")
    print("  baseline PARSEC R_V=3.1, stars above 8 Msun:")
    print(f"    {'subgroup':12s} {'probabilistic':>14s} {'thresholded':>12s} {'inflation':>10s}")
    for entry in record["member_side"]["baseline_PARSEC_rv3.1"]:
        print(f"    {entry['subgroup']:12s} "
              f"{entry['observed_above_8_probabilistic']:14.2f} "
              f"{entry['observed_above_8_thresholded']:12d} "
              f"{entry['eddington_inflation']:10.3f}")
    print(f"\n  Eddington inflation across the grid: "
          f"{record['member_side']['eddington_inflation_grid_range'][0]:.3f}"
          f"-{record['member_side']['eddington_inflation_grid_range'][1]:.3f}"
          "  (thresholding overcounts by this factor)")
    anchor = record["anchor_side"]
    print(f"\n  anchors: {anchor['countable_anchors']} countable, "
          f"{anchor['in_member_sample']} in the member sample, "
          f"{anchor['orphans']} orphaned")
    print(f"    orphans above 8 Msun in footprint: "
          f"{anchor['orphans_above_8_in_footprint']}  "
          f"(O {anchor['orphan_composition']['O']}, "
          f"B {anchor['orphan_composition']['B']}, "
          f"WR {anchor['orphan_composition']['WR']}, "
          f"unparsed {anchor['orphan_composition']['unparsed']})")
    print(f"    excluded: {anchor['orphans_outside_footprint']} outside the "
          f"{anchor['footprint_radius_deg']:.1f} deg footprint, "
          f"{anchor['orphans_unclassifiable']} with no spectral type")
    for entry in anchor["by_subgroup"]:
        print(f"    {entry['subgroup']:12s} {entry['orphan_anchors']:3d} orphans, "
              f"{entry['orphan_anchors_above_8']:3d} above 8 Msun, "
              f"{entry['unclassifiable']:2d} unclassifiable")
    print("\nwrote provenance/wp6_massive_census_execution.json")


def SkyCoordLite(ra_deg: np.ndarray, dec_deg: np.ndarray):
    """ICRS -> Galactic without pulling in astropy at import time."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic
    return coords.l.deg, coords.b.deg


if __name__ == "__main__":
    main()
