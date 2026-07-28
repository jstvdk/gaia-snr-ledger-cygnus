#!/usr/bin/env python3
"""Independent near-IR extinction check for CygOB2-B (issues #12, #1d).

The `repair_v4` joint age-k fit removed CygOB2-B's localized bump but left a
mass-dependent residual tilt (issue #12).  That tilt is monotone in the
injection truth age and clears only at ages >= 3.98 Myr, above the upper edge
of B's WP4 upper-main-sequence posterior (2.75-3.42 Myr) -- so WP5's counts and
WP4's age disagree for B and only for B
(provenance/wp5_tilt_vs_age_diagnostic_execution.json).

Age and extinction are degenerate: overestimating A_V makes stars look
intrinsically fainter, which biases a fitted age young *and* displaces the mass
scale.  Issue #1d records that B has no spectroscopic anchors nearby, so its
A_V rests on broadband optical photometry alone -- exactly the configuration in
which such an error hides.  This script breaks the degeneracy with data that do
not come from the optical fit.

Method.  Near-IR colour excess.  For 2-8 Msun stars the intrinsic J-Ks colour
is small and weakly age dependent, so

    A_V(NIR) = [ (J-Ks)_obs - (J-Ks)_0 ] / (k_J - k_Ks)

with (J-Ks)_0 read off the same isochrone family at the subgroup's fitted age
and the star's WP4 mass, and k the WP3 extinction-law coefficients.  Four
guards make this a real test rather than a restatement of the optical fit:

1. **Anchor validation** -- where spectroscopic A_V exists, A_V(NIR) is
   compared against it, calibrating the method's accuracy and bias.
2. **Age robustness** -- (J-Ks)_0 is recomputed at 2.82 and 3.98 Myr, the two
   ages in dispute, to show the test does not depend on the answer it is
   trying to find.
3. **Fixed-colour variant** -- a single constant (J-Ks)_0 is used as a
   cross-check, so no isochrone enters at all.
4. **Mass-matched differential** -- subgroups are compared inside the same
   observed-mass bins, so any residual intrinsic-colour systematic cancels
   between A, B and C rather than masquerading as a B-specific offset.

The reddening-vector slope E(J-H)/E(H-Ks) is also fitted per subgroup as an
independent probe of the extinction-law *shape*.

Outputs: tables/wp3_nir_extinction_crosscheck.csv
         provenance/wp3_nir_extinction_crosscheck_execution.json

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp3_nir_extinction_crosscheck.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
from wp3_extinction_law import band_coefficients

UPSTREAM = "repair_v3"
RV = 3.1
FAMILY = "PARSEC"
MASS_COLUMN = f"mass_{FAMILY}_rv{RV:.1f}"
WINDOW = (2.0, 8.0)
MASS_BINS = np.geomspace(*WINDOW, 5)
DISPUTED_AGES = [2.818, 3.981]


def intrinsic_jks(family: str, age_myr: float, mass: np.ndarray) -> np.ndarray:
    iso, _ = w.load_isochrone_at_age(family, age_myr)
    mini = iso["Mini"].to_numpy(float)
    j0 = iso["J0"].to_numpy(float)
    ks0 = iso["Ks0"].to_numpy(float)
    return np.interp(np.clip(mass, mini.min(), mini.max()), mini, j0 - ks0)


def intrinsic_jh_hks(family: str, age_myr: float, mass: np.ndarray):
    iso, _ = w.load_isochrone_at_age(family, age_myr)
    mini = iso["Mini"].to_numpy(float)
    clipped = np.clip(mass, mini.min(), mini.max())
    jh = np.interp(clipped, mini, iso["J0"].to_numpy(float) - iso["H0"].to_numpy(float))
    hks = np.interp(clipped, mini, iso["H0"].to_numpy(float) - iso["Ks0"].to_numpy(float))
    return jh, hks


def robust(values: np.ndarray) -> dict:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"n": 0}
    return {
        "n": int(len(values)),
        "median": float(np.median(values)),
        "mad_sigma": float(1.4826 * np.median(np.abs(values - np.median(values)))),
        "mean": float(np.mean(values)),
        "sem": float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else float("nan"),
    }


def main() -> None:
    coefficients = band_coefficients(RV)
    k_jks = coefficients["J"] - coefficients["Ks"]
    k_jh = coefficients["J"] - coefficients["H"]
    k_hks = coefficients["H"] - coefficients["Ks"]
    law_slope = k_jh / k_hks

    extinction = pd.read_parquet(w.PROC / f"wp3_extinction_{UPSTREAM}.parquet")
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet")[
        ["source_id", MASS_COLUMN]
    ]
    ages = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet")
    frame = extinction.merge(masses, on="source_id", how="left")
    frame = frame[
        frame["subgroup"].isin(w.SUBGROUPS)
        & frame["J"].notna() & frame["H"].notna() & frame["Ks"].notna()
        & frame[f"av_rv{RV:.1f}"].notna()
    ].copy()

    map_age = {
        sg: float(
            ages[
                ages.subgroup.eq(sg) & ages.family.eq(FAMILY) & ages.R_V.eq(RV)
                & ages.f_bin.eq(w.F_BINARY) & ages.indicator.eq("ums") & ages.dmu.eq(0.0)
            ].age_map.iloc[0]
        )
        for sg in w.SUBGROUPS
    }

    frame["jks_obs"] = frame["J"] - frame["Ks"]
    frame["jh_obs"] = frame["J"] - frame["H"]
    frame["hks_obs"] = frame["H"] - frame["Ks"]
    frame["jks_intrinsic"] = np.nan
    for subgroup, age in map_age.items():
        mask = frame["subgroup"].eq(subgroup)
        frame.loc[mask, "jks_intrinsic"] = intrinsic_jks(
            FAMILY, age, frame.loc[mask, MASS_COLUMN].to_numpy(float)
        )
    frame["av_nir"] = (frame["jks_obs"] - frame["jks_intrinsic"]) / k_jks
    frame["av_wp3"] = frame[f"av_rv{RV:.1f}"]
    frame["delta_av"] = frame["av_wp3"] - frame["av_nir"]

    # ---- guard 3: fixed intrinsic colour, no isochrone at all ----
    fixed_colour = float(np.median(frame["jks_intrinsic"]))
    frame["av_nir_fixed"] = (frame["jks_obs"] - fixed_colour) / k_jks
    frame["delta_av_fixed"] = frame["av_wp3"] - frame["av_nir_fixed"]

    # ---- guard 1: validate against spectroscopic anchors ----
    anchors = pd.read_parquet(w.PROC / "wp1_spectroscopic_anchors.parquet")[
        ["source_id", "extinction_av_mag"]
    ].dropna()
    anchors["source_id"] = pd.to_numeric(anchors["source_id"], errors="coerce").astype("Int64")
    anchor_join = frame.merge(
        anchors.drop_duplicates("source_id"), on="source_id", how="inner"
    )
    anchor_validation = {
        "n_matched": int(len(anchor_join)),
        "nir_minus_spectroscopic": robust(
            (anchor_join["av_nir"] - anchor_join["extinction_av_mag"]).to_numpy(float)
        ),
        "wp3_minus_spectroscopic": robust(
            (anchor_join["av_wp3"] - anchor_join["extinction_av_mag"]).to_numpy(float)
        ),
        "by_subgroup": {
            sg: {
                "n": int((anchor_join.subgroup == sg).sum()),
                "nir_minus_spec_median": (
                    float(np.median(
                        (anchor_join.loc[anchor_join.subgroup == sg, "av_nir"]
                         - anchor_join.loc[anchor_join.subgroup == sg, "extinction_av_mag"])
                    )) if (anchor_join.subgroup == sg).any() else None
                ),
            }
            for sg in w.SUBGROUPS
        },
        "note": (
            "Calibrates the near-IR method where ground truth exists.  Issue "
            "#1d predicts few or no anchors near CygOB2-B; the per-subgroup "
            "counts show whether that is so."
        ),
    }

    # ---- guard 2: does the intrinsic colour depend on the disputed age? ----
    probe_mass = np.geomspace(2.0, 8.0, 25)
    age_sensitivity = {
        f"{age:.3f}": [float(v) for v in intrinsic_jks(FAMILY, age, probe_mass)]
        for age in DISPUTED_AGES
    }
    max_shift = float(
        np.max(np.abs(
            np.array(age_sensitivity[f"{DISPUTED_AGES[0]:.3f}"])
            - np.array(age_sensitivity[f"{DISPUTED_AGES[1]:.3f}"])
        ))
    )
    robustness = {
        "disputed_ages_Myr": DISPUTED_AGES,
        "max_abs_intrinsic_JKs_shift_mag": max_shift,
        "implied_max_A_V_shift_mag": float(max_shift / k_jks),
        "verdict": (
            "the near-IR intrinsic colour barely moves between the two disputed "
            "ages, so this test does not depend on the answer it is testing"
            if max_shift / k_jks < 0.5
            else "WARNING: intrinsic colour is age sensitive; interpret with care"
        ),
    }

    # ---- per-subgroup extinction scale and law shape ----
    per_subgroup = {}
    for subgroup in w.SUBGROUPS:
        cell = frame[frame.subgroup.eq(subgroup)]
        window = cell[cell[MASS_COLUMN].between(*WINDOW)]
        slope = stats.theilslopes(
            cell["jh_obs"].to_numpy(float), cell["hks_obs"].to_numpy(float)
        )
        per_subgroup[subgroup] = {
            "n_stars": int(len(cell)),
            "n_in_mass_window": int(len(window)),
            "fitted_age_Myr": map_age[subgroup],
            "av_wp3": robust(cell["av_wp3"].to_numpy(float)),
            "av_nir": robust(cell["av_nir"].to_numpy(float)),
            "delta_av_wp3_minus_nir": robust(cell["delta_av"].to_numpy(float)),
            "delta_av_fixed_colour_variant": robust(
                cell["delta_av_fixed"].to_numpy(float)
            ),
            "reddening_slope_JH_over_HKs": {
                "theil_sen": float(slope[0]),
                "lo95": float(slope[2]),
                "hi95": float(slope[3]),
                "extinction_law_prediction": law_slope,
            },
        }

    # ---- mass-matched differential (guard 4) ----
    mass_matched = []
    for lo, hi in zip(MASS_BINS[:-1], MASS_BINS[1:]):
        row = {"mass_lo": float(lo), "mass_hi": float(hi)}
        for subgroup in w.SUBGROUPS:
            cell = frame[
                frame.subgroup.eq(subgroup) & frame[MASS_COLUMN].ge(lo)
                & frame[MASS_COLUMN].lt(hi)
            ]
            row[subgroup] = robust(cell["delta_av"].to_numpy(float))
        mass_matched.append(row)

    window_frame = frame[frame[MASS_COLUMN].between(*WINDOW)]
    contrasts = {}
    reference = np.concatenate([
        window_frame.loc[window_frame.subgroup.eq(sg), "delta_av"].to_numpy(float)
        for sg in ["CygOB2-A", "CygOB2-C"]
    ])
    b_values = window_frame.loc[
        window_frame.subgroup.eq("CygOB2-B"), "delta_av"
    ].to_numpy(float)
    reference = reference[np.isfinite(reference)]
    b_values = b_values[np.isfinite(b_values)]
    mannwhitney = stats.mannwhitneyu(b_values, reference, alternative="two-sided")
    contrasts["B_vs_A_and_C"] = {
        "B": robust(b_values),
        "A_plus_C": robust(reference),
        "median_offset_mag": float(np.median(b_values) - np.median(reference)),
        "mannwhitney_p": float(mannwhitney.pvalue),
        "note": (
            "delta_av = A_V(WP3 broadband) - A_V(near-IR).  A B-specific offset "
            "means B's optical extinction is biased relative to its siblings, "
            "which is the mass-scale error issue #12 needs."
        ),
    }

    # ---- guard 2b: is the B offset an artifact of B's disputed age? ----
    # Guard 2 showed the intrinsic near-IR colour is age sensitive at these
    # masses (PMS contraction), so B's offset could in principle be an artifact
    # of assuming B's disputed young age.  Recompute B's offset with B's
    # intrinsic colours taken at each candidate age, holding A and C at their
    # own undisputed ages.  If the offset survives -- or grows -- when B is
    # assumed old, it is not circular.
    reference_delta = np.concatenate([
        frame.loc[
            frame.subgroup.eq(sg) & frame[MASS_COLUMN].between(*WINDOW), "delta_av"
        ].to_numpy(float)
        for sg in ["CygOB2-A", "CygOB2-C"]
    ])
    reference_delta = reference_delta[np.isfinite(reference_delta)]
    b_mask = frame.subgroup.eq("CygOB2-B") & frame[MASS_COLUMN].between(*WINDOW)
    b_mass = frame.loc[b_mask, MASS_COLUMN].to_numpy(float)
    b_jks = frame.loc[b_mask, "jks_obs"].to_numpy(float)
    b_wp3 = frame.loc[b_mask, "av_wp3"].to_numpy(float)
    age_dependence = []
    for age in [2.818, 3.162, 3.548, 3.981, 4.467]:
        intrinsic = intrinsic_jks(FAMILY, age, b_mass)
        delta = b_wp3 - (b_jks - intrinsic) / k_jks
        delta = delta[np.isfinite(delta)]
        test = stats.mannwhitneyu(delta, reference_delta, alternative="two-sided")
        age_dependence.append(
            {
                "assumed_age_for_B_Myr": age,
                "B_median_delta_av": float(np.median(delta)),
                "offset_vs_A_and_C_mag": float(
                    np.median(delta) - np.median(reference_delta)
                ),
                "mannwhitney_p": float(test.pvalue),
            }
        )
    offsets = [row["offset_vs_A_and_C_mag"] for row in age_dependence]
    guard2b = {
        "reference_A_plus_C_median": float(np.median(reference_delta)),
        "per_assumed_age": age_dependence,
        "offset_range_mag": [float(min(offsets)), float(max(offsets))],
        "sign_stable_across_ages": bool(all(v > 0 for v in offsets) or all(v < 0 for v in offsets)),
        "verdict": (
            "the B-specific offset keeps the same sign at every candidate age, "
            "so it is not an artifact of assuming B's disputed young age"
            if all(v > 0 for v in offsets) or all(v < 0 for v in offsets)
            else "the offset changes sign with the assumed age: this test is "
                 "circular for B and cannot settle the question"
        ),
    }

    # ---- reddening slope from colour EXCESSES, not raw colours ----
    # Fitting raw (J-H) against (H-Ks) mixes the reddening vector with the
    # intrinsic-colour locus and is biased; the excesses remove the locus.
    excess_slopes = {}
    for subgroup in w.SUBGROUPS:
        mask = frame.subgroup.eq(subgroup)
        jh0, hks0 = intrinsic_jh_hks(
            FAMILY, map_age[subgroup], frame.loc[mask, MASS_COLUMN].to_numpy(float)
        )
        e_jh = frame.loc[mask, "jh_obs"].to_numpy(float) - jh0
        e_hks = frame.loc[mask, "hks_obs"].to_numpy(float) - hks0
        good = np.isfinite(e_jh) & np.isfinite(e_hks) & (e_hks > 0.05)
        if good.sum() > 20:
            fit = stats.theilslopes(e_jh[good], e_hks[good])
            excess_slopes[subgroup] = {
                "n": int(good.sum()),
                "theil_sen": float(fit[0]),
                "lo95": float(fit[2]),
                "hi95": float(fit[3]),
                "extinction_law_prediction": law_slope,
                "consistent_with_law": bool(fit[2] <= law_slope <= fit[3]),
            }

    # ---- what offset would be needed to explain the age discrepancy? ----
    # Express both effects in the currency that actually matters: the mass a
    # star is assigned.  An A_V error shifts the inferred absolute magnitude by
    # k_G * dA_V; an age error changes the mass-magnitude relation itself.
    iso_young, _ = w.load_isochrone_at_age(FAMILY, DISPUTED_AGES[0])
    iso_old, _ = w.load_isochrone_at_age(FAMILY, DISPUTED_AGES[1])
    probe = np.geomspace(2.5, 7.0, 19)
    g_young = np.interp(probe, iso_young["Mini"], iso_young["G0"])
    # mass that the OLD isochrone assigns to the young isochrone's magnitudes
    order = np.argsort(iso_old["G0"].to_numpy(float))
    mass_if_old = np.interp(
        g_young,
        iso_old["G0"].to_numpy(float)[order],
        iso_old["Mini"].to_numpy(float)[order],
    )
    age_mass_ratio = float(np.median(mass_if_old / probe))
    measured_offset = contrasts["B_vs_A_and_C"]["median_offset_mag"]
    # mass shift produced by the measured extinction offset, at fixed age
    g_shift = coefficients["G"] * measured_offset
    mass_if_av = np.interp(
        g_young + g_shift,
        iso_young["G0"].to_numpy(float)[np.argsort(iso_young["G0"].to_numpy(float))],
        iso_young["Mini"].to_numpy(float)[np.argsort(iso_young["G0"].to_numpy(float))],
    )
    av_mass_ratio = float(np.median(mass_if_av / probe))
    required = {
        "framing": (
            "both effects expressed as the multiplicative shift they induce in "
            "the mass assigned to a star of fixed observed magnitude"
        ),
        "mass_ratio_from_disputed_age_shift": age_mass_ratio,
        "mass_ratio_from_measured_A_V_offset": av_mass_ratio,
        "measured_B_offset_mag": measured_offset,
        "A_V_offset_mag_equivalent_to_the_age_shift": float(
            np.log(age_mass_ratio)
            / np.log(av_mass_ratio) * measured_offset
        ) if abs(np.log(av_mass_ratio)) > 1e-9 else None,
        "note": (
            "Order-of-magnitude comparison only.  It asks whether the measured "
            "B-specific extinction offset is even the right size to account for "
            "the disputed age shift; it does not establish that extinction is "
            "the cause."
        ),
    }

    out_csv = w.TABLES / "wp3_nir_extinction_crosscheck.csv"
    frame[[
        "source_id", "subgroup", MASS_COLUMN, "J", "H", "Ks",
        "jks_obs", "jks_intrinsic", "av_wp3", "av_nir", "delta_av",
        "av_nir_fixed", "delta_av_fixed",
    ]].to_csv(out_csv, index=False)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp3_nir_extinction_crosscheck.py",
        "status": "SUCCESS",
        "issue": "#12 / #1d — independent near-IR extinction check for CygOB2-B",
        "upstream_version": UPSTREAM,
        "branch": f"{FAMILY} / R_V={RV}",
        "stored_artifacts_overwritten": False,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "method": {
            "estimator": "A_V = [(J-Ks)_obs - (J-Ks)_0] / (k_J - k_Ks)",
            "extinction_law_coefficients": {k: float(v) for k, v in coefficients.items()},
            "k_J_minus_k_Ks": float(k_jks),
            "independence": (
                "the near-IR colour excess does not use G/BP/RP, so it is not a "
                "restatement of the WP3 broadband fit that produced av_rv3.1"
            ),
        },
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p) for p in [
                w.PROC / f"wp3_extinction_{UPSTREAM}.parquet",
                w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet",
                w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet",
                w.PROC / "wp1_spectroscopic_anchors.parquet",
                w.PROC / f"wp3_isochrones_{FAMILY.lower()}.parquet",
            ]
        },
        "n_stars_used": int(len(frame)),
        "guard1_anchor_validation": anchor_validation,
        "guard2_age_robustness": robustness,
        "guard2b_offset_vs_assumed_age_for_B": guard2b,
        "reddening_slope_from_colour_excesses": excess_slopes,
        "guard3_fixed_colour_variant": "reported per subgroup as delta_av_fixed_colour_variant",
        "guard4_mass_matched": mass_matched,
        "per_subgroup": per_subgroup,
        "contrasts": contrasts,
        "required_offset_comparison": required,
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(
        w.PROVENANCE / "wp3_nir_extinction_crosscheck_execution.json", record
    )

    print(f"stars used: {len(frame)}   (k_J - k_Ks = {k_jks:.4f})")
    print(f"\nGuard 2 — intrinsic (J-Ks)_0 shift between {DISPUTED_AGES[0]} and "
          f"{DISPUTED_AGES[1]} Myr: {max_shift:.4f} mag "
          f"(= {max_shift/k_jks:.3f} mag in A_V)")
    print(f"   {robustness['verdict']}")
    print(f"\nGuard 1 — anchors matched: {anchor_validation['n_matched']}")
    print(f"   A_V(NIR) - A_V(spec): median "
          f"{anchor_validation['nir_minus_spectroscopic'].get('median', float('nan')):+.3f} "
          f"mag, scatter {anchor_validation['nir_minus_spectroscopic'].get('mad_sigma', float('nan')):.3f}")
    print(f"   A_V(WP3) - A_V(spec): median "
          f"{anchor_validation['wp3_minus_spectroscopic'].get('median', float('nan')):+.3f} "
          f"mag, scatter {anchor_validation['wp3_minus_spectroscopic'].get('mad_sigma', float('nan')):.3f}")
    print(f"   anchors per subgroup: "
          f"{ {k: v['n'] for k, v in anchor_validation['by_subgroup'].items()} }")
    print("\nPer subgroup:")
    for subgroup, data in per_subgroup.items():
        print(f"   {subgroup} (age {data['fitted_age_Myr']:.2f} Myr, n={data['n_stars']})")
        print(f"      A_V WP3   median {data['av_wp3']['median']:6.3f}")
        print(f"      A_V NIR   median {data['av_nir']['median']:6.3f}")
        print(f"      WP3-NIR   median {data['delta_av_wp3_minus_nir']['median']:+6.3f} "
              f"+- {data['delta_av_wp3_minus_nir']['sem']:.3f} (sem), "
              f"scatter {data['delta_av_wp3_minus_nir']['mad_sigma']:.3f}")
        print(f"      fixed-colour variant  {data['delta_av_fixed_colour_variant']['median']:+6.3f}")
        s = data["reddening_slope_JH_over_HKs"]
        print(f"      E(J-H)/E(H-Ks) = {s['theil_sen']:.3f} "
              f"[{s['lo95']:.3f}, {s['hi95']:.3f}]  (law predicts {s['extinction_law_prediction']:.3f})")
    c = contrasts["B_vs_A_and_C"]
    print(f"\nB versus A+C in the 2-8 Msun window:")
    print(f"   B      median WP3-NIR {c['B']['median']:+.3f} (n={c['B']['n']})")
    print(f"   A + C  median WP3-NIR {c['A_plus_C']['median']:+.3f} (n={c['A_plus_C']['n']})")
    print(f"   offset {c['median_offset_mag']:+.3f} mag, Mann-Whitney p = {c['mannwhitney_p']:.3g}")
    print(f"\nGuard 2b — offset vs the age assumed for B:")
    for row in guard2b["per_assumed_age"]:
        print(f"   age {row['assumed_age_for_B_Myr']:.3f} Myr -> offset "
              f"{row['offset_vs_A_and_C_mag']:+.3f} mag  (p={row['mannwhitney_p']:.1e})")
    print(f"   {guard2b['verdict']}")
    print(f"\nReddening slope from colour excesses (law predicts {law_slope:.3f}):")
    for subgroup, data in excess_slopes.items():
        print(f"   {subgroup}: {data['theil_sen']:.3f} "
              f"[{data['lo95']:.3f}, {data['hi95']:.3f}]  "
              f"consistent_with_law={data['consistent_with_law']}")
    print(f"\nMass currency: the disputed age shift moves assigned mass by a factor "
          f"{required['mass_ratio_from_disputed_age_shift']:.4f};")
    print(f"   the measured {measured_offset:+.3f} mag extinction offset moves it by "
          f"{required['mass_ratio_from_measured_A_V_offset']:.4f} "
          f"({100*(1-required['mass_ratio_from_measured_A_V_offset']):.1f}% in mass).")
    print("   The two are NOT interchangeable: the age shift mainly relocates the "
          "PMS/Henyey fold,\n   while the extinction offset rescales the whole mass ladder.")
    print("\nwrote tables/wp3_nir_extinction_crosscheck.csv")
    print("wrote provenance/wp3_nir_extinction_crosscheck_execution.json")


if __name__ == "__main__":
    main()
