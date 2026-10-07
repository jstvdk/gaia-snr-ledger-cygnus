#!/usr/bin/env python3
"""Pre-registration of the repair_v9 sensitivity items S1 (IMF ceiling) and S2
(subgroup-assignment uncertainty), brief §4.5.

Written and committed while the repair_v9 chain was still running, i.e.
before any repair_v9 WP5 or WP7 result existed.  Owner delegation: decision O3
of provenance/decisions_2026_10_07_overnight.json (the owner pre-approved S1
and S2 and delegated their sign-off for that night).  Both run only if
repair_v9 is adopted (decision O4).

Run:  PYTHONPATH=scripts python3 scripts/repair_v9_s1_s2_prereg.py
"""
from __future__ import annotations

import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp6_mass_extension_decision import turnoff_mass

OUT = w.PROVENANCE / "repair_v9_s1_s2_prereg.json"
S1_TOL_REL = 0.30
S1_REPLAY_ABS = 0.10
S2_ASSOC_REL = 0.10
S2_K_REL = 0.15
S2_WEIGHT_CONSERVATION_REL = 1e-6
CODE = ["scripts/repair_v9_s1_imf_ceiling.py", "scripts/repair_v9_s2_label_uncertainty.py",
        "scripts/repair_v9_s1_s2_prereg.py", "scripts/wp7_ledger.py", "scripts/wp2_derive_subgroups.py",
        "scripts/wp5_fit_imf_joint.py", "scripts/wp5_joint_age_fit.py"]
INPUTS = ["data/processed/wp4_age_posteriors_repair_v9_headline.parquet",
          "data/processed/wp5_imf_normalization_repair_v8.parquet",
          "tables/wp2_subgroup_labels.parquet", "data/processed/wp2_members.parquet"]


def integral(lo: float, hi: float, a: float = 2.3) -> float:
    return max(0.0, (lo ** (1 - a) - hi ** (1 - a)) / (a - 1)) if lo < hi else 0.0


def s1_expectation() -> dict:
    """Coeval, MAP-age, analytic: N_sg ~ k_sg * integral[M_to(t_sg), m_max] M^-2.3,
    with k from repair_v8's baseline cell (the v9 k does not exist yet) and t
    the repair_v9 headline MAP.  The posterior widths smear the turnoff, which
    pulls the real ratios towards 1 relative to this point estimate."""
    head = pd.read_parquet(w.ROOT / INPUTS[0])
    n8 = pd.read_parquet(w.ROOT / INPUTS[1])
    out = {}
    for fam in w.FAMILIES:
        tot = {m: 0.0 for m in (100.0, 120.0, 150.0)}
        per = {}
        for sg in w.SUBGROUPS:
            t = float(head[head.subgroup.eq(sg) & head.family.eq(fam) & np.isclose(head.R_V, 3.1)
                           & np.isclose(head.f_bin, 0.4) & np.isclose(head.dmu, 0.0)].age_map.iloc[0])
            k = float(n8[n8.subgroup.eq(sg) & n8.family.eq(fam) & np.isclose(n8.R_V, 3.1)
                         & np.isclose(n8.alpha, 2.3)].k_median.iloc[0])
            mto = turnoff_mass(fam, t)
            per[sg] = {"age_Myr": t, "turnoff_Msun": round(mto, 1)}
            for m in tot:
                n = k * integral(min(mto, m), m)
                per[sg][f"N_{m:.0f}"] = round(n, 3)
                tot[m] += n
        out[fam] = {"subgroups": per, "association": {f"N_{m:.0f}": round(v, 3) for m, v in tot.items()},
                    "ratio_100_to_120": round(tot[100.0] / tot[120.0], 3),
                    "ratio_150_to_120": round(tot[150.0] / tot[120.0], 3)}
    return out


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; written once")
    for p in ("tables/repair_v9_s1_imf_ceiling.csv", "tables/repair_v9_s2_label_uncertainty.csv",
              "tables/repair_v9_s2_responsibilities.csv", "tables/wp7_ledger_repair_v9.csv"):
        if (w.ROOT / p).exists():
            raise SystemExit(f"{p} exists; S1/S2 can no longer be pre-registered blind")
    exp = s1_expectation()
    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_s1_s2_prereg.py",
        "python": sys.version.split()[0],
        "basis": "brief §4.5 (tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md); parent chain "
                 "provenance/repair_v9_prereg.json",
        "signoff": ("delegated by the owner for the night of 2026-10-07 (decision O3, "
                    "provenance/decisions_2026_10_07_overnight.json)"),
        "runs_only_if": "repair_v9 is adopted (I1-I5 pass), decision O4",
        "blindness": "written while the repair_v9 chain was running; no repair_v9 WP5/WP7 product existed",
        "S1": {
            "design": ("m_max in {100, 120, 150} Msun; WP7 engine unmodified (IMF_UPPER_LIMIT set in the "
                       "wp7_ledger namespace); the 18 alpha = 2.3 branches, each subgroup and the "
                       "association; all-explode; 2,000,000 iterations on the repair_v9 WP5 draws; the same "
                       "seed for every m_max (paired)"),
            "issue_14": ("150 Msun is an evolutionary turnoff in both grids (PARSEC 2.82 Myr; MIST between "
                         "2.00 and 2.25 Myr), below the table ceilings 300 / 210 Msun; TurnoffRelation "
                         "refuses to build if the grid does not bracket m_max"),
            "integrity": {"S1-I1": (f"at m_max = 120 the PARSEC baseline association mean reproduces "
                                    f"tables/wp7_ledger_repair_v9.csv within {S1_REPLAY_ABS} (different seeds)")},
            "analytic_expectation": exp,
            "analytic_note": s1_expectation.__doc__,
            "disclosure": ("the brief's 'analytic expectation on v8 inputs: 100 -> 7.0, 120 -> 8.4, 150 -> 9.7' "
                           "used the repair_v8 turnoffs (~58 Msun).  At the repair_v9 headline ages the turnoffs "
                           "are 74-96 Msun, close to the ceiling, so the ceiling matters much more; the "
                           "expectation above replaces the brief's numbers"),
            "predictions": {
                "S1-P1": "N_death(100) <= N_death(120) <= N_death(150) in every branch x subgroup row",
                "S1-P2": (f"PARSEC baseline association ratios N(100)/N(120) and N(150)/N(120) each within "
                          f"+-{S1_TOL_REL:.0%} (relative) of the analytic expectation above"),
            },
        },
        "S2": {
            "design": ("replay the frozen WP2 mixture (k = 3, full covariance, StandardScaler on l, b, pmra, "
                       "pmdec, the 50 frozen seeds, components named by wp2_derive_subgroups.name_components); "
                       "responsibilities = predict_proba averaged over seeds; every labelled star of the clean "
                       "WP2 set enters each subgroup with weight membership_probability x r_subgroup (the WP5 "
                       "fit's own counting rule); each star keeps its repair_v9 mass posterior; WP5 refitted on "
                       "the repair_v9 ages and responses; WP7 engine on both the repair_v9 (hard) and S2 (soft) "
                       "draws with identical seeds, 54 branches, 500,000 iterations"),
            "approximation": ("a relabelled star keeps the mass posterior computed at its hard subgroup's age; "
                              "the headline ages differ by <= 0.4 Myr, which moves masses near the turnoff, not "
                              "in WP5's calibration window below 8 Msun"),
            "integrity": {"S2-I1": (f"summed soft membership weight equals summed hard weight over the soft "
                                    f"stars within {S2_WEIGHT_CONSERVATION_REL} (responsibilities sum to 1)")},
            "disclosed_prior": ("the brief reports a median maximum responsibility of 0.90 and 51 % of stars "
                                "above 0.9, B clean, A and C mixing ~15 % each way"),
            "predictions": {
                "S2-P1": (f"PARSEC baseline association N_death (soft) within {S2_ASSOC_REL:.0%} of hard"),
                "S2-P2": (f"PARSEC baseline cell k_median (soft) within {S2_K_REL:.0%} of hard, each subgroup"),
                "S2-P3": "B's k_median changes least of the three subgroups (B is the cleanest)",
            },
        },
        "thresholds": {"S1_TOL_REL": S1_TOL_REL, "S1_REPLAY_ABS": S1_REPLAY_ABS,
                       "S2_ASSOC_REL": S2_ASSOC_REL, "S2_K_REL": S2_K_REL,
                       "S2_WEIGHT_CONSERVATION_REL": S2_WEIGHT_CONSERVATION_REL},
        "failed_predictions_rule": "scored as written; a failed prediction stays failed",
        "consumed_inputs": {p: w.sha256(w.ROOT / p) for p in INPUTS},
        "frozen_code": {p: w.sha256(w.ROOT / p) for p in CODE},
    })
    print(f"wrote {OUT.relative_to(w.ROOT)}")
    for fam, e in exp.items():
        print(fam, e["association"], "ratios", e["ratio_100_to_120"], e["ratio_150_to_120"])


if __name__ == "__main__":
    main()
