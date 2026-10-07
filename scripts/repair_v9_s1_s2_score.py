#!/usr/bin/env python3
"""Score S1 and S2 exactly as provenance/repair_v9_s1_s2_prereg.json states
(thresholds read from that record, never typed here).

Output:  provenance/repair_v9_s1_s2_outcome.json
Run:     PYTHONPATH=scripts python3 scripts/repair_v9_s1_s2_score.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

PREREG = w.PROVENANCE / "repair_v9_s1_s2_prereg.json"
OUT = w.PROVENANCE / "repair_v9_s1_s2_outcome.json"


def base(t: pd.DataFrame, fam="PARSEC") -> pd.DataFrame:
    return t[t.family.eq(fam) & np.isclose(t.R_V, 3.1) & np.isclose(t.sf_duration_Myr, 0.0)
             & (np.isclose(t.alpha, 2.3) if "alpha" in t else True)]


def main() -> None:
    pre = json.loads(PREREG.read_text())
    th = pre["thresholds"]
    out = {"created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "script": "scripts/repair_v9_s1_s2_score.py",
           "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)}}

    s1p = w.TABLES / "repair_v9_s1_imf_ceiling.csv"
    if s1p.exists():
        s1 = pd.read_csv(s1p)
        led = pd.read_csv(w.TABLES / "wp7_ledger_repair_v9.csv")
        ref = float(led[led.scope.eq("association") & led.family.eq("PARSEC") & np.isclose(led.R_V, 3.1)
                        & np.isclose(led.alpha, 2.3) & np.isclose(led.sf_duration_Myr, 0.0)
                        & led.explodability.eq("all_explode")].N_SN_mean.iloc[0])
        b = base(s1).set_index(["subgroup", "m_max_Msun"]).N_death_mean
        i1 = abs(b[("ALL", 120.0)] - ref) <= th["S1_REPLAY_ABS"]
        piv = s1.pivot_table(index=["family", "R_V", "sf_duration_Myr", "subgroup"],
                             columns="m_max_Msun", values="N_death_mean")
        p1 = bool(((piv[100.0] <= piv[120.0]) & (piv[120.0] <= piv[150.0])).all())
        exp = pre["S1"]["analytic_expectation"]["PARSEC"]
        r100, r150 = b[("ALL", 100.0)] / b[("ALL", 120.0)], b[("ALL", 150.0)] / b[("ALL", 120.0)]
        ok = lambda got, want: abs(got / want - 1.0) <= th["S1_TOL_REL"]  # noqa: E731
        p2 = ok(r100, exp["ratio_100_to_120"]) and ok(r150, exp["ratio_150_to_120"])
        out["S1"] = {"S1-I1": {"pass": bool(i1), "m120_baseline": float(b[("ALL", 120.0)]),
                               "repair_v9_ledger": ref},
                     "S1-P1": {"true": p1},
                     "S1-P2": {"true": bool(p2), "ratio_100_to_120": float(r100),
                               "ratio_150_to_120": float(r150), "expected": [exp["ratio_100_to_120"],
                                                                             exp["ratio_150_to_120"]]},
                     "baseline_PARSEC": {f"{sg}|{m:.0f}": float(v) for (sg, m), v in b.items()},
                     "baseline_MIST": {f"{sg}|{m:.0f}": float(v) for (sg, m), v in
                                       base(s1, "MIST").set_index(["subgroup", "m_max_Msun"]).N_death_mean.items()}}

    s2p = w.TABLES / "repair_v9_s2_label_uncertainty.csv"
    if s2p.exists():
        s2 = pd.read_csv(s2p)
        b = base(s2).set_index(["subgroup", "labels"])
        n_h, n_s = b.N_death_mean[("ALL", "hard")], b.N_death_mean[("ALL", "soft")]
        k = {sg: (b.k_median[(sg, "hard")], b.k_median[(sg, "soft")]) for sg in w.SUBGROUPS}
        rel = {sg: abs(s / h - 1.0) for sg, (h, s) in k.items()}
        hard_m = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v9.parquet")
        soft_m = pd.read_parquet(w.PROC / "wp4_mass_posteriors_repair_v9_s2soft.parquet")
        resp = pd.read_csv(w.TABLES / "repair_v9_s2_responsibilities.csv")
        ids = set(resp.source_id)
        hw = hard_m[hard_m.source_id.isin(ids) & hard_m.subgroup.isin(w.SUBGROUPS)].membership_probability.sum()
        sw = soft_m[soft_m.source_id.isin(ids) & soft_m.subgroup.isin(w.SUBGROUPS)].membership_probability.sum()
        out["S2"] = {"S2-I1": {"pass": bool(abs(sw / hw - 1) <= th["S2_WEIGHT_CONSERVATION_REL"]),
                               "hard_weight": float(hw), "soft_weight": float(sw)},
                     "S2-P1": {"true": bool(abs(n_s / n_h - 1) <= th["S2_ASSOC_REL"]),
                               "N_hard": float(n_h), "N_soft": float(n_s)},
                     "S2-P2": {"true": bool(max(rel.values()) <= th["S2_K_REL"]),
                               "k_hard_soft": {sg: [float(h), float(s)] for sg, (h, s) in k.items()},
                               "k_rel_change": {sg: float(v) for sg, v in rel.items()}},
                     "S2-P3": {"true": bool(min(rel, key=rel.get) == "CygOB2-B")},
                     "subgroup_N_hard_soft": {sg: [float(b.N_death_mean[(sg, "hard")]),
                                                   float(b.N_death_mean[(sg, "soft")])]
                                              for sg in w.SUBGROUPS},
                     "responsibilities": {"median_r_max": float(resp.r_max.median()),
                                          "fraction_r_max_gt_0p9": float((resp.r_max > 0.9).mean())}}
    w.write_json(OUT, out)
    print(json.dumps({k: v for k, v in out.items() if k in ("S1", "S2")}, indent=1)[:3000])


if __name__ == "__main__":
    main()
