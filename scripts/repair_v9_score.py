#!/usr/bin/env python3
"""repair_v9 -- score predictions V1-V5 and apply the adoption rule
(provenance/repair_v9_prereg.json).  Predictions are scored as written and
never gate adoption; adoption rests on I1-I5 alone.

Every comparison value is READ from the repair_v8 artifacts, never typed.

Outputs:
  provenance/repair_v9_outcome.json
  tables/repair_v9_before_after.csv

Run:  PYTHONPATH=scripts python3 scripts/repair_v9_score.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

PREREG = w.PROVENANCE / "repair_v9_prereg.json"
INTEGRITY = w.PROVENANCE / "repair_v9_integrity.json"
OUT = w.PROVENANCE / "repair_v9_outcome.json"
OUT_TABLE = w.TABLES / "repair_v9_before_after.csv"
V1B_RANGE = (4.0, 9.0)
V3_MIN = 0.5
V4_DROP_ABS, V4_DROP_REL, V4_ZERO_MAX, V4_YOUNG_MAX_MYR = 0.02, 0.02, 0.1, 2.53
V5_RATIO_MAX, V5_FROM_MYR = 1.25, 3.15
INTEGRITY_CHECKS = ("I1a", "I1b", "I1c", "I1d", "I1e", "I2", "I3", "I4", "I5")


def ledger(chain: str) -> pd.DataFrame:
    return pd.read_csv(w.TABLES / f"wp7_ledger_{chain}.csv")


def cell(l: pd.DataFrame, fam="PARSEC", rv=3.1, alpha=2.3, delta=0.0, expl="all_explode"):
    return l[l.family.eq(fam) & np.isclose(l.R_V, rv) & np.isclose(l.alpha, alpha)
             & np.isclose(l.sf_duration_Myr, delta) & l.explodability.eq(expl)]


def by_subgroup(sel: pd.DataFrame, col="N_SN_mean") -> dict:
    return {r.subgroup: float(getattr(r, col)) for r in sel.itertuples()}


def closure(chain: str, fam="PARSEC", rv=3.1, alpha=2.3) -> dict:
    c = pd.read_csv(w.TABLES / f"wp6_closure_{chain}.csv")
    c = c[c.family.eq(fam) & np.isclose(c.R_V, rv) & np.isclose(c.alpha, alpha)]
    return {r.subgroup: float(r.closure_ratio) for r in c.itertuples()}


def main() -> None:
    v8, v9 = ledger("repair_v8"), ledger("repair_v9")
    b8, b9 = by_subgroup(cell(v8)), by_subgroup(cell(v9))
    A, B, Cc, ALL = "CygOB2-A", "CygOB2-B", "CygOB2-C", "ALL"

    v1a = b9[ALL] < b8[ALL]
    v1b = V1B_RANGE[0] <= b9[ALL] <= V1B_RANGE[1]
    v1c = b9[A] < b8[A] and b9[B] < b8[B] and b9[Cc] > b8[Cc]

    c8, c9 = closure("repair_v8"), closure("repair_v9")
    c8m, c9m = closure("repair_v8", "MIST"), closure("repair_v9", "MIST")
    want = {A: "down", B: "down", Cc: "up"}
    moved = {sg: ("up" if c9[sg] > c8[sg] else "down") for sg in want}
    v2 = all(moved[sg] == want[sg] for sg in want)
    moved_mist = {sg: ("up" if c9m[sg] > c8m[sg] else "down") for sg in want}

    v3 = b9[Cc] > V3_MIN

    scan = pd.read_csv(w.TABLES / "repair_v9_age_scan.csv")
    viol, young = [], []
    for keys, s in scan.groupby(["family", "R_V", "alpha", "sf_duration_Myr", "subgroup"]):
        s = s.sort_values("age_Myr")
        m = s.N_death_mean.to_numpy()
        for j in range(1, len(m)):
            if m[j] < m[j - 1] - max(V4_DROP_ABS, V4_DROP_REL * m[j - 1]):
                viol.append({"cell": list(keys), "age_from": float(s.age_Myr.iloc[j - 1]),
                             "age_to": float(s.age_Myr.iloc[j]), "from": float(m[j - 1]), "to": float(m[j])})
        if np.isclose(keys[3], 0.0):
            y = s[s.age_Myr <= V4_YOUNG_MAX_MYR]
            for r in y.itertuples():
                if r.N_death_mean >= V4_ZERO_MAX:
                    young.append({"cell": list(keys), "age": r.age_Myr, "N_death_mean": r.N_death_mean})
    v4 = not viol and not young

    v5_bad = []
    for (fam, rv, idx), s in scan[np.isclose(scan.alpha, 2.3) & np.isclose(scan.sf_duration_Myr, 0.0)
                                  ].groupby(["family", "R_V", "scan_index"]):
        if s.age_Myr.iloc[0] < V5_FROM_MYR:
            continue
        lo, hi = s.N_death_mean.min(), s.N_death_mean.max()
        ratio = hi / lo if lo > 0 else np.inf
        if ratio > V5_RATIO_MAX:
            v5_bad.append({"family": fam, "R_V": rv, "age": float(s.age_Myr.iloc[0]),
                           "ratio_max_min": float(ratio),
                           "means": by_subgroup(s, "N_death_mean")})
    v5 = not v5_bad

    integ = json.loads(INTEGRITY.read_text())["checks"]
    missing = [c for c in INTEGRITY_CHECKS if c not in integ]
    passed = {c: bool(integ[c]["pass"]) for c in INTEGRITY_CHECKS if c in integ}
    adopt = not missing and all(passed.values())

    # before / after for every headline quantity of the ledger
    rows = []
    for name, sel8, sel9, col in [
        ("N_death baseline (PARSEC R_V 3.1 alpha 2.3 coeval)", cell(v8), cell(v9), "N_SN_mean"),
        ("P(>=1 death) baseline", cell(v8), cell(v9), "P_at_least_one"),
        ("P(last < 100 kyr) baseline", cell(v8), cell(v9), "P_last_SN_within_100kyr"),
    ]:
        x8, x9 = by_subgroup(sel8, col), by_subgroup(sel9, col)
        for sg in (A, B, Cc, ALL):
            rows.append({"quantity": name, "subgroup": sg, "repair_v8": x8[sg], "repair_v9": x9[sg]})
    head = lambda l: l[l.scope.eq("association") & np.isclose(l.alpha, 2.3)  # noqa: E731
                       & l.explodability.eq("all_explode")].N_SN_mean
    all54 = lambda l: l[l.scope.eq("association") & l.explodability.eq("all_explode")].N_SN_mean  # noqa: E731
    for name, f in [("headline range min (18 alpha 2.3 branches)", lambda l: head(l).min()),
                    ("headline range max (18 alpha 2.3 branches)", lambda l: head(l).max()),
                    ("full grid min (54 branches)", lambda l: all54(l).min()),
                    ("full grid max (54 branches)", lambda l: all54(l).max()),
                    ("full grid median (54 branches)", lambda l: all54(l).median())]:
        rows.append({"quantity": name, "subgroup": ALL, "repair_v8": float(f(v8)), "repair_v9": float(f(v9))})
    for sg in (A, B, Cc):
        rows.append({"quantity": "closure ratio baseline cell", "subgroup": sg,
                     "repair_v8": c8[sg], "repair_v9": c9[sg]})
    n8 = pd.read_parquet(w.PROC / "wp5_imf_normalization_repair_v8.parquet")
    n9 = pd.read_parquet(w.PROC / "wp5_imf_normalization_repair_v9.parquet")
    for sg in (A, B, Cc):
        k ={lab: float(n[n.subgroup.eq(sg) & n.family.eq("PARSEC") & np.isclose(n.R_V, 3.1)
                         & np.isclose(n.alpha, 2.3)].k_median.iloc[0]) for n, lab in ((n8, "8"), (n9, "9"))}
        t = {lab: float(n[n.subgroup.eq(sg) & n.family.eq("PARSEC") & np.isclose(n.R_V, 3.1)
                         & np.isclose(n.alpha, 2.3)].truth_age_posterior_mean_Myr.iloc[0])
             for n, lab in ((n8, "8"), (n9, "9"))}
        rows.append({"quantity": "WP5 k_median baseline cell", "subgroup": sg,
                     "repair_v8": k["8"], "repair_v9": k["9"]})
        rows.append({"quantity": "WP5 truth-age posterior mean (Myr)", "subgroup": sg,
                     "repair_v8": t["8"], "repair_v9": t["9"]})
    gate = lambda n: int(n.residual_gate_pass.astype(bool).sum())  # noqa: E731
    rows.append({"quantity": "WP5 residual-gate cells passing (of 54)", "subgroup": ALL,
                 "repair_v8": gate(n8), "repair_v9": gate(n9)})
    table = pd.DataFrame(rows)
    table["change"] = table.repair_v9 - table.repair_v8
    table.to_csv(OUT_TABLE, index=False)

    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_score.py",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "predictions": {
            "V1a": {"statement": "PARSEC baseline N_death < repair_v8", "repair_v8": b8[ALL],
                    "repair_v9": b9[ALL], "true": bool(v1a)},
            "V1b": {"statement": f"PARSEC baseline N_death in {list(V1B_RANGE)}", "repair_v9": b9[ALL],
                    "true": bool(v1b)},
            "V1c": {"statement": "A and B fall, C rises (baseline subgroup means)",
                    "repair_v8": {k: b8[k] for k in (A, B, Cc)},
                    "repair_v9": {k: b9[k] for k in (A, B, Cc)}, "true": bool(v1c)},
            "V2": {"statement": "baseline-cell closure ratio: A down, B down, C up (PARSEC); MIST reported",
                   "repair_v8": c8, "repair_v9": c9, "moved": moved, "true": bool(v2),
                   "report_only_MIST": {"repair_v8": c8m, "repair_v9": c9m, "moved": moved_mist}},
            "V3": {"statement": f"C baseline N_death > {V3_MIN}", "repair_v9": b9[Cc], "true": bool(v3)},
            "V4": {"statement": ("scan N_death(t) non-decreasing in every subgroup x branch (a drop beyond "
                                 f"max({V4_DROP_ABS}, {V4_DROP_REL:.0%}) is a violation), and < {V4_ZERO_MAX} "
                                 f"for t <= {V4_YOUNG_MAX_MYR} Myr on coeval branches"),
                   "monotonicity_violations": viol, "young_nonzero": young, "true": bool(v4)},
            "V5": {"statement": (f"subgroup curves agree within max/min <= {V5_RATIO_MAX} at every "
                                 f"t >= {V5_FROM_MYR} Myr (alpha 2.3, coeval, every family x R_V)"),
                   "violations": v5_bad, "true": bool(v5)},
        },
        "integrity": {"passed": passed, "missing": missing},
        "adoption": {"rule": "ADOPTED = repair_v9 iff I1-I5 all pass; predictions do not gate",
                     "adopt": bool(adopt)},
        "outputs": {str(OUT_TABLE.relative_to(w.ROOT)): w.sha256(OUT_TABLE)},
    })
    print(table.to_string(index=False))
    print("V1a", v1a, "V1b", v1b, "V1c", v1c, "V2", v2, moved, "V3", v3, "V4", v4, "V5", v5)
    print("integrity", passed, "missing", missing, "-> adopt", adopt)


if __name__ == "__main__":
    main()
