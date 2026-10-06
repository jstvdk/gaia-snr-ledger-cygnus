#!/usr/bin/env python3
"""Issue #20 Phase A' -- apply the frozen adoption rule (provenance/issue20b_prereg.json).

Outputs:
  tables/issue20b_adopted_ages.csv
  provenance/issue20b_outcome.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20b_score.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20b_common as B
import wp5_common as w

OUT_TABLE = w.TABLES / "issue20b_adopted_ages.csv"
OUT = w.PROVENANCE / "issue20b_outcome.json"
EPS = B.EPS_OUTLIER
KEYS = ("map", "lo68", "hi68", "lo95", "hi95")


def main() -> None:
    B.check_method_unchanged()
    run = json.loads((w.PROVENANCE / "issue20b_run_execution.json").read_text())
    tests = pd.read_csv(w.TABLES / "issue20b_age_tests.csv")
    curves = run["loglike_curves"]
    ages = {f: np.array(v) for f, v in run["native_ages"].items()}

    def result(test, sg, fam, rv):
        sel = tests[tests.test.eq(test) & tests.subgroup.eq(sg) & tests.family.eq(fam) & tests.R_V.eq(rv)]
        if test != "a":
            sel = sel[np.isclose(sel.eps, EPS)]
        r = sel.iloc[0]
        out = {"valid": bool(r.valid) if pd.notna(r.valid) else False, "n_stars": int(r.n_stars)}
        if out["n_stars"] > 0 and pd.notna(r.get("map")):
            out.update({k: float(r[k]) for k in KEYS})
            out["railed"] = bool(r.railed)
        return out

    rows, cells = [], {}
    for rv in B.R_V_BRANCHES:
        for fam in B.FAMILIES:
            for sg in B.SUBGROUPS:
                res = {t: result(t, sg, fam, rv) for t in ("a", "b", "c")}
                dec = B.adopt(res)
                rec = {"subgroup": sg, "family": fam, "R_V": rv, "status": dec["status"],
                       "agreeing_tests": ",".join(dec.get("tests", [])),
                       "valid_tests": ",".join(dec["valid_tests"]), "reason": dec.get("reason", "")}
                for t, r in res.items():
                    for k in ("map", "lo68", "hi68"):
                        rec[f"{t}_{k}"] = r.get(k, np.nan)
                    rec[f"{t}_valid"] = r["valid"]
                    rec[f"{t}_n"] = r["n_stars"]
                if dec["status"] == "ADOPTED":
                    lls = []
                    for t in dec["tests"]:
                        key = f"{t}|{'nan' if t == 'a' else EPS}|{sg}|{fam}|{rv}"
                        lls.append(np.array(curves[key]))
                    rec.update({f"adopted_{k}": v for k, v in
                                B.mixture_posterior(ages[fam], lls).items()})
                rows.append(rec)
                cells[(sg, fam, rv)] = rec
    table = pd.DataFrame(rows)
    table.to_csv(OUT_TABLE, index=False)

    # C versus A at R_V 3.1
    per_family = {}
    for fam in B.FAMILIES:
        a, c = cells[("CygOB2-A", fam, 3.1)], cells[("CygOB2-C", fam, 3.1)]
        if a["status"] != "ADOPTED" or c["status"] != "ADOPTED":
            per_family[fam] = "undetermined (A or C not adopted)"
        elif c["adopted_hi95"] < a["adopted_lo95"]:
            per_family[fam] = "C younger"
        elif not (c["adopted_hi68"] < a["adopted_lo68"] or a["adopted_hi68"] < c["adopted_lo68"]):
            per_family[fam] = "C coeval"
        else:
            per_family[fam] = "undetermined (intervals neither overlap at 68 % nor separate at 95 %)"
    vals = set(per_family.values())
    c_vs_a = ("C younger than A" if vals == {"C younger"} else
              "C coeval with A" if vals == {"C coeval"} else "C's age relative to A undetermined")

    robustness = {}
    for sg in B.SUBGROUPS:
        for fam in B.FAMILIES:
            ref = cells[(sg, fam, 3.1)]
            same = all(cells[(sg, fam, rv)]["status"] == ref["status"]
                       and cells[(sg, fam, rv)]["agreeing_tests"] == ref["agreeing_tests"]
                       for rv in (3.0, 3.5))
            robustness[f"{sg}|{fam}"] = "robust" if same else "R_V-fragile"

    mix = tests[tests.test.eq("c_mixture")]
    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20b_score.py",
        "preregistration": {"path": "provenance/issue20b_prereg.json",
                            "sha256": w.sha256(B.PREREG_PATH)},
        "C_versus_A": c_vs_a, "C_versus_A_per_family": per_family,
        "cells_R_V_3.1": {f"{k[0]}|{k[1]}": v for k, v in cells.items() if k[2] == 3.1},
        "robustness": robustness,
        "H3_guarded_mixture_report_only": mix.drop(columns=[c for c in mix.columns
                                                            if mix[c].isna().all()]).to_dict(orient="records"),
        "calibrations": run["calibrations"],
        "integrity": run["integrity"] if "integrity" in run else None,
        "outputs": {str(OUT_TABLE.relative_to(w.ROOT)): w.sha256(OUT_TABLE)},
    })
    print(f"C versus A: {c_vs_a}  {per_family}")
    cols = ["subgroup", "family", "R_V", "status", "agreeing_tests", "a_map", "b_map", "c_map",
            "adopted_map", "adopted_lo68", "adopted_hi68"]
    print(table[[c for c in cols if c in table]].round(2).to_string(index=False))
    print(robustness)


if __name__ == "__main__":
    main()
