#!/usr/bin/env python3
"""Issue #20 Phase A -- score the hypotheses by the frozen decision rule.

Reads the Phase A tables and execution records and applies
provenance/issue20_prereg.json "decision_rule" and "secondary_predictions"
mechanically.  Nothing is re-fitted here.

Output: provenance/issue20_phase_a_outcome.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20_score.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I
import wp5_common as w

OUT = w.PROVENANCE / "issue20_phase_a_outcome.json"
A, C = "CygOB2-A", "CygOB2-C"


def load(rel):
    return json.loads((w.ROOT / rel).read_text())


def family_outcome(post: pd.DataFrame, mix: pd.DataFrame) -> dict:
    """Per-family outcome on one settings cell: which hypotheses hold."""
    a = post[post.subgroup.eq(A)].iloc[0]
    c = post[post.subgroup.eq(C)].iloc[0]
    m = mix[mix.subgroup.eq(C)].iloc[0]
    railed = bool(a.railed or c.railed)
    h3 = bool(m.h3_criteria_met)
    h1 = (not railed) and c.hi95 < a.lo95
    h2 = (not railed) and not (c.hi68 < a.lo68 or a.hi68 < c.lo68)
    label = "H3" if h3 else ("H1" if h1 else ("H2" if h2 else "none"))
    single = "H1" if h1 else ("H2" if h2 else "none")
    return {"H1": bool(h1), "H2": bool(h2), "H3": h3, "railed": railed,
            "label": label, "single_age_label": single,
            "A": {k: float(a[k]) for k in ("map", "lo95", "lo68", "hi68", "hi95")},
            "C": {k: float(c[k]) for k in ("map", "lo95", "lo68", "hi68", "hi95")},
            "C_mixture": {k: (float(m[k]) if pd.notna(m[k]) else None) for k in
                          ("single_best_age", "age_young", "age_old", "phi_young",
                           "n_young", "n_old", "delta_BIC")}}


def combine(per_family: dict) -> str:
    labels = {f: o["label"] for f, o in per_family.items()}
    if all(o["H3"] for o in per_family.values()):
        return "H3"
    if all(o["H1"] for o in per_family.values()):
        return "H1"
    if all(o["H2"] for o in per_family.values()):
        return "H2"
    return "INCONCLUSIVE"


def single_age_combine(per_family: dict) -> str:
    if all(o["H1"] for o in per_family.values()):
        return "H1"
    if all(o["H2"] for o in per_family.values()):
        return "H2"
    return "INCONCLUSIVE"


def cell(post, mix, **sel):
    p, m = post, mix
    for k, v in sel.items():
        p = p[np.isclose(p[k], v)] if isinstance(v, float) else p[p[k].eq(v)]
        m = m[np.isclose(m[k], v)] if isinstance(v, float) else m[m[k].eq(v)]
    return p, m


def main() -> None:
    I.check_method_unchanged()
    pre = I.prereg()
    post = pd.read_csv(w.TABLES / "issue20_hrd_posteriors.csv")
    mix = pd.read_csv(w.TABLES / "issue20_hrd_mixture.csv")
    base = dict(sigma_mode="primary", cal=I.SIG_CAL_MAG, dmu=0.0, teff_mode="supergiant_scale")

    # ---------------------------------------------------------- the verdict
    decision, joint_decision = {}, {}
    for fam in I.FAMILIES:
        p, m = cell(post, mix, family=fam, R_V=3.1, f_bin=0.4, statistic="conditional", **base)
        decision[fam] = family_outcome(p, m)
        p, m = cell(post, mix, family=fam, R_V=3.1, f_bin=0.4, statistic="joint", **base)
        joint_decision[fam] = family_outcome(p, m)
    verdict = combine(decision)
    joint_verdict = combine(joint_decision)

    # robustness: every other R_V x f_bin cell of the registered grid
    robustness, dissent = [], []
    for fam in I.FAMILIES:
        for rv in I.R_V_BRANCHES:
            for fb in (0.3, 0.4, 0.5):
                if rv == 3.1 and fb == 0.4:
                    continue
                p, m = cell(post, mix, family=fam, R_V=rv, f_bin=fb,
                            statistic="conditional", **base)
                o = family_outcome(p, m)
                robustness.append({"family": fam, "R_V": rv, "f_bin": fb,
                                   "label": o["label"], "single_age_label": o["single_age_label"],
                                   "C_map": o["C"]["map"], "A_map": o["A"]["map"],
                                   "C_mix_dBIC": o["C_mixture"]["delta_BIC"]})
                if o["label"] != decision[fam]["label"]:
                    dissent.append(robustness[-1])

    # H3 specificity control: the same criteria on A
    a_mix = {}
    for fam in I.FAMILIES:
        _, m = cell(post, mix, family=fam, R_V=3.1, f_bin=0.4, statistic="conditional", **base)
        r = m[m.subgroup.eq(A)].iloc[0]
        a_mix[fam] = {"h3_criteria_met": bool(r.h3_criteria_met), "delta_BIC": float(r.delta_BIC),
                      "age_young": float(r.age_young), "age_old": float(r.age_old),
                      "phi_young": float(r.phi_young)}
    a_robust = []
    for fam in I.FAMILIES:
        for rv in I.R_V_BRANCHES:
            for fb in (0.3, 0.4, 0.5):
                _, m = cell(post, mix, family=fam, R_V=rv, f_bin=fb, statistic="conditional", **base)
                r = m[m.subgroup.eq(A)].iloc[0]
                a_robust.append({"family": fam, "R_V": rv, "f_bin": fb,
                                 "h3_criteria_met": bool(r.h3_criteria_met),
                                 "delta_BIC": float(r.delta_BIC)})
    not_specific = all(v["h3_criteria_met"] for v in a_mix.values())

    # oldest native age: PARSEC 10.00, MIST 10.12 Myr (report only; the frozen
    # H3 criteria do not test for it)
    grid_edge = {fam: {"C_old_component_Myr": decision[fam]["C_mixture"]["age_old"],
                       "A_old_component_Myr": a_mix[fam]["age_old"],
                       "both_at_oldest_grid_age": bool(
                           decision[fam]["C_mixture"]["age_old"] >= 9.99
                           and a_mix[fam]["age_old"] >= 9.99)}
                 for fam in I.FAMILIES}

    qualifiers = {
        "robust": not dissent,
        "dissenting_cells": dissent,
        "statistic_dependent": joint_verdict != verdict,
        "joint_statistic_verdict": joint_verdict,
        "H3_not_specific_to_C": bool(verdict == "H3" and not_specific),
    }

    # ------------------------------------------------ secondary predictions
    audit = load("provenance/issue20_star_audit_execution.json")
    frac = audit["S1_measured"]["fraction_nearest_C"]
    s1 = ("consistent with H2/H3" if frac >= 2 / 3 else
          "consistent with H1" if frac <= 1 / 3 else "neutral")
    refits = pd.read_csv(w.TABLES / "issue20_photometric_refits.csv")

    def refit(step, fam, rv=3.1, sub=C):
        r = refits[refits.step.eq(step) & refits.family.eq(fam) & refits.R_V.eq(rv)
                   & refits.subgroup.eq(sub)]
        return float(r.age_map.iloc[0])
    s2_rows = {fam: {"C_map_without_I_II": refit("A4b", fam),
                     "A_stored_map": refit("A_reference", fam, sub=A),
                     "young_persists": refit("A4b", fam) < refit("A_reference", fam, sub=A) - 0.5}
               for fam in I.FAMILIES}
    s2 = all(v["young_persists"] for v in s2_rows.values())
    s3_rows = {fam: {"C_hybrid_map": refit("A4a", fam), "old": refit("A4a", fam) >= 3.3}
               for fam in I.FAMILIES}
    s3 = all(v["old"] for v in s3_rows.values())
    clo = pd.read_csv(w.TABLES / "issue20_closure_by_c_age.csv")
    ages = [2.51, 3.16, 3.55, 3.98]
    mono = []
    for (fam, rv), blk in clo[clo.alpha.eq(2.3) & clo.C_age_Myr.notna()].groupby(["family", "R_V"]):
        seq = [float(blk[np.isclose(blk.C_age_Myr, a)].closure_ratio.iloc[0]) for a in ages]
        mono.append({"family": fam, "R_V": rv, "ratios": seq,
                     "non_decreasing": bool(np.all(np.diff(seq) >= -1e-12))})
    slopes = {a: float(clo[np.isclose(clo.C_age_Myr, a)].drop_duplicates(["family", "R_V"])
                       .closing_alpha.median()) for a in ages}
    s4 = all(r["non_decreasing"] for r in mono) and slopes[3.98] < slopes[2.51]

    secondary = {
        "S1": {"measured_fraction_nearest_C": frac,
               "n_class_I_II": len(audit["S1_measured"]["class_I_II_stars"]),
               "outcome": s1, "caveat": "partly circular (GMM labels in l, b, pm)",
               "stars": audit["S1_measured"]["class_I_II_stars"]},
        "S2": {"rows": s2_rows, "outcome": "young persists (consistent with H1)" if s2
               else "young does not persist"},
        "S3": {"rows": s3_rows, "outcome": "PASS (consistent with H2)" if s3 else
               "FAIL (hybrid not >= 3.3 Myr in both families)"},
        "S4": {"rows": mono, "median_closing_alpha_by_age": slopes,
               "outcome": "PASS" if s4 else "FAIL"},
    }

    # --------------------------------------------------------- integrity
    hrd = load("provenance/issue20_hrd_likelihood_execution.json")
    pho = load("provenance/issue20_photometric_execution.json")
    clx = load("provenance/issue20_closure_execution.json")
    integrity = {"I1": pho["integrity"]["I1"]["pass"], "I2": pho["integrity"]["I2"]["pass"],
                 "I3": hrd["integrity"]["I3"]["pass"], "I4": clx["integrity"]["I4"]["pass"],
                 "I5": hrd["integrity"]["I5"]["pass"],
                 "I6": "all inputs read through issue20_common.frozen; no FrozenInputMoved raised"}

    phase_b = {
        "H1": "ledger stands; Phase C opens",
        "H2": "repair_v9 (separately pre-registered, after the WP13 decision)",
        "H3": ("needs a user decision: the mixture preference is not specific to C"
               if qualifiers["H3_not_specific_to_C"] else
               "two-component C in repair_v9 or the formation-duration branch"),
        "INCONCLUSIVE": "carry C as a young/coeval branch; no novelty claim on C's age",
    }[verdict]

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_score.py",
        "preregistration": {"path": "provenance/issue20_prereg.json",
                            "created_utc": pre["created_utc"],
                            "sha256": w.sha256(I.PREREG_PATH)},
        "verdict": verdict,
        "verdict_rule_text": pre["decision_rule"],
        "qualifiers": qualifiers,
        "single_age_verdict_ignoring_H3": single_age_combine(decision),
        "decision_cells": decision,
        "joint_statistic_decision_cells": joint_decision,
        "robustness_cells": robustness,
        "H3_specificity_control_A": {"decision_cells": a_mix, "all_cells": a_robust},
        "mixture_old_component": grid_edge,
        "secondary_predictions": secondary,
        "integrity": integrity,
        "phase_b_action": phase_b,
        "consequence_for_n_death_registered": pre["consequences_for_n_death"].get(verdict),
    }
    w.write_json(OUT, record)
    print(f"VERDICT: {verdict}  (single-age, ignoring H3: {record['single_age_verdict_ignoring_H3']})")
    print(f"qualifiers: robust={qualifiers['robust']} statistic_dependent="
          f"{qualifiers['statistic_dependent']} H3_not_specific_to_C="
          f"{qualifiers['H3_not_specific_to_C']}")
    for fam, o in decision.items():
        print(f"  {fam}: A {o['A']['map']:.2f} [{o['A']['lo68']:.2f}, {o['A']['hi68']:.2f}]  "
              f"C {o['C']['map']:.2f} [{o['C']['lo68']:.2f}, {o['C']['hi68']:.2f}]  "
              f"H1 {o['H1']} H2 {o['H2']} H3 {o['H3']} (dBIC {o['C_mixture']['delta_BIC']:.1f}, "
              f"ages {o['C_mixture']['age_young']}/{o['C_mixture']['age_old']})")
    for k, v in secondary.items():
        print(f"  {k}: {v['outcome']}")
    print(f"  integrity: {integrity}")


if __name__ == "__main__":
    main()
