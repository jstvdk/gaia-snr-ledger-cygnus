#!/usr/bin/env python3
"""Issue #20 Phase A, step A6 -- C's closure ratio and closing slope at fixed
C ages 2.51, 3.16, 3.55 and 3.98 Myr.

Replicates the wp6_closure_test.py estimator for CygOB2-C and varies ONLY the
upper limit of the forward integral: cap = min(turnoff(family, t), 120 Msun)
applied to every one of C's stored node responses (provenance/issue20_prereg.json
"a6_method").  C's stored repair_v7 node responses (reused by repair_v8), their
WP4 node weights, C's repair_v8 k_median per alpha and C's repair_v8 observed
census are held fixed.  The age dependence of k and of the mass estimates of
C's observed stars is therefore NOT captured: that needs new injections at
3.2-4.0 Myr (repair_v9), which a read-only step cannot make.

I4: with the stored node-wise turnoff caps the replica reproduces C's rows of
tables/wp6_closure_repair_v8.csv.

Outputs:
  tables/issue20_closure_by_c_age.csv
  provenance/issue20_closure_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/issue20_closure_c_age.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

import issue20_common as I
import wp5_common as w
import wp5_joint_age_fit as J
from wp6_closure_attribution import closing_alpha
from wp6_closure_test import (INTEGRATION_FLOOR_MSUN, SN_THRESHOLD_MSUN,
                              observed_response, trapezoid_weights)
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass
from wp6_massive_injections import response_path as extension_path

OUT = w.TABLES / "issue20_closure_by_c_age.csv"
OUT_EXEC = w.PROVENANCE / "issue20_closure_execution.json"
SUBGROUP = "CygOB2-C"
RESPONSE_VERSION = "repair_v7"     # reused by the repair_v8 chain
C_AGES = (2.51, 3.16, 3.55, 3.98)


def node_responses(age_posterior, family, rv, native):
    prior = J.truth_age_nodes(age_posterior, SUBGROUP, family, rv, native,
                              snap=not J.uses_age_interpolation(RESPONSE_VERSION))
    nodes = []
    for age, weight in prior.items():
        base = J.node_response_path(SUBGROUP, family, rv, age, RESPONSE_VERSION)
        ext = extension_path(SUBGROUP, family, rv, age, RESPONSE_VERSION)
        combined = pd.concat([pd.read_parquet(base), pd.read_parquet(ext)],
                             ignore_index=True)
        draws = sorted(c for c in combined.columns if c.startswith("recovered_mass_draw_"))
        masses, above = observed_response(combined, SN_THRESHOLD_MSUN, draws)
        nodes.append({"age": float(age), "weight": float(weight),
                      "masses": masses, "above": above,
                      "files": [str(base.relative_to(w.ROOT)), str(ext.relative_to(w.ROOT))]})
    return nodes


def predicted_obs(nodes, k, alpha, cap_of_node):
    total = 0.0
    for node in nodes:
        cap = cap_of_node(node)
        window = (node["masses"] >= INTEGRATION_FLOOR_MSUN) & (node["masses"] <= cap)
        m, r = node["masses"][window], node["above"][window]
        total += node["weight"] * float(np.sum(trapezoid_weights(m) * m ** (-alpha) * r))
    return k * total


def main() -> None:
    method = I.check_method_unchanged()
    norm = pd.read_parquet(I.frozen("wp5_normalization_v8"))
    census = pd.read_csv(I.frozen("massive_census_v8"))
    stored = pd.read_csv(I.frozen("closure_v8"))
    age_posterior = pd.read_parquet(I.frozen("age_posteriors_v5"))
    native = {f: J.native_isochrone_ages(f) for f in w.FAMILIES}

    rows, i4, files = [], [], set()
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            nodes = node_responses(age_posterior, family, rv, native[family])
            for node in nodes:
                files.update(node["files"])
            observed = float(census[census.subgroup.eq(SUBGROUP) & census.family.eq(family)
                                    & census.R_V.eq(rv)].observed_above_8_probabilistic.iloc[0])
            variants = [("stored_node_prior", None)] + [(f"C_age_{t:.2f}", t) for t in C_AGES]
            for label, t in variants:
                if t is None:
                    def cap_of_node(node):
                        return min(turnoff_mass(family, node["age"]), IMF_UPPER_LIMIT)
                    cap_report = float(np.sum([n["weight"] * cap_of_node(n) for n in nodes]))
                else:
                    cap_t = min(turnoff_mass(family, t), IMF_UPPER_LIMIT)

                    def cap_of_node(node, cap_t=cap_t):
                        return cap_t
                    cap_report = cap_t
                ratios = {}
                for alpha in w.IMF_SLOPES:
                    cell = norm[norm.subgroup.eq(SUBGROUP) & norm.family.eq(family)
                                & norm.R_V.eq(rv) & norm.alpha.eq(alpha)]
                    k = float(cell.k_median.iloc[0])
                    pred = predicted_obs(nodes, k, alpha, cap_of_node)
                    ratios[alpha] = observed / pred
                    rows.append({"family": family, "R_V": rv, "variant": label,
                                 "C_age_Myr": t, "cap_Msun": cap_report, "alpha": alpha,
                                 "k_median": k, "observed_living": observed,
                                 "predicted_observed_living": pred,
                                 "closure_ratio": ratios[alpha]})
                    if t is None:
                        ref = stored[stored.subgroup.eq(SUBGROUP) & stored.family.eq(family)
                                     & stored.R_V.eq(rv) & stored.alpha.eq(alpha)]
                        ref_ratio = float(ref.closure_ratio.iloc[0])
                        i4.append({"family": family, "R_V": rv, "alpha": alpha,
                                   "replica": ratios[alpha], "stored": ref_ratio,
                                   "rel_diff": abs(ratios[alpha] / ref_ratio - 1.0)})
                slope = closing_alpha(np.array(list(ratios)), np.array(list(ratios.values())))
                for r in rows[-3:]:
                    r["closing_alpha"] = slope
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    i4_pass = bool(max(r["rel_diff"] for r in i4) <= 1e-9)

    summary = {}
    for label, block in table[table.alpha.eq(2.3)].groupby("variant"):
        slopes = table[table.variant.eq(label)].drop_duplicates(["family", "R_V"]).closing_alpha
        base = block[block.family.eq("PARSEC") & block.R_V.eq(3.1)]
        summary[label] = {"ratio_alpha2.3_PARSEC_rv3.1": float(base.closure_ratio.iloc[0]),
                          "closing_alpha_6cell_median": float(slopes.median())}
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue20_closure_c_age.py",
        "preregistration": "provenance/issue20_prereg.json",
        "method_check": method,
        "inputs": {n: I.prereg()["consumed_inputs"][n] for n in
                   ("wp5_normalization_v8", "massive_census_v8", "closure_v8",
                    "age_posteriors_v5")},
        "node_response_files": {f: w.sha256(w.ROOT / f) for f in sorted(files)},
        "held_fixed": "C's node responses and weights, k_median per alpha, observed census",
        "varies": "only the cap of the forward integral: min(turnoff(family, t), 120)",
        "summary": summary,
        "integrity": {"I4": {"pass": i4_pass, "rows": i4}},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    }
    w.write_json(OUT_EXEC, record)
    for label, s in summary.items():
        print(f"{label:18s} ratio(PARSEC 3.1, a=2.3) {s['ratio_alpha2.3_PARSEC_rv3.1']:.3f}  "
              f"median closing alpha {s['closing_alpha_6cell_median']:.3f}")
    print(f"I4 {'PASS' if i4_pass else 'FAIL'}")


if __name__ == "__main__":
    main()
