#!/usr/bin/env python3
"""repair_v9 sensitivity S1 -- the IMF upper limit m_max in {100, 120, 150} Msun
(brief §4.5; pre-registered in provenance/repair_v9_s1_s2_prereg.json).

IMF_UPPER_LIMIT = 120 is a convention, not a fit (wp6_mass_extension_decision).
It enters the ledger twice: the Poisson mean k * integral[M_min, m_max] M^-alpha
and the lifetime inversion tau(m), which must bracket m_max.  k itself is the
normalisation of dN/dM = k M^-alpha fitted below 8 Msun, so it does not depend
on m_max.  The WP7 engine is used unmodified: wp7_ledger reads IMF_UPPER_LIMIT
from its module namespace at call time, so this script sets that name before
building each TurnoffRelation and running each population.

Issue #14: above ~120 Msun the isochrone tables' maximum is a table ceiling at
young ages (PARSEC 300, MIST 210).  150 Msun is reached at 2.82 Myr (PARSEC)
and between 2.00 and 2.25 Myr (MIST), where the tabulated maximum still falls
with age, i.e. is an evolutionary turnoff; TurnoffRelation refuses to build if
the grid does not bracket m_max.

Branches: the 18 alpha = 2.3 branches (family x R_V x star-formation spread),
each subgroup and the association, all-explode, 2,000,000 iterations, on the
repair_v9 WP5 posterior draws.  Same seed for every m_max (paired draws).

Outputs:  tables/repair_v9_s1_imf_ceiling.csv
          provenance/repair_v9_s1_execution.json
Run:      PYTHONPATH=scripts python3 scripts/repair_v9_s1_imf_ceiling.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp7_ledger as L
from wp7_ledger_prereg import RECENT_WINDOW_MYR, SF_DURATIONS_MYR

M_MAX = (100.0, 120.0, 150.0)
ALPHA = 2.3
ITERATIONS = 2_000_000
WP5_VERSION = "repair_v9"
OUT = w.TABLES / "repair_v9_s1_imf_ceiling.csv"
OUT_JSON = w.PROVENANCE / "repair_v9_s1_execution.json"
PREREG = w.PROVENANCE / "repair_v9_s1_s2_prereg.json"


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; nothing is overwritten")
    draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{WP5_VERSION}.npz")
    original = L.IMF_UPPER_LIMIT
    rows = []
    try:
        for m_max in M_MAX:
            L.IMF_UPPER_LIMIT = m_max
            relations = {f: L.TurnoffRelation(f) for f in w.FAMILIES}
            for fi, fam in enumerate(w.FAMILIES):
                for ri, rv in enumerate(w.R_V_BRANCHES):
                    for di, delta in enumerate(SF_DURATIONS_MYR):
                        total = np.zeros(ITERATIONS, dtype=int)
                        last = np.full(ITERATIONS, np.inf)
                        for si, sg in enumerate(w.SUBGROUPS):
                            key = L.draw_key(sg, fam, rv, ALPHA)
                            k_all, age_all = draws[f"k__{key}"], draws[f"truth_age_draws__{key}"]
                            rng = np.random.default_rng([w.SEED, 1, fi, ri, di, si])  # same for every m_max
                            pick = rng.integers(0, k_all.size, ITERATIONS)
                            res = L.run_population(rng, k_all[pick], age_all[pick], ALPHA, delta,
                                                   relations[fam])
                            n, t = res["n_sn"]["all_explode"], res["t_last"]["all_explode"]
                            total += n
                            last = np.minimum(last, t)
                            rows.append(row(m_max, fam, rv, delta, sg, n, t))
                        rows.append(row(m_max, fam, rv, delta, "ALL", total, last))
            print(f"m_max {m_max:.0f} done", flush=True)
    finally:
        L.IMF_UPPER_LIMIT = original
    table = pd.DataFrame(rows)
    ref = table[np.isclose(table.m_max_Msun, 120.0)].set_index(
        ["family", "R_V", "sf_duration_Myr", "subgroup"]).N_death_mean
    table["ratio_to_120"] = [r.N_death_mean / ref[(r.family, r.R_V, r.sf_duration_Myr, r.subgroup)]
                             if ref[(r.family, r.R_V, r.sf_duration_Myr, r.subgroup)] > 0 else np.nan
                             for r in table.itertuples()]
    table.to_csv(OUT, index=False)
    w.write_json(OUT_JSON, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_s1_imf_ceiling.py",
        "preregistration": {"path": str(PREREG.relative_to(w.ROOT)), "sha256": w.sha256(PREREG)},
        "m_max_Msun": list(M_MAX), "alpha": ALPHA, "iterations": ITERATIONS,
        "inputs": {f"data/processed/wp5_imf_posterior_draws_{WP5_VERSION}.npz":
                   w.sha256(w.PROC / f"wp5_imf_posterior_draws_{WP5_VERSION}.npz")},
        "outputs": {str(OUT.relative_to(w.ROOT)): w.sha256(OUT)},
    })
    base = table[table.family.eq("PARSEC") & np.isclose(table.R_V, 3.1)
                 & np.isclose(table.sf_duration_Myr, 0.0)]
    print(base.pivot(index="subgroup", columns="m_max_Msun", values="N_death_mean").round(3))


def row(m_max, fam, rv, delta, sg, n, t) -> dict:
    return {"m_max_Msun": m_max, "family": fam, "R_V": rv, "alpha": ALPHA,
            "sf_duration_Myr": delta, "subgroup": sg,
            "N_death_mean": float(n.mean()), "N_death_median": float(np.median(n)),
            "N_death_p16": float(np.percentile(n, 16)), "N_death_p84": float(np.percentile(n, 84)),
            "P_at_least_one": float((n >= 1).mean()),
            "P_last_within_100kyr": float((t < RECENT_WINDOW_MYR).mean())}


if __name__ == "__main__":
    main()
