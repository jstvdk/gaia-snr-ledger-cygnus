#!/usr/bin/env python3
"""repair_v9 Stage 1 -- POST-HOC diagnostics of the real-data M1 bright-window fits.

NOT pre-registered; written after wp4v9_score.py had scored the stage.  Gates
nothing and changes no adopted row.  It explains why the real-data M1 ages are
unstable (reports/wp4v9_age_redesign.md, section 4):

  (1) per-star contributions P * [lnL(young) - lnL(old)] in the bright window,
      R_V 3.1, f_bin 0.4, both families;
  (2) two refits of M1 bright at every R_V (f_bin 0.4): with the outlier term
      (eps 0.05, density 0.2, as the M4 faint-end study) and without the
      composite spectroscopic members.

Outputs:
  tables/wp4v9_posthoc_star_contributions.csv
  tables/wp4v9_posthoc_refits.csv

Run:
  PYTHONPATH=scripts python3 scripts/wp4v9_posthoc_diagnostics.py
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import issue20_common as I20
import wp4v9_common as V
import wp5_common as w

OUT_STARS = w.TABLES / "wp4v9_posthoc_star_contributions.csv"
OUT_REFITS = w.TABLES / "wp4v9_posthoc_refits.csv"
# young / old comparison ages per subgroup (the real-data M1 MAPs vs the spectroscopic ages)
PAIRS = {"CygOB2-A": (1.2, 3.6), "CygOB2-B": (1.3, 4.5), "CygOB2-C": (1.0, 3.6)}


def _stars(rv):
    ext = pd.read_parquet(V.frozen("extinction_v5"))
    labels = pd.read_parquet(V.frozen("subgroup_labels"))
    npz = np.load(V.frozen("av_posterior_v5"))
    return V.load_stars(ext, labels, npz, rv)


def _bright(st, sg, exclude=()):
    s = V.in_window(st[st.subgroup.eq(sg)], "bright")
    s = s[~s.source_id.isin(exclude)].copy()
    s.attrs.update(st.attrs)
    return s


def run(args):
    fam, rv = args
    st = _stars(rv)
    anc = pd.read_parquet(V.frozen("anchor_hrd_v8"))
    composite = set(anc.source_id[[I20.parse_spectral_type(t)["composite"] for t in anc.spectral_type]])
    ages, plist = V.particles_by_age(pd.read_parquet(V.frozen(f"isochrones_{fam.lower()}")), 0.4)
    refits, contrib = [], []
    for sg in V.SUBGROUPS:
        s0 = _bright(st, sg)
        if rv == 3.1:
            ty, to = PAIRS[sg]
            jy, jo = int(np.argmin(abs(ages - ty))), int(np.argmin(abs(ages - to)))
            d = s0.P.values * (V.m1_loglike(s0, plist[jy], "bright") - V.m1_loglike(s0, plist[jo], "bright"))
            c = s0[["source_id", "is_anchor", "MG0_median", "c_obs", "av_median", "av_halfwidth68", "P"]].copy()
            c = c.assign(family=fam, subgroup=sg, age_young=float(ages[jy]), age_old=float(ages[jo]), dlnL=d)
            contrib.append(c.merge(anc[["source_id", "spectral_type"]], on="source_id", how="left"))
        for label, s, eps in (("outlier_eps0.05", s0, 0.05),
                              ("no_composite_anchors", _bright(st, sg, composite), 0.0)):
            tot = np.array([s.P.values @ V.m1_loglike(s, p, "bright", 0.0, eps, 0.2) for p in plist])
            p = V.posterior(ages, tot, len(s))
            refits.append({"family": fam, "R_V": rv, "f_bin": 0.4, "subgroup": sg, "variant": label,
                           "n_stars": len(s), **p})
    return refits, contrib


def main() -> None:
    V.check_method_unchanged()
    with ProcessPoolExecutor(6) as ex:
        res = list(ex.map(run, [(f, rv) for f in V.FAMILIES for rv in V.R_V_BRANCHES]))
    pd.DataFrame([r for a, _ in res for r in a]).to_csv(OUT_REFITS, index=False)
    pd.concat([c for _, b in res for c in b], ignore_index=True).to_csv(OUT_STARS, index=False)
    print("written", OUT_REFITS.name, OUT_STARS.name)


if __name__ == "__main__":
    main()
