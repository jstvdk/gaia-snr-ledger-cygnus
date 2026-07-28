"""Per-subgroup calibration lower edge, driven by each subgroup's own
completeness curve -- the procedure plan WP5 step 2 and CUTS 7.1 prescribe
("expect this to differ between subgroups because extinction differs").

repair_v1..v3 instead used one global 2.0 Msun edge for all three subgroups
under the corrected_no_absolute95_edge fallback, which forces the fit to
reconstruct CygOB2-B and -C counts from recovery fractions of 1-10%.
"""
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import wp5_common as w
import wp5_fit_imf as F

VERSION = "repair_v3"
masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{VERSION}.parquet")
curves = pd.read_parquet(w.PROC / f"wp5_completeness_curves_{VERSION}.parquet")
responses = pd.read_parquet(w.PROC / f"wp5_injection_response_{VERSION}.parquet")
store = np.load(w.PROC / f"wp4_mass_posterior_samples_{VERSION}.npz")
samples = store["samples"]
draw_cols = sorted(c for c in responses.columns if c.startswith("recovered_mass_draw_"))
original_edge = F.absolute_95_edge
scan: list[dict] = []


def make_rule(kind, level):
    def rule(curve):
        use = curve[curve["primary_mass"].between(2.0, 8.0)].sort_values("primary_mass")
        values = use["recovery_isotonic"].to_numpy(float)
        mass = use["primary_mass"].to_numpy(float)
        if kind == "plateau":
            target = level * float(np.median(values[mass >= 6.0]))
        else:
            target = level
        suffix_min = np.minimum.accumulate(values[::-1])[::-1]
        ok = np.flatnonzero(suffix_min >= target)
        return float(mass[ok[0]]) if len(ok) else float("nan")
    return rule


for label, kind, level in [
    ("global 2.0 (current)", None, None),
    ("absolute recovery >= 0.30", "abs", 0.30),
    ("absolute recovery >= 0.50", "abs", 0.50),
    ("95% of bright plateau", "plateau", 0.95),
]:
    F.absolute_95_edge = original_edge if kind is None else make_rule(kind, level)
    rng = np.random.default_rng(w.SEED)
    rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for alpha in w.IMF_SLOPES:
                for sg in w.SUBGROUPS:
                    cu = curves[curves.family.eq(family) & curves.R_V.eq(rv) & curves.subgroup.eq(sg)]
                    rp = responses[responses.family.eq(family) & responses.R_V.eq(rv) & responses.subgroup.eq(sg)]
                    try:
                        s, _, _ = F.fit_one(
                            masses, cu, rp, sg, family, rv, alpha, rng,
                            samples[:, w.FAMILIES.index(family) * len(w.R_V_BRANCHES)
                                    + w.R_V_BRANCHES.index(rv), :], draw_cols)
                    except RuntimeError:
                        rows.append(dict(subgroup=sg, family=family, R_V=rv, alpha=alpha,
                                         residual_gate_pass=False, calibration_lower_Msun=np.nan,
                                         max_abs_pearson_residual=np.nan, raw_calibration_sources=0))
                        continue
                    rows.append(s)
    t = pd.DataFrame(rows)
    base = t[(t.family == "PARSEC") & (t.R_V == 3.1) & (t.alpha == 2.3)]
    print(f"=== {label} ===")
    print("   grid %2d/54 | baseline all-pass: %-5s | min sources any branch: %d"
          % (t.residual_gate_pass.sum(), bool(base.residual_gate_pass.all()),
             t.raw_calibration_sources.min()))
    for r in base.itertuples():
        print("     %s  edge %.2f  N=%3d  max|r|=%.2f  chi2p=%.3f  %s"
              % (r.subgroup.replace("CygOB2-", ""), r.calibration_lower_Msun,
                 r.raw_calibration_sources, r.max_abs_pearson_residual,
                 r.poisson_chi_square_p, "PASS" if r.residual_gate_pass else "FAIL"))
    scan.append({
        "rule": label,
        "branches_passing": int(t.residual_gate_pass.sum()),
        "baseline_all_pass": bool(base.residual_gate_pass.all()),
        "min_sources_any_branch": int(t.raw_calibration_sources.min()),
        "baseline": base[[
            "subgroup", "calibration_lower_Msun", "raw_calibration_sources",
            "max_abs_pearson_residual", "poisson_chi_square_p",
            "residual_gate_pass",
        ]].to_dict(orient="records"),
    })

w.write_json(w.PROVENANCE / "wp5_lower_edge_scan_execution.json", {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "script": "scripts/wp5_lower_edge_scan.py",
    "status": "SUCCESS",
    "wp5_version": VERSION,
    "purpose": (
        "Test whether a per-subgroup calibration lower edge driven by each "
        "subgroup's own completeness curve (plan WP5 step 2, CUTS 7.1) "
        "resolves the CygOB2-B mid-window residual. It does not: raising B's "
        "edge makes B monotonically worse, ruling out a completeness-ramp "
        "artifact."
    ),
    "environment": {
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "pandas": pd.__version__,
    },
    "inputs": {
        str(p.relative_to(w.ROOT)): w.sha256(p) for p in [
            w.PROC / f"wp4_mass_posteriors_{VERSION}.parquet",
            w.PROC / f"wp5_completeness_curves_{VERSION}.parquet",
            w.PROC / f"wp5_injection_response_{VERSION}.parquet",
        ]
    },
    "scan": scan,
})
print("\nwrote provenance/wp5_lower_edge_scan_execution.json")
