#!/usr/bin/env python3
"""S3 -- the single authorized source of manuscript inputs, with a hard guard.

Issue #2 left `tables/wp5_imf_norm.csv` and `wp5_imf_norm.md` on disk as the
frozen pre-repair record in which 0 of 54 branches passed the mass-function
gate.  They are correct as history and catastrophic as manuscript inputs: they
contradict every accepted number in the chain.  The same hazard exists for every
WP5/WP6 product that has both an unversioned original and a `_repair_vN`
successor.

This module makes the safe path the only path.  WP10 resolves every input
through `resolve()`, which refuses anything on the forbidden list by raising,
and `audit()` greps the manuscript sources for forbidden references so a stray
hand-written path is caught too.  Import it; do not open manuscript inputs
directly.

Run standalone to write the manifest and run the audit:
  PYTHONPATH=scripts python3 scripts/wp10_inputs.py
"""
from __future__ import annotations

import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import wp5_common as w

WP5_VERSION = "repair_v7"
WP5_REPORT_VERSION = "repair_v6"   # the accepted gate record; v7 re-passed it
WP3_WP4_VERSION = "repair_v5"

# --------------------------------------------------------------- authorized
# Logical name -> path relative to the repository root.  Every WP5/WP6 entry is
# version-suffixed; nothing unversioned from those work packages appears here.
MANUSCRIPT_INPUTS: dict[str, str] = {
    # WP1-WP2 -- never repaired, so unversioned is correct for these
    "wp1_inventory": "tables/table1_wp1_inventory.md",
    "wp2_gate": "tables/table2_wp2_gate.md",
    "wp2_literature_recovery": "tables/table3_literature_recovery.md",
    "wp2_members": "data/processed/wp2_members.parquet",
    "wp2_subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    # WP3 -- repaired
    "wp3_extinction": f"data/processed/wp3_extinction_{WP3_WP4_VERSION}.parquet",
    # WP4 -- repaired
    "wp4_age_posteriors": f"data/processed/wp4_age_posteriors_{WP3_WP4_VERSION}.parquet",
    "wp4_masses": f"data/processed/wp4_mass_posteriors_{WP3_WP4_VERSION}.parquet",
    # WP5 -- repaired twice; the normalization consumed downstream is v7
    "wp5_normalization": f"data/processed/wp5_imf_normalization_{WP5_VERSION}.parquet",
    "wp5_posterior_draws": f"data/processed/wp5_imf_posterior_draws_{WP5_VERSION}.npz",
    "wp5_association_mass": f"data/processed/wp5_association_mass_{WP5_VERSION}.parquet",
    "wp5_imf_norm_table": f"tables/wp5_imf_norm_{WP5_REPORT_VERSION}.csv",
    "wp5_baseline_residuals": f"tables/wp5_baseline_residuals_{WP5_REPORT_VERSION}.csv",
    "wp5_gate_record": f"provenance/wp5_{WP5_REPORT_VERSION}_gate.json",
    "wp5_association_mass_reconciliation": "tables/wp5_association_mass_reconciliation.csv",
    # WP6 -- closure re-run under repair_v7
    "wp6_closure": "tables/wp6_closure_repair_v7.csv",
    "wp6_closure_attribution": "tables/wp6_closure_attribution_repair_v7.csv",
    "wp6_massive_census": "tables/wp6_massive_census.csv",
    "wp6_runaways": "tables/wp6_runaways.csv",
    "wp6_runaway_crossmatch": "tables/wp6_runaway_crossmatch.csv",
    "wp6_orphan_anchors": "tables/wp6_orphan_anchors.csv",
    "wp6_external_crosschecks": "tables/wp6_external_crosschecks.csv",
    # WP7-WP9 -- single versions, all computed on the repair_v7 chain
    "wp7_ledger": "tables/wp7_ledger.csv",
    "wp7_rsn_curves": "tables/wp7_rsn_curves.csv",
    "wp7_age_sensitivity": "tables/wp7_age_sensitivity.csv",
    "wp7_bh_threshold_scan": "tables/wp7_bh_threshold_scan.csv",
    "wp7_convergence": "tables/wp7_convergence.csv",
    "wp8_crosschecks": "tables/wp8_crosschecks.csv",
    "wp8_tension_list": "tables/wp8_tension_list.csv",
    "wp9_verdict": "tables/wp9_verdict.csv",
    "wp9_sensitivity": "tables/wp9_sensitivity.csv",
    # pre-WP10 work
    "age_reconciliation": "tables/wp4_wp5_age_reconciliation.csv",
    "alpha_headline_branch_sets": "tables/wp7_alpha_headline_branch_sets.csv",
    "binary_bound": "tables/wp7_binary_bound.csv",
    "binary_bound_branches": "tables/wp7_binary_bound_branches.csv",
    "binary_bound_harer_fig2": "tables/wp7_binary_bound_harer_fig2.csv",
}

# ---------------------------------------------------------------- forbidden
FORBIDDEN: dict[str, str] = {
    "tables/wp5_imf_norm.csv": (
        "issue #2 -- frozen pre-repair WP5 run, 0 of 54 branches passing; "
        "contradicts every accepted number.  Use wp5_imf_norm_repair_v6.csv "
        "or the repair_v7 normalization parquet."
    ),
    "tables/wp5_imf_norm.md": (
        "issue #2 -- markdown copy of the same 0/54 run"
    ),
    "wp5_imf_norm.md": (
        "issue #2 -- the pre-repair WP5 report, still headed 'BLOCKED AT THE "
        "VALIDATION GATE'"
    ),
    "wp5_completion_report.md": (
        "the pre-repair completion report; WP5 was accepted at repair_v6"
    ),
    "data/processed/wp5_imf_normalization.parquet": (
        "unversioned pre-repair normalization"
    ),
    "data/processed/wp5_association_mass.parquet": (
        "unversioned pre-repair association mass"
    ),
    "data/processed/wp5_imf_posterior_draws.npz": (
        "unversioned pre-repair posterior draws"
    ),
    "tables/wp5_association_mass.csv": (
        "unversioned pre-repair association mass"
    ),
    "tables/wp6_closure.csv": (
        "superseded by wp6_closure_repair_v7.csv; the pre-repair_v7 closure "
        "ratios were withdrawn (issues #16, #17)"
    ),
    "tables/wp6_closure_attribution.csv": (
        "superseded by wp6_closure_attribution_repair_v7.csv"
    ),
    "data/processed/wp2_members_failed_20260722.parquet": (
        "the rejected 2026-07-22 membership run"
    ),
}

# Files whose text is scanned for forbidden references by audit().
AUDITED_GLOBS = ("manuscript/**/*.tex", "manuscript/**/*.md", "manuscript/**/*.py")


class ForbiddenInput(RuntimeError):
    """Raised when the manuscript asks for a superseded artifact."""


def resolve(name: str) -> Path:
    """Return the absolute path of an authorized manuscript input."""
    if name in FORBIDDEN or name in {v for v in FORBIDDEN}:
        raise ForbiddenInput(f"{name}: {FORBIDDEN[name]}")
    if name not in MANUSCRIPT_INPUTS:
        raise KeyError(
            f"{name!r} is not a declared manuscript input.  Add it to "
            f"MANUSCRIPT_INPUTS in scripts/wp10_inputs.py -- deliberately, and "
            f"with its version suffix -- rather than opening a path directly."
        )
    relative = MANUSCRIPT_INPUTS[name]
    if relative in FORBIDDEN:
        raise ForbiddenInput(f"{relative}: {FORBIDDEN[relative]}")
    path = w.ROOT / relative
    if not path.exists():
        raise FileNotFoundError(f"declared manuscript input is missing: {relative}")
    return path


def audit() -> dict:
    """Check the declared inputs exist and no source references a forbidden one."""
    missing = [
        name for name, rel in MANUSCRIPT_INPUTS.items()
        if not (w.ROOT / rel).exists()
    ]
    overlap = sorted(set(MANUSCRIPT_INPUTS.values()) & set(FORBIDDEN))

    violations = []
    patterns = {
        rel: re.compile(re.escape(Path(rel).name) + r"(?![\w.-])")
        for rel in FORBIDDEN
    }
    scanned = []
    for pattern_glob in AUDITED_GLOBS:
        for path in w.ROOT.glob(pattern_glob):
            if not path.is_file():
                continue
            scanned.append(str(path.relative_to(w.ROOT)))
            text = path.read_text(errors="ignore")
            for line_number, line in enumerate(text.splitlines(), start=1):
                if line.lstrip().startswith("%") or "wp10_inputs" in line:
                    continue
                for rel, pattern in patterns.items():
                    # A versioned name contains the unversioned one as a
                    # prefix, so match the bare filename with no suffix after.
                    if pattern.search(line):
                        violations.append(
                            {
                                "file": str(path.relative_to(w.ROOT)),
                                "line": line_number,
                                "forbidden": rel,
                                "reason": FORBIDDEN[rel],
                                "text": line.strip()[:200],
                            }
                        )
    return {
        "declared_inputs": len(MANUSCRIPT_INPUTS),
        "missing_declared_inputs": missing,
        "forbidden_entries": len(FORBIDDEN),
        "authorized_forbidden_overlap": overlap,
        "files_scanned": scanned,
        "violations": violations,
        "pass": not missing and not overlap and not violations,
    }


def manifest() -> dict:
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp10_inputs.py",
        "item": "S3 of tasks/pre_wp10_assessment_brief.md",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {
            "wp3_wp4": WP3_WP4_VERSION,
            "wp5_consumed_downstream": WP5_VERSION,
            "wp5_gate_record": WP5_REPORT_VERSION,
        },
        "authorized_inputs": {
            name: {
                "path": rel,
                "sha256": w.sha256(w.ROOT / rel) if (w.ROOT / rel).exists() else None,
            }
            for name, rel in sorted(MANUSCRIPT_INPUTS.items())
        },
        "forbidden_inputs": FORBIDDEN,
        "audit": audit(),
    }


def main() -> None:
    record = manifest()
    w.write_json(w.PROVENANCE / "wp10_input_manifest.json", record)
    result = record["audit"]
    print(f"WP10 input manifest: {result['declared_inputs']} authorized, "
          f"{result['forbidden_entries']} forbidden")
    if result["missing_declared_inputs"]:
        print("  MISSING:", ", ".join(result["missing_declared_inputs"]))
    if result["authorized_forbidden_overlap"]:
        print("  OVERLAP:", ", ".join(result["authorized_forbidden_overlap"]))
    if result["violations"]:
        print(f"  {len(result['violations'])} forbidden references in sources:")
        for v in result["violations"][:20]:
            print(f"    {v['file']}:{v['line']}  ->  {v['forbidden']}")
    print(f"  scanned {len(result['files_scanned'])} manuscript source files")
    print(f"  audit: {'PASS' if result['pass'] else 'FAIL'}")
    print("wrote provenance/wp10_input_manifest.json")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
