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

import chain as C
import wp5_common as w

# Versions come from the active chain (scripts/chain.py).  repair_v7 is the
# chain whose WP6-WP12 products carry no suffix; on repair_v8 (issue #19) the
# same logical inputs resolve to their _repair_v8 siblings via C.tag.
WP5_VERSION = C.V["wp5"]
WP5_REPORT_VERSION = "repair_v6"   # the accepted gate record; v7 re-passed it
WP3_WP4_VERSION = C.V["wp4_ages"]
WP3_VERSION = C.V["wp3_extinction"]   # repair_v9: ages and extinction differ
WP4_MASS_VERSION = C.V["wp4_masses"]
T = C.tag_rel

# --------------------------------------------------------------- authorized
# Logical name -> path relative to the repository root.  Every WP5/WP6 entry is
# version-suffixed; nothing unversioned from those work packages appears here.
MANUSCRIPT_INPUTS: dict[str, str] = {
    # WP1-WP2 -- never repaired, so unversioned is correct for these
    "wp1_inventory": "tables/table1_wp1_inventory.md",
    "wp2_gate": "tables/table2_wp2_gate.md",
    "wp2_literature_recovery": "tables/table3_literature_recovery.md",
    "wp2_members": "data/processed/wp2_members.parquet",
    "wp2_control_members": "data/processed/wp2_control_members.parquet",
    "wp2_recovery_audit": "provenance/wp2_berlanas_recovery_audit.csv",
    "wp2_subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    # WP3 -- repaired
    "wp3_extinction": f"data/processed/wp3_extinction_{WP3_VERSION}.parquet",
    # WP4 -- repaired
    "wp4_age_posteriors": f"data/processed/wp4_age_posteriors_{WP3_WP4_VERSION}.parquet",
    "wp4_masses": f"data/processed/wp4_mass_posteriors_{WP4_MASS_VERSION}.parquet",
    # WP5 -- repaired twice; the normalization consumed downstream is v7
    "wp5_normalization": f"data/processed/wp5_imf_normalization_{WP5_VERSION}.parquet",
    "wp5_posterior_draws": f"data/processed/wp5_imf_posterior_draws_{WP5_VERSION}.npz",
    "wp5_association_mass": f"data/processed/wp5_association_mass_{WP5_VERSION}.parquet",
    "wp5_imf_norm_table": f"tables/wp5_imf_norm_{WP5_REPORT_VERSION}.csv",
    "wp5_baseline_residuals": f"tables/wp5_baseline_residuals_{WP5_REPORT_VERSION}.csv",
    "wp5_gate_record": f"provenance/wp5_{WP5_REPORT_VERSION}_gate.json",
    "wp5_association_mass_reconciliation": T("tables/wp5_association_mass_reconciliation.csv"),
    # WP6 -- closure re-run under repair_v7, and again under repair_v8
    "wp6_closure": f"tables/wp6_closure_{WP5_VERSION}.csv",
    "wp6_closure_attribution": f"tables/wp6_closure_attribution_{WP5_VERSION}.csv",
    "wp6_massive_census": T("tables/wp6_massive_census.csv"),
    "wp6_runaways": "tables/wp6_runaways.csv",
    "wp6_runaway_crossmatch": "tables/wp6_runaway_crossmatch.csv",
    "wp6_orphan_anchors": T("tables/wp6_orphan_anchors.csv"),
    "wp6_external_crosschecks": "tables/wp6_external_crosschecks.csv",
    # WP7-WP9 -- single versions, all computed on the repair_v7 chain
    "wp7_ledger": T("tables/wp7_ledger.csv"),
    "wp7_rsn_curves": T("tables/wp7_rsn_curves.csv"),
    "wp7_age_sensitivity": T("tables/wp7_age_sensitivity.csv"),
    "wp7_bh_threshold_scan": T("tables/wp7_bh_threshold_scan.csv"),
    "wp7_convergence": T("tables/wp7_convergence.csv"),
    "wp8_crosschecks": T("tables/wp8_crosschecks.csv"),
    "wp8_tension_list": T("tables/wp8_tension_list.csv"),
    "wp9_verdict": T("tables/wp9_verdict.csv"),
    "wp9_sensitivity": T("tables/wp9_sensitivity.csv"),
    # WP11 Part B -- post-hoc, pre-registered before scoring.  Registered here
    # deliberately: these are the ONLY post-WP1 comparisons in the manuscript
    # and the text must disclose that where it quotes them.
    "wp11_isotope_forecast": T("tables/wp11_isotope_forecast.csv"),
    "wp11_isotope_summary": T("tables/wp11_isotope_summary.csv"),
    # WP12 -- the manuscript-revision analysis.  Read-only over the frozen
    # repair_v7 chain; see provenance/wp12_revision_prereg.json.
    "wp12_gate_map": T("tables/wp12_wp5_gate_map.csv"),
    "wp12_combination_gate": T("tables/wp12_combination_gate.csv"),
    "wp12_branch_gate_table": T("tables/wp12_branch_gate_table.csv"),
    "wp12_closure_by_alpha": T("tables/wp12_closure_by_alpha.csv"),
    "wp12_closing_slopes": T("tables/wp12_closing_slopes.csv"),
    "wp12_mixed_slope_ledger": T("tables/wp12_mixed_slope_ledger.csv"),
    "wp12_c4_scan": T("tables/wp12_c4_scan.csv"),
    "wp12_c3_subtype": T("tables/wp12_c3_subtype.csv"),
    "wp12_scenario_score": T("tables/wp12_scenario_score.csv"),
    "wp12_neighbour_budget": T("tables/wp12_neighbour_budget.csv"),
    "wp12_cavity_share": T("tables/wp12_cavity_share.csv"),
    # pre-WP10 work
    "age_reconciliation": T("tables/wp4_wp5_age_reconciliation.csv"),
    "alpha_headline_branch_sets": T("tables/wp7_alpha_headline_branch_sets.csv"),
    "binary_bound": T("tables/wp7_binary_bound.csv"),
    "binary_bound_branches": T("tables/wp7_binary_bound_branches.csv"),
    "binary_bound_harer_fig2": T("tables/wp7_binary_bound_harer_fig2.csv"),
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
    # issue #19 -- the pre-repair WP4 run of 2026-07-23
    "data/processed/wp4_age_posteriors.parquet": (
        "issue #19 -- unversioned pre-repair WP4 age posterior; source of the "
        "withdrawn 2.25-5.67 Myr 'two-indicator' envelope.  Use "
        "wp4_age_posteriors_repair_v5.parquet"
    ),
    "data/processed/wp4_anchor_hrd.parquet": (
        "issue #19 -- anchor masses read at the pre-repair ages and copied onto "
        "every R_V branch.  Use wp4_anchor_hrd_repair_v8.parquet"
    ),
    "tables/wp4_ages_envelope.md": (
        "issue #19 -- pre-repair WP4 envelope table (2.25-5.67 Myr, PMS rows)"
    ),
    "tables/wp4_ages_table.md": (
        "issue #19 -- pre-repair WP4 age table"
    ),
}

# Once a later chain is active, the repair_v7 products it supersedes may not be
# quoted (CLAUDE.md rule 4).  They stay on disk; only the manuscript is barred.
SUPERSEDED_BY_LATER_CHAIN = [
    "data/processed/wp4_mass_posteriors_repair_v5.parquet",
    "data/processed/wp5_imf_normalization_repair_v7.parquet",
    "data/processed/wp5_imf_posterior_draws_repair_v7.npz",
    "data/processed/wp5_association_mass_repair_v7.parquet",
    "tables/wp6_closure_repair_v7.csv",
    "tables/wp6_closure_attribution_repair_v7.csv",
    "tables/wp6_massive_census.csv",
    "tables/wp6_orphan_anchors.csv",
    "tables/wp7_ledger.csv",
    "tables/wp7_rsn_curves.csv",
    "tables/wp7_age_sensitivity.csv",
    "tables/wp7_bh_threshold_scan.csv",
    "tables/wp7_convergence.csv",
    "tables/wp7_alpha_headline_branch_sets.csv",
    "tables/wp7_binary_bound.csv",
    "tables/wp7_binary_bound_branches.csv",
    "tables/wp8_crosschecks.csv",
    "tables/wp8_tension_list.csv",
    "tables/wp9_verdict.csv",
    "tables/wp9_sensitivity.csv",
    "tables/wp11_isotope_forecast.csv",
    "tables/wp11_isotope_summary.csv",
    "tables/wp4_wp5_age_reconciliation.csv",
]
if C.CHAIN != C.LEGACY:
    for _rel in SUPERSEDED_BY_LATER_CHAIN:
        FORBIDDEN[_rel] = (
            f"superseded by {C.CHAIN} (issue #19: anchor masses at pre-repair "
            "ages); preserved on disk, not quotable"
        )

# repair_v9 (issue #21, adopted 2026-10-07): the repair_v8 products -- built on
# the repair_v5 photometric ages -- are superseded in turn, with the repair_v5
# age table itself.  Generated from the same logical inputs, so nothing quoted
# by the manuscript on repair_v8 can be quoted on repair_v9.
SUPERSEDED_BY_REPAIR_V9 = [
    "data/processed/wp4_age_posteriors_repair_v5.parquet",
    "data/processed/wp4_mass_posteriors_repair_v8.parquet",
    "data/processed/wp4_anchor_hrd_repair_v8.parquet",
    "data/processed/wp5_imf_normalization_repair_v8.parquet",
    "data/processed/wp5_imf_posterior_draws_repair_v8.npz",
    "data/processed/wp5_association_mass_repair_v8.parquet",
    "tables/wp6_closure_repair_v8.csv",
    "tables/wp6_closure_attribution_repair_v8.csv",
] + [
    rel.replace(f"_{C.CHAIN}.", "_repair_v8.")
    for rel in MANUSCRIPT_INPUTS.values()
    if f"_{C.CHAIN}." in rel
]
if C.CHAIN == "repair_v9":
    for _rel in SUPERSEDED_BY_REPAIR_V9:
        FORBIDDEN[_rel] = (
            "superseded by repair_v9 (issue #21: ages from the Phase A' test-c "
            "spectroscopic HRD, not the repair_v5 photometric fit); preserved on "
            "disk, not quotable"
        )

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
