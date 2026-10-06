#!/usr/bin/env python3
"""Record the project owner's decision of 2026-10-06 on issue #21 (Phase A' report §4).

Written once; refuses to overwrite.  The decision is stored with the evidence
it rests on (paths and SHA-256), what it changes, and what it does NOT change.
Nothing is recomputed here.

Run:
  PYTHONPATH=scripts python3 scripts/decisions_2026_10_06.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import wp5_common as w

OUT = w.PROVENANCE / "decisions_2026_10_06.json"

EVIDENCE = [
    "reports/issue20b_phase_a_prime_report.md",
    "provenance/issue20b_outcome.json",
    "tables/issue20b_adopted_ages.csv",
    "tables/issue21_error_model_diagnostic.csv",
    "provenance/issue21_error_model_diagnostic_execution.json",
]

DECISIONS = {
    "5_issue21": {
        "question": "Issue #21: WP4 photometric ages rest on an independent-error approximation for "
                    "A_V (Phase A' report §4). Options: (a) follow the registered two-branch fallback; "
                    "(b) pre-register a WP4 age redesign first; (c) adopt spectroscopic HRD ages now.",
        "decision": "Option (b). Pre-register and run a WP4 age redesign (correct A_V error "
                    "propagation, a bright upper-MS indicator window, a separate study of the faint-end "
                    "mis-modelling, acceptance against the spectroscopic HRD), then rebuild the chain as "
                    "repair_v9 and only then run WP13.",
        "brief": "tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md",
        "changes": "work order only; no number moves until a pre-registered repair_v9 is adopted",
        "does_not_change": "repair_v8 stays the quoted chain until repair_v9 passes its adoption rule",
    },
}


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; a decision record is written once")
    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/decisions_2026_10_06.py",
        "decided_by": "project owner, 2026-10-06 ('Lets go b)')",
        "decisions": DECISIONS,
        "evidence": {p: w.sha256(w.ROOT / p) for p in EVIDENCE},
    })
    print(f"wrote {OUT.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
