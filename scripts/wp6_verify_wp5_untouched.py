#!/usr/bin/env python3
"""WP6 verification V1: the mass extension must not have moved WP5.

provenance/wp6_mass_extension_decision.json predicts under D1 that extending
the injection response above 18 Msun has NO effect on the accepted repair_v6
WP5 result, because wp5_common.MASS_GRID is unchanged and the extension is
written to a separate file prefix.  D1 is checkable, so it is checked rather
than asserted: every accepted repair_v6 artifact is re-hashed against the hash
the fit's own provenance recorded when it ran.

Output: provenance/wp6_verify_wp5_untouched.json

Run:
  PYTHONPATH=scripts python3 scripts/wp6_verify_wp5_untouched.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import wp5_common as w

ACCEPTED = "repair_v6"


def main() -> None:
    fit = json.loads(
        (w.PROVENANCE / f"wp5_imf_fit_execution_{ACCEPTED}.json").read_text(
            encoding="utf-8"
        )
    )
    recorded = fit.get("outputs", {})
    if not recorded:
        raise RuntimeError(
            f"wp5_imf_fit_execution_{ACCEPTED}.json carries no outputs block to "
            "verify against"
        )

    def expected_hash(entry):
        """The outputs block stores {"sha256": ..., "bytes": ...} per file;
        older records stored the bare digest.  Accept both."""
        if isinstance(entry, dict):
            return entry.get("sha256"), entry.get("bytes")
        return entry, None

    checks = []
    for relative, entry in recorded.items():
        expected, expected_bytes = expected_hash(entry)
        path = w.ROOT / relative
        if not path.exists():
            checks.append(
                {"path": relative, "status": "MISSING", "match": False}
            )
            continue
        actual = w.sha256(path)
        actual_bytes = path.stat().st_size
        match = bool(actual == expected)
        checks.append(
            {
                "path": relative,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "expected_bytes": expected_bytes,
                "actual_bytes": actual_bytes,
                "match": match,
                "status": "UNCHANGED" if match else "CHANGED",
            }
        )

    # The gate record is the artifact acceptance actually rests on.
    gate_path = w.PROVENANCE / f"wp5_{ACCEPTED}_gate.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate_state = {
        "accepted": gate.get("accepted"),
        "downstream_wp6_authorized": gate.get("downstream_wp6_authorized"),
        "blocking_reason": gate.get("blocking_reason"),
        "G3_pass": gate.get("gate_G3", {}).get("G3_pass"),
        "A_or_C_regressions": gate.get("gate_G3", {}).get("A_or_C_regressions"),
        "branches_passing": gate.get("gate_G3", {}).get("branches_passing"),
    }
    gate_intact = bool(
        gate_state["accepted"] is True
        and gate_state["downstream_wp6_authorized"] is True
        and gate_state["G3_pass"] is True
    )

    all_match = all(entry["match"] for entry in checks)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_verify_wp5_untouched.py",
        "status": "SUCCESS" if (all_match and gate_intact) else "FAILED",
        "verification": "V1 of provenance/wp6_mass_extension_decision.json",
        "claim_under_test": (
            "extending the injection response above 18 Msun leaves every "
            "accepted repair_v6 artifact byte-identical"
        ),
        "environment": {"python": sys.version, "platform": platform.platform()},
        "artifacts_checked": len(checks),
        "artifacts_unchanged": sum(1 for entry in checks if entry["match"]),
        "checks": checks,
        "gate_record": gate_state,
        "gate_intact": gate_intact,
        "V1_pass": bool(all_match and gate_intact),
    }
    w.write_json(w.PROVENANCE / "wp6_verify_wp5_untouched.json", record)

    for entry in checks:
        print(f"  {entry['status']:10s} {entry['path']}")
    print(f"\n  gate: accepted={gate_state['accepted']} "
          f"wp6_authorized={gate_state['downstream_wp6_authorized']} "
          f"G3={gate_state['G3_pass']} "
          f"regressions={gate_state['A_or_C_regressions']}")
    print(f"\nV1 {'PASS' if record['V1_pass'] else 'FAIL'} — "
          f"{record['artifacts_unchanged']}/{record['artifacts_checked']} "
          "artifacts byte-identical")
    print("wrote provenance/wp6_verify_wp5_untouched.json")


if __name__ == "__main__":
    main()
