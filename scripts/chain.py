#!/usr/bin/env python3
"""Which chain a downstream script reads and writes (issue #19).

repair_v7 wrote many WP6-WP12 products without a version suffix
(tables/wp7_ledger.csv, tables/wp6_massive_census.csv, ...).  Rule 3 forbids
overwriting them, so a later chain writes suffixed siblings instead:
``tag(path)`` is the identity on repair_v7 and inserts ``_<chain>`` before the
extension otherwise.

The chain is selected by the environment variable CYGOB2_CHAIN; unset, it is
ADOPTED -- the chain the manuscript quotes.  Each chain declares the versions
of the upstream products it consumes, so no downstream script hard-codes them.

    CYGOB2_CHAIN=repair_v8 PYTHONPATH=scripts python3 scripts/wp7_ledger.py
"""
from __future__ import annotations

import os
from pathlib import Path

LEGACY = "repair_v7"      # the chain whose products carry no suffix
ADOPTED = "repair_v9"     # adopted 2026-10-07: I1-I5 passed (repair_v9_integrity.json)
# previously repair_v8, adopted 2026-10-01: I1-I4 passed (issue19_repair_v8_integrity.json)

VERSIONS: dict[str, dict[str, str]] = {
    # repair_v7: WP3/WP4 from repair_v5, WP5/WP6 injections and fit at v7
    "repair_v7": {
        "wp3_extinction": "repair_v5",
        "wp4_ages": "repair_v5",
        "wp4_masses": "repair_v5",
        "wp5": "repair_v7",
        "responses": "repair_v7",
    },
    # repair_v8 (issue #19): anchor masses at the repair_v5 ages, per R_V.
    # Injection responses are reused from repair_v7 -- they never read the
    # anchor masses (provenance/issue19_repair_v8_prereg.json, "reuse").
    "repair_v8": {
        "wp3_extinction": "repair_v5",
        "wp4_ages": "repair_v5",
        "wp4_masses": "repair_v8",
        "wp5": "repair_v8",
        "responses": "repair_v7",
    },
    # repair_v9 (issue #21, decisions of 2026-10-07): the WP4 ages become the
    # headline table -- A and C their Phase A' test-c spectroscopic posteriors,
    # B the A + C product (not measured for B).  Everything else is unchanged;
    # the WP5 node and WP6 extension responses are regenerated at the new ages
    # (provenance/repair_v9_prereg.json).
    "repair_v9": {
        "wp3_extinction": "repair_v5",
        "wp4_ages": "repair_v9_headline",
        "wp4_masses": "repair_v9",
        "wp5": "repair_v9",
        "responses": "repair_v9",
    },
    # Integrity check I1: the repair_v9 code path fed the repair_v5 ages and the
    # repair_v7 responses must reproduce repair_v8.
    "repair_v9_replay": {
        "wp3_extinction": "repair_v5",
        "wp4_ages": "repair_v5",
        "wp4_masses": "repair_v9_replay",
        "wp5": "repair_v9_replay",
        "responses": "repair_v7",
    },
}

# repair_v9 age scan (decision 7 of 2026-10-07): eleven mini-chains, one per
# native isochrone age, every subgroup forced to that age.  Index i is the
# i-th native age of each family (PARSEC 2.00 ... 6.31, MIST 2.00 ... 6.37).
SCAN_POINTS = 11
for _i in range(SCAN_POINTS):
    _scan = f"repair_v9_scan{_i:02d}"
    VERSIONS[_scan] = {
        "wp3_extinction": "repair_v5",
        "wp4_ages": _scan,
        "wp4_masses": _scan,
        "wp5": _scan,
        "responses": _scan,
    }

CHAIN =os.environ.get("CYGOB2_CHAIN", ADOPTED)
if CHAIN not in VERSIONS:
    raise RuntimeError(f"unknown CYGOB2_CHAIN={CHAIN!r}; known: {sorted(VERSIONS)}")
V = VERSIONS[CHAIN]


def tag(path: Path | str) -> Path:
    """This chain's copy of a product that repair_v7 wrote unversioned."""
    path = Path(path)
    if CHAIN == LEGACY:
        return path
    return path.with_name(f"{path.stem}_{CHAIN}{path.suffix}")


def tag_rel(relative: str) -> str:
    """``tag`` for a repository-relative path string."""
    return str(tag(Path(relative)))
