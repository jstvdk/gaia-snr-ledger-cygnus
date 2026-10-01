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
ADOPTED = "repair_v8"     # adopted 2026-10-01: I1-I4 passed (issue19_repair_v8_integrity.json)

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
}

CHAIN = os.environ.get("CYGOB2_CHAIN", ADOPTED)
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
