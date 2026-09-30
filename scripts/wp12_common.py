#!/usr/bin/env python3
"""WP12 -- shared loaders and the read-only guard.

Every WP12 script goes through `frozen()` to open an input.  It checks the file
against the SHA-256 recorded in the preregistration and raises if it moved, so a
WP12 result can never be quietly computed from a re-run upstream product.
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import wp5_common as w

PREREG_PATH = w.PROVENANCE / "wp12_revision_prereg.json"

BASE = dict(family="PARSEC", R_V=3.1, alpha=2.3, sf_duration_Myr=0.0)
HEADLINE_ALPHAS = (2.0, 2.3)
SUPPORT_THRESHOLD = 0.5


class FrozenInputMoved(RuntimeError):
    """The artifact on disk is not the one WP12 preregistered."""


def prereg() -> dict:
    if not PREREG_PATH.exists():
        raise FileNotFoundError(
            "provenance/wp12_revision_prereg.json is missing.  Run "
            "scripts/wp12_prereg.py first -- WP12 results are only meaningful "
            "against a preregistration that predates them."
        )
    return json.loads(PREREG_PATH.read_text())


_PREREG: dict | None = None


def frozen(name: str) -> Path:
    """Resolve a preregistered input, verifying its hash."""
    global _PREREG
    if _PREREG is None:
        _PREREG = prereg()
    inputs = _PREREG["frozen_inputs"]
    if name not in inputs:
        raise KeyError(
            f"{name!r} is not a preregistered WP12 input.  WP12 may only read "
            "what wp12_prereg.py declared and hashed."
        )
    entry = inputs[name]
    path = w.ROOT / entry["path"]
    if not path.exists():
        raise FileNotFoundError(f"preregistered input is gone: {entry['path']}")
    digest = w.sha256(path)
    if digest != entry["sha256"]:
        raise FrozenInputMoved(
            f"{entry['path']} has changed since the WP12 preregistration.\n"
            f"  preregistered: {entry['sha256']}\n"
            f"  on disk      : {digest}\n"
            "WP12 consumes the frozen chain read-only.  If the upstream product "
            "genuinely changed, that is a new analysis and needs a new "
            "preregistration -- not a silent recomputation."
        )
    return path


def normalization() -> pd.DataFrame:
    """The repair_v7 per-cell mass-function fit, the version the ledger used."""
    return pd.read_parquet(frozen("wp5_normalization"))


def gate_map() -> pd.DataFrame:
    """One row per (subgroup, family, R_V, alpha) with its gate verdict."""
    columns = [
        "subgroup", "family", "R_V", "alpha", "k_median", "k_lo68", "k_hi68",
        "truth_age_posterior_mean_Myr", "poisson_chi_square_p",
        "residual_trend_p", "max_abs_pearson_residual", "residual_gate_pass",
    ]
    frame = normalization()[columns].copy()
    combo = ["family", "R_V", "alpha"]
    frame["all_subgroup_pass"] = frame.groupby(combo).residual_gate_pass.transform("all")
    frame["in_headline_set"] = frame.alpha.isin(HEADLINE_ALPHAS)
    return frame.sort_values(["family", "R_V", "alpha", "subgroup"]).reset_index(drop=True)


def combination_pass() -> pd.DataFrame:
    """One row per (family, R_V, alpha): does it pass in all three subgroups."""
    g = gate_map().groupby(["family", "R_V", "alpha"], as_index=False).agg(
        cells_passing=("residual_gate_pass", "sum"),
        all_subgroup_pass=("residual_gate_pass", "all"),
    )
    g["in_headline_set"] = g.alpha.isin(HEADLINE_ALPHAS)
    return g


def closure() -> pd.DataFrame:
    return pd.read_csv(frozen("wp6_closure"))


def ledger() -> pd.DataFrame:
    return pd.read_csv(frozen("wp7_ledger"))


def verdict() -> pd.DataFrame:
    return pd.read_csv(frozen("wp9_verdict"))


def branch_sets() -> pd.DataFrame:
    return pd.read_csv(frozen("wp7_branch_sets"))


def association_all_explode() -> pd.DataFrame:
    """The 54 all-explode association-scope ledger branches."""
    frame = ledger()
    return frame[
        frame.scope.eq("association") & frame.explodability.eq("all_explode")
    ].reset_index(drop=True)


def record(script: str, payload: dict, outputs: dict[str, Path]) -> dict:
    """Standard WP12 execution record: environment, prereg link, output hashes."""
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": script,
        "work_package": "WP12 -- manuscript revision analysis",
        "preregistration": "provenance/wp12_revision_prereg.json",
        "preregistration_sha256": w.sha256(PREREG_PATH),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        **payload,
        "outputs": {
            str(Path(p).relative_to(w.ROOT)): w.sha256(p)
            for p in outputs.values()
        },
    }


def score_prediction(pid: str, passed: bool, measured: dict) -> dict:
    spec = {p["id"]: p for p in prereg()["predictions"]}[pid]
    return {
        "id": pid,
        "statement": spec["statement"],
        "pass_if": spec["pass_if"],
        "outcome": "PASS" if passed else "FAIL",
        "measured": measured,
    }
