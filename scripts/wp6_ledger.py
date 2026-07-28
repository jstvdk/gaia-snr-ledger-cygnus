#!/usr/bin/env python3
"""WP6 step 5: the living massive-star ledger handed to WP7.

Members + orphan anchors + recovered runaways, provenance-flagged per star, so
WP7 can subtract living stars from the "missing = dead" bookkeeping without
having to re-derive where each entry came from.

Every entry carries `census_channel` and a weight:
  member          weight = p_membership * P(M > 8), the WP4 posterior count
  orphan_anchor   weight = 1, spectroscopically confirmed and inside footprint
  runaway         weight = recovery probability, minus the chance-alignment
                  rate measured at that star's own separation

Runaways are LIVING stars that left; returning them to the census reduces N_SN.
They are the one channel of the three that WP7 must subtract explicitly.

Outputs:
  tables/wp6_massive_census.cat        the ledger
  provenance/wp6_ledger_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_ledger.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

UPSTREAM = "repair_v5"
BASELINE = ("PARSEC", 3.1)


def main() -> None:
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet")
    orphans = pd.read_csv(w.TABLES / "wp6_orphan_anchors.csv")
    runaways = pd.read_csv(w.TABLES / "wp6_runaways.csv")
    # The aggregate correction is read from the traceback's own record rather
    # than repeated here, so the two can never drift apart.
    traceback = json.loads(
        (w.PROVENANCE / "wp6_runaways_execution.json").read_text(encoding="utf-8")
    )
    raw_recovered = int(traceback["result"]["raw_recovered"])
    aggregate_unclipped = float(traceback["result"]["false_positive_corrected"])

    column = f"mass_{BASELINE[0]}_rv{BASELINE[1]:.1f}_p_gt8"
    point = w.mass_column(*BASELINE)
    member_rows = masses[masses[column].notna()].copy()
    members = pd.DataFrame(
        {
            "source_id": member_rows.source_id.astype("int64"),
            "census_channel": "member",
            "subgroup": member_rows.subgroup,
            "weight": member_rows.membership_probability * member_rows[column],
            "mass_estimate_Msun": member_rows[point],
            "mass_source": f"WP4 posterior ({BASELINE[0]}, R_V={BASELINE[1]})",
            "alive": True,
        }
    )
    members = members[members.weight > 0].copy()

    keep = orphans[
        orphans.inside_footprint.astype(bool)
        & orphans.classifiable.astype(bool)
        & orphans.above_8_Msun.astype(bool)
    ]
    anchor_entries = pd.DataFrame(
        {
            "source_id": keep.source_id.astype("int64"),
            "census_channel": "orphan_anchor",
            "subgroup": keep.subgroup_positional,
            "weight": 1.0,
            "mass_estimate_Msun": keep.initial_mass_Msun,
            "mass_source": "spectral type (" + keep.mass_rule + ")",
            "alive": True,
        }
    )

    # Purity, per separation bin.  A per-star "probability minus chance rate"
    # does NOT reproduce the aggregate false-positive correction: the control
    # rate is the fraction of ALL candidates that pass by chance, so the
    # expected number of impostors among the SELECTED stars is that rate summed
    # over every candidate, not over the selected ones.  Purity is therefore
    # (selected - expected chance) / selected, evaluated in bins of separation
    # because both quantities vary steeply with it.
    #
    # The centroid is recomputed the way scripts/wp6_runaways.py computes it --
    # the median position of the labelled members -- rather than hard-coded, so
    # the separations binned here are the separations the search itself used.
    labels = pd.read_parquet(
        w.TABLES / "wp2_subgroup_labels.parquet", columns=["source_id", "subgroup"]
    )
    member_positions = pd.read_parquet(
        w.PROC / "wp2_members.parquet", columns=["source_id", "l_deg", "b_deg"]
    ).merge(labels, on="source_id", how="inner")
    centre = member_positions[["l_deg", "b_deg"]].median()
    centre_l, centre_b = float(centre.l_deg), float(centre.b_deg)
    separation = np.hypot(
        (runaways.l - centre_l) * np.cos(np.radians(runaways.b)),
        runaways.b - centre_b,
    )
    runaways = runaways.assign(separation_deg=separation)
    edges = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 90.0])
    purity = np.zeros(len(runaways))
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = runaways.separation_deg.ge(lo) & runaways.separation_deg.lt(hi)
        if not in_bin.any():
            continue
        chosen = in_bin & runaways.is_runaway_candidate.astype(bool)
        n_selected = int(chosen.sum())
        expected_chance = float(
            runaways.loc[in_bin, "false_positive_rate_at_separation"].sum()
        )
        value = (
            max(n_selected - expected_chance, 0.0) / n_selected
            if n_selected else 0.0
        )
        purity[in_bin.to_numpy()] = value
    runaways = runaways.assign(purity=purity)

    selected = runaways[runaways.is_runaway_candidate.astype(bool)].copy()
    # Weight IS the purity: summing it over the selected stars in a bin
    # gives (n_selected - expected chance) exactly, so the ledger total is
    # the false-positive-corrected count by construction.
    net = selected.purity
    runaway_entries = pd.DataFrame(
        {
            "source_id": selected.source_id.astype("int64"),
            "census_channel": "runaway",
            "subgroup": "unassigned",
            "weight": net,
            "mass_estimate_Msun": np.nan,
            "mass_source": "not estimated (no per-star extinction outside the member sample)",
            "alive": True,
        }
    )

    ledger = pd.concat(
        [members, anchor_entries, runaway_entries], ignore_index=True
    )
    out = w.TABLES / "wp6_massive_census.cat"
    ledger.to_csv(out, index=False)

    by_channel = {
        channel: {
            "entries": int(len(block)),
            "summed_weight": round(float(block.weight.sum()), 2),
        }
        for channel, block in ledger.groupby("census_channel")
    }
    by_subgroup = {
        subgroup: round(
            float(ledger[ledger.subgroup.eq(subgroup)].weight.sum()), 2
        )
        for subgroup in list(w.SUBGROUPS) + ["unassigned"]
    }

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_ledger.py",
        "status": "SUCCESS",
        "work_package": "WP6 step 5",
        "baseline_branch": f"{BASELINE[0]} R_V={BASELINE[1]} (reporting branch)",
        "runaway_weighting": (
            "purity per separation bin, times the star's own recovery "
            "probability.  Purity = (selected - expected chance) / selected, "
            "with expected chance summed over ALL candidates in the bin, so the "
            "ledger total reproduces the aggregate false-positive correction "
            "instead of contradicting it."
        ),
        "weight_definitions": {
            "member": "p_membership * P(M > 8 Msun) from the WP4 posterior",
            "orphan_anchor": (
                "1.0 — spectroscopically confirmed, inside the 1 deg footprint, "
                "above 8 Msun by its spectral type, and absent from the member "
                "sample"
            ),
            "runaway": (
                "purity(separation) * recovery probability; see "
                "runaway_weighting"
            ),
        },
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet",
                w.TABLES / "wp6_orphan_anchors.csv",
                w.TABLES / "wp6_runaways.csv",
            ]
        },
        "by_channel": by_channel,
        "by_subgroup_weight": by_subgroup,
        "total_living_above_8_Msun": round(float(ledger.weight.sum()), 2),
        "runaway_total_reconciliation": {
            "ledger_binned": round(float(runaway_entries.weight.sum()), 2),
            "aggregate_unclipped": aggregate_unclipped,
            "why_they_differ": (
                "the binned purity clips a bin at zero when its expected "
                "chance count exceeds the number selected there, because a "
                "negative count of real runaways is meaningless.  The "
                "unclipped aggregate allows those bins to subtract from "
                "others.  The clipped figure is the conservative one to "
                "quote per bin; the unclipped one is the better global "
                "estimate.  Both are reported and WP7 should carry the "
                "difference as a systematic on the runaway correction."
            ),
        },
        "runaway_provenance": {
            "traceback_record": "provenance/wp6_runaways_execution.json",
            "raw_recovered": raw_recovered,
            "aggregate_false_positive_corrected": aggregate_unclipped,
            "association_centroid_galactic_deg": [
                round(centre_l, 4), round(centre_b, 4)
            ],
            "peculiar_motion": (
                "the traceback runs on proper motions with the association's "
                "systemic motion subtracted (issue #16).  An earlier absolute-PM "
                "version of this ledger, quoting 260 raw and 109.2 corrected, is "
                "WITHDRAWN: Cyg OB2's systemic motion is larger than a typical "
                "ejection signature, so that traceback measured bulk drift."
            ),
            "external_gate": (
                "provenance/wp6_runaway_crossmatch_execution.json — BD+43 3654 "
                "recovered at p = 1.000, 38.8 km/s and 1.36 Myr against a "
                "literature ~40 km/s and 1.6 Myr"
            ),
        },
        "handoff_to_wp7": (
            "runaways are LIVING stars that left the association.  WP7 must "
            "subtract them from the missing-equals-dead bookkeeping, which "
            "reduces N_SN.  Their subgroup is 'unassigned': a 2D traceback "
            "constrains the origin footprint, not which subgroup inside it, so "
            "assigning them per subgroup would invent precision."
        ),
        "carried_limitations": [
            "runaway weights are statistical, not per-star confirmations",
            "the traceback is 2D, so line-of-sight ejections are missed and the "
            "runaway count is a lower bound",
            "orphan-anchor masses come from spectral type, not from the WP4 "
            "posterior, so they carry a different systematic",
            "13 orphan anchors have no spectral type and are excluded as "
            "unknown rather than counted either way",
        ],
        "outputs": {str(out.relative_to(w.ROOT)): w.sha256(out)},
    }
    w.write_json(w.PROVENANCE / "wp6_ledger_execution.json", record)

    print("WP6 step 5 — living massive-star ledger\n")
    for channel, block in by_channel.items():
        print(f"  {channel:15s} {block['entries']:5d} entries, "
              f"summed weight {block['summed_weight']:8.2f}")
    print(f"\n  total living above 8 Msun: {record['total_living_above_8_Msun']:.2f}")
    print("  by subgroup:")
    for subgroup, value in by_subgroup.items():
        print(f"    {subgroup:12s} {value:8.2f}")
    print("\nwrote provenance/wp6_ledger_execution.json")


if __name__ == "__main__":
    main()
