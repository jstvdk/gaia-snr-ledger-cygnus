#!/usr/bin/env python3
"""Pre-registration: WP8, external cross-checks of the supernova ledger.

Written BEFORE any comparison is scored.  Predictions and decision rules are
fixed here and are not amended afterwards.

WHAT WP8 IS
-----------
The ledger is confronted with every independent supernova marker in the field.
Nothing in the chain is tuned toward any of them -- a disagreement becomes an
issue in PROJECT_TRACE section 9, never a reason to move a number.  The markers
were frozen at WP1 on 2026-07-22, long before the ledger existed.

THE CENTRAL TEST
----------------
WP7 left one question dominating everything: does the supernova budget exist at
all?  On the all-explode branch the association has produced about 8 supernovae;
on any branch with a black-hole threshold at or below 40 Msun it has produced
EXACTLY ZERO, because the smallest turnoff on the grid is about 52 Msun.

PSR J2032+4127 is a NEUTRON STAR.  Neutron stars are made by successful
explosions, not by direct collapse.  Its existence is therefore a direct
observational discriminant between two branches the Gaia data cannot separate.

AND THE COMPANION IS OUR OWN STAR
---------------------------------
MT91 213, the pulsar's Be/B0V companion, is anchor
gaia_dr3:2067835682818358400 in this project's own catalogue: a spectroscopic
orphan anchor, 17 Msun by the B0V mass rule, inside the footprint, positionally
assigned to CygOB2-A at 0.115 deg from its centroid, and already counted in the
living census.  The pulsar is therefore not an external object being compared
against -- it is a compact remnant sitting next to a star the ledger counts.

That makes the test three-way rather than binary, and the pre-registration says
so before scoring:

  (a) massive stars above ~58 Msun DO make neutron stars     -> all-explode
  (b) the progenitor was less massive, so OLDER than CygOB2-A -> an older
      population exists that WP4 does not see
  (c) binary mass transfer stripped the progenitor, so a lower-mass star died
      early -> the effect WP7 section 10 already names as its largest
      unmodelled physics

Output: provenance/wp8_crosschecks_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp8_crosschecks_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

# ---- frozen WP1 marker values, quoted not recomputed -----------------------
PULSAR = {
    "name": "PSR J2032+4127",
    "P_s": 0.14324646628873475,
    "Pdot": 1.13062721074195e-14,
    "characteristic_age_yr": 200738.4,
    "pmra_masyr": -2.99,
    "pmdec_masyr": -0.74,
    "companion": "MT91 213 (B0V)",
    "companion_source_id": 2067835682818358400,
    "source": "ATNF PSRCAT v2.8.1, frozen 2026-07-22",
}
GAMMA_CYGNI = {
    "name": "G078.2+02.1 (gamma Cygni, DR4)",
    "age_kyr": (6.8, 10.0),
    "distance_kpc": (1.7, 2.6),
    "older_literature_distance_kpc": 1.5,
    "separation_from_centroid_deg": 2.27,
    "source": "Green 2024 VII/297; Leahy, Green & Ranasinghe 2013",
}
AL26 = {
    "complex_flux_ph_cm2_s": 3.9e-5,
    "complex_flux_error": 1.1e-5,
    "mean_lifetime_Myr": 1.05,
    "source": "Martin et al. 2009, A&A 506, 703",
}

# Braking-index and birth-period systematics on the characteristic age.
# tau = P / ((n-1) Pdot) * [1 - (P0/P)^(n-1)];  n = 3 and P0 << P give tau_c.
BRAKING_INDEX_RANGE = (2.0, 3.0)
BIRTH_PERIOD_FRACTION_RANGE = (0.0, 0.5)

# SNR visible lifetime in a normal-density ISM.  Deliberately generous at the
# top end so the non-detection argument is not manufactured by a short lifetime.
SNR_VISIBLE_LIFETIME_KYR = (20.0, 100.0)

CENTROID_L_DEG, CENTROID_B_DEG = 80.073, 0.816


def main() -> None:
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp8_crosschecks_prereg.py",
        "status": "PREREGISTERED",
        "work_package": "WP8",
        "objective": (
            "confront the WP7 ledger with every independent supernova marker in "
            "the field, in order of constraining power"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "validations_never_calibrations": (
            "no value in the chain is tuned toward any marker.  The markers were "
            "frozen at WP1 on 2026-07-22, before WP5, WP6 and WP7 existed.  A "
            "disagreement becomes an issue in PROJECT_TRACE section 9, never a "
            "reason to move a number."
        ),
        "inputs_frozen_at_wp1": {
            "pulsars": "data/processed/wp1_atnf_pulsars_wide.parquet (80 in box)",
            "snrs": "data/processed/wp1_green_snrs_wide.parquet (9 in box)",
            "al26": "wp1_sn_markers.md, Martin et al. 2009/2010",
            "ledger": "tables/wp7_ledger.csv, provenance/wp7_ledger_execution.json",
        },
        "markers": {
            "pulsar": PULSAR,
            "gamma_cygni": GAMMA_CYGNI,
            "al26": AL26,
        },
        "the_companion_is_our_own_star": {
            "anchor_uid": "gaia_dr3:2067835682818358400",
            "object_name": "[MT91] 213",
            "spectral_type": "B0V",
            "initial_mass_Msun": 17.0,
            "channel": "orphan anchor (one of the 27), counts_in_census = True",
            "subgroup_positional": "CygOB2-A",
            "centroid_separation_deg": 0.115,
            "why_not_a_wp2_member": (
                "it is a spectroscopic anchor that the astrometric membership "
                "selection did not recover.  Its Gaia astrometry is clean "
                "(RUWE 1.07, parallax error 0.016 mas), so this is a selection "
                "outcome and not a data-quality failure; the orphan-anchor "
                "channel exists precisely for such stars.  Reported, because a "
                "reader will ask whether the pulsar's companion is in the "
                "census -- it is."
            ),
            "consequence_for_the_test": (
                "the pulsar is not an external object being compared against.  "
                "It is a compact remnant beside a star the ledger already "
                "counts, in the footprint of the subgroup that produces most of "
                "the ledger's supernovae."
            ),
        },
        "check_1_pulsar": {
            "rank": "first -- the only marker that can settle a WP7 branch",
            "existence_argument": (
                "a neutron star requires a SUCCESSFUL explosion.  Direct "
                "collapse to a black hole produces no neutron star.  The "
                "islands branch predicts exactly zero supernovae in Cyg OB2, so "
                "the pulsar's existence is incompatible with it UNLESS "
                "explanation (b) or (c) below holds."
            ),
            "three_way_reading": {
                "a_massive_stars_explode": (
                    "the progenitor was above CygOB2-A's 57.9 Msun turnoff and "
                    "still made a neutron star -- the all-explode branch"
                ),
                "b_older_population": (
                    "the progenitor was less massive and therefore older than "
                    "CygOB2-A, implying a population WP4 does not resolve"
                ),
                "c_binary_stripping": (
                    "mass transfer stripped the progenitor's envelope so a "
                    "lower-mass star died early and left a neutron star.  WP7 "
                    "section 10 already names binary mass transfer as its "
                    "largest unmodelled effect, and WP6 measured multiplicity "
                    "above 8 Msun at f_bin ~ 0.7, so this is not a special plea."
                ),
            },
            "age_systematics_carried": {
                "braking_index_range": list(BRAKING_INDEX_RANGE),
                "birth_period_fraction_range": list(BIRTH_PERIOD_FRACTION_RANGE),
                "formula": "tau = P / ((n-1) Pdot) * [1 - (P0/P)^(n-1)]",
                "note": (
                    "the characteristic age is NOT an explosion-age "
                    "measurement, as wp1_sn_markers.md already states.  n = 3 "
                    "with P0 << P recovers tau_c = 200.7 kyr; smaller n "
                    "lengthens the true age, larger P0 shortens it.  The "
                    "comparison is made against the RANGE, not the point value."
                ),
                "binary_contamination": (
                    "J2032+4127 is in a highly eccentric decades-long orbit, so "
                    "the measured Pdot carries line-of-sight acceleration from "
                    "the companion.  This is acknowledged as a systematic on "
                    "tau_c that this project cannot remove."
                ),
            },
        },
        "check_2_gamma_cygni": {
            "rank": "second -- a two-sided check, and genuinely ambiguous",
            "two_sided": (
                "(i) does the ledger comfortably allow a supernova 7-10 kyr "
                "ago?  (ii) does our Gaia distance support physical association "
                "with Cyg OB2 at all?  Both must be answered; neither is "
                "allowed to be overclaimed."
            ),
            "distance_tension_is_real": (
                "our association distance is 1.62 kpc; Leahy+2013 infer "
                "1.7-2.6 kpc for the remnant from H I absorption.  These "
                "overlap only at the extreme low end of theirs.  The older "
                "~1.5 kpc class of estimates is carried as an explicit "
                "literature branch.  Association is UNSETTLED and WP8 will not "
                "settle it."
            ),
            "geometry_note": (
                f"the remnant sits "
                f"{GAMMA_CYGNI['separation_from_centroid_deg']} deg from the "
                f"association centroid, about 64 pc in projection at 1.62 kpc.  "
                f"That is outside the 1 deg footprint but well inside the "
                f"distance a runaway covers in a few Myr, so a runaway "
                f"progenitor is a live possibility rather than a contrivance."
            ),
        },
        "check_3_al26": {
            "rank": "third -- a consistency band, explicitly not a fit",
            "why_it_cannot_be_a_count": (
                "26Al is produced by Wolf-Rayet winds as well as by supernovae, "
                "and the Martin et al. flux is for the whole Cygnus complex, "
                "not for Cyg OB2 alone.  It constrains a COMBINATION.  Any "
                "attempt to invert it into a supernova count would be "
                "overreach, and is pre-emptively forbidden here."
            ),
            "permitted_claim": (
                "an order-of-magnitude consistency statement only: whether the "
                "ledger's recent activity is compatible with the observed "
                "1809 keV flux given standard yields"
            ),
        },
        "check_4_snr_absence": {
            "rank": "fourth -- and expected to be WEAK evidence",
            "visible_lifetime_kyr": list(SNR_VISIBLE_LIFETIME_KYR),
            "the_honest_framing": (
                "the Haerer scenario predicts an invisible remnant in a "
                "low-density cavity, and no catalogued SNR sits at the "
                "association's position.  But the ledger's rate is about 8 per "
                "Myr, so even with a GENEROUS 100 kyr visible lifetime and no "
                "cavity at all the expected number of visible remnants is "
                "below one.  A non-detection is therefore consistent with the "
                "cavity argument AND with ordinary Poisson luck, and cannot "
                "discriminate between them.  Saying so is the result."
            ),
        },
        "check_5_neighbours": {
            "rank": "fifth -- deliberately coarse",
            "scope": (
                "a literature-based bound on the supernova budget of Cyg OB1, "
                "OB9 and the surrounding field, to bound the probability that a "
                "cavity supernova came from somewhere other than Cyg OB2.  "
                "Labelled COARSE; no pipeline is run on them."
            ),
        },
        "predictions": [
            {
                "id": "X1",
                "statement": (
                    "the ledger assigns P(at least one supernova) > 0.5 on the "
                    "all-explode baseline branch, so the pulsar's existence is "
                    "comfortably accommodated rather than being a surprise"
                ),
                "falsifies_if": "P(>=1 SN) <= 0.5 on the baseline branch",
                "if_falsified": (
                    "gate-level tension: the ledger would be saying Cyg OB2 "
                    "probably produced no supernovae while a neutron star sits "
                    "inside it, which would require a documented resolution "
                    "before anything is reported"
                ),
            },
            {
                "id": "X2",
                "statement": (
                    "the islands branch assigns EXACTLY zero probability to any "
                    "supernova, so the pulsar formally excludes it under "
                    "reading (a) -- leaving (b) an older population or (c) "
                    "binary stripping as the only alternatives"
                ),
                "falsifies_if": "the islands branch assigns non-zero probability",
                "if_falsified": (
                    "WP7's L2 result would be contradicted and the explodability "
                    "branch would have to be re-derived"
                ),
            },
            {
                "id": "X3",
                "statement": (
                    "the ledger assigns P(last supernova within the pulsar's "
                    "systematics-widened age range) > 0.2, i.e. the pulsar's "
                    "age is a typical draw from the ledger rather than a tail "
                    "event"
                ),
                "falsifies_if": "the probability is <= 0.2",
                "if_falsified": (
                    "recorded as a mild tension and attributed between the "
                    "characteristic-age systematics and the ledger's rate, "
                    "without moving either"
                ),
            },
            {
                "id": "X4",
                "statement": (
                    "the expected number of CURRENTLY VISIBLE supernova "
                    "remnants implied by the ledger is below 1 even at the "
                    "generous 100 kyr visible lifetime and with no cavity, so "
                    "the observed absence is weak evidence"
                ),
                "falsifies_if": "the expected number is 1 or more",
                "if_falsified": (
                    "the non-detection becomes genuinely constraining and the "
                    "cavity argument gains real support, which would be a "
                    "stronger result than expected and must be reported as such"
                ),
            },
            {
                "id": "X5",
                "statement": (
                    "P(a supernova within gamma Cygni's 6.8-10 kyr age window) "
                    "lies between 0.03 and 0.12 on the baseline branch -- "
                    "non-negligible, so allowed, but not expected"
                ),
                "falsifies_if": "the probability falls outside [0.03, 0.12]",
                "if_falsified": (
                    "the rate argument underlying the ledger's recent history "
                    "is behaving differently from the flat R_SN(t) WP7 reports"
                ),
            },
        ],
        "gate_criteria": {
            "G8a": (
                "pulsar consistency RESOLVED -- either agreement, or a "
                "documented tension with candidate explanations.  This is the "
                "plan's stated gate for WP8."
            ),
            "G8b": (
                "the gamma Cygni distance ambiguity is reported as unsettled "
                "and not resolved in our favour"
            ),
            "G8c": (
                "the 26Al comparison stays a consistency band and is never "
                "inverted into a supernova count"
            ),
            "G8d": "a tension list is produced, even if empty",
        },
        "what_wp8_cannot_do": [
            "it cannot establish that gamma Cygni is associated with Cyg OB2; "
            "the distance evidence is genuinely ambiguous and will stay so",
            "it cannot convert the 26Al flux into a supernova count, because "
            "Wolf-Rayet winds contribute and the measurement is complex-wide",
            "it cannot distinguish reading (a) from (c) for the pulsar -- "
            "binary stripping and high-mass explodability both produce a "
            "neutron star, and no observable here separates them",
            "it cannot remove the binary-orbit contamination of the pulsar's "
            "Pdot, which is a timing problem beyond this project's data",
        ],
    }
    w.write_json(w.PROVENANCE / "wp8_crosschecks_prereg.json", record)

    print("WP8 — external cross-checks — pre-registered\n")
    print("  central test: PSR J2032+4127 is a NEUTRON STAR, so a successful")
    print("                explosion happened.  The islands branch predicts zero.")
    print(f"\n  and the companion is our own star: [MT91] 213, B0V, 17 Msun,")
    print(f"  orphan anchor in CygOB2-A, already counted in the living census.")
    print("\n  predictions:")
    for entry in record["predictions"]:
        head = entry["statement"].split(",")[0]
        print(f"    {entry['id']}  {head[:64]}")
    print("\nwrote provenance/wp8_crosschecks_prereg.json")


if __name__ == "__main__":
    main()
