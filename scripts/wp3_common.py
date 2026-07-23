#!/usr/bin/env python3
"""Shared constants/helpers for WP3 (per-star extinction & de-reddened CMDs)."""
import numpy as np

# --- distance posterior (single-population; WP2 membership manifest) ---
D_KPC = 1.6245          # provenance/wp2_membership_manifest.json one_component_mean_kpc
D_KPC_ERR = 0.045       # posterior sigma (task brief)
DIST_MODULUS = 5.0 * np.log10(D_KPC * 1e3) - 5.0        # ~11.05 mag
DIST_MODULUS_ERR = 5.0 / np.log(10) * (D_KPC_ERR / D_KPC)  # propagate sigma

# --- Gaia DR3 VEGAMAG zeropoints (Gaia DR3 documentation, G/BP/RP) ---
GAIA_ZP = {"G": 25.6874, "BP": 25.3385, "RP": 24.7479}
GAIA_MAG_FLOOR = 0.003  # calibration systematic floor added in quadrature (mag)
TMASS_MAG_FLOOR = 0.02  # 2MASS systematic floor

BANDS = ["G", "BP", "RP", "J", "H", "Ks"]

def gaia_mag_error(mag, flux_error, band):
    """Convert Gaia flux error to magnitude error via the band zeropoint.

    flux = 10**(-(mag - ZP)/2.5);  sigma_mag = 1.0857 * flux_error / flux.
    A calibration floor is added in quadrature.
    """
    mag = np.asarray(mag, float)
    flux_error = np.asarray(flux_error, float)
    zp = GAIA_ZP[band]
    flux = 10.0 ** (-(mag - zp) / 2.5)
    sig = 1.0857 * flux_error / flux
    return np.sqrt(sig ** 2 + GAIA_MAG_FLOOR ** 2)

def normalize_teff(teff):
    """Anchor teff column mixes K and kK (e.g. 33.6 meaning 33600 K). Values
    below 200 are interpreted as kilo-Kelvin."""
    t = np.asarray(teff, float)
    out = np.where(np.isfinite(t) & (t < 200), t * 1000.0, t)
    return out

# Bright-star parallax override: use individual (zero-point corrected) parallax
# instead of the group posterior when the parallax is well measured.
BRIGHT_G_MAX = 11.0        # "brightest members" threshold in apparent G
PARALLAX_SNR_MIN = 10.0    # and parallax S/N high enough to be defensible
