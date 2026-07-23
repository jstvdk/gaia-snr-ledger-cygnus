#!/usr/bin/env python3
"""WP3 extinction-law utilities: R_V-dependent band coefficients k_X = A_X / A_V.

Uses the Cardelli, Clayton & Mathis (1989) parametrization with the O'Donnell
(1994) optical update -- the same curve the PARSEC/CMD service uses for its
"constant" Av coefficients, so our reddening is consistent with the isochrone
service. Monochromatic evaluation at band effective wavelengths.

Branches: R_V = 3.0, 3.1 (baseline), 3.5  (per CUTS_AND_THRESHOLDS.md #17, WP3 plan).

Effective wavelengths (micron):
  Gaia EDR3 G / BP / RP : from the CMD 3.9 service (G2V), 0.63902 / 0.51826 / 0.78251
  2MASS J / H / Ks      : Cohen et al. (2003) isophotal, 1.235 / 1.662 / 2.159

Caveat (documented, not hidden): band coefficients are mildly SED-dependent for
broad optical bands (the Teff/extinction degeneracy). This is exactly why WP3
derives A_V for hot spectroscopic anchors from intrinsic colours rather than
broadband fitting. For the cooler broadband-fit majority the monochromatic
approximation at the G2V effective wavelength is adequate at the ~few-percent
level and is bracketed by the R_V branches.
"""
import numpy as np

# band effective wavelengths in micron
LAMBDA_UM = {
    "G":  0.639021,
    "BP": 0.518258,
    "RP": 0.782508,
    "J":  1.235,
    "H":  1.662,
    "Ks": 2.159,
}

def ccm_odonnell(wave_um, R_V=3.1):
    """A_lambda / A_V from CCM89 with O'Donnell (1994) optical coefficients.

    Valid across the optical-NIR (0.3-3.3 um^-1 in x = 1/lambda). Returns
    A_lambda/A_V for scalar or array wave_um.
    """
    x = 1.0 / np.atleast_1d(np.asarray(wave_um, dtype=float))  # inverse microns
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # Infrared: 0.3 <= x < 1.1
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61

    # Optical/NIR: 1.1 <= x < 3.3 (O'Donnell 1994 coefficients)
    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    a[opt] = (1 + 0.104 * y - 0.609 * y**2 + 0.701 * y**3 + 1.137 * y**4
              - 1.718 * y**5 - 0.827 * y**6 + 1.647 * y**7 - 0.505 * y**8)
    b[opt] = (1.952 * y + 2.908 * y**2 - 3.989 * y**3 - 7.985 * y**4
              + 11.102 * y**5 + 5.491 * y**6 - 10.805 * y**7 + 3.347 * y**8)

    # Near-UV: 3.3 <= x < 8 (CCM89), for completeness (BP edge only ~1.93)
    uv = (x >= 3.3) & (x < 8.0)
    xu = x[uv]
    fa = np.where(xu >= 5.9, -0.04473 * (xu - 5.9) ** 2 - 0.009779 * (xu - 5.9) ** 3, 0.0)
    fb = np.where(xu >= 5.9, 0.2130 * (xu - 5.9) ** 2 + 0.1207 * (xu - 5.9) ** 3, 0.0)
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb

    alav = a + b / R_V
    return alav if np.ndim(wave_um) else float(alav[0])

def band_coefficients(R_V=3.1):
    """Return dict band -> k_X = A_X/A_V for the six bands at this R_V."""
    return {b: ccm_odonnell(lam, R_V) for b, lam in LAMBDA_UM.items()}

R_V_BRANCHES = [3.0, 3.1, 3.5]

def coefficient_table():
    return {f"{rv:.1f}": band_coefficients(rv) for rv in R_V_BRANCHES}

if __name__ == "__main__":
    import json
    tab = coefficient_table()
    # cross-check against CMD 3.9 service (G2V, CCM89+O'Donnell, Rv=3.1):
    #   A_G/A_V = 0.83627, A_BP/A_V = 1.08337, A_RP/A_V = 0.63439
    ref = {"G": 0.83627, "BP": 1.08337, "RP": 0.63439}
    print("Band coefficients k_X = A_X/A_V:")
    for rv, d in tab.items():
        print(f"  R_V={rv}: " + ", ".join(f"{b}={v:.4f}" for b, v in d.items()))
    print("\nCross-check vs CMD 3.9 (Rv=3.1, G2V):")
    for b, r in ref.items():
        mine = tab["3.1"][b]
        print(f"  {b}: mine={mine:.4f}  CMD={r:.4f}  diff={mine-r:+.4f} ({100*(mine-r)/r:+.1f}%)")
    print(json.dumps(tab, indent=2))
