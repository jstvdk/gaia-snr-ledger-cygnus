# WP1 supernova markers

Frozen: 2026-07-22T10:29:58.255183+00:00

## Pulsars: ATNF PSRCAT

- Frozen release: ATNF Pulsar Catalogue v2.8.1; the exact public archive and database are stored under `data/raw/markers/`.
- Wide box l=72-88 deg, b=-5-8 deg: 80 pulsars in `data/processed/wp1_atnf_pulsars_wide.parquet`.
- PSR J2032+4127 / MT91 213 anchor:
  - P = 0.143246466289 s; P-dot = 1.13063e-14; characteristic age P/(2 P-dot) = 200.7 kyr.
  - Proper motion: mu_RA* = -2.99 mas/yr, mu_Dec = -0.74 mas/yr; ATNF parallax = 0.69 mas.
  - Association distance field = 1.33 kpc; catalogue association string: `GRS:4FGL_J2032.2+4127[aab+22],OPT:MT91_213(Cygnus_OB2)[lsk+15],XRS:[hnl+17]`.
  - The characteristic age is not an explosion-age measurement; braking index, birth period, and binary timing systematics must be carried in WP8.

## Supernova remnants

- Primary census: Green's 2024 October catalogue, VizieR VII/297 (310 confirmed Galactic SNRs). The wide box contains 9 entries, frozen in `data/processed/wp1_green_snrs_wide.parquet`.
- Gamma Cygni is Green G078.2+02.1 (DR4, gamma Cygni SNR), type `S`, angular size 60.0 arcmin.
- Physical anchor for Gamma Cygni: Leahy, Green & Ranasinghe (2013, MNRAS 436, 968; bibcode `2013MNRAS.436..968L`) infer d=1.7-2.6 kpc from H I absorption and a Sedov age of 6.8-10 kyr. Keep the older approximately 1.5 kpc class of estimates as an explicit literature branch; association with Cyg OB2 is unsettled.
- Chandra SNRcat snapshot: 6 entries in the same wide box. The Ferrand/Safi-Harb Manitoba SNRcat endpoint was not usable (the server returned its database-query error) on the freeze date; its raw response is retained rather than silently replacing missing high-energy fields.

## INTEGRAL radioactive-nuclide marker

- Martin et al. (2009/2010), A&A 506, 703 (arXiv:1001.1521; DOI `10.1051/0004-6361/200912178`): total Cygnus-region 1809 keV flux = (6.0 +/- 1.0) x 10^-5 ph cm^-2 s^-1; component attributed to the Cygnus complex = (3.9 +/- 1.1) x 10^-5 ph cm^-2 s^-1; 60Fe upper limit = 1.6 x 10^-5 ph cm^-2 s^-1. The inferred 26Al morphology is centred near Cyg OB2 and extends about 9 deg or more.
- Martin et al. (2010), A&A 511, A86 (arXiv:1001.1522; DOI `10.1051/0004-6361/200913385`) supplies the population-synthesis comparison. Both exact PDFs are frozen under `data/raw/markers/`.
- WP8 use: 26Al is a combined winds-plus-SN consistency constraint with a mean lifetime near 1 Myr, not a direct count of recent SNe and not a Cyg OB2-only measurement.
