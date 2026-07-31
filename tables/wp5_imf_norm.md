> ⚠️ **SUPERSEDED — FROZEN PRE-REPAIR RECORD. DO NOT QUOTE, DO NOT PROPAGATE.**
> This is the original WP5 run in which **0 of 54 branches passed** the
> mass-function residual gate (every `residual_gate_pass` below is `False`).
> It is retained only as the historical record of issue #2. The accepted
> normalization is [`wp5_imf_norm_repair_v6.md`](wp5_imf_norm_repair_v6.md) and
> the version consumed downstream is
> `data/processed/wp5_imf_normalization_repair_v7.parquet`.
> Banner added 2026-07-30 (item S3); the table below is byte-unchanged.
> Enforcement: [`scripts/wp10_inputs.py`](../scripts/wp10_inputs.py) refuses
> this file as a manuscript input.

# WP5 IMF-normalization branch table

| subgroup | family | R_V | alpha | k_median | k_lo68 | k_hi68 | poisson_chi_square_p | residual_gate_pass |
|---|---|---|---|---|---|---|---|---|
| CygOB2-A | PARSEC | 3 | 2 | 1.03e+03 | 975 | 1.1e+03 | 1.98e-15 | False |
| CygOB2-B | PARSEC | 3 | 2 | 1.02e+03 | 958 | 1.09e+03 | 1.52e-23 | False |
| CygOB2-C | PARSEC | 3 | 2 | 1.09e+03 | 1.01e+03 | 1.16e+03 | 1.61e-13 | False |
| CygOB2-A | PARSEC | 3 | 2.3 | 1.46e+03 | 1.37e+03 | 1.55e+03 | 3.04e-11 | False |
| CygOB2-B | PARSEC | 3 | 2.3 | 1.44e+03 | 1.34e+03 | 1.54e+03 | 1.67e-17 | False |
| CygOB2-C | PARSEC | 3 | 2.3 | 1.54e+03 | 1.44e+03 | 1.66e+03 | 1.51e-10 | False |
| CygOB2-A | PARSEC | 3 | 2.6 | 2e+03 | 1.87e+03 | 2.13e+03 | 3.97e-09 | False |
| CygOB2-B | PARSEC | 3 | 2.6 | 1.97e+03 | 1.83e+03 | 2.12e+03 | 3.84e-13 | False |
| CygOB2-C | PARSEC | 3 | 2.6 | 2.13e+03 | 1.96e+03 | 2.3e+03 | 2.14e-09 | False |
| CygOB2-A | PARSEC | 3.1 | 2 | 1.09e+03 | 1.03e+03 | 1.16e+03 | 1.44e-17 | False |
| CygOB2-B | PARSEC | 3.1 | 2 | 1.04e+03 | 971 | 1.1e+03 | 9.35e-25 | False |
| CygOB2-C | PARSEC | 3.1 | 2 | 1.08e+03 | 1e+03 | 1.16e+03 | 1.3e-16 | False |
| CygOB2-A | PARSEC | 3.1 | 2.3 | 1.53e+03 | 1.44e+03 | 1.63e+03 | 1.27e-13 | False |
| CygOB2-B | PARSEC | 3.1 | 2.3 | 1.46e+03 | 1.37e+03 | 1.56e+03 | 1.79e-18 | False |
| CygOB2-C | PARSEC | 3.1 | 2.3 | 1.53e+03 | 1.42e+03 | 1.64e+03 | 1.42e-13 | False |
| CygOB2-A | PARSEC | 3.1 | 2.6 | 2.1e+03 | 1.97e+03 | 2.23e+03 | 4.21e-12 | False |
| CygOB2-B | PARSEC | 3.1 | 2.6 | 2e+03 | 1.86e+03 | 2.15e+03 | 5.81e-14 | False |
| CygOB2-C | PARSEC | 3.1 | 2.6 | 2.11e+03 | 1.95e+03 | 2.28e+03 | 1.84e-12 | False |
| CygOB2-A | PARSEC | 3.5 | 2 | 1.22e+03 | 1.15e+03 | 1.29e+03 | 1.08e-36 | False |
| CygOB2-B | PARSEC | 3.5 | 2 | 1.15e+03 | 1.08e+03 | 1.23e+03 | 9.37e-42 | False |
| CygOB2-C | PARSEC | 3.5 | 2 | 1.17e+03 | 1.08e+03 | 1.25e+03 | 1.05e-10 | False |
| CygOB2-A | PARSEC | 3.5 | 2.3 | 1.72e+03 | 1.62e+03 | 1.82e+03 | 1.64e-29 | False |
| CygOB2-B | PARSEC | 3.5 | 2.3 | 1.63e+03 | 1.52e+03 | 1.74e+03 | 1.48e-33 | False |
| CygOB2-C | PARSEC | 3.5 | 2.3 | 1.68e+03 | 1.55e+03 | 1.81e+03 | 3.38e-09 | False |
| CygOB2-A | PARSEC | 3.5 | 2.6 | 2.37e+03 | 2.22e+03 | 2.52e+03 | 1.85e-25 | False |
| CygOB2-B | PARSEC | 3.5 | 2.6 | 2.25e+03 | 2.08e+03 | 2.41e+03 | 7.55e-28 | False |
| CygOB2-C | PARSEC | 3.5 | 2.6 | 2.34e+03 | 2.15e+03 | 2.55e+03 | 2.78e-09 | False |
| CygOB2-A | MIST | 3 | 2 | 976 | 918 | 1.04e+03 | 5.57e-12 | False |
| CygOB2-B | MIST | 3 | 2 | 1.05e+03 | 979 | 1.11e+03 | 9.53e-25 | False |
| CygOB2-C | MIST | 3 | 2 | 1.1e+03 | 1.03e+03 | 1.18e+03 | 1.16e-11 | False |
| CygOB2-A | MIST | 3 | 2.3 | 1.37e+03 | 1.28e+03 | 1.45e+03 | 3.29e-08 | False |
| CygOB2-B | MIST | 3 | 2.3 | 1.48e+03 | 1.38e+03 | 1.58e+03 | 6.45e-19 | False |
| CygOB2-C | MIST | 3 | 2.3 | 1.58e+03 | 1.47e+03 | 1.69e+03 | 1.09e-08 | False |
| CygOB2-A | MIST | 3 | 2.6 | 1.87e+03 | 1.75e+03 | 2e+03 | 3.44e-06 | False |
| CygOB2-B | MIST | 3 | 2.6 | 2.05e+03 | 1.91e+03 | 2.2e+03 | 9.95e-15 | False |
| CygOB2-C | MIST | 3 | 2.6 | 2.19e+03 | 2.03e+03 | 2.37e+03 | 1.83e-07 | False |
| CygOB2-A | MIST | 3.1 | 2 | 1.15e+03 | 1.08e+03 | 1.22e+03 | 1.33e-19 | False |
| CygOB2-B | MIST | 3.1 | 2 | 1.07e+03 | 999 | 1.13e+03 | 7.67e-28 | False |
| CygOB2-C | MIST | 3.1 | 2 | 1.1e+03 | 1.03e+03 | 1.18e+03 | 1.31e-09 | False |
| CygOB2-A | MIST | 3.1 | 2.3 | 1.63e+03 | 1.53e+03 | 1.72e+03 | 9.77e-15 | False |
| CygOB2-B | MIST | 3.1 | 2.3 | 1.51e+03 | 1.42e+03 | 1.61e+03 | 2.1e-21 | False |
| CygOB2-C | MIST | 3.1 | 2.3 | 1.57e+03 | 1.47e+03 | 1.69e+03 | 2.26e-07 | False |
| CygOB2-A | MIST | 3.1 | 2.6 | 2.25e+03 | 2.1e+03 | 2.39e+03 | 6.43e-12 | False |
| CygOB2-B | MIST | 3.1 | 2.6 | 2.09e+03 | 1.95e+03 | 2.24e+03 | 7.1e-17 | False |
| CygOB2-C | MIST | 3.1 | 2.6 | 2.19e+03 | 2.02e+03 | 2.36e+03 | 6.26e-07 | False |
| CygOB2-A | MIST | 3.5 | 2 | 1.14e+03 | 1.07e+03 | 1.21e+03 | 2.56e-30 | False |
| CygOB2-B | MIST | 3.5 | 2 | 1.09e+03 | 1.02e+03 | 1.17e+03 | 2.5e-41 | False |
| CygOB2-C | MIST | 3.5 | 2 | 1.13e+03 | 1.05e+03 | 1.21e+03 | 5.61e-14 | False |
| CygOB2-A | MIST | 3.5 | 2.3 | 1.6e+03 | 1.51e+03 | 1.7e+03 | 9.59e-24 | False |
| CygOB2-B | MIST | 3.5 | 2.3 | 1.56e+03 | 1.45e+03 | 1.66e+03 | 2.43e-33 | False |
| CygOB2-C | MIST | 3.5 | 2.3 | 1.62e+03 | 1.5e+03 | 1.75e+03 | 5.43e-12 | False |
| CygOB2-A | MIST | 3.5 | 2.6 | 2.21e+03 | 2.07e+03 | 2.35e+03 | 9.63e-20 | False |
| CygOB2-B | MIST | 3.5 | 2.6 | 2.16e+03 | 2.01e+03 | 2.32e+03 | 8.44e-28 | False |
| CygOB2-C | MIST | 3.5 | 2.6 | 2.26e+03 | 2.08e+03 | 2.46e+03 | 5.05e-12 | False |
