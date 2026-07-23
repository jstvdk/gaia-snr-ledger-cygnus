# wp1_gaia_wide

Content: exact Gaia DR3 wide box for WP6 runaway/traceback work, assembled from four non-overlapping Galactic-longitude tiles.

Selection: l=72-88 deg, b=-5-8 deg, parallax=0.35-1.10 mas, G<19. Raw parallax is retained; zero-point correction remains a downstream operation.

Frozen: 2026-07-22T12:35:59.907877+00:00; rows: 3133326; unique source_ids: 3133326; SHA-256: `afb54bb3e88c52c0c76546fe151a9abed107ed230516b36aa7ca046e6a4be7c6`.

Format: Parquet is the canonical analysis artifact. A FITS duplicate is omitted because this >3 million-row table is a WP6 support catalogue and Parquet preserves the nullable columns with substantially less storage.
