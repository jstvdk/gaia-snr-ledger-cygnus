# WP1 extinction references

Frozen: 2026-07-22T10:18:37.348940+00:00

## Baseline: Vergely, Lallement & Cox (2022)

- Product: Vergely, Lallement & Cox (2022), A&A 664, A174 (bibcode `2022A&A...664A.174V`; DOI `10.1051/0004-6361/202243319`).
- Frozen archive: VizieR `J/A+A/664/A174`, catalogue DOI `10.26093/cds/vizier.36640174`, last catalogue modification reported as 2024-02-15. The five official cube choices are frozen in `data/raw/extinction/vergely2022_extinction_cube_list.ecsv`.
- Frozen WP3 cube: `data/raw/extinction/explore_cube_density_values_025pc_v2.fits` (117034560 bytes; SHA-256 `e8212ac5730a5b318bdb4dd981b9e66c4e2af367a0c3a3a0a91938cfcb0d1ee3`). It spans 6 x 6 x 0.8 kpc with 10 pc sampling and approximately 25 pc correlation length, so Cyg OB2 at about 1.35-1.6 kpc is inside the volume. The 5 kpc/20 pc cube is a fallback sensitivity check only.
- Units: the cube is differential monochromatic extinction density at 550 nm in nanomagnitude per parsec. Integrate along the Sun-to-star ray to the sampled distance. Compare the resulting A0 with fitted A_V while retaining the paper's warning that A0 and A_V are close, not identical.

## Check map: Dharmawardena et al. (2022) Cygnus X

- Product: Dharmawardena et al. (2022), A&A 658, A166 (bibcode `2022A&A...658A.166D`; DOI `10.1051/0004-6361/202141298`), a non-negative Gaussian-process 3D dust reconstruction made specifically for Cygnus X.
- Frozen archive: all five Cygnus products from VizieR `J/A+A/658/A166/Cygnus` (median density, median cumulative extinction, and l/b/d grid boundaries). Their hashes are in `provenance/wp1_extinction_refs_execution.json`.
- Use the provided boundary arrays rather than assuming pixel centres. Compare its cumulative extinction to the integrated Vergely A0 and the WP3 fitted A_V; retain the different map methodology as a systematic branch.

## Optional third diagnostic: Bayestar19

- Product: Green et al. (2019), ApJ 887, 93 (bibcode `2019ApJ...887...93G`; DOI `10.3847/1538-4357/ab5362`; data DOI `10.7910/DVN/2EJ9TX`).
- Native-to-E(B-V) diagnostic conversion: `E(B-V) = 0.98 x 0.901 x E_B19 = 0.882980 x E_B19`.
- Access check at freeze time: failed without blocking WP3 because two local maps are frozen: HTTPError: 500 Server Error: INTERNAL SERVER ERROR for url: http://argonaut.skymaps.info/api/v2/bayestar2019/query.
- The 36-point attempted coverage grid is retained in `data/raw/extinction/wp1_bayestar19_coverage_grid.parquet`; successful rows with reliable-distance flags: 0.

## Binding use in WP3 step 3b

1. Query both maps at the same sampled member distance; never query only to infinity.
2. Compare the Vergely and Dharmawardena maps with fitted extinction spatially and as residual distributions per subgroup; do not replace per-star Gaia+2MASS/spectral-type fits with either map.
3. Treat angular-resolution mismatch explicitly: the maps cannot reproduce arcminute-scale cloud structure in the Cyg OB2 core.
4. These maps are external pipelines but not photometrically independent of Gaia/2MASS. Agreement is a consistency check, not independent confirmation.
5. Carry native quality flags and map choice as a labeled systematic branch; do not average the maps.
