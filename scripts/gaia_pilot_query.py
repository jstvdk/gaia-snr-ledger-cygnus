import pyvo

tap = pyvo.dal.TAPService(
    "https://gea.esac.esa.int/tap-server/tap"
)

query = """
SELECT TOP 100
    source_id,
    ra, dec,
    parallax, parallax_error,
    pmra, pmra_error,
    pmdec, pmdec_error,
    phot_g_mean_mag,
    phot_bp_mean_mag,
    phot_rp_mean_mag,
    ruwe,
    astrometric_params_solved,
    radial_velocity,
    radial_velocity_error
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', 308.3, 41.2, 5.0)
)
AND parallax BETWEEN 0.35 AND 1.10
AND phot_g_mean_mag < 19
"""

result = tap.search(query)
table = result.to_table()

print(len(table))
print(table.colnames)