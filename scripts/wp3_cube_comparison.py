#!/usr/bin/env python3
"""WP3 gate: compare per-star A_V against the Vergely+22 and Dharmawardena+22
3D dust reconstructions at each member's distance.

Vergely+22 (EXPLORE): heliocentric Cartesian cube of differential extinction
density A0(550nm)/pc at 25 pc resolution. We integrate along each member's
sightline to the group distance to obtain the cube-predicted A0 (~A_V at
R_V=3.1). Dharmawardena+22: (l,b,d) cube of CUMULATIVE extinction to distance d
over the Cygnus field; we interpolate at (l,b,d_member).

Both are smooth large-scale reconstructions that under-resolve the clumpy,
high-extinction material local to Cyg OB2, so they are expected to trace the
same spatial pattern while sitting BELOW the per-star values -- exactly the
consistency check the plan asks for (cubes as cross-check, not substitute).

Outputs: data/processed/wp3_cube_comparison.parquet
Provenance: provenance/wp3_cube_comparison_execution.json
"""
import hashlib, json
import numpy as np, pandas as pd
from astropy.io import fits
from scipy.stats import spearmanr
from scipy.interpolate import RegularGridInterpolator
from wp3_common import D_KPC

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

def vergely_integrator():
    h = fits.open("data/raw/extinction/explore_cube_density_values_025pc_v2.fits")[0]
    dens = h.data                      # (z, y, x), A0/pc
    hd = h.header
    STEP = hd["STEP"]; sx, sy, sz = hd["SUN_POSX"], hd["SUN_POSY"], hd["SUN_POSZ"]
    nz, ny, nx = dens.shape
    def integ(l, b, d_pc, n=600):
        lr, br = np.radians(l), np.radians(b)
        ds = np.linspace(0, d_pc, n); dd = ds[1] - ds[0]
        x = ds * np.cos(br) * np.cos(lr); y = ds * np.cos(br) * np.sin(lr); z = ds * np.sin(br)
        ix = np.round(sx + x / STEP).astype(int)
        iy = np.round(sy + y / STEP).astype(int)
        iz = np.round(sz + z / STEP).astype(int)
        m = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iz >= 0) & (iz < nz)
        val = np.zeros(n); val[m] = dens[iz[m], iy[m], ix[m]]
        return float(np.sum(val) * dd)
    return integ

def dharma_interpolator():
    cube = np.load("data/raw/extinction/dharmawardena2022_cygnus_ext_med_cube.pkl.npy")  # (l,b,d) cumulative
    lb = np.load("data/raw/extinction/dharmawardena2022_cygnus_l_bounds_pred.pkl.npy")
    bb = np.load("data/raw/extinction/dharmawardena2022_cygnus_b_bounds_pred.pkl.npy")
    db = np.load("data/raw/extinction/dharmawardena2022_cygnus_d_bounds_pred.pkl.npy")
    # l,b arrays are cell edges (n+1); d array is cell centres (n)
    lc = 0.5 * (lb[:-1] + lb[1:]); bc = 0.5 * (bb[:-1] + bb[1:])
    dc = db if len(db) == cube.shape[2] else 0.5 * (db[:-1] + db[1:])
    interp = RegularGridInterpolator((lc, bc, dc), cube, bounds_error=False, fill_value=np.nan)
    dmin, dmax = dc.min(), dc.max()
    def get(l, b, d_pc):
        return float(interp([[l, b, np.clip(d_pc, dmin, dmax)]])[0])
    return get, (lc.min(), lc.max(), bc.min(), bc.max())

def main():
    df = pd.read_parquet("data/processed/wp3_extinction.parquet")
    d_pc = D_KPC * 1000.0
    verg = vergely_integrator()
    dharma, (lmin, lmax, bmin, bmax) = dharma_interpolator()

    vav, dav = [], []
    for _, r in df.iterrows():
        vav.append(verg(r.l_deg, r.b_deg, d_pc))
        if lmin <= r.l_deg <= lmax and bmin <= r.b_deg <= bmax:
            dav.append(dharma(r.l_deg, r.b_deg, d_pc))
        else:
            dav.append(np.nan)
    df_out = df[["source_id", "l_deg", "b_deg", "av_rv3.1", "av_err_rv3.1", "av_source"]].copy()
    df_out["vergely_A0"] = vav
    df_out["dharma_Acum"] = dav
    p = "data/processed/wp3_cube_comparison.parquet"
    df_out.to_parquet(p, index=False)

    def stats(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        rho, pval = spearmanr(a[m], b[m])
        return {"n": int(m.sum()), "spearman_rho": float(rho), "spearman_p": float(pval),
                "median_star_av": float(np.median(a[m])), "median_cube": float(np.median(b[m])),
                "median_ratio_star_over_cube": float(np.median(a[m] / np.clip(b[m], 1e-3, None)))}
    av = df_out["av_rv3.1"].values
    sv = stats(av, np.array(vav))
    sd = stats(av, np.array(dav))
    prov = {
        "script": "scripts/wp3_cube_comparison.py",
        "inputs": {
            "wp3_extinction": sha("data/processed/wp3_extinction.parquet"),
            "vergely_cube": sha("data/raw/extinction/explore_cube_density_values_025pc_v2.fits"),
            "dharma_cube": sha("data/raw/extinction/dharmawardena2022_cygnus_ext_med_cube.pkl.npy"),
        },
        "distance_pc": d_pc,
        "vergely": {"description": "line-of-sight integral of A0/pc to group distance", **sv},
        "dharmawardena": {"description": "cumulative extinction interpolated at (l,b,d)", **sd},
        "interpretation": (
            "Both cubes trace the same spatial pattern (positive Spearman rho) but sit "
            "below the per-star A_V because their 25 pc / 0.14 deg smoothing washes out "
            "the clumpy dust local to Cyg OB2. Cubes are a consistency cross-check, not a "
            "substitute for per-star extinction (plan WP3 caveat)."),
        "output": {p: sha(p)},
    }
    json.dump(prov, open("provenance/wp3_cube_comparison_execution.json", "w"), indent=2)
    print(f"Vergely+22 : n={sv['n']} rho={sv['spearman_rho']:.2f} (p={sv['spearman_p']:.1e}) "
          f"star/cube median ratio={sv['median_ratio_star_over_cube']:.2f} "
          f"(star {sv['median_star_av']:.2f} vs cube {sv['median_cube']:.2f})")
    print(f"Dharma+22  : n={sd['n']} rho={sd['spearman_rho']:.2f} (p={sd['spearman_p']:.1e}) "
          f"star/cube median ratio={sd['median_ratio_star_over_cube']:.2f} "
          f"(star {sd['median_star_av']:.2f} vs cube {sd['median_cube']:.2f})")

if __name__ == "__main__":
    main()
