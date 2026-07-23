#!/usr/bin/env python3
"""WP3/WP4 shared input: fetch PARSEC (CMD 3.9) isochrones.

Solar metallicity (Z=0.0152), ages 1-10 Myr on a log grid (dlogage=0.05),
INTRINSIC (Av=0) synthetic photometry. Gaia EDR3 (= DR3 passbands) and 2MASS
are queried separately and merged on (logAge, Mini) because CMD has no single
Gaia+2MASS system. Reddening is applied downstream per R_V branch, so we only
need Av=0 absolute magnitudes here.

Output: data/processed/wp3_isochrones_parsec.parquet
Provenance appended by the caller (wp3_build_isochrones.py driver).
"""
import io, re, ssl, sys, time, urllib.request, urllib.parse
import numpy as np, pandas as pd

CMD_URL = "https://stev.oapd.inaf.it/cgi-bin/cmd_3.9"
BASE = "https://stev.oapd.inaf.it"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

PHOTSYS = {
    "gaia":  "YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat",
    "2mass": "YBC_tab_mag_odfnew/tab_mag_2mass_spitzer.dat",
}

def form(photsys_file):
    return {
        "cmd_version": "3.9",
        "track_parsec": "parsec_CAF09_v1.2S",
        "track_omegai": "0.00",
        "track_colibri": "parsec_CAF09_v1.2S_S_LMC_08_web",
        "track_postagb": "no",
        "n_inTPC": "10", "eta_reimers": "0.2",
        "kind_interp": "1", "kind_postagb": "-1",
        "photsys_file": photsys_file,
        "photsys_version": "YBCnewVega",
        "dust_sourceM": "dpmod60alox40", "dust_sourceC": "AMCSIC15",
        "kind_mag": "2", "kind_dust": "0",
        "extinction_av": "0.0",
        "extinction_coeff": "constant", "extinction_curve": "cardelli",
        "kind_LPV": "4",
        "imf_file": "tab_imf/imf_chabrier_lognormal_salpeter.dat",
        "isoc_isagelog": "1", "isoc_lagelow": "6.0", "isoc_lageupp": "7.005", "isoc_dlage": "0.05",
        "isoc_ismetlog": "0", "isoc_zlow": "0.0152", "isoc_zupp": "0.0152", "isoc_dz": "0.0",
        "isoc_metlow": "-2", "isoc_metupp": "0.3", "isoc_dmet": "0.0",
        "isoc_agelow": "1.0e9", "isoc_ageupp": "1.0e10", "isoc_dage": "0.0",
        "output_kind": "0", "output_evstage": "1",
        "lf_maginf": "-15", "lf_magsup": "20", "lf_deltamag": "0.5",
        "sim_mtot": "1.0e4", "output_gzip": "0",
        "submit_form": "Submit",
    }

def multipart(fields):
    boundary = "----wp3boundary1234567890"
    lines = []
    for k, v in fields.items():
        lines.append(f"--{boundary}")
        lines.append(f'Content-Disposition: form-data; name="{k}"')
        lines.append("")
        lines.append(str(v))
    lines.append(f"--{boundary}--")
    lines.append("")
    body = "\r\n".join(lines).encode()
    return body, boundary

def fetch(photsys_file):
    body, boundary = multipart(form(photsys_file))
    req = urllib.request.Request(CMD_URL, data=body, headers={
        "User-Agent": "Mozilla/5.0 wp3-pipeline",
        "Content-Type": f"multipart/form-data; boundary={boundary}"})
    html = urllib.request.urlopen(req, context=CTX, timeout=120).read().decode("utf-8", "ignore")
    m = re.search(r'href=["\']?(\.\./tmp/output\d+\.dat)', html)
    if not m:
        sys.stderr.write(html[:3000]); raise RuntimeError("no output link in CMD response")
    url = urllib.parse.urljoin("https://stev.oapd.inaf.it/cgi-bin/", m.group(1))
    txt = urllib.request.urlopen(url, context=CTX, timeout=120).read().decode("utf-8", "ignore")
    return txt, url

def parse_cmd(txt):
    header = None
    for line in txt.splitlines():
        if line.startswith("#") and ("Zini" in line or "logAge" in line):
            header = line.lstrip("#").split()
    df = pd.read_csv(io.StringIO(txt), sep=r"\s+", comment="#", header=None, engine="python")
    df.columns = header[:df.shape[1]]
    return df

def main():
    frames = {}
    prov = {}
    for tag, ps in PHOTSYS.items():
        for attempt in range(3):
            try:
                txt, url = fetch(ps); break
            except Exception as e:
                sys.stderr.write(f"[{tag}] attempt {attempt+1} failed: {e}\n"); time.sleep(5)
        else:
            raise RuntimeError(f"could not fetch {tag}")
        df = parse_cmd(txt)
        frames[tag] = df
        prov[tag] = {"photsys_file": ps, "output_url": url, "n_rows": int(len(df)),
                     "columns": list(df.columns)}
        print(f"[{tag}] {len(df)} rows, cols: {list(df.columns)[:12]}...")

    g = frames["gaia"]; t = frames["2mass"]
    key = ["Zini", "logAge", "Mini"]
    key = [k for k in key if k in g.columns and k in t.columns]
    # Gaia magnitude columns
    gcols = {c: c for c in g.columns if c in
             ("G_fSBmag", "G_BPmag", "G_RPmag", "Gmag", "G_BP_bmag", "G_BP_ftmag")}
    # standardized names
    gaia_map = {}
    for c in g.columns:
        cl = c.lower()
        if cl == "gmag": gaia_map[c] = "G0"
        elif "bp" in cl and "mag" in cl and "brmag" not in cl and gaia_map.get(c) is None and cl in ("g_bpmag","gbpmag","g_bp_mag"): gaia_map[c] = "BP0"
        elif "rp" in cl and "mag" in cl and cl in ("g_rpmag","grpmag","g_rp_mag"): gaia_map[c] = "RP0"
    tmass_map = {}
    for c in t.columns:
        cl = c.lower()
        if cl == "jmag": tmass_map[c] = "J0"
        elif cl == "hmag": tmass_map[c] = "H0"
        elif cl in ("ksmag", "kmag"): tmass_map[c] = "Ks0"
    print("gaia_map", gaia_map); print("tmass_map", tmass_map)

    extra = [c for c in ["label", "Mass", "logL", "logTe", "logg", "MH"] if c in g.columns]
    keep_g = key + extra + list(gaia_map)
    gg = g[keep_g].rename(columns=gaia_map)
    tt = t[key + list(tmass_map)].rename(columns=tmass_map)
    merged = gg.merge(tt, on=key, how="inner")
    merged["age_Myr"] = 10 ** merged["logAge"] / 1e6
    merged["family"] = "PARSEC"
    merged["track"] = "parsec_CAF09_v1.2S"
    out = "data/processed/wp3_isochrones_parsec.parquet"
    merged.to_parquet(out, index=False)
    print(f"WROTE {out}  {merged.shape}  ages {merged.age_Myr.min():.2f}-{merged.age_Myr.max():.2f} Myr "
          f"({merged.logAge.nunique()} ages)")
    import json
    prov["merged"] = {"rows": int(len(merged)), "columns": list(merged.columns),
                      "n_ages": int(merged.logAge.nunique()),
                      "mass_range": [float(merged["Mini"].min()), float(merged["Mini"].max())]}
    with open("provenance/wp3_parsec_fetch_execution.json", "w") as f:
        json.dump(prov, f, indent=2)

if __name__ == "__main__":
    main()
