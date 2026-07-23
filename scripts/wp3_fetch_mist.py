#!/usr/bin/env python3
"""WP3/WP4 shared input: fetch MIST v1.2 isochrones (mist.science interpolator).

Solar [Fe/H]=0, vvcrit=0.4, ages 1-10 Myr on a log grid (dlogage=0.05),
INTRINSIC (Av=0) synthetic photometry. The UBVRIplus bundle contains both
Gaia (EDR3) and 2MASS bands in a single file. Reddening is applied downstream
per R_V branch.

Output: data/processed/wp3_isochrones_mist.parquet
"""
import io, os, re, ssl, sys, json, time, zipfile, urllib.request
import numpy as np, pandas as pd

BASE = "https://mist.science"
FORM_URL = BASE + "/iso_form.php"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

FIELDS = {
    "version": "MIST1",
    "v_div_vcrit": "vvcrit0.4",
    "age_scale": "log10",
    "age_type": "range",
    "age_range_low": "6.0",
    "age_range_high": "7.005",
    "age_range_delta": "0.05",
    "age_value": "", "age_list": "",
    "FeH_value": "0",
    "alpha_value": "p0",
    "output_option": "photometry",
    "output": "UBVRIplus",
    "Av_value": "0",
}

def multipart(fields):
    boundary = "----wp3mistboundary987"
    lines = []
    for k, v in fields.items():
        lines += [f"--{boundary}", f'Content-Disposition: form-data; name="{k}"', "", str(v)]
    lines += [f"--{boundary}--", ""]
    return "\r\n".join(lines).encode(), boundary

def post():
    body, boundary = multipart(FIELDS)
    req = urllib.request.Request(FORM_URL, data=body, headers={
        "User-Agent": "Mozilla/5.0 wp3-pipeline",
        "Content-Type": f"multipart/form-data; boundary={boundary}"})
    r = urllib.request.urlopen(req, context=CTX, timeout=180)
    ctype = r.headers.get("Content-Type", "")
    raw = r.read()
    return raw, ctype, r.geturl()

def _from_zip(raw, src):
    zf = zipfile.ZipFile(io.BytesIO(raw))
    names = zf.namelist()
    # prefer the synthetic-photometry member (contains Gaia/2MASS columns)
    phot = [n for n in names if "UBVRIplus" in n or n.endswith(".iso.cmd") or ".cmd" in n]
    name = phot[0] if phot else names[0]
    return zf.read(name).decode("utf-8", "ignore"), f"{src} (zip:{name})"

def get_iso_text():
    raw, ctype, url = post()
    # Case 1: direct zip
    if raw[:2] == b"PK":
        return _from_zip(raw, url)
    txt = raw.decode("utf-8", "ignore")
    # Case 2: HTML with a link to the result file/zip
    m = re.search(r'href="([^"]+\.(?:zip|iso\.cmd|cmd|iso))"', txt)
    if m:
        link = m.group(1)
        if not link.startswith("http"):
            link = BASE + "/" + link.lstrip("/")
        data = urllib.request.urlopen(link, context=CTX, timeout=180).read()
        if data[:2] == b"PK":
            return _from_zip(data, link)
        return data.decode("utf-8", "ignore"), link
    # Case 3: response IS the cmd text
    if "log10_isochrone_age_yr" in txt or "isochrone_age" in txt or "EEP" in txt:
        return txt, url
    sys.stderr.write(f"ctype={ctype} len={len(raw)} head={txt[:800]}\n")
    raise RuntimeError("MIST: could not locate isochrone data")

def parse_mist(txt):
    lines = txt.splitlines()
    header = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("#") and ("EEP" in ln and "isochrone_age" in ln):
            header = ln.lstrip("#").split()
            data_start = i + 1
    if header is None:
        # header is the last comment line before data
        for i, ln in enumerate(lines):
            if ln.lstrip().startswith("# EEP") or (ln.startswith("#") and "initial_mass" in ln):
                header = ln.lstrip("#").split(); data_start = i + 1
    df = pd.read_csv(io.StringIO(txt), sep=r"\s+", comment="#", header=None, engine="python")
    df.columns = header[:df.shape[1]]
    return df

def main():
    for attempt in range(3):
        try:
            txt, src = get_iso_text(); break
        except Exception as e:
            sys.stderr.write(f"attempt {attempt+1}: {e}\n"); time.sleep(6)
    else:
        raise SystemExit("MIST fetch failed")
    open("scratch_mist_raw.txt", "w").write(txt[:5000])
    df = parse_mist(txt)
    print("MIST raw cols:", list(df.columns))
    # standardize
    def find(cols, *cands):
        for cand in cands:
            for c in cols:
                if c.lower() == cand.lower():
                    return c
        for cand in cands:
            for c in cols:
                if cand.lower() in c.lower():
                    return c
        return None
    cols = list(df.columns)
    ren = {
        find(cols, "log10_isochrone_age_yr"): "logAge",
        find(cols, "initial_mass"): "Mini",
        find(cols, "star_mass"): "Mass",
        find(cols, "log_Teff"): "logTe",
        find(cols, "log_g"): "logg",
        find(cols, "log_L"): "logL",
        find(cols, "phase"): "phase",
        find(cols, "Gaia_G_EDR3", "Gaia_G_DR2Rev", "Gaia_G"): "G0",
        find(cols, "Gaia_BP_EDR3", "Gaia_BP_DR2Rev", "Gaia_BP"): "BP0",
        find(cols, "Gaia_RP_EDR3", "Gaia_RP_DR2Rev", "Gaia_RP"): "RP0",
        find(cols, "2MASS_J"): "J0",
        find(cols, "2MASS_H"): "H0",
        find(cols, "2MASS_Ks"): "Ks0",
    }
    ren = {k: v for k, v in ren.items() if k is not None}
    df = df.rename(columns=ren)
    need = ["logAge", "Mini", "G0", "BP0", "RP0", "J0", "H0", "Ks0"]
    missing = [n for n in need if n not in df.columns]
    if missing:
        raise SystemExit(f"MIST missing columns {missing}; have {list(df.columns)}")
    keep = [c for c in ["logAge", "Mini", "Mass", "logTe", "logg", "logL", "phase",
                        "G0", "BP0", "RP0", "J0", "H0", "Ks0"] if c in df.columns]
    out = df[keep].copy()
    out["age_Myr"] = 10 ** out["logAge"] / 1e6
    out["family"] = "MIST"
    out["track"] = "MIST_v1.2_vvcrit0.4_feh0"
    # restrict to 1-10 Myr
    out = out[(out.age_Myr >= 0.99) & (out.age_Myr <= 10.5)]
    p = "data/processed/wp3_isochrones_mist.parquet"
    out.to_parquet(p, index=False)
    print(f"WROTE {p}  {out.shape}  ages {out.age_Myr.min():.2f}-{out.age_Myr.max():.2f} Myr "
          f"({out.logAge.nunique()} ages)")
    prov = {"source": src, "fields": FIELDS, "rows": int(len(out)),
            "n_ages": int(out.logAge.nunique()), "columns_renamed": ren,
            "mass_range": [float(out.Mini.min()), float(out.Mini.max())]}
    json.dump(prov, open("provenance/wp3_mist_fetch_execution.json", "w"), indent=2)

if __name__ == "__main__":
    main()
