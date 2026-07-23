"""
Provenance / consistency audit dump for gaia_snr_history_cygnus.

Run from the project root under the project environment:
    conda run -n cygob2-gaia --no-capture-output python audit.py

Writes AUDIT.txt so the state can be checked against the plan.

The project environment is mandatory.  A bare Python 3.14 installation has no
pandas and historically skipped the entire data-products section; this version
marks that condition as an audit failure instead of silently treating it as a
successful run.
"""

import hashlib
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path.cwd()
OUT = ROOT / "AUDIT.txt"
SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", ".venv", "venv", "node_modules"}
HASH_EXT = {".parquet", ".csv", ".fits", ".npy", ".npz", ".json", ".yaml", ".yml", ".md", ".py"}
MAX_HASH_BYTES = 500 * 1024 * 1024  # skip hashing anything above this

lines = []


def emit(s=""):
    lines.append(str(s))


def section(title):
    emit()
    emit("=" * 78)
    emit(title)
    emit("=" * 78)


def sh(cmd):
    try:
        return subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=60
        ).stdout.strip()
    except Exception as e:  # noqa: BLE001
        return f"<failed: {e}>"


def sha256(path):
    if path.stat().st_size > MAX_HASH_BYTES:
        return "<too large>"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def mtime(path):
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")


# ---------------------------------------------------------------- header
emit(f"AUDIT of {ROOT}")
emit(f"generated {datetime.now().isoformat(timespec='seconds')}")
emit(f"python {sys.version.split()[0]}")

# ---------------------------------------------------------------- git
section("GIT PROVENANCE")
if (ROOT / ".git").exists():
    emit("HEAD:")
    emit(sh("git log -1 --format='%H%n  author: %an%n  date:   %ai%n  msg:    %s'"))
    emit()
    emit("last 25 commits:")
    emit(sh("git log -25 --format='%h %ai %s'"))
    emit()
    emit("working tree status (uncommitted work = provenance gap):")
    st = sh("git status --porcelain")
    emit(st if st else "  clean")
    emit()
    emit("branch: " + sh("git rev-parse --abbrev-ref HEAD"))
else:
    emit("NO GIT REPOSITORY. Provenance cannot be reconstructed from history.")

# ---------------------------------------------------------------- inventory
section("FILE INVENTORY (mtime | size | sha256-16 | path)")
records = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
    for fn in sorted(filenames):
        p = Path(dirpath) / fn
        if p.name == OUT.name or p.is_symlink():
            continue
        try:
            stat = p.stat()
        except OSError:
            continue
        digest = sha256(p) if p.suffix.lower() in HASH_EXT else "-"
        rel = p.relative_to(ROOT)
        records.append((stat.st_mtime, mtime(p), stat.st_size, digest, str(rel)))

for _, mt, size, digest, rel in sorted(records, key=lambda r: r[4]):
    emit(f"{mt}  {size:>12,}  {digest:>16}  {rel}")
emit()
emit(f"total files: {len(records)}")

# ---------------------------------------------------------------- ordering
section("BUILD-ORDER CHECK (outputs older than inputs = STALE)")
emit("Files sorted newest-first. A WP-N product that predates a WP-(N-1)")
emit("product means the downstream stage was not re-run after the upstream")
emit("changed. Check the WP ordering here by eye.")
emit()
for _, mt, size, digest, rel in sorted(records, key=lambda r: -r[0])[:40]:
    emit(f"{mt}  {rel}")

# ---------------------------------------------------------------- tabular
section("DATA PRODUCTS: schema, row counts, key column distributions")
try:
    import pandas as pd
except ImportError:
    emit("*** AUDIT FAILURE: pandas not available; data products were NOT audited.")
    emit("Re-run with: conda run -n cygob2-gaia --no-capture-output python audit.py")
    pd = None

if pd is not None:
    tabular = [
        Path(r[4]) for r in records if Path(r[4]).suffix.lower() in {".parquet", ".csv"}
    ]
    for rel in tabular:
        p = ROOT / rel
        emit()
        emit("-" * 78)
        emit(f"FILE: {rel}")
        try:
            df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p)
        except Exception as e:  # noqa: BLE001
            emit(f"  <could not read: {e}>")
            continue
        emit(f"  rows: {len(df):,}   cols: {len(df.columns)}")
        emit(f"  columns: {list(df.columns)}")

        # membership / probability cuts
        for col in df.columns:
            lc = col.lower()
            if "probab" in lc or lc in {"p", "prob"}:
                try:
                    emit(f"  {col}: >0.5 -> {int((df[col] > 0.5).sum()):,}   "
                         f">0.7 -> {int((df[col] > 0.7).sum()):,}   "
                         f"median {df[col].median():.4f}")
                except Exception:  # noqa: BLE001
                    pass

        # label / flag columns
        for col in df.columns:
            lc = col.lower()
            if any(k in lc for k in ("label", "subgroup", "flag", "exempt", "class")):
                try:
                    vc = df[col].value_counts(dropna=False)
                    emit(f"  {col} value counts:")
                    for k, v in vc.head(15).items():
                        emit(f"      {str(k)[:60]:<62} {v:,}")
                    if len(vc) > 15:
                        emit(f"      ... ({len(vc)} distinct values)")
                except Exception:  # noqa: BLE001
                    pass

        # The subgroup schema has one canonical name and must never carry the
        # pre-closure distance-structure placeholder.
        subgroup_columns = [
            col for col in df.columns if col.lower() in {"subgroup", "subgroup_label"}
        ]
        frozen_failure = "failed" in p.name or "failure" in p.name
        if subgroup_columns and subgroup_columns != ["subgroup"] and not frozen_failure:
            emit(f"  *** SUBGROUP SCHEMA FAILURE: {subgroup_columns}")
        for col in ([] if frozen_failure else subgroup_columns):
            vals = df[col].astype(str)
            hits = vals.str.contains(
                "placeholder|distance_structure_unresolved|TODO|TBD|FIXME|dummy",
                case=False, na=False,
            )
            if hits.any():
                emit(f"  *** SUBGROUP PLACEHOLDER VALUES in '{col}': {int(hits.sum()):,} rows")

        # PLACEHOLDER DETECTION
        for col in df.select_dtypes(include=["object", "string"]).columns:
            try:
                vals = df[col].astype(str)
                hits = vals.str.contains(
                    "placeholder|TODO|TBD|FIXME|dummy", case=False, na=False
                )
                if hits.any():
                    emit(f"  *** PLACEHOLDER VALUES in '{col}': {int(hits.sum()):,} rows")
                    for v in vals[hits].unique()[:5]:
                        emit(f"        {v[:80]}")
            except Exception:  # noqa: BLE001
                pass

        # numeric sanity
        try:
            num = df.select_dtypes("number").replace(
                [float("inf"), float("-inf")], float("nan")
            )
            if len(num.columns):
                emit("  numeric summary (count / mean / std / min / max):")
                desc = num.describe().T[["count", "mean", "std", "min", "max"]]
                for name, row in desc.iterrows():
                    emit(f"      {str(name)[:28]:<30} {row['count']:>8.0f} "
                         f"{row['mean']:>12.4g} {row['std']:>12.4g} "
                         f"{row['min']:>12.4g} {row['max']:>12.4g}")
                nn = df.isna().sum()
                nn = nn[nn > 0]
                if len(nn):
                    emit("  columns with NaNs:")
                    for name, v in nn.items():
                        emit(f"      {str(name)[:40]:<42} {v:,}")
        except Exception:  # noqa: BLE001
            pass

# ---------------------------------------------------------------- docs
section("MARKDOWN DOCUMENTS: headings only")
for rel in sorted(Path(r[4]) for r in records if Path(r[4]).suffix.lower() == ".md"):
    p = ROOT / rel
    emit()
    emit(f"--- {rel}  (modified {mtime(p)}) ---")
    try:
        text = p.read_text(errors="replace")
    except Exception as e:  # noqa: BLE001
        emit(f"  <unreadable: {e}>")
        continue
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("#"):
            emit("  " + s[:110].rstrip())
        elif any(k in s.upper() for k in ("GATE", "PASS", "FAIL", "BLOCK", "TODO", "OPEN QUESTION")):
            emit("  > " + s[:110].rstrip())

# ---------------------------------------------------------------- figures
section("FIGURES")
figs = [
    r for r in records
    if Path(r[4]).suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".eps"}
]
if figs:
    for _, mt, size, _, rel in sorted(figs, key=lambda r: r[4]):
        emit(f"{mt}  {size:>10,}  {rel}")
else:
    emit("none found")

# ---------------------------------------------------------------- config
section("CONFIG / THRESHOLD FILES (full contents)")
for rel in sorted(Path(r[4]) for r in records):
    name = rel.name.lower()
    if (
        rel.suffix.lower() in {".yaml", ".yml", ".toml", ".ini", ".cfg"}
        or "threshold" in name
        or "config" in name
        or name in {"requirements.txt", "environment.yml", "pyproject.toml"}
    ):
        p = ROOT / rel
        emit()
        emit(f"--- {rel} ---")
        try:
            emit(p.read_text(errors="replace")[:6000])
        except Exception as e:  # noqa: BLE001
            emit(f"  <unreadable: {e}>")

OUT.write_text("\n".join(lines))
print(f"wrote {OUT}  ({len(lines)} lines)")
