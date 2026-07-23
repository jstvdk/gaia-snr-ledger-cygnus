#!/usr/bin/env python
"""Download frozen WP1 Gaia DR3 queries through TAP with resume/retry support."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyvo
import requests
from astropy.table import MaskedColumn
from pyvo.dal import AsyncTAPJob
from pyvo.dal.exceptions import DALServiceError
from urllib3.exceptions import HTTPError as Urllib3HTTPError


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TAP_URL = "https://gea.esac.esa.int/tap-server/tap"
TERMINAL_PHASES = {"COMPLETED", "ERROR", "ABORTED"}
RETRYABLE = (
    DALServiceError,
    requests.RequestException,
    Urllib3HTTPError,
    TimeoutError,
    OSError,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def fits_compatible_table(table):
    """Return a copy whose variable-length object strings are FITS-safe."""
    safe = table.copy()
    for name in safe.colnames:
        column = safe[name]
        if column.dtype.kind != "O":
            continue
        mask = np.ma.getmaskarray(column)
        values = np.asarray(
            ["" if is_masked else str(value) for value, is_masked in zip(column, mask)],
            dtype=str,
        )
        safe.replace_column(
            name,
            MaskedColumn(
                values,
                mask=mask,
                name=name,
                unit=getattr(column, "unit", None),
                description=getattr(column, "description", None),
                meta=getattr(column, "meta", None),
            ),
        )
    return safe


def retry_call(label: str, function, deadline: float, *, max_delay: int = 180):
    delay = 5
    attempt = 0
    while True:
        try:
            return function()
        except RETRYABLE as exc:
            attempt += 1
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise
            sleep_for = min(delay, remaining)
            print(
                f"{label}: transient connection error "
                f"({type(exc).__name__}); retry {attempt} in {sleep_for:.0f}s",
                flush=True,
            )
            time.sleep(sleep_for)
            delay = min(delay * 2, max_delay)


def run_query(name: str, timeout: int, new_job: bool, poll_interval: int) -> None:
    query_path = PROJECT_ROOT / "queries" / f"gaia_{name}.adql"
    output_dir = PROJECT_ROOT / "data" / "raw" / "gaia"
    provenance_dir = PROJECT_ROOT / "provenance"
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir.mkdir(parents=True, exist_ok=True)

    query_hash = sha256(query_path)
    state_path = provenance_dir / f"wp1_gaia_{name}_job.json"
    execution_path = provenance_dir / f"wp1_gaia_{name}_execution.json"
    fits_path = output_dir / f"wp1_gaia_{name}.fits"
    parquet_path = output_dir / f"wp1_gaia_{name}.parquet"
    deadline = time.monotonic() + timeout
    started = datetime.now(timezone.utc).isoformat()
    state = None

    if state_path.exists() and not new_job:
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("query_sha256") != query_hash:
            raise RuntimeError(
                f"Existing job state uses a different query. "
                f"Use --new-job only after reviewing {state_path}."
            )
        job = retry_call(
            "reopen TAP job",
            lambda: AsyncTAPJob(state["tap_job_url"], delete=False),
            deadline,
        )
        print(f"{name}: resuming TAP job {job.url}", flush=True)
    else:
        tap = pyvo.dal.TAPService(TAP_URL)
        job = tap.submit_job(query_path.read_text(encoding="utf-8"))
        # Keep the server-side job alive if this local process exits.
        job._delete_on_exit = False
        state = {
            "query_name": name,
            "tap_url": TAP_URL,
            "tap_job_url": job.url,
            "query_file": str(query_path.relative_to(PROJECT_ROOT)),
            "query_sha256": query_hash,
            "submitted_utc": started,
            "run_requested": False,
            "status": "SUBMITTED",
        }
        # Persist the URL before starting the job: this is the critical resume point.
        write_json(state_path, state)
        state["run_requested"] = True
        write_json(state_path, state)
        job.run()
        print(f"{name}: submitted TAP job {job.url}", flush=True)

    try:
        # A previously submitted but not-started job can still be started on resume.
        phase = retry_call("read TAP phase", lambda: job.phase, deadline)
        if phase == "PENDING":
            state["run_requested"] = True
            write_json(state_path, state)
            job.run()

        while True:
            phase = retry_call("poll TAP phase", lambda: job.phase, deadline)
            state["last_phase"] = phase
            state["last_poll_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(state_path, state)
            print(f"{name}: TAP phase={phase}", flush=True)
            if phase in TERMINAL_PHASES:
                break
            time.sleep(poll_interval)

        if phase != "COMPLETED":
            state["status"] = phase
            write_json(state_path, state)
            raise RuntimeError(f"TAP job finished in phase {phase}: {job.url}")

        results = retry_call("fetch TAP result", job.fetch_result, deadline)
        table = results.to_table()
        export_table = fits_compatible_table(table)
        export_table.write(fits_path, format="fits", overwrite=True)
        export_table.write(parquet_path, format="parquet", overwrite=True)

        state.update({
            "status": "COMPLETED",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "row_count": len(table),
            "column_count": len(table.colnames),
        })
        write_json(state_path, state)

        record = {
            "executed_utc": started,
            "query_name": name,
            "tap_url": TAP_URL,
            "tap_job_url": job.url,
            "tap_phase": phase,
            "query_file": str(query_path.relative_to(PROJECT_ROOT)),
            "query_sha256": query_hash,
            "row_count": len(table),
            "column_count": len(table.colnames),
            "columns": list(table.colnames),
            "outputs": {
                str(fits_path.relative_to(PROJECT_ROOT)): sha256(fits_path),
                str(parquet_path.relative_to(PROJECT_ROOT)): sha256(parquet_path),
            },
        }
        write_json(execution_path, record)
        print(json.dumps(record, indent=2), flush=True)
    except Exception as exc:
        state["status"] = "LOCAL_ERROR"
        state["error_type"] = type(exc).__name__
        state["error_message"] = str(exc)
        state["error_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(state_path, state)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="query name, e.g. narrow, wide, or narrow_tile01")
    parser.add_argument("--timeout", type=int, default=6 * 60 * 60)
    parser.add_argument("--poll-interval", type=int, default=30)
    parser.add_argument(
        "--new-job",
        action="store_true",
        help="submit a new job instead of resuming the saved job state",
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.name):
        parser.error("name may contain only letters, numbers, underscore, and hyphen")
    query_path = PROJECT_ROOT / "queries" / f"gaia_{args.name}.adql"
    if not query_path.is_file():
        parser.error(f"query file not found: {query_path}")
    run_query(args.name, args.timeout, args.new_job, args.poll_interval)


if __name__ == "__main__":
    main()
