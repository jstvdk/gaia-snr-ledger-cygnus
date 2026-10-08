# Transfer `data/` from the external SSD (instructions for an agent)

Written 2026-10-08. Follow these steps in order, from the repository root, after `git pull`. Do not improvise around a failure: **never regenerate data, and never delete anything except macOS `._*` metadata files.** If a step fails, stop and tell the owner.

## What is on the SSD

- **Drive:** "Crucial2TB" (exFAT).
- **Data:** folder `cygnus_data/` at the top level of the SSD. It is a copy of the project's `data/` tree (`raw/` and `processed/`, about 13 GB, 3,270 files), state of 2026-10-08, including the repair_v9 products.
- **Checksums:** `provenance/data_manifest_2026-10-08.sha256` in the repository lists every file with its SHA-256. A copy, `cygnus_data_manifest_2026-10-08.sha256`, sits at the top level of the SSD.
- **Verified:** the SSD copy was checked against this manifest after writing. All 3,270 files match.

## Steps

1. **Find the SSD.**
   - macOS: `/Volumes/Crucial2TB`.
   - Linux: `/run/media/$USER/Crucial2TB`, or check `lsblk -o NAME,LABEL,MOUNTPOINT`.
   - Set `SSD=<that path>` and confirm that `"$SSD/cygnus_data/processed"` exists.
2. **Check the repository.**
   - `git log --oneline -1` must show a commit at least as new as the one that added this file.
   - `provenance/data_manifest_2026-10-08.sha256` must exist. If it doesn't, use `"$SSD/cygnus_data_manifest_2026-10-08.sha256"` instead.
3. **Check free space:** you need at least 14 GB on the repository's disk (`df -h .`).
4. **Copy, without deleting anything.** The `data/` folder may already exist with older files; rsync only adds and updates.
   ```sh
   mkdir -p data
   rsync -rt --exclude '._*' --exclude '.DS_Store' "$SSD/cygnus_data/" data/
   ```
   - The trailing slashes matter. The result must be `data/raw` and `data/processed`, not `data/cygnus_data/...`.
   - Do **not** add `--delete`.
   - macOS's built-in rsync (openrsync) accepts these flags. Expect about 5–15 minutes over USB.
5. **Remove macOS metadata files only:**
   ```sh
   find data -name '._*' -size -5k -delete
   ```
6. **Verify every file.** No output followed by "OK" means every file matches.
   ```sh
   # macOS
   (cd data && shasum -a 256 -c ../provenance/data_manifest_2026-10-08.sha256 --quiet) && echo OK
   # Linux
   (cd data && sha256sum -c ../provenance/data_manifest_2026-10-08.sha256 --quiet) && echo OK
   ```
   - If files fail, re-copy them by content:
     ```sh
     rsync -rt --checksum --exclude '._*' "$SSD/cygnus_data/" data/
     ```
     Then verify again.
   - If anything still fails, stop and report the failing paths to the owner.
   - Files that exist in `data/` but are missing from the manifest are not an error.
7. **Check the project sees the data.** The environment comes from `INSTALL.md`.
   ```sh
   PYTHONPATH=scripts python scripts/wp10_inputs.py      # must end with: audit: PASS
   git diff --stat provenance/wp10_input_manifest.json   # timestamp/platform lines only
   git checkout provenance/wp10_input_manifest.json       # restore the committed record
   ```
8. **Report** to the owner: files copied, verification OK, audit PASS.
9. **Continue the project:** read `reports/repair_v9_completion_report.md` §0 ("Handoff: start here").

## Notes

- `data/` is gitignored and has no archive other than these copies. Keep the SSD copy until the owner says otherwise.
- **Run the analysis chain on only one machine at a time.** If you produce new data, copy it back to the SSD the same way, with source and destination reversed, and write a new manifest:
  ```sh
  (cd data && find . -type f ! -name '._*' ! -name '.DS_Store' -print0 | sort -z | xargs -0 shasum -a 256) > provenance/data_manifest_<date>.sha256
  ```
  On Linux, use `sha256sum` in place of `shasum -a 256`.
