# Setting up this project on a new machine

The repository holds code, tables, reports and provenance records. The data tree does not come with it. Three things are needed:

1. the code (`git clone`);
2. the Python environment (conda);
3. the data tree `data/` (about 8.5 GB, copied by hand).

If you only need to read or edit documents and write scripts, steps 1–2 are enough. **Nothing in the analysis chain runs without step 3.**

## 1. Get the code

```sh
git clone git@github.com:jstvdk/gaia-snr-ledger-cygnus.git gaia_snr_history_cygnus
cd gaia_snr_history_cygnus
git log --oneline -3          # expect c1218e4 or newer on main
```

Keep the folder name `gaia_snr_history_cygnus` (some notes and paths use it), but nothing in the code hard-codes the absolute path.

## 2. Create the environment

The project uses a conda environment called `cygob2-gaia` (Python 3.11). Install [Miniforge](https://github.com/conda-forge/miniforge) first if you have no conda.

**Same OS and CPU as the original machine (macOS, Apple silicon):** use the full pinned file. It records the exact builds, so the Monte Carlo results come out the same.

```sh
sed '/^prefix:/d' provenance/environment_cygob2-gaia.yml > /tmp/cygob2-gaia.yml
conda env create -n cygob2-gaia -f /tmp/cygob2-gaia.yml
```

**Linux, Windows or Intel Mac:** the pinned file contains platform-specific builds and will not solve. Use the declared dependencies instead:

```sh
sed '/^prefix:/d' provenance/environment_cygob2-gaia_from-history.yml > /tmp/cygob2-gaia.yml
conda env create -n cygob2-gaia -f /tmp/cygob2-gaia.yml
```

Versions in use on the original machine (from the injection provenance record): Python 3.11.15, numpy 2.4.6, pandas 3.0.3, scipy 1.17.1, scikit-learn 1.9.0. If the new machine resolves different versions, results can differ slightly at the Monte Carlo level. **Do not treat such differences as a finding, and never overwrite an existing product with them.**

Check it:

```sh
conda activate cygob2-gaia
export PYTHONPATH=scripts
python -c "import numpy, pandas, scipy, sklearn, pyarrow, astropy, healpy, dustmaps; print('ok')"
python scripts/wp10_inputs.py       # authorizes and audits the inputs; needs data/
```

Every script is run from the repository root with `PYTHONPATH=scripts`.

## 3. Move the data tree

`data/` is gitignored and has **no archive anywhere else**. The original machine's copy is the only one. Copy it, and do not regenerate it, because the archive queries are not guaranteed to return the same rows today.

| folder | size | content |
|---|---|---|
| `data/raw` | about 2 GB | Gaia, 2MASS, extinction, literature and spectroscopy downloads |
| `data/processed` | about 6.5 GB | every derived product, including the repair_v5, v8 and v9 age and mass products |

There are 1,970 files and none is larger than 1 GB, so a FAT32 drive would work, but use exFAT or APFS anyway. Quit any running job first.

**On the original machine, copy to the external SSD:**

```sh
rsync -a --info=progress2 data/ /Volumes/<SSD>/cygnus_data/
(cd data && find . -type f -print0 | sort -z | xargs -0 shasum -a 256) > /Volumes/<SSD>/cygnus_data.sha256
```

**On the new machine, from the repository root:**

```sh
mkdir -p data
rsync -a --info=progress2 /Volumes/<SSD>/cygnus_data/ data/
(cd data && find . -type f -print0 | sort -z | xargs -0 shasum -a 256) | diff - /Volumes/<SSD>/cygnus_data.sha256 && echo "data copy verified"
```

`diff` prints nothing and the message appears if every file matches. Use `sha256sum` instead of `shasum -a 256` on Linux; the line format is the same.

Also copy these if you want them, since they are not (or not entirely) in git:
- `manuscript/aa.cls` and `manuscript/aa.bst`: the A&A class files, deliberately not committed. Fetch them from <https://www.aanda.org/doc_journal/instructions/macro/aa/macro-latex-aa.zip> or copy them. They are only needed to compile the paper.
- `papers/`: the reference PDFs. Most are tracked in git; the folder is in `.gitignore`, so new files there are not.

Extinction queries use the Bayestar web service through `dustmaps`, so no local dust map needs copying. A re-run of the WP3 extinction step needs internet access.

## 4. What you can run

```sh
export PYTHONPATH=scripts
python scripts/wp10_inputs.py                     # input authorization and audit
python scripts/wp4v9_score.py                     # re-scores Stage 1 (needs data/ and the Stage 1 tables)
bash scripts/run_manuscript_chain.sh              # the whole manuscript chain
```

The manuscript needs [Tectonic](https://tectonic-typesetting.github.io/) (or `latexmk`) to compile: `brew install tectonic` on macOS. The recorded clean build used Tectonic 0.17.0 (`provenance/manuscript_build_execution.json`).

Long runs (hours) should be detached so they survive a closed terminal:

```sh
screen -dmS <name> caffeinate -i <script>     # macOS: caffeinate keeps it awake; on Linux omit it
tail -f <logfile>
```

## 5. Rules (also in CLAUDE.md)

- **Never open `AUDIT.txt`.** It is 16 MB of checksums; regenerate it with `audit.py` if needed.
- Nothing is overwritten and nothing is retuned. Superseded products stay on disk.
- Pre-registrations are written once, committed, and shown to the owner before any run.
- Commit only when the owner asks.

## 6. Where to start working

- Current task and state: [tasks/HANDOFF_repair_v9.md](tasks/HANDOFF_repair_v9.md).
- Project spine: `PROJECT_TRACE.md`. Rules and headline numbers: `CLAUDE.md`.

## 7. If you use Claude Code on the new machine

- Install it from <https://claude.com/claude-code>, open it in the repository root, and it reads `CLAUDE.md` automatically.
- Tell it: "Read tasks/HANDOFF_repair_v9.md and continue."
- Conversation history stays on the machine where it happened. The handoff file replaces it.
