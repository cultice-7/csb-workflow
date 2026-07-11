# Technological Prerequisites

- **Anaconda / Miniconda** — for environment management.

- **Python** — installed and managed via the `csb_workflow` conda environment (see [How to Set Up the Conda Environment](06-how-to-set-environment.md)).

- **Git Bash** (Windows) or any POSIX-compatible shell (macOS/Linux) — for running `git`/command-line steps.

- **Snakemake** — workflow orchestration; every processing step in this repository is a Snakemake rule.

- **Visual Studio Code** (or any editor) — used by the project authors to write and run the code; not strictly required, any Python-capable editor works.
- Core scientific/geospatial stack — all packages are listed in `envs/env.yml` and installed automatically when you create the conda environment.
