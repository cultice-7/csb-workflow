# How to Set Up the Conda Environment

From an Anaconda Prompt (or any shell with `conda` on the PATH), with your working directory set to the cloned `csb-workflow` folder:

```bash
# Navigate to the cloned repository folder
cd /path/to/your/csb-workflow

# Create the environment from the provided spec (first-time setup only)
conda env create -f envs/env.yml

# Activate the environment
conda activate csb_workflow
```

If you've already created the environment previously and `envs/env.yml` has since changed, update it instead of recreating it:

```bash
conda env update -f envs/env.yml --prune
```

The environment name is `csb_workflow` (defined in `envs/env.yml`). All Snakemake commands in this documentation assume this environment is active.
