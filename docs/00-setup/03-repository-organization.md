# Repository Organization

| Item | What it is |
| --- | --- |
| `data/` | Raw input data and all pipeline-generated outputs |
| `docs/` | This documentation |
| `envs/` | Conda environment specification (`env.yml`) needed to run the pipeline |
| `figures/` | Figures used in documentation |
| `scripts/` | All active pipeline scripts; `scripts/archive/` holds superseded versions kept for reference |
| | |
| `Snakefile` | Workflow definition — every pipeline rule with its inputs, outputs, parameters, and script path |
| `config.yml` | Global configuration: study states, file path templates, and algorithm parameters referenced throughout the Snakefile |
| `Scripts order to run the code.docx` | Memo listing the recommended execution order for pipeline stages |
| `LICENSE` | Project license |
