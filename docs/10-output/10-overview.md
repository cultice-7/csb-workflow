# Output Datasets

This subfolder describes the structure of all files produced by the pipeline and explains how to run Snakemake to generate them. The pipeline produces two parallel families of output — one built around Regrow fields, one built around CSB fields — plus supplementary data tables that are independently joined to each. Within each family, most merge operations produce their own standalone output file rather than one monolithic table, so users can pull in only the supplementary blocks they need.

Supplementary data blocks are referenced by numeric code (1–9) throughout both Regrow- and CSB-based outputs. The full code table appears in both [11 — Regrow-Based Output](11-regrow-based-output.md) and [12 — CSB-Based Output](12-csb-based-output.md).

## Contents

| File | Question it answers |
|---|---|
| [11 — Regrow-Based Output](11-regrow-based-output.md) | What output files does the Regrow branch produce, and what is in each? What are the supplementary data codes? |
| [12 — CSB-Based Output](12-csb-based-output.md) | What output files does the CSB branch produce, and what are the supplementary data codes? |
| [13 — Data Schema Links](13-data-schema-links.md) | Where are the codebooks, schema references, and dataset overview slides? |
| [14 — How to Obtain Output Datasets](14-how-to-obtain-output-datasets.md) | How do I run Snakemake to build all or selected outputs? What are the recommended execution chains? |
