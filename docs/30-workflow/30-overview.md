# Workflow: Merging Datasets

This subfolder explains how the raw Regrow and DISES datasets are cleaned and prepared, how Regrow and CSB fields are spatially joined with DISES farmer survey data and with the 9 supplementary farmland characteristic blocks. It covers the cleaning procedures, merge methodology, match-quality scoring, representative-field assignment, and how to combine the standalone output tables into a single analysis-ready dataset.

## Contents

### [31 — Cleaning and Processing Regrow Data](31-cleaning-regrow-data.md)

- [How is 2025 Regrow data combined with the 2014-2024 data?](31-cleaning-regrow-data.md#how-is-2025-regrow-data-combined-with-the-2014-2024-data)
- [How is Regrow field geometry extracted and cleaned?](31-cleaning-regrow-data.md#how-is-regrow-field-geometry-extracted-and-cleaned)
- [Transformation of the Regrow land management variables](31-cleaning-regrow-data.md#transformation-of-the-regrow-land-management-variables)
- [How are raw and processed Regrow datasets structured? What is a cultivation cycle?](31-cleaning-regrow-data.md#how-are-raw-and-processed-regrow-datasets-structured-what-is-a-cultivation-cycle)
- [How are the Regrow geometry and attribute tables joined?](31-cleaning-regrow-data.md#how-are-the-regrow-geometry-and-attribute-tables-joined)

### [32 — Cleaning DISES Data](32-cleaning-dises-data.md)

- [How is the raw DISES survey table cleaned?](32-cleaning-dises-data.md#how-is-the-raw-dises-survey-table-cleaned)
- [How is the raw DISES parcel shapefile cleaned and consolidated?](32-cleaning-dises-data.md#how-is-the-raw-dises-parcel-shapefile-cleaned-and-consolidated)

### [33 — Merging Regrow with DISES and Supplements](33-merging-regrow-dises-supplements.md)

- [How are Regrow and DISES spatially merged?](33-merging-regrow-dises-supplements.md#how-are-regrow-and-dises-spatially-merged)
- [What are the outcomes of the matching?](33-merging-regrow-dises-supplements.md#what-are-the-outcomes-of-the-matching)
- [How is matching quality evaluated?](33-merging-regrow-dises-supplements.md#how-is-matching-quality-evaluated)
- [What is a "representative field" attribute? Why do we need it?](33-merging-regrow-dises-supplements.md#what-is-a-representative-field-attribute-why-do-we-need-it)
- [What are the outcomes of the matching quality procedure?](33-merging-regrow-dises-supplements.md#what-are-the-outcomes-of-the-matching-quality-procedure)
- [How to generate Regrow sub-samples with available DISES data](33-merging-regrow-dises-supplements.md#how-to-generate-regrow-sub-samples-with-available-dises-data)
- [How to extract DISES-related variables](33-merging-regrow-dises-supplements.md#how-to-extract-dises-related-variables)
- [Merging Regrow with each supplementary dataset](33-merging-regrow-dises-supplements.md#merging-regrow-with-each-supplementary-dataset)
- [Output dataset structure and integration](33-merging-regrow-dises-supplements.md#output-dataset-structure-and-integration)

### [34 — Merging CSB with DISES and Supplements](34-merging-csb-dises-supplements.md)

- [How are CSB and DISES spatially merged?](34-merging-csb-dises-supplements.md#how-are-csb-and-dises-spatially-merged)
- [What are the outcomes of the matching?](34-merging-csb-dises-supplements.md#what-are-the-outcomes-of-the-matching)
- [How is matching quality evaluated?](34-merging-csb-dises-supplements.md#how-is-matching-quality-evaluated)
- [What is a representative-field attribute on the CSB side?](34-merging-csb-dises-supplements.md#what-is-a-representative-field-attribute-on-the-csb-side)
- [What are the outcomes of the matching quality procedure?](34-merging-csb-dises-supplements.md#what-are-the-outcomes-of-the-matching-quality-procedure)
- [How to generate CSB sub-samples with available DISES data](34-merging-csb-dises-supplements.md#how-to-generate-csb-sub-samples-with-available-dises-data)
- [How to extract DISES-related variables](34-merging-csb-dises-supplements.md#how-to-extract-dises-related-variables)
- [Merging CSB with each supplementary dataset](34-merging-csb-dises-supplements.md#merging-csb-with-each-supplementary-dataset)
- [Output dataset structure and integration](34-merging-csb-dises-supplements.md#output-dataset-structure-and-integration)
