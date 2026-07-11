# DISES Data

## What is the DISES dataset?

DISES (Dynamics of Integrated Socio-Environmental Systems) aggregates farmer survey responses collected in early 2024 about the **2023 growing season**. Respondents were onboarded via stratified random sampling of corn-soy owner-operators in the western Lake Erie Basin (Maumee watershed), across different farm sizes, restricted to farmers with at least one owned field (or contiguous owned fields) so their land could be spatially identified. **645 respondents** participated.

DISES is currently the project's only source of direct farmer characteristics. At the initial stage of our project, the team is using a smaller subset of Regrow/CSB fields with farmer-characteristic data to evaluate how predictive those characteristics are of land management decisions; results will inform whether/how to expand future survey rounds.

## Which area is covered by DISES?

DISES survey responses come specifically from the Maumee watershed / western Lake Erie Basin portion of Indiana, Michigan, and Ohio (`states_DISES` in `config.yml`) — not the entirety of those three states. This is a sub-region of the broader 7-state Regrow/CSB study area.

## What is a "DISES parcel"?

A DISES parcel is a tax-record-derived land unit. During dataset processing, all individual tax parcels belonging to the same survey respondent (`comp_id`) are **dissolved into one combined multi-part geometry** representing that farmer's entire landholding footprint — so the unit that actually gets spatially joined against Regrow/CSB fields is "one farmer's combined holdings," not "one tax parcel."

## What is the difference between a "Regrow field" and a "DISES parcel"?

A Regrow field is the output of a remote-sensing field-delineation *algorithm* based on observed crop-choice sequences; a DISES parcel is based on **tax/ownership records**, collected for a separate project with no direct association to this one. For the same physical land, the two may disagree in: **shape** (the polygon outline itself), **extent** (how much area is covered — a DISES holding may be larger or smaller than the corresponding Regrow fields), and **count** (one farmer's DISES holding may correspond to one Regrow field, several Regrow fields, or only partially overlap any Regrow field at all). This is exactly why the merging procedure in [Merging Regrow with DISES and Supplements](../30-workflow/31-merging-regrow-dises-supplements.md) needs explicit match-quality scoring rather than assuming a clean 1:1 correspondence. Because DISES was built independently, users should expect — and the match-quality variables are designed to surface — inconsistencies between the two datasets.

## What are the key attributes of the DISES dataset?

`[TODO — ADD INFORMATION]`

## Where can I find more information about DISES dataset attributes?

`[TODO — paste links]`

- DISES Codebook:
- DISES Survey Data Overview (slide deck):
