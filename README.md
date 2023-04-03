# flowpipe
Flow & mass cytometry analysis pipeline

## Preliminaries
The *flowpipe* R package is fairly easy to set up. In an R session:
```
install.packages("remotes") # If necessary.
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true") # https://github.com/r-lib/remotes#environment-variables
remotes::install_github("priscian/flowpipe")
library(flowpipe)

## Once the package has been installed as described above, all you need to use it is:
library(flowpipe)
```

Well, that's not quite *all* you need. Please continue reading for more details.

## Table of Contents
1. [Preliminaries](#preliminaries)
1. [Outline of analysis steps](#outline-of-analysis-steps)
1. [Notes](#notes)
1. [References](#references)

***

## Outline of analysis steps

It might be worthwhile here to describe briefly how our analysis pipeline works:

1. Check for common names & descriptions among the FCS files; suggest & perform renaming as necessary
1. Bead normalization (if necessary)
1. Debarcoding/deconvolution (if necessary)
1. IDs neg/dim/pos/bright events for each sample for each channel by silhouette-scanning [[Hu &al 2018]](#hu-et-al-2018)
    * The expression matrix for each FCS file is now accompanied by a "shadow matrix" of phenotypes
1. Create a single expression matrix (& shadow matrix) from all FCS files, IDed by sample
1. Using automatic biaxial gating, gate down to live cells (or further as required by the investigator)
1. Use PhenoGraph-based metaclustering [[Levine &al 2015]](#levine-et-al-2015) to mitigate batch effects & to assign each event to a cluster
    * Procedure: find PhenoGraph clusters for each sample; PhenoGraph-cluster on centroids of these clusters
1. Search for investigator-defined clusters (e.g. T cells, neutrophils, monocytes) & assign events to them
    * Two methods: find & merge PG clusters; biaxial gating down to population
1. Differential abundance analysis by cluster, based on user-supplied clinical metadata
    * Evaluates differential abundance in each cluster using a negative-binomial GLM (generalized linear model), comparing cell counts in each cluster for each condition relative to a control group & producing a log2 fold change
1. These steps are all accompanied by density plots, heatmaps, & UMAPs along the way, & a PDF report

***

## Notes

***

## References

1. <span id="hu-et-al-2018"/> Hu Z, Jujjavarapu C, Hughey JJ, et al. MetaCyto: A tool for automated meta-analysis of mass and flow cytometry data. Cell Rep 24:1377–88, 2018. [dx.doi.org/10.1016/j.celrep.2018.07.003](https://dx.doi.org/10.1016/j.celrep.2018.07.003)

2. <span id="levine-et-al-2015"/> Levine JH, Simonds EF, Bendall SC, et al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162(1):184–97, 2015. [dx.doi.org/10.1016/j.cell.2015.05.047](https://dx.doi.org/10.1016/j.cell.2015.05.047)
