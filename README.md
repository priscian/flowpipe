# flowpipe
Flow & mass cytometry analysis pipeline

## Preliminaries
The *flowpipe* R package is fairly easy to set up. In an R session:
```r
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
1. [Annotated example analysis](#annotated-example-analysis)
    * [Header material](#header-material)
    * [Gather FCS files and check](#gather-fcs-files-and-check)
1. [Notes](#notes)
1. [References](#references)

***

## Outline of analysis steps

It might be worthwhile here to describe briefly how our analysis pipeline works. Flow cytometry experiments often conclude with the production of [FCS computer files](https://en.wikipedia.org/wiki/Flow_Cytometry_Standard) containing investigators' raw results. After a flow cytometry experiment, it's not uncommon for investigators to import their FCS files into software such as FlowJo or FCS Express, and then proceed with an analysis "by hand," which can be time-consuming and somewhat subjective.

In the interest of saving time and limiting subjectiveness, the [URMC Flow Cytometry Resource](https://www.urmc.rochester.edu/research/flow-core.aspx) data-analysis team has developed a soup-to-nuts analysis pipeline for the R programming environment that we're calling *flowpipe*; it can semi-automatically handle most analytical tasks from pre-processing/pre-gating to phenotype clustering to differential-expression modeling. The length of a *flowpipe* analysis depends on the number and size of the FCS files provided to the software, but a typical run takes a few hours. Our software has aggregated a number of common techniques and algorithms (well-represented in the peer-reviewed literature) into a flexible parallel-processing framework that's meant to reduce investigators' analytical workload. As part of the *flowpipe* analysis process, investigators can provide "metadata" relevant to their flow-cytometry experiment: that is, for example, whether samples are cases or controls; a list of pre-defined phenotype gates for drilling down to interesting cell subsets; and patient data that can be incorporated into the differential-expression models.

Here's the outline of a typical *flowpipe* workflow:

1. Check for common names & descriptions among the FCS files; suggest & perform renaming as necessary
1. Bead normalization (if necessary)
1. Debarcoding/deconvolution (if necessary)
1. IDs neg/dim/pos/bright events for each sample for each channel by silhouette-scanning [[Hu &al 2018](#hu-et-al-2018)]
    * The expression matrix for each FCS file is now accompanied by a "shadow matrix" of phenotypes
1. Create a single expression matrix (& shadow matrix) from all FCS files, IDed by sample
1. Using automatic biaxial gating, gate down to live cells (or further as required by the investigator)
1. Use PhenoGraph-based metaclustering [[Levine &al 2015](#levine-et-al-2015)] to mitigate batch effects & to assign each event to a cluster
    * Procedure: find PhenoGraph clusters for each sample; PhenoGraph-cluster on centroids of these clusters
1. Search for investigator-defined clusters (e.g. T cells, neutrophils, monocytes) & assign events to them
    * Two methods: find & merge PG clusters; biaxial gating down to population
1. Differential abundance analysis by cluster, based on user-supplied clinical metadata
    * Evaluates differential abundance in each cluster using a negative-binomial GLM (generalized linear model), comparing cell counts in each cluster for each condition relative to a control group & producing a log₂ fold change
1. These steps are all accompanied by density plots, heatmaps, & UMAPs along the way, & a PDF report

***

## Annotated example analysis

This section will provide a detailed walkthrough of a *flowpipe* analysis of FCS files from a mass cytometry (CyTOF) experiment. The analyst should be familiar enough with the R programming environment to call functions, to understand the variety of R data types and objects, to make changes in this walkthrough appropriate to their own experiment, and sometimes to provide ad hoc code snippets for bridging parts of a *flowpipe* analysis together. These requirements shouldn't be onerous; *flowpipe* is designed robustly to do most of the common tedious and diffcult tasks associated with cytometry data analysis.

We'll use data stored on the [FlowRepository](https://flowrepository.org/) [[Spidlen &al 2012](#spidlen-et-al-2012)] from the very instructive paper "A Beginner's Guide To Analyzing and Visualizing Mass Cytometry Data" [[Kimball &al 2018](#kimball-et-al-2018)]. For convenience, we've provided this data in a single zipped archive [here](https://dl.dropboxusercontent.com/s/wd6g2ffstza8oc4/FlowRepository_FR-FCM-ZYDW_files.zip); if you wish to reproduce the walkthrough, download the archive, unzip it into a single directory[[*](#note-asterisk)], and keep track of where it's stored on your file system.

Let's step through the code and discuss it one section at a time; we'll provide the complete code file afterwards.

***

### Header material

This is a code that's common to all *flowpipe* analyses.

```r
## source("./kimball-&al-2018.R", keep.source = FALSE)
rm(list = ls(all.names = FALSE))

options(keystone_parallel = TRUE)

# setwd("your/preferred/working/directory")
data_dir <- "./data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
.cm <- memoise::cache_filesystem(path = data_dir)

if (!require("flowpipe")) {
  remotes::install_github("priscian/flowpipe")
  library(flowpipe)
}

report_dir <- "./report"
image_dir <- paste(report_dir, "images", sep = "/")
```

### Gather FCS files and check

```r
## Find paths of all FCS files in the given directory(-ies)
fcs_files <- keystone::list_files("C:/Users/priscian/Downloads/cytometry/FlowRepository_FR-FCM-ZYDW_files",
  "\\.fcs$", recursive = FALSE, ignore.case = TRUE, full.names = TRUE)

## Check channel names & descriptions among FCS files
## Add argument 'clear_1_cache = TRUE' to clear cache & rerun
cbs <- get_channels_by_sample(x = fcs_files, SUFFIX = "_cbs", clear_1_cache = FALSE)
```

***

## Notes

\*&nbsp; <span id="note-asterisk"> Though *flowpipe* is flexible enough not to require that all analysis FCS files be in a single directory, it's just easier here.

***

## References

1. <span id="hu-et-al-2018"/> Hu Z, Jujjavarapu C, Hughey JJ, et al. MetaCyto: A tool for automated meta-analysis of mass and flow cytometry data. Cell Rep 24:1377–88, 2018. [dx.doi.org/10.1016/j.celrep.2018.07.003](https://dx.doi.org/10.1016/j.celrep.2018.07.003)

1. <span id="levine-et-al-2015"/> Levine JH, Simonds EF, Bendall SC, et al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162(1):184–97, 2015. [dx.doi.org/10.1016/j.cell.2015.05.047](https://dx.doi.org/10.1016/j.cell.2015.05.047)

1. <span id="spidlen-et-al-2012"/> Spidlen J, Breuer K, Rosenberg C, Kotecha N, Brinkman RR. FlowRepository: A resource of annotated flow cytometry datasets associated with peer‐reviewed publications. Cytometry Part A 81(9):727–31, 2012. [dx.doi.org/10.1002/cyto.a.22106](https://dx.doi.org/10.1002/cyto.a.22106)

1. <span id="kimball-et-al-2018"/> Kimball AK, Oko LM, Bullock BL, Nemenoff RA, van Dyk LF, Clambey ET. A beginner's guide to analyzing and visualizing mass cytometry data. J Immunol 200(1):3–22, 2018. [dx.doi.org/10.4049/jimmunol.1701494](https://dx.doi.org/10.4049/jimmunol.1701494)
