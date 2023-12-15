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
    * [Gather FCS files and check them](#gather-fcs-files-and-check-them)
    * [Perform bead normalization and\/or debarcoding](#perform-bead-normalization-andor-debarcoding)
    * [Check channel names and descriptions among FCS files (again)](#check-channel-names-and-descriptions-among-fcs-files-again)
    * [Describe pre-gating channels](#describe-pre-gating-channels)
    * [Density plots](#density-plots)
    * [Find plus/minus events for each sample for each channel](#find-plusminus-events-for-each-sample-for-each-channel)
    * [Sample from and combine expression matrices](#sample-from-and-combine-expression-matrices)
    * [Metadata](#metadata)
    * [Plot sample cell counts](#plot-sample-cell-counts)
    * [Clustering](#clustering)
    * [Select channels to use in further analyses](#select-channels-to-use-in-further-analyses)
    * [Search for investigator-defined clusters](#search-for-investigator-defined-clusters)
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
    * Evaluates differential abundance in each cluster using a negative-binomial GLM (generalized linear model), comparing cell counts in each cluster for each condition relative to a control group & producing a [log₂ fold change](https://stackoverflow.com/questions/70696602/how-to-interpret-log-fold-change-log2fc-on-two-cases/70818010#70818010)
1. These steps are all accompanied by density plots, heatmaps, & UMAPs along the way, & a PDF report

***

## Annotated example analysis

This section will provide a detailed walkthrough of a *flowpipe* analysis of FCS files from a mass cytometry (CyTOF) experiment. The analyst should be familiar enough with the R programming environment to call functions, to understand the variety of R data types and objects, to make changes in this walkthrough appropriate to their own experiment, and sometimes to provide ad hoc code snippets for bridging parts of a *flowpipe* analysis together. These requirements shouldn't be onerous; *flowpipe* is designed robustly to do most of the common tedious and diffcult tasks associated with cytometry data analysis.

We'll use data stored on the [FlowRepository](https://flowrepository.org/) [[Spidlen &al 2012](#spidlen-et-al-2012)] from the very instructive paper "A Beginner's Guide To Analyzing and Visualizing Mass Cytometry Data" [[Kimball &al 2018](#kimball-et-al-2018)]. For convenience, we've provided this data in a single zipped archive [here](https://dl.dropboxusercontent.com/s/wd6g2ffstza8oc4/FlowRepository_FR-FCM-ZYDW_files.zip); if you wish to reproduce the walkthrough, download the archive, unzip it into a single directory[[*](#note-asterisk)], and keep track of where it's stored on your file system.

Let's step through the code and discuss it one section at a time; we'll provide the complete code file afterwards.

***

### Header material

This is code that's common to all *flowpipe* analyses.

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

### Gather FCS files and check them

```r
## Find paths of all FCS files in the given directory(-ies)
fcs_files <- keystone::list_files("C:/Users/priscian/Downloads/cytometry/FlowRepository_FR-FCM-ZYDW_files",
  "\\.fcs$", recursive = FALSE, ignore.case = TRUE, full.names = TRUE)

## Check channel names & descriptions among FCS files
## Add argument 'clear_1_cache = TRUE' to clear cache & rerun
cbs <- get_channels_by_sample(x = fcs_files, SUFFIX = "_cbs", clear_1_cache = FALSE)
```

The function `get_channels_by_sample()` compiles channel names and descriptions from the experiment's FCS files into a tabular `data.frame`; you can check this table to insure that channels are labeled correctly among all the samples. For the [Kimball &al 2018](#kimball-et-al-2018) data we're working on here, this is how the `cbs` variable looks:

```r
> cbs %>% tibble::as_tibble()
# A tibble: 62 × 10
   name         desc_01  desc_02  desc_03  desc_04  desc_05  desc_06  desc_07  desc_08  desc_09
   <chr>        <chr>    <chr>    <chr>    <chr>    <chr>    <chr>    <chr>    <chr>    <chr>
 1 Time         <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>
 2 Event_length <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>     <NA>
 3 Y89Di        89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45 89Y_CD45
 4 Pd102Di      102Pd    102Pd    102Pd    102Pd    102Pd    102Pd    102Pd    102Pd    102Pd
 5 Rh103Di      103Rh    103Rh    103Rh    103Rh    103Rh    103Rh    103Rh    103Rh    103Rh
 6 Pd104Di      104Pd    104Pd    104Pd    104Pd    104Pd    104Pd    104Pd    104Pd    104Pd
 7 Pd105Di      105Pd    105Pd    105Pd    105Pd    105Pd    105Pd    105Pd    105Pd    105Pd
 8 Pd106Di      106Pd    106Pd    106Pd    106Pd    106Pd    106Pd    106Pd    106Pd    106Pd
 9 Pd108Di      108Pd    108Pd    108Pd    108Pd    108Pd    108Pd    108Pd    108Pd    108Pd
10 Pd110Di      110Pd    110Pd    110Pd    110Pd    110Pd    110Pd    110Pd    110Pd    110Pd
# ℹ 52 more rows
# ℹ Use `print(n = ...)` to see more rows
```

The `clear_1_cache` argument, which you see frequently in these function calls, is available for functions whose return value—in addition to being returned to memory—is cached to disk. When the same function is called again, the return value will be retrieved from the disk cache instead of resulting from the function's operation; this caching can save a lot of time by skipping over time-consuming operations that have already been completed once. The usual value is `clear_1_cache = FALSE`; use `clear_1_cache = TRUE` to rerun the function and update the disk cache. (May be necessary if other function arguments have changed since the last call.)

### Perform bead normalization and\/or debarcoding

The function `bead_normalize()` provides a pipeline-friendly interface to `CATALYST::normCytof()` [[Crowell &al 2023](#crowell-et-al-2023)] function. (V. [Finck &al 2013](#finck-et-al-2013) for a discussion of bead-normalization techniques and algorithms in mass cytometry.) `bead_normalize()` returns a vector of file paths to the newly normalized set of FCS files.

```r
## 'fcs_files_01' is a vector of file paths to the normalized FCS files
fcs_files_01 <- bead_normalize(
  input_path = fcs_files,
  output_dir = paste(data_dir, "fcs/bead-normalized", sep = "/"),
  outfile_suffix = "-bead-normalized",
  ## "dvs" (for bead masses 140, 151, 153, 165, 175) or
  ##   "beta" (for masses 139, 141, 159, 169, 175) or numeric vector of masses:
  beads = "dvs",
  normCytof... = list(plot = FALSE),
  SUFFIX = "_beadnorm", clear_1_cache = FALSE
)
```

### Check channel names and descriptions among FCS files (again)

Run the function `get_channels_by_sample()` again after all normalization or debarcoding/deconvolution steps that produce new FCS files. This step insures that the channel names and descriptions are common among all the files. (If it throws a warning, use the function `rename_fcs_parameters_name_desc()` to fix the problem semi-automatically; type `help("rename_fcs_parameters_name_desc")` for details.)

```r
cbs_check <- get_channels_by_sample(x = fcs_files_01, SUFFIX = "_cbs-check", clear_1_cache = FALSE)

## Throw a warning if these parameters are still different among the files:
if (any_noncommon_descs(cbs_check))
  warning("FCS files may still not share common channel descriptions", immediate. = TRUE)
```

### Describe pre-gating channels

Variables called `𝘹𝘹𝘹𝘹_channels` in this walkthrough are so-called [named vectors](https://devtut.github.io/r/creating-vectors.html#creating-named-vectors) in R: the vector *values* are the channel names provided in each FCS file, and the vector *names* are the channel descriptions. For example, we're using here in the walkthrough a variable called `pregating_channels`, which looks like this—

```r
> pregating_channels
             <NA>              <NA>          89Y_CD45        141Pr_Gr-1       142Nd_CD11c
           "Time"    "Event_length"           "Y89Di"         "Pr141Di"         "Nd142Di"
       143Nd_GITR        144Nd_MHCI        145Nd_CD69         146Nd_CD8       148Nd_CD11b
        "Nd143Di"         "Nd144Di"         "Nd145Di"         "Nd146Di"         "Nd148Di"
    149Sm_p4E-BP1        150Nd_CD25         152Sm_CD3     154Sm_CTLA4-i        155Gd_IRF4
        "Sm149Di"         "Nd150Di"         "Sm152Di"         "Sm154Di"         "Gd155Di"
    156Gd_SigF-PE       158Gd_FoxP3        159Tb_PD-1  160Gd_KLRG1-FITC        161Dy_Tbet
        "Gd156Di"         "Gd158Di"         "Tb159Di"         "Gd160Di"         "Dy161Di"
       162Dy_Tim3    163Dy_LAG3-APC        164Dy_IkBa        166Er_CD19       167Er_NKp46
        "Dy162Di"         "Dy163Di"         "Dy164Di"         "Er166Di"         "Er167Di"
       168Er_Ki67        169Tm_Sca1 170Er_PDL2-Biotin        171Yb_CD44         172Yb_CD4
        "Er168Di"         "Tm169Di"         "Er170Di"         "Yb171Di"         "Yb172Di"
      173Yb_CD117       174Yb_MHCII        176Yb_ICOS
        "Yb173Di"         "Yb174Di"         "Yb176Di"
```

—and where the first line gives the vector names, and the second the vector values, all the way down. One way to create such an object is to concatenate name-value pairs "by hand":

```r
pregating_channels <- c("Time", "Event Length", 89Y_CD45 = "Y89Di", `141Pr_Gr-1` = "Pr141Di")
## Etc.
```

However, we've created `pregating_channels` through a more baroque method:

```r
pregating_channels <- cbs_check %>%
  dplyr::filter(
    (!is.na(desc_01) | stringr::str_detect(name, stringr::regex("time|event_length",
      ignore_case = TRUE))) &
    stringr::str_detect(attr(., "channels_by_sample_full_desc")$desc_01,
      ## The 1st part of the RE here, negative lookahead, insures that "Time" & channels containing "_" are kept:
      stringr::regex(
        "^(?!.*(_|time))|beads|_eq|barcode|bkg|bckg|background|intercalator|viability|live-dead",
        ignore_case = TRUE), negate = TRUE)) %>%
    (function(o) structure(o$name, .Names = o$desc_01))(.)
```

If that seems like too much effort to understand, you can review the variable `cbs_check` and build the vector "by hand" as described above. `pregating_channels` will be used to plot channel densities by sample, then a subset of the vector will be used for further analysis.

### Density plots

Per-sample per-channel density plots are diagnostic tools used primarily for insuring that no overwhelming batch effect informs the experiment data.

```r
## Density plots
plot_channel_densities_by_sample(
  x = fcs_files_01,
  image_dir = image_dir,
  current_image = 1,
  save_plot = TRUE,
  get_fcs_expression_subset... = list(
    b = 1/8,
    channels_subset = pregating_channels,
    channels_by_sample = cbs_check
  ),
  SUFFIX = "_channel-densities", clear_1_cache = FALSE
)
```

![Channel distribution for each sample (detail).](<inst/images/001_channel_dist-detail.png>)
**Figure.** Channel distribution for each sample (detail). The full image can be found as a PNG or zoomable PDF file in the results archive for this analysis, a large (~ 8 GB) ZIP file that can be downloaded from [here](https://www.dropbox.com/scl/fi/t2okywora7mszsccpmt2p/kimball-al-2018.zip).

### Find plus/minus events for each sample for each channel

This step identifies neg/dim/pos/bright events for each sample for each channel by silhouette-scanning [[Hu &al 2018](#hu-et-al-2018)]. The expression matrix for each FCS file is now accompanied by a "shadow matrix" of plus-minus (`-, d, +, ++`) phenotypes, which helps in setting predefined gates later.

```r
## Make augmented flowCore::flowFrame objects.
pmm_files <- prepare_augmented_fcs_data(
  x = fcs_files_01,
  b = 1/8,
  data_dir = paste(data_dir, "samples", sep = "/"),
  remove_outliers = TRUE,
  pmm_channels = cbs_check %>%
    dplyr::filter(stringr::str_detect(name, stringr::regex("time|event_length",
      ignore_case = TRUE)) |
    !is.na(desc_01)) %>% dplyr::pull(name),
  multisect... = list(max_sample = 5000),
  outfile_suffix = NULL,
  outfile_prefix = expression(outfile_prefix <- paste0(x %>% dirname %>% basename, "-")),
  SUFFIX = "_pmm-files", clear_1_cache = FALSE
)
```

### Sample from and combine expression matrices

Now we move on to the primary analyses, which benefit from aggregating all the FCS files/samples (or a downsampled subset of each files's events) into a single expression matrix. The gating-strategy list is used to remove unwanted events (e.g. dead cells, doublets) from each sample; this so-called pregating is summarized in PNG & PDF posters showing the strategies applied to all samples, [like so](inst/images/sample-pregating.png).

The function argument `sample_size = Inf` creates a complete aggregate expression matrix without downsampling any events. Use e.g. `sample_size = 5000` to randomly select a maximum of 5000 events from each sample.

```r
## Get stacked subsets of augmented expression matrices IDed by sample.
e <- get_expression_subset(
  x = pmm_files,
  gate... = list(
    strategy = list(
      `intact cells` = expression(
        x[, "Ir191Di"] %>% dplyr::between(quantile(., 0.25), quantile(., 0.95)) &
        x[, "Ir193Di"] %>% dplyr::between(quantile(., 0.25), quantile(., 0.95))
      ),
      `intact singlets` = expression(
        x[, "Event_length"] %>% { . < quantile(., 0.90) } &
        x[, "Ir191Di"] %>% { rep(TRUE, length(.)) }
      ),
      `live cells` = expression(
        x[, "Pt195Di"] %>% { . < quantile(., 0.90) } &
        x[, "Ir191Di"] %>% { rep(TRUE, length(.)) }
      )
    )
  ),
  #save_plot_fun = grDevices::x11, save_plot... = list(bg = "transparent"), # Plot to screen
  #save_plot... = NULL, # Don't plot to any device
  ## 'file' must be given to save PDF:
  save_plot... = list(file = paste(image_dir, "sample-pregating.pdf", sep = "/")),
  sample_size = Inf,
  SUFFIX = "_expression-subset", clear_1_cache = FALSE
)
```

### Metadata

A metadata data frame is used to identify samples as being in different experimental groups, and to provide additional sample-related covariates also used for grouping, or for adjusting statistical analyses so that the independent effect of the variables of interest is as precise as possible.

```r
## Create or import sample metadata here. Must minimally have 'id' & 'group' columns.
## This should be adapted specifically to each study, i.e. there's no default.
metadata <- read.csv(text = '
id,group
Clambey LO 110116 B6 lung1_01,B6
Clambey LO 110116 B6 lung 2_01,B6
Clambey LO 110116 B6 lung 3_01,B6
Clambey LO 11022016 B6 lung 4_01,B6
Clambey LO 11022016 B6 lung 5_01,B6
Clambey LO 110116 IL10KO lung 1_01,IL10KO
Clambey LO 110116 IL10KO lung 2_01,IL10KO
Calmbey LO 11022016 IL10 KO lung 3_01,IL10KO
Clambey LO 11022016 IL10 KO lung 4_01,IL10KO', comment.char = "%") %>%
  dplyr::mutate(group = factor(group)) -> # N.B. May have to reorder levels here, e.g. via package "forcats"
  metadata_i

## RE to 'stringr::str_extract()' concise sample ID from 'id_map' attribute in "pmm" object from 'get_expression_subset()'
sample_name_re <- "C(al|la)mbey.*?_\\d{2}"
## sort(stringr::str_extract(attr(e, "id_map") %>% names, sample_name_re)) == sort(metadata$id) # Should all be TRUE
```

Our metadata here is created by the code above, but it could just as easily be imported from a spreadsheet or CSV file at this point.

```r
> metadata
                                     id  group
1         Clambey LO 110116 B6 lung1_01     B6
2        Clambey LO 110116 B6 lung 2_01     B6
3        Clambey LO 110116 B6 lung 3_01     B6
4      Clambey LO 11022016 B6 lung 4_01     B6
5      Clambey LO 11022016 B6 lung 5_01     B6
6    Clambey LO 110116 IL10KO lung 1_01 IL10KO
7    Clambey LO 110116 IL10KO lung 2_01 IL10KO
8 Calmbey LO 11022016 IL10 KO lung 3_01 IL10KO
9 Clambey LO 11022016 IL10 KO lung 4_01 IL10KO
```

### Plot sample cell counts

Before and after pregating. The plots make use of the previously provided metadata.

![Starting event counts for each sample.](<inst/images/002_counts.png>)

![Event counts for each sample after pregating.](<inst/images/002_counts_gated.png>)

### Clustering

PhenoGraph [[Levine &al 2015](#levine-et-al-2015)] clusters expression data into subpopulations of events that are phenotypically similar in high-dimensional space. (More detailed, via Liu &al 2019 [[Liu &al 2019](#liu-et-al-2019)], PhenoGraph uses a k-nearest neighbors (KNN) technique to detect connectivity and density peaks among high-dimensional events.) In its first presentation in the literature [[Levine &al 2015](#levine-et-al-2015)], Phenograph is also used for batch correction; the procedure is to find PhenoGraph clusters for each sample, then to PhenoGraph-cluster again on the centroids of those clusters, to finally produce *metaclusters* to which each event can be assigned. This is an inexpensive batch-correction solution which doesn't require control samples, but which provides some measure of relief from batch effects.

We've come to prefer PhenoGraph clustering because a\) it's fairly conservative in creating clusters (i.e. you don't typically end up with, say, 400 of them); and b\) the metaclustering technique described in [Levine &al 2015](#levine-et-al-2015) provides low-cost batch correction. X-shift clustering [[Samusik &al 2016](#samusik-et-al-2016)] is good too, but less conservative and more difficult to use. In general, we prefer applying a clustering technique that assigns every event to a cluster (unlike e.g. Citrus [[Bruggner &al 2014](#bruggner-et-al-2014)]); FlowSOM [[Van Gassen &al 2015](#van-gassen-et-al-2015)] is perhaps the best tool in this regard, but its requirement that the user provide the final number of clusters perhaps limits its use in automatic analysis.

However … sometimes the PhenoGraph metaclustering results in too *few* clusters, as is the case with the relatively small number of samples (9) here, so we've instead asked FlowSOM to give us 25 clusters, and we'll proceed from there.

```r
## Clustering
clustering_channels <-
  pregating_channels[stringr::str_subset(pregating_channels %>% names, "_")]

cluster_id_mc <- make_metaclusters(
  x = e,
  channels = clustering_channels,
  make_clusters... = list(method = "Rphenograph", Rphenograph_k = 50),
  make_metaclusters... = list(method = "Rphenograph", Rphenograph_k = 15),
  #centroid_fun = mean,
  SUFFIX = "_cluster-id-mc", clear_1_cache = FALSE
)

cluster_id <- make_clusters(
  x = e,
  method = "FlowSOM", FlowSOM_k = 25, estimate_cluster_count = FALSE,
  #method = "Rphenograph", Rphenograph_k = 50,
  channels = clustering_channels,
  SUFFIX = "_cluster-id", clear_1_cache = FALSE
)

attr(e, "cluster_id") <- cluster_id
```

### Select channels to use in further analyses

After the clustering step, choose channels of primary interest for further analysis.

```r
## Change names of analysis channels to reflect antibodies (custom for each study)
analysis_channels <- clustering_channels %>%
  structure(names(.) %>% stringr::str_match(".*?_(.*)") %>% `[`(, 2), .Names = .) %>%
  `[<-`(is.na(.), names(.)[is.na(.)]) %>% c(Eu151Di = "CD64") # Include bead channel as CD64
## N.B. How to easily rename some column names w/ maps:
e_channels_abo <- structure(colnames(e), .Names = colnames(e))
colnames(e) <- e_channels_abo %>% `[<-`(names(analysis_channels), analysis_channels) %>%
  as.vector
```

### Search for investigator-defined clusters

[Recall](#find-plusminus-events-for-each-sample-for-each-channel) that events in the aggregate expression matrix have already been classified as one of the four plus-minus phenotypes (`-, d, +, ++`). These classifications allow investigators to define cell subtypes of interest based on channel description and plus-minus phenotype. Because there are four phenotypes, and because we often want a simpler binary classification, phenotypes can be combined into an inclusive binary *plus* as `+/++`, and *minus* as `-/d`; for example, CD45-*plus* is defined in code as `"cd45+/++"`, and CD64-*minus* as `"cd64-/d"`.

```r
## User-defined cell subtypes of interest
cell_subtypes <- list(
  `B cells` = c("cd45+/++", "cd19+/++"),
  `CD4+ T cells` = c("cd45+/++", "cd3+/++", "cd4+/++"),
  `CD8+ T cells` = c("cd45+/++", "cd3+/++", "cd8+/++"),
  `Natural killer cells` = c("cd45+/++", "nkp46+/++"),
  `CD11c+ Myeloid cells` = c("cd11c+/++"),
  `CD11b+ CD64+ Myeloid cells` = c("cd11b+/++", "cd64+/++"),
  `CD11b+ CD64- Myeloid cells` = c("cd11b+/++", "cd64-/d"),
  `Hematopoietic cells` = c("cd45+/++", "cd3-/d", "nkp46-/d", "cd19-/d", "cd4-/d",
    "cd8-/d", "cd11b-/d", "cd11c-/d", "cd64-/d"),
  `Non-hematopoietic cells` = c("cd45-/d")
)

## Search clusters for user-defined cell types & merge them if necessary
merged_clusters <- merge_clusters(
  x = e,#[, analysis_channels],
  clusters = cell_subtypes,
  channels = analysis_channels,
  label_threshold = 0.55, # Default is 0.90
  #which_cluster_set = NULL, # Search every event for cluster definitions
  SUFFIX = "_merged-clusters", clear_1_cache = FALSE
)
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

1. <span id="crowell-et-al-2023"/> Crowell H, Zanotelli V, Chevrier S, Robinson M. CATALYST: Cytometry dATa anALYSis Tools. R package version 1.26.0, 2023. [dx.doi.org/10.18129/B9.bioc.CATALYST](https://dx.doi.org/10.18129/B9.bioc.CATALYST)

1. <span id="finck-et-al-2013"/> Finck R, Simonds EF, Jager A, Krishnaswamy S, Sachs K, Fantl W, Pe'er D, Nolan GP, Bendall SC. Normalization of mass cytometry data with bead standards. Cytometry Part A 83(5):483–94, 2013. [dx.doi.org/10.1002/cyto.a.22271](https://dx.doi.org/10.1002/cyto.a.22271)

1. <span id="liu-et-al-2019"/> Liu X, Song W, Wong BY, et al. A comparison framework and guideline of clustering methods for mass cytometry data. Genome Biol 20(1):1–8, 2019. [dx.doi.org/10.1186/s13059-019-1917-7](https://dx.doi.org/10.1186/s13059-019-1917-7)

1. <span id="samusik-et-al-2016"/> Samusik N, Good Z, Spitzer MH, Davis KL, Nolan GP. Automated mapping of phenotype space with single-cell data. Nat Methods 13:493–96, 2016. [dx.doi.org/10.1038/nmeth.3863](https://dx.doi.org/10.1038/nmeth.3863)

1. <span id="bruggner-et-al-2014"/> Bruggner RV, Bodenmiller B, Dill DL, Tibshirani RJ, Nolan GP. Automated identification of stratifying signatures in cellular subpopulations. PNAS 111(26):E2770–E2777, 2014. [dx.doi.org/10.1073/pnas.1408792111](https://dx.doi.org/10.1073/pnas.1408792111)

1. <span id="van-gassen-et-al-2015"/> Van Gassen S, Callebaut B, Van Helden MJ, et al. FlowSOM: Using self‐organizing maps for visualization and interpretation of cytometry data. Cytometry Part A 87(7):636–45, 2015. [dx.doi.org/10.1002/cyto.a.22625](https://dx.doi.org/10.1002/cyto.a.22625)
