## http://cytoforum.stanford.edu/viewtopic.php?f=3&t=1026
##
## 1. Acquisition of raw data (hot off the machine, no additional processing)
## 2. Normalization. Either MATLAB/Finck style, or Fluidigm-style.
## 3. Debarcoding, if applicable.
## 4. Gating down to live intact singlets (removal of debris, beads, dead cells,
##   remaining cell-cell and bead-cell doublets).
## 5. Depending on your analysis pipeline, often export of new FCS files from the gated live intact singlets
##   (if using Cytobank or a FlowJo plug-in, you can usually specify the population to cluster without having to re-export).
## 6. Clustering.


### Apply silhouette cutoffs to flow data.

#' @export
find_plus_minus_by_channel <- function(
  x, # flowCore::flowFrame
  bins = 4, intensity_levels = c("-", "d", "+", "++"),
  zero_threshold = 0.90,
  excluded_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  multisect... = list(),
  verbose = TRUE
)
{
  temp <- flowCore::colnames(x); channels <- temp[stringr::str_detect(temp, excluded_channels_re, negate = TRUE)]

  e <- flowCore::exprs(x)
  colNames <- colnames(e)
  cutoffs <- list(); i <- 1
  pmm <- plyr::alply(e, 2,
    function(b)
    {
      sampleName <- tools::file_path_sans_ext(basename(flowCore::description(x)$FILENAME))
      if (verbose) {
        cat(sprintf("Processing sample %s, channel %s...", sampleName, colNames[i])); utils::flush.console()
      }

      if (colNames[i] %nin% channels) {
        rr <- rep(NA_integer_, length(b)); names(rr) <- colNames[i]
      } else {
        #plot(stats::density(b))

        ## If too many zeros, might be a blank channel.
        if (sum(b == 0) / length(b) < zero_threshold) {
          multisectArgs <- list(
            x = b,
            bins = bins,
            plot_cutoff = TRUE
          )
          multisectArgs <- utils::modifyList(multisectArgs, multisect...)
          cutoff <- do.call(multisect, multisectArgs)

          r <- Hmisc::cut2(b, cutoff)
        } else {
          cutoff <- rep(NA_real_, length(intensity_levels))

          r <- rep(NA_integer_, length(b))
        }
        length(intensity_levels) <- length(cutoff) + 1
        levels(r) <- intensity_levels
        cutoffs <<- c(cutoffs, list(cutoff))

        rr <- r; names(rr) <- colNames[i]
      }

      if (verbose) {
        cat(". Done.", fill = TRUE); utils::flush.console()
      }

      i <<- i + 1
      rr
    })

  colNames <- colnames(e)
  rm(e)
  plyr::l_ply(seq_along(pmm), function(i) { pmm[[i]] <<- as.data.frame(pmm[[i]]); names(pmm[[i]]) <<- colNames[i] })
  pmm <- purrr::reduce(pmm, dplyr::bind_cols)

  names(cutoffs) <- channels
  attr(pmm, "cutoffs") <- cutoffs

  pmm
}

## usage:
# pmm <- find_plus_minus_by_channel(ff, stringr::regex("time|event_length|sample_id", ignore_case = TRUE))


#' @export
get_channels_by_sample <- function(
  x,
  keep_sans_desc = NULL
)
{
  l0 <- sapply(seq_along(x),
    function(i)
    {
      ff <- flowCore::read.FCS(x[i], transformation = FALSE, truncate_max_range = FALSE)
      p <- flowCore::pData(flowCore::parameters(ff))

      p
    }, simplify = FALSE)

  l <- l0
  ## N.B. Use expression 'keep_sans_desc' to make changes to 'l' before continuing.
  plinth::poly_eval(keep_sans_desc)

  i <- 1
  d <- c(
    head(l, 1)[[1]] %>% dplyr::select(name, desc) %>% dplyr::rename(desc01 = desc) %>% list,
    tail(l, -1)
  ) %>%
  purrr::reduce(
    function(a, b)
    {
      i <<- i + 1
      newName <- sprintf("desc_%02d", i)

      dplyr::full_join(
        a,
        dplyr::select(b, name, desc) %>% dplyr::rename(!!newName := desc),
        by = "name")
    })

  channelNamesTables <- apply(dplyr::select(d, -name), 1, table)
  if ((channelNamesTables %>% sapply(length) > 1) %>% any)
    warning("Some channels may have multiple names among samples")

  structure(d, channel_names_tables = channelNamesTables)
}


## Create augmented 'flowCore::flowFrame' object files to use for analysis.
## N.B. Until event sampling is done, we'll be working w/ only 1 file at a time.
#' @export
prepare_fcs_data <- function(
  x, # Vector of file paths
  b = 1/150, # asinh transformation parameter: FCM = 1/150, CyTOF = 1/8 (v. MetaCyto vignette)
  get_channels_by_sample... = list(),
  manage_channels = NULL, # A function or 'rlang::expr({})' object
  excluded_transform_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  data_dir,
  remove_outliers = TRUE, flowCut... = list(),
  outfile_prefix = NULL, outfile_suffix = "_pmm", # &c, also possibly 'NULL'
  ...
)
{
  if (missing(data_dir))
    stop("Must specify a write directory for augmented FCS files")

  if (!missing(data_dir) && !dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)

  if (is.null(outfile_prefix))
    outfile_prefix <- ""
  else
    plinth::poly_eval(outfile_prefix)

  if (is.null(outfile_suffix))
    outfile_suffix <- ""
  else
    plinth::poly_eval(outfile_suffix)

  augmentedFcsFileNames <- sprintf(paste0("%s", basename(x), "%s.RData"), outfile_prefix, outfile_suffix)
  augmentedFcsFilePaths <- paste(data_dir, augmentedFcsFileNames, sep = "/")

  if (any(duplicated(augmentedFcsFileNames)))
    stop("Output filenames must be unique; use an 'outfile_prefix' expression to fix")

  #channels_by_sample <- get_channels_by_sample(x) # Or w/ more control:
  get_channels_by_sampleArgs <- list(
    x = x
  )
  get_channels_by_sampleArgs <- utils::modifyList(get_channels_by_sampleArgs, get_channels_by_sample..., keep.null = TRUE)
  channels_by_sample <- do.call(get_channels_by_sample, get_channels_by_sampleArgs)
  commonChannels <- structure(attr(channels_by_sample, "channel_names_tables") %>% sapply(names), .Names = channels_by_sample$name) %>%
    purrr::compact() %>% unlist

  ## Placeholder if 'poly_eval(manage_channels)' below invokes nothing:
  col_names <- commonChannels
  plinth::poly_eval(manage_channels)

  b <- rep(b, length.out = length(x))

  r <- sapply(seq_along(x),
    function(i)
    {
      ff <- flowCore::read.FCS(x[i], transformation = FALSE, truncate_max_range = FALSE)

      ## TODO: Compensation?
      ## accessing the spillover matrix
      # try(spillover(frame))
      # ff <- flowCore::compensate(ff, ff@description$SPILL)

      colNames <- col_names[!is.na(col_names)]
      if (!is_invalid(missingCols <- setdiff(names(colNames), flowCore::colnames(ff)))) {
        ## Add missing columns to flowFrame
        ff <- flowCore::fr_append_cols(
          ff,
          structure(
            matrix(rep(NA_real_, flowCore::nrow(ff) * length(missingCols)), ncol = length(missingCols)),
            .Dimnames = list(NULL, missingCols)
          )# %>% head
        )
      }

      ## Keep useful non-blank columns
      ff <- ff[, names(colNames)]; flowCore::colnames(ff) <- colNames

      ## Remove and/or rename channels (probably unnecessary now)
      p <- flowCore::pData(flowCore::parameters(ff)) # Use 'p' in expression or function

      if (remove_outliers) {
        flowCutArgs <- list(
          f = ff,
          Directory = paste(data_dir, "flowCut", sep = "/")
        )
        flowCutArgs <- utils::modifyList(flowCutArgs, flowCut..., keep.null = TRUE)

        flowCut_results <- do.call(flowCut::flowCut, flowCutArgs)
        ff <- flowCut_results$frame
        flowCut_results <- flowCut_results[names(flowCut_results) != "frame"]
      }

      temp <- flowCore::colnames(ff); transformChannels <- temp[stringr::str_detect(temp, excluded_transform_channels_re, negate = TRUE)]
      asinhTrans <- flowCore::arcsinhTransform(transformationId = "flowpipe-transformation", a = 1, b = b, c = 0)
      transList <- flowCore::transformList(transformChannels, asinhTrans)
      tff <- flowCore::transform(ff, transList)

      remove(ff)

      pmm <- find_plus_minus_by_channel(tff, ...)

      exprs_tff <- flowCore::exprs(tff)
      attr(exprs_tff, "plus_minus_matrix") <- pmm
      class(exprs_tff) <- c("pmm", class(exprs_tff))
      ## This bypasses the type checking on the 'exprs' slot; trick described here:
      ## https://stat.ethz.ch/R-manual/R-devel/library/methods/html/slot.html
      attr(tff, "exprs") <- exprs_tff

      if (exists("flowCut_results"))
        attr(tff, "flowCut_results") <- flowCut_results

      augmentedFcsFileName <- augmentedFcsFileNames[i]
      augmentedFcsFilePath <- augmentedFcsFilePaths[i]
      cat(sprintf("Saving augmented FCS file as %s...", augmentedFcsFileName)); utils::flush.console()
      save(tff, file = augmentedFcsFilePath)
      cat(". Done.", fill = TRUE)

      augmentedFcsFilePath
    }, simplify = TRUE)

  r
}


## Assure that 'exprs' subsetting includes the plus-minus matrix.
#' @export
`[.pmm` <- function(x, ...)
{
  y <- NextMethod("[") # Dispatch to generic

  attr(y, "plus_minus_matrix") <- attr(x, "plus_minus_matrix") %>% `[`(...)
  ## Correctly subset 'cluster_id' vector attribute if present
  if (!is.null(attr(x, "cluster_id"))) {
    cluster_id <- c()
    ## The following is hacky, but it appears to always work:
    extant_colname <- head(colnames(y), 1)
    if (!is_invalid(extant_colname)) {
      cim <- structure(rep(NA, prod(dim(x), na.rm = TRUE)),
        .Dim = dim(x), .Dimnames = dimnames(x)) %>% as.data.frame
      cim[, extant_colname] <- attr(x, "cluster_id")
      cluster_id <- cim %>% `[`(...) %>% `[`(, extant_colname)
    }
    attr(y, "cluster_id") <- cluster_id
  }
  class(y) <- class(x)

  y
}


#' @export
print.pmm <- function(x, n = formals(utils:::head.default)$n, ...)
{
  head(x, n) %>% print.default(...)
}


#' @export
as.matrix.pmm <- function(x, ...)
{
  unclass(x)
}


#' @export
`colnames<-` <- function(x, value)
  UseMethod("colnames<-")


#' @export
`colnames<-.default` <- function(x, value)
{
  base::`colnames<-`(x, value)
}


#' @export
`colnames<-.pmm` <- function(x, value)
{
  x <- `colnames<-.default`(x, value)
  colnames(attr(x,"plus_minus_matrix")) <- value

  x
}


#' @export
gate <- function(
  x, # Expression matrix of class "pmm"
  strategy = list(),
  visualize_channels... = list(),
  verbose = TRUE
)
{
  if (length(strategy) == 0) {
    warning("No gating strategy supplied; returning original expression matrix.")

    return (x)
  }

  plyr::l_ply(seq_along(strategy),
    function(a)
    {
      pmm <- attr(x, "plus_minus_matrix") %>% tibble::as_tibble()

      ## Plot gating strategy
      visualize_channelsArgs <- list(
        x = x,
        channels = strategy[a] # Named list w/ a single element
      )
      visualize_channelsArgs <-
        utils::modifyList(visualize_channelsArgs, visualize_channels..., keep.null = TRUE)

      do.call(visualize_channels, visualize_channelsArgs)

      x <<- x[with(pmm, plinth::poly_eval(strategy[[a]])), , drop = FALSE]
    })

  if (verbose && !is_invalid(names(strategy))) {
    cat("Gating strategies applied:", paste(names(strategy), collapse = ", "), fill = TRUE)
    utils::flush.console()
  }

  x
}


#' @export
get_expression_subset <- function(
  x, # Vector of file paths to augmented 'flowCore::flowFrame' objects
  gate... = list(), # Preprocess samples individually
  save_plot_fun = grDevices::pdf, save_plot... = list(),
  seed = 666,
  sample_size = 5000,
  callback = NULL # An expression
)
{
  if (!is_invalid(save_plot...)) { # Save gating sequences to PDF
    def_par <- par(no.readonly = TRUE)

    save_plotArgs <- list(
      width = 4.0 * length(gate...$strategy) + 1,
      height = 4.0 * length(x) + 1
    )
    save_plotArgs <- utils::modifyList(save_plotArgs, save_plot..., keep.null = TRUE)

    image_dir <- dirname(save_plotArgs$file)
    if (!dir.exists(image_dir))
      dir.create(image_dir, recursive = TRUE)

    do.call(save_plot_fun, save_plotArgs)

    ## TODO: This needs some serious work; I need to determine a max number of gates per page.
    # graphics::layout(
    #   matrix(seq(length(x) * length(gate...$strategy)), byrow = TRUE, ncol = length(gate...$strategy))
    # )
    par(mfrow = c(length(x), length(gate...$strategy)))
  }

  i <- 0
  ss <- sapply(x,
    function(a)
    {
      ## N.B. I may want to load into an environment if object naming isn't consistent.
      ## (Then use the only object in the environment.)
      load(a)

      exprs_tff <- flowCore::exprs(tff)
      sampleName <- tools::file_path_sans_ext(basename(flowCore::description(tff)$FILENAME))

      remove(tff)

      if (length(gate...) != 0L) {
        gateArgs <- list(
          x = exprs_tff,
          visualize_channels... =
            list(
              # N.B. Change the following line to use 'rlang::enquo()' leading to a call to 'mtext()':
              plot... = list(sub = sampleName),
              plot_end_callback = expression({
                graphics::title(main = sprintf("Gating strategy: %s", names(channels)))
              })
            )
        )
        gateArgs <- utils::modifyList(gateArgs, gate..., keep.null = TRUE)

        exprs_tff <- do.call(gate, gateArgs)
      }

      if (!is.null(seed))
        set.seed(seed)

      ## Create a new "pmm" object with a sample ID column
      e <- exprs_tff[sample(NROW(exprs_tff), pmin(sample_size, NROW(exprs_tff))), ]
      pmm <- attr(e, "plus_minus_matrix"); rownames(pmm) <- NULL

      i <<- i + 1
      id <- rep(i, NROW(e))
      r <- cbind(id = id, e)
      attr(r, "plus_minus_matrix") <- cbind(id = id, pmm)
      class(r) <- class(e)

      structure(r, total_gated_events = NROW(exprs_tff))
    }, simplify = FALSE)

  if (!is_invalid(save_plot...)) {
    dev.off()
    par(def_par)
  }

  if (!is.null(callback))
    plinth::poly_eval(callback) # E.g. standardize column names

  ## Assemble a new "pmm" object from all subsets.
  e <- purrr::reduce(ss, rbind)
  pmm <- purrr::reduce(sapply(ss, function(a) attr(a, "plus_minus_matrix"), simplify = FALSE), rbind); rownames(pmm) <- NULL
  attr(e, "plus_minus_matrix") <- pmm
  attr(e, "id_map") <- structure(seq(i), .Names = x)
  attr(e, "total_gated_events") <- sapply(ss, function(a) attr(a, "total_gated_events"))
  class(e) <- class(ss[[1]])

  e
}


#' @export
make_tsne_embedding <- function(
  x, # "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names otherwise
  seed = 666,
  ...
)
{
  set.seed(seed)
  tsne <- Rtsne::Rtsne(x[, channels, drop = FALSE], ...)

  tsne
}

## usage:
# cordon(make_tsne_embedding, perplexity = 30, envir = globalenv(), file_path = paste(data_dir, "flowpipe-test-tsne.RData", sep = "/"), variables = c("tsne"), timestamp... = list(use_seconds = TRUE), action = "save")


#' @export
make_umap_embedding <- function(
  x, # "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names otherwise
  n_neighbors = 15, min_dist = 0.2, metric = "euclidean",
  seed = 666,
  ...
)
{
  set.seed(seed)
  umap <- uwot::umap(
    x[, channels],
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    ...)

  umap
}

## usage:
# cordon(make_umap_embedding, n_neighbors = 15, envir = globalenv(), file_path = paste(data_dir, "flowpipe-test-umap.RData", sep = "/"), variables = c("umap"), timestamp... = list(use_seconds = TRUE), action = "save")


## Heavily borrowed from 'cytofkit:::FlowSOM_integrate2cytofkit()'
#' @export
simple_FlowSOM <- function (xdata, k, flow_seed = NULL, ...)
{
  xdata <- as.matrix(xdata)

  cat("Building SOM..."); utils::flush.console()

  ord <- tryCatch({
    map <- FlowSOM::SOM(xdata, silent = TRUE, ...)
    cat(". Done.", fill = TRUE)
    cat("Metaclustering to", k, "clusters..."); utils::flush.console()
    metaClusters <- suppressMessages(FlowSOM::metaClustering_consensus(map$codes,
      k = k, seed = flow_seed))
    cat(". Done.", fill = TRUE)
    cluster <- metaClusters[map$mapping[, 1]]
  }, error = function(e) { message("\nError: ", e$message); flush.console(); return (NULL) })

  if (is.null(ord)) {
    cluster <- NULL
  } else {
    if (length(ord) != NROW(xdata)) {
      message("\nError: FlowSOM failed.")

      return (NULL)
    }
    cluster <- ord
  }

  return (cluster)
}


#' @export
make_initial_clusters <- function(
  x, # Matrix, possibly the 'Y' object resulting from 'Rtsne::Rtsne()'
  channels = TRUE, # Vector of column names otherwise
  seed = 666,
  FlowSOM_k = NULL,
  VorteX_path = "./VorteX.jar",
  num_nearest_neighbors = 10,
  Xshift_command = "java -Xmx64G -cp \"%s\" standalone.Xshift -NUM_NEAREST_NEIGHBORS=%d",
  importConfig... = list(),
  tol = 1e-5
)
{
  set.seed(seed)

  x <- x[, channels, drop = FALSE]

  if (is.null(FlowSOM_k)) {
    ### Do X-Shift clustering

    ## X-Shift default arguments
    importConfigArgs <- list(
      CLUSTERING_COLUMNS = paste(seq(NCOL(x)), collapse = ","),
      LIMIT_EVENTS_PER_FILE = -1,
      TRANSFORMATION = "NONE",
      SCALING_FACTOR = 1,
      NOISE_THRESHOLD = 1.0,
      EUCLIDIAN_LENGTH_THRESHOLD = 0.0,
      RESCALE = "NONE",
      QUANTILE = 0.95,
      RESCALE_SEPARATELY = "false"
    )
    importConfigArgs <- utils::modifyList(importConfigArgs, importConfig..., keep.null = TRUE)

    ## N.B. Change this to a package file in directory "extdata" when the time comes:
    ic <- readLines("D:/Users/priscian/my_documents/urmc/2018/studies/flow/flowpipe/importConfig.txt")
    plyr::l_ply(names(importConfigArgs),
      function(a)
      {
        ic <<- stringr::str_replace_all(ic, paste0("%", a, "%"), importConfigArgs[[a]] %>% as.character)
      })

    ## Set up a temporary workspace
    d <- tempdir()
    p <- tempfile(tmpdir = d, fileext = ".fcs")

    ## Create temporary files
    p1 <- flowCore::write.FCS(flowCore::flowFrame(as.matrix(x)), p)
    writeLines(ic, paste(d, "importConfig.txt", sep = "/"))
    writeLines(normalizePath(p1, mustWork = FALSE), paste(d, "fcsFileList.txt", sep = "/"))

    XshiftCommand <- sprintf(Xshift_command, normalizePath(VorteX_path, mustWork = FALSE), num_nearest_neighbors)
    currentWorkingDir <- getwd()
    setwd(d)
    XshiftOutput <- system(XshiftCommand, intern = TRUE)
    setwd(currentWorkingDir)

    xx <- flowCore::read.FCS(paste(d, "out", basename(p1), sep = "/"))

    ## Are the original & X-Shift expression matrices the same except for some tolerance?
    if (!(dplyr::near(flowCore::exprs(xx)[, seq(NCOL(xx) - 1)], x, tol = tol) %>% all))
      warning("Input & X-Shift expression matrices don't match within tolerance")

    cluster_id <- flowCore::exprs(xx)[, "cluster_id"]

    ## This can be plotted;
    ## Examples here: rstudio-pubs-static.s3.amazonaws.com/362044_903076131972463e8fdfcc00885fc9a6.html
    cluster_graph <- igraph::read.graph(paste(d, "out", "mst.gml", sep = "/"), format = c("gml"))
  } else {
    ## Clusters using FlowSOM.
    cluster_id <- simple_FlowSOM(xdata = x[, , drop = FALSE], k = FlowSOM_k, flow_seed = seed)
  }

  cluster_id
}

## Also see:
# clusters_pg <- cytofkit::cytof_cluster(xdata = e[, channels, drop = FALSE], method = "Rphenograph")
## N.B. Might need to remove duplicate rows beforehand!


#' @export
summary.pmm <- function(
  x, # "pmm" object from 'get_expression_subset()'
  n = NULL, # Cluster numbers, NULL for all
  excluded_channels_re = stringr::regex("id|time|event_length", ignore_case = TRUE),
  overall_label_threshold = 0.95,
  label_threshold = 0.90,
  collapse = "",
  element_names = TRUE,
  as_list = FALSE
)
{
  if (is.null(attr(x, "cluster_id"))) {
    attr(x, "cluster_id") <- rep(1, NROW(x))
    warning("PMM object 'x' has no 'cluster_id' attribute")
  }

  channels <- colnames(x)[stringr::str_detect(colnames(x), excluded_channels_re, negate = TRUE)]
  if (is.null(n))
    n <- attr(x, "cluster_id") %>% unique %>% sort

  pmm <- attr(e, "plus_minus_matrix")[, channels]
  if (!is.null(overall_label_threshold)) {
    overall_channels <- (plyr::aaply(pmm, 2, table)/NROW(pmm) > overall_label_threshold) %>% apply(1, any) %>% `!`
    channels <- names(overall_channels)[overall_channels]
    if (any(!overall_channels))
      warning(sprintf("The following channels are overrepresented in all cells: %s",
        paste(names(overall_channels)[!overall_channels], collapse = " ")))
  }

  r <- sapply(n,
    function(i)
    {
      e <- x[attr(x, "cluster_id") %in% i, channels, drop = FALSE]
      pmm <- attr(e, "plus_minus_matrix")

      l <- sapply(colnames(pmm),
        function(a)
        {
          tab <- table(pmm[, a])/NROW(pmm)

          r <- ""
          if (any(tab >= label_threshold)) {
            r <- a %_% names(tab)[tab >= label_threshold]
          }
        }, simplify = FALSE) %>% unlist

      if (as_list)
        l
      else
        l %>% paste(collapse = collapse)
    }, simplify = ifelse(as_list, FALSE, TRUE))

  if (is.logical(element_names)) {
    if (element_names)
      names(r) <- as.character(n)
  } else if (is.character(element_names)) {
    names(r) <- element_names
  }

  r
}

## usage:
# summary(e, excluded_channels_re = excluded_channels_re, label_threshold = 0.90, as_list = FALSE)


#' @export
search <- function(x, ...)
  UseMethod("search")


#' @export
search.default <- function(x, ...)
{
  search(x, ...)
}


#' @export
search.pmm <- function(
  x, # "pmm" object from 'get_expression_subset()'
  query, # Vector of search terms based on channel names
  query_re = "^(%s)", # RegEx template for search
  summary... = list(), # Additional arguments to 'summary.pmm()'
  ids_only = TRUE
)
{
  summaryArgs <- list(
    x = x,
    as_list = TRUE
  )
  summaryArgs <- utils::modifyList(summaryArgs, summary..., keep.null = TRUE)

  sm <- do.call(summary, summaryArgs)

  r <- sapply(sm,
    function(a)
    {
      re <- stringr::regex(sprintf(query_re, paste(rex::escape(query), collapse = "|")), ignore_case = TRUE)
      d <- a[stringr::str_detect(a, re)]

      ## Were all the query terms found?
      if (length(table(d)) == length(query))
        return (d)

      NULL
    }, simplify = FALSE) %>% purrr::compact()

  if (is_invalid(r))
    return (NULL)

  if (ids_only)
    return (names(r))

  r
}

## usage:
# r <- search(e, c("cd45+", "cd3-"), summary... = list(excluded_channels_re = excluded_channels_re, label_threshold = 0.90))


#' @export
merge_clusters <- function(
  x, # "pmm" object from 'get_expression_subset()'
  clusters, # Named list of mutually exclusive cell subsets
  excluded_channels_re,
  label_threshold,
  search... = list(),
  verbose = TRUE,
  leftover_clusters = NULL
)
{
  searchArgs <- list(
    x = x
  )
  if (!missing(excluded_channels_re))
    searchArgs$summary...$excluded_channels_re <- excluded_channels_re
  searchArgs <- utils::modifyList(searchArgs, search..., keep.null = TRUE)

  label_thresholds <- structure(
    rep(formals(summary.pmm)$label_threshold, length(clusters)),
    .Names = names(clusters)
  )
  if (!missing(label_threshold)) {
    if (is_invalid(names(label_threshold)))
      names(label_threshold) <- rep("", length(label_threshold))

    namedThresholds <- label_threshold[names(label_threshold) != ""]
    if (!is_invalid(namedThresholds))
      label_thresholds <-
        replace(label_thresholds, names(namedThresholds), namedThresholds)

    unnamedThresholds <- label_threshold[names(label_threshold) == ""]
    if (!is_invalid(unnamedThresholds)) {
      indices <- names(label_thresholds) %nin% names(namedThresholds)
      label_thresholds[indices] <-
        rep(unnamedThresholds, length.out = length(label_thresholds[indices]))
    }
  }

  cc <- sapply(names(clusters),
    function(a)
    {
      searchArgs$query <- clusters[[a]]
      searchArgs$summary...$label_threshold <- label_thresholds[a]

      if (verbose) {
        cat(sprintf("Querying for '%s' clusters at %0.2f threshold...", a,
          searchArgs$summary...$label_threshold))
        utils::flush.console()
      }

      r <- do.call(search, searchArgs)

      if (verbose) {
        cat(". Done.", fill = TRUE); utils::flush.console()
      }

      r
    }, simplify = FALSE)

  cc0 <- cc[sapply(cc, is.null)]
  if (length(cc0) > 0)
    warning(sprintf("Clusters %s were not found", cc0 %>% names %>% sQuote %>% paste(collapse = ", ")))

  cc1 <- cc %>% purrr::compact()

  merged_clusters <- list(
    new_cluster = attr(x, "cluster_id"),
    orig_cluster = attr(x, "cluster_id")
  )

  if (is_invalid(cc1)) # No new clusters found
    return (merged_clusters)

  ## Create 'replace()' arguments
  replaceArgss <- sapply(names(cc1),
    function(a)
    {
      list(
        list = attr(x, "cluster_id") %in% cc1[[a]] %>% which,
        value = a
      )
    }, simplify = FALSE)

  new_cluster <- merged_clusters$new_cluster
  plyr::l_ply(replaceArgss,
    function(a) new_cluster <<- replace(new_cluster, a$list, a$value))

  if (!is.null(leftover_clusters)) local({
    i <- sapply(replaceArgss, function(a) a$list, simplify = FALSE) %>% unlist %>% as.vector
    leftovers_list <- rep(TRUE, length(new_cluster))
    leftovers_list[i] <- FALSE

    new_cluster <<- replace(new_cluster, leftovers_list, leftover_clusters)
  })

  merged_clusters$new_cluster <- new_cluster

  merged_clusters
}

### Differential expression

#' @export
do_differential_expression <- function(
  x, # "pmm" object from 'get_expression_subset()'
  m, # metadata data frame
  model_formula,
  id_map_re = ".*",
  metadata_id_var = "id"
)
{
  id_map <- structure(attr(x, "id_map"),
    .Names = names(attr(x, "id_map")) %>% basename %>% stringr::str_extract(id_map_re))
  ids <- names(id_map)[x[, "id"]]
  clusters <- attr(x, "cluster_id")
  cell_counts <- table(clusters, ids)
  ## Convert to percentages if samples are of uneven sizes
  cell_counts_percent <- (cell_counts /
    matrix(rep(table(ids), NROW(cell_counts)), nrow = NROW(cell_counts), byrow = TRUE)) * 100

  dge <- edgeR::DGEList(cell_counts_percent, lib.size = table(ids))
  design <- stats::model.matrix(model_formula,
    data = metadata %>% dplyr::full_join(structure(dataframe(rownames(dge$samples)), .Names = metadata_id_var), .))
  y <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  fit
}


#' @export
test_de_contrasts <- function(
  fit,
  contrasts,
  include_other_vars = FALSE
)
{
  if (is.logical(include_other_vars) && include_other_vars) {
    ov <-
      setdiff(colnames(fit$design)[stringr::str_detect(colnames(fit$design), "(Intercept)", negate = TRUE)], names(contrasts))
    contrasts <- utils::modifyList(
      contrasts,
      structure(as.list(ov), .Names = ov)
    )
  }

  res <- contrasts %>% plyr::llply(
    function(a, b)
    {
      glmQLFTestArgs <- utils::modifyList(
        list(glmfit = fit),
        structure(list(a), .Names = ifelse(is.matrix(a), "contrast", "coef")),
        keep.null = TRUE)
      do.call(edgeR::glmQLFTest, glmQLFTestArgs)
    })

  res_sig <- plyr::llply(res, function(a) a %>% edgeR::topTags(Inf) %>% as.data.frame %>% dplyr::filter(FDR < 0.05))

  res_sig
}
