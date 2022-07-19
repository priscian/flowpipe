## https://shenbaba.weebly.com/blog/how-to-use-the-pac-measure-in-consensus-clustering
## https://dx.doi.org/10.13140/RG.2.1.5075.8480
#' @export
get_optimal_cluster_count <- function(
  x, # Matrix, usu. "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names if not TRUE
  max_k = 28,
  seed = 666,
  ConsensusClusterPlus... = list()
)
{
  x <- x[, channels, drop = FALSE]

  ConsensusClusterPlusArgs <- list(
    d = as.matrix(x),
    maxK = max_k, # Max limit for 'clusterAlg = "hc"' is 28
    reps = 100,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    title = "consensus-clusters",
    innerLinkage = "complete",
    seed = seed,
    plot = NULL
  )
  ConsensusClusterPlusArgs <-
    utils::modifyList(ConsensusClusterPlusArgs, ConsensusClusterPlus..., keep.null = TRUE)

  ccp_res <- tryCatch({
    do.call(ConsensusClusterPlus::ConsensusClusterPlus, ConsensusClusterPlusArgs)
  }, error = function(e) { message("\nError: ", e$message); flush.console(); return (NULL) })

  if (is.null(ccp_res))
    return (NA_real_)

  ## PAC (proportion of ambiguous clustering) implementation
  Kvec <- 2:max_k
  x1 <- 0.1; x2 <- 0.9 # Threshold defining intermediate sub-interval
  PAC <- rep(NA, length(Kvec))
  names(PAC) <- paste0("K = ", Kvec) # 2:max_k
  for (i in Kvec) {
    M <- ccp_res[[i]]$consensusMatrix
    Fn <- stats::ecdf(M[lower.tri(M)])
    PAC[i - 1] <- Fn(x2) - Fn(x1)
  }

  ## Optimal K value
  optK <- Kvec[which.min(PAC)]

  optK
}


## Heavily borrowed from 'cytofkit:::FlowSOM_integrate2cytofkit()'
#' @export
simple_FlowSOM <- function(
  xdata,
  k,
  flow_seed = NULL,
  ...
)
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
make_clusters <- function(
  x, # Matrix, usu. "pmm" object from 'get_expression_subset()'
  method = c(
    "Rphenograph",
    "FlowSOM",
    "Xshift",
    "ClusterX",
    "DensVM"
  ),
  channels = TRUE, # Vector of column names if not TRUE
  seed = 666,
  ## cytofkit
  cytof_cluster... = list(), Rphenograph_k = 50,
  ## FlowSOM
  FlowSOM_k = 40, estimate_cluster_count = TRUE,
  ## X-shift
  VorteX_path = "./VorteX.jar",
  num_nearest_neighbors = 40,
  Xshift_command = "java -Xmx64G -cp \"%s\" standalone.Xshift -NUM_NEAREST_NEIGHBORS=%d",
  importConfig... = list(),
  tol = 1e-5
)
{
  ## This is necessary for using 'keystone::cordon()', because 1-arg version of 'match.arg()' fails -- TODO
  method <- match.arg(method, formals(make_clusters)$method %>% eval)

  x <- x[, channels, drop = FALSE]

  cluster_id <- switch(method,
    `FlowSOM` = (function() {
      ### Do FlowSOM clustering

      if (estimate_cluster_count) {
        opt_k <- get_optimal_cluster_count(x)
        ## 'get_optimal_cluster_count()' tends to fail for small data sets, so use minimum 'FlowSOM_k'
        if (is.na(opt_k) || opt_k < FlowSOM_k)
          FlowSOM_k <- max(opt_k, 3, na.rm = TRUE)
      }

      simple_FlowSOM(xdata = x, k = FlowSOM_k, flow_seed = seed)
    })(),

    `Xshift` = (function() {
      set.seed(seed)

      ### Do X-shift clustering

      ## For details v. file "X-shift_standalone_README.txt"
      ## X-shift default arguments
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
      ic <- readLines(system.file("inst/templates/importConfig.txt", package = "flowpipe"))
      plyr::l_ply(names(importConfigArgs),
        function(a)
        {
          ic <<- stringr::str_replace_all(ic, paste0("%", a, "%"), importConfigArgs[[a]] %>% as.character)
        })

      ## Set up a temporary workspace
      d <- tempdir(check = TRUE) %>% normalizePath(winslash = "/", mustWork = FALSE)
      p <- tempfile(tmpdir = d, fileext = ".fcs") %>% normalizePath(winslash = "/", mustWork = FALSE)
      ## N.B. "Output is written into an automatically created subdir within the current directory named 'out'."
      o <- paste(d, "out", sep = "/")

      ## Create temporary files
      p1 <- flowCore::write.FCS(flowCore::flowFrame(as.matrix(x)), p)
      writeLines(ic, paste(d, "importConfig.txt", sep = "/"))
      writeLines(normalizePath(p1, mustWork = FALSE), paste(d, "fcsFileList.txt", sep = "/"))

      XshiftCommand <- sprintf(Xshift_command, normalizePath(VorteX_path, mustWork = FALSE), num_nearest_neighbors)
      currentWorkingDir <- getwd()
      setwd(d)
      #if (!dir.exists(o)) dir.create(o, recursive = TRUE) # Make sure output directory exists
      XshiftOutput <- system(XshiftCommand, intern = TRUE)
      setwd(currentWorkingDir)

      ## On clustering failure, return single cluster
      #if (!file.exists(paste(o, basename(p1), sep = "/"))) browser()
      if (is.null(attr(XshiftOutput, "status"))) {
        xx <- flowCore::read.FCS(paste(o, basename(p1), sep = "/"))

        ## Are the original & X-shift expression matrices the same except for some tolerance?
        if (!(dplyr::near(flowCore::exprs(xx)[, seq(NCOL(xx) - 1)], x, tol = tol) %>% all))
          warning("Input & X-Shift expression matrices don't match within tolerance")

        cluster_id <- flowCore::exprs(xx)[, "cluster_id"]

        ## This can be plotted;
        ## Examples here: rstudio-pubs-static.s3.amazonaws.com/362044_903076131972463e8fdfcc00885fc9a6.html
        cluster_graph <- igraph::read.graph(paste(o, "mst.gml", sep = "/"), format = c("gml"))
      } else { # Clustering failed
        cluster_id <- rep(0, NROW(x))
      }

      return (structure(cluster_id, cluster_graph = cluster_graph, XshiftOutput = XshiftOutput))
    })(),

    (function() {
      ## N.B. This list is for 'cytofkit2::cytof_cluster()', which provides additional clustering methods:
      cytof_clusterArgs <- list(
        xdata = x,
        method  = method,
        Rphenograph_k = Rphenograph_k
      )
      ## N.B. This list is for 'Rphenograph::Rphenograph()', with only the one method:
      cytof_clusterArgs <- list(
        data = x,
        k = Rphenograph_k
      )
      cytof_clusterArgs <- utils::modifyList(cytof_clusterArgs, cytof_cluster..., keep.null = TRUE)

      tryCatch({
        # do.call(cytofkit2::cytof_cluster, cytof_clusterArgs)
        do.call(Rphenograph::Rphenograph, cytof_clusterArgs)
      }, error = function(e) { message("\nError: ", e$message); flush.console(); return (NULL) })
    })()
  )

  cluster_id
}

## Also see:
# clusters_pg <- cytofkit2::cytof_cluster(xdata = e[, channels, drop = FALSE], method = "Rphenograph")
## N.B. Might need to remove duplicate rows beforehand!


#' @export
make_metaclusters <- function(
  x, # "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names if not TRUE
  make_clusters... = list(),
  # make_clusters... = list( # Some useful defaults
  #   Rphenograph_k = 50,
  #   FlowSOM_k = 40, estimate_cluster_count = FALSE,
  #   num_nearest_neighbors = 30
  # ),
  make_metaclusters... = list(), # Passed to 'make_clusters()' for metaclustering step
  centroid_fun = median
)
{
  make_clustersArgs <- list(
    channels = channels
  )

  l <- x %>% as.data.frame %>%
    dplyr::group_by(id) %>%
    dplyr::group_map(
      .f = ~ (function(x, y)
      {
        ## Also see 'rlist::list.append()'
        mica <- utils::modifyList(c(make_clustersArgs, list(x = x %>% data.matrix)), make_clusters...,
          keep.null = TRUE)

        cluster_ids <- do.call(make_clusters, mica)
        ## If 'make_clusters()' fails, assign all events to single cluster
        if (is.null(cluster_ids))
          cluster_ids <- rep(1, NROW(x))
        cluster_ids %>% table(dnn = "sample " %_% y) %>% print; utils::flush.console()

        structure(cluster_ids, sample_id = dplyr::pull(y, id))
      })(.x, .y), .keep = TRUE)

  centroids <- x %>% as.data.frame %>%
    dplyr::select(all_of(c("id", channels %>% as.vector))) %>%
    dplyr::rename(sample_id = id) %>%
    dplyr::mutate(cluster_id = unlist(l)) %>%
    dplyr::relocate(cluster_id, .after = sample_id) %>%
    dplyr::group_by(sample_id, cluster_id) %>%
    dplyr::group_modify(
      .f = ~ (function(x, y)
      {
        d <- x %>% dplyr::select(-c("sample_id", "cluster_id"))

        ## Find centroid for each group
        dplyr::summarize(d, across(everything(), ~ centroid_fun(.x, na.rm = TRUE)))
      })(.x, .y), .keep = TRUE)

  make_metaclustersArgs <- make_clustersArgs

  mica <- utils::modifyList(c(make_metaclustersArgs,
    list(x = centroids %>% data.matrix)), make_metaclusters..., keep.null = TRUE)
  centroid_cluster_id <- do.call(make_clusters, mica)
  centroid_cluster_id %>% table(dnn = "metaclusters") %>% print

  ## Match metaclusters back to individual events
  centroids_clustered <- centroids %>%
    dplyr::ungroup() %>%
    dplyr::mutate(centroid_cluster_id = centroid_cluster_id) %>%
    dplyr::relocate(centroid_cluster_id, .after = cluster_id)

  centroid_cluster_df <- sapply(l,
    function(a) keystone::dataframe(sample_id = attr(a, "sample_id"), cluster_id = a),
      simplify = FALSE) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::left_join(centroids_clustered %>% dplyr::select(sample_id, cluster_id, centroid_cluster_id),
      by = c("sample_id", "cluster_id"))

  event_metacluster_id <- centroid_cluster_df$centroid_cluster_id
  if (is.na(event_metacluster_id) %>% any)
    warning("Some events are incorrectly unmatched to centroid clusters")
  event_metacluster_id %>%
    table(useNA = "always", dnn = "event metaclusters") %>% print

  event_metacluster_id
}


#' @export
summary.pmm <- function(
  x, # "pmm" object from 'get_expression_subset()'
  n = NULL, # Cluster numbers: NULL or TRUE for all, FALSE for every event
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  channels = colnames(x),
  merged_labels = list(
    `-/d` = c("-", "d"),
    `+/++` = c("+", "++")
    #`d/+` = c("d", "+")
  ),
  overall_merged_labels_index = 1:2,
  overall_label_threshold = Inf,
  label_threshold = 0.90,
  collapse = "", expression_level_sep = ",",
  element_names = TRUE,
  as_list = FALSE
)
{
  clusterId <- attr(x, "cluster_id")

  byEvent <- FALSE
  if (is.logical(n)) {
    if (!n) {
      byEvent <- TRUE
    }

    n <- NULL
  }

  if (is.null(attr(x, "cluster_id"))) {
    clusterId <- sprintf("%d", seq(NROW(x)))
    stop("PMM object 'x' has no 'cluster_id' attribute")
  }

  if (is.matrix(clusterId))
    clusterId <- clusterId[, which_cluster_set] %>% drop

  if (is.null(n)) {
    n <- clusterId %>% unique # But don't sort, because character-numbers don't stay in numeric order
    if (!byEvent)
      n <- n %>% sort
  }

  pmm <- attr(x, "plus_minus_matrix")[, channels]
  if (!is.null(overall_label_threshold)) {
    # comp <- (plyr::aaply(pmm, 2, table)/NROW(pmm)) %>%
    #   as.data.frame %>% tibble::rownames_to_column()
    comp <- sapply(pmm, table, simplify = FALSE) %>% purrr::compact() %>%
      { structure(dplyr::bind_rows(.) %>% as.data.frame, row.names = names(.)) } %>%
      data.matrix %>% `/`(NROW(pmm)) %>%
      as.data.frame %>% tibble::rownames_to_column()
    ## N.B. I don't want to include *all* merged labels here!
    plyr::l_ply(names(merged_labels)[overall_merged_labels_index],
      function(a)
      {
        comp <<- comp %@>% dplyr::rowwise() %@>% dplyr::mutate(
          !!a := sum(!!!rlang::syms(merged_labels[[a]]))
        )
      })
    comp <- comp %>% tibble::column_to_rownames() %>% data.matrix
    overall_channels <- (comp > overall_label_threshold) %>% apply(1, any) %>% `!`
    channels <- names(overall_channels)[overall_channels]
    if (any(!overall_channels))
      warning(sprintf("The following channels are overrepresented in all cells: %s",
        paste(names(overall_channels)[!overall_channels], collapse = " ")))
  }

  ## Here, 'comp' should look something like this:
  #                   -          d          +          ++       -/d       +/++
  # IL-23_p19 0.8236050 0.10851854 0.05529584 0.012580632 0.9321235 0.06787647
  # CD69      0.8493953 0.07843002 0.04751813 0.024656565 0.9278253 0.07217470
  # TGFb      0.8752095 0.08020844 0.02647643 0.018105611 0.9554180 0.04458204
  # IL-17A    0.8639330 0.07175749 0.04820795 0.016101551 0.9356905 0.06430950
  # IL-10     0.8733086 0.07338119 0.04515121 0.008158991 0.9466898 0.05331020
  # CCR7      0.8402868 0.09671154 0.04927721 0.013724493 0.9369983 0.06300171
  # [...]
  ##
  ## The non-'merged_labels' columns should add to 1, i.e. '(rowSums(comp[, 1:4]) == 1) %>% all' is TRUE.
  ## The 'merged_labels' columns should add to their component non-merged columns, e.g. "-/d" = "-" + "d".
  ## 'comp' summarizes the phenotypic composition of all clusters at once as the proportion of each label count
  ##   relative to all the events.

  if (byEvent) {
    e <- x[, channels, drop = FALSE]
    pmm <- attr(e, "plus_minus_matrix")

    mpmm <- as.matrix(pmm)
    merges <- sapply(names(merged_labels),
      function(a)
      {
        mpmm %in% merged_labels[[a]] %>% `dim<-`(dim(mpmm)) %>%
          `is.na<-`(. == FALSE) %>% `[<-`(., ., a)
      }, simplify = FALSE)

    allLabels <- c(list(as.list(t(mpmm))), sapply(merges, function(a) as.list(t(a))
      %>% `[<-`(is.na(.), list(NULL)), simplify = FALSE))

    r <- purrr::pmap(allLabels,
      function(...) { as.vector(c(...)) }) %>%
        `names<-`(rep(colnames(pmm), length.out = length(.))) %>%
        keystone::chunk(NCOL(pmm))

    if (!as_list) {
      r <- sapply(re,
        function(l) { sapply(names(l), function(a) a %_% paste(l[[a]], collapse = expression_level_sep)) %>%
          paste(collapse = collapse) }, simplify = FALSE)
    }
  } else {
    r <- sapply(n, # This doesn't appear to benefit if 'keystone::psapply()' is dropped in here -- it's worse, in fact!
      function(i)
      {
        e <- x[clusterId %in% i, channels, drop = FALSE]
        pmm <- attr(e, "plus_minus_matrix")

        l <- sapply(colnames(pmm),
          function(a)
          {
            comp <- (table(pmm[, a])/NROW(pmm)) %>%
              data.matrix %>% t %>% as.data.frame %>%
              `rownames<-`(a) %>% tibble::rownames_to_column()
            plyr::l_ply(names(merged_labels),
              function(a)
              {
                comp <<- comp %@>% dplyr::rowwise() %@>% dplyr::mutate(
                  !!a := sum(!!!rlang::syms(merged_labels[[a]]))
                )
              })
            comp <- comp %>% keystone::dataframe() %>% tibble::column_to_rownames() %>%
              data.matrix

            ## 'comp' should look something like this:
            #                   -          d         +         ++       -/d       +/++
            # IL-23_p19 0.9394749 0.03492733 0.0209564 0.00464135 0.9744023 0.02559775
            ##
            ## The names of all columns meeting 'label_threshold' (see below) are returned.

            rr <- ""
            if (any(comp >= label_threshold)) {
              rr <- colnames(comp)[comp >= label_threshold]
            }

            rr
          }, simplify = FALSE)

        if (as_list)
          l
        else
          sapply(names(l), function(a) a %_% paste(l[[a]], collapse = expression_level_sep)) %>%
            paste(collapse = collapse)
      }, simplify = ifelse(as_list, FALSE, TRUE))
  }

  if (is.logical(element_names)) {
    if (element_names)
      names(r) <- as.character(n)
  } else if (is.character(element_names)) {
    names(r) <- element_names
  }

  ## If 'as_list = TRUE', 'r' is a list the length of the unique cluster names in the current cluster set,
  ##   w/ elements named after the clusters; each element is a sub-list named after the channels/columns of 'x',
  ##   whose elements contain all the phenotype names (e.g. "-", "+", "+/++", &c) meeting
  ##   the proportion threshold for that channel & cluster. If no phenotype meets the threshold, "" is returned.
  ## If 'as_list = FALSE', 'r' is a list the length of the unique cluster names in the current cluster set,
  ##   each of whose elements is a single string displaying a full set of channel phenotypes separated
  ##   according to 'collapse' & 'expression_level_sep'.
  r
}

## usage:
# summary(e[, -1], label_threshold = 0.90, as_list = TRUE)


#' @export
search <- function(x, ...)
  UseMethod("search")


#' @export
search.default <- function(x, ...)
{
  search(x, ...)
}


## Search plus-minus matrix of "pmm" expression object for channel phenotypes given in 'query', e.g.
##   r <- search(e[, analysis_channels], query = c("cd45+/++", "cd3-/d"), summary... = list(which_cluster_set = 1, label_threshold = 0.55))
## Return Value: A vector of names of clusters whose 'query' channels all meet/exceed their 'label_threshold's,
##   i.e. each cluster returned is a hit for all the 'query' phenotypes.
#' @export
search.pmm <- function(
  x, # "pmm" object from 'get_expression_subset()'
  query, # Vector of search terms based on channel names
  query_re = "^(%s)$", # RegEx template for search
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
      test <- sapply(names(a), function(b) { if (all(a[[b]] == "")) return (NULL); paste0(b, a[[b]]) }, simplify = FALSE) %>%
        unlist(use.names = FALSE)
      ## This produces a list whose elements have >1-length vectors for each either-or query:
      baseQuery <- stringr::str_split(query, stringr::fixed("||", TRUE))
      re <- sapply(baseQuery, function(b) { stringr::regex(sprintf(query_re, paste(rex::escape(b %>% unlist), collapse = "|")), ignore_case = TRUE) }, simplify = FALSE)
      d <- sapply(re, function(b) stringr::str_subset(test, b), simplify = FALSE)

      ## Were all the query terms found?
      if (length(sapply(d, table, simplify = FALSE) %>% purrr::compact()) == length(baseQuery))
      ## If yes, return all those query terms that were found
        { return (d %>% unlist) }

      NULL
    }, simplify = FALSE) %>% purrr::compact()

  if (is_invalid(r))
    return (NULL)

  if (ids_only)
    return (names(r))

  r
}

## usage:
# r <- search(e[, analysis_channels], c("cd45+/++", "cd3-/d"), summary... = list(overall_label_threshold = Inf, label_threshold = 0.90))


#' @export
merge_clusters <- function(
  x, # "pmm" object from 'get_expression_subset()'
  clusters, # Named list of cell subsets
  channels,
  label_threshold,
  which_cluster_set = 1, # Column no. or name; NULL or FALSE to set off by-event search
  search... = list(),
  verbose = TRUE,
  leftover_clusters = NULL
)
{
  origClusterId <- attr(x, "cluster_id")

  byEvent <- FALSE
  if (is.null(which_cluster_set) || (is.logical(which_cluster_set) && !which_cluster_set)) {
    ## N.B. The "event" clusters must be run though 'sprintf()' to prevent exponentiation > 399999.
    attr(x, "cluster_id") <- sprintf("%d", seq(NROW(x)))
    which_cluster_set <- 1
    byEvent <- TRUE
  } else if (is.logical(which_cluster_set) && which_cluster_set) {
    which_cluster_set <- 1
  }

  searchArgs <- list(
    x = x,
    summary... = list(which_cluster_set = which_cluster_set)
  )
  if (!missing(channels))
    searchArgs$summary...$channels <- channels
  searchArgs <- utils::modifyList(searchArgs, search..., keep.null = TRUE)
  if (byEvent)
    searchArgs$summary... <- utils::modifyList(searchArgs$summary..., list(n = FALSE), keep.null = TRUE)

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

  tictoc::tic("Search clusters")

  #cc <- sapply(names(clusters),
  cc <- keystone::psapply(names(clusters),
    function(a)
    {
      searchArgs$query <- clusters[[a]]
      searchArgs$summary...$label_threshold <- label_thresholds[a]

      if (verbose) {
        if (!byEvent)
          cat(sprintf("Querying for '%s' clusters at %0.2f threshold...", a,
            searchArgs$summary...$label_threshold))
        else
          cat(sprintf("Querying for '%s' clusters at event level...", a))
        utils::flush.console()
      }

      r <- do.call(search, searchArgs)

      if (verbose) {
        cat(". Done.", fill = TRUE); utils::flush.console()
      }

      r
    }, simplify = FALSE)

  tictoc::toc()

  cc0 <- cc[sapply(cc, is.null)]
  if (length(cc0) > 0)
    warning(sprintf("Clusters %s were not found", cc0 %>% names %>% sQuote %>% paste(collapse = ", ")))

  cc1 <- cc %>% purrr::compact()

  clusterId <- attr(x, "cluster_id")
  if (is.matrix(clusterId))
    clusterId <- clusterId[, which_cluster_set]

  merged_clusters <- list(
    new_cluster_id = clusterId,
    orig_cluster_id = origClusterId
  )

  if (is_invalid(cc1)) { # No new clusters found
    if (byEvent)
      merged_clusters <- list(new_cluster_id = origClusterId, orig_cluster_id = origClusterId)

    return (merged_clusters)
  }

  ## Create 'replace()' arguments
  replaceArgss <- sapply(names(cc1),
    function(a)
    {
      list(
        list = clusterId %in% cc1[[a]],
        value = a
      )
    }, simplify = FALSE)

  newClusterId <- sapply(replaceArgss,
    function(a)
    {
      r <- replace(merged_clusters$new_cluster_id, a$list %>% which, a$value) %>%
        replace((!a$list) %>% which, NA_character_)

      r
    }, simplify = TRUE)

  ## Now collapse all the mutually exclusive columns together
  newMergedClusterId <- merge_mutually_exclusive_cols(newClusterId) %>%
    cbind(orig = merged_clusters$orig_cluster_id, .)
  ## N.B. If 'merged_clusters$orig_cluster_id' is already a matrix, the name "orig" is unused.

  merged_clusters$new_cluster_id <- newMergedClusterId

  merged_clusters
}


merge_mutually_exclusive_cols <- function(
  ..., # Combination of matrices/vectors having the same no. rows/length
  collapse = "|"
)
{
  d0 <- cbind(...); d <- rlang::duplicate(d0, shallow = FALSE)
  repeat {
    merge_comb <- utils::combn(seq(NCOL(d)), 2, simplify = FALSE)
    didMerge <- FALSE; startNcol <- NCOL(d)

    for (a in merge_comb) {
      print(a)
      colsAreMutuallyExclusive <-
        apply(d[, a], 1, function(b) (!is.na(b)) %>% sum, simplify = TRUE) %>% `==`(2) %>% any %>% `!`

      if (colsAreMutuallyExclusive) {
        ## Merge them into single column
        temp <- d[, a[1]] %>% `[<-`(!is.na(d[, a[2]]), d[, a[2]][!is.na(d[, a[2]])])

        d[, a[1]] <- temp
        ## Name of new column becomes a combination of both starting columns
        newColname <- paste(colnames(d[, a]), collapse = collapse)
        colnames(d)[a[1]] <- newColname
        d <- d[, -a[2], drop = TRUE]

        ## Don't ever finish right after a merge; check for mergeable columns at least once more.
        didMerge <- TRUE

        break
      }
    }

    ## Finish if no. cols are at minimum or are unchanged, & no merge just happened
    if ((NCOL(d) < 3 || NCOL(d) == startNcol) && !didMerge)
      break
  }

  #browser()
  d
}
#merge_mutually_exclusive_cols(orig = cluster_id, new_cluster)
#merge_mutually_exclusive_cols(new_cluster)
