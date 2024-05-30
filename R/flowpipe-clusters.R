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
    clusterAlg = "hc", distance = "pearson",
    # clusterAlg = "km", distance = "euclidean",
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
  VorteX_path = system.file("java/VorteX.jar", package = "flowpipe"),
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
      ic <- readLines(system.file("templates/importConfig.txt", package = "flowpipe"))
      plyr::l_ply(names(importConfigArgs),
        function(a)
        {
          ic <<- stringr::str_replace_all(ic, paste0("%", a, "%"), importConfigArgs[[a]] %>%
            as.character)
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

      XshiftCommand <- sprintf(Xshift_command, normalizePath(VorteX_path, mustWork = FALSE),
        num_nearest_neighbors)
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
        method = method,
        Rphenograph_k = Rphenograph_k
      )
      ## N.B. This list is for 'Rphenograph::Rphenograph()', with only the one method:
      # cytof_clusterArgs <- list(
      #   data = x,
      #   k = Rphenograph_k
      # )
      cytof_clusterArgs <- utils::modifyList(cytof_clusterArgs, cytof_cluster..., keep.null = TRUE)

      tryCatch({
        do.call(cytofkit2::cytof_cluster, cytof_clusterArgs)
        # do.call(Rphenograph::Rphenograph, cytof_clusterArgs)
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

        if (memoise::is.memoised(make_clusters))
          cluster_ids <- do.call(environment(make_clusters)$`_f`, mica)
        else
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
  if (memoise::is.memoised(make_clusters))
    centroid_cluster_id <- do.call(environment(make_clusters)$`_f`, mica)
  else
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

  #event_metacluster_id
  ## [11 Jan 2023] Make single relevant return value to move away from 'keystone::cordon()'.
  structure(event_metacluster_id, cluster_centroids = centroids_clustered)
}


#' @export
summary.pmm <- function(
  x, # "pmm" object from 'get_expression_subset()'
  n = NULL, # Cluster numbers: NULL or TRUE for all, FALSE for every event
  cluster_set,
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  channels = colnames(x),
  merged_labels = list(
    `-/d` = c("-", "d"),
    `+/++` = c("+", "++"),
    `d/+` = c("d", "+"),
    all = c("-", "d", "+", "++")
  ),
  overall_label_threshold = Inf,
  label_threshold = 0.90,
  collapse = ";", expression_level_sep = ",",
  element_names = TRUE,
  as_list = TRUE
)
{
  if (is_invalid(cluster_set)) {
    clusterId <- attr(x, "cluster_id")
    cluster_set <- NULL
  } else {
    clusterId <- cluster_set
  }

  byEvent <- FALSE
  if (is.logical(n)) {
    if (!n) {
      byEvent <- TRUE
    }

    n <- NULL
  }

  if (is.null(clusterId)) {
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

  pmm <- attr(x, "plus_minus_matrix")[, channels, drop = FALSE]
  if (!is.null(overall_label_threshold)) {
    # comp <- (plyr::aaply(pmm, 2, table)/NROW(pmm)) %>%
    #   as.data.frame %>% tibble::rownames_to_column()
    comp <- sapply(pmm, table, simplify = FALSE) %>% purrr::compact() %>%
      { structure(dplyr::bind_rows(.) %>% as.data.frame, row.names = names(.)) } %>%
      data.matrix %>% `/`(NROW(pmm)) %>%
      as.data.frame %>% tibble::rownames_to_column()
    plyr::l_ply(names(merged_labels),
      function(a)
      {
        comp <<- comp %@>% dplyr::rowwise() %@>% dplyr::mutate(
          !!a := sum(!!!rlang::syms(merged_labels[[a]]))
        )
      })
    comp <- comp %>% `rownames<-`(NULL) %>% tibble::column_to_rownames() %>% data.matrix
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
      r <- sapply(r,
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

            rr <- NULL # Default for channels that meet *none* of the label thresholds
            if (any(comp >= label_threshold)) {
              rr <- colnames(comp)[comp >= label_threshold]
            }

            rr
          }, simplify = FALSE) %>%
            purrr::compact() # Remove any channels that meet *none* of the label thresholds

        if (as_list)
          l
        else
          sapply(names(l), function(a) a %_% paste(l[[a]], collapse = expression_level_sep)) %>%
            paste(collapse = collapse)
      }, simplify = ifelse(as_list, FALSE, TRUE))
  }

  if (!byEvent) {
    if (is.logical(element_names)) {
      if (element_names)
        names(r) <- as.character(n)
    } else if (is.character(element_names)) {
      names(r) <- element_names
    }
  }

  ## If 'as_list = TRUE', 'r' is a list the length of the unique cluster names in the current cluster set,
  ##   w/ elements named after the clusters; each element is a sub-list named after the channels/columns of 'x',
  ##   whose elements contain all the phenotype names (e.g. "-", "+", "+/++", &c) meeting
  ##   the proportion threshold for that channel & cluster. If no phenotype meets the threshold, "" is returned.
  ## If 'as_list = FALSE', 'r' is a list the length of the unique cluster names in the current cluster set,
  ##   each of whose elements is a single string displaying a full set of channel phenotypes separated
  ##   according to 'collapse' & 'expression_level_sep'.
  structure(r, comp = comp) %>% keystone::add_class("summaryPmm")
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
  summary... = list(), # Additional arguments to 'summary.pmm()' or a "summaryPmm" object
  return_type = c("character", "logical", "grid")
)
{
  return_type <- match.arg(return_type)

  if (inherits(summary..., "summaryPmm")) {
    sm <- summary...
  } else {
    summaryArgs <- list(
      x = x,
      as_list = TRUE
    )
    summaryArgs <- utils::modifyList(summaryArgs, summary..., keep.null = TRUE)

    sm <- do.call(summary, summaryArgs)
  }

  comp <- attr(sm, "comp")
  ## Enumerate all possible event states as a named logical vector
  template <- expand.grid(rownames(comp), colnames(comp), stringsAsFactors = FALSE) %>%
    plyr::alply(1, function(a) { unlist(a, use.names = FALSE) %>%
    paste(collapse = "") }) %>% unlist(use.names = FALSE) %>%
    { structure(rep(FALSE, length(.)), .Names = .) }

  ## Multiple OR-conditional gates lead to multiple queries that need testing;
  ##   find all possible combinations and OR them to test for a hit.
  baseQuery <- stringr::str_split(query, "\\s*\\|\\|\\s*") ## Split query elements by "||"
  allQueries <- expand.grid(baseQuery, stringsAsFactors = FALSE) %>%
    plyr::alply(1, unlist, use.names = FALSE)
  mm <- lapply(allQueries,
    function(a)
    {
      sapply(a,
        function(b) { adist(b, names(template), fixed = TRUE) %>% which.min }, simplify = TRUE) %>%
        { names(template)[.] }
    })
  tests <- lapply(mm,
    function(a)
    {
      test <- template
      test[a] <- TRUE

      test
    })

  if (return_type == "grid") {
    r <- sapply(sm,
      function(a) {
        event <- template
        event[unlist(lapply(names(a), function(b) paste0(b, a[[b]])))] <- TRUE

        event
      }, simplify = TRUE)

    return (structure(t(r), gates = mm, query = query))
  }

  r <- lapply(sm,
    function(a)
    {
      event <- template
      event[unlist(lapply(names(a), function(b) paste0(b, a[[b]])))] <- TRUE

      ## Does this event/cluster include the same phenotypes as the query?
      Reduce(`||`, sapply(tests, function(b) sum(b & event) == sum(b)))
    }) %>% unlist(use.names = FALSE)

  if (is_invalid(r))
    return (NULL)

  if (return_type == "character")
    return (which(r) %>% as.character)

  r
}

## usage:
# r <- search(e[, analysis_channels], c("cd45+/++", "cd3-/d"), summary... = list(overall_label_threshold = Inf, label_threshold = 0.90))


#' @export
search_orig.pmm <- function(
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
      test <- sapply(names(a),
        function(b)
        {
          if (all(a[[b]] == "")) return (NULL); paste0(b, a[[b]])
        }, simplify = FALSE) %>% unlist(use.names = FALSE)
      ## This produces a list whose elements have >1-length vectors for each either-or query:
      baseQuery <- stringr::str_split(query, stringr::fixed("||", TRUE))
      re <- sapply(baseQuery,
        function(b)
        {
          stringr::regex(sprintf(query_re, paste(rex::escape(b %>% unlist), collapse = "|")),
            ignore_case = TRUE)
        }, simplify = FALSE)
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
  cluster_set,
  which_cluster_set = 1, # Column no. or name; NULL or FALSE to set off by-event search
  search... = list(),
  verbose = TRUE,
  leftover_clusters = NULL,
  make_gating_poster = FALSE, # Logical, or character path to directory for individual plots
  visualize_channels... = list(),
  devices = flowpipe:::graphics_devices,
  #save_plot_fun = grDevices::pdf, save_plot... = list(compress = FALSE)
  save_plot_fun = grDevices::cairo_pdf, save_plot... = list(onefile = TRUE)
)
{
  if (missing(cluster_set)) {
    clusterId <- attr(x, "cluster_id")
    cluster_set <- NULL
  } else {
    clusterId <- cluster_set
  }
  origClusterId <- clusterId

  byEvent <- FALSE
  if (is.null(which_cluster_set) || (is.logical(which_cluster_set) && !which_cluster_set)) {
    ## N.B. The "event" clusters must be run though 'sprintf()' to prevent exponentiation > 399999.
    #attr(x, "cluster_id") <- sprintf("%d", seq(NROW(x)))
    clusterId <- cluster_set <- sprintf("%d", seq(NROW(x)))
    which_cluster_set <- 1
    byEvent <- TRUE
  } else if (is.logical(which_cluster_set) && which_cluster_set) {
    which_cluster_set <- 1
  }

  searchArgs <- list(
    x = x,
    summary... = list(cluster_set = cluster_set, which_cluster_set = which_cluster_set)
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

  ### Create plots to visually follow a sequence of predefined gates

  gating_poster_dir <- NULL
  if (is.character(make_gating_poster)) {
    gating_poster_dir <- make_gating_poster
    make_gating_poster <- TRUE

    if (!dir.exists(gating_poster_dir))
      dir.create(gating_poster_dir, recursive = TRUE)
  }

  tictoc::tic("Search clusters")

  cc <- NULL
  if (make_gating_poster && byEvent) {
    ## This probably doesn't dispatch on 'summary' alone because of the name/position of the 1st argument
    sm <- do.call(summary.pmm, utils::modifyList(searchArgs$summary...,
      list(x = x), keep.null = TRUE))

    ## Prepare data set to proceed through & plot predefined gating sequences
    ## N.B. For size considerations, I might want to plot inside 'sapply()' & return NULL
    cc_grid <- keystone::psapply(names(clusters),
      function(a)
      {
        searchArgs$query <- clusters[[a]]
        searchArgs$summary...$label_threshold <- label_thresholds[a]
        searchArgs$return_type <- "grid"
        searchArgs$summary... <- sm

        if (verbose) {
          cat(sprintf("Querying for '%s' clusters at event level...", a))
          utils::flush.console()
        }

        r <- do.call(search, searchArgs)

        if (verbose) {
          cat(". Done.", fill = TRUE); utils::flush.console()
        }

        r
      }, simplify = FALSE)

    ## Ordering the colnames by decreasing length will prevent e.g. a match between
    ##   "CD4" & "CD45" before the regex search has gotten to "CD45".
    re <- stringr::regex(stringr::str_flatten(rex::escape(colnames(x)[colnames(x)
      %>% nchar %>% order(decreasing = TRUE)]), "|"))

    `cc+grobs` <- keystone::psapply(seq_along(cc_grid), # So 'a' can be used for numbering plots
      function(a)
      {
        chunks <- sapply(attr(cc_grid[[a]], "gates"), function(b) keystone::chunk(b, 2), simplify = FALSE)
        ## The first "chunk" will have the same no. of elements as all the others:
        tests <- sapply(seq_along(chunks[[1]]),
          function(b) sapply(chunks, function(g) g[[b]], simplify = FALSE) %>%
            unique, simplify = FALSE)

        query <- attr(cc_grid[[a]], "query")
        query_chunks <- keystone::chunk(query, 2)
        cat(sprintf("%s:", names(cc_grid)[a]), fill = TRUE); print(query); utils::flush.console()

        gated_events <- rep(TRUE, NROW(cc_grid[[a]]))
        flit <- sapply(seq_along(tests), # So 'b' can be used for numbering plots
          function(b)
          {
            ## 'NCOL(.)' handles the case where the test matrix has only one column:
            r <- Reduce(`|`, sapply(tests[[b]], function(g) { cc_grid[[a]][, g, drop = FALSE] %>% { `&`(.[, 1], .[, NCOL(.)]) } },
              simplify = FALSE), accumulate = TRUE)

            ## For each list element, create a biaxial plot
            grobs <- mapply(function(k, l)
            {
              plot_channels <- stringr::str_match_all(paste(k, collapse = " "), re)[[1]] %>% drop %>% unique
              ## N.B. Uncomment 'event_mask' just below to plot only events selected by the previous gate:
              visualize_channelsArgs <- list(
                x = x,
                channels = list(gated_events & r[[l]]),
                event_mask = gated_events,
                extract_gating_channels = function(...) plot_channels,
                points... = list(col = scales::alpha("red", 0.5)),
                plot_end_callback = function(...) { # A function will carry its environment along w/ itself
                  graphics::title(main = sprintf("Gate: %s", paste(query_chunks[[b]], collapse = " & ")), cex.main = 0.9, ...)
                  if (l > 1) graphics::mtext("(OR'd with previous gate)", cex = 0.9)
                  graphics::mtext(sprintf("Events: %d/%d", sum(gated_events & r[[l]]), sum(gated_events)),
                    side = 1, line = -1, cex = 0.8)
                }
              )
              visualize_channelsArgs <-
                utils::modifyList(visualize_channelsArgs, visualize_channels..., keep.null = TRUE)

              grobs <- list()
              if (!is.null(gating_poster_dir)) {
                # plyr::l_ply(seq_along(devices),
                #   function(d)
                #   {
                #     ext <- devices[[d]]$ext; devices[[d]]$ext <- NULL
                #     ## Reduce resolution for 'png()' etc. to a manageable value:
                #     if ("res" %in% names(formals(eval(parse(text = names(devices)[d]))))) devices[[d]]$res <- 150
                #     do.call(eval(parse(text = names(devices)[d])),
                #       modifyList(devices[[d]],
                #         list(
                #           width = 5, height = 5,
                #           file = sprintf("%s/%03d-%03d%s_gate-%s",
                #             gating_poster_dir, a, b, letters[l], paste(plot_channels, collapse = "&")) %_% paste0(".", ext)
                #         )
                #       )
                #     )
                #     dev.control(displaylist = "enable")

                #     do.call(visualize_channels, visualize_channelsArgs)

                #     if (d == length(devices)) {
                #       grobs <<- append(grobs, list(grDevices::recordPlot()))
                #     }
                #     dev.off()
                #   })

                gatePlotPath <- tempfile()
                grDevices::png(file = gatePlotPath, bg = "transparent")
                dev.control(displaylist = "enable")

                do.call(visualize_channels, visualize_channelsArgs)

                gatePlot <- grDevices::recordPlot()
                invisible(dev.off())
                unlink(gatePlotPath)

                grobs <- append(grobs, list(gatePlot))
              } else {
                do.call(visualize_channels, visualize_channelsArgs)
              }

              grobs
            }, tests[[b]], seq_along(r), USE.NAMES = TRUE, SIMPLIFY = FALSE)

            gated_events <<- gated_events & r[[length(r)]]
            print(table(gated_events)); utils::flush.console()

            grobs
          }, simplify = FALSE)

        list(gated_events = gated_events, grobs = flit %>% purrr::flatten())
      }, simplify = FALSE)

    grobs <- NULL
    if (!is.null(gating_poster_dir)) {
      grobs <- sapply(`cc+grobs`, function(a) a$grobs, simplify = FALSE) %>% `names<-`(names(cc_grid))
      ## Keep list of grobs for e.g. single plots, different image types:
      saveRDS(object = grobs, file = paste(data_dir, "gated-clusters-poster.rds", sep = "/"))
    }
    cc <- sapply(`cc+grobs`, function(a) a$gated_events, simplify = FALSE) %>%
      sapply(function(a) { as.vector(a) %>% which %>% as.character }, simplify = FALSE) %>% `names<-`(names(cc_grid))
    rm(`cc+grobs`)

    ## Finally, create full gating poster
    if (!is.null(grobs)) {
      max_gates <- sapply(grobs, length) %>% max
      grobs <- sapply(grobs, `length<-`, value = max_gates, simplify = FALSE)

      save_plotArgs <- list(
        width = min(5.0 * max_gates + 1, 200), # 200 in. is PDF maximum
        height = min(5.0 * length(grobs) + 1, 200), # 200 in. is PDF maximum
        #file = paste(gating_poster_dir, "gated-clusters-poster.pdf", sep = "/")
        filename = paste(gating_poster_dir, "gated-clusters-poster.pdf", sep = "/") # For 'grDevices::cairo_pdf()'
      )
      save_plotArgs <- utils::modifyList(save_plotArgs, save_plot..., keep.null = TRUE)

      do.call(save_plot_fun, save_plotArgs)

      ## Create a blank plot for empty grid cells (but not needed for 'cowplot::plot_grid()')
      if (FALSE) {
        blankPath <- tempfile()
        grDevices::png(file = blankPath, bg = "transparent")
        dev.control(displaylist = "enable")
        plot.new()
        blank <- grDevices::recordPlot()
        invisible(dev.off())
        unlink(blankPath)
      }

      cowplot::plot_grid(
        ## This creates a list of "recordedplot" objects:
        #plotlist = sapply(grobs %>% purrr::flatten(), function(a) if (is.null(a)) list(blank) else a),
        plotlist = sapply(grobs %>% purrr::flatten(), function(a) if (is.null(a)) list(NULL) else a),
        ncol = max_gates,
        hjust = 0, label_x = 0.01,
        labels = rep("", max_gates * length(grobs)) %>%
          `[<-`(seq(from = 1, by = max_gates, length.out = length(grobs)), names(grobs)),
        #label_colour = "darkgreen",
        label_size = 16
      ) %>% print

      dev.off()

      ## Convert PDF to PNG
      suppressWarnings(pdftools::pdf_convert(
        pdf = save_plotArgs$file,
        format = "png",
        dpi = 100,
        filenames = sprintf("%s.png", tools::file_path_sans_ext(save_plotArgs$file))
      ))
    }
  }

  if (is.null(cc)) {
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

        if (byEvent)
          searchArgs$summary... <- sm

        r <- do.call(search, searchArgs)

        if (verbose) {
          cat(". Done.", fill = TRUE); utils::flush.console()
        }

        r
      }, simplify = FALSE)
  }

  tictoc::toc()

  cc0 <- cc[sapply(cc, is.null)]
  if (length(cc0) > 0)
    warning(sprintf("Clusters %s were not found", cc0 %>% names %>% sQuote %>% paste(collapse = ", ")))

  cc1 <- cc %>% purrr::compact()

  #clusterId <- attr(x, "cluster_id")
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
  if (NCOL(d) < 3)
    return (d)

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
