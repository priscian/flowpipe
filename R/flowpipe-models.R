## Same as 'pwr::ES.h()'; used to find effect size between proportions.
## h = 0.20 small, h = 0.50 medium, h = 0.80 large
cohens_h <- function (p1, p2)
{
  2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
}


### Differential expression

do_differential_expression_single <- function(
  x, # "pmm" object from 'get_expression_subset()'
  m, # metadata data frame
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  model_formula,
  id_map_re = ".*",
  metadata_id_var = "id",
  glmQLFit... = list()
)
{
  # id_map <- structure(attr(x, "id_map"),
  #   .Names = names(attr(x, "id_map")) %>% basename %>% stringr::str_extract(id_map_re))
  id_map <- make_sample_id_map(x, id_map_re)
  ids <- names(id_map)[x[, "id"] %>% as.vector]

  clusterId <- attr(x, "cluster_id")
  if (is.matrix(clusterId))
    clusterId <- clusterId[, which_cluster_set] %>% drop

  ## Ignore events not associated w/ any cluster or whose sample isn't in the metadata:
  #keepEvents <- !is.na(clusterId) & ids %in% unique(m$id)
  keepEvents <- ids %in% unique(m$id)
  ids <- ids[keepEvents]
  clusterId <- clusterId[keepEvents]

  cell_counts_abo <- table(clusterId, ids)
  ## Zero counts for all clusters for any sample will cause the model to fail:
  if (all(colSums(cell_counts_abo) == 0)) {
    cell_counts <- NULL
  } else {
    cell_counts <- cell_counts_abo[, colSums(cell_counts_abo) != 0, drop = FALSE]
    if (any(colSums(cell_counts_abo) == 0)) {
      warning(sprintf("No counts for all clusters for all samples: %d samples -> %d",
        NCOL(cell_counts_abo), NCOL(cell_counts)), immediate. = TRUE)
    }
    ## Convert to percentages if samples are of uneven sizes (i.e. assume always)
    cell_counts_percent <- (cell_counts /
      matrix(rep(table(ids)[colnames(cell_counts)], NROW(cell_counts)), nrow = NROW(cell_counts), byrow = TRUE)) * 100
  }

  fit <- tryCatch({
    if (is.null(cell_counts))
      stop("No counts for all clusters for all samples")

    #dge <- edgeR::DGEList([cell_counts|cell_counts_percent], lib.size = table(ids))
    dge <- edgeR::DGEList(cell_counts_percent)
    ## N.B. 'left_join()' only here so that the design matrix matches up w/ the "cell_counts" table:
    design <- stats::model.matrix(model_formula,
      data = m %>% dplyr::left_join(structure(keystone::dataframe(rownames(dge$samples)), .Names = metadata_id_var), .))
    y <- edgeR::estimateDisp(y = dge, design = design)

    glmQLFitArgs <- list(
      y = y,
      design = design,
      ## N.B. Necessary to prevent occasional error in 'edgeR::glmQLFTest()' from result var.post = numeric(0):
      robust = FALSE
    )
    glmQLFitArgs <- utils::modifyList(glmQLFitArgs, glmQLFit..., keep.null = TRUE)

    fit <- structure(do.call(edgeR::glmQLFit, glmQLFitArgs),
      cell_counts = cell_counts, cell_counts_percent = cell_counts_percent)

    fit
  }, error = function(e)
    {
      message(sprintf("\nError modeling cluster set \"%s\": ", which_cluster_set), e$message)
      flush.console()

      ## Can't calculate Cohen's h if all events are in a single group:
      if (NCOL(cell_counts_percent) < 2) {
        message(sprintf("\nError calc'ing effect size for cluster set \"%s\"", which_cluster_set))
        flush.console()

        return (NULL)
      }

      group_combn <- t(utils::combn(colnames(cell_counts_percent), 2))
      ch <- apply(group_combn, 1,
        function(a)
        {
          plyr::alply(cell_counts_percent[, a, drop = FALSE], 1,
            function(b) { b <- b/100; cohens_h(b[1], b[2]) }) %>%
            unlist(use.names = FALSE) %>% `names<-`(rownames(cell_counts_percent))
        }, simplify = FALSE) %>%
          {
            structure(., .Names = plyr::alply(group_combn, 1,
              function(a) paste(keystone::backtick(a), collapse = "-")) %>%
              unlist(use.names = FALSE),
              class = "cohens_h")
          }

      return (ch)
    })

  fit
}


#' @export
do_differential_expression <- function(
  x,
  which_cluster_set = NULL,
  ...
)
{
  clusterId <- attr(x, "cluster_id")

  if (is.null(which_cluster_set)) {
    which_cluster_set <- 1
    if (is.matrix(clusterId))
      which_cluster_set <- colnames(clusterId)
  }

  fits <- sapply(which_cluster_set,
    function(a)
    {
      print(a)
      do_differential_expression_single(x, which_cluster_set = a, ...)
    }, simplify = FALSE) %>% purrr::compact()

  fits
}


test_contrasts_single <- function(
  fit,
  ## A helpful reminder of how to code contrasts: https://rcompanion.org/rcompanion/h_01.html
  ## https://rstudio-pubs-static.s3.amazonaws.com/79395_b07ae39ce8124a5c873bd46d6075c137.html
  ## For interactions, the latter w/ continuous variables, see:
  ##   http://bioconductor.org/packages/3.12/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  ##   https://support.bioconductor.org/p/87252/
  contrasts,
  include_other_vars = FALSE,
  alpha = 0.05,
  min_effect_size = 0.20,
  ...
)
{
  if (inherits(fit, "cohens_h")) {
    return (list(sig_results =
      sapply(fit, function(a) a[abs(a) >= min_effect_size], simplify = FALSE)))
  }

  if (is.logical(include_other_vars) && include_other_vars) {
    ov <-
      setdiff(colnames(fit$design)[stringr::str_detect(colnames(fit$design), "(Intercept)", negate = TRUE)],
        names(contrasts))
    contrasts <- utils::modifyList(
      contrasts,
      structure(as.list(ov), .Names = ov)
    )
  }

  res <- contrasts %>% plyr::llply(
    function(a)
    {
      glmQLFTestArgs <- utils::modifyList(
        list(glmfit = fit),
        structure(list(a), .Names = ifelse(is.matrix(a), "contrast", "coef")),
        keep.null = TRUE)
      do.call(edgeR::glmQLFTest, glmQLFTestArgs)
    })

  res_sig <- plyr::llply(res, function(a) a %>% edgeR::topTags(Inf, ...) %>% as.data.frame %>%
    dplyr::filter(dplyr::across(tidyselect::last_col(), ~ . < alpha)))

  inference <- list(results = res, sig_results = res_sig)

  inference
}


#' @export
test_contrasts <- function(
  fits,
  contrasts,
  ...
)
{
  inference <- sapply(fits,
    function(fit)
    {
      contrasts <- sapply(contrasts, function(b) keystone::poly_eval(b), simplify = FALSE)

      test_contrasts_single(fit, contrasts, ...)
    }, simplify = FALSE)

  inference
}
