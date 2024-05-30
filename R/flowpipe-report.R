#' @export
summarize_metadata <- function(
  dat,
  form,
  cap = "Summary of additional study data.",
  lab = "metadata",
  latex = TRUE,
  vu_summary... = list()
)
{
  vu_summaryArgs <- list(
    f = form,
    p = keystone::ptr(dat),
    summary_formula = summary,
    latex = latex,
    caption = ifelse(latex, Hmisc::latexTranslate(cap), cap),
    lab = lab,
    size = "smaller[1]",
    prn = TRUE,
    test = FALSE,
    prtest = c("P", "name"),
    exclude1 = FALSE,
    digits = 3L,
    verbose = FALSE,
    latex... = list(where = "!h")
  )
  vu_summaryArgs <- utils::modifyList(vu_summaryArgs, vu_summary..., keep.null = TRUE)

  do.call(keystone::vu_summary, vu_summaryArgs)

  nop()
}


summarize_all_clusters_single <- function(
  x, # "pmm" object from 'get_expression_subset()'
  cluster_set,
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  summary... = list()
)
{
  if (is_invalid(cluster_set))
    clusterId <- attr(x, "cluster_id")
  else
    clusterId <- cluster_set
  if (is.matrix(clusterId))
    clusterId <- clusterId[, which_cluster_set] %>% drop
  ## 'na.omit()' for when the cluster assignments incl. NAs:
  clusterId <- na.omit(clusterId) %>% as.vector %>% unique

  rr <- sapply(clusterId,
    function(a)
    {
      summaryArgs <- list(
        x = x,
        n = a,
        cluster_set = cluster_set,
        which_cluster_set = which_cluster_set,
        merged_labels = list(`-/d` = c("-", "d"), `+/++` = c("+", "++")),
        as_list = TRUE
      )
      summaryArgs <- utils::modifyList(summaryArgs, summary..., keep.null = TRUE)

      r <- do.call(summary, summaryArgs)

      r
    }, simplify = TRUE, USE.NAMES = FALSE)

  sac <- list()
  collapse <- utils::modifyList(formals(summary.pmm)["collapse"],
    summary...["collapse"], keep.null = FALSE)$collapse
  expression_level_sep <- utils::modifyList(formals(summary.pmm)["expression_level_sep"],
    summary...["expression_level_sep"], keep.null = FALSE)$expression_level_sep
  sac$list <- plyr::ldply(names(rr),
    function(a) { keystone::dataframe(cluster = a,
      phenotype = paste0(names(rr[[a]]), sapply(rr[[a]], paste, collapse = expression_level_sep), collapse = collapse)) }) %>%
    tibble::as_tibble()

  sac$table <- plyr::ldply(names(rr),
    function(a)
    {
      r <- sapply(names(rr[[a]]),
        function(b)
        {
          keystone::dataframe(antigen = b, expression = paste(rr[[a]][[b]], collapse = expression_level_sep))
        }, simplify = FALSE) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::bind_cols(cluster = a, .)

      r
    }) %>%
    tidyr::pivot_wider(id_cols = cluster, names_from = antigen, values_from = expression) %>%
    naniar::replace_with_na_if(.predicate = is.character, condition = ~ .x == "")

  sac
}


#' @export
summarize_all_clusters <- function(
  x,
  cluster_set,
  which_cluster_set = NULL,
  ...,
  callback = NULL
)
{
  if (missing(cluster_set)) {
    clusterId <- attr(x, "cluster_id")
    cluster_set <- NULL
  } else {
    clusterId <- cluster_set
  }

  if (is.null(which_cluster_set)) {
    which_cluster_set <- 1
    if (is.matrix(clusterId))
      which_cluster_set <- colnames(clusterId)
  }

  sac <- sapply(which_cluster_set,
    function(a)
    {
      #print(a)
      summarize_all_clusters_single(x, cluster_set = cluster_set, which_cluster_set = a, ...)
    }, simplify = FALSE)

  keystone::poly_eval(callback)

  sac
}

## usage:
# sac <- summarize_all_clusters(e, summary... = list(collapse = ";"))
# summarize_all_clusters_latex(sac, type = "table")
# summarize_all_clusters_latex(sac, type = "list", caption = "Blargh")


summarize_all_clusters_latex_single <- function(
  x, # list object from 'summarize_all_clusters()'
  type = c("list", "table"),
  caption = NULL, label = "label",
  longtabu = type == "list",
  latex... = list(),
  cat = FALSE
)
{
  type <- match.arg(type)
  xx <- x[[type]]

  xxx <- dplyr::mutate_all(xx, Hmisc::latexTranslate, greek = TRUE) %>%
    dplyr::rename_with(Hmisc::latexTranslate, greek = TRUE) %>%
    dplyr::mutate_at(-1, stringr::str_replace_all, pattern = "-", replacement = "--")
  i <- as.numeric(xxx$cluster)
  if (any(is.na(i))) i <- xxx$cluster
  xxx <- dplyr::arrange(xxx, i)

  latexArgs <- list(
    object = xxx,
    file = "",
    where = "!ht",
    size = "smaller[0]",
    ctable = FALSE, booktabs = TRUE,
    caption = caption, label = label,
    rowname = NULL,
    table.env = FALSE,
    center = "centering"
  )
  latexArgs <- utils::modifyList(latexArgs, latex..., keep.null = TRUE)

  l <- utils::capture.output(do.call(Hmisc::latex, latexArgs))

  if (longtabu) {
    longtabuTemplate <-
      "begin\\{longtabu\\} to \\\\textwidth {X[0.1]X}\n\\\\caption{%s\\\\label{%s}}\\\\\\\\"
    l <- l %>%
      stringr::str_replace("begin\\{tabular\\}.*$", sprintf(longtabuTemplate, caption, label)) %>%
      stringr::str_replace("\\{tabular\\}", "\\{longtabu\\}")
  }

  paste(l, collapse = "\n")
}


#' @export
summarize_all_clusters_latex <- Vectorize(summarize_all_clusters_latex_single, vectorize.args = "x", SIMPLIFY = FALSE)


## Create a document outside the main report
#' @export
make_external_latex_document <- function(
  x, # A potentially standalone LaTeX object (e.g. table, figure, etc.)
  file_path,
  docTemplate =
r"---{\documentclass{article}
\usepackage[paperwidth=50in,paperheight=50in]{geometry}
\usepackage{booktabs}
\usepackage{relsize}
\begin{document}
\pagestyle{empty}
%s
\end{document}}---"
)
{
  l <- sprintf(docTemplate, x)

  if (!dir.exists(dirname(file_path)))
    dir.create(dirname(file_path), recursive = TRUE)

  writeLines(l, con = file_path)
}
#make_external_latex_document(summarize_all_clusters_latex(sac, type = "table"), file_path = paste(report_dir, "cluster-phenotypes.tex", sep = "/"))


## Create spreadsheets of cluster cell counts by sample
#' @export
export_cluster_summary <- function(
  details, # Named list of objects from 'do_differential_expression()' & 'summarize_all_clusters()'
  spreadsheet_path,
  sep = "--",
  keep_pm_lists = FALSE,
  ...
)
{
  l1_0 <- sapply(names(details),
    function(a)
    {
      if (is.null(details[[a]][[1]]))
        return (NULL)

      sapply(details[[a]][[1]],
        function(b) { attr(b, "cell_counts") %>% unclass %>%
          as.data.frame %>% tibble::rownames_to_column(var = "cluster") },
        simplify = FALSE)
    }, simplify = FALSE)

  l1 <- sapply(names(l1_0),
    function(a)
    {
      if (is.null(l1_0[[a]]))
        return (NULL)

      r0 <- structure(l1_0[[a]], .Names = paste("events", a, names(l1_0[[a]]), sep = sep))
      rrapply::rrapply(r0,
        f = function(x)
        {
          i <- as.numeric(x$cluster)
          if (any(is.na(i))) i <- x$cluster

          dplyr::arrange(x, i)
        },
        class = "data.frame", how = "replace")
    }, simplify = FALSE) %>%
    purrr::flatten()

  l2 <- sapply(names(details), function(a)
    {
      if (is.null(details[[a]][[2]]))
        return (NULL)

      rrapply::rrapply(details[[a]][[2]],
        f = function(x, .xname)
        {
          i <- as.numeric(x$cluster)
          if (any(is.na(i))) i <- x$cluster

          if (.xname == "list" && !keep_pm_lists)
            return (NULL)

          dplyr::arrange(x, i)
        },
        class = "data.frame", how = "replace")
    }, simplify = FALSE, USE.NAMES = TRUE) %>% purrr::compact() %>%
    ## N.B. Note the braces around 'structure()' here -- they're necessary!
    {structure(
      .Data = rrapply::rrapply(., f = identity, class = "data.frame", how = "flatten"),
      .Names = rrapply::rrapply(.,
        f = function(x, .xname, .xparents) { paste(.xparents, collapse = sep) },
        class = "data.frame", how = "flatten") %>%
        as.vector %>% paste("PM", ., sep = sep)
    )}

  l3 <- sapply(names(details),
    function(a)
    {
      if (is.null(details[[a]][[3]]))
        return (NULL)

      rrapply::rrapply(details[[a]][[3]],
        f = function(x)
        {
          r <- cbind(cluster = rownames(x), t(scale(t(x))) %>% `[<-`(is.nan(.), NA) %>% as.data.frame)
          i <- as.numeric(r$cluster)
          if (any(is.na(i))) i <- r$cluster

          dplyr::arrange(r, i)
        },
        class = "matrix", how = "replace")
    }, simplify = FALSE, USE.NAMES = TRUE) %>% purrr::compact() %>%
    {structure(
      .Data = rrapply::rrapply(., f = identity, class = "data.frame", how = "flatten"),
      .Names = rrapply::rrapply(.,
        f = function(x, .xname, .xparents) { paste(.xparents, collapse = sep) },
        class = "data.frame", how = "flatten") %>%
        as.vector %>% paste("medians", ., sep = sep)
    )}

  l <- c(l1, l2, l3)

  sheetNames <- stringr::str_trunc(names(l), 31, "center")
  ## Rename sheets w/ duplicated names
  duplicateNames <- sheetNames %>% intersect(.[duplicated(.)])
  for (i in duplicateNames) {
    dupIndex <- which(sheetNames == i)
    ## Replace w/ sequential numbers:
    sheetNames[dupIndex] <-
      sapply(seq_along(dupIndex),
        function(j) sprintf("%s%01d", stringr::str_trunc(sheetNames[dupIndex[j]], 30, ellipsis = ""), j))
  }
  rio::export(l %>% `names<-`(sheetNames), spreadsheet_path, ...)

  l
}
