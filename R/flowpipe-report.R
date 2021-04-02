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
    p = plinth::ptr(dat),
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

  do.call(plinth::vu_summary, vu_summaryArgs)

  nop()
}


#' @export
summarize_all_clusters <- function(
  x, # "pmm" object from 'get_expression_subset()',
  summary... = list(),
  callback = NULL
)
{
  cluster_id <- attr(x, "cluster_id") %>% unique

  r <- sapply(cluster_id,
    function(a)
    {
      summaryArgs <- list(
        x = x,
        n = a,
        as_list = TRUE
      )
      summaryArgs <- utils::modifyList(summaryArgs, summary..., keep.null = TRUE)

      do.call(summary, summaryArgs)
    }, simplify = FALSE) %>%
    purrr::flatten()

  sac <- list()
  collapse <- utils::modifyList(formals(summary.pmm)["collapse"], summary...["collapse"], keep.null = FALSE)
  sac$list <- plyr::ldply(names(r), function(a) plinth::dataframe(cluster = a,
    phenotype = stringr::str_flatten(r[[a]], collapse = collapse))) %>%
    tibble::as_tibble()

  sac$table <- r %>% dplyr::bind_rows(.id = "cluster") %>%
    dplyr::mutate_at(-1, stringr::str_extract, pattern = "(-|\\+)$")

  plinth::poly_eval(callback)

  sac
}

## usage:
# sac <- summarize_all_clusters(e, summary... = list(collapse = ";"))
# summarize_all_clusters_latex(sac, type = "table")
# summarize_all_clusters_latex(sac, type = "list", caption = "Blargh")


#' @export
summarize_all_clusters_latex <- function(
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
    dplyr::mutate_at(-1, stringr::str_replace_all, pattern = "-", replacement = "--") %>%
    dplyr::arrange(as.numeric(cluster))
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


## Create a document outside the main report
#' @export
make_external_latex_document <- function(
  x, # A potentially standalone LaTeX object (e.g. table, figure, etc.)
  file_path,
  docTemplate =
r"---{\documentclass{article}
\usepackage[paperwidth=25in,paperheight=35in]{geometry}
\usepackage{booktabs}
\usepackage{relsize}
\begin{document}
\pagestyle{empty}
%s
\end{document}}---"
)
{
  l <- sprintf(docTemplate, x)

  writeLines(l, con = file_path)
}
#make_external_latex_document(summarize_all_clusters_latex(sac, type = "table"), file_path = paste(report_dir, "cluster-phenotypes.tex", sep = "/"))
