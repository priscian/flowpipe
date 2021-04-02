#' @export
find_cell_ontologies <- function
(
  x, # "pmm" object from 'get_expression_subset()'
  ontology_correction_map = c(), # E.g. ontology_correction_map = c(hCD45 = "CD45", CD235ab = "CD235a")
  flowCL_temp_dir = NULL, # Set to e.g. "D:/" to run 'flowCL::flowCL()'
  ...
)
{
  if (is.null(attr(x, "cluster_id"))) {
    attr(x, "cluster_id") <- rep(1, NROW(x))
    warning("PMM object 'x' has no 'cluster_id' attribute")
  }

  channels <- colnames(x)[stringr::str_detect(colnames(x), excluded_channels_re, negate = TRUE)]
  e <- x[, channels, drop = FALSE] %>%
    cbind(cluster = attributes(x)$cluster_id)

  mem <- MEM::MEM(e, ...)
  mem_scores <- mem$MEM_matrix[[1]]
  mem_scores <- mem_scores[order(as.numeric(rownames(mem_scores))), ]
  colnames(mem_scores) <- gsub("(\\+|\\-)", "", colnames(mem_scores))

  threshold <- median(abs(as.vector(mem_scores)), na.rm = TRUE)
  labels <- list(); ExpMrkrLst <- NULL
  plyr::l_ply(seq(NROW(mem_scores)),
  function(i)
  {
    m <- mem_scores[i, abs(mem_scores[i, ]) > 0]
    lab <- rep(NA_character_, length(m)); names(lab) <- names(m)
    #lab[m < -threshold] <- "-- "; is.na(m) <- m < -threshold
    lab[m < -threshold] <- "lo "; is.na(m) <- m < -threshold
    lab[m < 0] <- "- "; is.na(m) <- m < 0
    #lab[m > threshold] <- "++ "; is.na(m) <- m > threshold
    lab[m > threshold] <- "hi "; is.na(m) <- m > threshold
    lab[m > 0] <- "+ "; is.na(m) <- m > 0
    namesLab <- names(lab)
    plyr::l_ply(names(ontology_correction_map),
      function(a) { namesLab[grepl(a, namesLab, ignore.case = TRUE, perl = TRUE)] <<- ontology_correction_map[a]; NULL })
    names(lab) <- namesLab
    lab <- paste0(names(lab), trimws(lab))
    ExpMrkrLst <<- c(ExpMrkrLst, lab)

    labels[[rownames(mem_scores)[i]]] <<- trimws(paste(lab, collapse = ""))
  })

  flowCL_results <- NULL
  if (!is.null(flowCL_temp_dir)) {
    currentDir <- getwd()
    setwd(flowCL_temp_dir) # Change to this directory to avoid overlong path names.
    browser()
    flowCL_results <- flowCL::flowCL(MarkerList = labels, ExpMrkrLst = NULL, Verbose = FALSE, KeepArch = FALSE, VisualSkip = TRUE)
    setwd(currentDir)
  }

  list(mem = mem, mem_scores = mem_scores, flowCL_markerlist = labels, flowCL_results = flowCL_results)
}


#' @export
print_ontologies <- function()
{
  print(as.vector(ontology_results <- sapply(ontologies$Cell_Labels, function(x) x[[1]])))

  assign("ontology_results", ontology_results, envir = globalenv())
}
## N.B. This must go in the script file!
#print_ontologies()
