## Bead-normalization using the "CATALYST" functions
CATALYST_normCytof <- function(
  x, # flowCore::flowFrame
  beads,
  prepData... = list(), # Other arguments to 'CATALYST::prepData()'
  normCytof... = list(), # Other arguments to 'CATALYST::normCytof()'
  seed = 666 # Try changing this if 'CATALYST::normCytof()' fails
)
{
  prepDataArgs <- list(
    x = flowCore::flowSet(x),
    transform = TRUE,
    cofactor = 5,
    by_time = FALSE,
    FACS = FALSE
  )
  prepDataArgs <- utils::modifyList(prepDataArgs, prepData..., keep.null = TRUE)

  sce <- do.call(CATALYST::prepData, prepDataArgs)

  sce <- CATALYST::prepData(
    flowCore::flowSet(x),
    transform = TRUE,
    cofactor = prepDataArgs$cofactor,
    by_time = FALSE,
    FACS = FALSE
  )

  normCytofArgs <- list(
    x = sce,
    beads = beads,
    remove_beads = FALSE,
    overwrite = TRUE,
    transform = FALSE,
    cofactor = prepDataArgs$cofactor,
    plot = FALSE,
    verbose = TRUE
  )
  normCytofArgs <- utils::modifyList(normCytofArgs, normCytof..., keep.null = TRUE)

  if (!is.null(seed))
    set.seed(seed)
  res <- do.call(CATALYST::normCytof, normCytofArgs)

  ### Restore updated expression matrix to the flowFrame 'x'

  restored_colnames <- SummarizedExperiment:::rowData(res$data)$channel_name
  colData <- SingleCellExperiment::colData(res$data)
  new_exprs <- res$data@assays@data$exprs %>% t %>%
    rev_asinh(0, 1/prepDataArgs$cofactor) %>%
    `colnames<-`(restored_colnames)
  rm(res)
  ## Restore "non-mass" channels from 'SingleCellExperiment' object 'sce'
  new_exprs <- cbind(new_exprs,
    SingleCellExperiment:::int_colData(sce) %>% `[`(sapply(., is.numeric)) %>%
      data.matrix) %>% `[`(, flowCore::colnames(x))
  ## Remove bead & doublet events identified by 'CATALYST::normCytof()'
  keep_index <- !(colData$remove | colData$is_bead)
  new_exprs <- new_exprs[keep_index, ]

  flowCore::exprs(x) <- new_exprs

  x
}


#' @export
bead_normalize_single <- function(
  input_path,
  output_dir = ".", create_output_dir = TRUE,
  outfile_prefix = "", outfile_suffix = "", # also possibly 'NULL'
  ## "dvs" (for bead masses 140, 151, 153, 165, 175) or
  ##   "beta" (for masses 139, 141, 159, 169, 175) or numeric vector of masses:
  beads = c("dvs", "beta"),
  ...
)
{
  ff <- flowCore::read.FCS(input_path, transformation = FALSE, truncate_max_range = FALSE)

  ff <- CATALYST_normCytof(x = ff, beads = beads, ...)

  if (create_output_dir && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  fcsFileName <-
    sprintf(paste0("%s", basename(input_path), "%s.fcs"), outfile_prefix, outfile_suffix)
  fcsFilePath <- paste(output_dir, fcsFileName, sep = "/")
  flowCore::keyword(ff) <- list(`$FIL` = basename(fcsFilePath))
  cat(sprintf("Saving FCS file %s...", fcsFileName)); utils::flush.console()
  flowCore::write.FCS(ff, filename = fcsFilePath)
  cat(". Done.", fill = TRUE)

  fcsFilePath
}


#' @export
bead_normalize <- function(
  input_path, # Any vector of FCS files
  output_dir = ".", # Vector of output directories, recycled to 'length(input_path)',
  ... # Arguments passed on to 'bead_normalize_single()'
)
{
  outputDirs <- structure(rep(output_dir, length.out = length(input_path)), .Names = names(input_path))

  pp <- keystone::psapply(seq_along(input_path),
    function(a)
    {
      bead_normalize_single(input_path = input_path[a], output_dir = outputDirs[a], ...)
    }, simplify = FALSE)

  ret_val <- structure(pp %>% unlist %>% as.vector)

  ret_val
}
