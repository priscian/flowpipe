## Preliminary Tasks Possibly Needed before Further Processing
##
## 1. Compensation
## 2. Bead normalization (Use "premessa" interface)
## 3. Debarcoding/deconvolution (CyTOF)
## 4. Check channel names/descriptions
## 5. Introduce metadata
#' @export
preprocess <- function(
  x, # Vector of directory & file paths
  data_dir,
  queue = list(
    ## E.g.
    # keystone::expr_sub(prelim_compensate_expr, list(compensate... = list())),
    # keystone::expr_sub(prelim_bead_normalize_expr, list(bead_normalize... = list())),
    # keystone::expr_sub(prelim_debarcode_expr, list(debarcode... = list()))
  ),
  create_data_dir = TRUE,
  pattern = "(?i)\\.fcs$", list_files... = list(),
  result_fcs_path = paste(data_dir, "preprocessed-files.xlsx", sep = "/"), # Make NULL to skip
  ...
)
{
  if (create_data_dir && !dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)

  ## 'pp' will be a vector of FCS file paths to be preprocessed.
  pp0 <- NULL
  de <- dir.exists(x)
  if (any(de)) {
    pp0 <- c(pp0, sapply(x[de],
      function(a)
      {
        list_filesArgs <- list(
          path = a,
          pattern = pattern
        )
        list_filesArgs <- utils::modifyList(list_filesArgs, list_files..., keep.null = TRUE)

        do.call(keystone::list_files, list_filesArgs)
      }, simplify = FALSE)) %>% unlist(use.names = FALSE)

    pp0
  }
  if (any(!de)) {
    pp0 <- c(pp0, x[!de][file.exists(x[!de])])
  }

  pp <- structure(pp0, .Names = pp0)

  ### Run through preprocessing queue

  for (i in queue) {
    keystone::poly_eval(i)
  }

  if (!is.null(result_fcs_path))
    rio::export(dataframe(path = pp), result_fcs_path, overwrite = TRUE)

  pp
}


## Use to invoke function 'debarcode::debarcode()' during preprocessing.
#' @export
prelim_debarcode_expr <- expression({
  debarcodeArgs <- list(
    input_path = pp,
    output_dir = paste(data_dir, "deconvoluted", sep = "/")
  )
  debarcodeArgs <- utils::modifyList(debarcodeArgs, debarcode..., keep.null = TRUE)

  r <- do.call(debarcode::debarcode, debarcodeArgs)
  ## Don't include the "0" events files in this list:
  dbcFiles <- sapply(r, function(a) tail(sapply(a, tools::file_path_as_absolute), -1))
  pp <- pp %>% `[<-`(names(dbcFiles), dbcFiles) %>% unlist %>%
    structure(., .Names = .)
})
