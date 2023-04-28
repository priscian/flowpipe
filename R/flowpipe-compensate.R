#' @export
compensate <- function(
  x, # Vector of file paths
  output_dir = ".", create_output_dir = TRUE,
  outfile_prefix = "", outfile_suffix = "_compensated", # also possibly 'NULL'
  ...
)
{
  r <- keystone::psapply(pp,
    function(p)
    {
      ff <- flowCore::read.FCS(p, transformation = FALSE, truncate_max_range = FALSE)
      hasSpill <- tryCatch({ flowCore::spillover(ff); TRUE }, error = function(e) FALSE)
      if (hasSpill) {
        comp <- flowCore::compensation(flowCore::spillover(ff))
        ff <- flowCore::compensate(ff, comp)

        fcsFileName <-
          sprintf(paste0("%s", basename(tools::file_path_sans_ext(p)), "%s.fcs"), outfile_prefix, outfile_suffix)
        fcsFilePath <- paste(output_dir, fcsFileName, sep = "/")
        flowCore::keyword(ff) <- list(`$FIL` = basename(fcsFilePath))

        if (create_output_dir && !dir.exists(output_dir))
          dir.create(output_dir, recursive = TRUE)

        cat(sprintf("Saving FCS file %s...", fcsFileName)); utils::flush.console()
        flowCore::write.FCS(ff, filename = fcsFilePath)
        cat(". Done.", fill = TRUE)

        fcsFilePath
      } else { # I.e. if there's no spillover matrix ...
        p
      }
    }, simplify = FALSE)

  ## 'r' is a vector of new file paths named after the original paths
  r
}
