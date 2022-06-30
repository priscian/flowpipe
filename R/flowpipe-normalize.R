## Bead-normalization using the "premessa" functions
## A simplified drop-in replacement for 'premessa::normalize_folder()'
##   w/ additional flexibility in specifying files to process
#' @export
bead_normalize <- function(
  wd, # Directory path or vector of file paths
  output.dir.name, # For these args, see '?premessa::normalize_folder'
  ## N.B. Element names of 'beads.gates' must be full file paths. See '?premessa:::calculate_baseline'.
  beads.gates,
  beads.type, # Either "Fluidigm" or "Beta"
  baseline = NULL,
  pattern = "(?i)\\.fcs$", list_files... = list(),
  read.FCS... = list()
)
{
  ## 'pp' will be a vector of FCS file paths to be normalized to 'output.dir.name'.
  pp <- NULL
  de <- dir.exists(wd)
  if (any(de)) {
    pp <- c(pp, sapply(wd[de],
      function(a)
      {
        list_filesArgs <- list(
          path = a,
          pattern = pattern
        )
        list_filesArgs <- utils::modifyList(list_filesArgs, list_files..., keep.null = TRUE)

        do.call(keystone::list_files, list_filesArgs)
      }, simplify = FALSE)) %>% unlist(use.names = FALSE)

    pp
  }
  if (any(!de)) {
    pp <- c(pp, wd[!de][file.exists(wd[!de])])
  }

  baseline.data <- NULL
  if (is.null(baseline))
    baseline.data <- calculate_baseline(wd = pp, beads.type, files.type = "data", beads.gates)
  else
    baseline.data <- calculate_baseline(baseline, beads.type, files.type = "beads")

  if (!dir.exists(output.dir.name))
    dir.create(output.dir.name, recursive = TRUE)

  read.FCSArgs <- list(
    filename = NULL,
    transformation = "linearize",
    truncate_max_range = TRUE
  )

  ll <- lapply(names(beads.gates),
    function(f.name)
    {
      read.FCSArgs$filename <- f.name
      read.FCSArgs <- utils::modifyList(read.FCSArgs, read.FCS..., keep.null = TRUE)

      fcs <- do.call(read.FCS, ablineArgs)
      #fcs <- flowCore::read.FCS(f.name)

      beads.cols <- premessa:::find_bead_channels(fcs, beads.type)
      beads.cols.names <- premessa:::get_parameter_name(fcs, beads.cols)
      dna.col <- premessa:::find_dna_channel(fcs)
      m <- flowCore::exprs(fcs)
      beads.events <- premessa:::identify_beads(asinh(m/5), beads.gates[[f.name]], beads.cols.names, dna.col)
      beads.data <- m[beads.events, ]
      norm.res <- premessa:::correct_data_channels(m, beads.data, baseline.data, beads.cols.names)
      fcs <- flowCore::read.FCS(f.name, which.lines = 1)
      m <- NULL
      gc()
      m.normed <- norm.res$m.normed
      m.normed <- cbind(m.normed, beadDist = premessa:::get_mahalanobis_distance_from_beads(m.normed, beads.events, beads.cols.names))

      out.name <- paste(tools::file_path_sans_ext(basename(f.name)), "normalized.fcs", sep = "_")
      out.path <- paste(out.dir.path, out.name, sep = "/")
      flowCore::keyword(fcs) <- list(`$FIL` = out.name)
      premessa::write_flowFrame(premessa::as_flowFrame(m.normed, fcs), out.path)

      out.path
    }) %>% unlist(use.names = FALSE)

  ## Return value is a vector of new file paths named after the original paths
  structure(ll, .Names = names(beads.gates))
}


## A drop-in replacement for 'premessa::calculate_baseline()' that takes a vector of file paths
calculate_baseline <- function(
  wd, # Vector of file paths
  beads.type,
  files.type = c("data", "beads"),
  ## N.B. Element names of 'beads.gates' must be full file paths. See '?premessa:::calculate_baseline'.
  beads.gates = NULL
)
{
  files.type <- match.arg(files.type)
  files.list <- switch(files.type,
    data = names(beads.gates),
    beads = wd
  )
  ret <- lapply(files.list,
    function(f.name) {
      fcs <- flowCore::read.FCS(f.name)
      beads.cols.names <- premessa:::find_beads_channels_names(fcs, beads.type)
      dna.col <- premessa:::find_dna_channel(fcs)
      m <- flowCore::exprs(fcs)
      if (files.type == "data") {
        m.transformed <- asinh(m/5)
        sel <- premessa:::identify_beads(m.transformed, beads.gates[[f.name]], beads.cols.names, dna.col)
      }
      else {
        sel <- rep(TRUE, nrow(m))
      }

      return (m[sel, beads.cols.names])
    })
  ret <- Reduce("rbind", ret)
  ret <- apply(ret, 2, median)

  return (ret)
}


#' @export
prelim_bead_normalize_expr <- expression({
  bead_normalizeArgs <- list(
    wd = pp,
    output.dir.name = paste(data_dir, "normalized", sep = "/")
  )
  bead_normalizeArgs <- utils::modifyList(bead_normalizeArgs, bead_normalize..., keep.null = TRUE)

  r <- do.call(bead_normalize, bead_normalizeArgs)
  pp[names(r)] <- r
  pp <- structure(pp, .Names = pp %>% as.vector)
})
