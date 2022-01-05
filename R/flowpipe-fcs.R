#' @export
get_channels_by_sample <- function(
  x, # Vector of directory & file paths
  keep_sans_desc = expression({
    l <- sapply(l,
      function(a)
        ## When 'desc' parameter is missing, copy 'name' in:
        a %>% dplyr::mutate(desc = dplyr::case_when(is.na(desc) ~ name, TRUE ~ desc)), simplify = FALSE)
  })
)
{
  l0 <- sapply(x,
    function(a)
    {
      ff <- flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE)
      p <- flowCore::pData(flowCore::parameters(ff))

      p
    }, simplify = FALSE)

  l <- rlang::duplicate(l0, shallow = FALSE)
  ## N.B. Use expression 'keep_sans_desc' to make changes to 'l' before continuing.
  plinth::poly_eval(keep_sans_desc)

  d <- sapply(list(channels_by_sample = l0, channels_by_sample_full_desc = l),
    function(x) {
      i <- 1
      d <- c(
        head(x, 1)[[1]] %>% dplyr::select(name, desc) %>% dplyr::rename(desc_01 = desc) %>% list,
        tail(x, -1)
      ) %>%
      purrr::reduce(
        function(a, b)
        {
          i <<- i + 1
          newName <- sprintf("desc_%02d", i)

          dplyr::full_join(
            a,
            dplyr::select(b, name, desc) %>% dplyr::rename(!!newName := desc),
            by = "name")
        })

      d %>% dplyr::mutate(across(.cols = everything(), .fns = as.vector))
    }, simplify = FALSE)

  channelNamesCounts <- apply(dplyr::select(d$channels_by_sample_full_desc, -name), 1, table, simplify = FALSE) %>% unlist
  if (!all(channelNamesCounts == NCOL(d$channels_by_sample_full_desc) - 1))
    warning("Some channels may have multiple descriptions among samples. Run 'premessa::paneleditor_GUI()' to fix them.",
      immediate. = TRUE)

  desc_freq_by_channel <- apply(dplyr::select(d$channels_by_sample, -name), 1,
    function(a)
    {
      table(a, useNA = "ifany") %>%
        as.data.frame %>% `colnames<-`(c("desc", "freq")) %>%
        dplyr::mutate(desc = plinth::unfactor(desc)) %>%
        dplyr::arrange(desc(freq))
    }) %>% structure(.Names = d$channels_by_sample %>% dplyr::pull(name))

  desc_remap <- desc_freq_by_channel %>%
    sapply(function(a) { if (NROW(a) > 1) a$desc[1] else NULL }) %>%
    purrr::compact() %>% unlist

  channels_by_sample <- structure(d$channels_by_sample,
    #path_map = plinth::dataframe(desc = tail(names(d$channels_by_sample), -1), path = x),
    path_map = structure(x, .Names = tail(names(d$channels_by_sample), -1)),
    desc_freq_by_channel = desc_freq_by_channel,
    #channel_names_counts = channelNamesCounts,
    desc_remap = desc_remap,
    channels_by_sample_full_desc = d$channels_by_sample_full_desc
  )
  class(channels_by_sample) <- c("cbs", class(channels_by_sample)) %>% unique

  channels_by_sample
}


#' @export
any_noncommon_descs <- function(
  channels_by_sample # "cbs" object
)
{
  if (!inherits(channels_by_sample, "cbs"))
    stop("Input must be an object returned by 'get_channels_by_sample()'")

  channelNamesCounts <- apply(
    dplyr::select(attr(channels_by_sample, "channels_by_sample_full_desc"), -name),
    1, table, simplify = FALSE) %>% unlist

  structure(!all(channelNamesCounts ==
    NCOL(attr(channels_by_sample, "channels_by_sample_full_desc")) - 1),
    channel_names_counts = channelNamesCounts)
}


## Drop-in update of 'premessa::rename_fcs_parameters_name()' for regular 'flowCore::flowFrame's
#' @export
rename_fcs_parameters_name <- function(
  fcs, # flowCore::flowFrame
  names.map # Named vector of channels to rename
)
{
  ret <- fcs
  old.names <- flowCore::colnames(ret) %>% structure(., .Names = .)
  new.names <- old.names %>% `[<-`(names(names.map), names.map) %>% as.character
  stopifnot(!any(duplicated(new.names)))
  flowCore::colnames(ret) <- new.names

  spill.keyword <- grep("SPILL", names(flowCore::keyword(ret)), value = TRUE)
  if (length(spill.keyword) > 0) {
    m <- flowCore::keyword(ret)[[spill.keyword]]
    stopifnot(is.matrix(m))
    colnames(m) <- colnames(m) %>% structure(., .Names = .) %>%
      `[<-`(names(names.map), names.map) %>% as.character
    flowCore::keyword(ret) <- structure(list(m), .Names = spill.keyword)
  }

  return (ret)
}


append_fcs_blank_cols <- function(
  x, # flowCore::flowFrame
  new_cols,
  value = 0.0 # Or maybe 'NA_real_'
)
{
  flowCore::fr_append_cols(
    x,
    structure(
      matrix(rep(value, flowCore::nrow(x) * length(new_cols)), ncol = length(new_cols)),
      .Dimnames = list(NULL, new_cols)
    )
  )
}


## In-situ reshuffling of FCS columns to a common order
## (Slightly time-consuming, but keeps 'flowCore::read.flowSet' from failing)
#' @export
reorder_fcs_common_columns <- function(
  x, # Vector of file paths
  column_order = NULL, # If NULL, use the column ordering of first file in 'x'
  ...
)
{
  columnOrder <- column_order

  pp <- sapply(x,
    function(a)
    {
      ff <- flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE)

      if (is.null(columnOrder))
        columnOrder <<- flowCore::colnames(ff)

      if (all(columnOrder == flowCore::colnames(ff)))
        return (a)

      flowCore::write.FCS(ff[, columnOrder], a)

      a
    }, simplify = TRUE)

  pp
}


#' @export
rename_fcs_parameters_name_desc <- function(
  channels_by_sample,
  ## For non-automatic renaming, provide a named-vector argument to 'desc_remap':
  desc_remap = attr(channels_by_sample, "desc_remap"),
  name_remap = NULL, # Named vector of channels to rename, NA to remove channel
  output_dir = NULL,
  outfile_suffix = "_renamed",
  reorder_columns = TRUE # Put columns in each file in common order?
)
{
  cbs <- channels_by_sample
  path_map <- attr(cbs, "path_map")

  ## Rename descriptions according to map 'desc_remap'
  cbs_update <- cbs %>% tibble::as_tibble() %>% dplyr::mutate(
    across(.cols = starts_with("desc"),
      .fns = function(a)
      {
        structure(a, .Names = .$name) %>%
          `[<-`(names(desc_remap), desc_remap) %>%
          as.vector
      })
  ) %>% as.data.frame

  pp <- sapply(names(path_map),
    function(a)
    {
      descAbo <- structure(cbs[[a]], .Names = cbs$name)
      descNew <- structure(cbs_update[[a]], .Names = cbs_update$name)
      commonNames <- cbs_update$name %>% structure(., .Names = .)
      path <- path_map[a] %>% as.vector

      if (!identical(descAbo, descNew) || !is.null(name_remap)) { # Need to rewrite files w/ updated descriptions
        ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)

        if (!is.null(output_dir)) {
          if (!dir.exists(output_dir))
            dir.create(output_dir, recursive = TRUE)

          path <- sprintf("%s/%s%s.fcs", output_dir, tools::file_path_sans_ext(basename(path)), outfile_suffix)
        }

        ## Update channel descriptions
        p <- flowCore::pData(flowCore::parameters(ff))
        p$desc <- descNew[p$name] %>% as.vector
        flowCore::parameters(ff) <- as(p, "AnnotatedDataFrame")
        ## Might be able to do this w/ 'flowCore::markernames()'...?

        ## Rename/remove channels according to map 'name_remap'; NA means remove
        if (!is.null(name_remap)) {
          renameNameMap <- name_remap[!is.na(name_remap)]
          removeNameMap <- name_remap[is.na(name_remap)]

          if (!is_invalid(renameNameMap)) {
            if (!is_invalid(missingCols <- setdiff(names(renameNameMap), flowCore::colnames(ff))))
              ff <- append_fcs_blank_cols(ff, missingCols)

            flowCore::colnames(ff) <- flowCore::colnames(ff) %>% structure(., .Names = .) %>%
              `[<-`(names(renameNameMap), renameNameMap)

            commonNames <- commonNames %>% `[<-`(names(renameNameMap), renameNameMap) %>%
              as.vector %>% structure(., .Names = .)
          }

          if (length(removeNameMap) > 0) {
            ff <- ff[, setdiff(flowCore::colnames(ff), names(removeNameMap))]

            commonNames <- commonNames[setdiff(names(commonNames), names(removeNameMap))]
          }
        }

        if (!is_invalid(missingCols <- setdiff(commonNames %>% as.vector, flowCore::colnames(ff)))) {
          ## Add missing columns to flowFrame
          ff <- append_fcs_blank_cols(ff, missingCols)
        }

        ## Update channel descriptions again in case a blank channel nevertheless has a description
        p <- flowCore::pData(flowCore::parameters(ff))
        p$desc <- descNew[p$name] %>% as.vector
        flowCore::parameters(ff) <- as(p, "AnnotatedDataFrame")

        flowCore::keyword(ff) <- list(`$FIL` = basename(path))

        flowCore::write.FCS(ff, path)
      }

      path
    }, simplify = TRUE)


  if (reorder_columns) {
    reorder_fcs_common_columns(pp)
  }

  pp
}
