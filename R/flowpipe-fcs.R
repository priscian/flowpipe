#' @export
get_channels_by_sample <- function(
  x, # Vector of directory & file paths
  keep_sans_desc = expression({
    d$channels_by_sample_full_desc %<>%
      dplyr::mutate(across(starts_with("desc"), ~ dplyr::coalesce(.x, replace = .data$name)))
  }),
  verbose = TRUE
)
{
  l0 <- keystone::psapply(x,
    function(a)
    {
      ff <- flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE,
        which.lines = 1)
      p <- flowCore::pData(flowCore::parameters(ff))

      p
    }, simplify = FALSE)

  ## Are all the channels given in the same order? If no, 'flowCore::read.flowSet()' might fail.
  ## Decide yes or no by transitivity of equality.
  if (length(l0) > 1) {
    commonChannelOrder <- sapply(l0, function(a) { a$name }, simplify = FALSE) %>% {
      flit <- .
      seq(length(flit)) %>% utils::combn(m = 2) %>% keystone::dataframe() %>% as.list %>%
        sapply(function(b) { if (length(flit[[b[1]]]) != length(flit[[b[2]]])) return (FALSE);
          all(flit[[b[1]]] == flit[[b[2]]]) }) %>% all
    }
    if (!commonChannelOrder && verbose)
      warning("Some samples have their channels ordered differently from the other samples.",
        immediate. = TRUE)
  }

  #l <- rlang::duplicate(l0, shallow = FALSE)

  d <- sapply(list(channels_by_sample = l0),
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
    }, simplify = FALSE) %>% `$<-`("channels_by_sample_full_desc", .$channels_by_sample)

  ## N.B. Use expression 'keep_sans_desc' to make changes to
  ##   'd$channels_by_sample_full_desc' before continuing.
  keystone::poly_eval(keep_sans_desc)

  channelNamesCounts <- apply(dplyr::select(d$channels_by_sample_full_desc, -name), 1, table,
    simplify = FALSE) %>% unlist
  if (!all(channelNamesCounts == NCOL(d$channels_by_sample_full_desc) - 1) && verbose)
    warning("Some channels may have multiple descriptions among samples. Run 'premessa::paneleditor_GUI()' to fix them.",
      immediate. = TRUE)

  ## N.B. Should I randomly reorder multiple names for a single channel having the same frequency?
  desc_freq_by_channel <- apply(dplyr::select(d$channels_by_sample, -name), 1,
    function(a)
    {
      table(a, useNA = "ifany") %>%
        as.data.frame %>% `colnames<-`(c("desc", "freq")) %>%
        dplyr::mutate(desc = keystone::unfactor(desc)) %>%
        dplyr::arrange(desc(freq))
    }) %>% structure(.Names = d$channels_by_sample %>% dplyr::pull(name))

  desc_remap <- desc_freq_by_channel %>%
    sapply(function(a) { if (NROW(a) > 1) a$desc[1] else NULL }) %>%
    purrr::compact() %>% unlist

  channels_by_sample <- structure(d$channels_by_sample,
    #path_map = keystone::dataframe(desc = tail(names(d$channels_by_sample), -1), path = x),
    path_map = structure(x, .Names = tail(names(d$channels_by_sample), -1)),
    desc_freq_by_channel = desc_freq_by_channel,
    common_channels = commonChannelOrder,
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
## (Slightly time-consuming, but keeps 'flowCore::read.flowSet()' from failing)
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
  reorder_columns =  # Put columns in each file in common order?
    !attr(channels_by_sample, "common_channels")
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

  pp <- keystone::psapply(names(path_map),
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


#' @export
split_pmm_by_cluster <- function(
  x, # "pmm" object from 'get_expression_subset()'
  l, # Named list of cluster sets
  fcs_dir = NULL,
  export_id_map = FALSE,
  default_colname = "orig",
  ...
)
{
  lfs <- sapply(names(l),
    function(a)
    {
      ## Split "pmm" matrix into multiple 'flowCore::flowFrames' by cluster
      cid <- l[[a]]
      if (is_invalid(cid))
        stop("PMM object has no 'cluster_id' attribute")

      eff <- flowCore::flowFrame(exprs = x %>% unclass)
      #p <- flowCore::pData(flowCore::parameters(eff))

      if (!is.matrix(cid))
        cid <- structure(as.matrix(cid), .Dimnames = list(NULL, default_colname))

      fs <- keystone::psapply(colnames(cid),
        function(bc)
        {
          b <- cid[, bc]
          lff <- flowCore::split(eff[!is.na(b), ], b[!is.na(b)])

          if (!is_invalid(fcs_dir)) {
            plyr::l_ply(seq_along(lff),
              function(bb)
              {
                fcs_path <-
                  paste(fcs_dir, a,
                    fs::path_sanitize(paste0(paste(stringr::str_trunc(as.character(bc), 31), names(lff)[[bb]], sep = "#"),
                      ".fcs"), "_"), sep = "/")
                attr(lff[[bb]], "fcs_path") <<- fcs_path
                flowCore::keyword(lff[[bb]]) <<- list(`$FIL` = basename(fcs_path))
              })
          }
          ## Check:
          # sapply(lff, function(x) { attr(x, "fcs_path") })

          lff
        }, simplify = FALSE, .parallel = FALSE)

      if (export_id_map) {
        export_id_map(
          x = x,
          export_path = paste(fcs_dir, a, a %_% "-id-map.xlsx", sep = "/")
        )
      }

      fs
    }, simplify = FALSE)
  ## Check:
  # rrapply::rrapply(lfs, f = function(x, .xname, .xparents) { attr(x, "fcs_path") }, class = "flowFrame", how = "replace")

  if (!is.null(fcs_dir)) {
    fcs_paths <- rrapply::rrapply(lfs,
      f = function(x, .xname, .xparents)
      {
        fcs_path <- attr(x, "fcs_path")

        if (!dir.exists(dirname(fcs_path)))
          dir.create(dirname(fcs_path), recursive = TRUE)

        flowCore::write.FCS(x, fcs_path, ...)

        fcs_path
      }, class = "flowFrame", how = "flatten")
  }

  lfs
}


#' @export
export_id_map <- function(
  x, # "pmm" object from 'get_expression_subset()'
  export_path
)
{
  d <- keystone::dataframe(
    id = attr(x, "id_map") %>% as.vector,
    file = attr(x, "id_map") %>% names %>% basename %>% tools::file_path_sans_ext()
  )

  if (!missing(export_path)) {
    if (!dir.exists(dirname(export_path)))
      dir.create(dirname(export_path), recursive = TRUE)

    rio::export(d, export_path, overwrite = TRUE)
  }

  d
}


#' @export
count_events <- function(
  x, # Any vector of FCS or PMM file paths
  ...
)
{
  keystone::psapply(
    x,
    function(a)
    {
      rv <- NA_integer_
      fileExt <- tools::file_ext(a)

      if (stringr::str_detect(fileExt, stringr::regex("fcs$", ignore_case = TRUE))) {
        rv <- nrow(flowCore::read.FCS(a))
      } else if (stringr::str_detect(fileExt, stringr::regex("rdata$", ignore_case = TRUE))) {
        e <- new.env()
        load(a, envir = e)

        rv <- flowCore::exprs(e$tff) %>% NROW
      }

      rv
    }, ...)
}


#' @export
as_flowFrame <- function(x, ...)
  UseMethod("as_flowFrame")


#'@export
as_flowFrame.default <- function(x, ...)
{
  as_flowFrame.matrix(x, ...)
}


#'@export
as_flowFrame.matrix <- function(x, na.rm = FALSE, name, desc, ...)
{
  x <- data.matrix(x)

  if (missing(desc))
    desc = colnames(x)

  if (missing(name)) {
    name = colnames(x)
  } else {
    colnames(x) <- name
  }

  pData <- keystone::dataframe(
    name = name,
    desc = desc,
    range = apply(apply(x, 2, range, na.rm = na.rm), 2, keystone::diff_),
    minRange = apply(x, 2, min, na.rm = na.rm),
    maxRange = apply(x, 2, max, na.rm = na.rm)
  )

  ff <- new(getClass(structure("flowFrame", package = "flowCore")),
    exprs = x,
    parameters = Biobase::AnnotatedDataFrame(pData)
  )

  ff
}

## usage:
# d <- structure(rnorm(500), dim = c(50L, 10L), dimnames = list(NULL, LETTERS[1:10]))
# ff <- as_flowFrame(d)


#'@export
as_flowFrame.data.frame <- function(x, ...)
{
  as_flowFrame.matrix(x, ...)
}


#'@export
as_flowFrame.character <- function(x, ...)
{
  d <- rio::import(x)

  ff <- as_flowFrame.matrix(d, ...)

  flowCore::keyword(ff) <- list(`$FIL` = basename(x))

  ff
}
