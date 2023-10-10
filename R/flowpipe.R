#' @importFrom magrittr %>%
#' @importFrom tribe %@>% %<@>%
#' @importFrom keystone %_% cordon is_invalid dataframe %nin% nop
#' @importFrom magrittr %>%

## http://cytoforum.stanford.edu/viewtopic.php?f=3&t=1026
##
## 1. Acquisition of raw data (hot off the machine, no additional processing)
## 2. Normalization. Either MATLAB/Finck style, or Fluidigm-style.
## 3. Debarcoding, if applicable.
## 4. Gating down to live intact singlets (removal of debris, beads, dead cells,
##   remaining cell-cell and bead-cell doublets).
## 5. Depending on your analysis pipeline, often export of new FCS files from the gated live intact singlets
##   (if using Cytobank or a FlowJo plug-in, you can usually specify the population to cluster without having to re-export).
## 6. Clustering.


### Apply silhouette cutoffs to flow data.

#' @export
find_plus_minus_by_channel <- function(
  x, # flowCore::flowFrame
  bins = 4, intensity_levels = c("-", "d", "+", "++"),
  zero_threshold = 0.90,
  excluded_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  pmm_channels = NULL,
  multisect... = list(),
  resample_on_error = 15,
  verbose = TRUE
)
{
  if (!is.null(resample_on_error)) {
    if (is.logical(resample_on_error)) {
      if (resample_on_error) resample_on_error <- 5
      else resample_on_error <- 1
    }
  } else {
    resample_on_error <- 1
  }

  channels <- pmm_channels
  if (is.null(channels))
    channels <- stringr::str_subset(flowCore::colnames(x), excluded_channels_re, negate = TRUE)

  e <- flowCore::exprs(x)
  colNames <- colnames(e)
  cutoffs <- list(); starting_random_seed <- NULL
  i <- 1
  pmm <- plyr::alply(e, 2,
    function(b)
    {
      sampleName <- tools::file_path_sans_ext(basename(flowCore::description(x)$FILENAME))
      if (verbose) {
        cat(sprintf("Processing sample %s, channel %s...", sampleName, colNames[i])); utils::flush.console()
      }

      if (colNames[i] %nin% channels) {
        rr <- rep(NA_integer_, length(b)); names(rr) <- colNames[i]
      } else {
        ## If too many zeros etc., might be a blank channel.
        if (
          (sum(b == 0, na.rm = TRUE) / length(b) >= zero_threshold) ||
          all(is.na(b)) ||
          isTRUE(all.equal(sd(b), 0))
        ) {
          cutoff <- rep(NA_real_, bins - 1)

          r <- rep(NA_integer_, length(b))
        } else {
          multisectArgs <- list(
            x = b %>% `attr<-`("main_title", sprintf("%s|%s", sampleName, colNames[i])),
            bins = bins,
            plot_cutoff = TRUE,
            random_seed = 666
          )
          multisectArgs <- utils::modifyList(multisectArgs, multisect...)

          if (is.null(starting_random_seed) && !is.null(multisectArgs$random_seed))
            starting_random_seed <<- multisectArgs$random_seed

          cutoff <- rep(NA_real_, bins - 1) %>%
            `attr<-`("random_seed", multisectArgs$random_seed)

          j <- 0
          while (TRUE) {
            if (j >= resample_on_error) {
              cat("\nError: Too many restarts"); flush.console()
              cutoff <- rep(NA_real_, length(bins - 1)) %>%
                `attr<-`("random_seed", multisectArgs$random_seed)

              break
            }

            tryCatch({
              cutoff <- do.call(multisect, multisectArgs)

              ## Try to deal w/ errors resulting from quasi-successful 'multisect()' call.
              if (any(duplicated(cutoff))) {
                cat("\n  Warning: Duplicate cutoffs. Resampling..."); flush.console()

                j <- j + 1
                multisectArgs$random_seed <- attr(cutoff, "random_seed") + 1

                next
              }

              if (any(cutoff %in% range(b))) {
                cat("\n  Warning: Cutoff at data boundary. Resampling..."); flush.console()

                j <- j + 1
                multisectArgs$random_seed <- attr(cutoff, "random_seed") + 1

                next
              }

              break
            }, error =
              function(e)
              {
                message("\nError: ", e$message); flush.console()

                j <<- j + 1

                ## Try bumping up random seed & rerunning
                multisectArgs$random_seed <<- multisectArgs$random_seed + 1
              }
            )
          }

          ## On utter failure of 'multisect()', estimate cutoffs differently
          if (all(is.na(cutoff))) {
            ## Use mean as estimate of bisection points; cf. algorithm in file "flowpipe-silhouette.R"
            cat("\n  Warning: Falling back on mean estimates of cutpoints..."); flush.console()

            d <- stats::dist(b %>% `attributes<-`(NULL))

            cutoff1 <- mean(b, na.rm = TRUE)
            cutoff2 <- cutoff3 <- NULL

            if (bins != 2) {
              t1 <- mean(b[b > cutoff1], na.rm = TRUE)
              if (bins == 3) {
                cluster <- .bincode(b, sort(c(-Inf, cutoff1, t1, +Inf), decreasing = FALSE))
                ss <- cluster::silhouette(cluster, d) # Original R code (now in C++): cluster:::silhouette.default.R
                s1 <- mean(ss[, 3])
              }

              t2 <- mean(b[b <= cutoff1], na.rm = TRUE)
              if (bins == 3) {
                cluster <- .bincode(b, sort(c(-Inf, cutoff1, t2, +Inf), decreasing = FALSE))
                ss <- cluster::silhouette(cluster, d)
                s2 <- mean(ss[, 3])

                if (s2 > s1) cutoff2 <- t2
                else cutoff2 <- t1
              } else {
                cutoff2 <- t1; cutoff3 <- t2
              }

              cutoff <- sort(c(cutoff1, cutoff2, cutoff3), decreasing = FALSE)
              attr(cutoff, "random_seed") <- multisectArgs$random_seed
            }

            if (multisectArgs$plot_cutoff) {
              plot(stats::density(b), main = sprintf("%s|%s", sampleName, colNames[i]), cex.main = 0.8)
              keystone::vline(sprintf("%.2f", cutoff), abline... = list(col = "red"), text... = list(y = keystone::cp_coords()$y))
            }
          }

          r <- Hmisc::cut2(b, cutoff)
        }
        #length(intensity_levels) <- length(cutoff) + 1
        levels(r) <- intensity_levels
        cutoffs <<- c(cutoffs, list(cutoff))

        rr <- r; names(rr) <- colNames[i]
      }

      if (verbose) {
        cat(". Done.", fill = TRUE); utils::flush.console()
      }

      i <<- i + 1
      rr
    })

  colNames <- colnames(e)
  rm(e)
  plyr::l_ply(seq_along(pmm), function(i) { pmm[[i]] <<- as.data.frame(pmm[[i]]); names(pmm[[i]]) <<- colNames[i] })
  pmm <- purrr::reduce(pmm, dplyr::bind_cols)

  names(cutoffs) <- channels
  attr(pmm, "cutoffs") <- structure(cutoffs, starting_random_seed = starting_random_seed)

  pmm
}

## usage:
# pmm <- find_plus_minus_by_channel(ff, stringr::regex("time|event_length|sample_id", ignore_case = TRUE))


#' @export
get_fcs_expression_subset <- function(
  x, # Vector of file paths
  b = 1/150, # asinh transformation parameter: FCM = 1/150, CyTOF = 1/8 (NULL for no transformation)
  channels_subset = NULL,
  excluded_transform_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  sample_size = 10000, # 'Inf' for all events
  seed = 666,
  channels_by_sample = NULL # "cbs" object
)
{
  ## Transformation is 'asinh(a + b * x) + c'; a = shift about 0, b = scale factor, c = additive constant
  asinhTrans <- flowCore::arcsinhTransform(transformationId = "flowpipe-transformation", a = 1, b = b, c = 0)

  if (!is.null(seed))
    set.seed(seed)

  l <- keystone::psapply(x,
    function(a)
    {
      if (inherits(a, "flowFrame")) {
        ff <- a
        baseName <- a %>% flowCore::keyword() %>% `$`("FILE") %>% basename
      } else {
        ff <- flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE)
        baseName <- basename(a)
      }
      if (is.null(channels_subset))
        channels_subset <- flowCore::colnames(ff)
      ff <- ff[, channels_subset]
      tff <- ff
      if (!is.null(b)) {
        transformChannels <- flowCore::colnames(ff)[stringr::str_detect(flowCore::colnames(ff),
          excluded_transform_channels_re, negate = TRUE)]
        transList <- flowCore::transformList(transformChannels, asinhTrans)
        tff <- flowCore::transform(ff, transList)
      }
      e <- cbind(id = tools::file_path_sans_ext(baseName), flowCore::exprs(tff) %>%
        `[`(sample(NROW(.), pmin(sample_size, NROW(.))), ) %>% keystone::dataframe())

      e
    }, simplify = FALSE)

  exprs <- l %>% purrr::reduce(dplyr::bind_rows) %>%
    dplyr::mutate(id = as.factor(id)) %>%
    `attr<-`("sample_id_map",
      structure(.$id %>% keystone::unfactor() %>% unique, .Names = .$id %>% as.numeric %>% unique)) %>%
    dplyr::mutate(id = id %>% as.numeric)
  sample_id_map <- attr(exprs, "sample_id_map")
  exprs <- structure(data.matrix(exprs), sample_id_map = sample_id_map)

  ## Make a channels name-description map
  if (!is.null(channels_by_sample)) {
    attr(exprs, "channels_name_desc_map") <-
      attr(channels_by_sample, "channels_by_sample_full_desc") %>%
        (function(o) { structure(o$desc_01, .Names = o$name) })(.)
  }

  exprs
}


## Create augmented 'flowCore::flowFrame' object files to use for analysis.
#' @export
prepare_augmented_fcs_data <- function(
  x, # Vector of file paths
  b = 1/150, # asinh transformation parameter: FCM = 1/150, CyTOF = 1/8 (v. MetaCyto vignette)
  excluded_transform_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  channels_subset = NULL,
  data_dir,
  remove_outliers = TRUE, flowCut... = list(),
  outfile_prefix = expression(rep("", length(x))),
  outfile_suffix = "_pmm", # &c, also possibly 'NULL'
  overwrite = TRUE,
  ... # Additional arguments to 'find_plus_minus_by_channel()'
)
{
  if (missing(data_dir))
    stop("Must specify a write directory for augmented FCS files")

  if (!missing(data_dir) && !dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)

  if (is.null(outfile_prefix))
    outfile_prefix <- ""
  else
    outfile_prefix <- keystone::poly_eval(outfile_prefix)

  if (is.null(outfile_suffix))
    outfile_suffix <- ""
  else
    outfile_suffix <- keystone::poly_eval(outfile_suffix)

  asinhTrans <- flowCore::arcsinhTransform(transformationId = "flowpipe-transformation", a = 1, b = b, c = 0)

  ## Invoke parallel processing for no. files > 1
  sapply_fun <- sapply
  if (length(x) > 1L)
    sapply_fun <- keystone::psapply

  r <- sapply_fun(seq_along(x),
    function(i)
    {
      augmentedFcsFileName <-
        sprintf("%s%s%s.RData", outfile_prefix[i], basename(x[i]), outfile_suffix)
      augmentedFcsFilePath <- paste(data_dir, augmentedFcsFileName, sep = "/")

      if (!overwrite) {
        if (fs::file_exists(augmentedFcsFilePath)) {
          cat(sprintf("Skipping existing augmented FCS file %s.", augmentedFcsFileName), fill = TRUE); utils::flush.console()

          return (augmentedFcsFilePath)
        }
      }

      ff <- flowCore::read.FCS(x[i], transformation = FALSE, truncate_max_range = FALSE)
      if (is.null(channels_subset))
        channels_subset <- flowCore::colnames(ff)
      ff <- ff[, channels_subset]

      if (remove_outliers) {
        flowCutArgs <- list(
          f = ff,
          Directory = paste(data_dir, "flowCut", sep = "/")
        )
        flowCutArgs <- utils::modifyList(flowCutArgs, flowCut..., keep.null = TRUE)

        flowCut_results <- do.call(flowCut::flowCut, flowCutArgs)
        ff <- flowCut_results$frame
        flowCut_results <- flowCut_results[names(flowCut_results) != "frame"]
      }

      tff <- ff
      if (!is.null(b)) {
        transformChannels <- flowCore::colnames(ff)[stringr::str_detect(flowCore::colnames(ff),
          excluded_transform_channels_re, negate = TRUE)]
        transList <- flowCore::transformList(transformChannels, asinhTrans)
        tff <- flowCore::transform(ff, transList)
      }
      remove(ff)

      pmm <- find_plus_minus_by_channel(tff, ...)

      exprs_tff <- flowCore::exprs(tff)
      attr(exprs_tff, "plus_minus_matrix") <- pmm
      class(exprs_tff) <- c("pmm", class(exprs_tff))
      ## This bypasses the type checking on the 'exprs' slot; trick described here:
      ## https://stat.ethz.ch/R-manual/R-devel/library/methods/html/slot.html
      attr(tff, "exprs") <- exprs_tff

      if (exists("flowCut_results"))
        attr(tff, "flowCut_results") <- flowCut_results

      cat(sprintf("Saving augmented FCS file as %s...", augmentedFcsFileName)); utils::flush.console()
      save(tff, file = augmentedFcsFilePath)
      cat(". Done.", fill = TRUE)

      augmentedFcsFilePath
    }, simplify = TRUE)

  if (any(duplicated(r)))
    warning("Output file paths should be unique; use an 'outfile_prefix' expression to prevent overwriting", immediate. = TRUE)

  ret_val <- r

  ret_val
}


RBasicClasses <- NULL

## Assure that 'exprs' subsetting includes the plus-minus matrix.
#' @export
`[.pmm` <- function(x, ...)
{
  y <- NextMethod("[") # Dispatch to generic

  attr(y, "plus_minus_matrix") <- attr(x, "plus_minus_matrix") %>% `[`(...)
  ## Correctly subset 'cluster_id' vector attribute if present
  if (!is.null(attr(x, "cluster_id"))) {
    cluster_id <- c()
    ## The following is hacky, but it appears to always work:
    extant_colname <- head(colnames(y), 1)
    if (!is_invalid(extant_colname)) {
      cim <- structure(rep(NA, prod(dim(x), na.rm = TRUE)),
        .Dim = dim(x), .Dimnames = dimnames(x)) %>% as.data.frame
      cim[, extant_colname] <- attr(x, "cluster_id")
      cluster_id <- cim %>% `[`(...) %>% `[`(, extant_colname)
    }
    attr(y, "cluster_id") <- cluster_id
  }
  class(y) <- setdiff(class(x), RBasicClasses)
  if (is_invalid(dim(y))) {
    ## Resolves problems caused by not removing these classes from un'dim'med vectors
    class(y) <- setdiff(class(x), c("matrix", "array"))
  }

  ## Transfer unchanged attributes to 'y'
  uncommonAttributes <- setdiff(names(attributes(x)), names(attributes(y))) %>%
    `[`(. %nin% c("dim", "dimnames")) # In case of single-column subset
  plyr::l_ply(uncommonAttributes,
    function(a) { attr(y, a) <<- attributes(x)[[a]] })

  y
}


#' @export
print.pmm <- function(x, n = formals(utils:::head.default)$n, ...)
{
  head(x, n) %>% print.default(...)
}


#' @export
as.matrix.pmm <- function(x, ...)
{
  unclass(x)
}


#' @export
`colnames<-` <- function(x, value)
  UseMethod("colnames<-")


#' @export
`colnames<-.default` <- function(x, value)
{
  base::`colnames<-`(x, value)
}


#' @export
`colnames<-.pmm` <- function(x, value)
{
  x <- `colnames<-.default`(x, value)
  colnames(attr(x,"plus_minus_matrix")) <- value

  x
}


#' @export
gate <- function(
  x, # Expression matrix of class "pmm"
  strategy = list(),
  visualize_channels... = list(),
  verbose = TRUE
)
{
  if (length(strategy) == 0) {
    warning("No gating strategy supplied; returning original expression matrix.")

    return (x)
  }

  lstrategy <- keystone::psapply(strategy,
    function(a)
    {
      with(attr(x, "plus_minus_matrix") %>% tibble::as_tibble(), keystone::poly_eval(a))
    }, simplify = FALSE, .parallel = FALSE) %>%
    { Reduce(`&`, ., accumulate = TRUE) } %>%
    { c(list(rep(TRUE, NROW(x))), .) }

  gatePlots <- keystone::psapply(seq_along(strategy),
    function(a)
    {
      gatePlotPath <- tempfile()
      grDevices::png(file = gatePlotPath, bg = "transparent")
      dev.control(displaylist = "enable")

      ## Plot gating strategy
      visualize_channelsArgs <- list(
        #x = x,
        x = x[lstrategy[[a]], , drop = FALSE],
        channels = strategy[a] # Named list w/ a single element
      )
      visualize_channelsArgs <-
        utils::modifyList(visualize_channelsArgs, visualize_channels..., keep.null = TRUE)

      do.call(visualize_channels, visualize_channelsArgs)

      gatePlot <- grDevices::recordPlot()
      invisible(dev.off())
      unlink(gatePlotPath)

      gatePlot
    }, simplify = FALSE, .parallel = FALSE)

  x <- x[tail(lstrategy, 1)[[1]], , drop = FALSE]

  if (verbose && !is_invalid(names(strategy))) {
    cat("Gating strategies applied:", paste(names(strategy), collapse = ", "), fill = TRUE)
    utils::flush.console()
  }

  list(pmm = x, grobs = gatePlots)
}


#' @export
get_expression_subset <- function(
  x, # Vector of file paths to augmented 'flowCore::flowFrame' objects
  gate... = list(), # Preprocess samples individually
  save_plot_fun = grDevices::cairo_pdf, save_plot... = list(),
  seed = 666,
  sample_size = 10000,
  callback = NULL # An expression
)
{
  if (length(gate...) == 0L)
    save_plot... <- NULL

  ss <- keystone::psapply(seq_along(x),
    function(a)
    {
      ## N.B. I may want to load into an environment if object naming isn't consistent.
      ## (Then use the only object in the environment.)
      load(x[a])

      exprs_tff <- flowCore::exprs(tff)
      sampleName <- tools::file_path_sans_ext(basename(flowCore::description(tff)$FILENAME))

      remove(tff)

      grobs <- NULL
      if (length(gate...) != 0L) {
        gateArgs <- list(
          x = exprs_tff,
          visualize_channels... =
            list(
              # N.B. Change the following line to use 'rlang::enquo()' leading to a call to 'mtext()':
              #plot... = list(sub = sampleName),
              plot_end_callback = expression({
                graphics::title(main = sprintf("Gating strategy: %s", names(channels)))
                if (visualize_gates) {
                  graphics::mtext(sprintf("Events: %d/%d", NROW(xx), NROW(x)),
                    side = 1, line = -1, cex = 0.8)
                }
              })
            )
        )
        gateArgs <- utils::modifyList(gateArgs, gate..., keep.null = TRUE)

        flit <- do.call(gate, gateArgs)
        exprs_tff <- flit$pmm
        grobs <- flit$grobs
        rm(flit)
      }

      if (!is.null(seed))
        set.seed(seed)

      ## Create a new "pmm" object with a sample ID column
      e <- exprs_tff[sample(NROW(exprs_tff), pmin(sample_size, NROW(exprs_tff))), ]
      pmm <- attr(e, "plus_minus_matrix"); rownames(pmm) <- NULL

      id <- rep(a, NROW(e))
      r <- cbind(id = id, e)
      attr(r, "plus_minus_matrix") <- cbind(id = id, pmm)
      class(r) <- class(e)

      structure(r, total_gated_events = NROW(exprs_tff), grobs = grobs, sample_name = sampleName)
    }, simplify = FALSE) %>% `names<-`(x)

  ## Create pre-gating poster
  grobs <- structure(sapply(ss, function(a) attr(a, "grobs"), simplify = FALSE),
    .Names = sapply(ss, function(a) attr(a, "sample_name"), simplify = FALSE))
  anyGrobs <- !is_invalid(grobs %>% unlist %>% purrr::compact())
  if (anyGrobs) {
    max_gates <- sapply(grobs, length) %>% max
    grobs <- sapply(grobs, `length<-`, value = max_gates, simplify = FALSE)

    if (!is_invalid(save_plot...)) { # Save gating sequences to PDF
      save_plotArgs <- list(
        width = 5.0 * max_gates + 1,
        height = 5.0 * length(grobs) + 1
      )
      save_plotArgs <- utils::modifyList(save_plotArgs, save_plot..., keep.null = TRUE)

      if (!is.null(save_plotArgs$file)) {
        image_dir <- dirname(save_plotArgs$file)
        if (!dir.exists(image_dir))
          dir.create(image_dir, recursive = TRUE)
      }

      do.call(save_plot_fun, save_plotArgs)

      cowplot::plot_grid(
        ## This creates a list of "recordedplot" objects:
        #plotlist = sapply(grobs, function(a) if (is.null(a)) list(NULL) else a),
        plotlist = grobs %>% purrr::flatten(),
        ncol = max_gates,
        hjust = 0, label_x = 0.01,
        labels = rep("", max_gates * length(grobs)) %>%
          `[<-`(seq(from = 1, by = max_gates, length.out = length(grobs)), names(grobs)),
        #label_colour = "darkgreen",
        label_size = 16
      ) %>% print

      if (!is.null(save_plotArgs$file))
        dev.off()

      ## Convert PDF to PNG
      if (!is.null(save_plotArgs$file)) {
        suppressWarnings(pdftools::pdf_convert(
          pdf = save_plotArgs$file,
          format = "png",
          dpi = 100,
          filenames = sprintf("%s.png", tools::file_path_sans_ext(save_plotArgs$file))
        ))

        ## Resize very large PNGs to fit into LaTeX documents
        local({
          filePath <- sprintf("%s.png", tools::file_path_sans_ext(save_plotArgs$file))
          i <- magick::image_read(filePath)
          ii <- magick::image_info(i)[, c("width", "height")]
          if (any(ii %>% unlist > 16384)) {
            ni <- magick::image_scale(i, ifelse(ii$width > 16384, "", "x") %_% "16384")
            magick::image_write(ni, sprintf("%s-resized.%s", tools::file_path_sans_ext(filePath), tools::file_ext(filePath)))
          }
        })
      }
    }
  }

  if (!is.null(callback))
    keystone::poly_eval(callback) # E.g. standardize column names

  ## Assemble a new "pmm" object from all subsets.
  e <- purrr::reduce(ss, rbind)
  pmm <- purrr::reduce(sapply(ss, function(a) attr(a, "plus_minus_matrix"), simplify = FALSE), rbind); rownames(pmm) <- NULL
  attr(e, "plus_minus_matrix") <- pmm
  attr(e, "id_map") <- structure(pmm$id %>% unique, .Names = x[pmm$id %>% unique])
  attr(e, "total_gated_events") <- sapply(ss, function(a) attr(a, "total_gated_events"))
  class(e) <- class(ss[[1]])

  e
}


#' @export
rename_duplicates <- function(
  x, # Character vector
  fmt = "%s-%02d"
)
{
  if (any(duplicated(x))) {
    duplicateNames <- x %>% intersect(.[duplicated(.)])
    for(i in duplicateNames) {
      dupIndex <- which(x == i)
      ## Replace w/ sequential numbers:
      x[dupIndex] <- sapply(seq_along(dupIndex), function(j) sprintf(fmt, x[dupIndex[j]], j))
    }
  }

  x
}


#' @export
make_sample_id_map <- function(
  x, # Expression matrix of class "pmm"
  re, # Regular expression or vector of IDs the same length as 'attr(x, "id_map")'
  rename_dups = TRUE,
  ...
)
{
  sample_id_map <- attr(x, "id_map")
  if (length(re) == length(sample_id_map))
    mapNames <- re
  else
    mapNames <- names(sample_id_map) %>% basename %>% stringr::str_extract(re[1])

  if (rename_dups)
    mapNames <- rename_duplicates(mapNames, ...)

  names(sample_id_map) <- mapNames

  sample_id_map
}
