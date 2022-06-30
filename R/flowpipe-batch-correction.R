#' @export
batch_correct <- function(
  x, # Vector of file paths
  channels,
  batch = NULL, # Optional factor identifying batches
  output_dir = ".", output_suffix = "_batch-corr",
  batch_corr_fun = flowStats::gaussNorm, # https://rdrr.io/bioc/flowStats/man/gaussNorm.html
  batch_corr_fun... = flowpipe:::gaussNorm...,
  results_expr = flowpipe:::gaussNorm_results_expr,
  timer = TRUE,
  ...
)
{
  ## Other possible batch-correction functions (future expansion):
  ##   https://github.com/CUHIMSR/CytofBatchAdjust
  ##   https://rdrr.io/bioc/cydar/man/normalizeBatch.html

  if (timer) tictoc::tic("Batch correction")

  r <- do.call(batch_corr_fun, keystone::poly_eval(batch_corr_fun...))

  keystone::poly_eval(results_expr)

  if (timer) tictoc::toc()

  ret_val
}

### Arguments for 'flowStats::gaussNorm()'

gaussNorm... <- expression({
  list(
    flowset = flowCore::read.flowSet(x, transformation = FALSE, truncate_max_range = FALSE),
    channel.names = channels,
    max.lms = 1,
    debug = TRUE,
    fname = paste(output_dir, "batch-correction-sample", sep = "/")
  )
})

gaussNorm_results_expr <- expression({
  outFiles <- sprintf("%s%s.fcs", tools::file_path_sans_ext(basename(x)), outfile_suffix)
  flowCore::write.flowSet(
    x = r$flowset,
    outdir = output_dir,
    files = outFiles
  )
  batch_confidence = r$confidence
  ret_val <- paste(output_dir, outFiles, sep = "/")
})

### Arguments for 'cydar::normalizeBatch()'
normalizeBatch... <- expression({
  list(
    batch.x = sapply(x,
      function(i) { flowCore::read.FCS(i, transformation = FALSE, truncate_max_range = FALSE) %>% flowCore::exprs() },
      simplify = FALSE) %>% split(f = batch),
    batch.comp = NULL,
    mode = "warp",
    markers = channels
  )
})

normalizeBatch_results_expr <- expression({
  ## TBD, doesn't work anyway
  ## Probably need to turn 'exprs' matrices back into flowFrames & save as FCS
})

### Arguments for 'CytoNorm::CytoNorm.normalize()'

run_CytoNorm.normalize <- function(
  metadata,
  channels,
  transform, transform_reverse,
  FlowSOM.params... = list(),
  prepareFlowSOM... = list(),
  CytoNorm.train... = list(),
  CytoNorm.normalize... = list(),
  nCells = 1e6,
  quantile_prob = 0.99,
  seed = 666
)
{
  train_data <- dplyr::filter(metadata, Type == "Train")
  validation_data <- dplyr::filter(metadata, Type == "Validation")

  transformList <- flowCore::transformList(channels, transform)
  transformList.reverse <- flowCore::transformList(channels, transform_reverse)

  FlowSOM.params <- list(
    xdim = 15,
    ydim = 15,
    nClus = 10,
    scale = FALSE
  )
  FlowSOM.params <- utils::modifyList(FlowSOM.params, FlowSOM.params..., keep.null = TRUE)

  prepareFlowSOMArgs <- list(
    files = train_data$Path,
    colsToUse = channels,
    nCells = nCells,
    FlowSOM.params = FlowSOM.params,
    transformList = transformList,
    seed = seed
  )
  prepareFlowSOMArgs <- utils::modifyList(prepareFlowSOMArgs, prepareFlowSOM..., keep.null = TRUE)

  fsom <- do.call(CytoNorm::prepareFlowSOM, prepareFlowSOMArgs)

  cvs <- CytoNorm::testCV(fsom, cluster_values = seq(5, 50, by = 5))

  CytoNorm.trainArgs <- list(
    files = train_data$Path,
    labels = train_data$Batch,
    channels = channels,
    transformList = transformList,
    FlowSOM.params = FlowSOM.params %>% `$<-`("nCells", nCells),
    normMethod.train = CytoNorm::QuantileNorm.train,
    normParams = list(nQ = 101, goal = "mean"),
    seed = seed,
    verbose = TRUE
  )
  CytoNorm.trainArgs <- utils::modifyList(CytoNorm.trainArgs, CytoNorm.train..., keep.null = TRUE)

  model <- do.call(CytoNorm::CytoNorm.train, CytoNorm.trainArgs)

  CytoNorm.normalizeArgs <- list(
    model = model,
    files = validation_data$Path,
    labels = validation_data$Batch,
    transformList = transformList,
    transformList.reverse = transformList.reverse,
    normMethod.normalize = CytoNorm::QuantileNorm.normalize,
    clean = TRUE,
    verbose = TRUE
  )

  CytoNorm.normalizeArgs <- utils::modifyList(CytoNorm.normalizeArgs, CytoNorm.normalize..., keep.null = TRUE)

  do.call(CytoNorm::CytoNorm.normalize, CytoNorm.normalizeArgs)

  normalized_files <-
    keystone::list_files(CytoNorm.normalizeArgs$outputDir, "\\.fcs$",
      recursive = 1, ignore.case = TRUE, full.names = TRUE)

  ## This method produces wacky outlier values; fix them in situ!
  # plyr::l_ply(normalized_files,
  #   function(a)
  #   {
  #     ff <- flowCore::read.FCS(a, transformation = FALSE, truncate_max_range = FALSE)
  #     p0 <- flowCore::pData(flowCore::parameters(ff))
  #     e <- flowCore::exprs(ff) %>% tibble::as_tibble() %>% dplyr::mutate(
  #       dplyr::across(.cols = all_of(channels %>% as.vector),
  #         .fns = function(i) { q <- quantile(i, probs = quantile_prob); i[i > q] <- q; i[i < 0] <- 0; i })
  #     )
  #     ff_new <- flowCore::flowFrame(exprs = e %>% data.matrix)
  #     p <- flowCore::pData(flowCore::parameters(ff_new))
  #     p$desc <- p0$desc %>% as.vector
  #     flowCore::parameters(ff_new) <- as(p, "AnnotatedDataFrame")

  #     flowCore::write.FCS(ff_new, a)
  #   })

  ret_val <- list(
    fcs_files = c(metadata %>%
      dplyr::filter(Type == "Train") %>% dplyr::pull(Path), normalized_files) %>%
      sort,
    cvs = cvs
  )

  ret_val
}


run_CytoNorm.normalize... <- expression({
  list(
    metadata = batch,
    channels = channels,
    transform = flowpipe::cytof_transform, # Change for FCM
    transform_reverse = flowpipe::cytof_transform_reverse, # Change for FCM
    CytoNorm.normalize... = list(
      outputDir = output_dir,
      prefix = output_suffix
    )
  )
})

run_CytoNorm.normalize_results_expr <- expression({
  fcs_files <- r$fcs_files
  cvs <- r$cvs

  ret_val <- r
})

## À la CytoNorm
#' @export
fcm_transform <- flowCore::arcsinhTransform(transformationId = "fcmTransform", a = 0, b = 1/150, c = 0)

#' @export
fcm_transform_reverse <- function(x) { return (sinh(x) / (1/150)) }

#' @export
cytof_transform <- flowCore::arcsinhTransform(transformationId = "cytofTransform", a = 0, b = 1/5, c = 0)

#' @export
cytof_transform_reverse <- function(x) { return (sinh(x) / (1/5)) }
