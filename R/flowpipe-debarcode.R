#' @export
debarcode <- function(
  x, # flowCore::flowFrame
  barcoding_key
)
{
  if (is.null(barcoding_key))
    return (rep(basename(flowCore::description(x)$FILENAME), NROW(x)))

  plinth::poly_eval(barcoding_key$process_key)

  ## This shouldn't ever need tweaking, but we'll see:
  row_data <- flowCore::pData(flowCore::parameters(x)) %>%
    dplyr::select(name, desc) %>%
    dplyr::rename(channel_name = "name", marker_name = "desc") %>%
    dplyr::mutate(
      dplyr::across(.fns = as.character),
      marker_name = dplyr::case_when(is.na(marker_name) ~ channel_name, TRUE ~ marker_name)
    ) %>%
    textshape::column_to_rownames(loc = "marker_name")

  e <- flowCore::exprs(x)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = t(e) %>% `rownames<-`(as.vector(rownames(.)))),
    rowData = row_data)

  sce <- CATALYST::assignPrelim(sce, bc_key = bc_key, assay = "exprs", verbose = TRUE)
  sce <- CATALYST::estCutoffs(sce)
  sce <- CATALYST::applyCutoffs(sce)
  #table(sce$bc_id)

  sce$bc_id
}
