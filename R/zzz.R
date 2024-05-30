.onLoad <- function(...)
{
  ## For use in '`[.pmm`()'
  RBasicClasses <<- keystone::R_basic_classes() %>%
    sapply(class) %>% as.vector %>% unique #%>% c("matrix", "array")
}


memoise_package_function <- keystone::memoise_package_function_ %>%
  `formals<-`(value = formals(.) %>%
    `$<-`("env", "flowpipe")) %>%
  `environment<-`(asNamespace("flowpipe"))


#' @export
memoise_all_package_functions <- function(key, ...)
{
  mpf <- memoise_package_function %>%
    `formals<-`(value = formals(.) %>%
      `$<-`("cache", .cm))

  if (!missing(key)) {
    `_mpf` <- mpf
    mpf <- function(...) { `_mpf`(..., key = key) }
  }

  if (!exists(".cm", envir = globalenv())) {
    .cm <- memoise::cache_memory()
  }

  mpf("bead_normalize")
  mpf("get_channels_by_sample")
  mpf("rename_fcs_parameters_name_desc")
  #mpf("guess_channels")
  mpf("plot_channel_densities_by_sample")
  mpf("prepare_augmented_fcs_data")
  mpf("get_expression_subset")
  mpf("plot_cell_counts")
  mpf("make_clusters")
  mpf("make_metaclusters")
  mpf("merge_clusters")
  mpf("make_umap_embedding")
  mpf("summarize_all_clusters")
  mpf("plot_common_umap_viz")
  mpf("plot_heatmaps")
  mpf("plot_differential_abundance")
  mpf("export_cluster_summary")
  mpf("split_pmm_by_cluster")
}


.onAttach <- function(...)
{
  keystone:::.onAttach()

  ## Unlock '.[package]' variable to allow its modification:
  unlockBinding(".flowpipe", asNamespace("flowpipe"))

  ## Startup message
  msg <- flowpipeStartupMessage()
  if (!interactive())
    msg[1] <- paste("Package 'flowpipe' version", packageVersion("flowpipe"))

  packageStartupMessage(msg)

  invisible()
}


.onDetach <- function(...)
{
  keystone:::.onDetach()
}


flowpipeStartupMessage <- function()
{
  msg01 <- c(paste0(
r"---{                                __________
    ______                    .' ________ '.
   / __/ /___ _      ______  / /___  ___ `\ \
  / /_/ / __ \ | /| / / __ \/ / __ \/ _ \ [__]
 / __/ / /_/ / |/ |/ / /_/ / / /_/ /  __/  {{
/_/ /_/\____/|__/|__/ .___/_/ .___/\___/   }}
                   /_/     /_/             {{}---",
sprintf("\nVersion %-12s                       }}\n", as.character(packageVersion("flowpipe"))),
r"---{https://CART.urmc.edu                      {{
                                           }}
`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``'-.,

Type 'citation("flowpipe")' to acknowledge this R package in publications.

Use 'process_config_file()' to get started.
}---")
  )

  ## Alternative ASCII art:
  msg02 <- c(paste0(
r"---{    ______
   / __/ /___ _      ______  __ ___  ___
  / /_/ / __ \ | /| / / __ \/ / __ \/ _ \
 / __/ / /_/ / |/ |/ / /_/ / / /_/ /  __/____
/_/ /_/\____/|__/|__/ .___/_/ .___/\________ '.
                   /_/     /_/              `\ \}---",
sprintf("\nVersion %-12s                         [__]\n", as.character(packageVersion("flowpipe"))),
r"---{https://CART.urmc.edu                         }}
                                              {{
='``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``'-.,

Type 'citation("flowpipe")' to acknowledge this R package in publications.

Use 'process_config_file()' to get started.
}---")
  )

  return (msg01)
}
