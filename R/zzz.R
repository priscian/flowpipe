.onLoad <- function(...)
{
  if (!exists(".cm", envir = globalenv()))
    .cm <- memoise::cache_memory()

  bead_normalize <<- memoise::memoise(bead_normalize, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  get_channels_by_sample <<- memoise::memoise(get_channels_by_sample, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  rename_fcs_parameters_name_desc <<- memoise::memoise(rename_fcs_parameters_name_desc, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  plot_channel_densities_by_sample <<- memoise::memoise(plot_channel_densities_by_sample, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  prepare_augmented_fcs_data <<- memoise::memoise(prepare_augmented_fcs_data, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  get_expression_subset <<- memoise::memoise(get_expression_subset, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  plot_cell_counts <<- memoise::memoise(plot_cell_counts, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  make_clusters <<- memoise::memoise(make_clusters, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  make_metaclusters <<- memoise::memoise(make_metaclusters, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  make_umap_embedding <<- memoise::memoise(make_umap_embedding, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  merge_clusters <<- memoise::memoise(merge_clusters, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  summarize_all_clusters <<- memoise::memoise(summarize_all_clusters, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  plot_common_umap_viz <<- memoise::memoise(plot_common_umap_viz, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  plot_heatmaps <<- memoise::memoise(plot_heatmaps, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  plot_differential_abundance <<- memoise::memoise(plot_differential_abundance, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  export_cluster_summary <<- memoise::memoise(export_cluster_summary, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()

  split_pmm_by_cluster <<- memoise::memoise(split_pmm_by_cluster, cache = .cm) %>%
    keystone::patch_memoised_for_subcaching()
}


.onAttach <- function(...)
{
  keystone:::.onAttach()
}


.onDetach <- function(...)
{
  keystone:::.onDetach()
}
