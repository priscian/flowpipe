graphics_devices <- list(
  `grDevices::png` = list(height = 8.27, width = 11.69, units = "in", res = 600, type = "cairo", ext = "png"),
  `grDevices::pdf` = list(height = 8.27, width = 11.69, ext = "pdf")
)

## Density plots for all samples, all channels
#' @export
plot_channel_densities_by_sample <- function(
  x, # Vector of file paths
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  channel_dist_palette = randomcoloR::distinctColorPalette,
  devices = flowpipe:::graphics_devices,
  file_path_template = sprintf("%s/%03d%s", image_dir, current_image, "_channel_dist"),
  get_fcs_expression_subset... = list()
)
{
  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  get_fcs_expression_subsetArgs <- list(
    x = x,
    sample_size = 10000
  )
  get_fcs_expression_subsetArgs <- utils::modifyList(get_fcs_expression_subsetArgs, get_fcs_expression_subset..., keep.null = TRUE)

  e <- do.call(get_fcs_expression_subset, get_fcs_expression_subsetArgs)

  ## Identify channels by descriptions
  colnames(e) <- c(id = "id", attr(e, "channels_name_desc_map"))[colnames(e)] %>% as.vector

  num_col <- NCOL(e) - 1
  num_row <- e[, "id"] %>% as.vector %>% unique %>% length
  channel_colors <- channel_dist_palette(num_col)
  ## This one still works, but is slower:
  #ffn <- basename(x[e[, "id"] %>% as.vector]) %>% tools::file_path_sans_ext() %>% tools::file_path_sans_ext()
  ffn <- attr(e, "sample_id_map")[e[, "id"]] %>% as.vector
  density_plots <- sapply(seq(num_col), # First column is "id"
    function(i)
    {
      ggplot2::ggplot(as.data.frame(e) %>% dplyr::mutate(id = as.factor(ffn)),
        ggplot2::aes_string(x = plinth::backtick(colnames(e)[i + 1]), y = "id")) +
        ggridges::geom_density_ridges(color = channel_colors[i], fill = scales::alpha(channel_colors[i], 0.4)) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::theme_bw()
    }, simplify = FALSE)

  g <- cowplot::plot_grid(plotlist = density_plots, ncol = num_col)

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 50
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
        list(
          width = min(6.0 * num_col + 1, 200), # 200 in. is PDF maximum
          height = min(4.0 * num_row + 1, 200), # 200 in. is PDF maximum
          file = plinth::poly_eval(file_path_template) %_% paste0(".", ext)
        ))
      )

      print(g)

      dev.off()
    })
  } else {
    print(g)
  }
  par(def_par)

  plinth::nop()
}


#' @export
visualize_channels <- function(
  x, # Expression matrix of class "pmm"
  channels, # Character vector or gating expression a là 'gate()' as a 1-element named list
  extract_gating_channels = NULL, # Replace only if the default function doesn't work well
  plot... = list(),
  kde2d... = list(),
  contour... = list(),
  points... = list(),
  abline... = list(),
  plot_end_callback = NULL
)
{
  if (is.null(extract_gating_channels))
    extract_gating_channels <- function(m, s)
    {
      re <- stringr::regex(stringr::str_flatten(rex::escape(colnames(m)), "|"))
      stringr::str_match_all(as.character(s[[1]]), re)[[1]] %>% drop %>% unique
    }

  visualize_gates <- FALSE
  if (!is.character(channels)) {
    visualize_gates <- TRUE
    gating_channels <- extract_gating_channels(x, channels)
    cutoffs <- (x %>% attr("plus_minus_matrix") %>% attr("cutoffs"))[gating_channels]
  } else {
    gating_channels <- channels
  }

  if (is_invalid(x)) {
    plot.new()

    return (plinth::nop())
  }

  if (length(gating_channels) == 1) { # Density plot
    plot.densityArgs <- list(
      x = stats::density(x[, gating_channels[1]]),
      main = "",
      xlim = c(min(0, x[, gating_channels[1]]), max(x[, gating_channels[1]])),
      ylab = paste(gating_channels[1], "Density")
    )
    plot.densityArgs <- utils::modifyList(plot.densityArgs, plot..., keep.null = TRUE)

    do.call(plot, plot.densityArgs)

    if (visualize_gates) {
      plinth::vline(
        sprintf("%.2f", cutoffs[[gating_channels[1]]]),
        abline... = list(col = scales::alpha("red", 0.5), lty = "dashed"),
        text... = list(y = plinth::cp_coords()$y)
      )

      pmm <- attr(x, "plus_minus_matrix") %>% tibble::as_tibble()
      xx <- x[with(pmm, plinth::poly_eval(channels[[1]])), , drop = FALSE]

      p <- with(plot.densityArgs$x, plinth::dataframe(x = x, y = y)) %>%
        dplyr::filter(x >= min(xx[, gating_channels[1]], na.rm = TRUE) &
          x <= max(xx[, gating_channels[1]]), na.rm = TRUE) %>%
        dplyr::mutate(
          y = dplyr::case_when(x == min(x, na.rm = TRUE) | x == max(x, na.rm = TRUE) ~ 0.0, TRUE ~ y)
        ) %>% data.matrix

      polygon(p, col = scales::alpha("red", 0.5), border = NA)
    }
  } else { # Biaxial plot
    kde2dArgs <- list(
      x = x[, gating_channels[1]],
      y = x[, gating_channels[2]],
      n = 50,
      h = max(
        MASS::bandwidth.nrd(x[, gating_channels[1]]),
        MASS::bandwidth.nrd(x[, gating_channels[2]])
      )
    )
    kde2dArgs <- utils::modifyList(kde2dArgs, kde2d..., keep.null = TRUE)

    z <- do.call(MASS::kde2d, kde2dArgs)

    plotArgs <- list(
      x = x[, gating_channels],
      pch = ".",
      cex = 0.1,
      col = "gray",
      xlim = c(min(0, x[, gating_channels[1]]), max(x[, gating_channels[1]])),
      ylim = c(min(0, x[, gating_channels[2]]), max(x[, gating_channels[2]])),
      xlab = gating_channels[1], ylab = gating_channels[2]
    )
    plotArgs <- utils::modifyList(plotArgs, plot..., keep.null = TRUE)

    do.call(plot, plotArgs)

    contourArgs <- list(
      x = z,
      drawlabels = FALSE,
      nlevels = 11,
      add = TRUE
    )
    if (is.null(contourArgs$col))
      contourArgs$col <- rev(RColorBrewer::brewer.pal(contourArgs$nlevels, "RdYlBu"))
    contourArgs <- utils::modifyList(contourArgs, contour..., keep.null = TRUE)

    do.call(graphics::contour, contourArgs)

    if (visualize_gates) {
      pmm <- attr(x, "plus_minus_matrix") %>% tibble::as_tibble()
      xx <- x[with(pmm, plinth::poly_eval(channels[[1]])), , drop = FALSE]

      pointsArgs <- list(
        x = xx[, gating_channels],
        pch = plotArgs$pch,
        cex = 0.5,
        col = "red"
      )
      pointsArgs <- utils::modifyList(pointsArgs, points..., keep.null = TRUE)

      do.call(points, pointsArgs)

      ablineArgs <- list(
        col = scales::alpha("red", 0.5),
        lty = "dashed"
      )
      ablineArgs <- utils::modifyList(ablineArgs, abline..., keep.null = TRUE)

      do.call(graphics::abline, c(ablineArgs, list(v = cutoffs[[gating_channels[1]]])))
      do.call(graphics::abline, c(ablineArgs, list(h = cutoffs[[gating_channels[2]]])))
    }
  }

  plinth::poly_eval(plot_end_callback)

  plinth::nop()
}


#' @export
make_tsne_embedding <- function(
  x, # "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names otherwise
  seed = 666,
  ...
)
{
  set.seed(seed)
  tsne <- Rtsne::Rtsne(x[, channels, drop = FALSE], ...)

  tsne
}

## usage:
# cordon(make_tsne_embedding, perplexity = 30, envir = globalenv(), file_path = paste(data_dir, "flowpipe-test-tsne.RData", sep = "/"), variables = c("tsne"), timestamp... = list(use_seconds = TRUE), action = "save")


#' @export
make_umap_embedding <- function(
  x, # "pmm" object from 'get_expression_subset()'
  channels = TRUE, # Vector of column names otherwise
  n_neighbors = 15, min_dist = 0.2, metric = "euclidean",
  seed = 666,
  ...
)
{
  set.seed(seed)
  umap <- uwot::umap(
    x[, channels],
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    ...)

  umap
}

## usage:
# cordon(make_umap_embedding, n_neighbors = 15, envir = globalenv(), file_path = paste(data_dir, "flowpipe-test-umap.RData", sep = "/"), variables = c("umap"), timestamp... = list(use_seconds = TRUE), action = "run") # 'run' or 'save'


plot_common_umap_viz_single <- function(
  x, # "pmm" object from 'get_expression_subset()'
  which_cluster_set = 1, # Column no. or name
  channels = TRUE,
  m = NULL, # metadata; if NULL, plot will be skipped
  umap,
  sample_name_re = "^.*$",
  plot_palette = randomcoloR::distinctColorPalette,
  na.value = scales::alpha("grey95", 0.3), # Default is "grey50"; NA for transparent
  labels = ~ (function(x) { ifelse(is.na(x), "other", x) })(.x),
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  cluster_plot_only = FALSE,
  devices = flowpipe:::graphics_devices,
  file_path_template =
    sprintf("%s/%s/%03d%s-%s", image_dir,
      fs::path_sanitize(as.character(which_cluster_set), "_"),
      current_image, "_umap", fs::path_sanitize(as.character(which_cluster_set), "_")),
  seed = 666
)
{
  if (!is.null(seed))
    set.seed(seed)

  if (is.logical(channels) && channels)
    channels <- colnames(x)

  file_path_template <- plinth::poly_eval(file_path_template)

  if (save_plot && !dir.exists(dirname(file_path_template)))
    dir.create(dirname(file_path_template), recursive = TRUE)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ## Make plotting data set
  cluster_id <- attr(x, "cluster_id")
  if (is.matrix(cluster_id)) {
    cluster_id <- cluster_id[, which_cluster_set]
  }
  # sample_id_map <- structure(attr(x, "id_map"),
  #   .Names = names(attr(x, "id_map")) %>% basename %>% stringr::str_extract(sample_name_re))
  sample_id_map <- make_sample_id_map(x, sample_name_re)
  sample_id <- names(sample_id_map)[x[, "id"]]
  d <- plinth::dataframe(UMAP1 = umap[, 1], UMAP2 = umap[, 2], x[, channels],
    cluster_id = cluster_id %>% as.factor,
    sample_id = sample_id %>% as.factor
  )
  if (!is.null(m)) {
    id_group_map <- structure(m$group, .Names = m$id)
    group_id <- id_group_map[sample_id]
    d$group_id <- group_id
  } else {
    d$group_id <- "all"
    group_id <- d$group_id
  }

  ##### By cluster #####

  g3 <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, color = cluster_id)) +
    ggplot2::geom_point(size = 0.8) + # default: 'shape = 16'
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = plot_palette(length(unique(cluster_id))), na.value = na.value, labels = labels) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 3))

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 20, height = 15,
            file = file_path_template %_% "-clusters" %_% paste0(".", ext)
          )
        )
      )

      print(g3)

      dev.off()
    })
  } else {
    print(g3)
  }
  par(def_par)

  ## Facet per sample
  g3a <- g3 + ggplot2::facet_wrap(~ sample_id)

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 20, height = 15,
            file = file_path_template %_% "-clusters-by-sample" %_% paste0(".", ext)
          )
        )
      )

      print(g3a)

      dev.off()
    })
  } else {
    print(g3a)
  }
  par(def_par)

  ## Facet per group
  g3b <- g3 + ggplot2::facet_wrap(~ group_id)

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 20, height = 15,
            file = file_path_template %_% "-clusters-by-group" %_% paste0(".", ext)
          )
        )
      )

      print(g3b)

      dev.off()
    })
  } else {
    print(g3b)
  }
  par(def_par)

  if (cluster_plot_only)
    return (plinth::nop())

  ##### By sample #####

  g1 <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, color = sample_id)) +
    ggplot2::geom_point(size = 0.8) + # default: 'shape = 16'
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = plot_palette(length(unique(sample_id))), na.value = na.value) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 3))

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 20, height = 15,
            file = file_path_template %_% "-samples" %_% paste0(".", ext)
          )
        )
      )

      print(g1)

      dev.off()
    })
  } else {
    print(g1)
  }
  par(def_par)

  ##### By group #####

  g2 <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, color = group_id)) +
    ggplot2::geom_point(size = 0.8) + # default: 'shape = 16'
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = plot_palette(length(unique(group_id))), na.value = na.value) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 3))

  if (!is.null(m)) # Skip this plot if metadata arg is ignored
  {
    if (save_plot) {
      plyr::l_ply(names(devices),
      function(a)
      {
        ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
        ## Reduce resolution for 'png()' etc. to a manageable value:
        if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
        do.call(eval(parse(text = a)),
          modifyList(devices[[a]],
            list(
              width = 20, height = 15,
              file = file_path_template %_% "-groups" %_% paste0(".", ext)
            )
          )
        )

        print(g2)

        dev.off()
      })
    } else {
      print(g2)
    }
    par(def_par)
  }

  #browser()
  plinth::nop()
}


#' @export
plot_common_umap_viz <- function(
  x,
  which_cluster_set = NULL,
  ...
)
{
  clusterId <- attr(x, "cluster_id")

  if (is.null(which_cluster_set)) {
    which_cluster_set <- 1
    if (is.matrix(clusterId))
      which_cluster_set <- colnames(clusterId)
  }

  clusterPlotOnly <- FALSE
  plyr::l_ply(which_cluster_set,
    function(a)
    {
      plot_common_umap_viz_single(x, which_cluster_set = a, cluster_plot_only = clusterPlotOnly, ...)

      ## Only make group & sample plots once
      if (!clusterPlotOnly)
        clusterPlotOnly <<- TRUE

      plinth::nop()
    })
}


#' @export
plot_cell_counts <- function(
  x, # "pmm" object from 'get_expression_subset()'
  pmm_files, # Vector of file paths of "pmm" objects
  m, # metadata
  sample_name_re = "^.*$",
  plot_palette = randomcoloR::distinctColorPalette,
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  devices = flowpipe:::graphics_devices,
  file_path_template = sprintf("%s/%03d%s", image_dir, current_image, "_counts"),
  seed = 666
)
{
  if (!is.null(seed))
    set.seed(seed)

  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  file_path_template <- plinth::poly_eval(file_path_template)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ##### Cell counts #####

  cellCounts <- sapply(pmm_files,
    function(a)
    {
      e <- new.env()
      load(a, envir = e)

      flowCore::exprs(e$tff) %>% NROW
    }, simplify = TRUE)
  cell_counts <- cellCounts
  names(cell_counts) <- tools::file_path_sans_ext(basename(names(cell_counts))) %>%
    stringr::str_extract(sample_name_re) %>% rename_duplicates

  ggdf <- plinth::dataframe(sample_id = names(cell_counts), cell_counts = as.numeric(cell_counts))
  ## Add metadata info to 'ggdf'
  ggdf <- dplyr::left_join(ggdf, m %>% dplyr::select(id, group), by = c(sample_id = "id")) %>%
    dplyr::mutate(group = as.factor(group))
  tge <- structure(
    attr(x, "total_gated_events"),
    .Names = tools::file_path_sans_ext(basename(names(attr(x, "total_gated_events")))) %>%
      stringr::str_extract(sample_name_re) %>% rename_duplicates
  )
  ggdf <- dplyr::left_join(
    ggdf,
    structure(list(sample_id = names(tge), gated_cell_counts = as.vector(tge)), row.names = c(NA, length(tge)), class = "data.frame"),
    by = "sample_id"
  )

  ## Define colors for groups
  color_conditions <- plot_palette(nlevels(ggdf$group))
  names(color_conditions) <- levels(ggdf$group)

  plot_event_counts <- function(
    d,
    y
  )
  {
    g <- ggplot2::ggplot(d, ggplot2::aes_string(x = "sample_id", y = y, fill = "group")) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes_string(label = y), hjust = 0.5, vjust = -0.5, size = 2.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::scale_fill_manual(values = color_conditions, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE)

    g
  }

  g2 <- plot_event_counts(ggdf, "cell_counts") +
    ggplot2::labs(title = "Starting event counts")

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      #if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 100
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
        list(
          file = file_path_template %_% paste0(".", ext)
        ))
      )

      print(g2)

      dev.off()
    })
  } else {
    print(g2)
  }
  par(def_par)

  ##### Cell counts after initial gating #####

  g4 <- plot_event_counts(ggdf, "gated_cell_counts") +
    ggplot2::labs(title = "Event counts after gating")

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      #if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 100
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
        list(
          file = file_path_template %_% "_gated" %_% paste0(".", ext)
        ))
      )

      print(g4)

      dev.off()
    })
  } else {
    print(g4)
  }
  par(def_par)

  plinth::nop()
}


plot_heatmaps_single <- function(
  x, # "pmm" object from 'get_expression_subset()'
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  channels = TRUE,
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  devices = flowpipe:::graphics_devices,
  file_path_template =
    sprintf("%s/%03d%s-%s", image_dir, current_image, "_heatmap_channels-clusters", fs::path_sanitize(as.character(which_cluster_set), "_"))
)
{
  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  file_path_template <- plinth::poly_eval(file_path_template)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ##### Heatmap of channels vs. clusters #####

  if (is.logical(channels) && channels)
    channels <- colnames(x)

  ## 'cluster_matrix' contains the median expression for each channel for each cluster. Each row (each channel) is normalized by
  ## 'base::scale()' to allow relative comparison of normalized median expression among all clusters for each channel. This
  ## mostly disregards variance of the marker expression in each cluster, but offers an overview of each cluster's composition.
  sub_matrix <- x[, channels]

  cluster_id <- attr(x, "cluster_id")
  if (is.matrix(cluster_id))
    cluster_id <- cluster_id[, which_cluster_set] %>% drop

  cluster_matrix <- NULL
  for (i in sort(unique(cluster_id))) {
    cluster_matrix <- rbind(cluster_matrix, matrixStats::colMedians(sub_matrix[cluster_id == i, , drop = FALSE], na.rm = TRUE))
  }
  colnames(cluster_matrix) <- channels
  rownames(cluster_matrix) <- sort(unique(cluster_id))
  #par(mar = c(2, 2, 2, 2)) # Prob. doesn't do anything here

  if (any(dim(cluster_matrix) < 2))
    return (cluster_matrix)

  nr <- dim(cluster_matrix)[2]
  nc <- dim(cluster_matrix)[1]
  cexRow <- min(0.2 + 1/log10(nr), 1.0)
  cexCol <- min(0.2 + 1/log10(nc), 1.0)

  g5_expr <- expression({
    gplots::heatmap.2(
      t(cluster_matrix),
      col = gplots::bluered(100),
      trace = "none", density.info = "none",
      sepcolor = "white", sepwidth = c(0.001, 0.001),
      colsep = c(1:NCOL(t(cluster_matrix))),
      rowsep = c(1:NROW(t(cluster_matrix))),
      xlab = "cluster", ylab = "channel", scale = scale_,
      margins = c(15, 10), # Increase these to give more room to col & row labels, respectively
      cexRow = cexRow, cexCol = cexCol
    )
  })

  ## Center & scale heatmap values in the row (i.e. channel) direction
  scale_ <- "row"

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = file_path_template %_% "-row" %_% paste0(".", ext)
          )
        )
      )

      plinth::poly_eval(g5_expr)

      dev.off()
    })
  } else {
    plinth::poly_eval(g5_expr)
  }
  par(def_par)

  ## Center & scale heatmap values in the column (i.e. cluster) direction
  scale_ <- "column"

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = file_path_template %_% "-column" %_% paste0(".", ext)
          )
        )
      )

      plinth::poly_eval(g5_expr)

      dev.off()
    })
  } else {
    plinth::poly_eval(g5_expr)
  }
  par(def_par)

  cluster_matrix
}


#' @export
plot_heatmaps <- function(
  x,
  which_cluster_set = NULL,
  ...
)
{
  clusterId <- attr(x, "cluster_id")

  if (is.null(which_cluster_set)) {
    which_cluster_set <- 1
    if (is.matrix(clusterId))
      which_cluster_set <- colnames(clusterId)
  }

  cluster_matrices <- sapply(which_cluster_set,
    function(a)
    {
      plot_heatmaps_single(x, which_cluster_set = a, ...)
    }, simplify = FALSE)

  cluster_matrices
}


plot_differential_abundance_single <- function(
  x, # "pmm" object from 'get_expression_subset()'
  m, # metadata
  umap,
  fit, contrasts,
  which_cluster_set = 1, # If 'attr(x, "cluster_id")' is matrix, pick a column by name or number
  alpha = 0.05,
  sample_name_re = "^.*$",
  results_column = "logFC",
  plot_palette = randomcoloR::distinctColorPalette,
  na.value = scales::alpha("grey95", 0.3), # Default is "grey50"; NA for transparent
  labels = ~ (function(x) { ifelse(is.na(x), "other", x) })(.x),
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  devices = flowpipe:::graphics_devices,
  file_path_template =
    sprintf("%s/%03d%s-%s", image_dir, current_image, "_diff-abundance_significant-clusters", fs::path_sanitize(as.character(which_cluster_set), "_"))
)
{
  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ##### UMAP visualization of significant differential abundance between conditions #####

  ## Visualize the log2 fold change in clusters with a significant differential abundance between the two conditions

  cluster_id <- attr(x, "cluster_id")
  if (is.matrix(cluster_id))
    cluster_id <- cluster_id[, which_cluster_set] %>% drop
  cluster_id <- cluster_id  %>% as.character

  results_column <-
    structure(rep(results_column, length.out = length(contrasts)), .Names = names(contrasts))

  id_map <- make_sample_id_map(x, sample_name_re)
  ids <- names(id_map)[x[, "id"]]
  m1 <- m %>% dplyr::select(1, group) %>%
    dplyr::rename(id = 1)
  u <- sapply(names(contrasts),
    function(a)
    {
      contrasts_a <- plinth::poly_eval(contrasts[[a]])
      glmQLFTestArgs <- utils::modifyList(
        list(glmfit = fit),
        structure(list(contrasts_a), .Names = ifelse(is.matrix(contrasts_a), "contrast", "coef")),
        keep.null = TRUE)

      res_tags <- do.call(edgeR::glmQLFTest, glmQLFTestArgs) %>%
        edgeR::topTags(Inf) %>%
        as.data.frame %>%
        tibble::rownames_to_column("cluster_id")

      rc <- results_column[a] %>% as.vector
      if (is.numeric(rc))
        rc <- colnames(res_tags)[rc]

      diffexp_df <- structure(plinth::dataframe(cluster_id, umap),
        .Names = c("cluster_id", paste0("UMAP", seq(NCOL(umap))))) %>%
        dplyr::left_join(res_tags %>% dplyr::select(cluster_id, !!rc, FDR),
          by = "cluster_id") %>%
        dplyr::rename(cluster = "cluster_id", "logFC" := rc) %>%
        dplyr::mutate(id = ids) %>%
        dplyr::left_join(m1, by = "id") %>%
        dplyr::mutate(
          logFC = dplyr::case_when(FDR > alpha ~ 0.0, TRUE ~ logFC)
        )

      if (diffexp_df %>% dplyr::filter(FDR <= alpha) %>% is_invalid)
        return (NULL)

      col <- structure(
        rep(na.value, diffexp_df %>% dplyr::pull(cluster) %>% unique %>% length),
        .Names = diffexp_df %>% dplyr::pull(cluster) %>% unique
      )
      on <- diffexp_df %>% dplyr::filter(logFC != 0.0) %>% dplyr::pull(cluster) %>% unique
      off <- diffexp_df %>% dplyr::filter(logFC == 0.0) %>% dplyr::pull(cluster) %>% unique
      col[on] <-
        scales::alpha(plot_palette(on %>% length), 1.0)

      diffexp_df %>%
        dplyr::mutate(cluster = dplyr::case_when(cluster %nin% on ~ NA_character_, TRUE ~ cluster) %>% as.factor) %>%
      {
        ggplot2::ggplot(., ggplot2::aes(x = UMAP1, y = UMAP2, color = cluster)) +
        ## N.B. Use additional geometries & the 'breaks' arg to emphasize specific clusters
        ggplot2::scale_color_manual(values = plot_palette(length(unique(.$cluster))), na.value = na.value, labels = labels)
      } +
      ggplot2::geom_point(size = 0.8) + # default: 'shape = 16'
      cowplot::theme_cowplot() +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1), ncol = 2)) +
      ggplot2::ggtitle(sprintf("%s significant clusters", a), subtitle = sprintf("coef: %s", rc)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) -> p1

      p2 <- ggplot2::ggplot(diffexp_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = logFC)) +
        ggplot2::geom_point(alpha = 0.4, size = 0.8) + # default: 'shape = 16'
        ggplot2::scale_colour_gradient2(low = "blue", mid = "gray", high = "red", na.value = na.value) +
        ggplot2::ggtitle(sprintf("%s log2 fold change", a), subtitle = sprintf("coef: %s", rc)) +
        cowplot::theme_cowplot() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      list(p1, p2)
    }, simplify = FALSE) %>% purrr::compact() %>% purrr::flatten()

  if (length(u) == 0)
    return (NULL)

  g6 <- cowplot::plot_grid(plotlist = u, align = "v", ncol = 2)

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      ## Reduce resolution for 'png()' etc. to a manageable value:
      if ("res" %in% names(formals(eval(parse(text = a))))) devices[[a]]$res <- 150
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = plinth::poly_eval(file_path_template) %_% paste0(".", ext)
          )
        )
      )

      print(g6)

      dev.off()
    })
  } else {
    print(g6)
  }
  par(def_par)

  plinth::nop()
}


#' @export
plot_differential_abundance <- function(
  x,
  fits,
  which_cluster_set = NULL,
  ...
)
{
  clusterId <- attr(x, "cluster_id")

  if (is.null(which_cluster_set)) {
    which_cluster_set <- 1
    if (is.matrix(clusterId))
      which_cluster_set <- colnames(clusterId)
  }

  plyr::l_ply(which_cluster_set,
    function(a)
    {
      plyr::l_ply(fits,
        function(fit)
        {
          plot_differential_abundance_single(x, fit = fit, which_cluster_set = a, ...)
        })
    })

  plinth::nop()
}


#' @export
plot_xshift_network <- function(
  cluster_graph,
  image_dir = NULL,
  current_image = 0,
  save_plot = FALSE,
  devices = flowpipe:::graphics_devices,
  file_path_template = sprintf("%s/%03d%s", image_dir, current_image, "_xshift-clusters")
)
{
  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ##### X-Shift clusters #####

  # pdf("D:/Users/priscian/my_documents/urmc/2018/studies/flow/bandyopadhyay-lung/report/images/xshift-clusters.pdf", width = 11.0, height = 8.5)
  # igraph::plot.igraph(cluster_graph, vertex.size = 5, vertex.label.cex = 0.6, main = "Force-directed layout of X-shift clusters")
  # dev.off()

  g7_expr <- expression({
    igraph::plot.igraph(cluster_graph, vertex.size = 3, vertex.label.cex = 0.8, main = "Force-directed layout of X-shift clusters")
  })

  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
          width = 20, height = 20,
            file = plinth::poly_eval(file_path_template) %_% paste0(".", ext)
          )
        )
      )

      plinth::poly_eval(g7_expr)

      dev.off()
    })
  } else {
    plinth::poly_eval(g7_expr)
  }
  par(def_par)

  plinth::nop()
}
