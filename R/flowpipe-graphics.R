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

  if (length(gating_channels) == 1) { # Density plot
    plot.densityArgs <- list(
      x = stats::density(x[, gating_channels[1]]),
      main = "",
      ylab = paste(gating_channels[1], "Density")
    )
    plot.densityArgs <- utils::modifyList(plot.densityArgs, plot..., keep.null = TRUE)

    do.call(plot, plot.densityArgs)

    plinth::vline(
      sprintf("%.2f", cutoffs),
      abline... = list(col = scales::alpha("red", 0.5), lty = "dashed"),
      text... = list(y = plinth::cp_coords()$y)
    )
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

      do.call(graphics::abline, c(ablineArgs, c(v = cutoffs[[gating_channels[1]]])))
      do.call(graphics::abline, c(ablineArgs, c(h = cutoffs[[gating_channels[2]]])))
    }
  }

  plinth::poly_eval(plot_end_callback)

  plinth::nop()
}


#' @export
make_plots <- function(
  x, # "pmm" object from 'get_expression_subset()'
  pmm_files,
  image_dir = NULL,
  devices = list(
    `grDevices::png` = list(height = 8.27, width = 11.69, units = "in", res = 600, type = "cairo", ext = "png"),
    `grDevices::pdf` = list(height = 8.27, width = 11.69, ext = "pdf")
  ),
  save_plot = FALSE,
  breaks = 100,
  metadata,
  cell_counts_re = "",
  cell_counts_callback = NULL, # An expression to e.g. include metadata
  cell_counts_palette = randomcoloR::distinctColorPalette,
  get_expression_subset... = list(),
  channel_dist_palette = randomcoloR::distinctColorPalette,
  heatmap_channels = TRUE,
  umap,
  model_fit, model_coefs,
  metadata_group_var = "group",
  cluster_graph
)
{
  current_image <- 0

  if (save_plot && !dir.exists(image_dir))
    dir.create(image_dir, recursive = TRUE)

  def_par <- par(no.readonly = TRUE)#; dev.off()

  ##### Cell counts #####

  cell_counts <- sapply(pmm_files,
    function(a)
    {
      e <- new.env()
      load(a, envir = e)

      flowCore::exprs(e$tff) %>% NROW
    }, simplify = TRUE)
  names(cell_counts) <- tools::file_path_sans_ext(basename(names(cell_counts))) %>% stringr::str_extract(cell_counts_re)

  ggdf <- plinth::dataframe(sample_id = names(cell_counts), cell_counts = as.numeric(cell_counts))
  ## Use 'cell_counts_callback' to add metadata info to 'ggdf'.
  if (!is.null(cell_counts_callback))
    plinth::poly_eval(cell_counts_callback)

  ## Define colors for groups
  color_conditions <- cell_counts_palette(nlevels(ggdf$group))
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
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]], list(file = paste(image_dir, sprintf("%03d", current_image) %_% "_counts." %_% ext, sep = "/")))
      )

      print(g2)

      dev.off()
    })
  } else {
    print(g2)
  }
  par(def_par)

  ##### Density plots for all samples, all channels #####

  get_expression_subsetArgs <- list(
    x = pmm_files,
    sample_size = 10000
  )
  get_expression_subsetArgs <- utils::modifyList(get_expression_subsetArgs, get_expression_subset..., keep.null = TRUE)

  e <- do.call(get_expression_subset, get_expression_subsetArgs)

  num_col <- NCOL(e) - 1
  num_row <- e[, "id"] %>% as.vector %>% unique %>% length
  channel_colors <- channel_dist_palette(num_col)
  ffn <- basename(pmm_files[e[, "id"] %>% as.vector]) %>% tools::file_path_sans_ext() %>% tools::file_path_sans_ext()
  density_plots <- sapply(seq(num_col), # First column is "id"
    function(i)
    {
      ggplot2::ggplot(as.data.frame(e) %>% dplyr::mutate(id = as.factor(ffn)),
        ggplot2::aes_string(x = plinth::backtick(colnames(e)[i + 1]), y = "id")) +
        ggridges::geom_density_ridges(color = channel_colors[i], fill = scales::alpha(channel_colors[i], 0.4)) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::theme_bw()
    }, simplify = FALSE)

  g3 <- cowplot::plot_grid(plotlist = density_plots, ncol = num_col)

  current_image <- current_image + 1
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
          width = 6.0 * num_col + 1,
          height = 4.0 * num_row + 1,
          file = paste(image_dir, sprintf("%03d", current_image) %_% "_channel_dist." %_% ext, sep = "/")
        ))
      )

      print(g3)

      dev.off()
    })
  } else {
    print(g3)
  }
  par(def_par)

  ##### Cell counts after initial gating #####

  g4 <- plot_event_counts(ggdf, "gated_cell_counts") +
    ggplot2::labs(title = "Event counts after gating")

  current_image <- current_image + 1
  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]], list(file = paste(image_dir, sprintf("%03d", current_image) %_% "_gated_counts." %_% ext, sep = "/")))
      )

      print(g4)

      dev.off()
    })
  } else {
    print(g4)
  }
  par(def_par)

  ##### Heatmap of channels vs. clusters #####

  ## 'cluster_matrix' contains the median expression for each channel for each cluster. Each row (each channel) is normalized by
  ## 'base::scale()' to allow relative comparison of normalized median expression among all clusters for each channel. This
  ## mostly disregards variance of the marker expression in each cluster, but offers an overview of each cluster's composition.
  sub_matrix <- x[, heatmap_channels]
  cluster_id <- attr(x, "cluster_id")
  cluster_matrix <- NULL
  for(i in sort(unique(cluster_id))) {
    cluster_matrix <- rbind(cluster_matrix, matrixStats::colMedians(sub_matrix[cluster_id == i, ]))
  }
  colnames(cluster_matrix) <- heatmap_channels
  rownames(cluster_matrix) <- sort(unique(cluster_id))
  par(mar = c(2, 2, 2, 2))

  g5_expr <- expression({
    gplots::heatmap.2(
      t(cluster_matrix),
      col = gplots::bluered(100),
      trace = "none", density.info = "none",
      sepcolor = "white", sepwidth = c(0.001, 0.001),
      colsep = c(1:NCOL(t(cluster_matrix))),
      rowsep = c(1:NROW(t(cluster_matrix))),
      xlab = "cluster", ylab = "channel", scale = scale_
    )
  })

  ## Center & scale heatmap values in the row (i.e. channel) direction
  scale_ <- "row"

  current_image <- current_image + 1
  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = paste(image_dir, sprintf("%03d", current_image) %_% "_heatmap_channels-clusters-row." %_% ext, sep = "/")
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

  current_image <- current_image + 1
  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = paste(image_dir, sprintf("%03d", current_image) %_% "_heatmap_channels-clusters-column." %_% ext, sep = "/")
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

  ##### UMAP visualization of significant differential abundance between conditions #####

  ## Visualize the log2 fold change in clusters with a significant differential abundance between the two conditions

  cluster_id <- attr(x, "cluster_id") %>% as.character
  id_map <- structure(attr(x, "id_map"),
    .Names = names(attr(x, "id_map")) %>% basename %>% stringr::str_extract(cell_counts_re))
  ids <- names(id_map)[x[, "id"]]
  m <- metadata %>% dplyr::select(1, !!metadata_group_var) %>%
    dplyr::rename(id = 1)
  u <- sapply(model_coefs,
    function(a)
    {
      res_tags <- edgeR::glmQLFTest(fit, coef = a) %>%
        edgeR::topTags(Inf) %>%
        as.data.frame %>%
        tibble::rownames_to_column("cluster_id")

      diffexp_df <- structure(plinth::dataframe(cluster_id, umap),
        .Names = c("cluster_id", paste0("UMAP", seq(NCOL(umap))))) %>%
        dplyr::left_join(res_tags %>% dplyr::select(cluster_id, logFC, FDR),
          by = "cluster_id") %>%
        dplyr::rename(cluster = "cluster_id") %>%
        dplyr::mutate(id = ids) %>%
        dplyr::left_join(m, by = "id") %>%
        dplyr::mutate(
          logFC = dplyr::case_when(FDR > 0.05 ~ 0.0, TRUE ~ logFC)
        )

      if (diffexp_df %>% dplyr::filter(FDR <= 0.05) %>% is_invalid)
        return (NULL)

      col <- structure(
        rep(scales::alpha("lightgray", 0.3), diffexp_df %>% dplyr::pull(cluster) %>% unique %>% length),
        .Names = diffexp_df %>% dplyr::pull(cluster) %>% unique
      )
      on <- diffexp_df %>% dplyr::filter(logFC != 0.0) %>% dplyr::pull(cluster) %>% unique
      off <- diffexp_df %>% dplyr::filter(logFC == 0.0) %>% dplyr::pull(cluster) %>% unique
      col[on] <-
        scales::alpha(colorspace::qualitative_hcl(on %>% length, palette = "Dark 2"), 0.7)
      p1 <- ggplot2::ggplot(diffexp_df %>% dplyr::filter(cluster %in% off), ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(ggplot2::aes(color = cluster), alpha = 0.4, size = 0.5, show.legend = FALSE) +
        ggplot2::geom_point(data = diffexp_df %>% dplyr::filter(cluster %in% on),
          ggplot2::aes(color = cluster), alpha = 0.4, size = 0.5, show.legend = TRUE) +
        ## N.B. Use additional geometries & the 'breaks' arg to emphasize specific clusters
        ggplot2::scale_color_manual(values = col, breaks = on) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1), ncol = 2)) +
        ggplot2::ggtitle(sprintf("%s significant X-Shift clusters", a)) +
        cowplot::theme_cowplot() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


      p2 <- ggplot2::ggplot(diffexp_df, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(ggplot2::aes(color = logFC), alpha = 0.4, size = 0.5) +
        ggplot2::scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
        ggplot2::ggtitle(sprintf("%s log2 fold change", a)) +
        cowplot::theme_cowplot() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      list(p1, p2)
    }, simplify = FALSE) %>% purrr::compact() %>% purrr::flatten()

  g6 <- cowplot::plot_grid(plotlist = u, align = "v", ncol = 2)

  current_image <- current_image + 1
  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
            width = 11 * 2, height = 8.5 * 2,
            file = paste(image_dir, sprintf("%03d", current_image) %_% "_umap_significant-clusters." %_% ext, sep = "/")
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

  ##### X-Shift clusters #####

  # pdf("D:/Users/priscian/my_documents/urmc/2018/studies/flow/bandyopadhyay-lung/report/images/xshift-clusters.pdf", width = 11.0, height = 8.5)
  # igraph::plot.igraph(cluster_graph, vertex.size = 5, vertex.label.cex = 0.6, main = "Force-directed layout of X-shift clusters")
  # dev.off()

  g7_expr <- expression({
    igraph::plot.igraph(cluster_graph, vertex.size = 3, vertex.label.cex = 0.8, main = "Force-directed layout of X-shift clusters")
  })

  current_image <- current_image + 1
  if (save_plot) {
    plyr::l_ply(names(devices),
    function(a)
    {
      ext <- devices[[a]]$ext; devices[[a]]$ext <- NULL
      do.call(eval(parse(text = a)),
        modifyList(devices[[a]],
          list(
          width = 20, height = 20,
            file = paste(image_dir, sprintf("%03d", current_image) %_% "_xshift-clusters." %_% ext, sep = "/")
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

  #browser()
}
