#' @export
process_config_file <- function(
  path,
  storage_env = NULL, # globalenv(),
  get_initial_list_fcs_files... = list(),
  bead_normalize... = list(),
  list_files... = list(),
  rename_fcs_parameters_name_desc... = list(),
  debarcode... = list(),
  guess_channels... = list(
    channels_regex = list(
      flow = "",
      cytof =
        "^(?!.*(_))|^(?=.*(time|event_length))|beads|_eq|barcode|bkg|bckg|background|intercalator|viability|live-dead|cisplatin",
      spectral = ""
    )
  ),
  channel_sheet_path = NULL,
  skip_interactive_channel_sheet_exists = TRUE,
  prepare_metadata... = list(),
  metadata_sheet_path = NULL,
  metadata_required_columns = c("id", "group"),
  skip_interactive_metadata_sheet_exists = TRUE,
  impute_metadata = TRUE,
  b = NULL,
  plot_channel_densities_by_sample... = list(),
  prepare_augmented_fcs_data... = list(),
  get_expression_subset... = list(),
  plot_cell_counts... = list(),
  make_metaclusters... = list(),
  make_clusters... = list(),
  merge_clusters... = list(),
  do_differential_expression... = list(),
  test_contrasts... = list(),
  make_umap_embedding... = list(),
  summarize_all_clusters... = list(),
  plot_common_umap_viz... = list(),
  plot_heatmaps... = list(),
  plot_differential_abundance... = list(),
  export_cluster_summary... = list(),
  split_pmm_by_cluster... = list(),
  ...
)
{
  flowpipe_interactive_off <-
    !(is.null(getOption("flowpipe_interactive_off")) || !getOption("flowpipe_interactive_off"))

  if (missing(path) || is_invalid(path)) {
    if (interactive() && !flowpipe_interactive_off) {
      msg <-
r"---{
Choose the TOML configuration file used to describe this analysis. It should
contain a subsection "fcs_files" that gives paths to individual FCS files or
to a directory containing those files.
}---"
      message(msg); utils::flush.console()
      path <- svDialogs::dlg_open(title = "Open flowpipe configuration file",
        filters = c("TOML files (*.toml)", "*.toml"))$res

      if (is_invalid(path))
        stop("flowpipe is unable to find an experimental configuration file")
    } else {
      stop("A valid 'path' value must be provided to run this function non-interactively")
    }
  }

  #i <- RcppTOML::parseTOML(path, escape = FALSE)
  i <- blogdown::read_toml(path)
  if (!is_invalid(storage_env)) assign("config_toml", i, envir = storage_env)
  ii <- rlang::duplicate(i, shallow = FALSE)

  clear_1_cache <- FALSE

  ####################
  ### Setup
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "setup")

  type <- i$setup$type %>% stringr::str_trim() %>% tolower %>% substr(0, 1) %>%
    { c(f = "flow", c = "cytof", s = "spectral")[.] } %>% as.vector
  if (is_invalid(type)) {
    warning("Cytometry type is missing or misspecified; defaulting to \"flow\"", immediate. = TRUE)

    type <- "flow"
  }

  data_dir <- i$setup$data_dir %>% keystone::normalize_path()
  if (is_invalid(data_dir) || stringr::str_trim(data_dir) == "") {
    data_dir <- "./data"
    warning("Missing 'data_dir' has been set to the default directory", immediate. = TRUE)
  }
  if (!is_invalid(storage_env)) assign("data_dir", data_dir, envir = storage_env)
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  assign(".cm", memoise::cache_filesystem(path = data_dir), envir = globalenv())
  memoise_all_package_functions(key = as.character(i$setup$key))

  ## Directory for saving reports
  report_dir <- i$setup$report_dir
  if (is_invalid(report_dir) || stringr::str_trim(report_dir) == "") {
    report_dir <- "./report"
    warning("Missing 'report_dir' has been set to the default directory", immediate. = TRUE)
  }
  if (!is_invalid(storage_env)) assign("report_dir", report_dir, envir = storage_env)

  ## Directory for saving images
  image_dir <- i$setup$image_dir
  if (is_invalid(image_dir) || stringr::str_trim(image_dir) == "") {
    image_dir <- paste(report_dir, "images", sep = "/")
    warning("Missing 'image_dir' has been set to the default directory", immediate. = TRUE)
  }
  if (!is_invalid(storage_env)) assign("image_dir", image_dir, envir = storage_env)

  fcs_files <- readLines(textConnection(i$setup$fcs_files), skipNul = TRUE) %>%
    keystone::normalize_path() %>% stringr::str_trim() %>% stringi::stri_remove_empty()
  ## Are any of these paths directories? If so, find the FCS files in them.
  ## Arg 'list_files... = list(recursive = 𝘯)' or 'recursive = TRUE' to traverse subdirectories
  get_initial_list_fcs_files <- (function(...)
  {
    fcs_files %<>% sapply(
      function(a)
      {
        if (!fs::is_dir(a))
          return (a)

        list_filesArgs <- list(
          path = a,
          pattern = "\\.fcs$",
          ignore.case = TRUE,
          absolute = TRUE
        )
        list_filesArgs <- utils::modifyList(list_filesArgs, list_files..., keep.null = TRUE)

        do.call(keystone::list_files, list_filesArgs)
      }, simplify = TRUE, USE.NAMES = FALSE) %>% keystone::normalize_path() %>% unique

      fcs_files
  }) %>% memoise::memoise(cache = .cm) %>%
    keystone::patch_memoised_for_subcaching(name = "get_initial_list_fcs_files", key = as.character(i$setup$key))
  get_initial_list_fcs_filesArgs <- list(SUFFIX = "_initial-fcs", clear_1_cache = clear_1_cache)
  get_initial_list_fcs_filesArgs <-
    utils::modifyList(get_initial_list_fcs_filesArgs,
      get_initial_list_fcs_files..., keep.null = TRUE)

  do_stop_check(get_initial_list_fcs_filesArgs, stop_var = "STOP")
  fcs_files <- do.call(get_initial_list_fcs_files, get_initial_list_fcs_filesArgs)
  if (!is_invalid(storage_env)) assign("fcs_files", fcs_files, envir = storage_env)


  ## Check whether 'memoise()'d function should be rerun:
  get_caching_status <- function(args_list, status = clear_1_cache)
  {
    if (!is_invalid(args_list$clear_1_cache) && is.logical(args_list$clear_1_cache)) {
      status <- FALSE
      if (args_list$clear_1_cache)
        status <- TRUE
    }

    if (is_invalid(status))
      status <- FALSE

    status
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(get_initial_list_fcs_files...)

  ####################
  ### Data preparation
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "data_prep")

  fcs_output_dir <- i$data_prep$fcs_output_dir
  if (!is_invalid(storage_env)) assign("fcs_output_dir", fcs_output_dir, envir = storage_env)

  fcs_files_01 <- fcs_files

  ### Bead normalization for CyTOF (if necessary)

  norm_beads <- i$data_prep$norm_beads %>% stringr::str_trim() %>% tolower %>% substr(0, 1) %>%
    { c(d = "dvs", b = "beta")[.] } %>% as.vector
  if (!is_invalid(norm_beads)) {
    bead_normalizeArgs <- list(
      input_path = fcs_files,
      output_dir = paste(fcs_output_dir, "bead-normalized", sep = "/"),
      outfile_suffix = "-bead-normalized",
      ## "dvs" (for bead masses 140, 151, 153, 165, 175) or
      ##   "beta" (for masses 139, 141, 159, 169, 175) or numeric vector of masses:
      beads = "dvs",
      normCytof... = list(plot = FALSE),
      SUFFIX = "_beadnorm", clear_1_cache = clear_1_cache
    )
    bead_normalizeArgs <-
      utils::modifyList(bead_normalizeArgs, bead_normalize..., keep.null = TRUE)

    do_stop_check(bead_normalizeArgs, stop_var = "STOP")
    fcs_files_01 <- do.call(bead_normalize, bead_normalizeArgs)

    ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
    clear_1_cache <- get_caching_status(bead_normalize...)
  }
  if (!is_invalid(storage_env)) assign("fcs_files_01", fcs_files_01, envir = storage_env)

  ## Check channel names & descriptions among FCS files
  ## Add argument 'clear_1_cache = TRUE' to clear cache & rerun
  cbs <- get_channels_by_sample(x = fcs_files_01, SUFFIX = "_cbs", clear_1_cache = clear_1_cache)
  if (!is_invalid(storage_env)) assign("cbs", cbs, envir = storage_env)

  ## Create common channel names & descriptions among FCS files; save new FCS files
  name_remap <- sapply(i$data_prep$name_remap,
    function(a) { r <- a; if (stringr::str_trim(r) == "NA") r <- NA_character_; r })
  desc_remap <- c(attr(cbs, "desc_remap"), i$data_prep$description_remap %>% unlist)
  rename_fcs_parameters_name_descArgs <- list(
    channels_by_sample = cbs,
    desc_remap = desc_remap,
    name_remap = name_remap,
    output_dir = paste(fcs_output_dir, "renamed", sep = "/"),
    SUFFIX = "_rename-fcs", clear_1_cache = clear_1_cache
  )
  rename_fcs_parameters_name_descArgs <-
    utils::modifyList(rename_fcs_parameters_name_descArgs,
      rename_fcs_parameters_name_desc..., keep.null = TRUE)

  do_stop_check(rename_fcs_parameters_name_descArgs, stop_var = "STOP")
  fcs_files_01 <- do.call(rename_fcs_parameters_name_desc, rename_fcs_parameters_name_descArgs)
  if (!is_invalid(storage_env)) assign("fcs_files_01", fcs_files_01, envir = storage_env)

  if (interactive() && !flowpipe_interactive_off) {
    msg <-
r"---{
Here's a current, numbered list of this experiment's FCS files so far
(possibly after channel reorganizing & renaming). If debarcoding of some files
is necessary, use these numbers to match barcoding keys to the file names:
}---" %>%
      paste("\n", paste(seq_along(fcs_files_01) %_% ".", basename(fcs_files_01), collapse = "\n"),
        "\n", sep = "")
    message(msg); utils::flush.console()
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(rename_fcs_parameters_name_desc...)

  ## Check that renames/deletions worked; should produce no warnings
  cbs_check <-
    get_channels_by_sample(x = fcs_files_01, SUFFIX = "_cbs-check", clear_1_cache = clear_1_cache)
  if (!is_invalid(storage_env)) assign("cbs_check", cbs_check, envir = storage_env)

  if (any_noncommon_descs(cbs_check))
    warning("FCS files may still not share common channel descriptions", immediate. = TRUE)

  ### Debarcoding (if necessary)

  ## Create memoise()'d function that works in "debarcode" package namespace
  debarcode_fun <<- debarcode::debarcode %>%
    memoise::memoise(cache = .cm) %>%
    keystone::patch_memoised_for_subcaching(name = "debarcode")
  ## N.B. To wholly replace 'debarcode::debarcode' for this session, v.
  ## https://stackoverflow.com/questions/24331690/modify-package-function/58238931#58238931

  fcs_files_02 <- fcs_files_01

  if (!is_invalid(i$data_prep$barcoding_keys)) {
    fcsIndexRe <- "^\\s*\\[\\s*(\\d+)\\s*\\]\r*\n"
    flit <- i$data_prep$barcoding_keys %>% stringr::str_split_1("(\r*\n){2,}")
    fcs_index <- flit %>% stringr::str_extract(fcsIndexRe, group = 1) %>% readr::parse_number()
    barcoding_keys <- flit %>% stringr::str_remove(fcsIndexRe) %>%
      sapply(function(a) read.table(text = a, header = TRUE, check.names = FALSE), simplify = FALSE) %>%
      structure(.Names = fcs_files_01[fcs_index])

    debarcodeArgs <- list(
      input_path = fcs_files_01,
      key = barcoding_keys,
      output_dir = paste(fcs_output_dir, "deconvoluted", sep = "/"),
      SUFFIX = "_debarcoding", clear_1_cache = clear_1_cache
    )
    debarcodeArgs <- utils::modifyList(debarcodeArgs, debarcode..., keep.null = TRUE)

    do_stop_check(debarcodeArgs, stop_var = "STOP")
    debarcoding_details <- do.call(debarcode_fun, debarcodeArgs); gc()
    ## Debarcoding event counts:
    # attr(debarcoding_details, "details") %>% sapply(function(a) attr(a, "sample_id") %>% table)

    fcs_files_02 <- debarcoding_details %>% as.vector

    ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
    clear_1_cache <- get_caching_status(debarcode...)
  }
  if (!is_invalid(storage_env)) assign("fcs_files_02", fcs_files_02, envir = storage_env)

  ####################
  ### Channel selection
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "channels")

  ## Prepare a spreadsheet that selects channels for subsequent analyses
  guess_channelsArgs <- list(
    cbs = cbs_check,
    type = type
    #SUFFIX = "_channels", clear_1_cache = clear_1_cache
  )
  guess_channelsArgs <- utils::modifyList(guess_channelsArgs, guess_channels..., keep.null = TRUE)
  if (!is_invalid(i$channels$channels_regex))
    guess_channelsArgs$channels_regex[[type]] <- i$channels$channels_regex

  channel_sheet <- do.call(guess_channels, guess_channelsArgs)
  if (!is_invalid(storage_env)) assign("channel_sheet", channel_sheet, envir = storage_env)

  if (is_invalid(channel_sheet_path))
    channel_sheet_path <- paste(data_dir, "experiment-channels.xlsx", sep = "/")

  ## Have interactive user review 'channel_sheet':
  if (interactive() && !flowpipe_interactive_off &&
    !(skip_interactive_channel_sheet_exists && fs::file_exists(channel_sheet_path))) {
    msg <-
r"---{
------------------------------------------------------
The channels/markers sheet for this analysis is ready.
------------------------------------------------------
You can:

  • REVIEW or edit the new channels/markers sheet
  • SAVE this sheet to file system for review or external editing
  • CONTINUE flowpipe analysis using saved parameters sheet
  • QUIT flowpipe for now
}---"

    message(msg); utils::flush.console()
    repeat {
      choice <- keystone::ask_multiple_choice(c("Review", "Save", "Continue", "Quit"))

      ## React to each possible choice:
      switch(choice$lowercase,
        quit = {
          return (invisible(NULL))
        },

        review = {
          msg <-
r"---{
After making changes to the channels/markers sheet, click "synchronize" button,
then click "Done". If no changes, click "Done" and return to the R console.
}---"
          message(msg); utils::flush.console()
          keystone::press_enter_to_continue()

          channel_sheet <- DataEditR::data_edit(channel_sheet)
        },

        save = {
          # defaultFileName <- sprintf("experiment-channels_%s",
          #   keystone::make_current_timestamp(use_seconds = TRUE, seconds_sep = "+"))
          # params_path <- svDialogs::dlg_save(default = defaultFileName, title = "Save channels/markers sheet",
          #   filters = svDialogs::dlg_filters[c("xls", "csv"), ])$res

          rio::export(channel_sheet, channel_sheet_path)

          msg <- paste0(
r"---{
The channels/marker spreadsheet has been generated & can be found here:

}---",
          channel_sheet_path)
          message(msg); utils::flush.console()
        },

        continue = {
          if(!fs::file_exists(channel_sheet_path)) {
            rio::export(channel_sheet, channel_sheet_path)
            message("Channels/markers sheet was saved to file system for further use")

            msg <- paste0(
r"---{
The channels/marker spreadsheet has been generated & can be found here:

}---",
            channel_sheet_path)
            message(msg); utils::flush.console()
          } else {
            warning("Retrieving channels/markers sheet from file system",
              immediate. = TRUE)
          }

          break
        }
      )
    }
  } else if (!fs::file_exists(channel_sheet_path)) {
    rio::export(channel_sheet, channel_sheet_path)
  }

  channel_sheet <- rio::import(channel_sheet_path) %>%
    dplyr::mutate(
      use = stringr::str_trim(use),
      use =
        dplyr::case_when(use == "" | is.na(use) ~ NA_character_, TRUE ~ formals(guess_channels)$usage_marker)
    )
  if (!is_invalid(storage_env)) assign("channel_sheet", channel_sheet, envir = storage_env)

  ####################
  ### Prepare metadata
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "metadata")

  ## Prepare a spreadsheet that provides metadata for subsequent analyses
  if (is_invalid(metadata_sheet_path))
    metadata_sheet_path <- paste(data_dir, "experiment-metadata.xlsx", sep = "/")

  prepare_metadataArgs <- list(
    x = fcs_files_02
    #SUFFIX = "_metadata", clear_1_cache = clear_1_cache
  )
  prepare_metadataArgs <-
    utils::modifyList(prepare_metadataArgs, prepare_metadata..., keep.null = TRUE)
  if (!is_invalid(i$metadata$table)) {
    metadata_sheet <- read.table(text = i$metadata$table, header = TRUE, check.names = FALSE)

    if (!is_invalid(i$metadata$overwrite_spreadsheet) && i$metadata$overwrite_spreadsheet) {
      rio::export(metadata_sheet, metadata_sheet_path)
    }
  } else {
    metadata_sheet <- do.call(prepare_metadata, prepare_metadataArgs)
  }
  if (!is_invalid(storage_env)) assign("metadata_sheet", metadata_sheet, envir = storage_env)

  ## Warn that metadata "id" must contain a unique character sequence from each FCS file
  if (interactive() && !flowpipe_interactive_off) {
    msg <-
r"---{
Here's a current, numbered list of this experiment's FCS files so far. When
editing the metadata file, be sure that each "id" variable matched to a file
contains a UNIQUE character sequence from the corresponding file name:
}---" %>%
      paste("\n", paste(seq_along(fcs_files_02) %_% ".", basename(fcs_files_02), collapse = "\n"),
        "\n", sep = "")
    message(msg); utils::flush.console()
  }

  ## Have interactive user review 'channel_sheet':
  if (interactive() && !flowpipe_interactive_off &&
    !(skip_interactive_metadata_sheet_exists && fs::file_exists(metadata_sheet_path))) {
    msg <- c(paste0(
r"---{
------------------------------------------------------------------------------
The metadata sheet for this analysis is ready for review/editing. It MUST have}---",
sprintf("\ncolumns %s, but it can also include any additional\n",
  stringr::str_flatten_comma(dQuote(metadata_required_columns, q = FALSE))),
r"---{clinicopathological data suitable for differential-analysis modeling.
------------------------------------------------------------------------------
You can:

  • REVIEW or edit the new metadata sheet
  • SAVE this sheet to file system for review or external editing
  • CONTINUE flowpipe analysis using saved metadata sheet
  • QUIT flowpipe for now
}---"))

    message(msg); utils::flush.console()
    repeat {
      choice <- keystone::ask_multiple_choice(c("Review", "Save", "Continue", "Quit"))

      ## React to each possible choice:
      switch(choice$lowercase,
        quit = {
          return (invisible(NULL))
        },

        review = {
          msg <-
r"---{
After making changes to the metadata sheet, click "synchronize" button,
then click "Done". If no changes, click "Done" and return to the R console.
}---"
          message(msg); utils::flush.console()
          keystone::press_enter_to_continue()

          metadata_sheet <- DataEditR::data_edit(metadata_sheet)
        },

        save = {
          # defaultFileName <- sprintf("experiment-metadata_%s",
          #   keystone::make_current_timestamp(use_seconds = TRUE, seconds_sep = "+"))
          # params_path <- svDialogs::dlg_save(default = defaultFileName, title = "Save metadata sheet",
          #   filters = svDialogs::dlg_filters[c("xls", "csv"), ])$res

          rio::export(metadata_sheet, metadata_sheet_path)

          msg <- paste0(
r"---{
The metadata spreadsheet has been generated & can be found here:

}---",
          metadata_sheet_path)
          message(msg); utils::flush.console()
        },

        continue = {
          if(!fs::file_exists(metadata_sheet_path)) {
            rio::export(metadata_sheet, metadata_sheet_path)
            message("Metadata sheet was saved to file system for further use")

            msg <- paste0(
r"---{
The metadata spreadsheet has been generated & can be found here:

}---",
            metadata_sheet_path)
            message(msg); utils::flush.console()
          } else {
            warning("Retrieving metadata sheet from file system", immediate. = TRUE)
          }

          ## Does 'metadata_sheet' have all the requisite columns?
          if (!all(metadata_required_columns %in% names(metadata_sheet))) {
            warning(sprintf('Metadata does not contain required columns %s; please correct it',
              stringr::str_flatten_comma(dQuote(metadata_required_columns, q = FALSE))),
              immediate. = TRUE)

            next
          }

          break
        }
      )
    }
  } else if (!fs::file_exists(metadata_sheet_path)) {
    rio::export(metadata_sheet, metadata_sheet_path)
  }

  metadata_sheet <- rio::import(metadata_sheet_path) %>%
    dplyr::mutate(
      across(where(~ is.character(.x)) & !matches("^id$"), as.factor)
    )
  ## If variables' reference levels are provided in config sheet, set them:
  if (!is_invalid(i$metadata$reference_levels)) {
    plyr::l_ply(names(i$metadata$reference_levels),
      function(a)
      {
        metadata_sheet <<- metadata_sheet %>%
          dplyr::mutate(!!a := forcats::fct_relevel(.[[a]], i$metadata$reference_levels[[a]]))
      })
  }

  metadata <- metadata_sheet
  if (!is_invalid(storage_env)) assign("metadata", metadata, envir = storage_env)

  ####################
  ### Some "constants"
  ####################

  ## arcsinh scale parameter
  if (is_invalid(b)) {
    b <- switch(type,
      spectral = 1/150, # same as FCM, but this is a guess
      cytof = 1/8,
      1/150 # default, FCM
    )
  }
  if (!is_invalid(storage_env)) assign("b", b, envir = storage_env)

  ## Counter
  flowpipe_params(current_image_set = 0)

  ## Channels
  primary_channels <- channel_sheet %>% dplyr::filter(use == "X") %>%
    { structure(dplyr::pull(., name), .Names = dplyr::pull(., description)) }
  if (!is_invalid(storage_env)) assign("primary_channels", primary_channels, envir = storage_env)

  ### Density plots

  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  plot_channel_densities_by_sampleArgs <- list(
    x = fcs_files_02,
    image_dir = image_dir,
    current_image = flowpipe_params("current_image_set"),
    save_plot = TRUE,
    get_fcs_expression_subset... = list(
      b = b,
      channels_subset = primary_channels,
      channels_by_sample = cbs_check
    ),
    SUFFIX = "_channel-densities", clear_1_cache = clear_1_cache
  )
  plot_channel_densities_by_sampleArgs <-
    utils::modifyList(plot_channel_densities_by_sampleArgs,
      plot_channel_densities_by_sample..., keep.null = TRUE)

  do_stop_check(plot_channel_densities_by_sampleArgs, stop_var = "STOP")
  do.call(plot_channel_densities_by_sample, plot_channel_densities_by_sampleArgs)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(plot_channel_densities_by_sample...)

  ####################
  ### Prepare aggregate analysis data set
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "aggregate_data")

  remove_outliers <- i$aggregate_data$remove_outliers
  if (is_invalid(remove_outliers)) remove_outliers <- TRUE

  aggregate_files_parallel <- i$aggregate_data$parallel_processing
  options_keystone_parallel <- getOption("keystone_parallel")
  if (!is_invalid(aggregate_files_parallel) && !aggregate_files_parallel) {
    options(keystone_parallel = FALSE)
  }

  ## Make augmented flowCore::flowFrame objects.
  prepare_augmented_fcs_dataArgs <- list(
    x = fcs_files_02,
    b = b,
    data_dir = paste(data_dir, "samples", sep = "/"),
    remove_outliers = remove_outliers,
    pmm_channels = primary_channels,
    multisect... = list(max_sample = 5000),
    outfile_suffix = NULL,
    outfile_prefix = expression(outfile_prefix <- paste0(x %>% dirname %>% basename, "-")),
    overwrite = FALSE,
    SUFFIX = "_pmm-files", clear_1_cache = clear_1_cache
  )
  prepare_augmented_fcs_dataArgs <-
    utils::modifyList(prepare_augmented_fcs_dataArgs,
      prepare_augmented_fcs_data..., keep.null = TRUE)

  do_stop_check(prepare_augmented_fcs_dataArgs, stop_var = "STOP")
  pmm_files <-
    do.call(prepare_augmented_fcs_data, prepare_augmented_fcs_dataArgs,
      quote = any(sapply(prepare_augmented_fcs_dataArgs, is.expression))); gc()
  if (!is_invalid(storage_env)) assign("pmm_files", pmm_files, envir = storage_env)

  ## Restore state of parallel processing to before this function call:
  options(keystone_parallel = options_keystone_parallel)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(prepare_augmented_fcs_data...)


  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "pregating")
  i$pregating$stop <- NULL

  ## Create "pregating" expressions from the pseudocode in the config file.
  ## N.B. Currently supports 1 or 2 channels, & or | comparison, in-/equalities & "between" ranges.
  pg_channel_expr_template <- "x[, \"%s\"] %%>%% { %s }"
  pg_channel_expr_template_true <- "x[, \"%s\"] %%>%% { rep(%s, length(.)) }"
  pregating_strategies <- sapply(i$pregating,
    function(a)
    {
      compare <- a$compare; a$compare <- NULL
      strategy <- sapply(names(a),
        function(b)
        {
          template <- pg_channel_expr_template
          if (is.logical(a[[b]]))
            template <- pg_channel_expr_template_true
          ## Turn pseudocode into R code:
          r <- stringr::str_trim(a[[b]]) %>%
            stringr::str_replace_all("(\\()", "\\1., ") %>%
            stringr::str_replace_all("(between)\\s+(.*)", "dplyr::\\1(., \\2)") %>%
            stringr::str_replace_all("(\\<|\\>|\\!\\=|\\=\\=)", ". \\1")

          sprintf(template, b, r)
        }, simplify = FALSE)

      ## Turn pseudocode into valid R expression:
      parse(text = (strategy %>% unlist(use.names = FALSE) %>%
        paste(collapse = stringr::str_pad(compare[1L], 3, side = "both"))))
    }, simplify = FALSE)

  max_events_per_sample <- i$aggregate_data$max_events_per_sample
  if (is_invalid(max_events_per_sample)) max_events_per_sample <- Inf

  ## Get stacked subsets of augmented expression matrices IDed by sample.
  get_expression_subsetArgs <- list(
    x = pmm_files,
    gate... = list(strategy = pregating_strategies),
    save_plot... = list(file = paste(image_dir, "sample-pregating.pdf", sep = "/")), # 'file' must be given to save PDF
    sample_size = max_events_per_sample,
    SUFFIX = "_expression-subset", clear_1_cache = clear_1_cache
  )
  get_expression_subsetArgs <-
    utils::modifyList(get_expression_subsetArgs, get_expression_subset..., keep.null = TRUE)

  do_stop_check(get_expression_subsetArgs, stop_var = "STOP")
  e <- do.call(get_expression_subset, get_expression_subsetArgs)
  if (!is_invalid(storage_env)) assign("e", e, envir = storage_env)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(get_expression_subset...)


  ### Add sample IDs in attributes of aggregate 'e' matrix to metadata
  sampleId <- adist(names(attributes(e)$id_map), metadata$id) %>% apply(1, which.min)
  if (any(duplicated(sampleId)))
    error("Metadata contains non-unique IDs or IDs incorrectly matched to sample file names")

  metadata %<>% dplyr::mutate(sample_id = sampleId, .before = 1)
  metadata_i <- rlang::duplicate(metadata, shallow = FALSE)
  ## Impute metadata
  if (!is_invalid(i$metadata$impute) && i$metadata$impute) {
    metadata_i %<>%
      { mice::complete(mice::mice(., printFlag = FALSE),
        action = 1) } # 'action = n' returns nth imputation
  }

  if (!is_invalid(storage_env)) assign("metadata", metadata, envir = storage_env)
  if (!is_invalid(storage_env)) assign("metadata_i", metadata_i, envir = storage_env)


  ### Plot sample cell counts before & after gating for comparison at this point

  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  plot_cell_countsArgs <- list(
    x = e,
    ## If bead normalization, use starting FCS files before normalization:
    pmm_files = dplyr::case_when(!is_invalid(norm_beads) ~ fcs_files, .default = pmm_files),
    m = metadata,
    sample_name_re = NULL, # E.g. "D\\d{3}-.*-\\d{1,2}", or use 'm$id' when NULL
    image_dir = image_dir,
    current_image = flowpipe_params("current_image_set"),
    save_plot = TRUE,
    #devices = flowpipe:::graphics_devices["grDevices::pdf"],
    SUFFIX = "_cell-counts", clear_1_cache = clear_1_cache
  )
  plot_cell_countsArgs <-
    utils::modifyList(plot_cell_countsArgs, plot_cell_counts..., keep.null = TRUE)

  do_stop_check(plot_cell_countsArgs, stop_var = "STOP")
  do.call(plot_cell_counts, plot_cell_countsArgs)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(plot_cell_counts...)

  ####################
  ### Clustering
  ####################

  i$metaclustering$stop <- NULL

  do_metaclustering <- i$metaclustering$do_metaclustering
  i$metaclustering$do_metaclustering <- NULL

  i_metaclustering <- i$metaclustering[1]
  if (!is_invalid(i_metaclustering) &&
    !(is_invalid(names(i_metaclustering)) || stringr::str_trim(names(i_metaclustering)) == "")) {
    make_metaclusters_dots <-
      utils::modifyList(list(method = names(i_metaclustering)), i_metaclustering[[1]],
        keep.null = TRUE)
  } else {
    make_metaclusters_dots <- list(method = "Rphenograph", Rphenograph_k = 15)

    warning("Using default metaclustering parameters", immediate. = TRUE)
  }

  i$clustering$stop <- NULL

  label_threshold <- i$clustering$label_threshold
  i$clustering$label_threshold <- NULL

  i_clustering <- i$clustering[1]
  if (!is_invalid(i_clustering) &&
    !(is_invalid(names(i_clustering)) || stringr::str_trim(names(i_clustering)) == "")) {
    make_clusters_dots <-
      utils::modifyList(list(method = names(i_clustering)), i_clustering[[1]],
        keep.null = TRUE)
  } else {
    make_clusters_dots <- list(method = "Rphenograph", Rphenograph_k = 50)

    warning("Using default clustering parameters", immediate. = TRUE)
  }

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "metaclustering")

  ## Metaclustering
  cluster_id_mc <- NULL
  if (!is_invalid(do_metaclustering) && do_metaclustering) {
    make_metaclustersArgs <- list(
      x = e,
      channels = primary_channels,
      make_clusters... = make_clusters_dots,
      make_metaclusters... = make_metaclusters_dots,
      #centroid_fun = mean,
      SUFFIX = "_cluster-id-mc", clear_1_cache = clear_1_cache
    )
    make_metaclustersArgs <-
      utils::modifyList(make_metaclustersArgs, make_metaclusters..., keep.null = TRUE)

    do_stop_check(make_metaclustersArgs, stop_var = "STOP")
    cluster_id_mc <- do.call(make_metaclusters, make_metaclustersArgs)

    ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
    clear_1_cache <- get_caching_status(make_metaclusters...)
  }
  if (!is_invalid(storage_env)) assign("cluster_id_mc", cluster_id_mc, envir = storage_env)

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "clustering")

  ## Clustering
  make_clustersArgs <- list(
    x = e,
    channels = primary_channels,
    SUFFIX = "_cluster-id", clear_1_cache = FALSE
  )
  make_clustersArgs <-
    utils::modifyList(make_clustersArgs, make_clusters_dots, keep.null = TRUE)
  make_clustersArgs <-
    utils::modifyList(make_clustersArgs, make_metaclusters..., keep.null = TRUE)

  do_stop_check(make_clustersArgs, stop_var = "STOP")
  cluster_id <- do.call(make_clusters, make_clustersArgs)
  if (!is_invalid(storage_env)) assign("cluster_id", cluster_id, envir = storage_env)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(make_clusters...)

  ## Which set has more clusters? Assume that more is better:
  clusterDiff <- dplyr::n_distinct(cluster_id) - dplyr::n_distinct(cluster_id_mc)
  if (clusterDiff >= 0)
    attr(e, "cluster_id") <- cluster_id
  else
    attr(e, "cluster_id") <- cluster_id_mc


  ## Analysis channels
  e_channels_abo <- structure(colnames(e), .Names = colnames(e))
  analysis_channels <- channel_sheet %>% dplyr::filter(use == "X") %>%
    { structure(dplyr::pull(., description), .Names = dplyr::pull(., name)) }
  if (!is_invalid(storage_env)) assign("analysis_channels", analysis_channels, envir = storage_env)

  ## Rename columns of combined events to match markers
  colnames(e) <- e_channels_abo %>% `[<-`(names(analysis_channels), analysis_channels) %>%
    as.vector
  if (!is_invalid(storage_env)) assign("e", e, envir = storage_env)


  ### Clusters based on investigator-supplied cell-subtype definitions

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "cell_subtypes")
  i$cell_subtypes$stop <- NULL

  cluster_sets <- list(auto = attr(e, "cluster_id"))

  if (is_invalid(label_threshold))
    label_threshold <- 0.55
  cell_subtype_parallel_processing <- i$cell_subtypes$parallel_processing
  i$cell_subtypes$parallel_processing <- NULL

  cell_subtypes <- i$cell_subtypes

  if (!is_invalid(cell_subtypes)) {
    options_keystone_parallel <- getOption("keystone_parallel")
    if (!is_invalid(cell_subtype_parallel_processing) && !cell_subtype_parallel_processing) {
      options(keystone_parallel = FALSE)
    }

    ## Find event clusters, i.e. manually gated events
    merge_clustersArgs <- list(
      x = e,#[, analysis_channels],
      clusters = cell_subtypes,
      channels  = analysis_channels,
      label_threshold = label_threshold, # Default is 0.90
      which_cluster_set = NULL, # Search every event for cluster definitions
      make_gating_poster = paste(image_dir, "gated-clusters-poster", sep = "/"),
      SUFFIX = "_gated-clusters-poster", clear_1_cache = clear_1_cache
    )
    merge_clustersArgs <-
      utils::modifyList(merge_clustersArgs, merge_clusters..., keep.null = TRUE)

    do_stop_check(merge_clustersArgs, stop_var = "STOP")
    gated_clusters <- do.call(merge_clusters, merge_clustersArgs)
    if (!is_invalid(storage_env)) assign("gated_clusters", gated_clusters, envir = storage_env)


    ## Search clusters for user-defined cell types & merge them if necessary
    merge_clustersArgs$which_cluster_set <- 1
    merge_clustersArgs$make_gating_poster <- NULL
    merge_clustersArgs$SUFFIX = "_merged-clusters"

    merged_clusters <- do.call(merge_clusters, merge_clustersArgs)
    if (!is_invalid(storage_env)) assign("merged_clusters", merge_clusters, envir = storage_env)

    ## Restore state of parallel processing to before this function call:
    options(keystone_parallel = options_keystone_parallel)

    ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
    clear_1_cache <- get_caching_status(make_clusters...)

    cluster_sets %<>%
      utils::modifyList(list(gated = gated_clusters$new_cluster_id,
        merged = merged_clusters$new_cluster_id), keep.null = FALSE)# %>% purrr::compact()
  }
  ## Make 'cluster_sets$gated' the automatic-cluster labels if there are no user-defined clusters:
  if (is_invalid(cluster_sets$gated))
    cluster_sets$gated <- cluster_sets$auto
  if (!is_invalid(storage_env)) assign("cluster_sets", cluster_sets, envir = storage_env)

  ####################
  ### Differential expression
  ####################

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "differential_expression")

  model_formula <- i$differential_expression$model_formula
  if (is_invalid(model_formula))
    model_formula <- "~ group"
  model_formula %<>% as.formula
  if (!is_invalid(storage_env)) assign("model_formula", model_formula, envir = storage_env)

  do_differential_expressionArgs <- list(
    x = e,
    m = metadata_i,
    cluster_set = cluster_sets$gated,
    #id_map_re = sample_name_re,
    model_formula = model_formula
  )
  do_differential_expressionArgs <-
    utils::modifyList(do_differential_expressionArgs, do_differential_expression...,
      keep.null = TRUE)

  fits <- do.call(do_differential_expression, do_differential_expressionArgs)
  if (!is_invalid(storage_env)) assign("fits", fits, envir = storage_env)
  ## Follow-up testing e.g.:
  #sapply(fits, edgeR::gof, simplify = FALSE) %>% print
  ## N.B. The 'coef' arg of 'edgeR::glmQLFTest()' defaults to the last regressor of the design matrix, so name/number it explicitly.
  #sapply(fits, function(a) { edgeR::glmQLFTest(a, coef = 2) %>% edgeR::topTags(Inf) %>% as.data.frame %>% dplyr::filter(FDR < 0.05) }, simplify = FALSE)

  interesting_contrasts <- i$differential_expression$interesting_contrasts
  if (is_invalid(interesting_contrasts))
    stop("No contrasts supplied for differential analysis")
  interesting_contrasts %<>% sapply(
    function(a)
    {
      a %<>% stringr::str_trim()
      if (stringr::str_detect(a, "^makeContrasts\\(")) {
        a %<>% stringr::str_replace("(makeContrasts\\()", "\\1a = ") %>%
          stringr::str_replace("(\\))$", ", levels = fit$design\\1") %>%
          parse(text = .)
      }

      a
    }, simplify = FALSE)

  alpha <- i$differential_expression$alpha
  if (is_invalid(alpha))
    alpha <- 0.05

  test_contrastsArgs <- list(
    fit = fits,
    ## These should be group var + factor level of interest, or contrasts:
    contrasts = interesting_contrasts,
    include_other_vars = FALSE,
    alpha = alpha
  )
  test_contrastsArgs <-
    utils::modifyList(test_contrastsArgs, test_contrasts..., keep.null = TRUE)

  inference <- do.call(test_contrasts, test_contrastsArgs)
  if (!is_invalid(storage_env)) assign("inference", inference, envir = storage_env)
  #sapply(inference, function(a) a$sig_results, simplify = FALSE)
  #sapply(inference, function(a) plyr::llply(a$res, function(b) b %>% edgeR::topTags(Inf) %>% as.data.frame), simplify = FALSE)

  if (!is_invalid(cluster_sets$merged)) {
    do_differential_expressionArgs$cluster_set <- cluster_sets$merged
    do_differential_expressionArgs <-
      utils::modifyList(do_differential_expressionArgs, do_differential_expression...,
        keep.null = TRUE)

    fitsm <- do.call(do_differential_expression, do_differential_expressionArgs)
    if (!is_invalid(storage_env)) assign("fitsm", fitsm, envir = storage_env)
    ## Follow-up testing e.g.:
    #sapply(fitsm, edgeR::gof, simplify = FALSE) %>% print
    ## N.B. The 'coef' arg of 'edgeR::glmQLFTest()' defaults to the last regressor of the design matrix, so name/number it explicitly.
    #sapply(fitsm, function(a) { edgeR::glmQLFTest(a, coef = 2) %>% edgeR::topTags(Inf) %>% as.data.frame %>% dplyr::filter(FDR < 0.05) }, simplify = FALSE)

    test_contrastsArgs$fit <- fitsm
    test_contrastsArgs <-
      utils::modifyList(test_contrastsArgs, test_contrasts..., keep.null = TRUE)

    inferencem <- do.call(test_contrasts, test_contrastsArgs)
    if (!is_invalid(storage_env)) assign("inferencem", inferencem, envir = storage_env)
    #sapply(inferencem, function(a) a$sig_results, simplify = FALSE)
    #sapply(inferencem, function(a) plyr::llply(a$res, function(b) b %>% edgeR::topTags(Inf) %>% as.data.frame), simplify = FALSE)
  }

  ####################
  ### Plots & reporting
  ####################

  ## UMAP; results in variable 'umap'.
  make_umap_embeddingArgs <- list(
    x = e[, analysis_channels],
    #seed = 667,
    n_threads =
      ifelse(future::availableCores() > 1L, trunc(future::availableCores()/2), 1L),
    fast_sgd = TRUE,
    #n_sgd_threads = "auto", batch = TRUE,
    SUFFIX = "_umap", clear_1_cache = clear_1_cache
  )
  make_umap_embeddingArgs <-
    utils::modifyList(make_umap_embeddingArgs, make_umap_embedding..., keep.null = TRUE)

  do_stop_check(make_umap_embeddingArgs, stop_var = "STOP")
  umap <- do.call(make_umap_embedding, make_umap_embeddingArgs)
  if (!is_invalid(storage_env)) assign("umap", umap, envir = storage_env)
  #set.seed(666); plot(umap, col = randomcoloR::distinctColorPalette(cluster_id %>% unique %>% length)[cluster_id])

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(make_umap_embedding...)


  ## Summarize cluster phenotypes (gated clusters)
  summarize_all_clustersArgs <- list(
    x = e[, analysis_channels],
    cluster_set = cluster_sets$gated,
    summary... = list(label_threshold = 0.55, collapse = ";"),
    callback = expression({
      make_external_latex_document(
        summarize_all_clusters_latex(sac, type = "table") %>% paste(collapse = "\n\n"),
        file_path = paste(report_dir, "cluster-phenotypes.tex", sep = "/")
      )
    }),
    SUFFIX = "_sac", clear_1_cache = clear_1_cache
  )
  summarize_all_clustersArgs <-
    utils::modifyList(summarize_all_clustersArgs, summarize_all_clusters..., keep.null = TRUE)

  do_stop_check(summarize_all_clustersArgs, stop_var = "STOP")
  sac <- do.call(summarize_all_clusters, summarize_all_clustersArgs, quote = TRUE)
  if (!is_invalid(storage_env)) assign("sac", sac, envir = storage_env)
  ## Then:
  # summarize_all_clusters_latex(sac, type = "table")

  ## Summarize cluster phenotypes (merged clusters)
  if (!is_invalid(cluster_sets$merged)) {
    summarize_all_clustersArgs$cluster_set <- cluster_sets$merged
    summarize_all_clustersArgs$callback <- expression({
      make_external_latex_document(
        summarize_all_clusters_latex(sac, type = "table") %>% paste(collapse = "\n\n"),
        file_path = paste(report_dir, "cluster-phenotypes-merged.tex", sep = "/")
      )
    })
    summarize_all_clustersArgs$SUFFIX <- "_sacm"

    summarize_all_clustersArgs <-
      utils::modifyList(summarize_all_clustersArgs, summarize_all_clusters..., keep.null = TRUE)

    sacm <- do.call(summarize_all_clusters, summarize_all_clustersArgs, quote = TRUE)
    if (!is_invalid(storage_env)) assign("sacm", sacm, envir = storage_env)
    ## Then:
    # summarize_all_clusters_latex(sacm, type = "table")
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(summarize_all_clusters...)


  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## UMAP plots (gated clusters)
  plot_common_umap_vizArgs <- list(
    x = e,
    cluster_set = cluster_sets$gated,
    channels = analysis_channels,
    m = metadata,
    umap = umap,
    #sample_name_re = sample_name_re,
    label_clusters = TRUE,
    image_dir = paste(image_dir, "gated-clusters-umap", sep = "/"),
    current_image = flowpipe_params("current_image_set"),
    save_plot = TRUE,
    #devices = flowpipe:::graphics_devices["grDevices::pdf"],
    use_complete_centroids = TRUE,
    SUFFIX = "_umap-viz_gated", clear_1_cache = clear_1_cache
  )
  plot_common_umap_vizArgs <-
    utils::modifyList(plot_common_umap_vizArgs, plot_common_umap_viz..., keep.null = TRUE)

  do_stop_check(plot_common_umap_vizArgs, stop_var = "STOP")
  do.call(plot_common_umap_viz, plot_common_umap_vizArgs)

  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## UMAP plots (merged clusters)
  if (!is_invalid(cluster_sets$merged)) {
    plot_common_umap_vizArgs$cluster_set <- cluster_sets$merged
    plot_common_umap_vizArgs$current_image <- flowpipe_params("current_image_set")
    plot_common_umap_vizArgs$image_dir <- paste(image_dir, "merged-clusters-umap", sep = "/")
    plot_common_umap_vizArgs$SUFFIX <- "_umap-viz-merged"

    plot_common_umap_vizArgs <-
      utils::modifyList(plot_common_umap_vizArgs, plot_common_umap_viz..., keep.null = TRUE)

    do.call(plot_common_umap_viz, plot_common_umap_vizArgs)
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(plot_common_umap_viz...)


  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## Plot heatmaps (gated clusters)
  plot_heatmapsArgs <- list(
    x = e[, analysis_channels],
    m = metadata,
    cluster_set = cluster_sets$gated,
    image_dir = paste(image_dir, "heatmaps", sep = "/"),
    current_image = flowpipe_params("current_image_set"),
    save_plot = TRUE,
    #devices = flowpipe:::graphics_devices["grDevices::pdf"],
    file_path_template = expression(sprintf("%s/%03d%s-%s", image_dir, current_image, "_heatmap_channels-gated-clusters",
      fs::path_sanitize(stringr::str_trunc(as.character(which_cluster_set), 31), "_"))),
    SUFFIX = "_medians-gated", clear_1_cache = clear_1_cache
  )
  plot_heatmapsArgs <-
    utils::modifyList(plot_heatmapsArgs, plot_heatmaps..., keep.null = TRUE)

  do_stop_check(plot_heatmapsArgs, stop_var = "STOP")
  cluster_median_matrices <- do.call(plot_heatmaps, plot_heatmapsArgs, quote = TRUE)
  if (!is_invalid(storage_env)) assign("cluster_median_matrices", cluster_median_matrices,
    envir = storage_env)

  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## Plot heatmaps (merged clusters)
  if (!is_invalid(cluster_sets$merged)) {
    plot_heatmapsArgs$cluster_set <- cluster_sets$merged
    plot_heatmapsArgs$current_image <- flowpipe_params("current_image_set")
    plot_heatmapsArgs$file_path_template <- expression(sprintf("%s/%03d%s-%s", image_dir,
      current_image, "_heatmap_channels-merged-clusters",
      fs::path_sanitize(stringr::str_trunc(as.character(which_cluster_set), 31), "_")))
    plot_heatmapsArgs$SUFFIX <- "_medians-merged"

    plot_heatmapsArgs <-
      utils::modifyList(plot_heatmapsArgs, plot_heatmaps..., keep.null = TRUE)

    cluster_median_matrices_merged <- do.call(plot_heatmaps, plot_heatmapsArgs, quote = TRUE)
    if (!is_invalid(storage_env)) assign("cluster_median_matrices_merged",
      cluster_median_matrices_merged, envir = storage_env)
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(plot_heatmaps...)


  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## Plot differential abundance (gated clusters)
  plot_differential_abundanceArgs <- list(
    x = e,
    m = metadata,
    umap = umap,
    fit = fits,
    contrasts = interesting_contrasts,
    cluster_set = cluster_sets$gated,
    alpha = alpha,
    #sample_name_re = sample_name_re,
    image_dir = paste(image_dir, "diff-expression-umap/gated", sep = "/"),
    current_image = flowpipe_params("current_image_set"),
    save_plot = TRUE,
    #devices = flowpipe:::graphics_devices["grDevices::pdf"],
    SUFFIX = "_diff-abundance-gated", clear_1_cache = clear_1_cache
  )
  plot_differential_abundanceArgs <-
    utils::modifyList(plot_differential_abundanceArgs, plot_differential_abundance...,
      keep.null = TRUE)

  do_stop_check(plot_differential_abundanceArgs, stop_var = "STOP")
  do.call(plot_differential_abundance, plot_differential_abundanceArgs)

  ## Bump up image number:
  flowpipe_params(current_image_set = flowpipe_params("current_image_set") + 1)

  ## Plot differential abundance (merged clusters)
  if (!is_invalid(cluster_sets$merged)) {
    plot_differential_abundanceArgs$cluster_set <- cluster_sets$merged
    plot_differential_abundanceArgs$image_dir <-
      paste(image_dir, "diff-expression-umap/merged", sep = "/")
    plot_differential_abundanceArgs$current_image <- flowpipe_params("current_image_set")
    plot_differential_abundanceArgs$SUFFIX <- "_diff-abundance-merged"

    plot_differential_abundanceArgs <-
      utils::modifyList(plot_differential_abundanceArgs, plot_differential_abundance...,
        keep.null = TRUE)

    do.call(plot_differential_abundance, plot_differential_abundanceArgs)
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(plot_differential_abundance...)


  details <- list(gated = list(fits, sac, cluster_median_matrices))
  if (!is_invalid(cluster_sets$merged)) {
    details %<>%
      utils::modifyList(list(merged = list(fitsm, sacm, cluster_median_matrices_merged)),
        keep.null = TRUE)
  }

  export_cluster_summaryArgs <- list(
    details = details,
    spreadsheet_path = paste(report_dir, "cluster-summary-tables.xlsx", sep = "/"),
    keep_pm_lists = FALSE,
    overwrite = TRUE,
    SUFFIX = "_cluster-summary", clear_1_cache = clear_1_cache
  )
  export_cluster_summaryArgs <-
    utils::modifyList(export_cluster_summaryArgs, export_cluster_summary..., keep.null = TRUE)

  do_stop_check(export_cluster_summaryArgs, stop_var = "STOP")
  cluster_summary_tables <- do.call(export_cluster_summary, export_cluster_summaryArgs)
  if (!is_invalid(storage_env)) assign("cluster_summary_tables",
    cluster_summary_tables, envir = storage_env)

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(export_cluster_summary...)

  ## Stop pipeline here? (Errors & data can be evaluated from the command line.)
  do_stop_check(ii, "cluster_fcs")

  l <- list(gated = cluster_sets$gated)
  if (!is_invalid(cluster_sets$merged)) {
    l %<>% utils::modifyList(list(merged = cluster_sets$merged), keep.null = TRUE)
  }
  cluster_fcs_dir <- i$cluster_fcs$directory
  if (!is_invalid(cluster_fcs_dir)) {
    ## Break up "pmm" matrices into FCS files by cluster
    lfs <- split_pmm_by_clusterArgs <- list(
      x = e,
      l = l,
      fcs_dir = cluster_fcs_dir, # Set above near version check
      export_id_map = TRUE,
      SUFFIX = "_cluster-fcs", clear_1_cache = clear_1_cache
    )

    do_stop_check(split_pmm_by_clusterArgs, stop_var = "STOP")
    split_pmm_by_clusterArgs <-
      utils::modifyList(split_pmm_by_clusterArgs, split_pmm_by_cluster..., keep.null = TRUE)

    lfs <- do.call(split_pmm_by_cluster, split_pmm_by_clusterArgs)
    if (!is_invalid(storage_env)) assign("lfs", lfs, envir = storage_env)
  }

  ## If 'clear_1_cache = TRUE' for this function, rerun all subsequent cached functions
  clear_1_cache <- get_caching_status(split_pmm_by_cluster...)


  #browser()
}


guess_channels <- function(
  cbs, # "cbs" object from calling 'get_channels_by_sample()'
  type, # cytometry type
  channels_regex,
  extract_marker_regex = ".*?_(.*)",
  usage_marker = "X"
)
{
  cre <- channels_regex[[type]]
  cs <- cbs %>% dplyr::select(name = 1, description = 2)

  ## Assume 'cre' of length > 1 is a vector of channel names, not a RegEx
  if (length(cre) > 1) {
    cre %<>% stringr::str_trim() %>% { stringr::str_flatten(rex::escape(.), "|") }
    analysis_channels_index <-
      sapply(cs, stringr::str_detect, pattern = cre, simplify = FALSE) %>%
        purrr::reduce(`|`) %>% tidyr::replace_na(FALSE)
  } else {
    analysis_channels_index <-
      cbs %>% { stringr::str_detect(attr(., "channels_by_sample_full_desc")$desc_01,
        stringr::regex(cre, ignore_case = TRUE), negate = TRUE) }
  }

  cs %<>%
    dplyr::mutate(
      description =
        dplyr::case_when(analysis_channels_index ~ stringr::str_extract(description,
          extract_marker_regex, group = 1), TRUE ~ description),
      use = NA_character_,
      use =
        dplyr::case_when(analysis_channels_index ~ usage_marker, TRUE ~ use),
    )

  cs
}


prepare_metadata <- function(
  x, # vector of file paths
  ... # additional data to add to the metadata table
)
{
  d <- purrr::reduce(
    list(
      #list(sample_id = seq_along(x), id = tools::file_path_sans_ext(basename(x))),
      list(id = tools::file_path_sans_ext(basename(x))),
      list(...)
    ),
    dplyr::bind_cols)

  as.data.frame(d, optional = TRUE)
}


do_stop_check <- function(
  toml, # List object returned by 'blogdown::read_toml()' or other list
  table_name, # Name of top-level table to search in
  stop_var = "stop"
)
{
  st <- NULL
  msg <- "No error: program halted by user request."

  if (!is_invalid(table_name)) {
    if (!is_invalid(toml[[table_name]])) {
      st <- toml[[table_name]][[stop_var]]
      msg <- sprintf("No error: program halted by user request at step \"%s\".", table_name)
    }
  } else {
    st <- toml[[stop_var]]
    if (!is_invalid(toml$SUFFIX))
      msg <-
        sprintf("No error: program halted by user request at function w/ output suffix \"%s\"",
          toml$SUFFIX)
  }

  if (!is_invalid(st) && is.logical(st) && st)
    keystone::stop_no_error(msg)

  nop()
}
