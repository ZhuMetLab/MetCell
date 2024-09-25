#' @export
setMethod(
  "ExtractIndividualCell",
  signature = c("TimsData", "ExtractIndividualCellParam"),
  function(object, param) {
    # object <- tims_data
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Extracting peaks in individual cell...")

    # sub_dir <- file.path('group_peaks', 'density')
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'peak_in_cells')
    # files <- file.path(getwd(), object@experiment@tmp_dir, 'group_sc_peak', 'group_sc_peak_density')
    names(files) <- object@files


    if (.check_skip(files[1], param@rerun, show_message = FALSE)) {
      # browser()
      cat('  Using previouse results:\n    ')
      cat(paste0(files, collapse = '    \n    '))
      object@tmp_data_files$peak_in_cells <- files
      object@features <- readRDS(object@tmp_data_files$peak_in_cells)
      write.csv(object@features, file.path(object@experiment@res_dir, "01_feature_table.csv"), row.names = FALSE)
      object@tmp_data_files$detected_cell_number <- colnames(object@features)[grep(colnames(object@features), pattern = '@')]
      setwd(wd0)
      return(object)
    }
    ##### tmp data here #####
    param_list <- as.list(param)
    tmp <- .target_extract_individual_cell_data(object@tmp_data_files$tims_data_files[1], object@tmp_data_files$eim_peaks[1],
                                                param@mz_tol_ppm, param@mobility_tol, object@experiment@res_define_at,
                                                file.path(object@experiment@tmp_dir, 'tmp'), pool_size = param@pool_size)
    data_files <- tmp$tmp_data_files
    all_mobility <- tmp$all_mobility
    res_peak_df_single <- tmp$res_peak_df_single

    par_idx <- .gen_split_indexes(length(data_files), object@experiment@BPPARAM$workers)
    files_split <- .tmp_files(paste0('res', basename(data_files)), object@experiment@tmp_dir, 'peak_in_cells')



    # arg_list <- c(list("query_data_file" = unname(object@tmp_data_files$tims_data_files[1]),
    #                    # "sc_event_file" = unname(object@tmp_data_files$cell_superposition_frame),
    #                    "res_peak_df_single_file" = unname(object@tmp_data_files$eim_peaks[1]),
    #                    "res_define_at" = object@experiment@res_define_at),
    #               param_list)
    # peaks_in_cells <- .analysis_parser('.target_extract_individual_cell',
    #                                    arg_list,
    #                                    files)
    # browser()
    tmp_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(data_files[idxs], function(data_file) {
        # browser()
        list("data_file" = data_file,
          "res_peak_df_single" = data.frame(res_peak_df_single),
          "all_mobility" = as.array(all_mobility),
          "res_define_at" = object@experiment@res_define_at,
          "rerun" = param_list$rerun)
      })
      .parallel_parser(".target_extract_individual_cell", arg_list, files_split[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))

    peaks_in_cells <- do.call(cbind,lapply(tmp_files, readRDS))
    peaks_in_cells <- peaks_in_cells[, order(as.numeric(colnames(peaks_in_cells)), decreasing = FALSE)]
    is_keep <- apply(peaks_in_cells, 1, function(dr) { sum(!is.na(dr)) }) > param@minimal_cell_number
    peaks_in_cells <- peaks_in_cells[is_keep, , drop = FALSE]


    colnames(peaks_in_cells) <- paste0('cell#', colnames(peaks_in_cells), "@",
                                          basename(object@tmp_data_files$tims_data_files[1]))
    peaks_in_cells <- cbind(res_peak_df_single[, 1:10][is_keep, , drop = FALSE], peaks_in_cells)
    saveRDS(peaks_in_cells, file = files, version = 2)

    object@tmp_data_files$peaks_in_cells <- files
    object@features <- readRDS(object@tmp_data_files$peaks_in_cells)
    write.csv(object@features, file.path(object@experiment@res_dir, "01_feature_table.csv"), row.names = FALSE)
    object@tmp_data_files$detected_cell_number <- colnames(object@features)[grep(colnames(object@features), pattern = '@')]

    file.remove(data_files)
    file.remove(tmp_files)

    setwd(wd0)
    return(object)
  })
