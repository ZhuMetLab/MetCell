.target_extract_individual_cell_data <- function(query_data_file, res_peak_df_single_file, mz_tol_ppm, mobility_tol,
                                                 res_define_at, tmp_dir, pool_size = 100L) {
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  res_peak_df_single <- readRDS(res_peak_df_single_file)
  # res_peak_df_single <- res_peak_df_single[res_peak_df_single$charge == 1,]
  res_peak_df_single <- res_peak_df_single
  query_data <- readRDS(query_data_file)

  all_mobility <- query_data$all_mobility
  cell_frame <- query_data$all_frames
  cell_scan_list <- lapply(cell_frame, function(dt) {
    # setDT(dt)
    collapse::rsplit(dt, dt$scan)
  })
  rm(cell_frame)
  gc()

  if (mz_tol_ppm > 1) {
    mz_tol <- .ppm2dalton2(res_peak_df_single$mz, mz_tol_ppm, res_define_at = res_define_at)
  } else {
    mz_tol <- mz_tol_ppm
  }
  res_peak_df_single$mz_tol <- mz_tol
  res_peak_df_single$mz_min_ <- res_peak_df_single$mz - mz_tol
  res_peak_df_single$mz_max_ <- res_peak_df_single$mz + mz_tol
  res_peak_df_single$mobility_min_ <- res_peak_df_single$mobility - mobility_tol
  res_peak_df_single$mobility_max_ <- res_peak_df_single$mobility + mobility_tol
  # res_peak_df_single_ <- res_peak_df_single[, 2:ncol(res_peak_df_single)]
  # mobility_scan <- names(all_mobility)


  return(list("tmp_data_files" = .split_data_list(cell_scan_list, tmp_dir, split_method = "average", pool_size = pool_size),
              'res_peak_df_single' = res_peak_df_single,
              "all_mobility" = all_mobility))
  }


.target_extract_individual_cell <- function(data_file,
                                            res_peak_df_single,
                                            all_mobility,
                                            res_define_at,
                                            ...){

  res_peak_df_single_ <- res_peak_df_single[, 2:ncol(res_peak_df_single)]
  mobility_scan <- names(all_mobility)
  data_file <- readRDS(data_file)
  # data_file <- data_file[1:10]
  re_extact_in_cell <- lapply(seq(nrow(res_peak_df_single)), function(ir) {
    rf_peak <- unlist(res_peak_df_single_[ir,])
    # re_extact_in_cell <- t(pbapply::pbapply(res_peak_df_single[, 2:ncol(res_peak_df_single), drop = FALSE], 1, function(rf_peak) {
    scan_index <- mobility_scan[findRangeIndicesDecreasing(all_mobility, c(rf_peak["mobility_max_"], rf_peak["mobility_min_"]))]
    sapply(data_file, function(cell_scan) {
      eim <- data.table::rbindlist(lapply(cell_scan[scan_index], function(sc) {
        if (is.null(sc)) {
          return(NULL)
        }
        res <- sc[findRangeIndices(sc$mz, c(rf_peak["mz_min_"], rf_peak["mz_max_"])), , drop = FALSE]
        if (nrow(res) > 1) {
          res <- res[which.min(abs(res$mz - rf_peak["mz"])), , drop = FALSE]
        }
        res
      }))
      ifelse(nrow(eim) == 0, NA, sum(eim$intensity))
    })
  })
  re_extact_in_cell <- do.call(rbind, re_extact_in_cell)

  return(re_extact_in_cell)
  }
















