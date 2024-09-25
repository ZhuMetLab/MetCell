#### .query_sctims_data_segment ####
.query_sctims_data_segment <- function(
  data_file,
  time_limit_table,
  marker_mz,
  marker_mobility,
  mz_tolerance_ppm,
  mobility_range,
  intensity_abs_threshold_upper,
  intensity_abs_threshold_lower,
  marker_eic_peak_span,
  res_define_at = 200,
  pool_size = 2000L,
  overlaps = 2L,
  opentims_thread = 2L,
  ...
) {
  opentimsr::opentims_set_threads(opentims_thread)

  all_columns <- .set_opentims()
  D <- opentimsr::OpenTIMS(data_file)

  time_limit_table_load <- read.csv(time_limit_table,
                                    stringsAsFactors = FALSE)
  idx <- match(basename(data_file),
               time_limit_table_load$data_file)

  start_time <- time_limit_table_load$start_time[idx]
  end_time <- time_limit_table_load$end_time[idx]

  # precursor_info <- opentimsr::table2df(D, 'Precursors')$Precursors[, c("Id", "Parent")]
  precursor_info <- NA
  if ('Precursors' %in% opentimsr::tables_names(D)) {
    precursor_info <- opentimsr::table2df(D, 'Precursors')$Precursors[, c("Id", "Parent")]
  }
  frame_id <- opentimsr::table2df(D, 'Frames')$Frames
  ms1_frame_info <- frame_id[frame_id$MsMsType == 0, c("Id", "Time")]
  ms1_frame_info <- ms1_frame_info[ms1_frame_info$Time >= start_time &
                                     ms1_frame_info$Time <= end_time,]
  # browser()

  all_mobility <- opentimsr::get_inv_ion_mobilities(D, ms1_frame_info$Id[1])

  ids <- ms1_frame_info$Id
  frame_idxs <- .gen_split_indexes(nrow(ms1_frame_info), pool_size, overlaps = overlaps)

  # browser()

  # all_frames <- unlist(BiocParallel::bplapply(1:nrow(frame_idxs), function(irow) {
  all_frames <- unlist(pbapply::pbapply(frame_idxs, 1, function(dr) {
    dt <- opentimsr::query(D, frames = ids[dr[1]:dr[2]], columns = c("frame", "scan", "mz", "intensity"))
    frames <- collapse::rsplit(dt, dt$frame, cols = 2:4)

    .discover_sc_events(frames,
                        all_mobility,
                        marker_mz,
                        marker_mobility,
                        mz_tolerance_ppm,
                        mobility_range,
                        intensity_abs_threshold_upper,
                        intensity_abs_threshold_lower,
                        marker_eic_peak_span,
                        res_define_at)
  }), recursive = FALSE)
  # },BPPARAM = BiocParallel::bpparam()), recursive = FALSE)

  # all_frames <- pbapply::pblapply(ms1_frame_info$Id, function(frame) {
  #   opentimsr::query(D, frames = frame, columns = c("scan", "mz", "intensity"))
  # })
  # names(all_frames) <- ms1_frame_info$Id
  # all_mobility <- .query_scan_mobility(D, ms1_frame_info$Id)

  res <- list("ms1_frame_info" = ms1_frame_info,
              "all_frames" = all_frames,
              "all_mobility" = all_mobility,
              "precursor_info" = precursor_info)
  opentimsr::CloseTIMS(D)
  return(res)
}


#### .trim_tims_data ####
.trim_tims_data <- function(
  tims_data_file,
  time_limit_table,
  ...
) {
  query_data <- readRDS(tims_data_file)
  time_limit_table_load <- read.csv(time_limit_table,
                                    stringsAsFactors = FALSE)
  idx <- match(basename(tims_data_file),
               time_limit_table_load$data_file)

  start_time <- time_limit_table_load$start_time[idx]
  end_time <- time_limit_table_load$end_time[idx]

  ms1_frame_info <- query_data$ms1_frame_info[query_data$ms1_frame_info$Time >= start_time &
                                                query_data$ms1_frame_info$Time <= end_time,]
  all_frames <- query_data$all_frames[names(query_data$all_frames) %in% ms1_frame_info$Id]

  all_mobility <- query_data$all_mobility
  precursor_info <- query_data$precursor_info


  res <- list("ms1_frame_info" = ms1_frame_info,
              "all_frames" = all_frames,
              "all_mobility" = all_mobility,
              "precursor_info" = precursor_info)
  return(res)
}

#### .find_sc_events ####
# marker_mz = 760.5730
# marker_mobility = 1.38
# mz_tolerance_ppm = 20
# mobility_range = 0.05
# res_define_at = 200
# intensity_abs_threshold_upper = 600000
# intensity_abs_threshold_lower = 40000
# trim_sctims_data_files <- "results/tmp/trim_tims_data/01_J_4e5_1016_10M.d"
.discover_sc_events <- function(
  data,
  all_mobility,
  marker_mz,
  marker_mobility,
  mz_tolerance_ppm,
  mobility_range,
  intensity_abs_threshold_upper,
  intensity_abs_threshold_lower,
  marker_eic_peak_span,
  res_define_at = 200,
  ...
) {
  # integrate intensity of the marker in each frame #
  mz_tol <- .ppm2dalton(marker_mz, mz_tolerance_ppm, res_define_at)

  temp_intensity <- lapply(seq_along(data), function(i) {
    temp_ions <- .get_eims2(data[i],
                            all_mobility, marker_mz, mz_tol,
                            marker_mobility, mobility_range)
    data.frame(frame_index = names(data)[i],
               marker_intensity = sum(temp_ions[, 1]))
  })

  temp_intensity <- do.call(rbind, temp_intensity)

  # idx <- which(temp_intensity$marker_intensity >= intensity_abs_threshold_lower &
  #                temp_intensity$marker_intensity <= intensity_abs_threshold_upper)

  temp_intensity$sc <- FALSE
  temp_intensity$sc[ggpmisc:::find_peaks(temp_intensity$marker_intensity, span = marker_eic_peak_span) &
                      temp_intensity$marker_intensity >= intensity_abs_threshold_lower &
                      temp_intensity$marker_intensity <= intensity_abs_threshold_upper] <- TRUE
  temp_intensity$sc[1] <- FALSE
  temp_intensity$sc[nrow(temp_intensity)] <- FALSE
  ### check the local max of the single-cell event marker ###

  # browser()

  # temp_intensity$sc <- sapply(seq(nrow(temp_intensity)), function(i) {
  #   if (temp_intensity$sc[i]) {
  #     if (temp_intensity$marker_intensity[i] > signal_fold_change * temp_intensity$marker_intensity[i + 1] &
  #       temp_intensity$marker_intensity[i] > signal_fold_change * temp_intensity$marker_intensity[i - 1]) {
  #       return(TRUE)
  #     }else {
  #       return(FALSE)
  #     }
  #   }else {
  #     return(FALSE)
  #   }
  # })

  # table(temp_intensity$sc)
  return(data[temp_intensity$sc])
}


#### .unite_cell_superposition_frame #####
# trim_sctims_data_files = 'E:/04_single_cell/00_data_processing/20231129_k562/results/tmp/trim_tims_data/21E_K562_3.d'
# names(trim_sctims_data_files) <- '21E_K562_3.d'
# single_cell_events = 'E:/04_single_cell/00_data_processing/20231129_k562/results/tmp/single_cell_events/21E_K562_3.d'
# names(single_cell_events) <- '21E_K562_3.d'
# res_define_at = 200
# mz_tolerance_combine = 20
# n_thread = 9
.unite_cell_superposition_frame <- function(
  data_file,
  mz_tolerance_combine,
  res_define_at,
  ...
) {
  res <- lapply(readRDS(data_file), function(scan_data) {
    scan_data <- scan_data[order(mz)]
    if (mz_tolerance_combine > 1) {
      tol <- .ppm2dalton2(scan_data$mz, mz_tolerance_combine, res_define_at)
    } else {
      tol <- mz_tolerance_combine
    }
    dt_tol <- cbind(scan_data$mz - tol, scan_data$mz + tol)

    intensity_order <- order(scan_data$intensity, decreasing = TRUE)

    mz_bin_list <- list()
    for (ir in seq(nrow(scan_data))) {
      i <- intensity_order[ir]
      if (scan_data$intensity[i] == 0) {
        next
      }

      idx_in <- findRangeIndices(scan_data$mz, dt_tol[i,])
      idx_in <- idx_in[scan_data$intensity[idx_in] > 0]
      temp_interest <- scan_data[idx_in]
      data.table::set(temp_interest, j = "i", value = i)
      mz_bin_list <- append(mz_bin_list, list(temp_interest))
      data.table::set(scan_data, idx_in, "intensity", 0)
    }
    mz_bins <- data.table::rbindlist(mz_bin_list) %>%
      dplyr::group_by(i) %>%
      dplyr::summarise(scan = first(scan), mz = weighted.mean(mz, intensity), intensity = sum(intensity)) %>%
      dplyr::select(-1) %>%
      dplyr::arrange(mz) %>%
      setDT
    return(mz_bins)
  })
  return(res)
}

.unite_cell_superposition_data <- function(data_file, cell_number, tmp_dir, pool_size = 40L) {
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }
  query_data <- readRDS(data_file)

  cell_frame <- query_data$all_frames
  if(!is.null(cell_number)){
    max_cell <- min(cell_number, length(cell_frame))
  }else{
    max_cell <- length(cell_frame)
  }
  cell_frame <- cell_frame[1:max_cell]
  print(paste0('the cell number for cell superposition frame is ', max_cell))

  all_mobility <- query_data$all_mobility

  cell_scan_list <- lapply(cell_frame, function(dt) {
    collapse::rsplit(dt, dt$scan)
  })
  rm(cell_frame)
  rm(query_data)
  gc()

  all_scan <- unique(unlist(lapply(cell_scan_list, names)))
  all_scans <- pbapply::pbsapply(all_scan, function(sc) {
    data.table::rbindlist(lapply(cell_scan_list, function(scan_list) {
      if (sc %in% names(scan_list)) {
        scan_list[[sc]]
      }
    }))
  }, simplify = FALSE)
  return(list("tmp_data_files" = .split_data_list(all_scans, tmp_dir, split_method = "average", pool_size = pool_size),
              "all_mobility" = all_mobility))
}

.load_all_trim_data <- function(trim_sctims_data_files,
                                single_cell_events) {
  data_each_file <- lapply(names(trim_sctims_data_files), function(nm) {
    trim_data <- trim_sctims_data_files[nm]
    sc_events <- single_cell_events[nm]

    query_data <- readRDS(trim_data)
    sc_events <- readRDS(sc_events)
    sc_index <- as.character(sc_events$frame_index[which(sc_events$sc)])
    cell_frame <- query_data$all_frames[sc_index]
    return(cell_frame)
  })

  data_all_file <- do.call(c, data_each_file)
  names(data_all_file) <- seq_along(data_all_file)
  return(data_all_file)
}

#### search peak targets in the cell superposition frame ####
.search_peak_target <- function(cell_superposition_frame_file,
                                mz_tol_ppm = 20,
                                res_define_at = 200,
                                n_skip = 1,
                                peak_target_length = 15,
                                ...) {
  superposition_frame <- readRDS(cell_superposition_frame_file)
  # make sure all scans are sorted by mz
  if (!all(superposition_frame$all_frames$`1` %>%
             group_by(scan) %>%
             summarise(unique(sign(diff(mz)))) %>%
             select(2) == 1)) {
    stop("Please make sure scans are sorted by m/z!")
  }

  nr_cumsum <- c(0, cumsum(table(superposition_frame$all_frames$`1`$scan)) - 1)
  v <- superposition_frame$all_frames$`1`$intensity
  ov <- order(v, decreasing = TRUE) - 1
  mz <- superposition_frame$all_frames$`1`$mz
  mzd <- .ppm2dalton2(mz, mz_tol_ppm, res_define_at)
  si <- nr_cumsum

  peak_target_res <- find_peak_targets(v, mz, mzd, si, ov,
                                       num = peak_target_length, thr = 0, n_skip = n_skip)
  return(peak_target_res)
}

#### detect peaks in the peak targets ####
.detect_peak_cell_superposition <- function(cell_superposition_frame_file,
                                            res_pktg_file,
                                            smooth_window = 10,
                                            peak_span_eim_detection = 21,
                                            peak_span_eim_integration = 27,
                                            signal_sd_threshold = 0.35,
                                            single_charge_line_slope = 0.001,
                                            single_charge_line_intercept = 0.35,
                                            skip_invalid_peaks = TRUE,
                                            ...) {
  # load the data ####
  cell_superposition_frame <- readRDS(cell_superposition_frame_file)
  res_pktg <- readRDS(res_pktg_file)

  # create a template ####
  all_mobility <- cell_superposition_frame$all_mobility
  temp_plate <- rep(0, length(all_mobility))
  names(temp_plate) <- names(all_mobility)

  # calculate the boundary and the length ####
  scan_min <- 1
  scan_max <- length(all_mobility)
  peak_span_eim_integration_one_side <- as.integer((peak_span_eim_integration - 1) / 2)
  loess_control <- loess.control(cell = 0.3)
  # peak detection in the cell superposition frame  ####
  res_peak <- pbapply::pblapply(seq_along(res_pktg), function(i) {
    # res_peak <- lapply(seq_along(res_pktg), function(i) {
    # cat(i, '\t')
    # i <- 43
    pktg <- res_pktg[[i]]

    # smooth the peak ####
    temp_pktg <- temp_plate
    temp_pktg[seq(pktg$rg[1], pktg$rg[2])] <- pktg$ints
    temp_pktg_smooth <- .smooth_loess(data = temp_pktg,
                                      degree = 1,
                                      window = smooth_window,
                                      control = loess_control)
    temp_mz <- temp_plate
    temp_mz[seq(pktg$rg[1], pktg$rg[2])] <- pktg$mzs

    # apex ####
    apex <- .find_apex_full_range(data_smooth = temp_pktg_smooth,
                                  peak_target_range = pktg$rg,
                                  min_intensity = 0,
                                  span = peak_span_eim_detection,
                                  ignore_threshold = 0.03)
    if (is.null(apex)) {
      return(NULL)
    }

    # using the fixed peak boundary(27 points, and 13 points each side) ####
    peak_bd_1 <- max(apex - peak_span_eim_integration_one_side, scan_min)
    peak_bd_2 <- min(apex + peak_span_eim_integration_one_side, scan_max)

    # calculate the peak sd ####
    peak_sd <- .get_smooth_sd(temp_pktg[peak_bd_1:peak_bd_2],
                              temp_pktg_smooth[peak_bd_1:peak_bd_2])

    if (peak_sd > signal_sd_threshold && skip_invalid_peaks) {
      return(NULL)
    }

    # return the result ####
    new_mz <- .get_weighted_mz(mzs = temp_mz[peak_bd_1:peak_bd_2],
                               intensities = temp_pktg[peak_bd_1:peak_bd_2])
    new_mobility <- all_mobility[apex]
    names(new_mobility) <- NULL
    new_int <- sum(temp_pktg[peak_bd_1:peak_bd_2])

    return(data.frame(mz = new_mz,
                      mobility = new_mobility,
                      int = new_int))

  })

  # the final data is an data frame ####
  res_peak_df <- do.call(rbind, res_peak)

  # calculate the charge number for each ion ####
  res_peak_df$charge <- 1
  res_peak_df$predict_mobility <- single_charge_line_slope * res_peak_df$mz + single_charge_line_intercept
  idx_single_charge <- which(res_peak_df$mobility > res_peak_df$predict_mobility)
  res_peak_df$charge[-idx_single_charge] <- 2

  # name features ####
  res_peak_df$ccs <- .mobility2ccs(k0 = res_peak_df$mobility,
                                   mz = res_peak_df$mz,
                                   z = res_peak_df$charge)

  res_peak_df <- cbind(data.frame(name = .get_peak_name2(res_peak_df[, c("mz", "ccs")])),
                       res_peak_df)
  # browser()
  all_col <- colnames(res_peak_df)
  res_peak_df <- res_peak_df[, c(c('name', 'mz', 'ccs'), all_col[-which(all_col %in% c('name', 'mz', 'ccs'))])]
  # res_peak_df <- select(res_peak_df,
  #                       name, mz, ccs, mobility, everything())
  return(res_peak_df)
}

.find_apex_full_range <- function(data_smooth = NULL,
                                  peak_target_range,
                                  min_intensity = 0,
                                  span = 7,
                                  ignore_threshold = 0.01) {
  data_smooth <- round(data_smooth, 10)
  apex <- which(ggpmisc:::find_peaks(data_smooth, span = span, ignore_threshold = ignore_threshold))
  apex <- apex[apex >= peak_target_range[1] & apex <= peak_target_range[2]]

  if (length(apex) >= 1) {
    apex <- apex[which.max(data_smooth[apex])]
  }else {
    apex <- NULL
  }
  return(apex)
}


#### target extraction in individual cells and generate the final result ####
# .target_extract_individual_cell <- function(res_peak_df_single_file,
#                                             query_data_file,
#                                             mz_tol_ppm = 20,
#                                             res_define_at = 200,
#                                             mobility_tol = 0.04,
#                                             minimal_cell_number = 10,
#                                             ...) {
#   res_peak_df_single <- readRDS(res_peak_df_single_file)
#   res_peak_df_single <- res_peak_df_single[res_peak_df_single$charge == 1,]
#   query_data <- readRDS(query_data_file)
#
#   all_mobility <- query_data$all_mobility
#   cell_frame <- query_data$all_frames
#   cell_scan_list <- lapply(cell_frame, function(dt) {
#     # setDT(dt)
#     collapse::rsplit(dt, dt$scan)
#   })
#   rm(cell_frame)
#   gc()
#
#   if (mz_tol_ppm > 1) {
#     mz_tol <- .ppm2dalton2(res_peak_df_single$mz, mz_tol_ppm, res_define_at = res_define_at)
#   } else {
#     mz_tol <- mz_tol_ppm
#   }
#   res_peak_df_single$mz_tol <- mz_tol
#   res_peak_df_single$mz_min_ <- res_peak_df_single$mz - mz_tol
#   res_peak_df_single$mz_max_ <- res_peak_df_single$mz + mz_tol
#   res_peak_df_single$mobility_min_ <- res_peak_df_single$mobility - mobility_tol
#   res_peak_df_single$mobility_max_ <- res_peak_df_single$mobility + mobility_tol
#   res_peak_df_single_ <- res_peak_df_single[, 2:ncol(res_peak_df_single)]
#   mobility_scan <- names(all_mobility)
#   # re_extact_in_cell <- data.table::transpose(data.table::rbindlist(BiocParallel::bplapply(seq(nrow(res_peak_df_single)), function(ir) {
#   re_extact_in_cell <- do.call(rbind, BiocParallel::bplapply(seq(nrow(res_peak_df_single)), function(ir) {
#     rf_peak <- unlist(res_peak_df_single_[ir,])
#   # re_extact_in_cell <- t(pbapply::pbapply(res_peak_df_single[, 2:ncol(res_peak_df_single), drop = FALSE], 1, function(rf_peak) {
#     scan_index <- mobility_scan[findRangeIndicesDecreasing(all_mobility, c(rf_peak["mobility_max_"], rf_peak["mobility_min_"]))]
#     sapply(cell_scan_list, function(cell_scan) {
#       eim <- data.table::rbindlist(lapply(cell_scan[scan_index], function(sc) {
#         if (is.null(sc)) {
#           return(NULL)
#         }
#         res <- sc[findRangeIndices(sc$mz, c(rf_peak["mz_min_"], rf_peak["mz_max_"])), , drop = FALSE]
#         if (nrow(res) > 1) {
#           res <- res[which.min(abs(res$mz - rf_peak["mz"])), , drop = FALSE]
#         }
#         res
#       }))
#       ifelse(nrow(eim) == 0, NA, sum(eim$intensity))
#     })
#   }, BPPARAM = BiocParallel::bpparam()))
#
#   # re_extact_in_cell <- pbapply::pblapply(seq(nrow(res_peak_df_single))[1:10], function(i) {
#   #   # i <- 160
#   #   # cat(i, '\t')
#   #   # i <- 254
#   #   rf_peak <- res_peak_df_single[i,]
#   #   rf_peak_in_cell <- lapply(seq_along(cell_frame), function(fr) {
#   #
#   #     temp_eim <- .get_eims3(cell_frame[fr],
#   #                            query_data$all_mobility,
#   #                            mz = rf_peak$mz,
#   #                            mz_tol = rf_peak$mz_tol,
#   #                            mobility = rf_peak$mobility,
#   #                            mobility_tol = mobility_tol)
#   #
#   #     if (is.null(temp_eim)) {
#   #       return(NA)
#   #     }
#   #
#   #     return(sum(temp_eim$intensity))
#   #   })
#   #
#   #   rf_peak_in_cell <- do.call(rbind, rf_peak_in_cell)
#   #
#   #   return(rf_peak_in_cell)
#   # })
#
#   # idx <- sapply(re_extact_in_cell, function(k) {
#   #   sum(!is.na(k))
#   # })
#   #
#   # idx <- which(idx >= minimal_cell_number)
#   #
#   # xxx <- lapply(re_extact_in_cell[idx], function(ppp) {
#   #   t(cbind(ppp))
#   # })
#   #
#   # re_extact_in_cell_df <- do.call(rbind, xxx)
#   # re_extact_in_cell_df <- as.data.frame(re_extact_in_cell_df)
#   #
#   # colnames(re_extact_in_cell_df) <- paste0('cell#', names(cell_scan_list), "@",
#   #                                          basename(query_data_file))
#   is_keep <- apply(re_extact_in_cell, 1, function(dr) { sum(!is.na(dr)) }) > minimal_cell_number
#   re_extact_in_cell <- re_extact_in_cell[is_keep, , drop = FALSE]
#   colnames(re_extact_in_cell) <- paste0('cell#', names(cell_scan_list), "@",
#                                         basename(query_data_file))
#   final_table <- cbind(res_peak_df_single[, 1:10][is_keep, , drop = FALSE], re_extact_in_cell)
#
#   return(final_table)
#
# }


#### isotope annotation ####
.annotate_isotope <- function(temp_feature,
                              mz_tol_ppm = 20,
                              res_define_at = 200,
                              mobility_tol = 0.01,
                              isotope_delta = 1.003355,
                              isotope_max_num = 4,
                              isotope_int_ratio_check = TRUE,
                              isotope_int_ratio_cutoff = 500,
                              ...) {

  temp_feature$base_peak <- NA
  temp_feature$isotope <- NA

  for (i in seq(nrow(temp_feature))) {
    # i <- 1
    cat(i, '\t')
    this_ion <- temp_feature[i,]

    if (!is.na(this_ion$base_peak)) {
      next
    }

    this_ion$base_peak <- this_ion$name
    this_ion$isotope <- "[M]"
    mz_isotope_list <- this_ion$mz + seq(0, isotope_max_num - 1) * isotope_delta
    isotope_matrix <- .get_iso_mass_range(mz = mz_isotope_list,
                                          mz_tol_ppm = mz_tol_ppm,
                                          res_define_at = res_define_at)

    idx <- temp_feature$mz > this_ion$mz &
      abs(this_ion$mobility - temp_feature$mobility) < mobility_tol &
      is.na(temp_feature$base_peak)
    # temp_peak_list <- rbind(this_ion, temp_feature[idx, ])
    temp_peak_list <- temp_feature[idx,]

    temp_res_iso <- lapply(2:length(mz_isotope_list), function(p) {
      mz_min <- isotope_matrix[p, 1]
      mz_max <- isotope_matrix[p, 2]

      label <- paste0("[M+", p - 1, "]")

      idx_iso <- which(temp_peak_list$mz > mz_min & temp_peak_list$mz < mz_max)

      if (length(idx_iso) > 0) {
        res_iso <- temp_peak_list[idx_iso,]
        res_iso$base_peak <- this_ion$name
        res_iso$isotope <- label
        res_iso$theo_mz <- mz_isotope_list[p]
        return(res_iso)
      }else {
        return(NULL)
      }
    })

    this_ion$theo_mz <- this_ion$mz

    temp_res_iso <- rbind(this_ion,
                          do.call(rbind, temp_res_iso))

    temp_res_iso$mz_diff <- abs(temp_res_iso$mz - temp_res_iso$theo_mz)
    temp_res_iso$int_ratio <- temp_res_iso$int / temp_res_iso$int[1]

    #### check the intensity ratio ####
    simulate_iso_df <- .simulate_theo_isotope(this_ion$mz,
                                              isotope_max_num = isotope_max_num)
    isp_match <- match(temp_res_iso$isotope, simulate_iso_df$label)
    temp_res_iso$int_theo <- simulate_iso_df$int_theo[isp_match]
    temp_res_iso$int_error <- abs(temp_res_iso$int_ratio - temp_res_iso$int_theo) / temp_res_iso$int_theo * 100

    if (isotope_int_ratio_check) {
      temp_res_iso <- temp_res_iso[which(temp_res_iso$int_error <= isotope_int_ratio_cutoff),]
    }


    #### remove those replicate ####

    temp_res_iso <- temp_res_iso[order(temp_res_iso$mz, temp_res_iso$mz, decreasing = c(FALSE, FALSE)),]
    temp_res_iso <- temp_res_iso[!duplicated(temp_res_iso$isotope),]

    #### output the final data ####
    final_match <- match(temp_res_iso$name,
                         temp_feature$name)

    temp_feature$base_peak[final_match] <- temp_res_iso$base_peak
    temp_feature$isotope[final_match] <- temp_res_iso$isotope

  }

  return(temp_feature)
}


.get_iso_mass_range <- function(mz,
                                mz_tol_ppm = 10,
                                res_define_at = 500) {
  result <- sapply(mz, function(x) {
    if (x >= res_define_at) {
      x * (1 + c(-1, 1) * mz_tol_ppm * 1e-6)
    } else {
      temp1 <- x + res_define_at * c(-1, 1) * mz_tol_ppm * 1e-6
    }
  })
  t(result)
}

.simulate_theo_isotope <- function(mz,
                                   isotope_max_num = 4) {
  # mz = 760.5732
  # isotope_max_num = 4
  # calculate simulated carbon number with alkane (CnH2n+2)
  num_carbon <- (mz - 1.0078 * 2) %/% 14.0156
  simulate_alkane_formula <- paste0('C', num_carbon, 'H', 2 * num_carbon + 2)

  options(readr.num_columns = 0)

  mol <- Rdisop::getMolecule(simulate_alkane_formula, z = 1)
  simulate_isotope <- as.data.frame(t(Rdisop::getIsotope(mol, index = seq(isotope_max_num))
  )
  )
  colnames(simulate_isotope) <- c('mz', 'int_theo')

  if (isotope_max_num > 1) {
    label <- c('[M]', paste0("[M+", seq(isotope_max_num - 1), "]"))
  } else {
    label <- '[M]'
  }

  # simulate_isotope <- simulate_isotope %>%
  #   dplyr::mutate(label = label) %>%
  #   dplyr::mutate(int_theo = int_theo/int_theo[1])

  simulate_isotope$label <- label
  simulate_isotope$int_theo <- simulate_isotope$int_theo / simulate_isotope$int_theo[1]

  return(simulate_isotope)

}

#### match with the library using m/z and ccs #####
.match_metabolite_info <- function(exp_data,
                                   lib_data,
                                   adduct_table,
                                   toleranceCCS,
                                   tolerancemz,
                                   typeCCS,
                                   res_define_at,
                                   info_col,
                                   other_col) {
  # info_col <- colnames(object@features)[which(!colnames(object@features) %in% row.names(object@sample_groups))]
  # other_col <- colnames(object@features)[which(colnames(object@features) %in% row.names(object@sample_groups))]

  # colnames(lib_data)[7:11] <- c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H]-', '[M+HCOO]-')
  colnames(lib_data)[19:25] <- c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H2O+H]+', '[M-H]-', '[M+Na-2H]-', '[M+HCOO]-')

  adducts <- adduct_table$name

  df_for_match <- data.frame(id = rep(lib_data$id, length(adducts)),
                             exact_mass = rep(lib_data$monoisotope_mass, length(adducts)))
  df_for_match$adduct <- rep(adduct_table$name, each = nrow(lib_data))
  df_for_match$massdiff <- rep(adduct_table$massdiff, each = nrow(lib_data))
  idx <- match(adducts, colnames(lib_data))
  # browser()
  if(length(idx) == 1){
    df_for_match$ccs <- c(lib_data[, idx])
  }else{
    df_for_match$ccs <- do.call(c, lib_data[, idx])
  }

  df_for_match$mz <- df_for_match$exact_mass + df_for_match$massdiff


  if (tolerancemz > 1) {
    df_for_match$mz_tol <- sapply(df_for_match$mz, function(mz) .ppm2dalton(mz, tolerancemz, res_define_at))
  } else {
    df_for_match$mz_tol <- tolerancemz
  }

  if (typeCCS == 'percentage') {
    df_for_match$ccs_tol <- df_for_match$ccs * toleranceCCS / 100
  } else {
    df_for_match$ccs_tol <- toleranceCCS
  }

  df_for_match$mz_up <- df_for_match$mz + df_for_match$mz_tol
  df_for_match$mz_low <- df_for_match$mz - df_for_match$mz_tol

  df_for_match$ccs_up <- df_for_match$ccs + df_for_match$ccs_tol
  df_for_match$ccs_low <- df_for_match$ccs - df_for_match$ccs_tol


  res_for_each_feature <- lapply(seq(nrow(exp_data)), function(i) {
    temp_row <- exp_data[i,]

    if (!temp_row$isotope == "[M]") {
      final_idx <- NULL
    }else {
      idx_mz <- which(temp_row$mz > df_for_match$mz_low &
                        temp_row$mz < df_for_match$mz_up)
      idx_ccs <- which(temp_row$ccs > df_for_match$ccs_low &
                         temp_row$ccs < df_for_match$ccs_up)
      final_idx <- intersect(idx_mz,
                             idx_ccs)
    }

    if (length(final_idx) == 0) {
      return(data.frame(id = NA,
                        compound = NA,
                        adduct = NA,
                        mz_error = NA,
                        ccs_error = NA,
                        formula = NA,
                        smiles = NA,
                        raw_id = NA))
    }else {
      idx_id <- match(df_for_match$id[final_idx],
                      lib_data$id)
      return(data.frame(id = paste0(lib_data$id[idx_id], collapse = ';'),
                        compound = paste0(lib_data$name[idx_id], collapse = ';'),
                        adduct = paste0(df_for_match$adduct[final_idx], collapse = ';'),
                        mz_error = paste0(round(abs(df_for_match$mz[final_idx] - temp_row$mz) / df_for_match$mz[final_idx] * 10^6, digits = 0), collapse = ';'),
                        ccs_error = paste0(round(abs(df_for_match$ccs[final_idx] - temp_row$ccs) / df_for_match$ccs[final_idx] * 100, digits = 0), collapse = ';'),
                        formula = paste0(lib_data$formula[idx_id], collapse = ';'),
                        smiles = paste0(lib_data$smiles[idx_id], collapse = ';'),
                        raw_id = paste0(lib_data$raw_id[idx_id], collapse = ';')))


    }
  })
  res_for_each_feature <- do.call(rbind, res_for_each_feature)

  output_table <- cbind(exp_data[, info_col],
                        res_for_each_feature,
                        exp_data[, other_col])
  return(output_table)
}

###############################################################################
#### previous code and function #####
###############################################################################
# cell_superposition_frame = 'E:/04_single_cell/00_package/07_cell_superposition/results/tmp/cell_superposition_frame/cell_superposition_frame'
# mz_tol = 20
# res_define_at = 200
# mobility_tol = 0.075
# min_intensity = 30
# min_points = 15
# min_points_peak = 20
# n_skip = 1
# smooth_window = 10
# peak_span_eim = 11
# snthreshold = 3
# skip_invalid_peaks = TRUE


# .detect_reference_peaks <- function(cell_superposition_frame,
#                                     mz_tol = 20,
#                                     res_define_at = 200,
#                                     mobility_tol = 0.075,
#                                     min_intensity = 30,
#                                     min_points = 15,
#                                     min_points_peak = 20,
#                                     n_skip = 1 ,
#                                     smooth_window = 10,
#                                     peak_span_eim = 11,
#                                     snthreshold = 3,
#                                     skip_invalid_peaks = TRUE,
#                                     ...){
#   cell_superposition_frame <- readRDS(cell_superposition_frame)
#   temp_frame <- cell_superposition_frame$all_frames
#
#   temp_frame[[1]]$num_ion <- paste0('#', seq(nrow(temp_frame[[1]])))
#   temp_frame[[1]] <- temp_frame[[1]][order(temp_frame[[1]]$intensity, decreasing = TRUE), ]
#
#   if (mz_tol > 1) {
#     temp_frame[[1]]$mz_tol <- sapply(temp_frame[[1]]$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
#   } else {
#     temp_frame[[1]]$mz_tol <- mz_tol
#   }
#
#
#   # temp_frame[[1]]$mz_up <- temp_frame[[1]]$mz + temp_frame[[1]]$mz_tol
#   # temp_frame[[1]]$mz_low <- temp_frame[[1]]$mz - temp_frame[[1]]$mz_tol
#   #
#   # first_rm <- sapply(seq(nrow(temp_frame[[1]])), function(i){
#   #   # cat(i, '\t')
#   #   # length(which(temp_frame[[1]]$mz_up >= temp_frame[[1]]$mz[i] &
#   #   #                temp_frame[[1]]$mz_low <= temp_frame[[1]]$mz[i] &
#   #   #                abs(temp_frame[[1]]$scan - temp_frame[[1]]$scan[i] < 50))) >= min_points/2
#   #   length(which(temp_frame[[1]]$mz_up >= temp_frame[[1]]$mz[i] &
#   #                  temp_frame[[1]]$mz_low <= temp_frame[[1]]$mz[i])) >= min_points/2
#   # },
#   # simplify = TRUE)
#   # temp_frame[[1]] <- temp_frame[[1]][which(first_rm), ]
#
#   roi_list <- list()
#   while (nrow(temp_frame[[1]]) > 1) {
#     cat(nrow(temp_frame[[1]]), '\t')
#     this_ion <- temp_frame[[1]][1, ]
#     temp_eim <- .get_eims3(temp_frame, cell_superposition_frame$all_mobility,
#                            mz = this_ion$mz,
#                            mz_tol = this_ion$mz_tol,
#                            mobility = cell_superposition_frame$all_mobility[as.character(this_ion$scan)],
#                            mobility_tol = mobility_tol)
#     temp_roi <- .find_roi(temp_eim$intensity,
#                           min_intensity = min_intensity,
#                           min_points = min_points,
#                           n_skip = n_skip,
#                           ref = match(this_ion$num_ion, temp_eim$num_ion))
#     if(!is.null(temp_roi)){
#       for(t in 1:nrow(temp_roi)){
#         temp_eim$roi <- FALSE
#         temp_eim$roi[temp_roi[t, 1]: temp_roi[t, 2]] <- TRUE
#         roi_list[[length(roi_list) + 1]] <- temp_eim
#         rm_idx <- temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion[which(!is.na(temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion))]
#         idx <- match(rm_idx, temp_frame[[1]]$num_ion)
#         temp_frame[[1]] <- temp_frame[[1]][-idx, ]
#       }
#     }else{
#       temp_frame[[1]] <- temp_frame[[1]][-1, ]
#     }
#   }
#
#   res <- list()
#   res[[1]] <- roi_list
#   names(res) <- '1'
#   # browser()
#
#
#   res_peaks_in_all_cells <- lapply(seq_along(res), function(cell_index){
#     # cell_index <- 1
#     one_cell <- res[[cell_index]]
#     res_in_one_cell <- lapply(seq_along(one_cell), function(roi_index){
#       # roi_index <- 1
#       # cat(roi_index, '\t')
#       temp <- one_cell[[roi_index]]
#       temp$intensity_smooth <- .smooth_loess(data = temp$intensity,
#                                              degree = 1,
#                                              window = smooth_window)
#
#
#       apex_eim <- temp$k0[.find_apex(temp$intensity,
#                                      temp$intensity_smooth,
#                                      span = peak_span_eim,
#                                      min_points = min_points_peak,
#                                      min_intensity = min_intensity,
#                                      n_skip = n_skip,
#                                      ref_index = which.max(temp$intensity_smooth),
#                                      find_roi = FALSE)]
#
#       if (length(apex_eim) == 0) {
#         return(NULL)
#       }
#
#       peak_bd <- temp$k0[.findLocalMin(temp$intensity_smooth)]
#
#       peak_bd_1 <- peak_bd[which(apex_eim > peak_bd)]
#       if(length(peak_bd_1) == 1){
#         peak_bd_1 <- peak_bd_1
#       } else if(length(peak_bd_1) == 0){
#         peak_bd_1 <- apex_eim - 0.015
#       } else if(length(peak_bd_1) > 1){
#         peak_bd_1 <- peak_bd_1[which.min(abs(apex_eim - peak_bd_1 - 0.015))]
#       }
#
#       peak_bd_2 <- peak_bd[which(apex_eim < peak_bd)]
#       if(length(peak_bd_2) == 1){
#         peak_bd_2 <- peak_bd_2
#       } else if(length(peak_bd_2) == 0){
#         peak_bd_2 <- apex_eim + 0.015
#       } else if(length(peak_bd_2) > 1){
#         peak_bd_2 <- peak_bd_2[which.min(abs(peak_bd_2 - apex_eim - 0.015))]
#       }
#
#       peak_range <- temp[temp$k0 >= peak_bd_1 & temp$k0 <= peak_bd_2, ]
#       peak_range <- peak_range[!is.na(peak_range$mz),]
#       new_mz <- .get_weighted_mz(peak_range$mz, peak_range$intensity)
#       peak_intensity <- sum(peak_range$intensity)
#
#
#       peak_sd <- .get_smooth_sd(temp$intensity, temp$intensity_smooth)
#
#       peak_quality <-.get_peak_quality(temp$intensity,
#                                       which.max(temp$intensity_smooth),
#                                       min_intensity = min_intensity,
#                                       min_points = min_points,
#                                       n_skip = n_skip,
#                                       snthreshold = snthreshold,
#                                       skip_invalid_peaks = skip_invalid_peaks)
#
#
#       if(peak_sd > 0.35 && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       if(is.null(peak_quality) && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       final_df <- data.frame(mobility = apex_eim,
#                              mobility_min = peak_bd_1,
#                              mobility_max = peak_bd_2,
#                              mz = new_mz,
#                              into = peak_intensity,
#                              sn = peak_quality$snr,
#                              baseline = peak_quality$baseline)
#       return(final_df)
#
#     })
#     return(res_in_one_cell)
#   })
#
#   names(res_peaks_in_all_cells) <- '1'
#
#   res_all_cell_peaks_df <- lapply(seq_along(res_peaks_in_all_cells), function(ccc){
#     # browser()
#     one_cell <- res_peaks_in_all_cells[[ccc]]
#     file_name <- 'cell_superposition_frame'
#     cell_name <- paste0(file_name, '@', names(res_peaks_in_all_cells)[ccc])
#     res_for_peaks <- lapply(one_cell, function(i){
#       if(!is.null(i)){
#         final_df <- i
#         final_df$sample <- cell_name
#         final_df$raw_file <- file_name
#         return(final_df)
#       }else{
#         return(NULL)
#       }
#     })
#     res_for_peaks <- do.call(rbind, res_for_peaks)
#     return(res_for_peaks)
#   })
#   res_all_cell_peaks_df <- do.call(rbind, res_all_cell_peaks_df)
#   row.names(res_all_cell_peaks_df) <- paste0('#', seq(nrow(res_all_cell_peaks_df)))
#   return(res_all_cell_peaks_df)
# }


# .detect_reference_peaks <- function(cell_superposition_frame,
#                                     mz_tol = 20,
#                                     res_define_at = 200,
#                                     mobility_tol = 0.075,
#                                     min_intensity = 30,
#                                     min_points = 15,
#                                     min_points_peak = 20,
#                                     n_skip = 1 ,
#                                     smooth_window = 10,
#                                     peak_span_eim = 11,
#                                     snthreshold = 3,
#                                     skip_invalid_peaks = TRUE,
#                                     ...){
#   cell_superposition_frame = 'E:/04_single_cell/00_package/17_1020cell_superposition/results/tmp/cell_superposition_frame/cell_superposition_frame'
#   mz_tol = 20
#   res_define_at = 200
#   mobility_tol = 0.075
#   min_intensity = 30
#   min_points = 15
#   min_points_peak = 20
#   n_skip = 1
#   smooth_window = 10
#   peak_span_eim = 11
#   snthreshold = 3
#   skip_invalid_peaks = TRUE
#
#   cell_superposition_frame <- readRDS(cell_superposition_frame)
#   temp_frame <- cell_superposition_frame$all_frames
#
#   temp_frame[[1]]$num_ion <- paste0('#', seq(nrow(temp_frame[[1]])))
#   row.names(temp_frame[[1]]) <- temp_frame[[1]]$num_ion
#   temp_frame[[1]] <- temp_frame[[1]][order(temp_frame[[1]]$intensity, decreasing = TRUE), ]
#
#   if (mz_tol > 1) {
#     temp_frame[[1]]$mz_tol <- sapply(temp_frame[[1]]$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
#   } else {
#     temp_frame[[1]]$mz_tol <- mz_tol
#   }
#
#   temp_frame_list <- split(temp_frame$`1`, temp_frame$`1`$scan)
#
#   roi_list <- list()
#   p_ion <- 1
#   # pb <- txtProgressBar(min = 1, max = nrow(temp_frame[[1]]), style = 3)
#   while (!all(temp_frame[[1]]$intensity == 0)) {
#     cat(p_ion, '\t')
#
#     if(!temp_frame[[1]]$intensity[p_ion] == 0){
#
#       this_ion <- temp_frame[[1]][p_ion, ]
#       temp_eim <- .get_eims3(temp_frame, cell_superposition_frame$all_mobility,
#                              mz = this_ion$mz,
#                              mz_tol = this_ion$mz_tol,
#                              mobility = cell_superposition_frame$all_mobility[as.character(this_ion$scan)],
#                              mobility_tol = mobility_tol)
#       temp_roi <- .find_roi(temp_eim$intensity,
#                             min_intensity = min_intensity,
#                             min_points = min_points,
#                             n_skip = n_skip,
#                             ref = match(this_ion$num_ion, temp_eim$num_ion))
#       if(!is.null(temp_roi)){
#         for(t in 1:nrow(temp_roi)){
#
#           roi_list[[length(roi_list) + 1]] <- temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion
#           temp_frame[[1]][na.omit(temp_eim$num_ion[temp_roi[t, 1]: temp_roi[t, 2]]), ]$intensity <- 0
#         }
#       }else{
#         temp_frame[[1]]$intensity[p_ion] <- 0
#       }
#     }
#     # setTxtProgressBar(pb, p_ion)
#     p_ion <- p_ion + 1
#
#   }
#   # close(pb)
#   res <- list()
#   res[[1]] <- roi_list
#   names(res) <- '1'
#   # browser()
#
#
#   res_peaks_in_all_cells <- lapply(seq_along(res), function(cell_index){
#     # cell_index <- 1
#     one_cell <- res[[cell_index]]
#     res_in_one_cell <- lapply(seq_along(one_cell), function(roi_index){
#       # roi_index <- 1
#       # cat(roi_index, '\t')
#       temp <- one_cell[[roi_index]]
#       temp$intensity_smooth <- .smooth_loess(data = temp$intensity,
#                                              degree = 1,
#                                              window = smooth_window)
#
#
#       apex_eim <- temp$k0[.find_apex(temp$intensity,
#                                      temp$intensity_smooth,
#                                      span = peak_span_eim,
#                                      min_points = min_points_peak,
#                                      min_intensity = min_intensity,
#                                      n_skip = n_skip,
#                                      ref_index = which.max(temp$intensity_smooth),
#                                      find_roi = FALSE)]
#
#       if (length(apex_eim) == 0) {
#         return(NULL)
#       }
#
#       peak_bd <- temp$k0[.findLocalMin(temp$intensity_smooth)]
#
#       peak_bd_1 <- peak_bd[which(apex_eim > peak_bd)]
#       if(length(peak_bd_1) == 1){
#         peak_bd_1 <- peak_bd_1
#       } else if(length(peak_bd_1) == 0){
#         peak_bd_1 <- apex_eim - 0.015
#       } else if(length(peak_bd_1) > 1){
#         peak_bd_1 <- peak_bd_1[which.min(abs(apex_eim - peak_bd_1 - 0.015))]
#       }
#
#       peak_bd_2 <- peak_bd[which(apex_eim < peak_bd)]
#       if(length(peak_bd_2) == 1){
#         peak_bd_2 <- peak_bd_2
#       } else if(length(peak_bd_2) == 0){
#         peak_bd_2 <- apex_eim + 0.015
#       } else if(length(peak_bd_2) > 1){
#         peak_bd_2 <- peak_bd_2[which.min(abs(peak_bd_2 - apex_eim - 0.015))]
#       }
#
#       peak_range <- temp[temp$k0 >= peak_bd_1 & temp$k0 <= peak_bd_2, ]
#       peak_range <- peak_range[!is.na(peak_range$mz),]
#       new_mz <- .get_weighted_mz(peak_range$mz, peak_range$intensity)
#       peak_intensity <- sum(peak_range$intensity)
#
#
#       peak_sd <- .get_smooth_sd(temp$intensity, temp$intensity_smooth)
#
#       peak_quality <-.get_peak_quality(temp$intensity,
#                                        which.max(temp$intensity_smooth),
#                                        min_intensity = min_intensity,
#                                        min_points = min_points,
#                                        n_skip = n_skip,
#                                        snthreshold = snthreshold,
#                                        skip_invalid_peaks = skip_invalid_peaks)
#
#
#       if(peak_sd > 0.35 && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       if(is.null(peak_quality) && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       final_df <- data.frame(mobility = apex_eim,
#                              mobility_min = peak_bd_1,
#                              mobility_max = peak_bd_2,
#                              mz = new_mz,
#                              into = peak_intensity,
#                              sn = peak_quality$snr,
#                              baseline = peak_quality$baseline)
#       return(final_df)
#
#     })
#     return(res_in_one_cell)
#   })
#
#   names(res_peaks_in_all_cells) <- '1'
#
#   res_all_cell_peaks_df <- lapply(seq_along(res_peaks_in_all_cells), function(ccc){
#     # browser()
#     one_cell <- res_peaks_in_all_cells[[ccc]]
#     file_name <- 'cell_superposition_frame'
#     cell_name <- paste0(file_name, '@', names(res_peaks_in_all_cells)[ccc])
#     res_for_peaks <- lapply(one_cell, function(i){
#       if(!is.null(i)){
#         final_df <- i
#         final_df$sample <- cell_name
#         final_df$raw_file <- file_name
#         return(final_df)
#       }else{
#         return(NULL)
#       }
#     })
#     res_for_peaks <- do.call(rbind, res_for_peaks)
#     return(res_for_peaks)
#   })
#   res_all_cell_peaks_df <- do.call(rbind, res_all_cell_peaks_df)
#   row.names(res_all_cell_peaks_df) <- paste0('#', seq(nrow(res_all_cell_peaks_df)))
#   return(res_all_cell_peaks_df)
# }
#### 20231024 slice the whole data frame to increase the speed of finding roi ####
# .detect_reference_peaks <- function(cell_superposition_frame,
#                                     mz_tol = 20,
#                                     res_define_at = 200,
#                                     mobility_tol = 0.075,
#                                     min_intensity = 30,
#                                     min_points = 15,
#                                     min_points_peak = 20,
#                                     n_skip = 1 ,
#                                     smooth_window = 10,
#                                     peak_span_eim = 11,
#                                     snthreshold = 3,
#                                     skip_invalid_peaks = TRUE,
#                                     ...){
#   cell_superposition_frame = 'E:/04_single_cell/00_package/17_1020cell_superposition/results/tmp/cell_superposition_frame/cell_superposition_frame'
#   mz_tol = 20
#   res_define_at = 200
#   mobility_tol = 0.075
#   min_intensity = 100
#   min_points = 15
#   min_points_peak = 20
#   n_skip = 0
#   smooth_window = 10
#   peak_span_eim = 11
#   snthreshold = 3
#   skip_invalid_peaks = TRUE
#
#   cell_superposition_frame <- readRDS(cell_superposition_frame)
#   # cell_superposition_frame$all_frames$`1` <- cell_superposition_frame$all_frames$`1`[which(cell_superposition_frame$all_frames$`1`$intensity > min_intensity), ]
#   temp_frame <- cell_superposition_frame$all_frames
#   all_mobility <- cell_superposition_frame$all_mobility
#
#   temp_frame[[1]]$num_ion <- paste0('#', seq(nrow(temp_frame[[1]])))
#   row.names(temp_frame[[1]]) <- temp_frame[[1]]$num_ion
#   # all_intensity <- temp_frame$`1`$intensity
#   # names(all_intensity) <- row.names(temp_frame[[1]])
#
#   temp_frame[[1]] <- temp_frame[[1]][order(temp_frame[[1]]$intensity, decreasing = TRUE), ]
#
#   if (mz_tol > 1) {
#     temp_frame[[1]]$mz_tol <- sapply(temp_frame[[1]]$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
#   } else {
#     temp_frame[[1]]$mz_tol <- mz_tol
#   }
#
#   #### slice the whole data frame into the list according to the scan #####
#   temp_frame_list <- split(temp_frame$`1`, temp_frame$`1`$scan)
#   roi_list <- list()
#
#   #### do loop for all the ions in the data frame ####
#   for(p_ion in seq(nrow(temp_frame[[1]]))){
#     # p_ion <- 1
#     cat(p_ion, '\t')
#     ###### if the intensity of the ion was 0, next #####
#     if(temp_frame[[1]]$intensity[p_ion] == 0){
#       next
#     }
#
#     ###### calculate the scan range ######
#     temp_mobility <- all_mobility[as.character(temp_frame[[1]]$scan[p_ion])]
#     temp_mz <- temp_frame[[1]]$mz[p_ion]
#     temp_mz_tol <- temp_frame[[1]]$mz_tol[p_ion]
#     target_scan <- names(all_mobility)[abs(all_mobility - temp_mobility) <= mobility_tol]
#
#     ####### extract the data and determine the roi #####
#     temp_eim <- .get_eims4(scan_list = temp_frame_list[target_scan],
#                           mz = temp_mz,
#                           mz_tol = temp_mz_tol)
#
#     temp_eim$intensity <- temp_frame$`1`[temp_eim$num_ion, ]$intensity
#     ###### switch the value of roi into 0 ######
#     temp_roi <- .find_roi(temp_eim$intensity,
#                           min_intensity = min_intensity,
#                           min_points = min_points,
#                           n_skip = n_skip,
#                           ref = match(names(temp_mobility), row.names(temp_eim)))
#     # temp_roi <- .find_roi(temp_eim$intensity,
#     #                       min_intensity = min_intensity,
#     #                       min_points = min_points,
#     #                       n_skip = n_skip,
#     #                       ref = NULL)
#
#     ###### output the ion number of the roi into the roi list #####
#     if(!is.null(temp_roi)){
#       for(t in 1:nrow(temp_roi)){
#         roi_list[[length(roi_list) + 1]] <- temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion
#         temp_frame[[1]][na.omit(temp_eim$num_ion[temp_roi[t, 1]: temp_roi[t, 2]]), ]$intensity <- 0
#       }
#     }else{
#       temp_frame[[1]]$intensity[p_ion] <- 0
#     }
#     # temp_frame[[1]]$intensity[p_ion] <- 0
#   }
#
#
#   res <- list()
#   res[[1]] <- roi_list
#   names(res) <- '1'
#   # browser()
#
#
#   res_peaks_in_all_cells <- lapply(seq_along(res), function(cell_index){
#     # cell_index <- 1
#     one_cell <- res[[cell_index]]
#     res_in_one_cell <- lapply(seq_along(one_cell), function(roi_index){
#       # roi_index <- 1
#       # cat(roi_index, '\t')
#       temp <- one_cell[[roi_index]]
#       temp$intensity_smooth <- .smooth_loess(data = temp$intensity,
#                                              degree = 1,
#                                              window = smooth_window)
#
#
#       apex_eim <- temp$k0[.find_apex(temp$intensity,
#                                      temp$intensity_smooth,
#                                      span = peak_span_eim,
#                                      min_points = min_points_peak,
#                                      min_intensity = min_intensity,
#                                      n_skip = n_skip,
#                                      ref_index = which.max(temp$intensity_smooth),
#                                      find_roi = FALSE)]
#
#       if (length(apex_eim) == 0) {
#         return(NULL)
#       }
#
#       peak_bd <- temp$k0[.findLocalMin(temp$intensity_smooth)]
#
#       peak_bd_1 <- peak_bd[which(apex_eim > peak_bd)]
#       if(length(peak_bd_1) == 1){
#         peak_bd_1 <- peak_bd_1
#       } else if(length(peak_bd_1) == 0){
#         peak_bd_1 <- apex_eim - 0.015
#       } else if(length(peak_bd_1) > 1){
#         peak_bd_1 <- peak_bd_1[which.min(abs(apex_eim - peak_bd_1 - 0.015))]
#       }
#
#       peak_bd_2 <- peak_bd[which(apex_eim < peak_bd)]
#       if(length(peak_bd_2) == 1){
#         peak_bd_2 <- peak_bd_2
#       } else if(length(peak_bd_2) == 0){
#         peak_bd_2 <- apex_eim + 0.015
#       } else if(length(peak_bd_2) > 1){
#         peak_bd_2 <- peak_bd_2[which.min(abs(peak_bd_2 - apex_eim - 0.015))]
#       }
#
#       peak_range <- temp[temp$k0 >= peak_bd_1 & temp$k0 <= peak_bd_2, ]
#       peak_range <- peak_range[!is.na(peak_range$mz),]
#       new_mz <- .get_weighted_mz(peak_range$mz, peak_range$intensity)
#       peak_intensity <- sum(peak_range$intensity)
#
#
#       peak_sd <- .get_smooth_sd(temp$intensity, temp$intensity_smooth)
#
#       peak_quality <-.get_peak_quality(temp$intensity,
#                                        which.max(temp$intensity_smooth),
#                                        min_intensity = min_intensity,
#                                        min_points = min_points,
#                                        n_skip = n_skip,
#                                        snthreshold = snthreshold,
#                                        skip_invalid_peaks = skip_invalid_peaks)
#
#
#       if(peak_sd > 0.35 && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       if(is.null(peak_quality) && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       final_df <- data.frame(mobility = apex_eim,
#                              mobility_min = peak_bd_1,
#                              mobility_max = peak_bd_2,
#                              mz = new_mz,
#                              into = peak_intensity,
#                              sn = peak_quality$snr,
#                              baseline = peak_quality$baseline)
#       return(final_df)
#
#     })
#     return(res_in_one_cell)
#   })
#
#   names(res_peaks_in_all_cells) <- '1'
#
#   res_all_cell_peaks_df <- lapply(seq_along(res_peaks_in_all_cells), function(ccc){
#     # browser()
#     one_cell <- res_peaks_in_all_cells[[ccc]]
#     file_name <- 'cell_superposition_frame'
#     cell_name <- paste0(file_name, '@', names(res_peaks_in_all_cells)[ccc])
#     res_for_peaks <- lapply(one_cell, function(i){
#       if(!is.null(i)){
#         final_df <- i
#         final_df$sample <- cell_name
#         final_df$raw_file <- file_name
#         return(final_df)
#       }else{
#         return(NULL)
#       }
#     })
#     res_for_peaks <- do.call(rbind, res_for_peaks)
#     return(res_for_peaks)
#   })
#   res_all_cell_peaks_df <- do.call(rbind, res_all_cell_peaks_df)
#   row.names(res_all_cell_peaks_df) <- paste0('#', seq(nrow(res_all_cell_peaks_df)))
#   return(res_all_cell_peaks_df)
# }
#


# trim_sctims_data_files = 'E:/04_single_cell/00_package/07_cell_superposition/results/tmp/trim_tims_data/metabolomics_3.d'
# single_cell_events = 'E:/04_single_cell/00_package/07_cell_superposition/results/tmp/single_cell_events/metabolomics_3.d'
# reference_peaks = 'E:/04_single_cell/00_package/07_cell_superposition/results/tmp/reference_peaks/reference_peaks'
# res_define_at = 200
# mz_tol = 20
# mobility_tol = 0.05
# mobility_tol_ref = 0.01
# mobility_tol_int = 0.01
# smooth_window = 10
# peak_span_eim = 11
# min_points = 5
# snthreshold = 3
# n_skip = 1
# min_intensity = 0
# skip_invalid_peaks = TRUE

# .extract_cell_peaks <- function(trim_sctims_data_files,
#                                 single_cell_events,
#                                 reference_peaks,
#                                 res_define_at,
#                                 mz_tol = 20,
#                                 mobility_tol = 0.05,
#                                 mobility_tol_ref = 0.01,
#                                 mobility_tol_int = 0.01,
#                                 smooth_window = 10,
#                                 peak_span_eim = 11,
#                                 min_points = 5,
#                                 snthreshold = 3,
#                                 n_skip = 1,
#                                 min_intensity = 0,
#                                 skip_invalid_peaks = TRUE,
#                                 ...
#                                 ){
#   res_all_cell_peaks_df <- readRDS(reference_peaks)
#   query_data <- readRDS(trim_sctims_data_files)
#   sc_event <- readRDS(single_cell_events)
#   sc_index <- as.character(sc_event$frame_index[which(sc_event$sc)])
#   all_mobility <- query_data$all_mobility
#   cell_frame <- query_data$all_frames[sc_index]
#
#   if (mz_tol > 1) {
#     res_all_cell_peaks_df$mz_tol <- sapply(res_all_cell_peaks_df$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
#   } else {
#     res_all_cell_peaks_df$mz_tol <- mz_tol
#   }
#
#   re_extact_in_cell <- lapply(seq(nrow(res_all_cell_peaks_df)), function(i){
#     # i <- 1
#     # cat(i, '\t')
#     # i <- 254
#     rf_peak <- res_all_cell_peaks_df[i, ]
#     rf_peak_in_cell <- lapply(seq_along(cell_frame), function(fr){
#       # fr <- 14
#       # cat(fr, '\t')
#       temp_eim <- .get_eims3(cell_frame[fr],
#                              query_data$all_mobility,
#                              mz = rf_peak$mz,
#                              mz_tol = rf_peak$mz_tol,
#                              mobility = rf_peak$mobility,
#                              mobility_tol = mobility_tol)
#       if(is.null(temp_eim)){
#         return(NULL)
#       }
#
#       temp_eim$intensity_smooth <- .smooth_loess(data = temp_eim$intensity,
#                                                  degree = 1,
#                                                  window = smooth_window)
#       # apex_idx <- .find_apex(temp_eim$intensity,
#       #                        temp_eim$intensity_smooth,
#       #                        span = peak_span_eim,
#       #                        min_points = min_points,
#       #                        min_intensity = min_intensity,
#       #                        n_skip = n_skip,
#       #                        ref_index = which.min(abs(rf_peak$mobility - temp_eim$k0)),
#       #                        find_roi = TRUE)
#       # apex_eim <- temp_eim$k0[apex_idx]
#
#       # if (length(apex_eim) == 0) {
#       #   return(NULL)
#       # }
#       #
#       # if(abs(apex_eim - rf_peak$mobility) > mobility_tol_ref){
#       #   return(NULL)
#       # }
#       apex_idx <- which.min(abs(rf_peak$mobility - temp_eim$k0))
#       apex_eim <- temp_eim$k0[apex_idx]
#       peak_sd <- .get_smooth_sd(temp_eim$intensity, temp_eim$intensity_smooth)
#
#       peak_quanlity <-.get_peak_quality(temp_eim$intensity,
#                                         apex_idx,
#                                         min_intensity = min_intensity,
#                                         min_points = min_points,
#                                         n_skip = n_skip,
#                                         snthreshold = snthreshold,
#                                         skip_invalid_peaks = skip_invalid_peaks)
#
#
#       if(peak_sd > 0.35 && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       if(is.null(peak_quanlity) && skip_invalid_peaks){
#         return(NULL)
#       }
#
#       peak_range <- temp_eim[which(temp_eim$k0 > apex_eim - mobility_tol_int & temp_eim$k0 < apex_eim + mobility_tol_int),
#       ]
#       peak_range <- peak_range[!is.na(peak_range$mz),]
#       # range_idx <- which(temp_eim$k0 > apex_eim - mobility_tol_int & temp_eim$k0 < apex_eim + mobility_tol_int)
#       peak_intensity <- sum(peak_range$intensity)
#       # range_idx <- range_idx
#       final_mz <- .get_weighted_mz(mzs = peak_range$mz,
#                                    intensities = peak_range$intensity)
#
#       final_df <- data.frame(mobility = apex_eim,
#                              mobility_min = apex_eim - mobility_tol_int,
#                              mobility_max = apex_eim + mobility_tol_int,
#                              mz = final_mz,
#                              into = peak_intensity,
#                              sn = peak_quanlity$snr,
#                              baseline = peak_quanlity$baseline,
#                              sample = paste0(basename(trim_sctims_data_files), '@', sc_index[fr]),
#                              raw_file = basename(trim_sctims_data_files),
#                              peak_groups = row.names(rf_peak))
#
#       return(final_df)
#     })
#
#     rf_peak_in_cell <- do.call(rbind, rf_peak_in_cell)
#
#     return(rf_peak_in_cell)
#   })
#
#   re_extact_in_cell_df <- do.call(rbind, re_extact_in_cell)
#   return(re_extact_in_cell_df)
# }
# peaks <- tims_data@peaks
# sampleGroups <- sample_groups
# minSamples = 1
# minFraction = 0.1
# .generate_sc_feature <- function(peaks,
#                                  sampleGroups,
#                                  minSamples = 1,
#                                  minFraction = 0.1,
#                                  ...){
#   .reqCols <- c("mz", "mobility", "sample", "into")
#   if (!all(.reqCols %in% colnames(peaks))){
#     stop("Required columns ",
#          paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
#                 collapse = ", "), " not found in 'peaks' parameter")
#   }
#
#   sample_name <- row.names(sampleGroups)
#   sampleGroups <- as.character(sampleGroups$class)
#   names(sampleGroups) <- sample_name
#   sampleGroupNames <- unique(sampleGroups)
#   sampleGroupTable <- table(sampleGroups)
#   nSampleGroups <- length(sampleGroupTable)
#
#   row.names(peaks) <- paste0('#', seq(nrow(peaks)))
#
#   unique_pg <- unique(peaks$peak_groups)
#   res <- lapply(unique_pg, function(pg){
#     # pg <- unique_pg[100]
#     idx <- which(peaks$peak_groups == pg)
#     temp_peaks <- peaks[idx, ]
#
#     tt <- table(sampleGroups[unique(temp_peaks[, "sample"])])
#     if (!any(tt / sampleGroupTable[names(tt)] >= minFraction &
#              tt >= minSamples)){
#       return(NULL)
#     }
#
#     gcount <- rep(0, length(sampleGroupNames))
#     names(gcount) <- sampleGroupNames
#     gcount[names(tt)] <- as.numeric(tt)
#
#
#     # res_mat <- rbind(res_mat,
#     #                  c(median(temp_peaks[, "mz"]),
#     #                    range(temp_peaks[, "mz"]),
#     #                    median(temp_peaks[, "mobility"]),
#     #                    range(temp_peaks[, "mobility"]),
#     #                    length(idx),
#     #                    gcount)
#     # )
#     res_ql <- as.data.frame(matrix(c(median(temp_peaks[, "mz"]),
#                             range(temp_peaks[, "mz"]),
#                             median(temp_peaks[, "mobility"]),
#                             range(temp_peaks[, "mobility"]),
#                             length(idx),
#                             gcount), nrow = 1))
#     colnames(res_ql) <- c("mz", "mzmin", "mzmax", "mobility", "mobilitymin", "mobilitymax",
#                           "npeaks", sampleGroupNames)
#     res_qt <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(sample_name)))
#     colnames(res_qt) <- sample_name
#
#     idx_ss <- match(temp_peaks$sample,
#                     colnames(res_qt))
#
#     res_qt[1, idx_ss] <- temp_peaks$into
#
#     res_this_pg <- cbind(res_ql, res_qt)
#     return(res_this_pg)
#   })
#
#   res_df <- do.call(rbind, res)
#   return(res_df)
# }





