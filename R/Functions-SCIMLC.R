################################################################################
# trim_sctims_data_files = 'E:/04_single_cell/00_package/08_pressure_cell/results/tmp/trim_tims_data/sc06_4e4_450psi.d'
# # single_cell_events = 'E:/04_single_cell/00_package/08_pressure_cell/results/tmp/single_cell_events/sc06_4e4_450psi.d'
# trim_sctims_data_files = 'E:/04_single_cell/00_package/05_sc_local_max_long/results/tmp/trim_tims_data/metabolomics_3.d'
# single_cell_events = 'E:/04_single_cell/00_package/05_sc_local_max_long/results/tmp/single_cell_events/metabolomics_3.d'
# # trim_sctims_data_files = 'E:/04_single_cell/00_package/10_pressure/results/tmp/trim_tims_data/03_sc.d'
# # single_cell_events = 'E:/04_single_cell/00_package/10_pressure/results/tmp/single_cell_events/03_sc.d'
# # trim_sctims_data_files = 'E:/04_single_cell/00_package/11_pressure_2/results/tmp/trim_tims_data/02_sc.d'
# # single_cell_events = 'E:/04_single_cell/00_package/11_pressure_2/results/tmp/single_cell_events/02_sc.d'
# # tims_data <- readRDS('E:/04_single_cell/00_package/11_pressure_2/tims_data')
# # bpparam <- tims_data@experiment@BPPARAM
# mz_tol = 20
# res_define_at = 200
# mobility_tol = 0.05 # the range to search the roi around the selected ion
# min_intensity = 0 # absolute intensity threshold
# min_points = 10 # minimal point to generate a ROI
# n_skip = 1 # the maximum point to skip when finding a roi
# smooth_window = 10
# peak_span_eim = 11
# snthreshold = 3
# skip_invalid_peaks = TRUE
# rerun = FALSE

.extract_single_cell_im_data_local_max <- function(
  # about loading the data
  trim_sctims_data_files,
  single_cell_events,
  
  # about find roi
  mz_tol, # ppm 
  res_define_at, 
  mobility_tol = 0.05, # the range to search the roi around the selected ion  
  min_intensity = 0, # absolute intensity threshold
  min_points = 10, # minimal point to generate a ROI
  n_skip = 1, # the maximum point to skip when finding a roi
  
  # about peak detection with local maximum
  smooth_window = 10, 
  peak_span_eim = 11, 
  snthreshold = 3, 
  skip_invalid_peaks = TRUE,
  ...
) {
  # browser()
  query_data <- readRDS(trim_sctims_data_files)
  sc_event <- readRDS(single_cell_events)
  sc_index <- as.character(sc_event$frame_index[which(sc_event$sc)])
  # sc_index <- sc_index[1:10]
  # names(query_data$all_frames)
  query_data$all_frames <- query_data$all_frames[which(names(query_data$all_frames) %in% sc_index)]
  all_mobility <- query_data$all_mobility
  
  # all_st <- Sys.time()
  cl <- parallel::makeCluster(10L)
  res <- parallel::parLapply(cl = cl, sc_index, function(sci, 
                                                         query_data, 
                                                         mz_tol, 
                                                         min_points, 
                                                         min_intensity, 
                                                         n_skip, 
                                                         res_define_at, 
                                                         mobility_tol, 
                                                         .ppm2dalton, 
                                                         .get_eims3, 
                                                         .find_roi){
  # res <- parallel::parLapply(cl = cl, sc_index, function(sci, query_data){
  # res <- pbapply::pblapply(sc_index, function(sci){
  # res <- lapply(sc_index, function(sci){
    # cat(sci, '\t')
    # sci <- '2826'
    temp_frame <- query_data$all_frames[sci]
    temp_frame[[1]]$num_ion <- paste0('#', seq(nrow(temp_frame[[1]])))
    temp_frame[[1]] <- temp_frame[[1]][order(temp_frame[[1]]$intensity, decreasing = TRUE), ]
    
    
    if (mz_tol > 1) {
      temp_frame[[1]]$mz_tol <- sapply(temp_frame[[1]]$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
    } else {
      temp_frame[[1]]$mz_tol <- mz_tol
    }
    
    temp_frame[[1]]$mz_up <- temp_frame[[1]]$mz + temp_frame[[1]]$mz_tol
    temp_frame[[1]]$mz_low <- temp_frame[[1]]$mz - temp_frame[[1]]$mz_tol
    
    
    # st_time <- Sys.time()
    first_rm <- sapply(seq(nrow(temp_frame[[1]])), function(i){
      # cat(i, '\t')
      # length(which(temp_frame[[1]]$mz_up >= temp_frame[[1]]$mz[i] &
      #                temp_frame[[1]]$mz_low <= temp_frame[[1]]$mz[i] &
      #                abs(temp_frame[[1]]$scan - temp_frame[[1]]$scan[i] < 50))) >= min_points/2
      length(which(temp_frame[[1]]$mz_up >= temp_frame[[1]]$mz[i] &
                     temp_frame[[1]]$mz_low <= temp_frame[[1]]$mz[i])) >= min_points/2
      },
      simplify = TRUE)

    temp_frame[[1]] <- temp_frame[[1]][which(first_rm), ]
    # temp_frame[[1]] <- temp_frame[[1]][which(temp_frame[[1]]$intensity > 16), ]
    roi_list <- list()
    
    while (nrow(temp_frame[[1]]) > 1) {
      # browser()
      # cat(nrow(temp_frame[[1]]), '\t')
      this_ion <- temp_frame[[1]][1, ]
      # if(length(which(temp_frame[[1]]$mz_up >= this_ion$mz &
      #                 temp_frame[[1]]$mz_low <= this_ion$mz)) < min_points/2){
      #   temp_frame[[1]] <- temp_frame[[1]][-1, ]
      #   next
      # }
      # 
      temp_eim <- .get_eims3(temp_frame, query_data$all_mobility,
                             mz = this_ion$mz,
                             mz_tol = this_ion$mz_tol,
                             mobility = query_data$all_mobility[as.character(this_ion$scan)],
                             mobility_tol = mobility_tol)
      # if(length(which(!is.na(temp_eim$mz))) < min_points/2){
      #     temp_frame[[1]] <- temp_frame[[1]][-1, ]
      #     next
      # }
      
      temp_roi <- .find_roi(temp_eim$intensity, min_intensity = min_intensity,
                            min_points = min_points,
                            n_skip = n_skip,
                            ref = match(this_ion$num_ion, temp_eim$num_ion))
      if(!is.null(temp_roi)){
        for(t in 1:nrow(temp_roi)){
          temp_eim$roi <- FALSE
          temp_eim$roi[temp_roi[t, 1]: temp_roi[t, 2]] <- TRUE
          roi_list[[length(roi_list) + 1]] <- temp_eim
          rm_idx <- temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion[which(!is.na(temp_eim[temp_roi[t, 1]: temp_roi[t, 2], ]$num_ion))]
          idx <- match(rm_idx, temp_frame[[1]]$num_ion)
          temp_frame[[1]] <- temp_frame[[1]][-idx, ]
        }
      }else{
        temp_frame[[1]] <- temp_frame[[1]][-1, ]
      }
    }
    # ed_time <- Sys.time()
    # cat(ed_time - st_time, '\n')
    return(roi_list)
    gc()
  # })
  }, 
  query_data = query_data, 
  mz_tol = mz_tol, 
  min_points = min_points, 
  min_intensity = min_intensity, 
  n_skip = n_skip, 
  res_define_at = res_define_at, 
  mobility_tol = mobility_tol, 
  .ppm2dalton = .ppm2dalton, 
  .get_eims3 = .get_eims3, 
  .find_roi = .find_roi)
  
  parallel::stopCluster(cl)
  # all_ed <- Sys.time()
  # all_ed - all_st
  names(res) <- sc_index
  
  res_peaks_in_all_cells <- lapply(seq_along(res), function(cell_index){
    # cell_index <- 1
    one_cell <- res[[cell_index]]
    res_in_one_cell <- lapply(seq_along(one_cell), function(roi_index){
      # roi_index <- 1
      temp <- one_cell[[roi_index]]
      temp$intensity_smooth <- .smooth_loess(data = temp$intensity, 
                                             degree = 1, 
                                             window = smooth_window)
      ### ref_index should be the maximun in the roi ###
      apex_idx <- .find_apex(temp$intensity, 
                             temp$intensity_smooth,
                             span = peak_span_eim,
                             min_points = min_points,
                             min_intensity = min_intensity,
                             n_skip = n_skip,
                             ref_index = which.max(temp$intensity_smooth),
                             find_roi = FALSE)
      apex_eim <- temp$k0[apex_idx]
      
      if (length(apex_eim) == 0) {
        return(NULL)
      }
      
      peak_bd <- temp$k0[.findLocalMin(temp$intensity_smooth)]
      
      peak_bd_1 <- peak_bd[which(apex_eim > peak_bd)]
      if(length(peak_bd_1) == 1){
        peak_bd_1 <- peak_bd_1
      } else if(length(peak_bd_1) == 0){
        peak_bd_1 <- apex_eim - 0.015
      } else if(length(peak_bd_1) > 1){
        peak_bd_1 <- peak_bd_1[which.min(abs(apex_eim - peak_bd_1 - 0.015))]
      }
      
      peak_bd_2 <- peak_bd[which(apex_eim < peak_bd)]
      if(length(peak_bd_2) == 1){
        peak_bd_2 <- peak_bd_2
      } else if(length(peak_bd_2) == 0){
        peak_bd_2 <- apex_eim + 0.015
      } else if(length(peak_bd_2) > 1){
        peak_bd_2 <- peak_bd_2[which.min(abs(peak_bd_2 - apex_eim - 0.015))]
      }
      
      peak_range <- temp[temp$k0 >= peak_bd_1 & temp$k0 <= peak_bd_2, ]
      peak_range <- peak_range[!is.na(peak_range$mz),]
      new_mz <- .get_weighted_mz(peak_range$mz, peak_range$intensity)
      peak_intensity <- sum(peak_range$intensity)
      
      
      peak_sd <- .get_smooth_sd(temp$intensity, temp$intensity_smooth)
      
      peak_quanlity <-.get_peak_quality(temp$intensity,
                                        apex_idx,
                                        min_intensity = min_intensity,
                                        min_points = min_points,
                                        n_skip = n_skip,
                                        snthreshold = 3,
                                        skip_invalid_peaks = skip_invalid_peaks)
      
      
      if(peak_sd > 0.35 && skip_invalid_peaks){
        return(NULL)
      }
      
      if(is.null(peak_quanlity) && skip_invalid_peaks){
        return(NULL)
      }
      
      final_df <- data.frame(mobility = apex_eim, 
                             mobility_min = peak_bd_1, 
                             mobility_max = peak_bd_2, 
                             mz = new_mz, 
                             into = peak_intensity, 
                             sn = peak_quanlity$snr, 
                             baseline = peak_quanlity$baseline)
      return(final_df)
      
      })
    return(res_in_one_cell)
  })
  
  names(res_peaks_in_all_cells) <- sc_index
  
  res_all_cell_peaks_df <- lapply(seq_along(res_peaks_in_all_cells), function(ccc){
    # browser()
    one_cell <- res_peaks_in_all_cells[[ccc]]
    file_name <- basename(trim_sctims_data_files)
    cell_name <- paste0(file_name, '@', names(res_peaks_in_all_cells)[ccc])
    res_for_peaks <- lapply(one_cell, function(i){
      if(!is.null(i)){
        final_df <- i
        final_df$sample <- cell_name
        final_df$raw_file <- file_name
        return(final_df)
      }else{
        return(NULL)
      }
    })
    res_for_peaks <- do.call(rbind, res_for_peaks)
    return(res_for_peaks)
  })
  res_all_cell_peaks_df <- do.call(rbind, res_all_cell_peaks_df)

}



.findLocalMin <- function(x, m = 10){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape > 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] >= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  unname(pks)
}
# res <- readRDS('./20230807_intermedian/roi_detected_res.rds')
# pdf('test_local_minimum_10pts_10min.pdf')
# 
# for(i in seq_along(res[["3824"]])){
#   temp <- res[["3824"]][[i]]
#   temp$intensity_smooth <- .smooth_loess(data = temp$intensity, 
#                                          degree = 1, 
#                                          window = 10)
#   plot(temp$k0, temp$intensity, type = 'o')
#   lines(temp$k0, temp$intensity_smooth, type = 'l', col = 'red')
#   
#   .get_smooth_sd(temp$intensity, temp$intensity_smooth)
#   
#   .get_peak_quality(temp$intensity,
#                     which.max(temp$intensity_smooth),
#                     min_intensity = 0,
#                     min_points = 7,
#                     n_skip = 1,
#                     snthreshold = 3,
#                     skip_invalid_peaks = TRUE)
#   abline(v = temp$k0[which.max(temp$intensity_smooth)], col = 'red')
#   abline(v = temp$k0[findLocalMin(temp$intensity_smooth)], col = 'blue')
# }
# 
# dev.off()


################################################################################

.group_sc_peaks_density <- function(
    peaks,
    sampleGroups, 
    bw = 0.02,
    mz_bin_size = 0.01,
    max_features = 50,
    res_define_at = 200,
    minFraction = 0.2, 
    minSamples = 10, 
    plot_density = FALSE,
    plot_dir = '.',
    ...
) {
  # modify from xcms
  # sampleGroups <- peaks$cell_index
  .reqCols <- c("mz", "mobility", "sample", "into")
  if (!all(.reqCols %in% colnames(peaks))){
    stop("Required columns ",
         paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
                collapse = ", "), " not found in 'peaks' parameter")
  }
  
  sample_name <- row.names(sampleGroups)
  sampleGroups <- as.character(sampleGroups$class)
  names(sampleGroups) <- sample_name
  sampleGroupNames <- unique(sampleGroups)
  sampleGroupTable <- table(sampleGroups)
  nSampleGroups <- length(sampleGroupTable)
  
  # if (max(peaks[, "sample"]) > length(sampleGroups)){
  #   stop("Sample indices in 'peaks' are larger than there are sample",
  #        " groups specified with 'sampleGroups'!")
  # }
  
  
  peaks <- cbind(peaks[, .reqCols, drop = FALSE],
                 index = row.names(peaks))
  
  ## Order peaks matrix by mz
  peaks <- peaks[order(peaks[, "mz"]), , drop = FALSE]
  rownames(peaks) <- NULL
  MolRange <- range(peaks[, "mobility"])
  
  
  ## Define the mass slices and the index in the peaks matrix with an mz
  ## value >= mass[i].
  mass <- seq(peaks[1, "mz"], peaks[nrow(peaks), "mz"] + mz_bin_size,
              by = mz_bin_size / 2)
  masspos <- find_greater_equal_than(peaks[, "mz"], mass)
  
  densFrom <- MolRange[1] - 3 * bw
  densTo <- MolRange[2] + 3 * bw
  ## Increase the number of sampling points for the density distribution.
  densN <- max(512, 2 * 2^(ceiling(log2(diff(MolRange) / (bw / 2)))))
  endIdx <- 0
  message("Processing ", length(mass) - 1, " mz slices ... ",
          appendLF = FALSE)
  resL <- vector("list", (length(mass) - 2))
  for (i in seq_len(length(mass)-2)) {
    # cat(i, '\t')
    ## That's identifying overlapping mz slices.
    startIdx <- masspos[i]
    endIdx <- masspos[i + 2] - 1
    if (endIdx - startIdx < 0)
      next
    resL[[i]] <- .group_sc_density(peaks[startIdx:endIdx, , drop = FALSE],
                                   bw = bw, densFrom = densFrom,
                                   densTo = densTo, densN = densN,
                                   sampleGroups = sampleGroups,
                                   sampleGroupTable = sampleGroupTable,
                                   minFraction = minFraction,
                                   minSamples = minSamples,
                                   maxFeatures = max_features)
  }
  message("OK")
  res <- do.call(rbind, resL)
  
  # if (nrow(res)) {
  #   ## Remove groups that overlap with more "well-behaved" groups
  #   numsamp <- rowSums(
  #     as.matrix(res[, (match("npeaks", colnames(res)) +1):(ncol(res) -1),
  #                   drop = FALSE]))
  #   uorder <- order(-numsamp, res[, "npeaks"])
  #   
  #   uindex <- rect_unique(
  #     as.matrix(res[, c("mzmin", "mzmax", "mobilitymin", "mobilitymin"),
  #                   drop = FALSE]), uorder-1, c(0,0))
  #   res <- res[uindex, , drop = FALSE]
  #   rownames(res) <- NULL
  # }
  
  if (nrow(res)) {
    ## Remove groups that overlap with more "well-behaved" groups
    uorder <- order(res[, "npeaks"])
    
    uindex <- rect_unique(
      as.matrix(res[, c("mzmin", "mzmax", "mobilitymin", "mobilitymax"),
                    drop = FALSE]), uorder-1, c(0, 0))
    res <- res[which(uindex == 1), , drop = FALSE]
    # peaks2 <- peaks[!peaks$peak_idx %in% do.call(c, res_u[, 'peakidx']), , drop = FALSE]
  }
  
  res
}




.group_sc_density <- function(x, bw, densFrom, densTo, densN, sampleGroups,
                              sampleGroupTable, minFraction,
                              minSamples, maxFeatures) {
  # modify from xcms
  # browser()
  den <- density(x[, "mobility"], bw = bw, from = densFrom, to = densTo,
                 n = densN)
  maxden <- max(den$y)
  deny <- den$y
  sampleGroupNames <- names(sampleGroupTable)
  nSampleGroups <- length(sampleGroupNames)
  pk_index <- x$index
  col_nms <- c("mzmed", "mzmin", "mzmax", "mobilitymed", "mobilitymin", "mobilitymax",
               "npeaks", sampleGroupNames)
  res_mat <- matrix(nrow = 0, ncol = length(col_nms),
                    dimnames = list(character(), col_nms))
  res_idx <- list()
  while (deny[maxy <- which.max(deny)] > maxden / 20 && nrow(res_mat) <
         maxFeatures) {
    grange <- find_dense_min(deny, maxy - 1)
    deny[grange[1]:grange[2]] <- 0
    gidx <- which(x[,"mobility"] >= den$x[grange[1]] &
                    x[,"mobility"] <= den$x[grange[2]] &
                    !is.na(pk_index))
    ## Determine the sample group of the samples in which the peaks
    ## were detected and check if they correspond to the required limits.
    tt <- table(sampleGroups[unique(x[gidx, "sample"])])
    if (!any(tt / sampleGroupTable[names(tt)] >= minFraction &
             tt >= minSamples)){
      next
    }
    gcount <- rep(0, length(sampleGroupNames))
    names(gcount) <- sampleGroupNames
    gcount[names(tt)] <- as.numeric(tt)
    res_mat <- rbind(res_mat,
                     c(median(x[gidx, "mz"]),
                       range(x[gidx, "mz"]),
                       median(x[gidx, "mobility"]),
                       range(x[gidx, "mobility"]),
                       length(gidx),
                       gcount)
    )
    res_idx <- c(res_idx, list(unname(sort(x[gidx, "index"]))))
    pk_index[gidx] <- NA
  }
  # if (sleep > 0) {
  #   ## Plot the density
  #   plot(den, main = paste(round(min(x[,"mz"]), 2), "-",
  #                          round(max(x[,"mz"]), 2)))
  #   ## Highlight peaks per sample group.
  #   for (j in seq_len(nSampleGroups)) {
  #     ## Which peaks belong to this sample group.
  #     cur_group_samples <- which(sampleGroups == sampleGroupNames[j])
  #     idx <- x[, "sample"] %in% cur_group_samples
  #     points(x[idx, "mobility"], x[idx, "into"] /
  #              max(x[, "into"]) * maxden,
  #            col = j, pch=20)
  #   }
  #   for (j in seq_len(nrow(res_mat)))
  #     abline(v = res_mat[j, 5:6], lty = "dashed", col = j)
  #   Sys.sleep(sleep)
  # }
  res <- as.data.frame(res_mat)
  res$peakidx <- res_idx
  res
}





.fill_single_cell_peaks <- function(
    info, # import as the feature table 
    tims_data_file,
    mz_tol = 20,
    frame_integration_range = 1,
    mobility_intgration_range = 0.02,
    res_define_at = 200,
    data_file = NULL,
    ...
) {
  # browser()
  query_data <- readRDS(tims_data_file)
  num_frames <- nrow(query_data$ms1_frame_info)
  
  info <- do.call(rbind, info)
  
  if (mz_tol > 1) {
    info$mz_tol <- sapply(info$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    info$mz_tol <- mz_tol
  }
  
  mobility_col <- "mobility"
  
  
  info$extracted_intensity <- apply(info[, c("mz", mobility_col, "mz_tol", "target_frame")], 1, function(dr) {
    # browser()
    mz <- as.numeric(dr["mz"])
    mobility <- as.numeric(dr[mobility_col])
    mz_tol <- as.numeric(dr["mz_tol"])
    target_frame <- as.character(dr["target_frame"])
    
    # target_idx <- which(query_data$ms1_frame_info$Id == target_frame)
    # frame_start <- max(0, target_idx - frame_integration_range)
    # extract_frames <- query_data$ms1_frame_info$Id[frame_start:min(num_frames, target_idx + frame_integration_range)]
    # num_extract_frames <- length(extract_frames)
    
    # get eim from the target frame
    eim <- .get_eims2(query_data$all_frames[as.character(target_frame)], query_data$all_mobility,
                      mz, mz_tol, mobility, mobility_intgration_range)
    sum(eim[, 1])
  })
  # browser()
  if (!is.null(data_file)) {
    saveRDS(info, file = data_file, version = 2)
  } else {
    return(info)
  }
}


