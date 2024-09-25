#' @export
setMethod(
  "ExtractSCIMData",
  "signature" = c("TimsData", "ExtractSCIMDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Extracting single-cell IM data...")
    # browser()
    param_list <- as.list(param)
    # par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    
    # files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'single_cell_EIM')
    names(files) <- object@files
    
    # sc_im_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
    #   idxs <- seq(dr[1], dr[2])
    #   cat('processing files: ', dr[1], '-', dr[2], '\n')
    #   # browser()
    #   arg_list <- lapply(object@files[idxs], function(data_file) {
    #     c(list("trim_sctims_data_files" = unname(object@tmp_data_files$trim_sctims_data_files[data_file]),
    #            "single_cell_events" = unname(object@tmp_data_files$single_cell_events[data_file]),
    #            "res_define_at" = object@experiment@res_define_at),
    #       param_list)
    #   })
    #   # browser()
    #   # .parallel_parser(".extract_single_cell_im_data", arg_list, files[idxs], object@experiment@BPPARAM)
    #   .parallel_parser(".extract_single_cell_im_data_local_max", arg_list, files[idxs], object@experiment@BPPARAM)
    #   # .analysis_parser(".extract_single_cell_im_data_local_max", arg_list[idxs], files[idxs])
    # }, simplify = FALSE))
    
    #### use .analysis_parser to replace .parallel_parser ####
    sc_im_data_files <- do.call(c, lapply(object@files, function(file) {
      # idxs <- seq(dr[1], dr[2])
      cat('processing files: ', file, '\n')
      # browser()
      # arg_list <- lapply(object@files[idxs], function(data_file) {
      #   c(list("trim_sctims_data_files" = unname(object@tmp_data_files$trim_sctims_data_files[data_file]),
      #          "single_cell_events" = unname(object@tmp_data_files$single_cell_events[data_file]),
      #          "res_define_at" = object@experiment@res_define_at),
      #     param_list)
      # })
      
      arg_list <- c(list("trim_sctims_data_files" = unname(object@tmp_data_files$trim_sctims_data_files[file]),
                         "single_cell_events" = unname(object@tmp_data_files$single_cell_events[file]),
                         "res_define_at" = object@experiment@res_define_at),
                    param_list)
      
      .analysis_parser(".extract_single_cell_im_data_local_max", arg_list, files[file])
    }))
    
    names(sc_im_data_files) <- object@files
    object@tmp_data_files$sc_im_data_files <- sc_im_data_files
    
    object@peaks <- do.call(rbind, lapply(sc_im_data_files, function(peak_file) {
      peaks <- readRDS(peak_file)
      rownames(peaks) <- NULL
      peaks
    }))
    rownames(object@peaks) <- .gen_indexes(object@peaks)
    
    # object <- tims_data
    # object@peaks$raw_file <- sapply(object@peaks$sample, function(sss){
    #   paste0(strsplit(sss, '_')[[1]][1:2], collapse = '_')
    # })
    # browser()
    old_sample_groups <- object@sample_groups
    temp_sg <- object@peaks[, c('sample', 'raw_file')]
    temp_sg <- temp_sg[!duplicated(temp_sg$sample), ]
    idx <- match(temp_sg$raw_file, row.names(old_sample_groups))
    sample_groups <- data.frame(class = old_sample_groups$class[idx])
    row.names(sample_groups) <- as.character(temp_sg$sample)
    object@sample_groups <- sample_groups
    
    
    
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "GroupSCPeaks",
  signature = c("TimsData", "GroupSCDensityParam"),
  function(object, param) {
    # object <- tims_data
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Grouping SC peaks using density method...")
    
    # sub_dir <- file.path('group_peaks', 'density')
    files <- .tmp_files('group_data', object@experiment@tmp_dir, 'group_peaks')
    # files <- file.path(getwd(), object@experiment@tmp_dir, 'group_sc_peak', 'group_sc_peak_density')
    names(files) <- 'group_sc_peak'
    arglist <- c(list("peaks" = object@peaks,
                      "sampleGroups" = object@sample_groups,
                      "plot_dir" = dirname(files)),
                 as.list(param))
    
    res <- .analysis_parser('.group_sc_peaks_density',
                            arglist,
                            files, TRUE)
    
    # generate the feature table with quantification information with res above and peaks in the objects.
    # empty_df <- data.frame(matrix(data = NA, ncol = nrow(object@sample_groups), nrow = nrow(res)))
    all_sample <- row.names(object@sample_groups)
    # browser()
    ff_feature <- lapply(seq(nrow(res)), function(i){
      # browser()
      # i <- 3
      temp_idx <- res$peakidx[[i]]
      temp_peaks <- object@peaks[row.names(object@peaks) %in% temp_idx,]
      # sample_idx <- match(all_sample, temp_peaks$sample)
      
      sample_idx <- sapply(all_sample, function(nms){
        # browser()
        this_feature_this_sample <- which(temp_peaks$sample == nms)
        if(length(this_feature_this_sample) > 1){
          min_mobility_idx <- which.min(abs(res$mobilitymed[i] - temp_peaks$mobility[this_feature_this_sample]))
          return(this_feature_this_sample[min_mobility_idx])
        }else if(length(this_feature_this_sample) == 0){
          return(NA)
        }else if(length(this_feature_this_sample) == 1){
          return(this_feature_this_sample)
        }
      }, simplify = TRUE)
      
      this_feature_df <- data.frame(matrix(data = temp_peaks$into[sample_idx], ncol = length(sample_idx)))
      colnames(this_feature_df) <- all_sample
      return(this_feature_df)
      # return(sample_idx)
    })
    
    ff_feature <- do.call(rbind, ff_feature)
    
    # col_features <- colnames(res$features)
    # pk_nms <- rownames(object@peaks)
    # res$peak_groups <- merged_indexes[rownames(res$features)[rownames(res$features) %in% names(merged_indexes)]]
    # res$features <- as.data.frame(t(sapply(res$peak_groups, function(idx) {
    #   peaks <- object@peaks[fastmatch::fmatch(idx, pk_nms), , drop = FALSE]
    #   c(median(peaks[, "mz"]),
    #     range(peaks[, "mz"]),
    #     median(peaks[, "mobility"]),
    #     range(peaks[, "mobility"]),
    #     median(peaks[, "rt_align"]),
    #     range(peaks[, "rt_align"]),
    #     nrow(peaks),
    #     colMedians(peaks[, c("ccs", "target_intensity", "height", "height_fit", "area"), drop = FALSE]))
    # })), stringsAsFactors = FALSE)
    # colnames(res$features) <- col_features
    
    # object@.processHistory
    res_temp <- res[, 1:ncol(res) - 1]
    res_temp$ccs <- .mobility2ccs(res_temp$mobilitymed, res_temp$mzmed)
    final_df <- cbind(res_temp, ff_feature)
    final_df <- cbind(data.frame(name = .get_peak_name2(final_df[, c("mzmed", "ccs")])), final_df)
    object@tmp_data_files$align_file <- files
    object@features <- final_df
    object@peak_groups <- res
    
    write.csv(object@features, file.path(object@experiment@res_dir, "features.csv"), row.names = FALSE)
    
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "FillSCPeaks",
  signature = c("TimsData", "FillSCPeakParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    # object <- tims_data
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Filling single-cell peaks...")
    files <- .tmp_files(object@files, object@experiment@tmp_dir, "fill_sc_peaks")
    names(files) <- object@files
    
    param_list <- as.list(param)
    # rev_models <- readRDS(object@tmp_data_files$align_file)$rev_models
    # rt_correct_files <- object@tmp_data_files$correct_rt_files
    smp_names <- rownames(object@sample_groups)
    raw_data_file <- sapply(smp_names, function(smp){
      strsplit(smp, '@')[[1]][1]
    })
    
    
    data_index <- data.frame(smp_names = smp_names, 
                             raw_data_file = raw_data_file)
    
    
    data_index$is_filled_sample <- sapply(smp_names, function(i){
      any(is.na(object@features[, i]))
    })
    
    
    
    par_idx <- .gen_split_indexes(length(object@files), object@experiment@BPPARAM$workers)
    fill_scpeak_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        # browser()
        temp_raw_data <- basename(data_file)
        match_sample <- data_index$smp_names[which(data_index$raw_data_file == temp_raw_data)]
        
        # browser()
        info <- lapply(match_sample, function(sp){
          is_fill <- is.na(object@features[, sp])
          res <- object@features[is_fill, c('mzmed', 'mobilitymed', 'name'), drop = FALSE]
          colnames(res) <- c('mz', 'mobility', 'name')
          res$target_frame <- strsplit(sp, '@')[[1]][2]
          res$sample <- sp
          return(res)
        })
        
        names(info) <- match_sample
        
        # is_fill <- is.na(object@features[, smp_names[idx]])
        # info <- object@features[is_fill, , drop = FALSE]
        
        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[data_file]),
               "info" = info,
               "res_define_at" = object@experiment@res_define_at, 
               'data_file' = unname(files[data_file])),
          param_list)
        # browser()
      })
      # browser()
      .parallel_parser(".fill_single_cell_peaks", arg_list, files[idxs],
                       object@experiment@BPPARAM, save_in_analysis = TRUE)
    }, simplify = FALSE))
    names(fill_scpeak_files) <- object@files
    object@tmp_data_files$fill_scpeak_files <- fill_scpeak_files
    
    filled_peaks <- do.call(rbind, lapply(fill_scpeak_files, readRDS))
    row.names(filled_peaks) <- NULL
    
    object@filled_peaks <- filled_peaks
    
    ### file the peak area to the feature table 
    for (smp in unique(filled_peaks$sample)) {
      # smp <- 'metabolomics_3.d@3824'
      temp_filled <- filled_peaks[filled_peaks$sample == smp, ]
      
      idx <- match(temp_filled$name, object@features$name)
      object@features[idx, smp] <- temp_filled$extracted_intensity
      # is_fill <- is.na(object@features[, smp_names[idx]])
      # object@features[is_fill, smp_names[idx]] <- readRDS(fill_peak_files[idx])
    }
    write.csv(object@features, file.path(object@experiment@res_dir, "features_filled.csv"), row.names = FALSE)
    setwd(wd0)
    return(object)
  })
