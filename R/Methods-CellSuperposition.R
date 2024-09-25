#' @export
setMethod(
  "QueryTimsDataSegment",
  signature = c("TimsData", "QueryTimsDataSegmentParam", "DiscoverSCEventsParam"),
  function(object, param, sc_param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Querying tims data...")

    param_list <- c(as.list(param), as.list(sc_param))
    par_idx <- .gen_split_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'tims_data')
    names(files) <- object@files

    tims_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c("data_file" = data_file,
          param_list)
      })
      .parallel_parser(".query_sctims_data_segment", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files <- tims_data_files

    setwd(wd0)
    return(object)
  })


#' @export
setMethod(
  "TrimSCTimsData",
  signature = c("TimsData", "TrimSCTimsDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    # browser()
    message("Trimming single-cell metabolomics tims data...")

    param_list <- as.list(param)
    par_idx <- .gen_split_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'trim_tims_data')
    names(files) <- object@files

    trim_sctims_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        # browser()
        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[data_file])),
          param_list)
      })
      # browser()
      .parallel_parser(".trim_tims_data", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(trim_sctims_data_files) <- object@files
    object@tmp_data_files$trim_sctims_data_files <- trim_sctims_data_files

    setwd(wd0)
    return(object)
  })


#' @export
setMethod(
  "DiscoverSCEvents",
  signature = c("TimsData", "DiscoverSCEventsParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    # browser()
    message("Discovering single-cell events...")

    param_list <- as.list(param)
    par_idx <- .gen_split_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'single_cell_events')
    names(files) <- object@files

    single_cell_events <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        # browser()
        c(list("trim_sctims_data_files" = unname(object@tmp_data_files$trim_sctims_data_files[data_file])),
          param_list)
      })
      .parallel_parser(".discover_sc_events", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(single_cell_events) <- object@files
    object@tmp_data_files$single_cell_events <- single_cell_events

    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "UniteCellSuperpositionFrame",
  signature = c("TimsData", "UniteCellSuperpositionFrameParam"),
  function(object, param) {
    # object <- tims_data
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Uniting cell superposition frame...")

    files <- .tmp_files('cell_superposition_frame', object@experiment@tmp_dir, 'cell_superposition_frame')
    names(files) <- 'cell_superposition_frame'

    if (.check_skip(files[1], param@rerun, show_message = FALSE)) {
      cat('  Using previouse results:\n    ')
      cat(paste0(files, collapse = '    \n    '))
      object@tmp_data_files$cell_superposition_frame <- files
      setwd(wd0)
      return(object)
    }
    tmp <- .unite_cell_superposition_data(object@tmp_data_files$tims_data_files[1], param@cell_number, file.path(object@experiment@tmp_dir, 'tmp'), pool_size = param@pool_size)
    data_files <- tmp$tmp_data_files
    all_mobility <- tmp$all_mobility

    param_list <- as.list(param)
    par_idx <- .gen_split_indexes(length(data_files), object@experiment@BPPARAM$workers)
    files_split <- .tmp_files(paste0('res', basename(data_files)), object@experiment@tmp_dir, 'cell_superposition_frame')

    tmp_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(data_files[idxs], function(data_file) {
        c("data_file" = data_file,
          "res_define_at" = object@experiment@res_define_at,
          param_list)
      })
      .parallel_parser(".unite_cell_superposition_frame", arg_list, files_split[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))


    res_data <- unlist(lapply(tmp_files, readRDS), recursive = FALSE)
    res_frame <- list()
    res_frame$all_frames <- list()
    res_frame$all_frames[[1]] <- data.table::rbindlist(res_data[as.character(sort(as.numeric(names(res_data))))])
    names(res_frame$all_frames)[1] <- '1'
    res_frame$all_mobility <- all_mobility
    scan_range <- range(res_frame$all_frames$`1`$scan)

    res_frame$all_mobility <- all_mobility[as.numeric(names(res_frame$all_mobility)) >= scan_range[1]
                                           & as.numeric(names(res_frame$all_mobility)) <= scan_range[2]]
    saveRDS(res_frame, file = files)

    file.remove(data_files)
    file.remove(tmp_files)
    object@tmp_data_files$cell_superposition_frame <- files
    setwd(wd0)
    return(object)
  })


#' @export
setMethod(
  "SearchPeakTargets",
  signature = c("TimsData", "SearchPeakTargetParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    ############################################################################
    message("Searching peak targets in the cell superposition frame...")
    files <- .tmp_files('peak_targets', object@experiment@tmp_dir, 'peak_targets')
    names(files) <- 'peak_targets'

    arglist_pt <- c(list("cell_superposition_frame_file" = object@tmp_data_files$cell_superposition_frame,
                         "res_define_at" = object@experiment@res_define_at),
                    as.list(param))

    res_pt <- .analysis_parser('.search_peak_target', arglist_pt, files)

    object@tmp_data_files$peak_targets <- files

    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "DetectEIMPeaks",
  signature = c("TimsData", "DetectEIMPeaksParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    ############################################################################
    message("Detecting EIM peaks in the cell superposition frame...")
    # sub_dir <- file.path('group_peaks', 'density')
    files <- .tmp_files('eim_peaks', object@experiment@tmp_dir, 'eim_peaks')
    # files <- file.path(getwd(), object@experiment@tmp_dir, 'group_sc_peak', 'group_sc_peak_density')
    names(files) <- 'eim_peaks'
    arglist <- c(list("cell_superposition_frame_file" = object@tmp_data_files$cell_superposition_frame,
                      "res_pktg_file" = object@tmp_data_files$peak_targets),
                 as.list(param))

    res <- .analysis_parser('.detect_peak_cell_superposition',
                            arglist,
                            files, TRUE)

    object@tmp_data_files$eim_peaks <- files
    ############################################################################

    setwd(wd0)
    return(object)
  })




#' @export
setMethod(
  "IdentifySCPeaks",
  signature = c("TimsData", "IsotopeParam", "SCMatchParam"),
  function(object, isoparam, param) {
    # browser()
    wd0 <- getwd()
    # object <- tims_data
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, isoparam)
    object@.processHistory <- c(object@.processHistory, param)

    ############################################################################
    message("Annotating isotope")
    iso_file <- .tmp_files("isotope_files", object@experiment@tmp_dir, "isotope_files")
    names(iso_file) <- "isotope_files"

    info_col_isotope <- colnames(object@features)[!colnames(object@features) %in% object@tmp_data_files$detected_cell_number]
    sample_col_isotope <- object@tmp_data_files$detected_cell_number

    temp_feature <- object@features[order(object@features$int, decreasing = TRUE), info_col_isotope]

    iso_arglist <- c(list("temp_feature" = temp_feature[temp_feature$charge == 1, ],
                          "res_define_at" = object@experiment@res_define_at),
                     as.list(isoparam))

    res <- .analysis_parser('.annotate_isotope',
                            iso_arglist,
                            iso_file, TRUE)

    object@tmp_data_files$isotope_files <- iso_file


    ############################################################################
    message("Identifying single-cell peaks...")

    iso_res <- readRDS(object@tmp_data_files$isotope_files)
    id_feature <- object@features
    id_feature <- id_feature[id_feature$charge == 1, ]
    match_iso <- match(id_feature$name,
                       iso_res$name)
    id_feature$base_peak <- iso_res$base_peak[match_iso]
    id_feature$isotope <- iso_res$isotope[match_iso]

    feature_with_iso <- id_feature[, c(info_col_isotope, 'base_peak', 'isotope', sample_col_isotope)]

    write.csv(feature_with_iso, file.path(object@experiment@res_dir, "02_isotope_annotation_table.csv"), row.names = FALSE)
    # files <- .tmp_files("identify_peaks", object@experiment@tmp_dir, "identify_peaks")
    # names(files) <- "identify_peaks"

    pkg <- getPackageName()
    if (pkg == ".GlobalEnv") {
      pkg <- "MetCell"
    }
    library <- param@library
    # load the library
    if (is.null(library)) {
      lib_data <- readRDS(system.file(package = pkg, 'library', 'ms1info'))
    }else {
      lib_data <- read.csv(lib,
                           stringsAsFactors = FALSE)
    }

    # load the experimental data
    # exp_data <- object@features

    adduct_table <- read.csv(system.file('adducts',
                                         paste0("adducts_", object@experiment@ion_mode, ".csv"),
                                         package = pkg), stringsAsFactors = FALSE)

    if (!is.null(param@adduct_for_id)) {

      adm <- match(param@adduct_for_id,
                   adduct_table$name)
      adduct_table <- adduct_table[adm,]

    }


    toleranceCCS <- param@toleranceCCS[1]
    typeCCS <- param@typeCCS
    tolerancemz <- param@tolerancemz


    info_col <- colnames(id_feature)[!colnames(id_feature) %in% object@tmp_data_files$detected_cell_number]
    other_col <- object@tmp_data_files$detected_cell_number

    # info_col <- colnames(object@features)[which(!colnames(object@features) %in% row.names(object@sample_groups))]
    # other_col <- colnames(object@features)[which(colnames(object@features) %in% row.names(object@sample_groups))]

    match_res <- .match_metabolite_info(id_feature,
                                        lib_data,
                                        adduct_table,
                                        toleranceCCS,
                                        tolerancemz,
                                        typeCCS,
                                        object@experiment@res_define_at,
                                        info_col,
                                        other_col)

    # final_result <- cbind(object@features,
    #                       match_res)


    write.csv(match_res, file.path(object@experiment@res_dir, "03_metabolite_annotation_table.csv"), row.names = FALSE)
    setwd(wd0)
    return(object)
  })


#' #' @export
#' setMethod(
#'   "GenerateSCFeature",
#'   signature = c("TimsData", "GenerateSCFeatureParam"),
#'   function(object, param) {
#'     # object <- tims_data
#'     wd0 <- getwd()
#'     setwd(object@experiment@wd)
#'     object@.processHistory <- c(object@.processHistory, param)
#'
#'     message("Getting single-cell features...")
#'
#'     # sub_dir <- file.path('group_peaks', 'density')
#'     files <- .tmp_files('single_cell_feature', object@experiment@tmp_dir, 'single_cell_feature')
#'     # files <- file.path(getwd(), object@experiment@tmp_dir, 'group_sc_peak', 'group_sc_peak_density')
#'     names(files) <- 'single_cell_feature'
#'     # browser()
#'     old_sample_groups <- object@sample_groups
#'     temp_sg <- object@peaks[, c('sample', 'raw_file')]
#'     temp_sg <- temp_sg[!duplicated(temp_sg$sample), ]
#'     idx <- match(temp_sg$raw_file, row.names(old_sample_groups))
#'     sample_groups <- data.frame(class = old_sample_groups$class[idx])
#'     row.names(sample_groups) <- as.character(temp_sg$sample)
#'     sampleGroups <- sample_groups
#'     arglist <- c(list("peaks" = object@peaks,
#'                       "sampleGroups" = sampleGroups),
#'                  as.list(param))
#'
#'     res <- .analysis_parser('.generate_sc_feature',
#'                             arglist,
#'                             files, TRUE)
#'
#'     object@tmp_data_files$single_cell_feature <- files
#'
#'     final_df <- readRDS(object@tmp_data_files$single_cell_feature)
#'
#'     object@features <- final_df
#'     # object@peak_groups <- res
#'
#'     write.csv(object@features, file.path(object@experiment@res_dir, "features.csv"), row.names = FALSE)
#'     setwd(wd0)
#'     return(object)
#'   })


