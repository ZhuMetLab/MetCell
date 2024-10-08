# project: MetCell
# File name: AllClasses.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:50
# Copyright (c) 2022- ZhuMSLab ALL right reserved

# setClass("Met4DXParam", contains = "VIRTUAL")
setClass("MetCellParam", contains = "VIRTUAL")

setClassUnion("nullOrCharacter", c("NULL", "character"))
setClassUnion("nullOrNumeric", c("NULL", "numeric"))
# setClassUnion("data.frameOrCharacter", c("data.frame", "character"))

setClass("SmoothParam", contains = "MetCellParam")

setClass("GaussianSmoothParam",
         slots = c(
           window = "numeric",
           alpha = "numeric"
         ),
         contains = "SmoothParam")

setClass("LOESSSmoothParam",
         slots = c(
           span = "numeric",
           degree = "numeric",
           window = "nullOrNumeric"
         ),
         contains = "SmoothParam")

setClass("FillPeakParam",
         slots = c(
           mz_tol = "numeric",
           frame_integration_range = "numeric",
           mobility_intgration_range = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("FinalizeFeatureParam",
         slots = c(
           min_fraction = "numeric",
           msms_assign_method = "character",
           min_num_samples = "numeric",
           valid_eic_peak = "logical",
           valid_eim_peak = "logical",
           snthreshold_eic = "nullOrNumeric",
           snthreshold_eim = "nullOrNumeric",
           quant_method = "character",
           col_max = "character",
           col_quant = "character",
           output_MS2 = 'logical',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("MatchBetweenRunParam",
         slots = c(
           mz_tol = "numeric",
           frame_range = "numeric",
           frame_integration_range = "numeric",
           mobility_range = "numeric",
           mobility_intgration_range = "numeric",
           min_points = "numeric",
           min_intensity = "numeric",
           n_skip = "numeric",
           interpolate_method = "nullOrCharacter",
           smooth_method = "character",
           peak_span_eim = "numeric",
           peak_span_eic = "numeric",
           smooth_window_eim = "numeric",
           smooth_window_eic = "numeric",
           keep_profile = "logical",
           filter_outlier_peaks = "logical",
           allowed_mobility_shift = "numeric",
           allowed_rt_shift = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("GroupParam", contains = "MetCellParam")

setClass("GroupDensityParam",
         slots = c(
           bw = "numeric",
           mz_bin_size = "numeric",
           mobility_bin_size = "numeric",
           max_features = "numeric",
           plot_density = "logical",
           rerun = "logical"
         ),
         contains = "GroupParam"
)

setClass("GroupLandmarksParam",
         slots = c(
           rerun = "logical"
         ),
         contains = "GroupParam"
)

setClass("CorrectRTParam", contains = "MetCellParam")

setClass("CorrectRTLandmarksParam",
         slots = c(
           mz_tol = "numeric",
           strict_rt_constrains = "logical",
           rt_tol_landmark = "numeric",
           rt_tol_match = "numeric",
           ccs_tol = "numeric",
           cutoff = "numeric",
           min_num_fragments = "numeric",
           method_match = "character",
           method_score = "character",
           weight_intensity = "numeric",
           weight_mz = "numeric",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           include_precursor = "logical",
           mz_tol_ms1 = "numeric",
           mz_tol_ms2 = "numeric",
           rerun = "logical"
         ),
         contains = "CorrectRTParam"
)

setClass("AlignPeakParam",
         slots = c(
           mz_tol = "numeric",
           strict_rt_constrains = "logical",
           rt_tol_landmark = "numeric",
           rt_tol_match = "numeric",
           ccs_tol = "numeric",
           cutoff = "numeric",
           min_num_fragments = "numeric",
           method_match = "character",
           method_score = "character",
           weight_intensity = "numeric",
           weight_mz = "numeric",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           include_precursor = "logical",
           mz_tol_ms1 = "numeric",
           mz_tol_ms2 = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam")

setClass("AlignPeakParamRaw",
         slots = c(
           mz_tol = "numeric",
           strict_rt_constrains = "logical",
           rt_tol_landmark = "numeric",
           rt_tol_match = "numeric",
           ccs_tol = "numeric",
           cutoff = "numeric",
           min_num_fragments = "numeric",
           method_match = "character",
           method_score = "character",
           weight_intensity = "numeric",
           weight_mz = "numeric",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           include_precursor = "logical",
           mz_tol_ms1 = "numeric",
           mz_tol_ms2 = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam")

setClass("DereplicatePeaksParam",
         slots = c(
           order_column = "character",
           mz_tol = "numeric",
           rt_tol = "numeric",
           mobility_tol = "numeric",
           match_msms = "logical",
           msms_cutoff = "numeric",
           method_match = "character",
           method_score = "character",
           weight_intensity = "numeric",
           weight_mz = "numeric",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           include_precursor = "logical",
           mz_tol_ms1 = "numeric",
           mz_tol_ms2 = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam")

setClass("ExtractIMDataParam",
         slots = c(
           order_column = "character",
           mz_tol = "numeric",
           frame_range = "numeric",
           frame_integration_range = "numeric",
           mobility_range = "numeric",
           mobility_intgration_range = "numeric",
           min_points = "numeric",
           min_intensity = "numeric",
           n_skip = "numeric",
           interpolate_method = "nullOrCharacter",
           smooth_method = "character",
           snthreshold = "numeric",
           peak_span_eim = "numeric",
           peak_span_eic = "numeric",
           smooth_window_eim = "numeric",
           smooth_window_eic = "numeric",
           skip_invalid_eic_peaks = "logical",
           skip_invalid_eim_peaks = "logical",
           keep_profile = "logical",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("QueryTimsDataParam",
         slots = c(
           opentims_thread = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("ReadSpectraParam",
         slots = c(
           intensity_from = "character",
           include_precursor = "logical",
           mz_tol = "numeric",
           denoise = "logical",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass('BinPrecursorParam',
         slots = c(
           mz_tol = "numeric",
           rt_tol = "numeric",
           mobility_tol = "numeric",
           distance_cutoff = "numeric",
           weight_rt = "numeric",
           weight_mobility = "numeric",
           weight_msms = "numeric",
           method_match = "character",
           method_score = "character",
           weight_intensity = "numeric",
           weight_mz = "numeric",
           int_threshold_abs = "numeric",
           int_threshold_rel = "numeric",
           include_precursor = "logical",
           mz_tol_ms1 = "numeric",
           mz_tol_ms2 = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("SearchParam",
         slots = c(ppm = "numeric",
                   scoreRT = "logical",
                   toleranceRT = "numeric",
                   scoreCCS = "logical",
                   toleranceCCS = "numeric",
                   typeCCS = "character",
                   adductIncluded = "nullOrCharacter",
                   adductExcluded = "nullOrCharacter",
                   adductFile = "nullOrCharacter",
                   classIncluded = "nullOrCharacter",
                   classExcluded = "nullOrCharacter",
                   useMS1ResDefine = 'logical',
                   updateRefMZ = 'logical',
                   resDefineAt = "numeric"
         ),
         contains = c("ParamSpectraTools")
)


setClass("MatchParam",
         slots = c(ppm = "numeric",
                   cutoff = "numeric",
                   methodMatch = "character",
                   methodScore = "character",
                   weightMZ = "numeric",
                   weightIntensity = "numeric",
                   includePrecursor = "logical",
                   ppmPrecursorFilter = "numeric",
                   ms2range = "nullOrNumeric",
                   thrIntensityAbs = "nullOrNumeric",
                   thrIntensityRel = "nullOrNumeric",
                   intensityExpNormed = "logical",
                   intensityLibNormed = "logical",
                   tuneLibSpectra = "logical",
                   useMS1ResDefine = 'logical',
                   resDefineAt = "numeric",
                   normIntensity = 'logical',
                   intensityNormedMethod = 'character'
         ),
         contains = c("ParamSpectraTools")
)


setClass('CombineParam',
         slots = c(
           cutoff = "numeric",
           weightRT = "numeric",
           weightCCS = "numeric",
           weightMSMS = "numeric",
           scoreMSMS = "character"
         ),
         contains = "MetCellParam"
)

setClass("QueryTimsDataSegmentParam",
         slots = c(
           opentims_thread = "numeric",
           time_limit_table = 'character',
           pool_size = "integer",
           overlaps = "integer",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("TrimSCTimsDataParam",
         slots = c(
           time_limit_table = 'character',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("DiscoverSCEventsParam",
         slots = c(
           marker_mz = 'numeric',
           marker_mobility = 'numeric',
           mz_tolerance_ppm = 'numeric',
           mobility_range = 'numeric',
           intensity_abs_threshold_upper = 'nullOrNumeric',
           intensity_abs_threshold_lower = 'nullOrNumeric',
           marker_eic_peak_span = 'numeric'
         ),
         contains = "MetCellParam"
)

# setClass("ExtractSCIMDataParam",
#          slots = c(
#            mz_tol = 'numeric',
#            mobility_tol = 'numeric',
#            min_intensity = 'numeric',
#            min_points = 'numeric',
#            n_skip = 'numeric',
#            methodRoi = 'character',
#            methodBaseline = 'character',
#            methodSmooth = 'character',
#            peakwidth = 'numeric',
#            snthr = 'numeric',
#            prefilter = 'numeric',
#            fitgauss = 'logical',
#            rerun = "logical"
#          ),
#          contains = "Met4DXParam"
# )


setClass("ExtractSCIMDataParam",
         slots = c(
           mz_tol = 'numeric',
           mobility_tol = 'numeric',
           min_intensity = 'numeric',
           min_points = 'numeric',
           n_skip = 'numeric',
           smooth_window = 'numeric',
           peak_span_eim = 'numeric',
           snthreshold = 'numeric',
           skip_invalid_peaks = 'logical',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("GroupSCDensityParam",
         slots = c(
           bw = 'numeric',
           mz_bin_size = 'numeric',
           max_features = 'numeric',
           minFraction = 'numeric',
           minSamples = 'numeric',
           plot_density = 'logical',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


# setClass("FillSCPeakParam",
#          slots = c(
#            mz_tol = "numeric",
#            frame_integration_range = "numeric",
#            mobility_intgration_range = "numeric",
#            rerun = "logical"
#          ),
#          contains = "MetCellParam"
# )


setClass("SCMatchParam",
         slots = c(
           typeCCS = "character",
           toleranceCCS = "numeric",
           tolerancemz = "numeric",
           adduct_for_id = "nullOrCharacter",
           library = "nullOrCharacter",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

setClass("UniteCellSuperpositionFrameParam",
         slots = c(
           mz_tolerance_combine = "numeric",
           cell_number = 'nullOrNumeric',
           pool_size = "integer",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("SearchPeakTargetParam",
         slots = c(
           mz_tol_ppm = "numeric",
           n_skip = 'numeric',
           peak_target_length = 'numeric',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("DetectEIMPeaksParam",
         slots = c(
           smooth_window = "numeric",
           peak_span_eim_detection = "numeric",
           peak_span_eim_integration = "numeric",
           signal_sd_threshold = "numeric",
           single_charge_line_slope = "numeric",
           single_charge_line_intercept = "numeric",
           skip_invalid_peaks = "logical",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)

# IsotopeParam

setClass("IsotopeParam",
         slots = c(
           mz_tol_ppm = "numeric",
           mobility_tol = "numeric",
           isotope_delta = "numeric",
           isotope_max_num = "numeric",
           isotope_int_ratio_check = "logical",
           isotope_int_ratio_cutoff = "numeric",
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


setClass("ExtractIndividualCellParam",
         slots = c(
           mz_tol_ppm = "numeric",
           mobility_tol = 'numeric',
           minimal_cell_number = 'numeric',
           pool_size = 'numeric',
           rerun = "logical"
         ),
         contains = "MetCellParam"
)


# setClass("GenerateSCFeatureParam",
#          slots = c(
#            minSamples = 'numeric',
#            minFraction = 'numeric',
#            rerun = "logical"
#          ),
#          contains = "MetCellParam"
# )


########################################################################################################################
## Data Classes
########################################################################################################################

#' @export
setClass("Experiment",
         slots = c(
           wd = "character",
           injection_order = 'nullOrCharacter',
           ion_mode = "character",
           res_dir = "character",
           tmp_dir = "character",
           ce = "character",
           ms1range = "nullOrNumeric",
           ms2range = "nullOrNumeric",
           rt_range = "numeric",
           lc_column = "nullOrCharacter",
           lc_method = "nullOrCharacter",
           res_define_at = "numeric",
           BPPARAM = "BiocParallelParam"
         ),
         contains = "VIRTUAL"
)

#' @export
setClass("DDAExperiment",
         contains = "Experiment"
)

#' @export
setClass("DIAExperiment",
         contains = "Experiment"
)

#' @export
setClass("SingleCellExperiment",
         contains = "Experiment"
)


#' @export
setClass("TimsData",
         slots = c(
           experiment = "Experiment",
           files = "character",
           features = "data.frame",
           peaks = "data.frame",
           filled_peaks = "data.frame",
           peak_groups = "list",
           spectra = "list",
           sample_groups = "data.frame",
           tmp_data_files = "list",
           .processHistory = "list"
         ))
#' @export
setClass("IMData",
         slots = c(
           info = "numeric",
           peak_quality = "list",
           eic = "list",
           eic_rt = "numeric",
           eic_mz = "numeric",
           eim = "list",
           eim_mobility = "numeric",
           eim_mobility_interpolate = "nullOrNumeric",
           profile_data = "data.frame"
         ))
