# project: MetCell
# File name: AllGenerics.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:50
# Copyright (c) 2022- ZhuMSLab ALL right reserved

setGeneric("TimsData",
           function(experiment, ...)
             standardGeneric("TimsData"))

setGeneric("ReadSpectraData",
           function(object, param, ...)
             standardGeneric("ReadSpectraData"))

setGeneric("BinPrecursors",
           function(object, param, ...)
             standardGeneric("BinPrecursors"))

setGeneric("QueryTimsData",
           function(object, param, ...)
             standardGeneric("QueryTimsData"))

setGeneric("ExtractIMData",
           function(object, param, ...)
             standardGeneric("ExtractIMData"))

setGeneric("DereplicatePeaks",
           function(object, param, ...)
             standardGeneric("DereplicatePeaks"))

setGeneric("AlignPeaks",
           function(object, param, ref_index, ref_sample, num_pools, ...)
             standardGeneric("AlignPeaks"))

setGeneric("CorrectRT",
           function(object, param, ref_index, ref_sample, ...)
             standardGeneric("CorrectRT"))

setGeneric("GroupPeaks",
           function(object, param, ...)
             standardGeneric("GroupPeaks"))

setGeneric("MatchBetweenRuns",
           function(object, param, ...)
             standardGeneric("MatchBetweenRuns"))

setGeneric("FinalizeFeatures",
           function(object, param, ...)
             standardGeneric("FinalizeFeatures"))

setGeneric("FillPeaks",
           function(object, param, ...)
             standardGeneric("FillPeaks"))

setGeneric("SmoothData",
           function(object, data_slot, param, ...)
             standardGeneric("SmoothData"))

setGeneric("SetData",
           function(object, data_slot, data_range, ...)
             standardGeneric("SetData"))

setGeneric("SetInfo",
           function(object, apex_eic, apex_eim, frame_integration_range, ...)
             standardGeneric("SetInfo"))

setGeneric("InterpolateEIM",
           function(object, mobility_range, interpolate_method, target_intensity, ...)
             standardGeneric("InterpolateEIM"))

setGeneric("DropProfile",
           function(object, interpolation_method, ...)
             standardGeneric("DropProfile"))

setGeneric("GetPeakInfo",
           function(object, ...)
             standardGeneric("GetPeakInfo"))

setGeneric("IdentifyPeaks",
           function(object, search_param, match_param, combine_param, ...)
             standardGeneric("IdentifyPeaks"))

setGeneric("QueryTimsDataSegment",
           function(object, param, sc_param, ...)
             standardGeneric("QueryTimsDataSegment"))

setGeneric("TrimSCTimsData",
           function(object, param, ...)
             standardGeneric("TrimSCTimsData"))

setGeneric("DiscoverSCEvents",
           function(object, param, ...)
             standardGeneric("DiscoverSCEvents"))

setGeneric("ExtractSCIMData",
           function(object, param, ...)
             standardGeneric("ExtractSCIMData"))

setGeneric("GroupSCPeaks",
           function(object, param, ...)
             standardGeneric("GroupSCPeaks"))

setGeneric("FillSCPeaks",
           function(object, param, ...)
             standardGeneric("FillSCPeaks"))

setGeneric("IdentifySCPeaks",
           function(object, isoparam, param, ...)
             standardGeneric("IdentifySCPeaks"))

setGeneric("UniteCellSuperpositionFrame",
           function(object, param, ...)
             standardGeneric("UniteCellSuperpositionFrame"))

setGeneric("SearchPeakTargets",
           function(object, param,...)
             standardGeneric("SearchPeakTargets"))

setGeneric("DetectEIMPeaks",
           function(object, param,...)
             standardGeneric("DetectEIMPeaks"))

setGeneric("ExtractIndividualCell",
           function(object, param, ...)
             standardGeneric("ExtractIndividualCell"))

# setGeneric("GenerateSCFeature",
#            function(object, param, ...)
#              standardGeneric("GenerateSCFeature"))
