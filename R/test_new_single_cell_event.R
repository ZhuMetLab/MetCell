# arker_mz = 734.5736
# marker_mobility = 1.395
# mz_tolerance_ppm = 20
# mobility_range = 0.05
# intensity_abs_threshold_upper = 800000
# intensity_abs_threshold_lower = 50000
# signal_fold_change = 2
#
# temp_intensity <- readRDS('E:/04_single_cell/00_data_processing/20240130_am/20240130_lsec/results/tmp/single_cell_events/16J_LSEC_1.d')
# plot(temp_intensity$frame_index, temp_intensity$marker_intensity, type = 'l')
# temp_intensity$sc <- FALSE
#
#
#
#
# apex <- which(ggpmisc:::find_peaks(temp_intensity$marker_intensity, span = marker_eic_peak_span) &
#                 temp_intensity$marker_intensity >= intensity_abs_threshold_lower &
#                 temp_intensity$marker_intensity <= intensity_abs_threshold_upper)
# length(apex)
# # points()
#
# temp_intensity$sc[apex] <- TRUE
# temp_intensity$sc[1] <- FALSE
# temp_intensity$sc[nrow(temp_intensity)] <- FALSE
#
#
# # ppp <- readRDS(tims_data@tmp_data_files$single_cell_events)
# # table(ppp$sc)
# # plot(ppp$frame_index, ppp$marker_intensity, type = 'l')
# points(temp_intensity$frame_index[which(temp_intensity$sc)],
#        temp_intensity$marker_intensity[which(temp_intensity$sc)],
#        col = 'red')
#
