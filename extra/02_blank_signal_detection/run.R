library(MetCell)
wd <- '.' # the wd should be the same as this
setwd(wd)


#### Experiment setup ##########################################################
exp <- Experiment(wd = wd,
                  nSlaves = 7, # Number of threads for data processing
                  ion_mode = 'positive', # ion mode: "positive" or "negative"
                  experiment_type = 'SingleCell')
tims_data <- TimsData(exp)


#### Single-cell frame detection ###############################################
param <- QueryTimsDataSegmentParam(time_limit_table = 'time_limit_table.csv',
                                   pool_size = 2000L, # how many frame to be processed in one thread
                                   rerun = FALSE)
sc_param <- DiscoverSCEventsParam(marker_mz = 158.1532, # m/z of a blank ion in your blank sample (Da)
                                  marker_mobility = 0.673, # Mobility value of the blank ion in your blank sample (V·s/cm2)
                                  mz_tolerance_ppm = 20, # m/z tolerance to extract and integrate blank ion in each frame (ppm)
                                  mobility_range = 0.05, # Ion mobility tolerance to extract and integrate blank ion in each frame (V·s/cm2)
                                  intensity_abs_threshold_upper = 70000, # The upper intensity limit of detected pulse EIM peaks (counts)
                                  intensity_abs_threshold_lower = 30000, # The lower intensity limit of detected pulse EIM peak (counts)
                                  marker_eic_peak_span = 1 # Peak span to detect pulse EIM peaks of the marker (point)
)
tims_data <- QueryTimsDataSegment(tims_data, param, sc_param)



#### Cell superposition ########################################################
param <- UniteCellSuperpositionFrameParam(mz_tolerance_combine = 20, # m/z tolerance to group and aggregate ions in the same scans across cells (ppm)
                                          cell_number = 600, # number of frame for cell superposition, please keep the same of that use in cell samples (default: 600)
                                          pool_size = 40L, # how many scans to processed in one thread
                                          rerun =F)
tims_data <- UniteCellSuperpositionFrame(tims_data, param)


#### Peak detection ############################################################
param <- SearchPeakTargetParam(mz_tol_ppm = 20, # m/z tolerance to assemble ion mobilogram in cell superposition frame (ppm)
                               n_skip = 0, # The maximal number of skipped MS1 data point allowed for ion mobilogram assemble
                               peak_target_length = 25, # The minimal number of an ion mobilogram (point)
                               rerun = FALSE)
tims_data <- SearchPeakTargets(tims_data, param)

param <- DetectEIMPeaksParam(smooth_window = 5, # Smooth window of LOESS applied in ion mobilogram during peak detection (point)
                             peak_span_eim_detection = 21, # Peak span to detect EIM apex in ion mobilogram (point)
                             peak_span_eim_integration = 27, # Peak span to EIM integration (point)
                             signal_sd_threshold = 0.1, # The threshold of normalized standard noise
                             skip_invalid_peaks = TRUE,
                             single_charge_line_slope = 0.0009, # slope of the line to distinguish singly- and multiply-charged ions
                             single_charge_line_intercept = 0.35, # intercept of the line to distinguish singly- and multiply-charged ions
                             rerun = FALSE)
tims_data <- DetectEIMPeaks(tims_data, param)
