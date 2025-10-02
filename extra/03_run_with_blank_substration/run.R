library(MetCell)
wd <- '.' # the wd should be the same as this
setwd(wd)


#### Experiment setup ##########################################################
exp <- Experiment(wd = wd,
                  nSlaves = 3, # Number of threads for data processing
                  ion_mode = 'positive', # ion mode: "positive" or "negative"
                  experiment_type = 'SingleCell') # data type: "SingleCell" for metabolomics or "SingleCellLipidomics" for lipidomics
tims_data <- TimsData(exp)


#### Single-cell frame detection ###############################################
param <- QueryTimsDataSegmentParam(time_limit_table = 'time_limit_table.csv',
                                   pool_size = 2000L, # how many frames to be processed in one thread
                                   rerun = FALSE)
sc_param <- DiscoverSCEventsParam(marker_mz = 760.5845, # m/z of the maker ion (Da)
                                  marker_mobility = 1.402, # Mobility value of the marker ion (V·s/cm2)
                                  mz_tolerance_ppm = 20, # m/z tolerance to extract and integrate marker ion in each frame (ppm)
                                  mobility_range = 0.05, # Ion mobility tolerance to extract and integrate marker ion in each frame (V·s/cm2)
                                  intensity_abs_threshold_upper = 600000, # The upper intensity limit of detected pulse EIM peaks (counts)
                                  intensity_abs_threshold_lower = 60000, # The lower intensity limit of detected pulse EIM peak (counts)
                                  marker_eic_peak_span = 5 # Peak span to detect pulse EIM peaks of the marker (point)
                                  )
tims_data <- QueryTimsDataSegment(tims_data, param, sc_param)



#### Cell superposition ########################################################
param <- UniteCellSuperpositionFrameParam(mz_tolerance_combine = 20, # m/z tolerance to group and aggregate ions in the same scans across cells (ppm)
                                          cell_number = 600, # Cell number for cell superposition (number)
                                          pool_size = 40L, # how many scans to processed in one thread
                                          rerun =FALSE)
tims_data <- UniteCellSuperpositionFrame(tims_data, param)


#### Peak detection ############################################################
param <- SearchPeakTargetParam(mz_tol_ppm = 20, # m/z tolerance to assemble ion mobilogram in cell superposition frame (ppm)
                               n_skip = 1, # The maximal number of skipped MS1 data point allowed for ion mobilogram assemble
                               peak_target_length = 15, # The minimal number of an ion mobilogram (point)
                               rerun = FALSE)
tims_data <- SearchPeakTargets(tims_data, param)

param <- DetectEIMPeaksParam(smooth_window = 10, # Smooth window of LOESS applied in ion mobilogram during peak detection (point)
                             peak_span_eim_detection = 21, # Peak span to detect EIM apex in ion mobilogram (point)
                             peak_span_eim_integration = 27, # Peak span to EIM integration (point)
                             signal_sd_threshold = 0.35, # The threshold of normalized standard noise
                             skip_invalid_peaks = TRUE,
                             single_charge_line_slope = 0.0009, # slope of the line to distinguish singly- and multiply-charged ions
                             single_charge_line_intercept = 0.35, # intercept of the line to distinguish singly- and multiply-charged ions
                             blank_peak_file = 'eim_peaks', # path to the eim_peaks file of the blank sample
                             mz_tol_match_blank = 5, # m/z tolerance to match peaks between single-cell data and the blank sample data (ppm)
                             res_define_at = 200,
                             ccs_tol_match_blank = 0.01, # ccs tolerance to match peaks between single-cell data and the blank sample data (%)
                             intensity_threshold = 3, # the peak intensity ratio between the single-cell data and the blank sample data, peaks with ratios greater than this values would be kept
                             rerun = FALSE)
tims_data <- DetectEIMPeaks(tims_data, param)


#### Peak quantification #######################################################
param <- ExtractIndividualCellParam(mz_tol_ppm = 20, # m/z tolerance to extract MS1 data points in each single-cell frame (ppm)
                                    mobility_tol = 0.04, # Ion mobility to extract MS1 data points in each single-cell frame (V·s/cm2)
                                    minimal_cell_number = 10, # The minimal cell number of a peak to be kept
                                    pool_size = 150, # how many cells to be extracted in one thread
                                    rerun = FALSE)
tims_data <- ExtractIndividualCell(tims_data, param)


#### Metabolite annotation #####################################################
iso_parame <- IsotopeParam(mz_tol_ppm = 20, # m/z tolerance to annotate isotope peaks (ppm)
                           mobility_tol = 0.01, # Ion mobility tolerance to annotate isotope peaks (V·s/cm2)
                           isotope_delta = 1.003355, # m/z difference to search a isotopic peak (Da)
                           isotope_max_num = 4, # Maximun number of isotopic peak to be annotated: [M] --> [M+3]
                           isotope_int_ratio_check = TRUE,
                           isotope_int_ratio_cutoff = 500,
                           rerun = FALSE)
param <- SCMatchParam(typeCCS = "percentage",
                      toleranceCCS = 6, # CCS tolerance to match with metabolite library (%)
                      tolerancemz = 10, # m/z tolerance to match with metabolite library (ppm)
                      adduct_for_id = c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H2O+H]+'), # adduct forms for putative metabolite annotation: positive mode to select: '[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H2O+H]+'; negative mode to select: '[M-H]-', '[M+Na-2H]-', '[M+HCOO]-'
                      rerun = FALSE)
tims_data <- IdentifySCPeaks(tims_data, iso_parame, param)
