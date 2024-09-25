# project: MetCell
# File name: zzz.R.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:51
# Copyright (c) 2022- ZhuMSLab ALL right reserved



#' Processing MS data
#'
#' This package processes MS data for peak detection and MSMS spectra assignment
#' to generate feature table and corresponding MSMS spectra.
#' @docType _PACKAGE
#' @author Mingdu Luo(luomd@sioc.ac.cn), Yandong Yin (\email{yinyandong@@sioc.ac.cn})
#' @import SpectraTools opentimsr BiocParallel ggpmisc smoother Rcpp splus2R fastmatch mzR parallel dplyr data.table RcppProgress collapse
#' @importFrom data.table ":="
#' @importFrom Rcpp evalCpp
#' @useDynLib MetCell
#' @name MetCell
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\n========== WELCOME to ",
                        getPackageName(), " v", packageVersion(getPackageName()),
                        " ===========",
                        "\nDeveloped by zhulab for data processing for  IM-resolved mass cytometry",
                        "\nImported packages:",
                        "\n  SpectraTools v", packageVersion("SpectraTools"),
                        "\n  OpenTIMSR v", packageVersion("opentimsr"),
                        "\n================================================\n",
                        "Version 1.0.18 (20230918) \n",
                        "o Package released \n",
                        "o Metabolite library contained 135638 compounds",
                        "\n================================================\n")
}

.onLoad <- function(libname, pkgname) {
  # require(methods)
  .set_package_options(pkgname, 'production')
  # .set_package_options(pkgname, 'development')
  getOption('BioC')[[pkgname]]
}
