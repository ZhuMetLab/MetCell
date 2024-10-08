# project: MetCell
# File name: utils.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:59
# Copyright (c) 2022- ZhuMSLab ALL right reserved

.get_num_cores <- function() {
  opt <- getOption("BioC")[[getPackageName()]]
  if (opt$run_env == "development") {
    return(1)
  }
  if ("num_cores" %in% names(opt)) {
    return(opt$num_cores)
  }
  return(parallel::detectCores() - 2)
}

.gen_bpparam <- function() {
  num_cores <- .get_num_cores()
  os <- Sys.info()["sysname"]
  switch(os,
         "Darwin" = {
           BiocParallel::register(BiocParallel::MulticoreParam(workers = num_cores, progressbar = TRUE))
         },
         "Linux" = {
           BiocParallel::register(BiocParallel::SnowParam(workers = num_cores, progressbar = TRUE))
         },
         "Windows" = {
           BiocParallel::register(BiocParallel::SnowParam(workers = num_cores, progressbar = TRUE))
         }
  )
  return(BiocParallel::bpparam())
}

.gen_split_indexes <- function(num_data, split_size, overlaps = 0) {
  split_size <- split_size - overlaps * 2
  index_from <- seq(from = 1, to = num_data, by = split_size)
  index_to <- c(index_from[-1] - 1, num_data)
  res <- cbind(index_from, index_to)
  if ((num_row <- nrow(res)) > 1) {
    res[2:num_row, 1] <- res[2:num_row, 1] - overlaps
    res[1:(num_row - 1), 2] <- res[1:(num_row - 1), 2] + overlaps
  }
  return(res)
}

.get_ref_lib <- function() {
  sys_info <- Sys.info()
  sys_name <- sys_info["sysname"]
  machine <- sys_info["machine"]
  pkg <- getPackageName()
  if (pkg == ".GlobalEnv") {
    pkg <- "MetCell"
  }
  if (sys_name == "Darwin" || sys_name == "Linux" && machine == "x86_64") {
    ref_lib <- system.file(package = pkg, "lib", "linux64", "libtimsdata.so")
  } else if (sys_name == "Windows") {
    if (machine == "x86-64") {
      ref_lib <- system.file(package = pkg, "lib", "win64", "timsdata.dll")
    } else if (machine == "i686") {
      ref_lib <- system.file(package = pkg, "lib", "win32", "timsdata.dll")
    }
  } else {
    stop("Unsupported OS")
  }
  return(ref_lib)
}

.set_opentims <- function(accept_bruker_eula = TRUE) {
  if (accept_bruker_eula) {
    opentimsr::setup_bruker_so(.get_ref_lib())
    all_columns <- c("frame", "scan", "tof", "intensity", "mz", "inv_ion_mobility", "retention_time")
  } else {
    all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
  }
  return(all_columns)
}

.show_info <- function(info) {
  paste0("*****", info, "*****\n")
}

.show_seperate_lines <- function(info) {
  sepLine <- paste(rep("=", floor(40 - nchar(info) / 2)), collapse = "")
  message(sepLine, info, sepLine)
}

.check_file <- function(fn, dirpath = ".") {
  if (!file.exists(fn)) {
    warnTxt <- paste0("file '", fn,
                      "' does not found! Probaboly the previous data was moved.")
    fn <- list.files(dirpath, pattern = basename(fn), recursive = TRUE)
    message(paste0(warnTxt, "\nUsing '", fn, "' in current working dir instead."))
    if (length(fn) == 0) {
      stop(warnTxt)
    }
  }
  return(fn)
}

.check_skip <- function(fn, rerun, show_message = TRUE) {
  # check if use the existed results. if rerun = FALSE, ignore this check
  a <- file.exists(fn) & (!rerun)
  if (a & show_message) {
    cat("  Using previous results: ", fn, "\n")
  }
  a
}

.load_data <- function(file, keep.name = FALSE, env) {
  if (missing(env)) env <- new.env()
  b <- load(file, envir = env)
  if (keep.name | length(b) > 1) {
    r <- lapply(b, function(b1) env[[ b1]])
    names(r) <- b
    r
  } else {
    env[[b]]
  }
}

.gen_indexes <- function(data, prefix = "#") {
  len = 0
  if (is.data.frame(data) || is.matrix(x)) {
    len <- nrow(data)
  } else if (is.vector(x)) {
    len <- length(data)
  } else {
    stop("Unsupported data type")
  }
  return(paste0(prefix, seq(len)))
}

.parallel_parser <- function(analysis, arg_list, files,
                             BPPARAM = BiocParallel::BiocParallelParam(),
                             save_in_analysis = FALSE, return_data = FALSE) {
  is_find <- sapply(seq_along(arg_list), function(idx) {
    .check_skip(files[idx], arg_list[[idx]]$rerun, show_message = FALSE)
  })

  is_dev <- getOption('BioC')[[getPackageName()]]$run_env == 'development'
  if (all(is_find)) {
    cat('  Using previouse results:\n    ')
    cat(paste0(files, collapse = '    \n    '))
    res <- splus2R::ifelse1(return_data, lapply(files, readRDS), as.list(files))
    cat('\n')
  } else {
    if (is_dev) {
    res <- lapply(seq_along(arg_list), function(idx) {
      data_file <- files[idx]
      arg <- arg_list[[idx]]
      if (!.check_skip(data_file, arg$rerun)) {
        if (save_in_analysis) {
          arg$data_file <- data_file
          res_data <- do.call(analysis, arg)
        } else {
          res_data <- do.call(analysis, arg)
          saveRDS(res_data, file = data_file, version = 2)
        }
      } else if (return_data) {
        res_data <- readRDS(data_file)
      }
      return(splus2R::ifelse1(return_data, res_data, data_file))
    })
    } else {
      res <- BiocParallel::bplapply(seq_along(arg_list), function(idx, files, save_in_analysis) {
        data_file2 <- files[idx]
        arg <- arg_list[[idx]]
        if (!.check_skip(data_file2, arg$rerun)) {
          if (save_in_analysis) {
            arg$data_file2 <- data_file2
            res_data <- do.call(analysis, arg)
          } else {
            res_data <- do.call(analysis, arg)
            saveRDS(res_data, file = data_file2, version = 2)
          }
        } else if(return_data) {
          res_data <- readRDS(data_file2)
        }
        return(splus2R::ifelse1(return_data, res_data, data_file2))
      }, BPPARAM = BPPARAM, files = files, save_in_analysis = save_in_analysis)
    }
  }
  gc()
  if (return_data) {
    names(res) <- names(arg_list)
    return(res)
  } else {
    return(do.call(c, res))
  }
}

.analysis_parser <- function(analysis, arg, data_file, return_data = FALSE) {
  if (!.check_skip(data_file, arg$rerun)) {
    res_data <- do.call(analysis, arg)
    saveRDS(res_data, data_file, version = 2)
  } else if (return_data) {
    res_data <- readRDS(data_file)
  }
  return(splus2R::ifelse1(return_data, res_data, data_file))
}

.tmp_files <- function(files, tmp_dir, sub_dir = NULL) {
  tmp_path <- ifelse(is.null(sub_dir), tmp_dir, file.path(tmp_dir, sub_dir))
  if (!dir.exists(tmp_path)) {
    dir.create(tmp_path, recursive = TRUE)
  }
  file.path(tmp_path, basename(files))
}

.ppm2dalton <- function(mz, ppm, res_define_at) {
  prod(max(res_define_at, mz), ppm, 1e-6)
}

.ppm2dalton2 <- function(mz, ppm, res_define_at) {
  mz[mz <= res_define_at] <- res_define_at
  mz * ppm * 1e-6
}


.percentage2value <- function(value, percentage) {
  value * percentage / 100
}

.get_weighted_mz <- function(mzs, intensities) {
  return(sum(mzs * intensities) / sum(intensities))
}

.get_weight_mz_sint <- function(mzs, intensities) {
  sint <- sum(intensities)
  return(sum(mzs * intensities) / sint, sint)
}

.interpolate_data <- function(data, x_range = NULL, interpolate_method = NULL, xcol = "k0", ycol = "intensity") {
  if (is.null(interpolate_method)) {
    res <- data[!is.na(data[, xcol]), , drop = FALSE]
  } else {
    res <- switch(interpolate_method,
                  "linear" = .interpolate_data_linear(data, x_range, xcol, ycol),
                  "scans" = .interpolate_data_scans(data, xcol, ycol)
    )
  }

  return(res)
}

.interpolate_data_scans <- function(data, xcol = "k0", ycol = "intensity") {
  col_names <- c(xcol, ycol)
  inter <- median(na.omit(diff(data[, xcol])))
  scan_range <- range(as.numeric(rownames(data)))
  scans <- seq(scan_range[1], scan_range[2])
  num_scans <- length(scans)
  res <- data.frame(rep(NA, num_scans), rep(0, num_scans))
  rownames(res) <- scans
  colnames(res) <- col_names
  res[na.omit(match(rownames(data), rownames(res))),] <- data[, col_names]
  for (idx in which(is.na(res[, xcol]))) {
    res[idx, xcol] <- res[idx - 1, xcol] + inter
  }
  return(res)
}

.interpolate_data_linear <- function(data, x_range, xcol, ycol) {
  data <- data[data[, ycol] > 0, , drop = FALSE]
  new_data <- data.frame(seq(x_range[1], x_range[2], by = 0.002), 0)
  colnames(new_data) <- c(xcol, ycol)
  res <- rbind(new_data, data)
  res <- res[order(res[, 1], decreasing = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  return(res)
}

.split_data_list <- function(data_list, data_folder = NULL, split_method = "sequential", pool_size = NULL) {
  len_data <- length(data_list)
  if (is.null(pool_size)) {
    pool_size <- ceiling(len_data /.get_num_cores())
  }
  res <- switch(split_method,
                "sequential" = {
                  apply(.gen_split_indexes(len_data, pool_size), 1, function(dr) {
                    data_list[dr[1]:dr[2]]
                  }, simplify = FALSE)
                },
                "average" = {
                  seq_step <- ceiling(len_data / pool_size)
                  lapply(1:seq_step, function(i) {
                    data_list[seq(i, len_data, seq_step)]
                  })
                })
  if (!is.null(data_folder)) {
    file_path <- file.path(data_folder, paste0("tmp_data", seq_along(res), ".rds"))
    for (i in seq_along(res)) {
      saveRDS(res[[i]], file_path[i], version = 2)
    }
    return(file_path)
  } else {
    return(res)
  }
}
