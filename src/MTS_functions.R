# Helper functions for the initial processing step of the MTS pipeline

#---- Load Packages ----
library(tidyverse)
library(magrittr)
library(data.table)
library(tidyr)
library(scam)
library(dr4pl)
library(sva)
library(readr)
library(stats)
library(reshape2)
library(dplyr)
library(hdf5r)

#---- Reading ----
# HDF5 file reader
read_hdf5 <- function(filename, index = NULL) {
  fun_call <- match.call()
  hdf5_obj <- hdf5r::H5File$new(filename, mode = "r+")
  hdf5_attributes <- hdf5r::h5attributes(hdf5_obj)
  matrix_dims <- hdf5_obj[["0/DATA/0/matrix"]][["dims"]]
  if (is.null(index)) {
    index <- list(1:matrix_dims[1], 1:matrix_dims[2])
  }
  data_matrix <- hdf5_obj[["0/DATA/0/matrix"]][index[[1]],
                                               index[[2]]]
  if (is.null(dim(data_matrix))) {
    data_matrix %<>% matrix(nrow = length(index[[1]]),
                            ncol = length(index[[2]]))
  }
  data_matrix %<>%
    magrittr::set_rownames(hdf5_obj[["0/META/ROW/id"]][index[[1]]] %>%
                             gsub(" *$", "", .)) %>%
    magrittr::set_colnames(hdf5_obj[["0/META/COL/id"]][index[[2]]] %>%
                             gsub(" *$", "", .))
  hdf5_obj$close_all()
  return(data_matrix)
}

# convert a wide platemap to log form
make_long_map <- function(df) {
  pert1 <- df %>%
    dplyr::select(!contains("2"))
  pert2 <- df %>%
    dplyr::select(contains("2"), pert_plate, pert_well, pert_time,
                  prism_replicate, is_well_failure, profile_id)

  colnames(pert2) <- sapply(colnames(pert2), FUN = function(x) rename_col(x))

  if (ncol(pert2) > 6) {
    new_map <- dplyr::bind_rows(pert1, pert2)
  } else {
    pert1  %<>%
      dplyr::filter(pert_iname != "Untrt") %>%
      dplyr::mutate(pert_type = ifelse(pert_iname %in% c("PBS", "DMSO"), "ctl_vehicle", pert_type)) %>%
      dplyr::rename(pert_name = pert_iname)

    if (!("pert_mfc_id") %in% colnames(pert1)){
      pert1 %<>% dplyr::mutate(pert_mfc_id = pert_id)
    }

    return(pert1)
  }

  new_map %<>%
    dplyr::filter(pert_iname != "Untrt") %>%
    dplyr::select(intersect(colnames(.), colnames(pert2))) %>%
    dplyr::mutate(pert_type = ifelse(pert_iname %in% c("PBS", "DMSO"), "ctl_vehicle", pert_type))

  overview <- new_map %>%
    dplyr::group_by(pert_well, pert_plate, prism_replicate, profile_id) %>%
    dplyr::summarize(pert_types = paste(unique(pert_type), collapse = fixed("+")),
                     pert_names = paste(unique(pert_iname), collapse = fixed("+")),
                     n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(new_map) %>%
    dplyr::filter(!(pert_type == "ctl_vehicle" & str_detect(pert_types, "trt")))  %>%
    dplyr::mutate(pert_iname = ifelse(pert_types == "ctl_vehicle", pert_names, pert_iname),
                  pert_vehicle = ifelse(pert_types == "ctl_vehicle", pert_names, pert_vehicle),
                  pert_id = ifelse(pert_types == "ctl_vehicle", pert_names, pert_id),
                  pert_mfc_id = ifelse(pert_types == "ctl_vehicle", pert_names, pert_mfc_id)) %>%
    dplyr::select(-pert_types, -pert_names, -n) %>%
    dplyr::distinct() %>%
    dplyr::rename(pert_name = pert_iname)

  return(overview)
}

#---- Normalization ----
# calculate control barcode medians
control_medians <- function(X) {
  ref <- X %>%
    dplyr::filter(pert_type == "ctl_vehicle") %>%  # look at controls
    dplyr::group_by(prism_replicate, pert_well) %>%
    dplyr::mutate(mLMFI = median(logMFI)) %>%  # median for rep
    dplyr::group_by(prism_replicate, rid) %>%
    dplyr::mutate(mmLMFI = logMFI - mLMFI + median(mLMFI)) %>%  # normalized value for rep
    dplyr::summarize(rLMFI = median(mmLMFI)) %>%  # median normalized value across reps
    dplyr::left_join(X)

  return(ref)
}

# fit scam to control barcode profiles and normalize other data
normalize <- function(X, barcodes) {
  normalized <- X %>%
    dplyr::group_by(prism_replicate, pert_well) %>%
    dplyr::mutate(LMFI = scam(y ~ s(x, bs = "micv", k = 5),
                              data = tibble(
                                y = rLMFI[rid %in% barcodes$rid],
                                x = logMFI[rid %in% barcodes$rid])) %>%
                    predict(newdata = tibble(x = logMFI))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-logMFI)

  return(normalized)
}

#---- SSMD calculations ----
# calculate SSMD and NNMD
calc_ssmd <- function(X) {
  SSMD_table <- X %>%
    # look at controls
    dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon")) %>%
    dplyr::distinct(ccle_name, rid, pert_type, prism_replicate, LMFI, profile_id,
                    pert_time, pool_id, culture) %>%
    # group common controls
    dplyr::group_by(pert_type, prism_replicate, pert_time, ccle_name, rid,
                    pool_id, culture) %>%
    # take median and mad of results
    dplyr::summarize(med = median(LMFI, na.rm = TRUE),
                     mad = mad(LMFI, na.rm = TRUE)) %>%
    # add to table
    dplyr::mutate(pert_type_md = paste0(pert_type, '_md'),
                  pert_type_mad = paste0(pert_type, '_mad')) %>%
    # spread to columns
    tidyr::spread(key = pert_type_md, value = med, fill = 0) %>%
    tidyr::spread(key = pert_type_mad, value = mad, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pert_type) %>%
    # give each control all values (median and mad for vehicle and poscon)
    dplyr::group_by(prism_replicate, ccle_name, pert_time, rid, pool_id, culture) %>%
    dplyr::summarize_all(sum) %>%
    # calculate SSMD and NNMD
    dplyr::mutate(ssmd = (ctl_vehicle_md - trt_poscon_md) /
                    sqrt(ctl_vehicle_mad^2 +trt_poscon_mad^2),
                  nnmd = (ctl_vehicle_md - trt_poscon_md) / ctl_vehicle_mad)

  return(SSMD_table)
}

#---- Dose-Response Parameters ----
# area under curve given dose-response parameters
compute_auc <- function(l, u, ec50, h, md, MD) {
  f1 = function(x) pmax(pmin((l + (u - l)/(1 + (2^x/ec50)^h)),1), 0 )
  return(tryCatch(integrate(f1, log2(md),log2(MD))$value/(log2(MD/md)),
                  error = function(x) NA))
}

# log IC50 from given dose-response parameters
compute_log_ic50 <- function(l, u, ec50, h, md, MD) {
  if((l >= 0.5) | (u <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (l + (u - l)/(1 + (2^x/ec50)^h) - 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,
                    error = function(x) NA))
  }
}

# area over curve given dose-response parameters
compute_aoc <- function(l, u, ec50, h, md, MD) {
  f1 = function(x) 1 - pmax(pmin((l + (u - l)/(1 + (2^x/ec50)^h)), 1), -1)
  return(tryCatch(integrate(f1, log2(md),log2(MD))$value/(log2(MD/md)),
                  error = function(x) NA))
}

# log gr50 from given dose-response parameters
compute_log_gr50 <- function(l, u, ec50, h, md, MD) {
  if((l >= 0.5) | (u <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (l + (u - l)/(1 + (2^x/ec50)^h) - 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,
                    error = function(x) NA))
  }
}

#---- Batch Correction ----
# corrects for pool effects using ComBat
apply_combat <- function(Y) {

  # create "condition" column to be used as "batches"
  df <- Y %>%
    dplyr::distinct(ccle_name, prism_replicate, LFC, culture, pool_id, pert_well) %>%
    tidyr::unite(cond, culture, pool_id, prism_replicate, sep = "::") %>%
    dplyr::filter(is.finite(LFC))

  # calculate means and sd's of each condition
  batch <- df$cond
  m <- rbind(df$LFC,
             rnorm(length(df$LFC),
                   mean =  mean(df$LFC, na.rm = TRUE),
                   sd = sd(df$LFC, na.rm = TRUE)))

  # use ComBat to align means and sd's of conditions
  combat <- sva::ComBat(dat = m, batch = batch) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(ccle_name = df$ccle_name, cond = df$cond, pert_well = df$pert_well) %>%
    dplyr::rename(LFC.cb = V1) %>%
    dplyr::mutate(culture = stringr::word(cond, 1, sep = stringr::fixed("::")),
                  pool_id = stringr::word(cond, 2, sep = stringr::fixed("::")),
                  prism_replicate = stringr::word(cond, 3, sep = stringr::fixed("::"))) %>%
    dplyr::select(-cond, -V2)

  combat_corrected <- Y %>%
    dplyr::left_join(combat) %>%
    .$LFC.cb

  return(combat_corrected)
}
