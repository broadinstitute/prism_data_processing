# Helper functions for the initial processing step of the MTS pipeline

#---- Load Packages ----
library(plyr)
library(tidyverse)
library(magrittr)
library(data.table)
library(tidyr)
library(scam)
library(dr4pl)
library(sva)
library(stats)
library(reshape2)
library(dplyr)

#---- Normalization ----
# calculate control barcode medians
control_medians <- function(X, control) {
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
    dplyr::group_by(prism_replicate, pert_well)%>%
    dplyr::mutate(LMFI = scam(y ~ s(x, bs = "mpi"),
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
    dplyr::distinct(ccle_name, pert_type, prism_replicate,LMFI, profile_id,
                    pool_id, culture) %>%
    # group common controls
    dplyr::group_by(pert_type, prism_replicate, ccle_name, pool_id, culture) %>%
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
    dplyr::group_by(prism_replicate, ccle_name, pool_id, culture) %>%
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
  df <- Y %>%
    dplyr::distinct(ccle_name, prism_replicate, LFC, culture, pool_id) %>%
    tidyr::unite(cond, culture, pool_id, prism_replicate, sep = "::") %>%
    dplyr::filter(is.finite(LFC))
  
  batch <- df$cond
  m <- rbind(df$LFC,
             rnorm(length(df$LFC),
                   mean =  mean(df$LFC, na.rm = TRUE),
                   sd = sd(df$LFC, na.rm = TRUE)))
  
  combat <- sva::ComBat(dat = m, batch = batch) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(ccle_name = df$ccle_name, cond = df$cond) %>%
    dplyr::rename(LFC.cb = V1) %>%
    dplyr::mutate(culture = stringr::word(cond, 1, sep = stringr::fixed("::")),
                  pool_id = stringr::word(cond, 2, sep = stringr::fixed("::")),
                  prism_replicate = stringr::word(cond, 3, sep = stringr::fixed("::"))) %>%
    dplyr::select(-cond, -V2)
  
  Y %>%
    dplyr::left_join(combat) %>%
    .$LFC.cb
}
