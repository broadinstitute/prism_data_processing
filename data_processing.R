# Script to go from normalized logMFI values to LFC and DRC calculations.

# import necessary libraries and functions using MTS_functions.R
suppressMessages(source("./src/MTS_functions.R"))

#---- Read arguments ----
script_args <- commandArgs(trailingOnly = TRUE)
if (length(script_args) != 3) {
  stop("Please supply necessary arguments",
       call. = FALSE)
}
base_dir <- script_args[1]
out_dir <- script_args[2]
compound <- script_args[3]

write_name <- stringr::str_replace_all(compound, "[[:punct:]\\s]+", "-")

comp_dir <- paste(out_dir, write_name, sep = fixed("/"))
if (!dir.exists(comp_dir)) {dir.create(comp_dir, recursive = T)}

#---- Load the data ----
print("Loading data and pre-processing")
logMFI_normalized <- data.table::fread(paste0(base_dir, "/logMFI_NORMALIZED.csv"))

# ensure single pert_mfc_id for compound
mfc_ids <- logMFI_normalized %>%
  dplyr::filter(pert_name == compound) %>%
  dplyr::distinct(pert_mfc_id) %>%
  .$pert_mfc_id
if (length(mfc_ids) > 1) {
  mfc_id <- mfc_ids[which.max(nchar(mfc_ids))]
  logMFI_normalized %<>%
    dplyr::mutate(pert_mfc_id = ifelse(pert_name == compound, mfc_id, pert_mfc_id))
}

# filter out base plate (early timepoint)
base_normalized <- logMFI_normalized %>%
  dplyr::filter(str_detect(prism_replicate, "BASE"))
logMFI_normalized %<>%
  dplyr::filter(!str_detect(prism_replicate, "BASE"))

# filter to relevant replicates
plates <- logMFI_normalized %>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::distinct(prism_replicate)
logMFI_normalized %<>% dplyr::filter(prism_replicate %in% plates$prism_replicate)

# read in QC results
SSMD_TABLE <- data.table::fread(paste0(base_dir, "/SSMD_TABLE.csv"))

#---- Compute log-fold changes ----
print("Calculating log-fold changes")
LFC_TABLE <- logMFI_normalized %>%
  # join with SSMD (to filter bad lines)
  dplyr::inner_join(SSMD_TABLE %>%
                      dplyr::distinct(ccle_name, prism_replicate, culture, pass),
                    by = c("prism_replicate", "ccle_name", "culture")) %>%
  dplyr::filter(pass) %>%
  dplyr::group_by(prism_replicate, ccle_name, culture) %>%
  # calculate LFC (LMFI - median(LMFIcontrol))
  dplyr::mutate(LFC = LMFI - median(LMFI[pert_type == "ctl_vehicle"])) %>%
  dplyr::distinct(pert_mfc_id, pert_name, prism_replicate, culture, rid, LFC,
                  pert_type, ccle_name, pert_dose, pert_well, pool_id, pert_time,
                  profile_id, pert_idose) %>%
  dplyr::ungroup()

#---- Correct for pool effects ----
print("ComBat correcting")
LFC_TABLE %<>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::mutate(compound_plate = stringr::word(prism_replicate, 1,
                                               sep = stringr::fixed("_"))) %>%
  tidyr::unite(col = "condition", pert_name, pert_dose, compound_plate, pert_time,
               sep = "::", remove = FALSE) %>%
  split(.$condition) %>%
  purrr::map_dfr(~dplyr::mutate(.x, LFC.cb = apply_combat(.))) %>%
  dplyr::select(-condition)

#---- Compute growth rates ----
# control (base) and DMSO
print("Calculating GR metrics")
CONTROL_GR <- tryCatch(expr = {base_normalized %>%
    dplyr::group_by(pert_time, ccle_name, rid, pool_id, culture) %>%  # no compound to group by
    dplyr::summarize(mLMFI.c = median(LMFI),
                     n.c = n(),
                     var.c = (mad(LMFI)^2/n.c) * pi/2,
                     .groups = "drop") %>%  # n = replicates (~300)
    dplyr::select(-n.c) %>%
    dplyr::ungroup() %>%
    dplyr::rename(pert_base_time = pert_time) %>%
    # join with DMSO
    dplyr::inner_join(logMFI_normalized %>%
                        dplyr::filter(pert_type == "ctl_vehicle") %>%
                        dplyr::mutate(assay_length = as.numeric(str_sub(pert_time, 1, -2))/24) %>%
                        dplyr::group_by(ccle_name, rid, pool_id, culture, pert_time, assay_length) %>%
                        dplyr::summarize(mLMFI.d = median(LMFI),
                                         n.d = n(),
                                         var.d = (mad(LMFI)^2/n.d) * pi/2,
                                         .groups = "drop") %>%
                        dplyr::select(-n.d),
                      by = c("ccle_name", "rid", "pool_id", "culture")) %>%
    dplyr::mutate(base_day_num = as.numeric(str_sub(pert_base_time, 1, -2))/24)
}, error = function(e) {
  return(NA)
})
# treatment
GR_TABLE <- tryCatch(expr = {logMFI_normalized %>%
    # now group by compound
    dplyr::group_by(pert_mfc_id, pert_type, pert_name, pert_dose, pert_idose,
                    ccle_name, rid, pool_id, culture, pert_time) %>%
    dplyr::summarize(mLMFI.t = median(LMFI),
                     n.t = n(),
                     var.t = (mad(LMFI)^2/n.t) * pi/2,
                     .groups = "drop") %>%  # n.t = 3 (replicates)
    dplyr::select(-n.t) %>%
    dplyr::inner_join(CONTROL_GR,
                      by = c("ccle_name", "rid", "pool_id", "culture", "pert_time")) %>%
    dplyr::ungroup()
}, error = function(e) {
  return(tibble())
})

# combined
GR_TABLE <- tryCatch(expr = {GR_TABLE %>%
    # calc control change (DMSO - base)/(t - base day),
    # treatment change (treatment - DMSO)/t - control,
    # use to calc Z (treatment/control) and GR (2^Z - 1)
    dplyr::mutate(control_lfc = (mLMFI.d - mLMFI.c)/(assay_length - base_day_num),
                  treatment_control_lfc = (mLMFI.t - mLMFI.d)/(assay_length),
                  treatment_lfc = treatment_control_lfc + control_lfc,
                  Z = treatment_lfc/control_lfc,
                  var.treatment = (var.t/assay_length^2) + (var.d/(assay_length - base_day_num)^2) +
                    (var.c*(1/(assay_length-base_day_num) - 1/assay_length)^2),
                  var.control = (var.c + var.d)/(assay_length - base_day_num)^2,
                  GR = (2^Z) - 1) %>%
    dplyr::distinct(pert_mfc_id, pert_type, pert_name, pert_dose, pert_idose,
                    rid, ccle_name, culture, pool_id, pert_time, assay_length, base_day_num,
                    control_lfc, treatment_lfc, Z, var.treatment, var.control, GR) %>%
    dplyr::rename(base_day = base_day_num)
}, error = function(e) {
  return(tibble())
})

#---- Compute dose-response parameters ----
# table with each compound cell line combo and number of doses
DRC_TABLE_cb <- LFC_TABLE %>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::distinct(ccle_name, culture, pert_mfc_id, pert_name, pert_dose, pert_time) %>%
  dplyr::count(ccle_name, culture, pert_mfc_id, pert_name, pert_time) %>%
  dplyr::filter(n > 3)  # only fit curves with 4+ doses

if (nrow(DRC_TABLE_cb > 0)) {
  print("Fitting dose-response curves")
  DRC_TABLE_cb %<>% dplyr::mutate(ix = 1:n())
  DRC_cb <- list()  # empty tibble to track results

  # loop through compound cell line combos fitting curves
  for (jx in 1:nrow(DRC_TABLE_cb)) {
    d = DRC_TABLE_cb %>%
      dplyr::filter(ix == jx) %>%
      dplyr::left_join(LFC_TABLE, by = c("ccle_name", "culture", "pert_mfc_id",
                                         "pert_name", "pert_time"))

    # fit curve
    a = tryCatch(dr4pl(dose = d$pert_dose, response = 2^d$LFC.cb,
                       method.init = "logistic", trend = "decreasing"),
                 error = function(e) return(NA))
    # get parameters
    param <- tryCatch(summary(a)$coefficients$Estimate, error = function(e) return(NA))
    if (!is.na(param)) {
      d %<>%
        dplyr::mutate(pred = dr4pl::MeanResponse(pert_dose, param))
      d %<>%
        dplyr::mutate(e = (2^LFC.cb - pred)^2)  # prediction residuals

      mse <- mean(d$e)
      R2 <- 1 - (sum(d$e)/(nrow(d) * var(d$LFC.cb)))

      x <- tibble(ix = jx,
                  min_dose = min(d$pert_dose),
                  max_dose = max(d$pert_dose),
                  upper_limit = param[1],
                  ec50 = param[2],
                  slope = -param[3],
                  lower_limit = param[4],
                  convergence = a$convergence) %>%
        dplyr::mutate(auc = compute_auc(lower_limit,
                                        upper_limit,
                                        ec50,
                                        slope,
                                        min(d$pert_dose),
                                        max(d$pert_dose)),
                      log2.ic50 = compute_log_ic50(lower_limit,
                                                   upper_limit,
                                                   ec50,
                                                   slope,
                                                   min(d$pert_dose),
                                                   max(d$pert_dose)),
                      mse = mse,
                      R2 = R2)
      DRC_cb[[jx]] <- x
    }
  }

  if (length(DRC_cb) > 0) {
    DRC_TABLE_cb <- DRC_cb %>%
      dplyr::bind_rows() %>%
      dplyr::filter(convergence) %>%
      dplyr::left_join(DRC_TABLE_cb, by = c("ix")) %>%
      dplyr::select(-ix, -convergence, -n)
  } else {
    print("Unable to fit any dose-response curves in LFC space")
    DRC_TABLE_cb <- tibble()
  }
} else {
  DRC_TABLE_cb <- tibble()
}

# GROWTH RATE DOSE-RESPONSE
if (nrow(GR_TABLE) > 0) {
  print("Fitting growth dose-response curves")
  DRC_TABLE_growth <- GR_TABLE %>%
    dplyr::distinct(ccle_name, culture, pert_mfc_id, pert_name, pert_dose, pert_time) %>%
    dplyr::count(ccle_name, culture, pert_mfc_id, pert_name, pert_time) %>%
    dplyr::filter(n > 3)

  if (nrow(DRC_TABLE_growth > 0)) {
    DRC_TABLE_growth %<>% dplyr::mutate(ix = 1:n())
    DRC_gr <- list()

    for (jx in 1:nrow(DRC_TABLE_growth)) {
      d = DRC_TABLE_growth %>%
        dplyr::filter(ix == jx) %>%
        dplyr::left_join(GR_TABLE, by = c("ccle_name", "culture", "pert_mfc_id",
                                          "pert_name", "pert_time"))

      a = tryCatch(dr4pl(dose = d$pert_dose, response = d$GR,
                         method.init = "logistic", trend = "decreasing"),
                   error = function(e) return(NA))
      param <- tryCatch(summary(a)$coefficients$Estimate, error = function(e) return(NA))
      if (!is.na(param)) {
        d %<>%
          dplyr::mutate(pred = dr4pl::MeanResponse(pert_dose, param))
        d %<>%
          dplyr::mutate(e = (GR - pred)^2)  # prediction residuals

        mse <- mean(d$e)
        R2 <- 1 - (sum(d$e)/(nrow(d) * var(d$GR)))

        x <- tibble(ix = jx,
                    min_dose = min(d$pert_dose),
                    max_dose = max(d$pert_dose),
                    upper_limit = param[1],
                    ec50 = param[2],
                    slope = -param[3],
                    lower_limit = param[4],
                    convergence = a$convergence) %>%
          dplyr::mutate(aoc = compute_aoc(lower_limit,
                                          upper_limit,
                                          ec50,
                                          slope,
                                          min(d$pert_dose),
                                          max(d$pert_dose)),
                        log2.gr50 = compute_log_gr50(lower_limit,
                                                     upper_limit,
                                                     ec50,
                                                     slope,
                                                     min(d$pert_dose),
                                                     max(d$pert_dose)),
                        mse = mse,
                        R2 = R2)
        DRC_gr[[jx]] <- x
      }
    }

    if (length(DRC_gr) > 0) {
      DRC_TABLE_growth <- DRC_gr %>%
        dplyr::bind_rows() %>%
        dplyr::filter(convergence) %>%
        dplyr::left_join(DRC_TABLE_growth, by = c("ix")) %>%
        dplyr::select(-ix, -convergence, -n)
    } else {
      print("Unable to fit any dose-response curves in GR space")
      DRC_TABLE_growth <- tibble()
    }
  } else {
    DRC_TABLE_growth <- tibble()
  }
}

#---- Make collapsed LFC table ----

LFC_COLLAPSED_TABLE <- LFC_TABLE %>%
  dplyr::mutate(compound_plate = stringr::word(prism_replicate, 1,
                                               sep = stringr::fixed("_"))) %>%
  dplyr::group_by(ccle_name, culture, pool_id, pert_name, pert_mfc_id,
                  pert_dose, pert_idose, compound_plate, pert_time) %>%
  # LFC and LFC.cb values will be medians across replicates
  dplyr::summarize(LFC = median(LFC, na.rm = TRUE),
                   LFC.cb = median(LFC.cb, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::ungroup()


#---- Write to .csv ----
# compound data (logMFI, DRC, LFC)
logMFI_normalized %>%
  dplyr::bind_rows(base_normalized) %>%
  readr::write_csv(., paste0(comp_dir, "/logMFI_normalized.csv"))
readr::write_csv(LFC_TABLE, paste0(comp_dir, "/LFC_TABLE.csv"))
readr::write_csv(LFC_COLLAPSED_TABLE, paste0(comp_dir, "/LFC_COLLAPSED_TABLE.csv"))

if(nrow(DRC_TABLE_cb) > 0)  {
  readr::write_csv(DRC_TABLE_cb, paste0(comp_dir, "/DRC_TABLE.csv"))
}

# GR data if it exists
if(nrow(GR_TABLE) > 0) {
  readr::write_csv(GR_TABLE, paste0(comp_dir, "/GR_TABLE.csv"))
  if(nrow(DRC_TABLE_growth) > 0) {
    readr::write_csv(DRC_TABLE_growth, paste0(comp_dir, "/DRC_TABLE_GR.csv"))
  }
}


#---- Generate DRC plots ----
if(nrow(DRC_TABLE_cb) > 0) {
  print("Generating plot PDFs")
  # dose response data
  DRC_TABLE_cb %<>% dplyr::arrange(auc)

  # create .pdf
  pdf(paste0(comp_dir, "/", toupper(write_name), "_DRCfigures.pdf"))

  # loop through each cell line treated by compound and plot DRC
  conditions <- DRC_TABLE_cb %>% dplyr::distinct(ccle_name, culture, pert_time)
  for(j in 1:nrow(conditions)) {
    condition <- conditions[j,]
    assay_time <- condition$pert_time
    cell_line <- condition$ccle_name
    culture <- condition$culture

    d <- DRC_TABLE_cb %>% dplyr::inner_join(condition,
                                            by = c("ccle_name", "culture", "pert_time"))
    d_cult_line <- LFC_TABLE %>% dplyr::inner_join(condition,
                                                   by = c("ccle_name", "culture", "pert_time"))

    # DRC curve function
    f1 = function(x) {
      d$lower_limit + (d$upper_limit - d$lower_limit)/
        (1 + (2^x/d$ec50)^d$slope)
    }

    # sequence for plotting curve
    xx = seq(min(log2(d_cult_line$pert_dose)),
             max(log2(d_cult_line$pert_dose)),
             length.out = 1000)

    # plot individual data points and DRC fit line
    p = d_cult_line %>%
      ggplot() +
      geom_point(aes(x = log2(pert_dose),
                     color = prism_replicate, y = 2^LFC.cb)) +
      geom_line(data = tibble(x = xx, y = f1(xx)),
                aes(x = x, y = y, group = 1),  lwd =1 ) +
      ylim(0,2) + theme_bw() +
      labs(x = 'log2(Dose) (uM)', y = 'Viability', color = "",
           title = paste0(compound, "-", assay_time, "\n", cell_line,' - ', culture,
                          "\nAUC:", round(d$auc, 2),
                          " - IC50:", round(2^d$log2.ic50, 2)))
    # outputs to .pdf
    print(p)
  }
  # closes .pdf
  invisible(dev.off())
}
