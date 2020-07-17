# MTS data processing pipeline (from logMFI to dose-response and viability)
# See MTS_Pipeline.md for details
# Replace data_path with appropriate file to run

# import necessary libraries and functions using MTS_functions.R
suppressMessages(source("./src/MTS_functions.R"))

#---- Load Data ----
data_path <- # INSERT path to logMFI file
unprocessed <- data.table::fread(data_path)

print("Loaded the data")
assay_length <- 5  # change to length of experiment
base_day <- 1  # change to day of early time point measurements

# split into 300 and 500
PR300 <- unprocessed %>%
  dplyr::filter(culture == "PR300")
PR500 <- unprocessed %>%
  dplyr::filter(culture == "PR500")

# create barcode tables
PR300_barcodes <- PR300 %>%
  dplyr::filter(pool_id == "CTLBC")
PR500_barcodes <- PR500 %>%
  dplyr::filter(pool_id == "CTLBC")

# filter out base plate (day 0/1, for use later)
PR300_base <- PR300 %>%
  dplyr::filter(str_detect(prism_replicate, "BASE"))
PR300 %<>%
  dplyr::filter(!str_detect(prism_replicate, "BASE"))
PR500_base <- PR500 %>%
  dplyr::filter(str_detect(prism_replicate, "BASE"))
PR500 %<>%
  dplyr::filter(!str_detect(prism_replicate, "BASE"))


#---- Basic QC and Normalization ----
# compute control barcode median of medians for normalization
PR300_control_medians <- control_medians(PR300)

# fit curve to controls and predict test conditions
PR300_normalized <- normalize(PR300_control_medians, PR300_barcodes)

if(nrow(PR300_base) > 0) {
  # generate reference profile to normalize base data
  PR300_profile <- PR300_normalized %>%
    dplyr::filter(rid %in% PR300_barcodes$rid) %>%
    dplyr::group_by(rid) %>%
    dplyr::mutate(rLMFI = mean(rLMFI)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(rid, rLMFI)
  
  PR300_base_normalized <- PR300_base %>%
    dplyr::left_join(PR300_profile) %>%
    normalize(., PR300_barcodes)
}
# join with other info (LMFI is normalized, logMFI is not)
PR300_normalized %<>%
  dplyr::left_join(PR300) %>%
  dplyr::select(-logMFI)

# REPEAT with 500
PR500_control_medians <- control_medians(PR500)

PR500_normalized <- normalize(PR500_control_medians, PR500_barcodes)

if(nrow(PR500_base) > 0) {
  PR500_profile <- PR500_normalized %>%
    dplyr::filter(rid %in% PR500_barcodes$rid) %>%
    dplyr::group_by(rid) %>%
    dplyr::mutate(rLMFI = mean(rLMFI)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(rid, rLMFI)
  
  PR500_base_normalized <- PR500_base %>%
    dplyr::left_join(PR500_profile) %>%
    normalize(., PR500_barcodes)
}
PR500_normalized %<>%
  dplyr::left_join(PR500) %>%
  dplyr::select(-logMFI)

#---- Generate SSMD table ----
# calculate SSMD and NNMD
SSMD_table_300 <- calc_ssmd(PR300_normalized)

# calculate error rate of normalized table (based on threshold classifier)
PR300_error <- PR300_normalized %>%
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon"),
                is.finite(LMFI)) %>%
  dplyr::group_by(rid, ccle_name ,prism_replicate) %>%
  dplyr::summarize(error_rate =
                     min(PRROC::roc.curve(LMFI[pert_type == "ctl_vehicle"],
                                          LMFI[pert_type == "trt_poscon"],
                                          curve = TRUE )$curve[,1] + 1 -
                           PRROC::roc.curve(LMFI[pert_type == "ctl_vehicle"],
                                            LMFI[pert_type == "trt_poscon"],
                                            curve = TRUE )$curve[,2])/2)

# join with SSMD table
SSMD_table_300 %<>%
  dplyr::left_join(PR300_error)

# REPEAT with 500
SSMD_table_500 <- calc_ssmd(PR500_normalized)

PR500_error <- PR500_normalized %>%
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon"),
                is.finite(LMFI)) %>%
  dplyr::group_by(rid, ccle_name ,prism_replicate) %>%
  dplyr::summarize(error_rate =
                     min(PRROC::roc.curve(LMFI[pert_type == "ctl_vehicle"],
                                          LMFI[pert_type == "trt_poscon"],
                                          curve = TRUE )$curve[,1] + 1 -
                           PRROC::roc.curve(LMFI[pert_type == "ctl_vehicle"],
                                            LMFI[pert_type == "trt_poscon"],
                                            curve = TRUE )$curve[,2])/2)

SSMD_table_500 %<>%
  dplyr::left_join(PR500_error)

# combine 300 and 500 tables
SSMD_TABLE <- dplyr::bind_rows(SSMD_table_500, SSMD_table_300) %>%
  # if error rate <= .05 then pass
  dplyr::mutate(pass = error_rate <= 0.05,
                compound_plate =  stringr::word(prism_replicate, 1,2,
                                                sep = stringr::fixed("_"))) %>%
  dplyr::filter(pool_id != "CTLBC") %>%
  dplyr::group_by(ccle_name, culture, compound_plate) %>%
  # sum how many passes per cell line
  dplyr::mutate(n.rep = sum(pass, na.rm = TRUE)) %>%
  dplyr::ungroup()


#---- Compute log-fold changes ----
LFC_TABLE <- PR300_normalized %>%
  # combine tables
  dplyr::bind_rows(PR500_normalized) %>%
  # join with SSMD
  dplyr::inner_join(SSMD_TABLE %>%
                      dplyr::distinct(ccle_name, prism_replicate,
                                      culture, pass, n.rep)) %>%
  # more than 1 passing sample
  dplyr::filter(pass, n.rep > 1) %>%
  dplyr::group_by(prism_replicate, ccle_name, culture) %>%
  # calculate LFC (LMFI - median(LMFIcontrol))
  dplyr::mutate(LFC = LMFI - median(LMFI[pert_type == "ctl_vehicle"])) %>%
  dplyr::distinct(pert_mfc_id, pert_name, prism_replicate, culture, rid, LFC,
                  pert_type, ccle_name, pert_dose, pert_well, pool_id,
                  profile_id, pert_idose) %>%
  dplyr::ungroup()


#---- Correct for pool effects ----
LFC_TABLE %<>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::mutate(compound_plate = stringr::word(prism_replicate, 1,2,
                                               sep = stringr::fixed("_"))) %>%
  tidyr::unite(col = "condition", pert_mfc_id, pert_dose, compound_plate,
               sep = "::", remove = FALSE) %>%
  split(.$condition) %>%
  purrr::map_dfr(~dplyr::mutate(.x, LFC.cb = apply_combat(.))) %>%
  dplyr::select(-condition)


#---- Compute growth rates ----
# control (base) and DMSO
CONTROL_GR_300 <- tryCatch(expr = {PR300_base_normalized %>%
    dplyr::group_by(ccle_name, rid, pool_id, culture) %>%  # no compound to group by
    dplyr::summarize(mLMFI.c = median(LMFI),
                     n.c = n(),
                     var.c = (mad(LMFI)^2/n.c) * pi/2) %>%  # n = replicates (~300)
    dplyr::select(-n.c) %>%
    # join with DMSO
    dplyr::inner_join(PR300_normalized %>%
                        dplyr::filter(pert_type == "ctl_vehicle") %>%
                        dplyr::group_by(ccle_name, rid, pool_id, culture) %>%
                        dplyr::summarize(mLMFI.d = median(LMFI),
                                         n.d = n(),
                                         var.d = (mad(LMFI)^2/n.d) * pi/2) %>%
                        dplyr::select(-n.d)) %>%
    dplyr::ungroup()
}, error = function(e) {
  return(NA)
})
# repeat with PR500
CONTROL_GR_500 <- tryCatch(expr = {PR500_base_normalized %>%
    dplyr::group_by(ccle_name, rid, pool_id, culture) %>%
    dplyr::summarize(mLMFI.c = median(LMFI),
                     n.c = n(),
                     var.c = (mad(LMFI)^2/n.c) * pi/2) %>%
    dplyr::select(-n.c) %>%
    dplyr::inner_join(PR500_normalized %>%
                        dplyr::filter(pert_type == "ctl_vehicle") %>%
                        dplyr::group_by(ccle_name, rid, pool_id, culture) %>%
                        dplyr::summarize(mLMFI.d = median(LMFI),
                                         n.d = n(),
                                         var.d = (mad(LMFI)^2/n.d) * pi/2) %>%
                        dplyr::select(-n.d)) %>%
    dplyr::ungroup()
}, error = function(e) {
  return(NA)
})
# treatment
GR_300 <- tryCatch(expr = {PR300_normalized %>%
    dplyr::filter(pool_id != "CTLBC") %>% # no control barcodes
    # now group by compound
    dplyr::group_by(pert_mfc_id, pert_type, pert_name, pert_dose, pert_idose,
                    ccle_name, rid, pool_id, culture) %>%
    dplyr::summarize(mLMFI.t = median(LMFI),
                     n.t = n(),
                     var.t = (mad(LMFI)^2/n.t) * pi/2) %>%  # n.t = 3 (replicates)
    dplyr::select(-n.t) %>%
    dplyr::inner_join(CONTROL_GR_300) %>%
    dplyr::ungroup()
}, error = function(e) {
  return(tibble())
})
GR_500 <- tryCatch(expr = {PR500_normalized %>%
    dplyr::filter(pool_id != "CTLBC") %>%
    dplyr::group_by(pert_mfc_id, pert_type, pert_name, pert_dose, pert_idose,
                    ccle_name, rid, pool_id, culture) %>%
    dplyr::summarize(mLMFI.t = median(LMFI),
                     n.t = n(),
                     var.t = (mad(LMFI)^2/n.t) * pi/2) %>%
    dplyr::select(-n.t) %>%
    dplyr::inner_join(CONTROL_GR_500) %>%
    dplyr::ungroup()
}, error = function(e) {
  return(tibble())
})
# combined
GR_TABLE <- tryCatch(expr = {dplyr::bind_rows(GR_300, GR_500) %>%
    # calc control change (DMSO - base)/(t - base day),
    # treatment change (treatment - DMSO)/t - control,
    # use to calc Z (treatment/control) and GR (2^Z - 1)
    dplyr::mutate(control_lfc = (mLMFI.d - mLMFI.c)/(assay_length - base_day),
                  treatment_control_lfc = (mLMFI.t - mLMFI.d)/(assay_length),
                  treatment_lfc = treatment_control_lfc + control_lfc,
                  Z = treatment_lfc/control_lfc,
                  var.treatment = (var.t/assay_length^2) + (var.d/(assay_length - base_day)^2) +
                    (var.c*(1/(assay_length-base_day) - 1/assay_length)^2),
                  var.control = (var.c + var.d)/(assay_length - base_day)^2,
                  GR = (2^Z) - 1) %>%
    dplyr::distinct(pert_mfc_id, pert_type, pert_name, pert_dose, pert_idose, ccle_name, culture, pool_id,
                    control_lfc, treatment_lfc, Z, var.treatment, var.control, GR) %>%
    dplyr::filter(control_lfc > 0)
}, error = function(e) {
  return(tibble())
})

#---- Compute dose-response parameters ----
# table with each compound cell line combo and number of doses
DRC_TABLE_cb <- LFC_TABLE %>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::distinct(ccle_name, culture, pert_mfc_id, pert_name, pert_dose) %>%
  dplyr::count(ccle_name, culture, pert_mfc_id, pert_name) %>%
  dplyr::filter(n > 4) %>%  # only fit curves with 4+ doses
  dplyr::mutate(ix = 1:n())

DRC_cb <- tibble()  # empty tibble to track results

# loop through compound cell line combos fitting curves
for(jx in 1:nrow(DRC_TABLE_cb)) {
  d = DRC_TABLE_cb %>%
    dplyr::filter(ix == jx) %>%
    dplyr::left_join(LFC_TABLE)
  
  # fit curve
  a = tryCatch(dr4pl(dose = d$pert_dose, response = 2^d$LFC.cb,
                     method.init = "logistic", trend = "decreasing"),
               error = function(e) NA)
  
  # get parameters
  param <- tryCatch(summary(a)$coefficients$Estimate, error = function(e) NA)
  if(!is.na(param)) {
    d %<>%
      dplyr::mutate(pred = dr4pl::MeanResponse(pert_dose, param),
                    e = (2^LFC.cb - pred)^2)
    mse <- mean(d$e)
    R2 <- 1 - (sum(d$e)/(nrow(d) * var(d$LFC.cb)))
    auc <- compute_auc(param[4], param[1], param[2], -param[3],
                       min(d$pert_dose),max(d$pert_dose))
    log2.ic50 <- compute_log_ic50(param[4], param[1], param[2], -param[3],
                                  min(d$pert_dose), max(d$pert_dose))
    
    x <- tibble(ix = jx, 
                min_dose = min(d$pert_dose), max_dose = max(d$pert_dose),
                upper_limit = param[1], lower_limit = param[4],
                ec50 = param[2], slope = -param[3],
                convergence = a$convergence,
                auc = auc, log2.ic50 = log2.ic50, mse = mse, R2 = R2)
    DRC_cb %<>% dplyr::bind_rows(x)
  }
}

DRC_TABLE_cb <- DRC_cb %>%
  dplyr::filter(convergence) %>%
  dplyr::left_join(DRC_TABLE_cb) %>%
  dplyr::select(-ix, -convergence, -n)

#---- Make collapsed LFC table ----

LFC_COLLAPSED_TABLE <- LFC_TABLE %>%
  dplyr::mutate(compound_plate = stringr::word(prism_replicate, 1,2,
                                               sep = stringr::fixed("_"))) %>%
  dplyr::group_by(ccle_name, culture, pert_name, pert_mfc_id,
                  pert_dose, pert_idose, compound_plate) %>%
  # LFC and LFC.cb values will be medains across replicates
  dplyr::summarize(LFC = median(LFC, na.rm = TRUE),
                   LFC.cb = median(LFC.cb, na.rm = TRUE))


#---- Write to .csv ----
# SSMD table
readr::write_csv(SSMD_TABLE, paste0(dirname(data_path), "/SSMD_TABLE.csv"))

# normalized controls
logMFI_normalized <- dplyr::bind_rows(PR500_normalized, PR300_normalized) %>%
  dplyr::filter(pool_id != "CTLBC") %>%
  dplyr::select(-rLMFI, -rid)

readr::write_csv(logMFI_normalized,
                 paste0(dirname(data_path), "/logMFI_NORMALIZED.csv"))

# compound data (DRC, LFC)
compounds <- dplyr::distinct(LFC_TABLE, pert_name, pert_mfc_id)
for(i in 1:nrow(compounds)) {
  id <- compounds[[i, "pert_mfc_id"]]  # Broad ID (unique)
  name <- compounds[[i, "pert_name"]]  # name (human readable)
  name <- stringr::str_replace_all(name, "[[:punct:]\\s]+", "-")
  
  # output directory
  path <- paste0(dirname(data_path), "/", name)
  if(!dir.exists(path)) {
    dir.create(path)
  }
  
  lfc <- dplyr::filter(LFC_TABLE, pert_mfc_id == id)
  drc <- dplyr::filter(DRC_TABLE_cb, pert_mfc_id == id)
  lfc_coll <-
    dplyr::filter(LFC_COLLAPSED_TABLE, pert_mfc_id == id)
  
  readr::write_csv(lfc, paste0(path, "/LFC_TABLE.csv"))
  readr::write_csv(lfc_coll, paste0(path, "/LFC_COLLAPSED_TABLE.csv"))
  readr::write_csv(drc, paste0(path, "/DRC_TABLE.csv"))
  
  # GR data if it exists
  if(nrow(GR_TABLE) > 0) {
    gr <- dplyr::filter(GR_TABLE, pert_mfc_id == id)
    # drc_gr <- dplyr::filter(DRC_TABLE_growth, pert_mfc_id == id)
    readr::write_csv(gr, paste0(path, "/GR_TABLE.csv"))
    # readr::write_csv(drc_gr, paste0(path, "/DRC_TABLE_GR.csv"))
  }
}

#---- Generate DRC plots ----
# generate a .pdf of graphs for each compound
for(i in 1:nrow(compounds)) {
  id <- compounds[[i, "pert_mfc_id"]]
  name <- compounds[[i, "pert_name"]]
  name <- stringr::str_replace_all(name, "[[:punct:]\\s]+", "-")
  
  # filter to just see that compound
  compound_DRC <- DRC_TABLE_cb %>%
    dplyr::filter(pert_mfc_id == id) %>%
    dplyr::arrange(desc(auc))
  
  # tracks LFC info
  compound_LFC <- LFC_TABLE %>%
    dplyr::filter(pert_mfc_id == id)
  
  # create .pdf
  pdf(paste0(dirname(data_path), "/", name, "/",
             toupper(name), "_DRCfigures.pdf"))
  
  # loop through each cell line treated by compound and plot DRC
  cell_lines <- compound_DRC$ccle_name %>% unique()
  for(cell_line in cell_lines) {
    d <- compound_DRC %>%
      dplyr::filter(ccle_name == cell_line)
    cultures <- d$culture %>% unique()
    # for each culture generate a graph
    for(cult in cultures){
      d_cult <- dplyr::filter(d, culture == cult)
      d_cult_line <- dplyr::filter(compound_LFC, culture == cult,
                                   ccle_name == cell_line)
      # DRC curve function
      f1 = function(x) {
        d_cult$lower_limit + (d_cult$upper_limit - d_cult$lower_limit)/
          (1 + (2^x/d_cult$ec50)^d_cult$slope)
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
        ylim(0,2) +
        labs(x = 'log2(Dose) (uM)', y = 'Viability', color = "",
             title = paste0(toupper(id), " - ",
                            d_cult$pert_name, "\n", cell_line,' - ', cult,
                            "\nAUC:", round(d_cult$auc,2),
                            " - IC50:", round(2^d_cult$log2.ic50,2)))
      # outputs to .pdf
      print(p)
    }
  }
  # closes .pdf
  dev.off()
}
