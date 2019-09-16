# Script to run initial processing on MTS data (logMFI.csv files)

# import necessary libraries and functions using MTS_functions.R
source("./MTS_functions.R")


#---- Load Data ----
data_path <- # INSERT PATH TO DATA HERE: should be "~/.../logMFI.csv"
unprocessed <- data.table::fread(data_path)

# check for duplicates and append "_2" to pert_mfc_id
dups <- unprocessed %>%
  dplyr::filter(pert_type == "trt_cp") %>%
  dplyr::group_by(ccle_name, pert_mfc_id, pert_dose, culture,
                  prism_replicate, pool_id) %>%
  dplyr::mutate(count = dplyr::n_distinct(pert_well)) %>%
  dplyr::filter(count > 1) %>%
  dplyr::ungroup()

if (nrow(dups) > 0) {
  dups %<>%
    dplyr::distinct(ccle_name, rid, pert_mfc_id,pert_dose, culture,
                    prism_replicate, pool_id) %>%
    plyr::join(., unprocessed, match = "first")
  
  unprocessed %<>%
    dplyr::anti_join(dups)
  
  dups %<>%
    dplyr::mutate(pert_mfc_id = paste(pert_mfc_id, 2, sep = "_"))
  
  unprocessed %<>%
    dplyr::bind_rows(dups)
}

# split into 300 and 500
PR300 <- unprocessed %>%
  filter(culture == "PR300")

PR500 <- unprocessed %>%
  filter(culture == "PR500")

# create barcode tables
PR300_barcodes <- PR300 %>%
  filter(pool_id == "CTLBC")

PR500_barcodes <- PR500 %>%
  filter(pool_id == "CTLBC")


#---- Basic QC and Normalization ----

# compute median of medians for normalization
PR300_normalized <- PR300 %>%
  dplyr::filter(pert_type == "ctl_vehicle") %>%
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::mutate(mLMFI = median(logMFI)) %>%
  dplyr::group_by(prism_replicate, rid) %>%
  dplyr::mutate(mmLMFI = logMFI - mLMFI + median(mLMFI)) %>%
  dplyr::summarize(rLMFI = median(mmLMFI)) %>%
  dplyr::left_join(PR300)

# fit curve to controls and predict test conditions
PR300_normalized <- PR300_normalized %>%
  dplyr::group_by(prism_replicate, pert_well)%>%
  dplyr::mutate(LMFI = scam(y ~ s(x, bs = "mpi"),
                          data = tibble(
                            y = rLMFI[rid %in% PR300_barcodes$rid],
                            x = logMFI[rid %in% PR300_barcodes$rid])) %>% 
                  predict(newdata = tibble(x = logMFI))) %>%
  dplyr::ungroup() %>%
  dplyr::select(-logMFI)

# join with other info
PR300_normalized <- PR300_normalized %>%
  left_join(PR300) %>%
  select(-logMFI)

# REPEAT with 500
PR500_normalized <- PR500 %>%
  dplyr::filter(pert_type == "ctl_vehicle") %>%
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::mutate(mLMFI = median(logMFI)) %>%
  dplyr::group_by(prism_replicate, rid) %>%
  dplyr::mutate(mmLMFI = logMFI - mLMFI + median(mLMFI)) %>%
  dplyr::summarize(rLMFI = median(mmLMFI)) %>%
  dplyr::left_join(PR500)

PR500_normalized <- PR500_normalized %>%
  dplyr::group_by(prism_replicate, pert_well)%>%
  dplyr::mutate(LMFI = scam(y ~ s(x, bs = "mpi"),
                          data = tibble(
                            y = rLMFI[rid %in% PR500_barcodes$rid],
                            x = logMFI[rid %in% PR500_barcodes$rid])) %>% 
                  predict(newdata = tibble(x = logMFI))) %>%
  dplyr::ungroup() %>%
  dplyr::select(-logMFI)

PR500_normalized <- PR500_normalized %>%
  left_join(PR500) %>%
  select(-logMFI)

#---- Generate SSMD table ----

SSMD_table_300 <- PR300_normalized %>%
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

# calculate error rate of normalized table
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
SSMD_table_300 <- SSMD_table_300 %>%
  dplyr::left_join(PR300_error)

# REPEAT with 500
SSMD_table_500 <- PR500_normalized %>%
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon")) %>%
  dplyr::distinct(ccle_name, pert_type, prism_replicate,LMFI, profile_id,
                  pool_id, culture) %>%
  dplyr::group_by(pert_type, prism_replicate, ccle_name, pool_id, culture) %>%
  dplyr::summarize(med = median(LMFI, na.rm = TRUE), 
                   mad = mad(LMFI, na.rm = TRUE)) %>%
  dplyr::mutate(pert_type_md = paste0(pert_type, '_md'),
                pert_type_mad = paste0(pert_type, '_mad')) %>%
  tidyr::spread(key = pert_type_md, value = med, fill = 0) %>% 
  tidyr::spread(key = pert_type_mad, value = mad, fill = 0) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-pert_type) %>%
  dplyr::group_by(prism_replicate, ccle_name, pool_id, culture) %>%
  dplyr::summarise_all(sum) %>%
  dplyr::mutate(ssmd = (ctl_vehicle_md - trt_poscon_md) / 
                  sqrt(ctl_vehicle_mad^2 +trt_poscon_mad^2),
                nnmd = (ctl_vehicle_md - trt_poscon_md) / ctl_vehicle_mad)

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

SSMD_table_500 <- SSMD_table_500 %>%
  dplyr::left_join(PR500_error)

# combine 300 and 500 tables
SSMD_TABLE <- dplyr::bind_rows(SSMD_table_500, SSMD_table_300) %>%
  # if error rate <= .05 then pass
  dplyr::mutate(pass = error_rate <= 0.05,
                compound_plate =  word(prism_replicate, 1,2, 
                                       sep = fixed("_"))) %>%
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
  dplyr::mutate(compound_plate = word(prism_replicate, 1,2,
                                      sep = fixed("_"))) %>% 
  tidyr::unite(col = "condition", pert_mfc_id, pert_dose, compound_plate,
               sep = "::", remove = FALSE) %>%
  split(.$condition) %>%
  purrr::map_dfr(~dplyr::mutate(.x, LFC.cb = apply_combat(.))) %>%
  dplyr::select(-condition)

#---- Compute dose-response parameters ----

# table of all cultures with >4 doses (necessary number to fit curve)
# with uncorrected values
DRC_TABLE <- LFC_TABLE %>%
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(ccle_name, culture, pert_mfc_id, pert_dose) %>%
  dplyr::count(ccle_name, culture, pert_mfc_id) %>%
  dplyr::filter(n > 4) %>%
  dplyr::mutate(ix = 1:n())

# empty tibble to be filled with DRC parameters
DRC <- tibble()

for(jx in 1:nrow(DRC_TABLE)){
  d = DRC_TABLE %>%
    # select one row (cell line x compund) and join with corresponding LFC data
    dplyr::filter(ix == jx) %>%
    dplyr::left_join(LFC_TABLE)
  
  # fit decreasing logistic function to data
  a = tryCatch(dr4pl(dose = d$pert_dose,
                     response = 2^d$LFC,
                     method.init = "logistic",
                     trend = "decreasing"), 
               error = function(e) NA)
  # store parameters
  param <- tryCatch(summary(a)$coefficients$Estimate, error = function(e) NA)
  # if successfully generated parameters make dataframe of them
  if(!is.na(param)){
    # extract mean squared error
    error_val <- tryCatch(as.double(gof.dr4pl(a)[[2,3]]),
                          error = function(e) NA)
    # calcuclate our own MSE (capping predictions and values at 1)
    mse_capped <- 0
    for(z in 1:nrow(d)){
      e <- (max(0, min(1, 2^d$LFC[[z]]))
            - max(0, min(1, dr4pl::MeanResponse(d$pert_dose[[z]], param))))^2
      mse_capped <- mse_capped + e
    }
    mse_capped <- mse_capped/nrow(d)
    
    x <- tibble(ix = jx,
               upper_limit = param[1],
               ec50 = param[2],
               slope = -param[3],
               lower_limit = param[4],
               convergence = a$convergence) %>%
      # compute area under DRC curve
      dplyr::mutate(auc = compute_auc(lower_limit,
                                      upper_limit,
                                      ec50,
                                      slope, 
                                      min(d$pert_dose),
                                      max(d$pert_dose)),
                    # compute IC50
                    log2.ic50 = compute_log_ic50(lower_limit, 
                                                 upper_limit, 
                                                 ec50, slope, 
                                                 min(d$pert_dose), 
                                                 max(d$pert_dose)),
                    mse = error_val,
                    mse_cap = mse_capped)
    # add to dataframe
    DRC %<>% 
      dplyr::bind_rows(x)
  }
}

# join with full table (by ix)
DRC_TABLE <- DRC %>%
  dplyr::filter(convergence) %>% 
  dplyr::left_join(DRC_TABLE) %>% 
  dplyr::select(-ix, -convergence, -n)


# REPEAT with ComBat corrected values
DRC_TABLE_cb <- LFC_TABLE %>%
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(ccle_name, culture, pert_mfc_id, pert_dose) %>%
  dplyr::count(ccle_name, culture, pert_mfc_id) %>%
  dplyr::filter(n > 4) %>%
  dplyr::mutate(ix = 1:n())

DRC_cb <- tibble()

for(jx in 1:nrow(DRC_TABLE_cb)){
  d = DRC_TABLE_cb %>%
    dplyr::filter(ix == jx) %>%
    dplyr::left_join(LFC_TABLE)
  
  a = tryCatch(dr4pl(dose = d$pert_dose,
                     response = 2^d$LFC.cb,
                     method.init = "logistic",
                     trend = "decreasing"), 
               error = function(e) NA)
  param <- tryCatch(summary(a)$coefficients$Estimate, error = function(e) NA)
  if(!is.na(param)){
    error_val <- tryCatch(as.double(gof.dr4pl(a)[[2,3]]),
                          error = function(e) NA)
    mse_capped <- 0
    for(z in 1:nrow(d)){
      e <- (max(0, min(1, 2^d$LFC[[z]]))
            - max(0, min(1, dr4pl::MeanResponse(d$pert_dose[[z]], param))))^2
      mse_capped <- mse_capped + e
    }
    mse_capped <- mse_capped/nrow(d)
    x <- tibble(ix = jx, 
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
                    mse = error_val,
                    mse_cap = mse_capped)
    DRC_cb %<>% dplyr::bind_rows(x)
  }
}

DRC_TABLE_cb <- DRC_cb %>%
  dplyr::filter(convergence) %>% 
  dplyr::left_join(DRC_TABLE_cb) %>% 
  dplyr::select(-ix, -convergence, -n)


# add back pert_name to DRC tables
pert_names <- dplyr::select(LFC_TABLE, pert_mfc_id, pert_name) %>%
  unique()
DRC_TABLE %<>%
  dplyr::left_join(pert_names)
DRC_TABLE_cb %<>%
  dplyr::left_join(pert_names)


#---- Make collapsed LFC table ----

LFC_COLLAPSED_TABLE <- LFC_TABLE %>% 
  dplyr::mutate(compound_plate = word(prism_replicate, 1,2, 
                                      sep = fixed("_"))) %>% 
  dplyr::group_by(ccle_name, culture, pert_name, pert_mfc_id,
                  pert_dose, pert_idose, compound_plate) %>%
  # LFC and LFC.cb values will be medains across replicates
  dplyr::summarize(LFC = median(LFC, na.rm = TRUE),
                   LFC.cb = median(LFC.cb, na.rm = TRUE))


#---- Write to .csv ----

readr::write_csv(LFC_TABLE, 
                 paste0(dirname(data_path), "/LFC_TABLE.csv"))
readr::write_csv(LFC_COLLAPSED_TABLE, 
                 paste0(dirname(data_path), "/LFC_COLLAPSED_TABLE.csv"))
readr::write_csv(SSMD_TABLE, paste0(dirname(data_path), "/SSMD_TABLE.csv"))
readr::write_csv(DRC_TABLE, paste0(dirname(data_path), "/DRC_TABLE.csv"))
readr::write_csv(DRC_TABLE_cb, paste0(dirname(data_path), "/DRC_TABLE_cb.csv"))


#---- Generate DRC plots ----

# generate a .pdf of graphs for each compound
for(compound in pert_names$pert_mfc_id){
  
  # filter to just see that compound
  compound_DRC <- DRC_TABLE_cb %>%
    dplyr::filter(pert_mfc_id == compound) %>%
    dplyr::arrange(desc(auc))
  
  # tracks LFC info
  compound_LFC <- LFC_TABLE %>%
    dplyr::filter(pert_mfc_id == compound)
  
  # create .pdf
  pdf(paste0(dirname(data_path), "/", toupper(compound), "_DRCfigures.pdf"))
  
  # loop through each cell line treated by compound and plot DRC
  cell_lines <- compound_DRC$ccle_name %>% unique()
  for(cell_line in cell_lines){
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
                       color = prism_replicate, y = 2^LFC)) + 
        geom_line(data = tibble(x = xx, y = f1(xx)),
                  aes(x = x, y = y, group = 1),  lwd =1 ) +
        ylim(0,2) + 
        labs(x = 'log2(dose)', y = 'viability', color = "",
             title = paste0(toupper(compound), " - ", 
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

