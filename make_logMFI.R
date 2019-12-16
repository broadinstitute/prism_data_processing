# script to take raw data and make logMFI csv

library(dplyr)
library(data.table)
library(tidyverse)
library(readr)
library(dr4pl)
library(magrittr)
library(hdf5r)


# script takes the name of the directory where data is stored as arg
script_args <- commandArgs(trailingOnly = TRUE)
if (length(script_args) == 0) {
  stop("Please supply path to data directory", call. = FALSE)
}

# directory with data (also where new folders will be created)
base_dir <- script_args[1]

# paths to data (make sure directory of data has these files)
path_500 <- paste0(base_dir,"/PR500_LMFI.gctx")
path_300 <- paste0(base_dir,"/PR300_LMFI.gctx")
path_cell_info_500 <- paste0(base_dir,"/PR500_cell_info.txt")
path_cell_info_300 <- paste0(base_dir,"/PR300_cell_info.txt")
path_inst_info_500 <- paste0(base_dir,"/PR500_inst_info.txt")
path_inst_info_300 <- paste0(base_dir,"/PR300_inst_info.txt")
path_skipped <- paste0(base_dir, "/skipped_wells.csv")  # if this data exists


# HDF5 file reader
read_hdf5 = function (filename, index = NULL) {
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


#---- Load the data ----

# data table linking drugs to projects (collaborators)
key_table <- data.table::fread(path_key)

# data table of skipped wells (make sure pool ids are 3 digits)
skipped <- data.table::fread(path_skipped) %>%
  dplyr::rename(pert_mfc_id = broad_sample,
                pert_well = `Wells Skipped`,
                pert_dose = mmoles_per_liter,
                rep = Replicate,
                pool_id = Pool,
                plate = pert_plate_src) %>%
  dplyr::distinct(plate, pert_mfc_id, pert_well, rep, pool_id)

# read in logMFI data
PR500 <- read_hdf5(path_500)
PR300 <- read_hdf5(path_300)
rownames(PR500) = paste0(rownames(PR500), "_", "PR500")
rownames(PR300) = paste0(rownames(PR300), "_", "PR300")

# read in cell info
cell_info_500 <- data.table::fread(path_cell_info_500) %>%
  dplyr::distinct(rid, ccle_name, pool_id) %>%
  dplyr::mutate(culture = "PR500") %>%
  dplyr::mutate(rid = paste0(rid, "_", culture)) %>%
  dplyr::mutate(pool_id = ifelse(pool_id == "", "CTLBC", pool_id))

cell_info_300 <- data.table::fread(path_cell_info_300) %>%
  dplyr::distinct(rid, ccle_name, pool_id) %>%
  dplyr::mutate(culture = "PR300") %>%
  dplyr::mutate(rid = paste0(rid, "_", culture)) %>%
  dplyr::mutate(pool_id = ifelse(pool_id == "", "CTLBC", pool_id))

# read in inst info
inst_info_500 <- data.table::fread(path_inst_info_500) %>%
  dplyr::filter(!is_well_failure) %>%
  dplyr::distinct(profile_id, pert_dose, pert_idose, pert_iname, pert_mfc_id,
                  pert_type, pert_well, prism_replicate) %>%
  dplyr::mutate(compound_plate =  word(prism_replicate, 1,2, sep = fixed("_")))

inst_info_300 <- data.table::fread(path_inst_info_300) %>%
  dplyr::filter(!is_well_failure) %>%
  dplyr::distinct(profile_id, pert_dose, pert_idose, pert_iname, pert_mfc_id,
                  pert_type, pert_well, prism_replicate) %>%
  dplyr::mutate(compound_plate =  word(prism_replicate, 1,2, sep = fixed("_")))

# ensure unique profile IDs
PR500 <- PR500[, inst_info_500$profile_id %>% unique()]
PR300 <- PR300[, inst_info_300$profile_id %>% unique()]

# melt matrices into data tables and join with inst and cell info
PR500_molten <- log2(PR500) %>%
  reshape2::melt(varnames = c("rid", "profile_id"), value.name = "logMFI") %>%
  dplyr::left_join(cell_info_500) %>%
  dplyr::left_join(inst_info_500)

PR300_molten <- log2(PR300) %>%
  reshape2::melt(varnames = c("rid", "profile_id"), value.name = "logMFI") %>%
  dplyr::left_join(cell_info_300) %>%
  dplyr::left_join(inst_info_300)

# bind tables together (reorder columns)
master_logMFI <- PR500_molten %>%
  dplyr::bind_rows(PR300_molten) %>%
  dplyr::mutate(pert_name = pert_iname) %>%
  dplyr::select(profile_id, rid, ccle_name, pool_id, culture, prism_replicate,
                pert_type, pert_dose, pert_idose, pert_mfc_id, pert_name, pert_well,
                logMFI) %>%
  dplyr::filter(prism_replicate != "PMTS.BASE001_PR500.C25_0H_X1_P5")

# change validation (.es) to treatment for processing
master_logMFI$pert_type[which(master_logMFI$pert_type == "trt_poscon.es")] <-
  "trt_cp"

# remove pools that were skipped
master_logMFI %<>%
  tidyr::separate(profile_id,
                  c("plate", "ignore", "ignore2", "rep", "ignore3"), "_",
                  remove = FALSE) %>%
  dplyr::select(-ignore, -ignore2, -ignore3)

master_logMFI$rep <- substr(master_logMFI$rep, 1, 2)

master_logMFI %<>%
  dplyr::anti_join(skipped) %>%
  dplyr::select(-plate, -rep)

readr::write_csv(master_logMFI, paste0(base_dir, "/logMFI.csv"))
