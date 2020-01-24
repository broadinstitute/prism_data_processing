# Script to run MTS biomarker analyses
# Replace expression_dir and response_dir with appropriate directories to run

# source functions from functions script (also loads libraries)
suppressMessages(source("./src/analysis_functions.R"))

# base directories
expression_dir <- # INSERT path to folder with general expression data
response_dir <- # INSERT path to folder with project specific results

#---- Load data ----
print("Loading the data")
GE <- data.table::fread(paste0(expression_dir, "/ge.csv")) %>%  # gene exp.
  as.matrix(., rownames = "V1")
CNA <- data.table::fread(paste0(expression_dir, "/cna.csv")) %>%  # copy num
  as.matrix(., rownames = "V1")
MET <- data.table::fread(paste0(expression_dir, "/met.csv")) %>%  # metablomics
  as.matrix(., rownames = "V1")
miRNA <- data.table::fread(paste0(expression_dir, "/mirna.csv")) %>%  # micro RNA
  as.matrix(., rownames = "V1")
MUT <- data.table::fread(paste0(expression_dir, "/mut.csv")) %>%  # mutations
  as.matrix(., rownames = "V1")
PROT <- data.table::fread(paste0(expression_dir, "/prot.csv")) %>%  # proteomics
  as.matrix(., rownames = "V1")
XPR <- data.table::fread(paste0(expression_dir, "/xpr.csv")) %>%  # CRISPR KO
  as.matrix(., rownames = "V1")
LIN <- data.table::fread(paste0(expression_dir, "/lin.csv")) %>%  # lineage
  as.matrix(., rownames = "V1")
shRNA <- data.table::fread(paste0(expression_dir, "/shrna.csv")) %>%  # lineage
  as.matrix(., rownames = "V1")
repurposing <- data.table::fread(paste0(expression_dir, "/rep.csv")) %>%  # repurposing
  as.matrix(., rownames = "V1")

continuous_features <- list(GE, CNA, MET, miRNA, PROT, XPR, shRNA, repurposing)
continuous_names <- c("GE", "CNA", "MET", "miRNA", "PROT", "XPR", "shRNA", "REP")

discrete_features <- list(LIN, MUT)
discrete_names <- c("LIN", "MUT")

# combinations of features for multivariate
X.all <- data.table::fread(paste0(expression_dir, "/x-all.csv")) %>%
  subset(select=which(!duplicated(names(.)))) %>%
  unique() %>%
  as.matrix(., rownames = "V1")
X.ccle <- data.table::fread(paste0(expression_dir, "/x-ccle.csv")) %>%
  subset(select=which(!duplicated(names(.)))) %>%
  unique() %>%
  as.matrix(., rownames = "V1")

# DRC table
DRC <- data.table::fread(paste0(response_dir, "/DRC_TABLE.csv"))
max_dose <- max(DRC$max_dose)
DRC %<>%
  dplyr::distinct(ccle_name, culture, pert_name, pert_mfc_id, auc, log2.ic50) %>%
  dplyr::mutate(log2.ic50 = ifelse((is.finite(auc) & is.na(log2.ic50)),
                                   3 * max_dose, log2.ic50),
                log2.auc = log2(auc)) %>%
  tidyr::gather(key = "dose", value = "response", log2.auc, log2.ic50) %>%
  dplyr::filter(is.finite(response))

# LFC table
LFC <- data.table::fread(paste0(response_dir, "/LFC_TABLE.csv")) %>%
  dplyr::distinct(pert_name, ccle_name, culture, pert_idose, pert_mfc_id, LFC.cb) %>%
  dplyr::rename(response = LFC.cb) %>%
  dplyr::rename(dose = pert_idose) %>%
  dplyr::filter(is.finite(response))

# combined
Yall <- dplyr::bind_rows(DRC, LFC)

print("Loaded the data")


#---- RUN ANALYSIS ----

# each run is a compound and a dose (dose can be actual dose or AUC/IC50)
runs <- Yall %>%
  dplyr::distinct(pert_mfc_id, pert_name, dose)

# empty tibbles to store results
continuous <- tibble()
discrete <- tibble()
rf_results <- tibble()
# enet_results <- tibble()
multi_models <- tibble()

print("Starting analysis")

# loop through runs and do analysis
for (i in 1:nrow(runs)) {

  run <- runs[i,]  # current run, print for debugging/sanity/progress bar
  print(paste(run$pert_name, run$dose))

  # select relevant responses and convert to vector
  Y <- Yall %>%
    dplyr::filter(pert_mfc_id == run$pert_mfc_id, dose == run$dose)
  y <- Y$response; names(y) <- Y$ccle_name

  if (all(is.na(y))) {
    next
  }

  #---- MUTLIVARIATE ANALYSIS ----

  if (length(intersect(rownames(X.all), names(y))) >= 20) {
    multi_results_all <- multivariate(X.all, y, k = 10, n = 250)

    # get results and add to aggregation table
    if (nrow(multi_results_all[[1]]) > 0) {
      result_rf_all <- multi_results_all[[1]] %>%
        dplyr::arrange(desc(RF.imp.mean)) %>%
        dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                      pert_name = run$pert_name,
                      pert_idose = run$dose,
                      rank = 1:n(),
                      model = "all")
    }

    if (nrow(multi_results_all[[2]]) > 0) {
      result_enet_all <- multi_results_all[[2]] %>%
        dplyr::arrange(desc(enet.coef.mean)) %>%
        dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                      pert_name = run$pert_name,
                      pert_idose = run$dose,
                      rank = 1:n(),
                      model = "all")
    }

    result_models_all <- multi_results_all[[3]] %>%
      dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                    pert_name = run$pert_name,
                    pert_idose = run$dose,
                    model = "all")

  } else {
    result_rf_all <- tibble()
    result_enet_all <- tibble()
    result_models_all <- tibble()
  }

  # repeat with CCLE data set
  if (length(intersect(rownames(X.ccle), names(y))) >= 20) {
    multi_results_ccle <- multivariate(X.ccle, y, k = 10, n = 250)
    if (nrow(multi_results_ccle[[1]]) > 0) {
      result_rf_ccle <- multi_results_ccle[[1]] %>%
        dplyr::arrange(desc(RF.imp.mean)) %>%
        dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                      pert_name = run$pert_name,
                      pert_idose = run$dose,
                      rank = 1:n(),
                      model = "ccle")
    }
    if (nrow(multi_results_ccle[[2]])) {
      result_enet_ccle <- multi_results_ccle[[2]] %>%
        dplyr::arrange(desc(abs(enet.coef.mean))) %>%
        dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                      pert_name = run$pert_name,
                      pert_idose = run$dose,
                      rank = 1:n(),
                      model = "ccle")
    }
    result_models_ccle <- multi_results_ccle[[3]] %>%
      dplyr::mutate(pert_mfc_id = run$pert_mfc_id,
                    pert_name = run$pert_name,
                    pert_idose = run$dose,
                    model = "ccle")
  } else {
    result_rf_ccle <- tibble()
    result_enet_ccle <- tibble()
    result_models_ccle <- tibble()
  }

  rf_results %<>%
    dplyr::bind_rows(result_rf_all) %>%
    dplyr::bind_rows(result_rf_ccle)
  # enet_results %<>%
  #   dplyr::bind_rows(result_enet_all) %>%
  #   dplyr::bind_rows(result_enet_ccle)
  multi_models %<>%
    dplyr::bind_rows(result_models_all) %>%
    dplyr::bind_rows(result_models_ccle)

  #---- UNIVARIATE CONTINUOUS ANALYSIS ----
  for (j in 1:length(continuous_features)) {
    feature_table <- continuous_features[[j]]
    feature_type <- continuous_names[j]

    # run correlation
    corr_table <- correlate(feature_table, y, 250)

    corr_table %<>%
      dplyr::mutate(pert_mfc_id = run$pert_mfc_id) %>%
      dplyr::mutate(pert_name = run$pert_name) %>%
      dplyr::mutate(dose = run$dose) %>%
      dplyr::mutate(feature_type = feature_type)

    continuous %<>%
      dplyr::bind_rows(corr_table)
  }

  #---- UNIVARIATE DISCRETE ANALYSIS ----
  if(run$dose != "log2.ic50") {
    for(z in 1:length(discrete_features)) {
      feature_table <- discrete_features[[z]]
      feature_type <- discrete_names[z]

      # run t-test
      t_table <- discrete_test(feature_table, y)

      t_table %<>%
        dplyr::mutate(pert_mfc_id = run$pert_mfc_id) %>%
        dplyr::mutate(pert_name = run$pert_name) %>%
        dplyr::mutate(dose = run$dose) %>%
        dplyr::mutate(feature_type = feature_type)

      # only keep top 100 mutations
      if (feature_type == "MUT" & nrow(t_table) > 0) {
        t_table %<>%
          dplyr::arrange(q) %>%
          dplyr::mutate(rank = 1:n()) %>%
          dplyr::filter(rank <= 100) %>%
          dplyr::select(-rank)
      }

      discrete %<>%
        dplyr::bind_rows(t_table)
    }
  }

  #---- GE CORRELATION (-lineage) ----
  ge_tab <- correlate_regressed(GE, y, LIN, 250) %>%
    dplyr::mutate(pert_mfc_id = run$pert_mfc_id) %>%
    dplyr::mutate(pert_name = run$pert_name) %>%
    dplyr::mutate(dose = run$dose) %>%
    dplyr::mutate(feature_type = "GE_noLIN")

  continuous %<>%
    dplyr::bind_rows(ge_tab)

  print(paste0(i, " of ", nrow(runs), " complete"))
}

# write results to project folder
readr::write_csv(continuous, paste0(response_dir, "/continuous_associations.csv"))
readr::write_csv(discrete, paste0(response_dir, "/discrete_associations.csv"))
readr::write_csv(rf_results, paste0(response_dir, "/RF_table.csv"))
# readr::write_csv(enet_results, paste0(response_dir, "/ENet_table.csv"))
readr::write_csv(multi_models, paste0(response_dir, "/Model_table.csv"))
