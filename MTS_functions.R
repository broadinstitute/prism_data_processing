library(useful)
library(magrittr)
library(tidyverse)
library(glmnet)
library(scam)
library(dr4pl)
library(readr)
library(sva)
library(sandwich)
library(ranger)


#---- General functions

# area under curve given dose response parameters
compute_auc <- function(l, u, ec50, h, md, MD) {
  #if( l > 1) l = 1
  #if( l < 0) l = 0
  f1 = function(x) pmax(pmin((l + (u - l)/(1 + (2^x/ec50)^h)),1), 0 )
  integrate(f1, log2(md),log2(MD))$value/(log2(MD/md))
}

# log IC50 from given DRC parameters
compute_log_ic50 <- function(l, u, ec50, h, md, MD) {
  if((l >= 0.5) | (u <= 0.5)){
    return(NA)
  }else{
    f1 = function(x) (l + (u - l)/(1 + (2^x/ec50)^h) - 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root, 
                    error = function(x) NA))
  }
}

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
    dplyr::mutate(culture = word(cond, 1, sep = fixed("::")),
                  pool_id = word(cond, 2, sep = fixed("::")),
                  prism_replicate = word(cond, 3, sep = fixed("::"))) %>%
    dplyr::select(-cond, -V2)
  
  Y %>% 
    dplyr::left_join(combat) %>% 
    .$LFC.cb
}


#---- Conitnuous variables biomarkers ----
# function to correlate a vector y with a matrix X using only the top features
correlate <- function(X, y, max_features) {
  
  # shared cell lines (feature matrix and results)
  overlap <- intersect(rownames(X), names(y))
  
  # correlate and filter top n by effect size
  corr_mat <- cor(X[overlap,], y[overlap], use = "pairwise.complete.obs")
  top_n <- scale(X[overlap, rank(-abs(corr_mat)) <= max_features])
  top_n <- top_n[, apply(top_n, 2, var) > 0]
  
  # calculate p values and make tibble
  corr_table <- tibble(feature = colnames(top_n),
                       coeff = apply(top_n, 2, function(x) {
                         m = lm(y ~ x, data = data.frame(y = y[overlap],
                                                         x = x))
                         coef(m)[['x']]
                       }),
                       p.val = apply(top_n, 2, function(x) {
                         m = lm(y ~ x, data = data.frame(y = y[overlap],
                                                         x = x))
                         l = lmtest::coeftest(m, 
                                              vcov = 
                                                sandwich::vcovHC(m, 
                                                                 type = "HC1"))
                         return(l[2, 4])
                       })
  )
  
  # calculate q values
  corr_table %<>%
    dplyr::arrange(desc(abs(coeff))) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::mutate(q.val = p.adjust(c(p.val,rep(1, ncol(X) - ncol(top_n))),
                                   method = "BH")[1:ncol(top_n)])
  
  return(corr_table)
}

# funtion to correlate continous variables with viability metrics
correlate_all <- function(X, Y, comps, max_features, feature_type){
  
  out_table <- tibble()
  
  # for each compound
  for (i in 1:nrow(comps)) {
    comp <- comps[i,]
    
    # select relevant data
    d <- dplyr::filter(Y, pert_mfc_id == comp$pert_mfc_id) %>%
      dplyr::group_by(ccle_name) %>%
      dplyr::summarize(response = mean(response, na.rm = TRUE)) %>%
      dplyr::filter(is.finite(response))
    
    # convert to response matrix
    d_mat <- d$response; names(d_mat) <- d$ccle_name
    
    corr_table <- correlate(X, d_mat, max_features)
    
    # calculate q values
    corr_table %<>%
      dplyr::mutate(pert_mfc_id = comp$pert_mfc_id) %>%
      dplyr::mutate(dose = Y$dose[1]) %>%
      dplyr::mutate(feature_type = feature_type)
    
    out_table %<>%
      dplyr::bind_rows(corr_table)
  }
  
  return(out_table)
}


#---- Discrete variables biomarkers ----
# funtion to run t-test for discrete variables based on viability metrics
discrete_test <- function(X, y) {
  
  out_table <- tibble()  # tracks output
  feats <- colnames(X)  # features to run tests on
  overlap <- dplyr::intersect(rownames(X), names(y))
  y <- y[overlap]
  X <- X[overlap,]
  
  # loop through each feature (e.g. lineage)
  for (feat in feats) {
    has_feat <- rownames(X[X[,feat] == 1,])
    group <- y[has_feat]
    others <- dplyr::setdiff(y, group)
    
    # need more than one data point to generate a group mean and test
    if (length(group) > 1) {
      
      # get t test results (two-sample unpaired)
      t <- stats::t.test(group, others)
      p_val <- t$p.value
      t_stat <- t$statistic["t"]
      effect_size <- mean(group) - mean(others)
      
      result <- tibble(feature = feat,
                       effect_size = effect_size,
                       t = t_stat,
                       p = p_val
      )
      out_table %<>%
        dplyr::bind_rows(result)
    }
  }
  # calculate q values (BH procedure)
  out_table$q <- p.adjust(c(out_table$p,rep(1, length(feats) - nrow(out_table))),
                          method = "BH")[1:nrow(out_table)]
  
  return(out_table)
}
