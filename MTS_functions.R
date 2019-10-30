library(plyr)
library(tidyverse)
library(magrittr)
library(data.table)
library(tidyr)
library(glmnet)
library(ranger)
library(scam)
library(dr4pl)
library(glmnet)
library(sva)
library(sandwich)
library(lmtest)
library(stats)
library(PRROC)
library(reshape2)
library(dplyr)


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
    dplyr::mutate(culture = stringr::word(cond, 1, sep = stringr::fixed("::")),
                  pool_id = stringr::word(cond, 2, sep = stringr::fixed("::")),
                  prism_replicate = stringr::word(cond, 3, sep = stringr::fixed("::"))) %>%
    dplyr::select(-cond, -V2)

  Y %>%
    dplyr::left_join(combat) %>%
    .$LFC.cb
}
