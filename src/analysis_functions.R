# Helper functions for the analysis step of the MTS pipeline

#---- Load Packages ----
library(plyr)
library(tidyverse)
library(magrittr)
library(tidyr)
library(tibble)
library(glmnet)
library(ranger)
library(glmnet)
library(sandwich)
library(stats)
library(reshape2)
library(dplyr)
library(lmtest)


#---- Conitnuous variables biomarkers ----
# function to correlate a vector y with a matrix X using only the top features
correlate <- function(X, y, max_features) {
  
  # shared cell lines (feature matrix and results)
  overlap <- intersect(rownames(X), names(y))
  
  if (length(overlap) <= 10 & !all(is.na(y[overlap]))) {
    return(tibble())
  }
  
  # correlate and filter top n by effect size
  corr_mat <- cor(X[overlap,], y[overlap], use = "pairwise.complete.obs")
  top_n <- scale(X[overlap, rank(-abs(corr_mat)) <= max_features])
  top_n <- top_n[, apply(top_n, 2, function(x) var(x, na.rm = TRUE)) > 0]
  y <- scale(y)[,1]
  
  corr_table <- tryCatch(expr = {
    tibble(feature = colnames(top_n),
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
  }, error = function(e) {
    print(paste("Unable to correlate: ", e))
    return(tibble())
  })
  
  if (length(corr_table) < 1) {
    return(corr_table)
  }
  
  # calculate q values
  corr_table %<>%
    dplyr::arrange(desc(abs(coeff))) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::mutate(q.val = p.adjust(c(p.val,rep(1, ncol(X) - ncol(top_n))),
                                   method = "BH")[1:ncol(top_n)])
  
  return(corr_table)
}


#---- Discrete variables biomarkers ----
# funtion to run t-test for discrete variables based on viability metrics
discrete_test <- function(X, y) {
  
  out_table <- tibble()  # tracks output
  feats <- colnames(X)  # features to run tests on
  overlap <- dplyr::intersect(rownames(X), names(y))
  
  if (length(overlap) < 10) {
    return(out_table)
  }
  
  y <- y[overlap]
  X <- X[overlap,]
  
  # loop through each feature (e.g. lineage)
  for (feat in feats) {
    has_feat <- rownames(X[X[,feat] == 1,])
    group <- y[has_feat]
    others <- dplyr::setdiff(y, group)
    
    # need more than one data point to generate a group mean and test
    if (length(group) >= 5 & length(others) >= 5) {
      
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
  
  if (nrow(out_table) < 1) {
    return(out_table)
  }
  
  # calculate q values (BH procedure)
  out_table$q <- p.adjust(c(out_table$p,
                            rep(1, length(feats) - nrow(out_table))),
                          method = "BH")[1:nrow(out_table)]
  
  return(out_table)
}


# function to do multivariate analysis on MTS data
# takes named vector of responses (y) and feature matrix (X)
# k = number of cross validation cycles
# n = number of features to test for each run (e.g. top 250)
multivariate <- function(X, y, k = 10, n = 250){
  y <- y[is.finite(y)]  # only finite values
  cl <- sample(intersect(rownames(X), names(y)))  # overlapping rows
  X <- scale(X[cl,]); y = y[cl]; N = floor(length(cl)/k)
  
  colnames(X) %<>% make.names()
  yhat_rf <- y; yhat_enet <- y
  
  SS = tibble()  # output
  
  # run cross validation k times
  for (kx in 1:k) {
    test <- cl[((kx - 1) * N + 1):(kx * N)]  # select test set
    train <- setdiff(cl, test)  # everything else is training
    X_train <- X[train,]
    
    # assumes no NAs in X
    X_train <- X_train[,apply(X_train,2,var) > 0]
    X_test <- X[test,]
    
    # select top n correlated features in X
    X_train <- X_train[,rank(-abs(cor(X_train,y[train])))<= n]
    
    # makes RF and ENet models
    rf <- ranger::ranger(y ~ .,
                         data = as.data.frame(cbind(X_train, y = y[train])),
                         importance = "impurity")
    
    # # enet errors out under certain conditions
    # enet <- tryCatch(expr = {
    #   glmnet::cv.glmnet(x = X_train,
    #                     y = y[train],
    #                     alpha = .5)
    # }, error = function(e) {
    #   print(paste0("Warning, there was an error: ", e))
    #   return(NULL)
    # })
    
    # # generate predictions and get variable importances
    # if(!is.null(enet)) {
    #   yhat_enet[test] <-
    #     predict(enet, newx = X_test[, colnames(X_train)])
    #   a <- coef(enet)
    # } else {
    #   yhat_enet[test] <- NA
    #   a <- NULL
    # }
    
    enet <- NULL
    
    yhat_rf[test] <-
      predict(rf, data = as.data.frame(X_test[, colnames(X_train)]))$predictions
    
    # get variable importance metrics for rf
    ss <- tibble(feature = names(rf$variable.importance),
                 RF.imp = rf$variable.importance / sum(rf$variable.importance),
                 RF.mse = mean((yhat_rf[test] - y[test])^2)
    )
    
    if(!is.null(enet)) {
      ss %<>%
        dplyr::full_join(tibble(feature = a@Dimnames[[1]][a@i+1],
                                enet.coef = a@x), by = "feature") %>%
        dplyr::mutate(enet.coef = ifelse(is.na(enet.coef), 0, enet.coef),
                      enet.mse = mean((yhat_enet[test] - y[test])^2, na.rm = T))
    }
    
    ss %<>%
      dplyr::mutate(fold = kx)
    
    # add this model to the output
    SS %<>%
      dplyr::bind_rows(ss)
  }
  
  # separate into RF and ENet tables, and one large table of model stats
  model_table <- tibble()
  
  RF.importances <- SS %>%
    dplyr::distinct(feature, RF.imp, fold) %>%
    reshape2::acast(feature ~ fold, value.var = "RF.imp")
  
  # average values across validation runs
  RF.table <- tibble(feature = rownames(RF.importances),
                     RF.imp.mean = RF.importances %>%
                       apply(1, function(x) mean(x, na.rm = T)),
                     RF.imp.sd = RF.importances %>%
                       apply(1, function(x) sd(x, na.rm = T)),
                     RF.imp.stability = RF.importances %>%
                       apply(1, function(x) mean(!is.na(x)))) %>%
    # only keep features used in > half the models
    dplyr::filter((RF.imp.stability > 0.5), feature != "(Intercept)")
  
  # table ofrf model level statistics
  rf_mod <- RF.table %>%
    dplyr::mutate(MSE = mean(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                             na.rm =T),
                  MSE.se = sd(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                              na.rm = T),
                  R2 = 1 - MSE/var(y[1:(k*N)], na.rm = T),
                  PearsonScore = cor(y[1:(k*N)], yhat_rf[1:(k*N)],
                                     use = "pairwise.complete.obs")) %>%
    dplyr::distinct(MSE, MSE.se, R2, PearsonScore) %>%
    dplyr::mutate(type = "RF")
  
  model_table %<>% dplyr::bind_rows(rf_mod)
  
  if ("enet.coef" %in% colnames(SS)) {
    enet.importances <- SS %>%
      dplyr::distinct(feature, enet.coef, fold) %>%
      reshape2::acast(feature ~ fold, value.var = "enet.coef")
    
    enet.table <- tibble(feature = rownames(enet.importances),
                         enet.coef.mean = enet.importances %>%
                           apply(1, function(x) mean(x[x!= 0], na.rm = T)),
                         enet.coef.sd =   enet.importances %>%
                           apply(1, function(x) sd(x[x!= 0], na.rm = T)),
                         enet.coef.stability = enet.importances %>%
                           apply(1, function(x) mean(x != 0, na.rm = T))) %>%
      dplyr::filter((enet.coef.stability >= 0.7), feature != "(Intercept)")
    
    enet_mod <- enet.table %>%
      dplyr::mutate(MSE = mean(dplyr::distinct(SS, enet.mse, fold)$enet.mse,
                               na.rm = T),
                    MSE.se = sd(dplyr::distinct(SS, enet.mse, fold)$enet.mse,
                                na.rm = T),
                    R2 = 1 - MSE/var(y[1:(k*N)], na.rm = T),
                    PearsonScore = cor(y[1:(k*N)], yhat_enet[1:(k*N)],
                                       use = "pairwise.complete.obs")) %>%
      dplyr::distinct(MSE, MSE.se, R2, PearsonScore) %>%
      dplyr::mutate(type = "ENet")
  } else {
    enet.importances <- tibble()
    enet.table <- tibble()
    enet_mod <- tibble()
  }
  
  
  model_table %<>% dplyr::bind_rows(enet_mod)
  
  
  return(list(RF.table, enet.table, model_table, yhat_rf, yhat_enet))
}

# function to correlate a vector y with a matrix X using only the top features
correlate_regressed <- function(X, y, regressors, max_features) {
  y <- y[is.finite(y)]  # only finite values
  
  # shared cell lines (feature matrix and results)
  overlap <- intersect(rownames(X), names(y))
  overlap <- intersect(overlap, rownames(regressors))
  
  if (length(overlap) <= 10) {
    return(tibble())
  }
  
  # correlate and filter top n by effect size
  corr_mat <- cor(X[overlap,], y[overlap], use = "pairwise.complete.obs")
  top_n <- scale(X[overlap, rank(-abs(corr_mat)) <= 2 * max_features])
  top_n <- top_n[, apply(top_n, 2, function(x) var(x, na.rm = TRUE)) > 0]
  y <- scale(y)[,1]
  R <- regressors[overlap,]
  
  corr_table <- tryCatch(expr = {
    tibble(feature = colnames(top_n),
           coeff = apply(top_n, 2, function(x) {
             m = lm(y[overlap] ~ x + R)
             coef(m)[['x']]
           }),
           p.val = apply(top_n, 2, function(x) {
             m = lm(y[overlap] ~ x + R)
             l = lmtest::coeftest(m,
                                  vcov =
                                    sandwich::vcovHC(m,
                                                     type = "HC1"))
             return(l[2, 4])
           }))
    
  }, error = function(e) {
    print(paste("Unable to correlate: ", e))
    return(tibble())
  })
  
  if (length(corr_table) < 1) {
    return(corr_table)
  }
  
  # calculate q values
  corr_table %<>%
    dplyr::arrange(desc(abs(coeff))) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::mutate(q.val = p.adjust(c(p.val,rep(1, ncol(X) - ncol(top_n))),
                                   method = "BH")[1:ncol(top_n)]) %>%
    dplyr::filter(rank <= max_features)
  
  return(corr_table)
}

