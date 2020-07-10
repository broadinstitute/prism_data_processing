# Helper functions for the analysis step of the MTS pipeline

#---- Load Packages ----
library(plyr)
library(tidyverse)
library(magrittr)
library(ranger)
library(sandwich)
library(reshape2)
library(lmtest)
library(limma)
library(corpcor)
library(ashr)
library(gmodels)
library(WGCNA)

#---- Continuous ----
# correlation function
lin_associations <- function (X, Y, W=NULL, n.min=11, shrinkage=T, alpha=0) {
  # load necessary packages
  
  # NA handling
  X.NA <- !is.finite(X); Y.NA <- !is.finite(Y)
  N <- (t(!X.NA) %*% (!Y.NA)) - 2
  
  # if there are provided confounders, calculate P matrix and use to control
  if (!is.null(W)) {
    P <- W %*% corpcor::pseudoinverse(t(W) %*% W) %*% t(W)
    X[X.NA] <- 0; X <- X - (P %*% X); X[X.NA] <- NA
    Y[Y.NA] <- 0; Y <- Y - (P %*% Y); Y[Y.NA] <- NA
    N <- N - dim(W)[2]
  }
  
  # scale matrices
  X <- scale(X); sx <- attributes(X)$`scaled:scale`
  Y <- scale(Y); sy <- attributes(Y)$`scaled:scale`
  
  # calculate correlations
  rho <- WGCNA::cor(X, Y, use = "pairwise.complete.obs")
  
  # correct correlation estimates
  beta <- t(t(rho / sx) * sy)
  beta.se <- t(t(sqrt(1 - rho ^ 2) / sx) * sy) / sqrt(N)
  
  # generate p and q values (without adaptive shrinkage)
  p.val <- 2 * pt(-abs(beta / beta.se), N)
  p.val[N < n.min] <- NA
  p.val[(sx == 0) | !is.finite(sx), ] <- NA
  p.val[, (sy == 0) | !is.finite(sy)]  <- NA
  q.val = apply(p.val, 2, function(x) p.adjust(x, method = "BH"))
  
  # if shrinkage wanted, generate res.table with adaptive shrinkage results
  if(shrinkage){
    res.table <- list()
    for (ix in 1:dim(Y)[2]) {
      fin <- is.finite(p.val[, ix])
      res <- tryCatch(ashr::ash(beta[fin, ix], beta.se[fin, ix],
                                mixcompdist = "halfuniform", alpha = alpha)$result,
                      error = function(e) NULL
      )
      if (!is.null(res)) {
        res$dep.var <- colnames(Y)[ix]
        res$ind.var <- rownames(res)
        res <- cbind(res, rho=rho[rownames(res),], q.val=q.val[rownames(res),])
        res.table[[ix]] <- res
      }
    }
    res.table <- dplyr::bind_rows(res.table)
  } else{
    res.table <- NULL
  }
  
  return(list(
    N = N,
    rho = rho,
    beta = beta,
    beta.se = beta.se,
    p.val = p.val,
    q.val = q.val,
    res.table = res.table
  ))
}

#---- Discrete variables biomarkers ----
# function to run limma on a matrix  (helper function for discrete test)
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL,
                               target_type = 'Gene', limma_trend = FALSE) {
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = p.adjust(p.left, method = 'BH'),
                             q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}

# funtion to run t-test for discrete variables based on viability metrics
discrete_test <- function(X, y) {
  
  y <- y[is.finite(y)]  # only finite values
  X <- X[, apply(X, 2, function(x) !any(is.na(x)))]
  
  overlap <- dplyr::intersect(rownames(X), names(y))
  
  if (length(overlap) < 10) {
    stop("Not enough overlapping data")
  }
  
  y <- y[overlap]
  X <- X[overlap,]
  
  out_table <- run_lm_stats_limma(X, y, target_type = "feature")
  out_table %<>%
    dplyr::select(feature, EffectSize, t_stat, p.value, q.value) %>%
    dplyr::rename(effect_size = EffectSize)
  
  return(out_table)
}


#---- Multivariate ----
# function to fit a random forest analysis to MTS data
# takes named vector of responses (y) and feature matrix (X)
# k = number of cross validation cycles
# n = number of features to test for each run (e.g. top 250)
random_forest <- function(X, y, k = 10, n = 500){
  y <- y[is.finite(y)]  # only finite values
  cl <- sample(intersect(rownames(X), names(y)))  # overlapping rows
  y = y[cl]; N = floor(length(cl)/k)
  
  colnames(X) %<>% make.names()
  yhat_rf <- y; yhat_enet <- y
  
  SS = tibble()  # output
  
  # run cross validation k times
  for (kx in 1:k) {
    test <- cl[((kx - 1) * N + 1):(kx * N)]  # select test set
    train <- setdiff(cl, test)  # everything else is training
    X_train <- X[train,]
    
    # assumes no NAs in X
    X_train <- X_train[,apply(X_train, 2, var) > 0]
    X_test <- X[test,]
    
    # select top n correlated features in X
    X_train <- X_train[,rank(-abs(cor(X_train, y[train])))<= n]
    
    # makes RF and ENet models
    rf <- ranger::ranger(y ~ .,
                         data = as.data.frame(cbind(X_train, y = y[train])),
                         importance = "impurity")
    
    yhat_rf[test] <-
      predict(rf, data = as.data.frame(X_test[, colnames(X_train)]))$predictions
    
    # get variable importance metrics for rf
    ss <- tibble(feature = names(rf$variable.importance),
                 RF.imp = rf$variable.importance / sum(rf$variable.importance),
                 fold = kx)
    
    # add this model to the output
    SS %<>%
      dplyr::bind_rows(ss)
  }
  
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
  
  # table of model level statistics
  mse <- mean((yhat_rf - y)^2)
  mse.se <- sqrt(var((yhat_rf - y)^2))/sqrt(length(y))
  r2 <- 1 - (mse/var(y))
  ps <- cor(y, yhat_rf, use = "pairwise.complete.obs")
  model_table = tibble(MSE = mse,
                       MSE.se = mse.se,
                       R2 = r2,
                       PearsonScore = ps)
  
  return(list(rf.fit = RF.table, mod.tab = model_table, preds = yhat_rf))
}

# modified random forest function to incorporate sample weights
# weights passed as a named vector
random_forest_weighted <- function(X, y, weights = NULL, k = 10, n = 500){
  y <- y[is.finite(y)]  # only finite values
  cl <- sample(intersect(rownames(X), names(y)))  # overlapping rows
  X <- scale(X[cl,]); y = y[cl]; N = floor(length(cl)/k)
  
  # make weight vector
  if(is.null(weights)) {
    weights <- rep(1, length(cl))
    names(weights) <- cl
  }
  W <- weights[cl]
  
  colnames(X) %<>% make.names()
  yhat_rf <- y;
  
  SS = tibble()  # output
  
  # run cross validation k times
  for (kx in 1:k) {
    test <- cl[((kx - 1) * N + 1):(kx * N)]  # select test set
    train <- setdiff(cl, test)  # everything else is training
    X_train <- X[train,]
    W_train <- W[train]
    
    # assumes no NAs in X
    X_train <- X_train[,apply(X_train,2,var) > 0]
    X_test <- X[test,]
    
    # select top n correlated features in X
    X_train <- X_train[,rank(-abs(cor(X_train,y[train])))<= n]
    
    # makes RF and ENet models
    rf <- ranger::ranger(y ~ .,
                         data = as.data.frame(cbind(X_train, y = y[train])),
                         case.weights = W_train,
                         importance = "impurity")
    
    yhat_rf[test] <-
      predict(rf, data = as.data.frame(X_test[, colnames(X_train)]))$predictions
    
    # get variable importance metrics for rf
    ss <- tibble(feature = names(rf$variable.importance),
                 RF.imp = (rf$variable.importance / sum(rf$variable.importance)),
                 fold = kx)
    
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
  mse <- sum((yhat_rf - y)^2 * W/sum(W))
  r2 <- 1 - (mse/sum((y - mean(y))^2 * W/sum(W)))
  ps <- cov.wt(cbind(y, yhat_rf), wt = W, cor = TRUE)$cor[1,2]
  rf_mod <- RF.table %>%
    dplyr::mutate(MSE = mse, R2 = r2, PearsonScore = ps) %>%
    dplyr::distinct(MSE, R2, PearsonScore) %>%
    dplyr::mutate(type = "RF")
  
  model_table %<>% dplyr::bind_rows(rf_mod)
  
  
  return(list(rf.fit = RF.table, mod.tab = model_table, preds = yhat_rf))
}
