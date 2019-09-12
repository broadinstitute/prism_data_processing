# Script to run multivariate analysis on data using random forest and enet

# import necessary libraries and functions using MTS_functions.R
source("./MTS_functions.R")


# function to do multivariate analysis on MTS data
# takes named vector of responses (y) and feature matrix (X)
# k = number of cross validation cycles
# n = number of features to test for each run (e.g. top 250)
biomarkers_cv <- function(y, X, k = 10, n = 250){
  y <- y[is.finite(y)]  # only finite values
  cl <- sample(intersect(rownames(X), names(y)))  # overlapping rows
  X <- scale(X[cl,]); y = y[cl]; N = floor(length(cl)/k)
  colnames(X) %<>% make.names()
  yhat_rf <- y; yhat_enet <- y
  
  SS = tibble()  # output
  
  # run cross validation k times
  for (kx in 1:k) {
    print(paste0( (kx - 1) * N, " - ", kx * N))  # for debugging (loading bar)
    test <- cl[((kx - 1) * N + 1):(kx * N)]  # select test set
    train <- setdiff(cl, test)  # everything else is training
    X_train <- X[train,]
    
    # assumes no NAs in X
    X_train <- X_train[,apply(X_train,2,var) > 0]
    X_test <- X[test,]
    
    # select top n correlated features in X
    X_train <- X_train[,rank(-abs(cor(X_train,y[train])))<= n]
    
    # makes RF and ENet models
    rf <- ranger(y ~ .,
                data = as.data.frame(cbind(X_train, y = y[train])),
                importance = "impurity")
    
    enet <- cv.glmnet(x = X_train,
                     y = y[train], 
                     alpha = .5)
    
    # generate predictions
    yhat_enet[test] <- 
      predict.cv.glmnet(enet, newx = X_test[, colnames(X_train)])
    yhat_rf[test] <-
      predict(rf, data = as.data.frame(X_test[, colnames(X_train)]))$predictions
    
    # get variable importance metrics
    a <- coef(enet)
    ss <- tibble(feature = names(rf$variable.importance),
                RF.imp = rf$variable.importance / sum(rf$variable.importance),
                RF.mse = mean((yhat_rf[test] - y[test])^2)
    ) %>%
      dplyr::full_join(tibble(feature = a@Dimnames[[1]][a@i+1],
                              enet.coef = a@x)) %>%
      dplyr::mutate(enet.coef = ifelse(is.na(enet.coef), 0, enet.coef),
                    enet.mse = mean((yhat_enet[test] - y[test])^2, na.rm = T),
                    fold = kx)
    
    # add this model to the output
    SS %<>%
      bind_rows(ss)
  }
  
  # separate into RF and ENet tables
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
    dplyr::filter((RF.imp.stability > 0.5), feature != "(Intercept)") %>%
    dplyr::mutate(RF.MSE = mean(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                                na.rm =T),
                  RF.MSE.se = sd(dplyr::distinct(SS, RF.mse, fold)$RF.mse,
                                 na.rm = T),
                  RF.R2 = 1 - RF.MSE/var(y[1:(k*N)], na.rm = T),
                  RF.PearsonScore = cor(y[1:(k*N)], yhat_rf[1:(k*N)],
                                        use = "pairwise.complete.obs"))
  
  
  enet.importances <- SS %>% 
    dplyr::distinct(feature, enet.coef, fold) %>%
    reshape2::acast(feature ~ fold, value.var = "enet.coef")
  
  enet.table <- tibble(feature = rownames(enet.importances),
                      enet.coef.mean = enet.importances %>%
                        apply(1, function(x) mean(x[x!= 0], na.rm = T)),
                      enet.coef.sd =   enet.importances %>%
                        apply(1, function(x) sd(x[x!= 0], na.rm = T)),
                      enet.coef.stability =   enet.importances %>%
                        apply(1, function(x) mean(x != 0, na.rm = T))) %>%
    dplyr::filter((enet.coef.stability >= 0.7), feature != "(Intercept)") %>%
    dplyr::mutate(enet.MSE = mean(dplyr::distinct(SS, enet.mse, fold)$enet.mse, 
                                  na.rm = T),
                  enet.MSE.se = sd(dplyr::distinct(SS, enet.mse, fold)$enet.mse, 
                                   na.rm = T),
                  enet.R2 = 1 - enet.MSE/var(y[1:(k*N)], na.rm = T),
                  enet.PearsonScore = cor(y[1:(k*N)], yhat_enet[1:(k*N)], 
                                          use = "pairwise.complete.obs"))
  
  return(list(RF.table, enet.table, yhat_rf, yhat_enet))
}

X <- data.table::fread()  # INSERT path to feature table
df <- data.table::fread()  # INSERT path to response table

y <- df$response; names(y) <- df$ccle_name

results <- biomarkers_cv(y, X)
names(results) <- c("RFImp", "ENetImp", "RFPred", "ENetPred")
