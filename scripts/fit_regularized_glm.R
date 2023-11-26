#!/usr/bin/env Rscript
# $title fit regularized glm to molecular features
# adapted from supplementary material of Horvath_ Genome Biol_2013 paper
# $input predictor_matrix, response_vector
# $output coefficient of predictor, test MSE
# $author giuseppe
# $date Apr 2020

library(glmnet)
library(doMC)
registerDoMC(cores=4)

# define regularized glm function
reg_glm_wrap <- function(X.m = X.m,
                        y.v = y.v,
                        prop_train = .7,
                        family = "gaussian", # "lasso", "ridge"
                        penalty = "elastic.net",
                        standardidze = F,
                        nfolds = 10,
                        nlambda = 100,
                        directory = directory,
                        name.plt = name.plt) {
  # match obs of X.m & y.v
  rownames(X.m) <- sample_name(rownames(X.m), 4)
  names(y.v) <- sample_name(names(y.v), 4)
  X.m <- X.m[rownames(X.m) %in% names(y.v),]
  y.v <- y.v[names(y.v) %in% rownames(X.m)]
  X.m <- X.m[match(names(y.v), rownames(X.m)),]
  # remove obs with missing response value
  del.idx <- is.na(y.v)
  y.v <- y.v[!del.idx]
  X.m <- X.m[!del.idx,]
  
  # set penalty type
  alpha <- switch(penalty, elastic.net = .5, lasso = 1, ridge = 0)
  # split obs into training set & test set
  train <- sample(1:nrow(X.m), nrow(X.m)*prop_train)
  test <- -train
  y_test <- y.v[test]
  
  # fit the model to the training data
  fit_train <- glmnet(X.m[train, ], y.v[train], family = family, 
                      alpha = alpha, nlambda = nlambda, standardidze = standardidze)
  # estimate lambda parameter with CV
  fit_cv <- cv.glmnet(X.m[train, ], y.v[train], family = family, alpha = alpha, parallel = T)
  pdf(file = paste0(directory,"scatterplot_MSE_logLambda_CV_",penalty,"_",name.plt,".pdf"),
      width = 8, height = 8, pointsize = 14)
  plot(fit_cv)
  dev.off()
  bestlam <- fit_cv$lambda.min
  # get coefficient of predictor
  coeff.x <- coef(fit_train, s = bestlam)
  # compute test MSE
  pred_test <- predict(fit_train, s = bestlam, newx = X.m[test,])
  test_MSE <- mean((pred_test-y_test)^2)
  res <- list(fit_train, fit_cv, coeff.x, test_MSE)
  names(res) <- c("fit_train", "fit_cv", "coeff.x", "test_MSE")
  return(res)
}
