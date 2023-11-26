#!/usr/bin/env Rscript
# $title fit regression tree to molecular features
# $input predictor_matrix, response_vector
# $output coefficient of predictor, test MSE
# $author giuseppe
# $date Apr 2020

library(randomForest)

# define random forest regression function
regression.tree_wrap <- function(X.m = X.m,
                         y.v = y.v,
                         ylab = ylab,
                         prop_train = .7,
                         tree.type = "randomForest",  # bagging
                         ntree = 1e3,
                         importance = T,
                         n.var = 10,
                         width = 8,
                         height = 8,
                         point.size_plt = 14,
                         directory = directory,
                         name.plt = name.plt) {
  colnames(X.m) <- gsub("-","_", colnames(X.m), fixed = T)
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

  # split obs into training set & test set
  train <- sample(1:nrow(X.m), nrow(X.m)*prop_train)
  test <- -train
  y_test <- y.v[test]
  
  # fit the model to the training data
  #mtry <- ifelse(tree.type == "randomForest", ncol(X.m)/3, ncol(X.m))
  mtry <- ifelse(tree.type == "randomForest", round(sqrt(ncol(X.m))), ncol(X.m))
  fit_train <- randomForest(X.m[train,], y.v[train], mtry = mtry, 
                            ntree = ntree, importance = importance)
  # compute test MSE
  yhat_test <- predict(fit_train, newdata=X.m[-train ,])
  test_MSE <- mean((yhat_test-y_test)^2)
  # retrieve importance of predictors
  importance <- importance(fit_train)
  res <- list(fit_train, importance, test_MSE)
  names(res) <- c("fit_train", "importance_predictor", "test_MSE")
  # plot importance measures
  pdf(file = paste0(directory,"plot_importance_predictor_",tree.type,"_",name.plt,".pdf"),
      width = width, height = height, pointsize = point.size_plt)
  varImpPlot(fit_train, n.var = n.var, main = "Variable importance plot")
  dev.off()
  return(res)
}
