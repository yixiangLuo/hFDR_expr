#### modified from https://github.com/lihualei71/dbhPaper/blob/master/R/dBH_lm_HIV_expr.R

library(here)
library(glmnet)   # lasso
library(latex2exp)

if (!file.exists(here("data", "HIV", "HIV_data.RData"))){
  source(here("R", "HIV", "HIV_preprocess.R"))
}
load(here("data", "HIV", "HIV_data.RData"))

get_X_y <- function(drugClass_X_index = 1, drug_y_index = 1){
  X <- data[[drugClass_X_index]]$X
  Y <- data[[drugClass_X_index]]$Y
  y <- Y[[drug_y_index]]
  
  
  ## Log-transform the drug resistance measurements.
  y <- log(y)
  
  ## Remove patients with missing measurements.
  missing <- is.na(y)
  y <- y[!missing]
  X <- X[!missing,]
  
  ## Remove predictors that appear less than 3 times.
  X <- X[,colSums(X) >= 3]
  
  ## Remove duplicate predictors.
  X <- X[, colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  X <- scale(X, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2) / (NROW(X)-1)))
  y <- y - mean(y)
  
  ## Get names
  genes <- colnames(X)
  drug_class <- names(data)[drugClass_X_index]
  drug_name <- names(Y)[drug_y_index]
  
  return(list(X = X, y = y,
              genes = genes, drug_class = drug_class, drug_name = drug_name))
}

get_FDP <- function(rejs, genes, drug_class, stat_signals = NULL){
  signals <- signal_genes[[paste0(drug_class, 1)]]
  signals <- c(signals, genes[stat_signals])
  rejs <- genes[rejs]
  false_discoveries <- sum(!(rejs %in% signals))
  FDP <- false_discoveries / max(1, length(rejs))
  
  return(FDP)
}

gene_lambda <- function(X, y, n_lambda = 100){
  n <- NROW(X)
  p <- NCOL(X)
  
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
  
  lambda_max <- max(abs(t(X) %*% y)) / n
  sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p-1))
  lambda_min <- sigma / n
  lambda_seq <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
  
  return(lambda_seq)
}


