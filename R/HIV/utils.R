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


plot.ROC <- function(ROC, log_cv = T){
  
  method_colors <- c("#e41a1c", "orange1", "#377eb8")
  method_names <- c("LARS", "elastic", "FS")
  legend_names <- c("LARS", "Elastic net", "Forward stepwise")
  
  par(mar = c(4,4,3,3), pty = "s")
  
  if(log_cv) plot(x = ROC[[1]]$hFDR, y = ROC[[1]]$MSE_CV, xlab = "", ylab = "", type = "n", log = "y")
  else plot(x = ROC[[1]]$hFDR, y = ROC[[1]]$MSE_CV, xlab = "", ylab = "", type = "n")
  
  leg <- c()
  lty <- c()
  lwd <- c()
  col <- c()
  
  for(i in 1:length(ROC)){
    method <- method_names[i]
    errors <- ROC[[method]]
    lines(x = errors$hFDR, y = errors$MSE_CV, col = method_colors[i], lwd = 2)
    leg <- c(leg, legend_names[i])
    lty <- c(lty, 1)
    lwd <- c(lwd, 2)
    col <- c(col, method_colors[i])
  }
  
  title(xlab = TeX("\\widehat{FDR}$"), ylab = "CV MSE", line = 2.5, cex.lab = 1)
  
  legend("topright", inset = 0.05,
         legend = leg,
         lty = lty, lwd = lwd,
         col = col)
  
  invisible()
}


# ## Getting the positions of the gene list and the TSM gene list for each drug
# signal_genes <- list()
# for (drug_class in c("PI", "NRTI", "NNRTI")){
#   base_url <- 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
#   tsm_url <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')
#   tsm_df <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
#   signal_genes[[drug_class]] <- rep(list(tsm_df[, 1]),
#                                     length(data[[drug_class]]$Y))
# }
# signal_genes <- do.call(c, signal_genes)
# 
# get_position <- function(x){
#   sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)
# }
# 
# get_FDP <- function(rejs, genes, drug_class){
#   signals <- signal_genes[[paste0(drug_class, 1)]]
#   # rejs <- unique(get_position(genes[rejs]))
#   # false_discoveries <- length(setdiff(rejs, signals))
#   rejs <- (get_position(genes[rejs]))
#   false_discoveries <- sum(!(rejs %in% signals))
#   FDP <- false_discoveries / max(1, length(rejs))
#   
#   return(FDP)
# }




