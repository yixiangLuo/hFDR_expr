#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(tictoc)
library(hFDR)

source(here("R", "methods.R"))
source(here("R", "HIV", "utils.R"))
source(here("R", "plot.R"))

# read data
# x1, y1-7; x2, y1-6, x3, y1-3
# x1 y2; x3
indeces <- list("1" = 1:7, "2" = 1:6, "3" = 1:3)
# indeces <- list("1" = 7, "2" = 1:6, "3" = 1:3)
n_cores <- 14
method <- "lasso"  # LARS, lasso

## Getting the positions of the gene list and the TSM gene list for each drug
signal_genes <- list()
for (drug_class in c("PI", "NRTI", "NNRTI")){
  base_url <- 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
  tsm_url <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')
  tsm_df <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
  signals <- unlist(sapply(1:NROW(tsm_df), function(row){
    paste(tsm_df[row, 1], unlist(strsplit(tsm_df[row, 2], " ")), sep = ".")
  }))
  signal_genes[[drug_class]] <- rep(list(paste0("P", signals)),
                                    length(data[[drug_class]]$Y))
}
signal_genes <- do.call(c, signal_genes)

for(ind in 1:length(indeces)){
  drugClass_X_index <- as.integer(names(indeces[ind]))
  for(drug_y_index in indeces[[ind]]){
    print(paste("X", drugClass_X_index))
    print(paste("y", drug_y_index))
    load(here("data", "HIV", "HIV_data.RData"))
    data <- get_X_y(drugClass_X_index, drug_y_index)
    X <- data$X
    y <- data$y


    n <- NROW(X)
    p <- NCOL(X)

    if(method == "LARS"){
      tune_seq <- gene_lambda(X, y)
    } else if(method == "lasso"){
      tune_seq <- gene_lambda(X, y)
    } else stop()
    
    select <- function(X, y, lambda){
      select_variables(X, y, lambda, method = method)
    }
    pred_fit <- function(X.new, X, y, lambda){
      model_predict(X.new, X, y, lambda, method = method, use_mle = F)
    }
    
    cv.obj <- cv.model(X, y, tune_seq, pred_fit, select, nfold = 10)
    
    lambda.hFDR <- subset_lambda(cv.obj, n_tune = 40)
    tic()
    hFDR.obj <- hFDR(X, y, model = "gausslinear", select = method, lambda = lambda.hFDR,
                     psi = "pval", se = T, n_sample.se = 10, n_cores = 14)
    toc()

    tvals <- lm_to_t(y, X)
    pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")
    stat_signals <- which(pvals < 0.01 / p)
    
    selected <- select_variables(X, y, tune_seq, method)

    FDP <- sapply(1:length(tune_seq), function(tune_i){
      selected <- which(selected[, tune_i])
      get_FDP(selected, data$genes, data$drug_class, stat_signals)
    })
    errors <- list(FDP = FDP)

    info <- list(m1 = length(signal_genes[[paste0(data$drug_class, 1)]]), dim_X = dim(X),
                 drug_name = data$drug_name)
    # 
    # 
    expr_name <- paste0(drugClass_X_index, "_", drug_y_index, "-", data$drug_name)

    save(hFDR.obj, cv.obj, errors, method, info,
         file = here("data", "HIV", "expr", paste0(expr_name, ".RData")))
    # load(here("data", "HIV", "expr", paste0(expr_name, ".RData")))
    
    print(info)
    
    hFDR.obj$lambda <- hFDR.obj$lambda * n
    cv.obj$lambda <- cv.obj$lambda * n
    
    # plot for paper
    pdf(file = here("figs", "HIV",
                    paste0(expr_name, ".pdf")), width = 5, height = 3.8)  # width = 7, height = 5.5
    par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
    plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("fs"))-1, # ylab = "type I/II error",
         show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, log_x = T, xlab = "-log(lambda)", FDP.leg = expression(paste("TSM ", widehat(FDP))))
    dev.off()
    
    # ### plot for slides
    # 
    # halfslide.aspect.ratio <- 1.4
    # halfslide.width <- 4.5
    # halfslide.mar <- c(4,4,2,0.5)
    # halfslide.legend <- FALSE
    # halfslide.axes <- FALSE
    # 
    # pdf(file = here("figs","talk","hiv",paste0(expr_name, ".pdf")), width = halfslide.width, height = halfslide.width/halfslide.aspect.ratio)
    # par(mar = halfslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
    # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("fs"))-1,
    #      show_cv = T, log_cv = F, log_x = T, xlab = "-log(lambda)", 
    #      show_legend=halfslide.legend, show_extra_axes = halfslide.axes, main = data$drug_name)
    # dev.off()
    
  }
}


