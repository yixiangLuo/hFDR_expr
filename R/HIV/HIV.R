#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(tictoc)

source(here("R", "lasso_hFDR.R"))
source(here("R", "hFDR.R"))
source(here("R", "HIV", "utils.R"))

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
    # 
    # 
    # n <- NROW(X)
    # p <- NCOL(X)
    # 
    # if(method == "LARS"){
    #   tune_seq <- gene_lambda(X, y)
    # } else if(method == "lasso"){
    #   tune_seq <- gene_lambda(X, y)
    # } else stop()
    # 
    # cv.obj <- cv.predict(X, y, method, tune_seq,
    #                      nfold = 10, rescale = F, use_mle = F, fit_model = T)
    # 
    # tic()
    # hFDR.obj <- hFDR_on_cv(cv.obj, model = "fixed-X", n_tune = 40, n_cores = n_cores)
    # toc()
    # 
    # tvals <- lm_to_t(y, X)
    # pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")
    # stat_signals <- which(pvals < 0.01 / p)
    # 
    # FDP <- sapply(1:length(tune_seq), function(tune_i){
    #   selected <- which(cv.obj$selected[, tune_i])
    #   get_FDP(selected, data$genes, data$drug_class, stat_signals)
    # })
    # errors <- list(FDP = FDP)
    # 
    info <- list(m1 = length(signal_genes[[paste0(data$drug_class, 1)]]), dim_X = dim(X),
                 drug_name = data$drug_name)
    # 
    # 
    expr_name <- paste0(drugClass_X_index, "_", drug_y_index, "-", data$drug_name)
    # 
    # save(hFDR.obj, cv.obj, errors, method, info,
    #      file = here("data", "HIV", "expr", paste0(expr_name, ".RData")))
    load(here("data", "HIV", "expr", paste0(expr_name, ".RData")))
    
    print(info)
    
    # plot for paper
    pdf(file = here("figs", "HIV",
                    paste0(expr_name, ".pdf")), width = 5, height = 3.8)  # width = 7, height = 5.5
    par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
    plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1, # ylab = "type I/II error",
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
    # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
    #      show_cv = T, log_cv = F, 
    #      show_legend=halfslide.legend, show_extra_axes = halfslide.axes, main = data$drug_name)
    # dev.off()
    
  }
}

# expr_name
# load(here("data", "HIV", "expr", method, paste0("1_1-APV.RData")))
# colSums(cv.obj$selected)
# rownames(pvals)[which(cv.obj$selected[,5])]
# tvals$tvals[which(cv.obj$selected[,5])]
# pvals
# 
# for(ind in 1:length(indeces)){
#   drugClass_X_index <- as.integer(names(indeces[ind]))
#   for(drug_y_index in indeces[[ind]]){
#     print(paste("X", drugClass_X_index))
#     print(paste("y", drug_y_index))
#     load(here("data", "HIV", "HIV_data.RData"))
#     data <- get_X_y(drugClass_X_index, drug_y_index)
#     X <- data$X
#     y <- data$y
#     
#     
#     n <- NROW(X)
#     p <- NCOL(X)
#     
#     guess_null_by = "p-value"
#     
#     ## lasso selection + hFDR at certain lambda
#     sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p))
#     lambda <- 2*sigma / n
#     
#     # lasso
#     lasso.fit <- glmnet::glmnet(X, y, lambda = lambda,
#                                 intercept = F, standardize = F,
#                                 standardize.response = F, family = "gaussian")
#     
#     lasso.selected <- which(as.matrix(lasso.fit$beta) != 0)
#     
#     # # hFDR
#     # FDR.est <- hFDR.fixedX(X, y, lambda, guess_null_by,
#     #                        method = "lasso")
#     # 
#     # FDP.true <- get_FDP(lasso.selected, data$genes, data$drug_class)
#     # 
#     # print(paste0("Estiamted FDR: ", FDR.est))
#     # print(paste0("'True' FDP: ", FDP.true))
#     
#     
#     
#     ## glasso selection + hFDR at various rho
#     # Compute hFDR
#     lambda_max <- max(abs(t(X) %*% y)) / n
#     lambda_min <- sigma / n
#     n_lambda <- 20
#     lambda_seq <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
#     
#     lasso.fit <- glmnet::glmnet(X, y, lambda = lambda_seq,
#                                 intercept = F, standardize = F,
#                                 standardize.response = F, family = "gaussian")
#     
#     
#     # hFDR <- hFDR.fixedX(X, y, lambda_seq, guess_null_by,
#     #                     method = "lasso", n_cores = n_cores)
#     
#     glmnet.pack <- pack_glmnet(X, y, lasso.fit)
#     data.pack <- process_data(X, y)
#     hFDR.res <- calc_lasso_hFDR(glmnet.pack, data.pack, guess_null_by, couple = T, n_cores = n_cores)
#     hFDR <- hFDR.res$hFDR_lambda
#     DRj_lambda <- hFDR.res$DRj_lambda
#     
#     
#     glmnet.pack <- pack_glmnet(X, y, lasso.fit)
#     data.pack <- process_data(X, y)
#     dev <- est_lasso_hFDR_std(glmnet.pack, data.pack,
#                               guess_null_by,
#                               method = "Bootstrap.lasso_ols", # Bootstrap.lasso_ols, Bootstrap.nonpara
#                               measure = "quantile", couple = F, n_cores = n_cores)$deviate.est
#     # std <- 0
#     
#     
#     lasso.cv <- cv.glmnet(X, y, lambda = lambda_seq, alpha = 1,
#                           nfolds = 10, intercept = F, standardize = F,
#                           standardize.response = F, family = "gaussian")
#     
#     
#     hFDR.obj <- structure(list(X = X,
#                                y = y,
#                                glmnet.cv = lasso.cv,
#                                lambda = lambda_seq,
#                                hFDR = hFDR,
#                                hFDR.std = NA,
#                                hFDR.low = dev$low,
#                                hFDR.up = dev$up,
#                                name = "estimated FDR"),
#                           class = 'lasso.cv.hFDR')
#     
#     tvals <- data.pack$vjy_obs / sqrt(data.pack$RSS_X / (n-p))
#     pvals <- pvals_t(tvals, n-p, side = "two")
#     stat_signals <- which(pvals < 0.01 / p)
#     
#     plot.lasso.cv.hFDR <- function(hFDR.obj, sign.lambda = -1,...){
#       xlab <- expression(Log(lambda))
#       if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
#       
#       nlambda <- length(hFDR.obj$glmnet.cv$lambda)
#       lasso.FDP <- sapply(1:nlambda, function(lambda_index){
#         lasso.selected <- which(lasso.cv$glmnet.fit$beta[, lambda_index] != 0)
#         get_FDP(lasso.selected, data$genes, data$drug_class, stat_signals)
#       })
#       
#       
#       plot.range <- range(hFDR.obj$hFDR, hFDR.obj$hFDR.low, hFDR.obj$hFDR.up)
#       plot.range[1] <- max(plot.range[1], 0)
#       plot.range[2] <- min(plot.range[2], 1)
#       # plot.range <- range(hFDR.obj$hFDR)
#       plot.args <- list(x = sign.lambda * log(hFDR.obj$lambda),
#                         y = hFDR.obj$hFDR,
#                         xlim = range(sign.lambda*log(hFDR.obj$glmnet.cv$lambda)),
#                         ylim = plot.range,
#                         xlab = xlab, ylab = "estimated FDR and scaled CV MSE", type = "n")
#       new.args <- list(...)
#       if(length(new.args)) plot.args[names(new.args)] <- new.args
#       
#       cv.range <- abs(diff(range(c(hFDR.obj$glmnet.cv$cvm, hFDR.obj$glmnet.cv$cvlo, hFDR.obj$glmnet.cv$cvup))))
#       do.call("plot", plot.args)
#       lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
#             y = (hFDR.obj$glmnet.cv$cvm - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
#             col = "dodgerblue3")
#       # lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
#       #       y = (hFDR.obj$glmnet.cv$cvlo - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
#       #       col = "dodgerblue3", lty = 2)
#       # lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
#       #       y = (hFDR.obj$glmnet.cv$cvup - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
#       #       col = "dodgerblue3", lty = 2)
#       error.bars(sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
#                  (hFDR.obj$glmnet.cv$cvup - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
#                  (hFDR.obj$glmnet.cv$cvlo - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
#                  width = 0.01, col = alpha("dodgerblue3", 0.1))
#       lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
#             y = lasso.FDP,
#             col = "black", lty = 2)
#       error.bars(sign.lambda * log(hFDR.obj$lambda),
#                  hFDR.obj$hFDR.up, hFDR.obj$hFDR.low,
#                  width = 0.01, alpha("red", 0.3))
#       points(sign.lambda*log(hFDR.obj$lambda), hFDR.obj$hFDR,
#              pch = 20, col = "red")
#       axis(side = 3, at = sign.lambda*log(hFDR.obj$glmnet.cv$lambda),
#            labels = paste(hFDR.obj$glmnet.cv$nz), tick = FALSE, line = 0)
#       abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.min), lty = 3)
#       abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.1se), lty = 3)
#       legend("bottomright", inset = 0.05,
#              legend = c("CV MSE", "FDP", "hFDR"),
#              lty = c(1, 2, 1), lwd = c(1, 1, 0),
#              pch = c(26, 26, 19), col = c("dodgerblue3", "black", "red"))
#       invisible()
#     }
#     
#     
#     # plot
#     pdf(file = here("figs", "HIV",
#                     paste0(data$drug_class, drugClass_X_index, "-Y",
#                            drug_y_index, ".pdf")),
#         width = 7,height = 5.5)
#     plot.lasso.cv.hFDR(hFDR.obj)
#     dev.off()
#     
#     DRj_lambda.data <- data.frame(hFDR = NULL, j = NULL, log_lambda = NULL)
#     for(j in 1:p){
#       for(lambda_i in 1:n_lambda){
#         DRj_lambda.data <- rbind(DRj_lambda.data,
#                                  data.frame(hFDR = sum(DRj_lambda[1:j, lambda_i]), j = as.character(p-j), log_lambda = log(lambda_seq[lambda_i])))
#       }
#     }
#     plot <- ggplot(DRj_lambda.data) +
#       geom_ribbon(aes(x = -log_lambda, ymax = hFDR, fill = j), ymin = 0, alpha = 0.2) +
#       geom_line(aes(x = -log_lambda, y = hFDR, group = j), alpha = 1, color = "black", size = 0.2) +
#       guides(fill = F, color = F) +
#       theme_bw()
#     
#     ggsave(filename = here("figs", "HIV",
#                            paste0(data$drug_class, drugClass_X_index, "-Y",
#                                   drug_y_index, "-detail.pdf")),
#            plot, width = 7, height = 5.5)
#     
#     m1 <- length(signal_genes[[paste0(data$drug_class, 1)]])
#     print(m1)
#     print(dim(X))
#   }
# }



### Some code Will was using to audit answers for specific experiments

# tic()
# hFDR.obj <- hFDR_on_cv(cv.obj, model = "fixed-X", n_tune = 40, n_cores = n_cores, decompose=TRUE, skip_se = TRUE)
# toc()
# matplot(t(hFDR.obj$hFDR),type="l")
# abline(v=8)
# 
# selected <- cv.obj$selected[,cv.obj$tune_seq %in% hFDR.obj$tune_seq]
# selected
# (selected <- which(selected[,8]))
# sort(hFDR.obj$hFDR[,7], decreasing=TRUE)
# order(hFDR.obj$hFDR[,7], decreasing=TRUE)
# pvals[order(hFDR.obj$hFDR[,7], decreasing=TRUE)[1:5],]
# signal_genes[[paste0(data$drug_class,1)]]
# hFDR.obj
# 
# (top.vars <- order(hFDR.obj$hFDR[,7], decreasing=TRUE)[1:3])
# 
# cor(data$X[,top.vars])
# mod <- lm(data$X[,148] ~ data$X[,-148])
# cor(mod$fitted.values, data$X[,148])
# plot(mod$fitted.values, data$X[,148])
# 
# mod <- lm(data$X[,147] ~ data$X[,-147])
# cor(mod$fitted.values, data$X[,147])
# plot(mod$fitted.values, data$X[,147])
# 
# mod <- lm(data$X[,275] ~ data$X[,-275])
# cor(mod$fitted.values, data$X[,275])
# plot(mod$fitted.values, data$X[,275])
# 
# summary(cor(data$X)[upper.tri(cor(data$X))])
# 
# apply(data$X[,order(hFDR.obj$hFDR[,7], decreasing=TRUE)[1:3]],2,table)
# 
# data$X


