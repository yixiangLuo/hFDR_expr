#!/usr/bin/env Rscript

library(here)
library(tictoc)
library(tidyverse)

source(here("R", "lasso_hFDR.R"))
source(here("R", "hFDR.R"))
source(here("R", "HIV", "utils.R"))

# read data
# x1, y1-7; x2, y1-6, x3, y1-3
# x1 y2; x3


indeces <- list("1" = 1:7, "2" = 1:6, "3" = 1:3)

for(ind in 1:length(indeces)){
  drugClass_X_index <- as.integer(names(indeces[ind]))
  for(drug_y_index in indeces[[ind]]){
  
    print(paste("X", drugClass_X_index))
    print(paste("y", drug_y_index))
    
    res <- get_X_y(drugClass_X_index, drug_y_index)
    X <- res$X
    y <- res$y
    drug_name <- res$drug_name
  
    n <- NROW(X)
    p <- NCOL(X)
  
    guess_null_by = "p-value"
  
    ## glasso selection + hFDR at various rho
    # Compute hFDR
    lambda_seq <- gene_lambda(X, y, n_lambda = 100)
  
    result <- data.frame()
    ROC <- list()
  
    # method <- "lasso"
    # print(method)
    # hFDR <- hFDR.fixedX(X, y, tune_seq = lambda_seq, guess_null_by, method, n_cores = 14)
    # cv.obj <- cv.predict(X, y, method, tune_seq = lambda_seq, nfold = 10, use_mle = F)
    # result <- rbind(result, data.frame(hFDR = hFDR, cv.obj = cv.obj$cv.mean,
    #                                    method = method, drug_name = drug_name))
    
    method <- "LARS"
    print(method)
    tic()
    hFDR <- hFDR.fixedX(X, y, tune_seq = lambda_seq, guess_null_by, method, n_cores = 14)
    cv.obj <- cv.predict(X, y, method, tune_seq = lambda_seq, nfold = 10, use_mle = F)
    ROC[[method]] <- list(hFDR = hFDR, MSE_CV = cv.obj$cv.mean)
    toc()
    
    method <- "elastic"
    print(method)
    tic()
    hFDR <- hFDR.fixedX(X, y, tune_seq = lambda_seq, guess_null_by, method, n_cores = 14)
    cv.obj <- cv.predict(X, y, method, tune_seq = lambda_seq, nfold = 10, use_mle = F)
    ROC[[method]] <- list(hFDR = hFDR, MSE_CV = cv.obj$cv.mean)
    toc()
    
    method <- "FS"
    print(method)
    tic()
    hFDR <- hFDR.fixedX(X, y, tune_seq = 1:round(0.9*p), guess_null_by, method, n_cores = 1)
    cv.obj <- cv.predict(X, y, method, tune_seq = 1:round(0.9*p), nfold = 10, use_mle = F)
    ROC[[method]] <- list(hFDR = hFDR, MSE_CV = cv.obj$cv.mean)
    toc()
    
    ROC <- structure(ROC, class = 'ROC')
  
    expr_name <- paste0(drugClass_X_index, "_", drug_y_index, "-", drug_name)
    
    save(ROC,
         file = here("data", "HIV", "ROC", paste0(expr_name, ".RData")))
    # load(here("data", "HIV", "ROC", paste0(expr_name, ".RData")))
    
    # plot
    pdf(file = here("figs", "HIV", "ROC",
                    paste0("ROC-", expr_name, ".pdf")), width = 5, height = 5)  # width = 7, height = 5.5
    par(mar = c(4,4,3,3), pty = "s")
    plot(ROC)
    dev.off()
  }
}

# method_colors <- c("#e41a1c", "orange1", "#377eb8")
# method_names <- c("LARS", "elastic", "FS")
# legend_names <- c("LARS", "Elastic net", "Forward stepwise")
# names(method_colors) <- method_names
# 
# plot <- ggplot(result, aes(x = hFDR, y = log(cv.obj), color = method)) +
#   # geom_line() +
#   geom_path() +
#   scale_color_manual(values = method_colors, labels = legend_names, breaks = method_names) +
#   facet_grid(NULL, vars(drug_name), scales = "free") +
#   theme_bw() +
#   theme(aspect.ratio = 1,
#         panel.grid = element_blank(),
#         strip.text = element_text(size = 15),
#         axis.title = element_text(size = 13),
#         axis.text = element_text(size = 10),
#         legend.position = "right",
#         legend.title=element_text(size=9),
#         legend.text=element_text(size=9)) +
#   labs(x = TeX("\\widehat{FDR}$"), y = "CV MSE")
# 
# ggsave(filename = here("figs", "HIV", "ROC",
#                        paste0("ROC-", expr_name, ".pdf")),
#        plot, width = 5, height = 4)



