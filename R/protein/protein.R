#!/usr/bin/env Rscript

library(here)
library(readxl)
library(igraph)
library(tidyverse)
library(tictoc)

source(here("R", "hFDR.R"))
source(here("R", "protein", "utils.R"))

method <- "glasso"
n_cores <- 14

# plot DAG and moralized graph
pdf(file = here("figs", "protein", "DAG.pdf"), width = 6.5, height = 6)  # width = 7, height = 5.5
plot.graph(adjmat.DAG, "directed")
dev.off()
pdf(file = here("figs", "protein", "UDG.pdf"), width = 6.5, height = 6)  # width = 7, height = 5.5
plot.graph(adjmat.true, "undirected")
dev.off()

for(expr in 0:14){
  print(expr)
  # read data
  if(expr > 0){
    X <- read_data(exprs = expr)
  } else{
    X <- do.call(rbind, lapply(1:14, function(i){
      read_data(exprs = i)
    }))
  }
  
  expr_name <- paste0("expr_", expr)
  
  n <- NROW(X)
  p <- NCOL(X)
  pair_num <- p*(p-1)/2

  tune_seq <- gene_rho(X, n_tune = 100)

  cv.obj <- cv.ggraph(X, method, tune_seq, nfold = 10, use_mle = F, fit_model = T)

  tic()
  hFDR.obj <- hFDR_on_cv(cv.obj, model = "graphGauss", n_tune = 40, n_cores = n_cores)
  toc()

  pvals <- guess_null(X, y = NA, tune_seq, Xcov.true = NA, method = method, model = "graphGauss")$pvals
  stat_signals <- which(pvals < 0.01 / pair_num)
  nulls <- null_pairs[!(null_pairs %in% stat_signals)]
  H0 <- (1:pair_num) %in% nulls

  FDP <- calc_sel_FDP(cv.obj$selected, H0)
  errors <- list(FDP = FDP)

  info <- list(dim_X = dim(X))

  print(info)
  
  
  save(hFDR.obj, cv.obj, errors, method, info,
       file = here("data", "protein", "expr", paste0(expr_name, ".RData")))
  # load(here("data", "protein", "expr", paste0(expr_name, ".RData")))
  
  # plot
  pdf(file = here("figs", "protein", paste0(expr_name, ".pdf")), width = 5, height = 3.8)  # width = 7, height = 5.5
  par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
  plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1, # ylab = "type I/II error",
       show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, log_x = T, xlab = "-log(lambda)", FDP.leg = expression(paste("DAG ", widehat(FDP))))
  dev.off()
  
  
  
  # fullslide.aspect.ratio <- 1.6
  # fullslide.width <- 6
  # fullslide.mar <- c(4,4,3,4)
  # 
  # pdf(file = here("figs","talk","protein",paste0(expr_name, "_noFDR.pdf")), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_FDR = F, show_hFDR = F, show_FPR = F, log_cv = F, leg_bty="n",
  #      FDP.leg = expression(paste("Valid'n. ", widehat(FDP))))
  # dev.off()
  # 
  # pdf(file = here("figs","talk","protein",paste0(expr_name, "_FDR_nohFDR.pdf")), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = F, show_FPR = F, log_cv = F, leg_bty="n",
  #      FDP.leg = expression(paste("Valid'n. ", widehat(FDP))))
  # dev.off()
  # 
  # pdf(file = here("figs","talk","protein",paste0(expr_name, "_hFDR.pdf")), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, leg_bty="n",
  #      FDP.leg = expression(paste("Valid'n. ", widehat(FDP))))
  # dev.off()
  
}


# # glasso selection + hFDR at certain rho
# rho <- 0.078
# 
# # glasso
# S <- var(X)
# glasso.res <- glasso(S, rho)
# Theta <- glasso.res$wi
# Sigma <- glasso.res$w
# glasso.selected <- which(Theta[upper.tri(Theta)] != 0)
# 
# plot.graph(Theta)
# print(sum(abs(Theta[upper.tri(Theta)])))
# 
# # hFDR
# FDR.est <- hFDR.glasso(X, rho, guess_null_by = "p-value",
#                        method = "glasso", mc_size = 100)[1]
# 
# FDP.true <- FDP(glasso.selected)
# 
# print(paste0("Estiamted FDR: ", FDR.est))
# print(paste0("'True' FDP: ", FDP.true))
