
library(here)

library(tidyverse)
library(glmnet)
library(tictoc)

source(here("R", "lasso_hFDR.R"))
source(here("R", "hFDR.R"))
source(here("R", "utils.R"))
source(here("R", "methods.R"))


### scenario = "Independent X" ## "Independent X", "Correlated X", "Forward Stepwise"

for(scenario in c("Independent X", "Correlated X", "Forward Stepwise")) {
  
  generate_data <- function(scenario, way = "linear", setting_seed = 1, noise_seed = -1){
    if(way == "linear"){
      
      p <- 200
      n <- 3*p
      pi1 <- 20 / p
      mu1 <- 6 / sqrt(n)
      
      set.seed(setting_seed)
      
      ### change design matrix settings here. see utils.R > gene_X function for more available settings
      X <- gene_X(ifelse(scenario == "Correlated X", "X_AR_Strong", "IID_Normal"), 
                  n, p, setting_seed, scaling = F, model_X = T)$X  # "IID_Normal", "MCC_Block", "X_AR_Strong" 
      beta <- genmu(p, pi1, mu1, "random", 1) * (1 + rexp(p, 1))/2
      H0 <- beta == 0
      
      set.seed(noise_seed)
      y <- X %*% beta + rnorm(n)
    } else if(way == "hbeta"){ # generate data by sample OLS \hat\beta first and then recover X, y
      set.seed(setting_seed)
      
      block.size <- 10
      n.blocks <- 10
      m <- block.size * n.blocks
      m1 <- 5
      nn.indices <- sample(1:m, size = m1, replace=FALSE)
      nn.pattern <- rep(0,m)
      nn.pattern[nn.indices] <- 1
      
      mu <- 1.5 * nn.pattern # * (1 + rexp(p, 1))/2
      
      rho <- 0.8
      
      set.seed(noise_seed)
      z.block <- rnorm(n.blocks)
      z <- mu + rho * rep(z.block, each=block.size) + sqrt(1-rho^2) * rnorm(m)
      
      p <- length(z)
      n <- 3*p
      H0 <- nn.pattern == 0
      
      blockSigma <- matrix(rho^2, block.size, block.size)
      diag(blockSigma) <- 1
      beta_cov <- as.matrix(diag(n.blocks) %x% blockSigma)
      cov_mat <- solve(beta_cov)
      
      R <- chol(cov_mat)
      basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
      X <- basis %*% R
      
      proj_y <- X %*% z
      
      y <- rnorm(n, sd = 1)
      y <- proj_y + y - lm(formula = y ~ X + 0)$fitted.values
    }
    
    return(list(X = X, y = y, H0 = H0, beta = beta))
  }
  
  # generate lamabda sequence for lasso
  gene_lambda <- function(X, y, scaler = 1){
    n <- NROW(X)
    p <- NCOL(X)
    
    y <- scale(y, center = T, scale = F)
    X <- scale(X, center = T, scale = F)
    X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
    
    lambda_max <- max(abs(t(X) %*% y)) / n * scaler
    sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p-1))
    lambda_min <- sigma / n
    n_lambda <- 100
    lambda_seq <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
    
    return(lambda_seq)
  }
  
  setting_seed <- 1
  data <- generate_data(scenario = scenario, way = "linear", setting_seed, noise_seed = -1)
  X <- data$X
  y <- data$y
  H0 <- data$H0
  beta <- data$beta
  
  n <- NROW(X)
  p <- NCOL(X)
  
  
  ### change selection method here, for Gaussian linear model, currently support "lasso", "elastic", "FS", "LARS", see hFDR.R
  if(scenario == "Forward Stepwise") {
    method <- "FS"
    tune_seq <- 1:round(0.75*p)   # number of selections
  } else {
    method <- "lasso"
    tune_seq <- gene_lambda(X, y, scaler = 1.4)
  }
  
  cv.obj <- cv.predict(X, y, method, tune_seq,
                       nfold = 10, rescale = F, use_mle = F, fit_model = T)
  
  
  FDP <- calc_sel_FDP(cv.obj$selected, H0)
  FPP <- 1-calc_sel_power(cv.obj$selected, H0)
  error.samples <- sapply(1:100, function(mc_i){
    data.mc <- generate_data(scenario = scenario, way = "linear", setting_seed, noise_seed = -mc_i)
    X.mc <- data.mc$X
    y.mc <- data.mc$y
    H0.mc <- data.mc$H0
    select.res <- select_variables(X.mc, y.mc, tune_seq, method)
    FDP <- calc_sel_FDP(select.res, H0.mc)
    FPP <- 1-calc_sel_power(select.res, H0.mc)
    return(c(FDP, FPP))
  })
  n_tune <- length(tune_seq)
  FDR <- rowMeans(error.samples[1:n_tune, ])
  FPR <- rowMeans(error.samples[(n_tune+1):(2*n_tune), ])
  
  errors <- list(FDP = FDP, FPP = FPP, FDR = FDR, FPR = FPR)
  
  tic()
  hFDR.obj <- hFDR_on_cv(cv.obj, model = "fixed-X", n_tune = ifelse(method=="FS",60,40), n_cores = 16)
  toc()
  
  tic()
  hFDR.decomp <- hFDR_on_cv(cv.obj, model = "fixed-X", n_tune = ifelse(method=="FS",60,40), n_cores = 16, decompose = TRUE, skip_se = TRUE)
  toc()
  
  #### These should match
  colSums(hFDR.decomp$hFDR) - hFDR.obj$hFDR
  
  
  #### draw figures for paper
  fullcolumn.aspect.ratio <- 1.6
  fullcolumn.width <- 6
  fullcolumn.mar <- c(4,4,3,4)
  
  halfcolumn.aspect.ratio <- 1.4
  halfcolumn.width <- 4.5
  halfcolumn.mar <- c(4,4,3,4)
  halfcolumn.legend <- FALSE
  halfcolumn.axes <- FALSE
  
  pdf(file = here("figs","intro", paste0(gsub(" ","_",scenario), ".pdf")), width = fullcolumn.width, height = fullcolumn.width/fullcolumn.aspect.ratio)
  par(mar = fullcolumn.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
       show_cv = T, log_cv = F)  # leg_bty="n"
  dev.off()
  
  
  print(sum(cv.obj$selected[, cv.obj$ind.min]))
  print(sum(cv.obj$selected[, cv.obj$ind.1se]))
  
  
  # ### draw figures for slides. use show_cv, show_hFDR, show_FPR to control whether to show cv curve, hFDR, FDR/FDP if provided, and FPR/FPP. See hFDR.R > plot.hFDR
  # 
  # fullslide.aspect.ratio <- 1.6
  # fullslide.width <- 6
  # fullslide.mar <- c(4,4,3,4)
  # 
  # halfslide.aspect.ratio <- 1.4
  # halfslide.width <- 4.5
  # halfslide.mar <- c(4,4,2,0.5)
  # halfslide.legend <- FALSE
  # halfslide.axes <- FALSE
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"fullslide_noFDR.pdf"), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = F, show_FPR = F, show_FDR = F, log_cv = F, leg_bty="n")
  # dev.off()
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"fullslide_FDR_nohFDR.pdf"), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = F, show_FPR = F, log_cv = F, leg_bty="n")
  # dev.off()
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"fullslide_hFDR.pdf"), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, leg_bty="n")
  # dev.off()
  # 
  # # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"fullslide_hFDR_decomp.pdf"), width = fullslide.width, height = fullslide.width/fullslide.aspect.ratio)
  # # par(mar = fullslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  # #      show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, leg_bty="n")
  # # dev.off()
  # 
  # 
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"halfslide_noFDR.pdf"), width = halfslide.width, height = halfslide.width/halfslide.aspect.ratio)
  # par(mar = halfslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = F, show_FPR = F, show_FDR = F, log_cv = F, 
  #      show_legend=halfslide.legend, show_extra_axes = halfslide.axes, main = scenario)
  # dev.off()
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"halfslide_FDR_nohFDR.pdf"), width = halfslide.width, height = halfslide.width/halfslide.aspect.ratio)
  # par(mar = halfslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = F, show_FPR = F, log_cv = F, leg_bty="n", 
  #      show_legend=halfslide.legend, show_extra_axes = halfslide.axes, main = scenario)
  # dev.off()
  # 
  # pdf(file = here("figs","talk","intro",gsub(" ","_",scenario),"halfslide_hFDR.pdf"), width = halfslide.width, height = halfslide.width/halfslide.aspect.ratio)
  # par(mar = halfslide.mar, cex = 1, cex.axis = 1, cex.lab = 1)
  # plot(hFDR.obj, cv.obj, errors, sign.tune = 2*(method %in% c("FS"))-1,
  #      show_cv = T, show_hFDR = T, show_FPR = F, log_cv = F, leg_bty="n", 
  #      show_legend=halfslide.legend, show_extra_axes = halfslide.axes, main = scenario)
  # dev.off()
  
}
