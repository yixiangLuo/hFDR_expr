library(glmnet)
source(here("R", "homotopy_lasso.R"))

lasso.cv.hFDR <- function(X, y, lasso.cv, n_cores = 1){
  if(!all(class(lasso.cv) == "cv.glmnet")){
    stop()
  }
  
  n <- NROW(X)
  p <- NCOL(X)
  data.pack <- process_data(X, y)
  
  sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
  
  nlambda <- length(lasso.cv$lambda)
  nsample <- 20
  lambda0.index <- which.min(abs(lasso.cv$lambda * n - 2*sigma_hat))
  lambda.end <- max(lambda0.index, lasso.cv$index + nlambda / nsample)
  lambda.end <- min(nlambda, round(lambda.end))
  
  ind.sample <- round(seq(from = 1, to = lambda.end, length.out = nsample))
  ind.sample <- sort(unique(c(ind.sample, lasso.cv$index)))
  lambda <- lasso.cv$lambda[ind.sample]
  lasso.fit <- list(lambda = lambda, beta = lasso.cv$glmnet.fit$beta[, ind.sample])
  
  glmnet.pack <- pack_glmnet(X, y, lasso.fit)
  data.pack <- process_data(X, y)
  
  hFDR <- calc_lasso_hFDR(glmnet.pack, data.pack, guess_null_by = "p-value", couple = T, n_cores = n_cores)$hFDR_lambda
  deviate.est <- est_lasso_hFDR_std(glmnet.pack, data.pack,
                                    guess_null_by = "p-value", couple = F, n_cores = n_cores)$deviate.est
  
  
  
  if(deviate.est$type == "add"){
    hFDR.low <- hFDR - deviate.est$low
    hFDR.up <- hFDR - deviate.est$up
  } else if(deviate.est$type == "mult"){
    hFDR.low <- hFDR * deviate.est$low
    hFDR.up <- hFDR * deviate.est$up
  } else if(deviate.est$type == "val"){
    hFDR.low <- deviate.est$low
    hFDR.up <- deviate.est$up
  }
  
  hFDR.res <- structure(list(call = match.call(),
                             X = X,
                             y = y,
                             glmnet.cv = lasso.cv,
                             lambda = lambda,
                             hFDR = hFDR,
                             hFDR.std = deviate.est,
                             hFDR.low = hFDR.low,
                             hFDR.up = hFDR.up,
                             name = "estimated FDR"),
                        class = 'lasso.cv.hFDR')
  
  return(hFDR.res)
}


pack_glmnet <- function(X, y, lasso.fit, tol = 1e-7){
  n <- NROW(X)
  p <- NCOL(X)
  
  Sigma <- t(X) %*% X
  Xy <- c(t(y) %*% X)
  
  lambdas <- lasso.fit$lambda * n
  nlambda <- length(lambdas)
  
  selected.list <- list()
  Sigma.selected.inv.list <- list()
  pre_selected <- NA
  for(lambda_i in 1:nlambda){
    selected <- sort(which(abs(lasso.fit$beta[, lambda_i]) > tol))
    
    if(!identical(pre_selected, selected)){
      Sigma.selected.inv <- if(length(selected) > 0){
        list(solve(Sigma[selected, selected]))
      } else{ list(NULL) }
    }
    
    Sigma.selected.inv.list[lambda_i] <- Sigma.selected.inv
    selected.list[[lambda_i]] <- selected
    
    pre_selected <- selected
  }
  
  return(list(Sigma = Sigma, Xy = Xy,
              lambdas = lambdas, betas = lasso.fit$beta,
              selected.list = selected.list,
              Sigma.selected.inv.list = Sigma.selected.inv.list,
              tol = tol))
}


calc_lasso_hFDR <- function(glmnet.pack, data.pack, guess_null_by, couple = T, n_cores = 1){
  X <- data.pack$X
  y <- data.pack$y
  tol <- glmnet.pack$tol

  n <- NROW(X)
  p <- NCOL(X)

  lambdas <- glmnet.pack$lambdas
  nlambda <- length(lambdas)

  # DRj_lambda <- matrix(NA, nrow = p, ncol = nlambda)

  null.guessed.res <- guess_null.lasso(data.pack, method = guess_null_by, tune_seq = lambdas)
  null.guessed <- null.guessed.res$nulls
  unselect.prob <- null.guessed.res$probs

  parallel <- process_parallel(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector

  debug_df <- data.frame()
  
  DRj_lambda <- forall(j = 1:p, .options.multicore = list(preschedule = F)) %exec% {
    DRj <- rep(NA, nlambda)
    for(lambda_i in 1:nlambda){

      selected <- glmnet.pack$selected.list[[lambda_i]]
      lasso.homopath <- NULL

      lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
                         beta_lasso = glmnet.pack$betas[, lambda_i],
                         selected = selected,
                         Sigma = glmnet.pack$Sigma,
                         Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
                         lambda = lambdas[lambda_i],
                         tol = tol)
      Xjy_range <- c(data.pack$Xy_bound[, j])
      
      
      # y_Sj <- y - data.pack$vjy_obs[j] * data.pack$vj_mat[, j]
      # lasso_Sj <- glmnet::glmnet(X, y_Sj, lambda = lambdas[lambda_i] / n,
      #                               intercept = F, standardize = F,
      #                               standardize.response = F, family = "gaussian")
      # Rj <- max(1, sum(abs(lasso_Sj$beta) > 1e-10))
      # threshold <- 1 / (p / Rj + 1)
      # null.guessed.res <- guess_null(data.pack, threshold = threshold)
      # null.guessed <- null.guessed.res$nulls
      # unselect.prob <- null.guessed.res$probs
      
      if(!null.guessed[j, lambda_i]){
        # DRj_lambda[j, lambda_i] <- 0
        DRj[lambda_i] <- 0
      } else{
        if(couple == T){
          # start.time <- Sys.time()
          
          lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
          res <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack, tol = tol)
          
          # end.time <- Sys.time()
          # time.taken <- as.numeric(difftime(end.time, start.time, units = "secs"))

          DRj[lambda_i] <- res$DRj
          prob_j_unselect <- 1 - res$Evj
          
          # DRj[lambda_i] <- if(DRj[lambda_i] > 0) length(lasso.homopath$Xjy_nodes) else -length(lasso.homopath$Xjy_nodes)

        } else if(couple == F){
          prob_j_unselect <- calc_prob_unselect(lasso.pack, j, data.pack)
          DRj[lambda_i] <- (1 - prob_j_unselect) / max(length(selected), 1)
        }
      }
      
      DRj[lambda_i] <- DRj[lambda_i] / unselect.prob[j, lambda_i]
      
      
      
      # if(lambda_i == 0){
      #   tvals <- data.pack$vjy_obs / sqrt(data.pack$RSS_X / (n-p))
      #   pvals <- pvals_t(tvals, n-p, side = "two")
      #   if(!(j %in% null.guessed) || couple == F){
      #     lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
      #     res <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack, tol = tol)
      #   }
      #   Xj_proj <- X[, j] - data.pack$vj_mat[, j] * sum(data.pack$vj_mat[, j] * X[, j])
      #   Xj_proj_normed <- Xj_proj / sqrt(sum(Xj_proj^2))
      # 
      #   debug_df <- rbind(debug_df,
      #                     data.frame(j = j, DRj = res$DRj, prob_select = res$Evj, FDR_j = DRj[lambda_i],
      #                                threshold = threshold, Rj = Rj,
      #                                pval = pvals[j], guess_null = (pvals[j] >= threshold),
      #                                lasso_selected = (j %in% selected),
      #                                vjy = data.pack$vjy_obs[j], Xjy_Sj = sum(y * Xj_proj_normed),
      #                                Xj_R2 = sum(Xj_proj^2) / sum(X[, j]^2)))
      # }
      
    }
    return(DRj)
  }

  DRj_lambda <- do.call(rbind, DRj_lambda)
  
  # sigma <- sqrt(sum(lm(y ~ X + 1)$residuals^2) / (n-p-1))
  # hbeta <- solve(t(X)%*%X)%*%t(X)%*%y
  # std <- sqrt(diag(solve(t(X)%*%X))) * sigma
  # tvals <- lm_to_t(y, X)$tvals
  # pvals <- pvals_t(tvals, n-p-1)
  # 
  # data <- data.frame(beta=BETA, hbeta = hbeta, hbeta_std = std, tvals = tvals, pvals = pvals, sel_prob = DRj_lambda)
  # data <- round(data, digits = 4)
  # 
  # browser()
  
  # zero_contrib_nodes=DRj_lambda
  # zero_contrib_nodes[zero_contrib_nodes >= 0]=0
  # zero_contrib_nodes=abs(colSums(zero_contrib_nodes))
  # pos_contrib_nodes=DRj_lambda
  # pos_contrib_nodes[pos_contrib_nodes >= 0]=0
  # pos_contrib_nodes=DRj_lambda
  # pos_contrib_nodes[pos_contrib_nodes <= 0]=0
  # pos_contrib_nodes=abs(colSums(pos_contrib_nodes))
  # n_0_contrib=DRj_lambda
  # n_0_contrib[n_0_contrib>=0]=0
  # n_0_contrib[n_0_contrib<0]=1
  # n_0_contrib=colSums(n_0_contrib)
  # n_kill=colSums(DRj_lambda==0)
  # data = data.frame(lambda=lambdas/n, n_0_contrib = n_0_contrib, n_kill=n_kill, pos_contrib_nodes = pos_contrib_nodes, zero_contrib_nodes = zero_contrib_nodes, avg_zero_contrib_nodes = zero_contrib_nodes/n_0_contrib, avg_pos_contrib_nodes=pos_contrib_nodes/(p-n_0_contrib-n_kill))
  # 
  # browser()

  # DRj_lambda <- DRj_lambda / unselect.prob

  hFDR_lambda <- colSums(DRj_lambda)

  return(list(hFDR_lambda = hFDR_lambda, DRj_lambda = DRj_lambda))
}

# est_lasso_hFDR_std <- function(glmnet.pack, data.pack, guess_null_by,
#                                method = "Bootstrap.lasso_ols", measure = "std",
#                                couple = F, n_cores = 1){
#   X <- data.pack$X
#   y <- data.pack$y
#   tol <- glmnet.pack$tol
#   
#   n <- NROW(X)
#   p <- NCOL(X)
#   
#   lambdas <- glmnet.pack$lambdas
#   nlambda <- length(lambdas)
#   
#   
#   if(method == "Poincare"){
#     
#     indicator_nabla <- matrix(0, nrow = n, ncol = nlambda)
#     R_nabla_Xy <- matrix(0, nrow = p, ncol = nlambda)
#     indicators <- rep(0, nlambda)
#     
#     null.guessed <- guess_null.lasso(data.pack, method = guess_null_by, tune_seq = lambdas)$nulls
#     
#     for(j in 1:p){
#       for(lambda_i in 1:nlambda){
#         if(null.guessed[j, lambda_i]){
#           selected <- glmnet.pack$selected.list[[lambda_i]]
#           lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
#                              beta_lasso = glmnet.pack$betas[, lambda_i],
#                              selected = selected,
#                              Sigma = glmnet.pack$Sigma,
#                              Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
#                              lambda = lambdas[lambda_i],
#                              tol = tol)
#           Xjy_range <- c(data.pack$Xy_bound[, j])
#           
#           lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
#           
#           nablas <- calc_sensitivity(lasso.homopath, j, data.pack)
#           indicator_nabla[, lambda_i] <- indicator_nabla[, lambda_i] + nablas$indicator_nabla
#           R_nabla_Xy[j, lambda_i] <- nablas$R_nabla_Xjy
#           indicators[lambda_i] <- indicators[lambda_i] + res$Evj
#         }
#       }
#     }
#     
#     R_nabla <- t(X %*% R_nabla_Xy)
#     indicator_nabla <- t(indicator_nabla)
#     R <- sapply(1:nlambda, function(lambda_i){
#       max(length(glmnet.pack$selected.list[[lambda_i]]), 1)
#     })
#     hFDR_nabla <- (indicator_nabla * R - R_nabla * indicators) / R^2 * 2
#     
#     std.est <- sqrt(rowSums(hFDR_nabla^2)) * sigma_hat
#     deviate.est <- list(low = -std.est, up = std.est, type = "add")
#   }
#   
#   if(method == "Bootstrap.Bayes"){
#     lm.Bayes <- suppressWarnings(stan_lm(y ~ X + 0, prior = R2(0.5), refresh = 0)) # R2(0.5)
#     Bayes.samples <- t(rstan::extract(lm.Bayes$stanfit)$beta[1:10, 1, ])
#     
#     n_samples <- NCOL(Bayes.samples)
#     
#     var.est.res <- sapply(1:n_samples, function(iter){
#       beta_hat <- Bayes.samples[, iter]
#       
#       X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
#                      level_score = data.pack$level_score)
#       
#       hFDR.samples <- matrix(NA, nlambda, mc_size)
#       
#       for(boot_i in 1:mc_size){
#         X.boot <- X
#         y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
#         
#         data.pack.boot <- process_data(X.boot, y.boot)
#         
#         lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
#                                          intercept = T, standardize = F,
#                                          standardize.response = F, family = "gaussian")
#         glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
#         
#         hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, guess_null_by, couple = F, n_cores = n_cores)$hFDR_lambda
#         
#         hFDR.samples[, boot_i] <- hFDR.boot
#       }
#       
#       return(rowSums((hFDR.samples - rowMeans(hFDR.samples))^2) / (mc_size-1))
#     })
#     
#     std.est <- sqrt(rowMeans(var.est.res))
#     deviate.est <- list(low = -std.est, up = std.est, type = "add")
#   }
#   
#   if(method %in% c("Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara", "Efron_Stein")) {
#     mc_size <- 20
#     
#     if(method == "Bootstrap.ols"){
#       beta_hat <- data.pack$vjy_obs * sqrt(data.pack$level_score)
#     }
#     if(method %in% c("Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se")){
#       lasso.cv <- cv.glmnet(X, y, nfolds = 5, intercept = F, standardize = F,
#                             standardize.response = F, family = "gaussian")
#       beta_hat.lasso <- lasso.cv$glmnet.fit$beta[, lasso.cv$index[1]]
#       
#       beta_hat <- rep(0, p)
#       if(any(beta_hat.lasso != 0)){
#         beta_hat[beta_hat.lasso != 0] <- lm(y ~ X[, beta_hat.lasso != 0] + 0)$coefficients
#       }
#     }
#     
#     X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
#                    level_score = data.pack$level_score)
#     
#     sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
#     
#     hFDR.samples <- matrix(NA, nlambda, mc_size)
#     FDP.samples <- matrix(NA, nlambda, mc_size)
#     
#     # pvals.samples <- matrix(NA, p, mc_size)
#     # lasso_selected.samples <- matrix(NA, p, mc_size)
#     
#     for(boot_i in 1:mc_size){
#       
#       if(method %in% c("Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se")){
#         X.boot <- X
#         y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
#       }
#       if(method == "Efron_Stein"){
#         leave_out_num <- n / mc_size
#         leave_out_index <- (((boot_i-1)*leave_out_num+1) : (boot_i*leave_out_num))
#         X.boot <- X[-leave_out_index, ]
#         X.boot <- scale_X(X.boot, F)
#         y.boot <- y[-leave_out_index]
#       }
#       if(method == "Bootstrap.nonpara"){
#         bs.samples <- sample(1:n, n, replace = T)
#         X.boot <- X[bs.samples, ]
#         y.boot <- y[bs.samples]
#         X.boot <- scale_X(X.boot, F)
#       }
#       data.pack.boot <- process_data(X.boot, y.boot)
#       
#       lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
#                                        intercept = F, standardize = F,
#                                        standardize.response = F, family = "gaussian")
#       glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
#       
#       hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, guess_null_by, couple = T, n_cores = n_cores)$hFDR_lambda
#       
#       hFDR.samples[, boot_i] <- hFDR.boot
#       FDP.samples[, boot_i] <- calc_lasso_FDP(lasso.fit.boot, beta_hat == 0)
#       
#       # pvals.samples[, boot_i] <- pvals_t(lm_to_t(y.boot, X.boot)$tvals, n-p)
#       # lambda_i = 1
#       # lasso_selected.samples[, boot_i] <- (abs(lasso.fit.boot$beta[, lambda_i]) > 1e-8)
#     }
#     FDR.samples <- rowMeans(FDP.samples)
#     
#     # pvals <- pvals_t(lm_to_t(y, X)$tvals, n-p)
#     # lambda_i = 1
#     # lasso.fit <- glmnet::glmnet(X, y, lambda = lambdas / n,
#     #                                  intercept = F, standardize = F,
#     #                                  standardize.response = F, family = "gaussian")
#     # lasso_selected <- (abs(lasso.fit$beta[, lambda_i]) > 1e-8)
#     # beta_sigma <- sqrt(diag(solve(t(X)%*%X))) * sigma_hat
#     # data <- data.frame(pvals = pvals, lasso_selected = lasso_selected, Xjy = abs(t(X) %*% y),
#     #                    beta_hat = beta_hat, beta_sigma = beta_sigma,
#     #                    beta_hat_pval = 2 * pnorm(abs(beta_hat), sd = beta_sigma, lower.tail = FALSE),
#     #                    pvals_left = 2 * pnorm(abs(beta_hat - 2*beta_sigma), sd = beta_sigma, lower.tail = FALSE),
#     #                    pvals_right = 2 * pnorm(abs(beta_hat + 2*beta_sigma), sd = beta_sigma, lower.tail = FALSE))
#     # row.names(data) <- 1:p
#     # threshold <- 1 / (p / (sum(data$lasso_selected)) + 1)
#     # browser()
#     
#     if(method %in% c("Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara")){
#       deviate.est <- calc_deviate(hFDR.samples, FDR.samples, measure = measure)
#     } else if(method %in% c("Efron_Stein")){
#       hFDR_lambda <- calc_lasso_hFDR(glmnet.pack, data.pack, guess_null_by, couple = T, n_cores = n_cores)$hFDR_lambda
#       std.est <- sqrt(rowSums((hFDR_lambda - hFDR.samples)^2))
#       deviate.est <- list(low = -std.est, up = std.est, type = "add")
#     }
#   }
#   
#   return(list(deviate.est = deviate.est))
# }




# bias_var_iid <- function(X, beta, sigma, lambda, eta){
#   H0 <- beta == 0
#   X_norm <- sqrt(colSums(X^2))
#   Xy_mean <- X_norm^2 * beta
#   Xy_pval_thres <- abs(qnorm(eta/2, 0, X_norm[i] * sigma))
#   
#   prob <- sapply(1:length(beta), function(i){
#     prob_sel <- pnorm(-lambda, Xy_mean[i], X_norm[i] * sigma) + (1 - pnorm(lambda, Xy_mean[i], X_norm[i] * sigma))
#     prob_pval <- pnorm(Xy_pval_thres, Xy_mean[i], X_norm[i] * sigma) - pnorm(-Xy_pval_thres, Xy_mean[i], X_norm[i] * sigma)
#     prob_join <- if(Xy_pval_thres > lambda){
#       prob_sel - (1-prob_pval)
#     } else{
#       0
#     }
#     return(c(prob_sel, prob_pval, prob_join))
#   })
#   prob_sel <- prob[1, ]
#   prob_pval <- prob[2, ]
#   prob_join <- prob[3, ]
#   
#   ER <- max(1, sum(prob_sel))
#   bias <- (prob_join * (!H0)) / ER / (1 - eta)
# }

guess_null.lasso <- function(data.pack, method = "p-value", tune_seq){
  
  n <- data.pack$n
  p <- data.pack$p
  nlambda <- length(tune_seq)
  
  X <- data.pack$X
  y <- data.pack$y
  
  if(method == "p-value"){
    thresholds <- matrix(NA, nrow = p, ncol = nlambda)
    nulls <- matrix(NA, nrow = p, ncol = nlambda)
    probs <- matrix(NA, nrow = p, ncol = nlambda)
    
    tvals <- data.pack$vjy_obs / sqrt(data.pack$RSS_X / (n-p))
    pvals <- pvals_t(tvals, n-p, side = "two")
    
    for(j in 1:p){
      for(lambda_i in 1:nlambda){
        # y_Sj <- y - data.pack$vjy_obs[j] * data.pack$vj_mat[, j]
        # lasso_Sj <- glmnet::glmnet(X, y_Sj, lambda = tune_seq[lambda_i] / n,
        #                            intercept = F, standardize = F,
        #                            standardize.response = F, family = "gaussian")
        # Rj <- max(1, sum(abs(lasso_Sj$beta) > 1e-10))
        # threshold <- 1 / (p / Rj + 1)
        threshold <- 0.1
        thresholds[j, lambda_i] <- threshold
        nulls[j, lambda_i] <- (pvals[j] > threshold)
        probs[j, lambda_i] <- 1 - threshold
      }
    }
  } else if(method == "selection"){
    sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
    lambda0 <- sigma_hat / n
    lasso.fit <- glmnet::glmnet(X, y, lambda = lambda0,
                                intercept = F, standardize = F,
                                standardize.response = F, family = "gaussian")
    glmnet.pack <- pack_glmnet(X, y, lasso.fit)
    
    selected <- glmnet.pack$selected.list[[1]]
    # nulls <- which(!((1:p) %in% selected))
    # probs <- rep(1, p)
    
    lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
                       beta_lasso = glmnet.pack$betas[, 1],
                       selected = selected,
                       Sigma = glmnet.pack$Sigma,
                       Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[1]],
                       lambda = lambda0,
                       tol = glmnet.pack$tol)
    # for(j in nulls){
    #   probs[j] <- calc_prob_unselect(lasso.pack, j, data.pack)
    # }
    thresholds <- matrix(NA, nrow = p, ncol = nlambda)
    nulls <- matrix(NA, nrow = p, ncol = nlambda)
    probs <- matrix(NA, nrow = p, ncol = nlambda)
    for(j in 1:p){
      nulls[j, ] <- !(j %in% selected)
      probs[j, ] <- calc_prob_unselect(lasso.pack, j, data.pack)
    }
  }
  
  return(list(nulls = nulls, probs = probs, thresholds = thresholds))
}

calc_DRj <- function(lasso.homopath, j, Xjy_range, data.pack, tol = 1e-7){
  res_norm2 <- data.pack$RSS_Xnoj[j]
  df <- data.pack$n - data.pack$p
  trans <- data.pack$trans
  
  n_nodes <- length(lasso.homopath$Xjy_nodes)
  mid_beta <- (lasso.homopath$beta_at_nodes[, 1:(n_nodes-1)] + 
    lasso.homopath$beta_at_nodes[, 2:n_nodes]) / 2
  mid_beta <- matrix(mid_beta, ncol = n_nodes-1)
  DPj <- (abs(mid_beta[j, ]) > tol) / pmax(colSums(abs(mid_beta) > tol), 1)
  
  trunc.low <- sum(lasso.homopath$Xjy_nodes <= min(Xjy_range))
  trunc.up <- sum(lasso.homopath$Xjy_nodes >= max(Xjy_range))
  if(trunc.low < 1 || trunc.up < 1) browser()
  main_trunk <- if(trunc.low+1 <= n_nodes-trunc.up){
    lasso.homopath$Xjy_nodes[(trunc.low+1):(n_nodes-trunc.up)]
  } else { NULL }
  # int_nodes <- c(min(Xjy_range), main_trunk, max(Xjy_range))
  int_nodes <- main_trunk
  int_nodes <- Xjy_to_vjy(int_nodes, trans, j)
  DPj <- DPj[trunc.low:(n_nodes-trunc.up)]
  
  CDF <- c(trans$Xv[j] < 0, vjy_CDF(int_nodes, res_norm2, df), trans$Xv[j] > 0)
  DRj <- sum(abs(diff(CDF)) * DPj)
  Evj <- sum(abs(diff(CDF)) * (DPj > 0))
  
  return(list(DRj = DRj, Evj = Evj))
}

calc_prob_unselect <- function(lasso.pack, j, data.pack){
  p <- data.pack$p
  n <- data.pack$n
  df <- n - p
  
  Xjy_range <- c(data.pack$Xy_bound[, j])
  
  if(!(j %in% lasso.pack$selected)){
    Xjy_nodes <- c()
    for(direction in c(1, -1)){
      res_at_node <- lasso_homotopy_step(lasso.pack, j, Xjy_range, direction)
      Xjy_nodes <- c(Xjy_nodes, res_at_node$Xjy)
    }
    Xjy_nodes[Xjy_nodes <= min(Xjy_range)] <- min(Xjy_range)
    Xjy_nodes[Xjy_nodes >= max(Xjy_range)] <- max(Xjy_range)
    int_nodes <- sort(Xjy_to_vjy(Xjy_nodes, data.pack$trans, j))
    prob_unselect <- abs(diff(vjy_CDF(int_nodes, data.pack$RSS_Xnoj[j], df)))
  } else{
    lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
    Evj <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack,
                    tol = lasso.pack$tol)$Evj
    prob_unselect <- 1 - Evj
  }
  
  return(prob_unselect)
}

calc_sensitivity <- function(lasso.homopath, j, data.pack){
  p <- data.pack$p
  n <- data.pack$n
  df <- n - p
  
  Xjy_range <- c(data.pack$Xy_bound[, j])
  sen_coef <- lasso.homopath$sen_coef
  
  Xjy_nodes <- NULL
  for(i in 1:length(lasso.homopath$Xjy_nodes)){
    selected <- which(abs(lasso.homopath$beta_at_nodes[, i]) > 1e-7)
    if(!(j %in% selected)) Xjy_nodes <- c(Xjy_nodes, lasso.homopath$Xjy_nodes[i])
  }
  if(length(Xjy_nodes) != 2){
    warning("always selected.")
    return(list(indicator_nabla = 0, R_nabla_Xjy = 0))
  }
  
  permute <- 1e-4
  
  Xjy_nodes[Xjy_nodes <= min(Xjy_range) + permute] <- min(Xjy_range) + permute
  Xjy_nodes[Xjy_nodes >= max(Xjy_range) - permute] <- max(Xjy_range) - permute
  Xjy_nodes <- sort(Xjy_nodes[1:2])
  
  vjy_nodes <- Xjy_to_vjy(Xjy_nodes, data.pack$trans, j)
  vjy_permute <- Xjy_to_vjy(Xjy_nodes + permute, data.pack$trans, j)
  derivative_Xy <- (abs(diff(vjy_CDF(vjy_permute, data.pack$RSS_Xnoj[j], df))) - 
    abs(diff(vjy_CDF(vjy_nodes, data.pack$RSS_Xnoj[j], df)))) / permute
  
  derivative_y2 <- (abs(diff(vjy_CDF(vjy_nodes, data.pack$RSS_Xnoj[j] + permute, df))) - 
    abs(diff(vjy_CDF(vjy_nodes, data.pack$RSS_Xnoj[j], df)))) / permute
  
  indicator_nabla <- derivative_Xy * sen_coef + derivative_y2 * 2 * data.pack$y
  
  Xjy_nodes <- lasso.homopath$Xjy_nodes
  Xjy_nodes[Xjy_nodes >= max(Xjy_range) + 10 * abs(diff(Xjy_range))] <- max(Xjy_range) + 10 * abs(diff(Xjy_range))
  Xjy_nodes[Xjy_nodes <= min(Xjy_range) - 10 * abs(diff(Xjy_range))] <- min(Xjy_range) - 10 * abs(diff(Xjy_range))
  if(length(Xjy_nodes) <= 2){
    R_nabla_Xjy <- 0
  } else{
    Xjy <- sum(data.pack$X[, j] * data.pack$y)
    Xjy <- min(Xjy, max(Xjy_range))
    Xjy <- max(Xjy, min(Xjy_range))
    Xjy.index <- min(max(which(Xjy >= Xjy_nodes)), length(Xjy_nodes) - 1)
    
    mid.node <- (Xjy_nodes[Xjy.index] + Xjy_nodes[Xjy.index+1]) / 2
    end.node <- Xjy_nodes[Xjy.index + (Xjy >= mid.node)]
    mid.R <- length(lasso.homopath$selected_at_nodes_right[[Xjy.index]])
    end.R <- if(Xjy >= mid.node){
      length(lasso.homopath$selected_at_nodes_right[[Xjy.index+1]])
    } else{
      ifelse(Xjy.index == 1, mid.R,
             length(lasso.homopath$selected_at_nodes_right[[Xjy.index-1]]))
    }
    end.R <- (mid.R + end.R) / 2
    poly.coef <- solve(rbind(c(mid.node^2, mid.node, 1),
                             c(end.node^2, end.node, 1),
                             c(2*mid.node, 1, 0)),
                       c(mid.R, end.R, 0))
    
    R_nabla_Xjy <- 2*poly.coef[1]*Xjy + poly.coef[2]
    
    # R_nabla_Xjy <- (end.R - mid.R) / (end.node - mid.node)
    
    # R <- sapply(lasso.homopath$selected_at_nodes_right, length)
    # Xjy_nodes <- lasso.homopath$Xjy_nodes
    # R_nabla_Xjy <- lm(R ~ Xjy_nodes + 1)$coefficients[2]
    # browser()
    
    # max_try <- 10
    # mc_size <- 50
    # for(power_i in 1:max_try){
    #   width <- abs(diff(Xjy_range)) / 2 / 2^power_i
    #   cond_stats <- sapply(runif(mc_size, min(Xjy_range), max(Xjy_range)), function(Xjy.sample){
    #     calc_cond_exp_var(Xjy.sample, width, j, data.pack, lasso.homopath)
    #   })
    #   noise_ratio <- mean(cond_stats[2, ]) / max(var(cond_stats[1, ]), 1e-7)
    #   if(noise_ratio <= 1) break
    # }
    # Xjy_nodes.fit <- seq(max(min(Xjy_range), Xjy - width),
    #                      min(max(Xjy_range), Xjy + width), length.out = 100)
    # R.fit <- sapply(Xjy_nodes.fit, function(Xjy){
    #   max(length(homotopy_fit(Xjy, lasso.homopath)$selected), 1)
    # })
    # R_nabla_Xjy <- lm(R.fit ~ Xjy_nodes.fit + 1)$coefficients[2]
    R_nabla_Xjy <- 0
  }
  
  return(list(indicator_nabla = indicator_nabla, R_nabla_Xjy = R_nabla_Xjy))
}

calc_cond_exp_var <- function(midpoint, width, j, data.pack, lasso.homopath){
  p <- data.pack$p
  n <- data.pack$n
  df <- n - p
  
  Xjy_range <- sort(c(data.pack$Xy_bound[, j]))
  Xjy_range <- c(max(Xjy_range[1], midpoint - width),
                 min(Xjy_range[2], midpoint + width))
  
  Xjy_nodes <- lasso.homopath$Xjy_nodes
  nodes.in_range.index <- which((Xjy_nodes > Xjy_range[1]) & (Xjy_nodes < Xjy_range[2]))
  range_intv.left_node.index <- c(max(which(Xjy_nodes <= Xjy_range[1])), nodes.in_range.index)
  
  intv.nodes.Xjy <- c(Xjy_range[1], Xjy_nodes[nodes.in_range.index], Xjy_range[2])
  intv.nodes.vjy <- Xjy_to_vjy(intv.nodes.Xjy, data.pack$trans, j)
  intv.mass <- abs(diff(vjy_CDF(intv.nodes.vjy, data.pack$RSS_Xnoj[j], df)))
  intv.mass <- intv.mass / sum(intv.mass)
  
  intv.R <- sapply(lasso.homopath$selected_at_nodes_right, length)[range_intv.left_node.index]
  intv.R <- pmax(intv.R, 1)
  
  cond_mean <- sum(1/intv.R * intv.mass)
  cond_var <- sum(1/intv.R^2 * intv.mass) - cond_mean^2
  
  return(c(cond_mean, cond_var))
}


plot.lasso.cv.hFDR <- function(hFDR.obj, H0, lasso.FDR, lasso.TPR, sign.lambda = -1,
                               show_cv = T, show_hFDR = T, show_TPR = F, ...){
  xlab <- expression(Log(lambda))
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  
  nlambda <- length(hFDR.obj$glmnet.cv$lambda)
  lasso.FDP <- sapply(1:nlambda, function(lambda_index){
    lasso.selected <- which(lasso.cv$glmnet.fit$beta[, lambda_index] != 0)
    false_discover <- length(intersect(lasso.selected, which(H0)))
    false_discover / max(length(lasso.selected), 1)
  })
  lasso.TPP <- sapply(1:nlambda, function(lambda_index){
    lasso.selected <- which(lasso.cv$glmnet.fit$beta[, lambda_index] != 0)
    true_discover <- length(intersect(lasso.selected, which(!H0)))
    true_discover / max(sum(!H0), 1)
  })
  
  
  plot.range <- range(hFDR.obj$hFDR, hFDR.obj$hFDR.low, hFDR.obj$hFDR.up)
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  # plot.range <- range(hFDR.obj$hFDR)
  plot.args <- list(x = sign.lambda * log(hFDR.obj$lambda),
                    y = hFDR.obj$hFDR,
                    xlim = range(sign.lambda*log(hFDR.obj$glmnet.cv$lambda)),
                    ylim = plot.range,
                    xlab = xlab, ylab = "values", type = "n")
  
  # "estimated FDR and scaled CV MSE"
  # new.args <- list(...)
  # if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  cv.range <- abs(diff(range(c(hFDR.obj$glmnet.cv$cvm, hFDR.obj$glmnet.cv$cvlo, hFDR.obj$glmnet.cv$cvup))))
  
  leg <- c()
  lty <- c()
  lwd <- c()
  pch <- c()
  col <- c()
  
  do.call("plot", plot.args)
  
  if(show_cv){
    lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
          y = (hFDR.obj$glmnet.cv$cvm - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
          col = "dodgerblue3")
    error.bars(sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
               (hFDR.obj$glmnet.cv$cvup - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
               (hFDR.obj$glmnet.cv$cvlo - min(hFDR.obj$glmnet.cv$cvlo)) / cv.range * abs(diff(plot.range)) + min(plot.range),
               width = 0.01, col = alpha("dodgerblue3", 0.1))
    
    abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.min), lty = 3)
    abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.1se), lty = 3)
    
    leg <- c(leg, "CV MSE")
    lty <- c(lty, 1)
    lwd <- c(lwd, 1)
    pch <- c(pch, 26)
    col <- c(col, "dodgerblue3")
  }
  
  if(show_TPR){
    lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
          y = 1-lasso.TPP,
          col = "orange3", lty = 2)
    lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
          y = 1-lasso.TPR,
          col = "orange3")
    
    leg <- c(leg, "Type II error", "Type II error rate")
    lty <- c(lty, 2, 1)
    lwd <- c(lwd, 1, 1)
    pch <- c(pch, 26, 26)
    col <- c(col, "orange3", "orange3")
  }
  
  lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
        y = lasso.FDP,
        col = "black", lty = 2)
  lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
        y = lasso.FDR,
        col = "black")
  
  leg <- c(leg, "FDP", "FDR")
  lty <- c(lty, 2, 1)
  lwd <- c(lwd, 1, 1)
  pch <- c(pch, 26, 26)
  col <- c(col, "black", "black")
  
  
  if(show_hFDR){
    error.bars(sign.lambda * log(hFDR.obj$lambda),
               hFDR.obj$hFDR.up, hFDR.obj$hFDR.low,
               width = 0.01, alpha("red", 0.3))
    points(sign.lambda*log(hFDR.obj$lambda), hFDR.obj$hFDR,
           pch = 20, col = "red")
    
    leg <- c(leg, "hFDR")
    lty <- c(lty, 1)
    lwd <- c(lwd, 0)
    pch <- c(pch, 19)
    col <- c(col, "red")
  }
  
  axis(side = 3, at = sign.lambda*log(hFDR.obj$glmnet.cv$lambda),
       labels = paste(hFDR.obj$glmnet.cv$nz), tick = FALSE, line = 0)
  
  legend("bottomright", inset = 0.05,
         legend = leg,
         lty = lty, lwd = lwd,
         pch = pch, col = col)
  invisible()
}


plot_path <- function(hFDR.obj, sign.lambda = -1,...){
  xlab <- expression(Log(lambda))
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  
  nbeta <- NCOL(hFDR.obj$glmnet.cv$glmnet.fit$beta)
  plot.args <- list(x = sign.lambda * log(hFDR.obj$lambda),
                    xlim = range(sign.lambda*log(hFDR.obj$glmnet.cv$lambda)),
                    ylim = range(as.vector(hFDR.obj$glmnet.cv$glmnet.fit$beta)),
                    xlab = xlab, ylab = "coefficients", type = "n")
  new.args <- list(...)
  if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  do.call("plot", plot.args)
  for(beta_i in 1:nbeta){
    lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
          y = as.vector(hFDR.obj$glmnet.cv$glmnet.fit$beta[beta_i, ]),
          col = "grey")
  }
  invisible()
}















process_X <- function(X){
  n <- NROW(X)
  p <- NCOL(X)
  
  QR <- qr(X)
  
  pivot_back <- sort.list(QR$pivot)
  
  Q_X <- qr.Q(QR, complete = F)
  R_X <- qr.R(QR, complete = F)[, pivot_back]
  
  # compute the a basic matrix needed by cknockoff
  # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
  vj_mat <- sapply(1:p, function(j){
    # Q_X[, j] = X_j orthogonal projected onto X_{1:j-1}
    # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
    #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
    # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
    coord <- forwardsolve(t(R_X[j:p,j:p]), c(1, rep(0, p-j)))
    vj <- Q_X[, j:p] %*% matrix(coord, nrow = p-j+1)
    vj <- vj / sqrt(sum(vj^2))
  })
  
  return(list(X = X, Q_X = Q_X, vj_mat = vj_mat, level_score = diag(solve(t(X)%*%X))))
}

process_y <- function(X.pack, y, ep = 1e-4){
  X <- X.pack$X
  Q_X <- X.pack$Q_X
  vj_mat <- X.pack$vj_mat
  
  n <- NROW(X)
  p <- NCOL(X)
  
  vjy_obs <- c(matrix(y, nrow=1) %*% vj_mat)
  RSS_X <- abs(sum(y^2) - sum((matrix(y, nrow=1) %*% Q_X)^2))
  RSS_Xnoj <- RSS_X + vjy_obs^2
  
  Xv <- colSums(X * vj_mat)
  X_perpv_y <- c(t(y) %*% X) - Xv * vjy_obs
  
  trans <- list(Xv = Xv, X_perpv_y = X_perpv_y)
  
  vy_bound <- pmax(abs(vjy_quantile(ep, RSS_Xnoj, df = n-p)), abs(vjy_obs))
  Xy_bound <- sapply(1:p, function(j){
    sort(c(vjy_to_Xjy(-vy_bound[j], trans, j),
           vjy_to_Xjy(vy_bound[j], trans, j)))
  })
  
  return(list(X = X, y = y, Q_X = Q_X, vj_mat = vj_mat, level_score = X.pack$level_score,
              RSS_X = RSS_X, RSS_Xnoj = RSS_Xnoj, trans = trans,
              Xy_bound = Xy_bound, vjy_obs = vjy_obs,
              n = n, p = p))
}

process_data <- function(X, y, ep = 1e-4){
  X.pack <- process_X(X)
  data.pack <- process_y(X.pack, y, ep)
  
  return(data.pack)
}

Xjy_to_vjy <- function(Xjy, trans, j){
  (Xjy - trans$X_perpv_y[j]) / trans$Xv[j]
}

vjy_to_Xjy <- function(vjy, trans, j){
  trans$Xv[j] * vjy + trans$X_perpv_y[j]
}

# the distribution function of Xj^perp y conditional on Sj. See 11.summary Thm 3.1
vjy_CDF <- function(vjy, res_norm2, df){
  return(pt(vjy * sqrt(df) / sqrt(pmax(res_norm2 - vjy^2, 0)), df = df))
}

vjy_quantile <- function(prob, res_norm2, df){
  t_quantile <- qt(prob, df)
  if(abs(t_quantile) < Inf){
    vjy <- sign(t_quantile) * sqrt(t_quantile^2 * res_norm2 / (t_quantile^2 + df))
  } else{
    vjy <- sign(t_quantile) * sqrt(res_norm2)
  }
  return(vjy)
}




# calc_lasso_hFDR <- function(glmnet.pack, data.pack, couple = T, var.est){
#   X <- data.pack$X
#   y <- data.pack$y
#   tol <- glmnet.pack$tol
#   
#   n <- NROW(X)
#   p <- NCOL(X)
#   
#   sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
#   
#   lambdas <- glmnet.pack$lambdas
#   nlambda <- length(lambdas)
#   
#   DRj_lambda <- matrix(NA, nrow = p, ncol = nlambda)
#   DRj_lambda.org <- matrix(NA, nrow = p, ncol = nlambda)
#   
#   indicator_nabla <- matrix(0, nrow = n, ncol = nlambda)
#   R_nabla_Xy <- matrix(0, nrow = p, ncol = nlambda)
#   indicators <- rep(0, nlambda)
#   
#   null.guessed.res <- var.est$guess_null(data.pack)
#   null.guessed <- null.guessed.res$nulls
#   unselect.prob <- null.guessed.res$probs
#   
#   for(j in 1:p){
#     
#     for(lambda_i in 1:nlambda){
#       selected <- glmnet.pack$selected.list[[lambda_i]]
#       lasso.homopath <- NULL
#       
#       lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
#                          beta_lasso = glmnet.pack$betas[, lambda_i],
#                          selected = selected,
#                          Sigma = glmnet.pack$Sigma,
#                          Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
#                          lambda = lambdas[lambda_i],
#                          tol = tol)
#       Xjy_range <- c(data.pack$Xy_bound[, j])
#       
#       if(!(j %in% null.guessed)){
#         DRj_lambda[j, lambda_i] <- 0
#       }
#       else{
#         if(couple == T){
#           lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
#           res <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack, tol = tol)
#           
#           DRj_lambda[j, lambda_i] <- res$DRj
#           prob_j_unselect <- 1 - res$Evj
#           
#         } else if(couple == F){
#           prob_j_unselect <- calc_prob_unselect(lasso.pack, j, data.pack)
#           DRj_lambda[j, lambda_i] <- (1 - prob_j_unselect) / max(length(selected), 1)
#         }
#         
#         if(var.est$do & var.est$method == "Poincare"){
#           if(couple == F) lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
#           
#           nablas <- calc_sensitivity(lasso.homopath, j, data.pack)
#           indicator_nabla[, lambda_i] <- indicator_nabla[, lambda_i] + nablas$indicator_nabla
#           R_nabla_Xy[j, lambda_i] <- nablas$R_nabla_Xjy
#           indicators[lambda_i] <- indicators[lambda_i] + res$Evj
#         }
#       }
#     }
#     
#   }
#   
#   DRj_lambda <- DRj_lambda / unselect.prob
#   
#   hFDR_lambda <- colSums(DRj_lambda)
#   
#   if(var.est$do){
#     if(var.est$method == "Poincare"){
#       R_nabla <- t(X %*% R_nabla_Xy)
#       indicator_nabla <- t(indicator_nabla)
#       R <- sapply(1:nlambda, function(lambda_i){
#         max(length(glmnet.pack$selected.list[[lambda_i]]), 1)
#       })
#       hFDR_nabla <- (indicator_nabla * R - R_nabla * indicators) / R^2 * 2
#     } else if(var.est$method == "Bootstrap.Bayes"){
#       n_samples <- NCOL(var.est$Bayes.samples)
#       
#       var.est.res <- sapply(1:n_samples, function(iter){
#         beta_hat <- var.est$Bayes.samples[, iter]
#         
#         X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
#                        level_score = data.pack$level_score)
#         var.est.boot <- var.est
#         var.est.boot$do <- F
#         mc_size <- var.est$mc_size
#         
#         hFDR.samples <- matrix(NA, nlambda, mc_size)
#         
#         for(boot_i in 1:mc_size){
#           X.boot <- X
#           y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
#           
#           data.pack.boot <- process_data(X.boot, y.boot)
#           
#           lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
#                                            intercept = F, standardize = F,
#                                            standardize.response = F, family = "gaussian")
#           glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
#           
#           hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, couple = F,
#                                        var.est.boot)$hFDR_lambda
#           
#           hFDR.samples[, boot_i] <- hFDR.boot
#         }
#         
#         return(rowSums((hFDR.samples - rowMeans(hFDR.samples))^2) / (mc_size-1))
#       })
#       
#     }
#     else{
#       if(var.est$method == "Bootstrap.ols"){
#         beta_hat <- data.pack$vjy_obs * sqrt(data.pack$level_score)
#       } else if(var.est$method %in% c("Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se")){
#         beta_hat <- rep(0, p)
#         if(any(var.est$beta_hat.lasso != 0)){
#           beta_hat[var.est$beta_hat.lasso != 0] <- lm(y ~ X[, var.est$beta_hat.lasso != 0] + 0)$coefficients
#         }
#       }
#       
#       X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
#                      level_score = data.pack$level_score)
#       var.est.boot <- var.est
#       var.est.boot$do <- F
#       mc_size <- var.est$mc_size
#       
#       hFDR.samples <- matrix(NA, nlambda, mc_size)
#       FDP.samples <- matrix(NA, nlambda, mc_size)
#       
#       for(boot_i in 1:mc_size){
#         if(var.est$method == "Efron_Stein.oracle"){
#           leave_out_num <- n / mc_size
#           leave_out_index <- (((boot_i-1)*leave_out_num+1) : (boot_i*leave_out_num))
#           
#           hFDR.sample <- sapply(1:10, function(iter){
#             X.boot <- X
#             y.boot <- y
#             y.boot[leave_out_index] <- X[leave_out_index, ] %*% var.est$beta + rnorm(leave_out_num, sd = 1)
#             data.pack.boot <- process_data(X.boot, y.boot)
#             lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
#                                              intercept = F, standardize = F,
#                                              standardize.response = F, family = "gaussian")
#             glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
#             
#             hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, couple = F,
#                                          var.est.boot)$hFDR_lambda
#           })
#           hFDR.samples[, boot_i] <- rowMeans(hFDR.sample)
#           next
#         }
#         
#         if(var.est$method %in% c("Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se")){
#           X.boot <- X
#           y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
#         } else if(var.est$method == "Efron_Stein"){
#           leave_out_num <- n / mc_size
#           leave_out_index <- (((boot_i-1)*leave_out_num+1) : (boot_i*leave_out_num))
#           X.boot <- X[-leave_out_index, ]
#           X.boot <- scale_X(X.boot, F)
#           y.boot <- y[-leave_out_index]
#         } else if(var.est$method == "Bootstrap.nonpara"){
#           bs.samples <- sample(1:n, n, replace = T)
#           X.boot <- X[bs.samples, ]
#           y.boot <- y[bs.samples]
#           X.boot <- scale_X(X.boot, F)
#         }
#         data.pack.boot <- process_data(X.boot, y.boot)
#         
#         lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
#                                          intercept = F, standardize = F,
#                                          standardize.response = F, family = "gaussian")
#         glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
#         
#         hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, couple = F,
#                                      var.est.boot)$hFDR_lambda
#         
#         hFDR.samples[, boot_i] <- hFDR.boot
#         FDP.samples[, boot_i] <- calc_lasso_FDP(lasso.fit.boot, beta_hat == 0)
#       }
#       FDR.samples <- rowMeans(FDP.samples)
#     }
#     
#     if(var.est$method == "Poincare"){
#       std.est <- sqrt(rowSums(hFDR_nabla^2)) * sigma_hat
#       deviate.est <- list(low = -std.est, up = std.est, type = "add")
#     } else if(var.est$method %in% c("Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara")){
#       deviate.est <- calc_deviate(hFDR.samples, FDR.samples, measure = var.est$measure)
#       # sqrt(rowSums((hFDR.samples - rowMeans(hFDR.samples))^2) / (mc_size-1))
#     } else if(var.est$method %in% c("Efron_Stein", "Efron_Stein.oracle")){
#       std.est <- sqrt(rowSums((hFDR_lambda - hFDR.samples)^2))
#       deviate.est <- list(low = -std.est, up = std.est, type = "add")
#     } else if(var.est$method %in% c("Bootstrap.Bayes")){
#       std.est <- sqrt(rowMeans(var.est.res))
#       deviate.est <- list(low = -std.est, up = std.est, type = "add")
#     }
#     
#   } else{
#     deviate.est <- NA
#   }
#   
#   # return(list(hFDR_lambda = hFDR_lambda, DRj_lambda = DRj_lambda.org,
#   #             deviate.est = deviate.est,
#   #             beta_hat.ols = beta_hat.ols, beta_hat.lasso = beta_hat.lasso,
#   #             beta_hat = beta_hat, sigma_hat = sigma_hat, pvals = pvals,
#   #             couple_diff = couple_diff, phi = 1/lambda0_prob_unselect))
#   
#   return(list(hFDR_lambda = hFDR_lambda, deviate.est = deviate.est))
# }


# lasso.FDR.select <- function(X, y, lasso.fit, alpha){
#   lambdas <- lasso.fit$lambda
#   
#   hFDR.obj <- lasso.hFDR(X, y, lasso.fit)
#   
#   lambda.hat.index <- max(which(hFDR.obj$hFDR <= alpha))
#   if(abs(lambda.hat.index) == Inf) stop("lambda range too narrow")
#   
#   lambda.cand.index <- max(which(hFDR.obj$hFDR <= 1.5 * alpha))
#   candidates <- which(lasso.fit$beta[, lambda.cand.index] != 0)
#   
#   first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
#   time_of_entry <- apply(lasso.fit$beta, 1, first_nonzero)
#   names(time_of_entry) <- NULL
#   time_of_entry[is.na(time_of_entry)] <- nlambda + 1
#   
#   mc_size <- 50
#   print(candidates)
#   preselected <- sapply(candidates, function(j){
#     print(j)
#     y_samples <- y_condj_sample(y, X, j, mc_size)
#     
#     bj.mc <- rep(NA, mc_size)
#     DPj.mc <- rep(NA, mc_size)
#     for(mc_i in 1:mc_size){
#       lasso.fit.mc <- glmnet::glmnet(X, y_samples[, mc_i], lambda = lambdas,
#                                      intercept = F, standardize = F,
#                                      standardize.response = F, family = "gaussian")
#       hFDR.obj.mc <- lasso.hFDR(X, y_samples[, mc_i], lasso.fit.mc)
#       
#       lambda.hat.mc.index <- max(which(hFDR.obj.mc$hFDR <= alpha))
#       if(abs(lambda.hat.mc.index) == Inf) stop("lambda range too narrow")
#       
#       tvals <- lm_to_t(y_samples[, mc_i], X)
#       pvals <- pvals_t(tvals$tvals, n-p, side = "two")
#       
#       bj.mc[mc_i] <- (pvals[j] > 0.5) *2 * hFDR.obj$hFDRj[j, lambda.hat.mc.index]
#       
#       time_of_entry.mc <- apply(lasso.fit.mc$beta, 1, first_nonzero)
#       names(time_of_entry.mc) <- NULL
#       time_of_entry.mc[is.na(time_of_entry.mc)] <- nlambda + 1
#       
#       DPj.mc[mc_i] <- (time_of_entry.mc[j] <= time_of_entry[j]) / 
#         length(union(which(lasso.fit.mc$beta[, lambda.hat.mc.index] != 0), j))
#     }
#     # browser()
#     
#     return(mean(DPj.mc) <= mean(bj.mc))
#   })
#   preselected <- candidates[preselected]
#   
#   R.hat <- which(lasso.fit$beta[, lambda.hat.index] != 0)
#   Rj.hat <- sapply(preselected, function(j){
#     length(union(R.hat, j))
#   })
#   rand.prune <- (length(preselected) > 0) && any(length(preselected) < Rj.hat)
#   if(rand.prune){
#     u <- runif(length(preselected))
#     R <- max(which(sapply(1:length(preselected), function(r){
#       sum(u <= r / Rj.hat) >= r
#     })))
#     selected <- preselected[u <= R / Rj.hat]
#   } else{
#     selected <- preselected
#   }
#   
#   return(list(naive.selected = R.hat, pre.selected = preselected, 
#               rand.prune = rand.prune, selected = selected))
# }

