
FS_hFDR.fixedX <- function(X, y, n_sels, guess_null_by = "p-value", decompose = FALSE){
  guess.res <- guess_null(X, y, n_sels, Xcov.true = NA, method = "FS",
                          model = "fixed-X")
  hNulls <- guess.res$nulls
  filter_prob_unselect <- guess.res$probs
  
  data.pack <- process_FS_data(X, y)
  
  X <- data.pack$X
  y <- data.pack$y
  

  n <- data.pack$n
  p <- data.pack$p
  max_step <- max(n_sels)

  candidates <- 1:p
  thresholds_up <- matrix(NA, nrow = p, ncol = max_step)
  thresholds_low <- matrix(NA, nrow = p, ncol = max_step)
  
  for(step in 1:max_step){
    res <- move_forward(X, y, candidates, move = step < max_step)
    selected <- candidates[res$cand_sel]
    for(cand in 1:length(candidates)){
      j <- candidates[cand]
      vjy_obs <- data.pack$vjy_obs[j]
      vj <- data.pack$vj_mat[, j]
      if(j != selected){
        y_proj_j <- res$y_proj[cand]
        vj_to_Xj <- sum(vj * X[, j])
        thresh_1 <- (res$y_proj[res$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
        thresh_2 <- (-res$y_proj[res$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
        thresholds_up[j, step] <- max(thresh_1, thresh_2)
        thresholds_low[j, step] <- min(thresh_1, thresh_2)
        
        # mc_size=10000
        # H = X[, -j]%*%solve(t(X[, -j])%*%X[, -j])%*%t(X[, -j])
        # Xjj=H %*% X[, j]
        # Xjj = Xjj / sqrt(sum(Xjj^2))
        # Ho = diag(n) - H
        # vj = Ho %*% X[, j]
        # vj = vj / sqrt(sum(vj^2))
        # y.samples <- sampler.fixedX(X, y, j, mc_size)
        # # print(t(Xjj) %*% y.samples)
        # # print(mean((t(X[, j]) %*% y.samples)))
        # print(sum(abs(t(X[, j]) %*% y.samples) >= res$y_proj[res$cand_sel])/mc_size)
        # print(sum((t(vj) %*% y.samples) >= thresh_1)/mc_size)
        # print(vjy_CDF(max(thresholds_low[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p) +
        #         1 - vjy_CDF(min(thresholds_up[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p))
        # browser()
      } else if(sum(hNulls[j, ]) > 0 && max(n_sels[hNulls[j, ]]) >= step){
        X_v <- X
        y_v <- y
        candidates_v <- candidates[-cand]
        for(step_v in step:min(max_step, max(n_sels[hNulls[j, ]]))){
          if(step_v == p){
            thresholds_up[j, step_v] <- 0
            thresholds_low[j, step_v] <- 0
            break
          }
          res_v <- move_forward(X_v, y_v, candidates_v, aux = j, move = step_v < max_step)
          y_proj_j <- sum(X_v[, j] * y_v)
          vj_to_Xj <- sum(vj * X_v[, j])
          thresh_1 <- (res_v$y_proj[res_v$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
          thresh_2 <- (-res_v$y_proj[res_v$cand_sel] - (y_proj_j - vjy_obs*vj_to_Xj)) / vj_to_Xj
          thresholds_up[j, step_v] <- max(thresh_1, thresh_2)
          thresholds_low[j, step_v] <- min(thresh_1, thresh_2)
          
          X_v <- res_v$X
          y_v <- res_v$y
          candidates_v <- res_v$candidates
        }
      }
    }
    X <- res$X
    y <- res$y
    candidates <- res$candidates
  }
  
  DRj_tunes <- matrix(0, nrow = p, ncol = length(n_sels))
  for(j in 1:p){
    for(tune_i in 1:length(n_sels)){
      step <- n_sels[tune_i]
      if(hNulls[j, tune_i]){
        select_prob_low <- vjy_CDF(max(thresholds_low[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p)
        select_prob_up <- 1 - vjy_CDF(min(thresholds_up[j, 1:step]), res_norm2 = data.pack$RSS_Xnoj[j], df = n-p)
        DRj_tunes[j, tune_i] <- min(1, select_prob_low + select_prob_up) / step / filter_prob_unselect[j, tune_i]
      }
    }
  }
  if(decompose == T) {
    hFDR_tunes <- DRj_tunes
  } else {
    hFDR_tunes <- colSums(DRj_tunes)
  }
  return(hFDR_tunes)
}



FS_hFDR.modelX <- function(X, y, Xcov.true, n_sels, guess_null_by = "p-value", mc_size = 50){
  guess.res <- guess_null(X, y, n_sels, Xcov.true, method = "FS",
                          model = "model-X")
  hNulls <- guess.res$nulls
  filter_prob_unselect <- guess.res$probs
  
  n <- NROW(X)
  p <- NCOL(X)
  max_step <- max(n_sels)
  
  candidates <- 1:p
  X.samples <- lapply(1:p, function(j){
    if(sum(hNulls[j, ]) > 0){
      Xj.samples <- sampler.modelX(X, Xcov.true, j, mc_size)
      
      Xj.samples <- scale(Xj.samples, center = T, scale = F)
      Xj.samples <- scale(Xj.samples, center = F, scale = sqrt(colSums(Xj.samples^2)))
    } else{
      NA
    }
  })
  
  X <- scale(X, center = T, scale = F)
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  
  nsamples_selected <- matrix(0, nrow = p, ncol = max_step)
  
  for(step in 1:max_step){
    res <- move_forward(X, y, candidates, move = step < max_step)
    selected <- candidates[res$cand_sel]
    for(cand in 1:length(candidates)){
      j <- candidates[cand]
      if(sum(hNulls[j, ]) > 0 && max(n_sels[hNulls[j, ]]) >= step){
        if(j != selected){
          threshold <- abs(res$y_proj[res$cand_sel])
          
          samples_selected <- (abs(matrix(y, nrow = 1) %*% X.samples[[j]]) >= threshold)
          nsamples_selected[j, step] <- sum(samples_selected)
          if(step < max(n_sels[hNulls[j, ]])){
            X.samples[[j]] <- X.samples[[j]][, !samples_selected]
            
            X.samples[[j]] <- X.samples[[j]] - (t(X[, selected]) %*% X.samples[[j]]) %x% X[, selected]
            X.samples[[j]] <- scale(X.samples[[j]], center = F, scale = sqrt(colSums(X.samples[[j]]^2)))
          }
        } else{
          X_v <- X
          y_v <- y
          candidates_v <- candidates[-cand]
          for(step_v in step:min(max_step, max(n_sels[hNulls[j, ]]))){
            if(step_v == p){
              nsamples_selected[j, step_v] <- NCOL(X.samples[[j]])
              break
            }
            res_v <- move_forward(X_v, y_v, candidates_v, move = step_v < max_step)
            threshold <- abs(res_v$y_proj[res_v$cand_sel])
            
            samples_selected <- (abs(matrix(y_v, nrow = 1) %*% X.samples[[j]]) >= threshold)
            nsamples_selected[j, step_v] <- sum(samples_selected)
            if(step_v < max(n_sels[hNulls[j, ]])){
              X.samples[[j]] <- X.samples[[j]][, !samples_selected]
              
              X.samples[[j]] <- X.samples[[j]] - (t(X_v[, selected]) %*% X.samples[[j]]) %x% X_v[, selected]
              X.samples[[j]] <- scale(X.samples[[j]], center = F, scale = sqrt(colSums(X.samples[[j]]^2)))
            }
            
            X_v <- res_v$X
            y_v <- res_v$y
            candidates_v <- res_v$candidates
          }
        }
      }
      
    }
    X <- res$X
    y <- res$y
    candidates <- res$candidates
  }
  
  DRj_tunes <- matrix(0, nrow = p, ncol = length(n_sels))
  for(j in 1:p){
    for(tune_i in 1:length(n_sels)){
      step <- n_sels[tune_i]
      if(hNulls[j, tune_i]){
        DRj_tunes[j, tune_i] <- sum(nsamples_selected[j, 1:step]) / mc_size / step / filter_prob_unselect[j, tune_i]
      }
    }
  }
  hFDR_tunes <- colSums(DRj_tunes)
  if(is.na(sum(hFDR_tunes))) browser()
  
  return(hFDR_tunes)
}



move_forward <- function(X, y, candidates, aux = NULL, move = T){
  y_proj <- as.vector(matrix(y, nrow=1) %*% X[, candidates])
  sel <- which.max(abs(y_proj))
  sel_j <- candidates[sel]
  
  if(move){
    candidates <- candidates[-sel]
    X_sel <- X[, sel_j]
    y <- y - sum(y * X_sel) * X_sel
    
    updates <- c(candidates, aux)
    # X[, updates] <- X[, updates] - (t(X_sel) %*% X[, updates]) %x% X_sel
    # X[, updates] <- scale(X[, updates], center = F, scale = sqrt(colSums(X[, updates]^2)))
    for(j in c(candidates, aux)){
      X[, j] <- X[, j] - sum(X[, j] * X_sel) * X_sel
      X[, j] <- X[, j] / sqrt(sum(X[, j]^2))
    }
  }
  
  return(list(X = X, y = y, candidates = candidates, cand_sel = sel, y_proj = y_proj))
}


# guess_null.FS <- function(X, y, n_sels, guess_null_by, method = "p-value",
#                           model){
#   
#   n <- data.pack$n
#   p <- data.pack$p
#   nlambda <- length(tune_seq)
#   
#   X <- data.pack$X
#   y <- data.pack$y
#   
#   if(method == "p-value"){
#     thresholds <- matrix(NA, nrow = p, ncol = nlambda)
#     nulls <- matrix(NA, nrow = p, ncol = nlambda)
#     probs <- matrix(NA, nrow = p, ncol = nlambda)
#     
#     tvals <- data.pack$vjy_obs / sqrt(data.pack$RSS_X / (n-p))
#     pvals <- pvals_t(tvals, n-p, side = "two")
#     
#     for(j in 1:p){
#       for(lambda_i in 1:nlambda){
#         y_Sj <- y - data.pack$vjy_obs[j] * data.pack$vj_mat[, j]
#         lasso_Sj <- glmnet::glmnet(X, y_Sj, lambda = tune_seq[lambda_i] / n,
#                                    intercept = F, standardize = F,
#                                    standardize.response = F, family = "gaussian")
#         Rj <- max(1, sum(abs(lasso_Sj$beta) > 1e-10))
#         threshold <- 1 / (p / Rj + 1)
#         thresholds[j, lambda_i] <- threshold
#         nulls[j, lambda_i] <- (pvals[j] > threshold)
#         probs[j, lambda_i] <- 1 - threshold
#       }
#     }
#   } else if(method == "selection"){
#     sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
#     lambda0 <- sigma_hat / n
#     lasso.fit <- glmnet::glmnet(X, y, lambda = lambda0,
#                                 intercept = F, standardize = F,
#                                 standardize.response = F, family = "gaussian")
#     glmnet.pack <- pack_glmnet(X, y, lasso.fit)
#     
#     selected <- glmnet.pack$selected.list[[1]]
#     # nulls <- which(!((1:p) %in% selected))
#     # probs <- rep(1, p)
#     
#     lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
#                        beta_lasso = glmnet.pack$betas[, 1],
#                        selected = selected,
#                        Sigma = glmnet.pack$Sigma,
#                        Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[1]],
#                        lambda = lambda0,
#                        tol = glmnet.pack$tol)
#     # for(j in nulls){
#     #   probs[j] <- calc_prob_unselect(lasso.pack, j, data.pack)
#     # }
#     thresholds <- matrix(NA, nrow = p, ncol = nlambda)
#     nulls <- matrix(NA, nrow = p, ncol = nlambda)
#     probs <- matrix(NA, nrow = p, ncol = nlambda)
#     for(j in 1:p){
#       nulls[j, ] <- !(j %in% selected)
#       probs[j, ] <- calc_prob_unselect(lasso.pack, j, data.pack)
#     }
#   }
#   
#   return(list(nulls = nulls, probs = probs, thresholds = thresholds))
# }












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

process_FS_data <- function(X, y){
  n <- NROW(X)
  p <- NCOL(X)
  
  X <- scale(X, center = T, scale = F)
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  
  QR <- qr(X)
  
  pivot_back <- sort.list(QR$pivot)
  
  Q_X <- qr.Q(QR, complete = F)
  R_X <- qr.R(QR, complete = F)[, pivot_back]
  
  RSS_X <- abs(sum(y^2) - sum((matrix(y, nrow=1) %*% Q_X)^2))
  
  # y <- as.vector(matrix(y, nrow=1) %*% Q_X)
  
  # X <- R_X
  # X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
  
  
  # Q_X <- diag(p)
  
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
  
  vjy_obs <- c(matrix(y, nrow=1) %*% vj_mat)
  
  RSS_Xnoj <- RSS_X + vjy_obs^2
  
  data.pack <- list(X = X, y = y, n = n, p = p,
                    vj_mat = vj_mat, vjy_obs = vjy_obs,
                    RSS_X = RSS_X, RSS_Xnoj = RSS_Xnoj)
  
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


# library(here)
# source(here("R", "hFDR.R"))
# source(here("R", "utils.R"))
# set.seed(1)
# 
# p <- 100
# n <- 3*p
# pi1 <- 5 / p
# mu1 <- 1
# 
# # mcc5, iid6
# X.res <- gene_X("IID_Normal", n, p, 2, model_X = T)  # "IID_Normal", "MCC_Block", "X_AR_Strong"
# X <- X.res$X
# Xcov.true <- X.res$Xcov.true
# beta <- genmu(p, pi1, mu1, "random", 1)
# H0 <- beta == 0
# y <- X %*% beta + rnorm(n)
# n_sels <- 1:10
# mc_size <- 100
# # FS_hFDR.modelX(X, y, Xcov.true, n_sels, guess_null_by = "p-value", mc_size)
# # hFDR.modelX(X, y, Xcov.true, n_sels, guess_null_by = "p-value",
# #             method = "FS", mc_size, n_cores = 14)
# 
# FS_hFDR.fixedX(X, y, n_sels)
# hFDR.fixedX(X, y, n_sels, guess_null_by = "p-value", method = "FS",
#             mc_size = 100, n_cores = 14)


