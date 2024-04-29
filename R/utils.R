library(foreach)
library(doParallel)
library(leaps)

gene_X <- function(X_type = "IID_Normal", n, p, X_seed = NA,
                   scaling = T, model_X = F,
                   signal = 1, nonnulls = NA, cov_seed = NA, mis_model = F){
  if(is.na(X_seed)){
    X_seed <- NULL
  }
  if(is.na(cov_seed)){
    cov_seed <- NULL
  }
  
  with_seed(X_seed, {
    if(X_type == "IID_Normal"){
      cov_mat <- diag(p)
      
      # X <- matrix(rnorm(n*p), n)
    } else if(X_type == "Coef_AR"){
      rho <- 0.5
      cov_mat <- solve(rho^(abs(outer(1:p, 1:p, "-"))))
    } else if(X_type == "X_AR"){
      rho <- 0.5
      cov_mat <- rho^(abs(outer(1:p, 1:p, "-")))
    } else if(X_type == "X_AR_Strong"){
      rho <- 0.8
      cov_mat <- rho^(abs(outer(1:p, 1:p, "-")))
    } else if(X_type == "Homo_Block"){
      rho <- 0.5
      block_size <- 10
      
      blockSigma <- matrix(rho, block_size, block_size)
      diag(blockSigma) <- 1
      
      cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
    } else if(X_type == "MCC"){
      cov_mat <- matrix(-1/p, p, p)
      diag(cov_mat) <- 1
    } else if(X_type == "MCC_Block"){
      block_size <- 5
      
      blockSigma <- matrix(-1/block_size, block_size, block_size)
      diag(blockSigma) <- 1
      
      cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
    } else if(substr(X_type, 1, 3) == "PEC"){
      rho <- as.numeric(substr(X_type, 5, nchar(X_type))) / 10^(nchar(X_type)-4)
      beta_cov <- matrix(rho, p, p)
      diag(beta_cov) <- 1
      cov_mat <- solve(beta_cov)
    } else if(substr(X_type, 1, 3) == "NEC"){
      rho <- as.numeric(substr(X_type, 5, nchar(X_type))) / 10^(nchar(X_type)-4)
      cov_mat <- matrix(rho, p, p)
      diag(cov_mat) <- 1
    } else if(X_type == "Inv_Sparse"){
      beta_cov <- diag(1, nrow = p, ncol = p)
      lower_tri <- lower.tri(beta_cov)
      
      with_seed(cov_seed, {
        beta_cov[lower_tri] <- 1 * signal / sqrt(n) * (1:sum(lower_tri)) %in% sample(1:sum(lower_tri), nonnulls) * (1 + rexp(sum(lower_tri), 1)) / 2
      })
      
      beta_cov <- (t(beta_cov) + beta_cov) / 2
      diag(beta_cov) <- pmax(diag(beta_cov), 1.1 * (rowSums(beta_cov) - diag(beta_cov)))
      
      cov_mat <- solve(beta_cov)
    } else if(X_type == "Sparse"){
      sparsity <- 0.05
      # X <- diag(1, nrow = n, ncol = p)
      # lower_tri <- lower.tri(X)
      # X[lower_tri] <- replicate(sum(lower_tri), rbinom(1, 1, sparsity))
      for(n_try in 1:10){
        X <- matrix(rbinom(n*p, 1, sparsity), n, p)
        if(all(eigen(t(X)%*%X)$values > 0)) break
      }
      if(n_try == 10) stop("cannot generate sparse X.")
    }
    
    if(!(X_type %in% c("Sparse"))){
      scaler <- diag(1/sqrt(diag(cov_mat)))
      cov_mat <- scaler %*% cov_mat %*% scaler
      
      R <- chol(cov_mat)
      basis <- if(model_X){
        if(mis_model == 1){
          matrix(rexp(n*p) * (rbinom(n*p, 1, 0.5)*2-1), n)
        } else{
          matrix(rnorm(n*p), n)
        }
      } else{
        qr.Q(qr(matrix(rnorm(n*p), n)))
      }
      X <- basis %*% R
      if(mis_model == 2){
        X <- X * matrix(runif(n*p, 0, 2), n, p)
        # dof <- 1
        # X <- sqrt(dof / rchisq(n, dof)) * X
        # X <- rmvl(n, 0, (cov_mat+t(cov_mat))/2)
      }
    }
  })
  
  if(scaling) X <- scale_X(X)
  if(!exists("cov_mat")) cov_mat <- NA
  if(!exists("beta_cov")) beta_cov <- NA
  
  return(list(X = X, Xcov.true = cov_mat, beta_cov = beta_cov))
}

rt.mult <- function(n, p, df, cov_mat){
  R <- chol(cov_mat)
  basis <- matrix(rnorm(n*p), n)
  X <- basis %*% R
  X <- sqrt(df / rchisq(n, df)) * X
  return(X)
}

glasso.invSigma <- function(p, nonnulls){
  sparsity <- nonnulls / (p*(p-1)/2)
  invSigma <- diag(1, nrow = p, ncol = p)
  lower_tri <- lower.tri(invSigma)
  invSigma[lower_tri] <- 0.9 * rbinom(sum(lower_tri), 1, sparsity)
  invSigma <- (t(invSigma) + invSigma) / 2
  diag(invSigma) <- pmax(diag(invSigma), 1.1 * (rowSums(invSigma) - diag(invSigma)))
  
  return(invSigma)
}

scale_X <- function(X, center = F){
  X <- scale(X, center = center, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2) / (NROW(X)-1)))
}

beta_to_sparse <- function(beta, beta.sample, Sigma, sigma, mah.tol = 0.1, m1){
  p <- length(beta)
  mah.target <- p

  noise.sample <- c(beta.sample) - beta
  mah.sample <- c(t(noise.sample) %*% Sigma %*% noise.sample / sigma^2)
  # print(mah.sample)
  
  for(try_i in 1:100){
    noise.sample <- noise.sample * sqrt(mah.target / mah.sample)
    beta.sample <- beta + noise.sample
    
    threshold <- sort(abs(beta.sample), partial = p-(m1-1))[p-(m1-1)]
    beta.sample[abs(beta.sample) < threshold] <- 0
    
    noise.sample <- beta.sample - beta
    mah.sample <- c(t(noise.sample) %*% Sigma %*% noise.sample / sigma^2)
    # print(mah.sample)
    
    if(abs(mah.sample - mah.target) / mah.target <= mah.tol) break
    if(try_i == 100) warning("not convergent")
  }

  return(list(beta.sample = beta.sample, mah.sample = mah.sample))
}

sample_beta <- function(beta, Sigma, Sigma.inv.sqrt, sigma,
                        quantile_num = 9, sample_size = 100){
  if(sample_size < quantile_num+1) stop()
  
  p <- NCOL(Sigma)
  sample_trials <- 5000
  
  noise.sample <- Sigma.inv.sqrt %*% matrix(rnorm(sample_trials * p, sd = sigma), p, sample_trials)
  mahalanobis.sample <- sapply(1:sample_trials, function(sample_i){
    t(noise.sample[, sample_i]) %*% Sigma %*% noise.sample[, sample_i] / sigma^2
  })
  mahalanobis.quantiles <- unname(quantile(mahalanobis.sample, (1:quantile_num)/(quantile_num+1)))
  
  quantile_size <- abs(diff(c(0, mahalanobis.quantiles)))
  quantile_size <- quantile_size^(log(3, base = max(quantile_size)/min(quantile_size)))
  quantile_size <- ceiling(quantile_size / sum(quantile_size) * sample_size)
  sample_exceed <- sum(quantile_size) - (sample_size - 1)
  quantile_index <- quantile_num
  while(sample_exceed > 0){
    if(quantile_size[quantile_index] > 1){
      quantile_size[quantile_index] <- quantile_size[quantile_index] - 1
      sample_exceed <- sample_exceed - 1
    }
    quantile_index <- quantile_index - 1
    if(quantile_index == 0) quantile_index <- quantile_num
  }
  
  sample_indeces <- c()
  mahalanobis.quantiles <- c(0, mahalanobis.quantiles)
  for(quantile_index in 1:quantile_num){
    candidates <- which(mahalanobis.sample > mahalanobis.quantiles[quantile_index] &
                          mahalanobis.sample <= mahalanobis.quantiles[quantile_index+1])
    sample_selected <- sample(candidates, quantile_size[quantile_index],
                              replace = quantile_size[quantile_index] > length(candidates))
    sample_indeces <- c(sample_indeces, sample_selected)
  }
  
  beta.samples <- beta + cbind(rep(0, p), noise.sample[, sample_indeces])
  mah_dis.cdf <- c(0, sapply(mahalanobis.sample[sample_indeces], function(mah_dis){
    sum(mahalanobis.sample <= mah_dis) / sample_trials
  }))
  
  return(list(beta = beta.samples, mah_dis.cdf = mah_dis.cdf,
              mah_dis = mahalanobis.sample[sample_indeces], 
              mah.quantiles = mahalanobis.quantiles))
}


makeup_vectors <- function(...){
  vars <- list(...)
  max_length <- max(sapply(vars, length))
  
  env <- parent.frame()
  for(var_i in 1:length(vars)){
    var_name <- names(vars[var_i])
    var_length <- length(vars[[var_i]])
    
    if(var_length < max_length){
      var_list <- list()
      for(i in 1:max_length) var_list <- c(var_list, vars[[var_i]])
      
      env[[var_name]] <- var_list
    } else{
      env[[var_name]] <- as.list(vars[[var_i]])
    }
  }
  
  invisible()
}

# naively sampling y conditional on Sj
y_condj_sample <- function(y, X, j, sample_size, data.pack = NULL){
  
  n <- NROW(X)
  
  if(is.null(data.pack)){
    projj_y <- lm(formula = y ~ X[, -j] + 0)$fitted.values
  } else{
    proj_y <- matrix(y, nrow = 1) %*% data.pack$Q_X
    proj_y <- c(data.pack$Q_X %*% matrix(proj_y, ncol = 1))
    projj_y <- proj_y - data.pack$vjy_obs[j] * data.pack$vj_mat[, j]
  }
  
  radius <- sqrt(sum(y^2) - sum(projj_y^2))
  
  y_sample <- matrix(rnorm(n * sample_size), nrow = n)
  
  if(is.null(data.pack)){
    y_sample_projj <-lm(formula = y_sample ~ X[, -j] + 0)$fitted.values
  } else{
    y_sample_proj <- data.pack$Q_X %*% (t(data.pack$Q_X) %*% y_sample)
    y_sample_projj <- y_sample_proj - data.pack$vj_mat[, j] %*% (t(data.pack$vj_mat[, j]) %*% y_sample)
  }
  y_sample <- y_sample - y_sample_projj
  y_sample <- scale(y_sample, center = FALSE, scale = sqrt(colSums(y_sample^2)) / radius)
  y_sample <- projj_y + y_sample
  
  if(!is.null(data.pack)){
    Xjy.sample <- t(X[, j]) %*% y_sample
    in_range <- (Xjy.sample >= min(data.pack$Xy_bound[, j])) & (Xjy.sample <= max(data.pack$Xy_bound[, j]))
    if(sum(in_range) < sample_size){
      bootstrap_samples <- sample(which(in_range), sum(!in_range))
      y_sample[, !in_range] <- y_sample[, bootstrap_samples]
    }
  }
  
  return(y_sample)
}

y_from_beta <- function(X, X_proj, beta, sigma){
  n <- NROW(X)
  
  proj_y <- X %*% beta
  
  y <- rnorm(n, sd = sigma)
  y <- proj_y + y - X_proj %*% y
  
  return(y)
}

dir_trunc <- function(x, direction){
  if(length(x) == 0) return(NULL)
  
  # x <- sapply(x, function(e){
  #   if(direction * e > 0){e} else{direction * Inf}
  # })
  for(i in 1:length(x)){
    if(direction * x[i] <= 0){
      x[i] <- direction * Inf
    }
  }
  return(x)
}

solve_mat22 <- function(mat){
  a <- mat[1, 1]
  b <- mat[1, 2]
  c <- mat[2, 1]
  d <- mat[2, 2]
  det <- a*d-b*c
  inv <- matrix(c(d, -c, -b, a), ncol = 2) / det
  return(inv)
}

# gene_lambda_seq <- function(Sigma, beta, nlambda, sigma, n, scaling = T){
#   p <- NCOL(Sigma)
# 
#   lambda_max <- max(abs(Sigma %*% beta)) / n * ifelse(scaling, 1, 1.5*n)
#   lambda_min <- min(lambda_max * 0.9, 1 * sigma / n * ifelse(scaling, 1, 1.5*sqrt(n)))
#   k <- (0:(nlambda-1)) / nlambda
#   lambda_seq <- lambda_max * (lambda_min/lambda_max)^k
#   
#   return(lambda_seq)
# }
# 
# gene_glasso_rho_seq <- function(Sigma, n, ntune = 10, lower.tail = 0.95){
#   rho_max <- max(abs((Sigma)[upper.tri(Sigma)])) * n / 4
#   
#   X.null <- gene_X(X_type = "Inv_Sparse", n, p = NCOL(Sigma), 1,
#                    scaling = F, model_X = T, signal = 0, nonnulls = 0)$X
#   Sigma.null <- t(X.null) %*% X.null
#   rho_min <- unname(quantile(abs((Sigma.null)[upper.tri(Sigma.null)]), lower.tail))
#   
#   tune_seq <- exp(seq(log(rho_max), log(rho_min), length.out = ntune)) / n
#   
#   return(tune_seq)
# }
# 
# gene_FS_tune_seq <- function(X_type, n, p, scaling_X, model_X,
#                              pi1, mu1, posit_type, ntune){
#   FDR <- sapply(1:10, function(seed){
#     X <- gene_X(X_type, n, p, seed, scaling_X, model_X)$X
#     beta <- genmu(p, pi1, mu1, posit_type, 1)
#     H0 <- beta == 0
#     y <- X %*% beta + rnorm(n)
#     FS.select <- select_variables(X, y, tune_seq = 1:min(round(5*p*pi1), p, n), "FS")
#     FDP <- calc_sel_FDP(FS.select, H0)
#     return(FDP)
#   })
#   FDR <- rowMeans(FDR)
#   
#   start <- max(1, min(which(FDR > 0.01)))
#   end <- which.min(abs(FDR - 0.75))
#   end <- max(end, round(p * 0.4))
#   distance <- end - start
#   if(distance + 1 <= ntune){
#     tune_seq <- (start + min(0, p - (start + ntune - 1))) : min(p, start + ntune - 1)
#   } else{
#     tune_seq <- (start:end)[round(seq(1, distance + 1, length.out = ntune))]
#   }
#   return(tune_seq)
# }


gene_tune_seq <- function(set, sel_method, mu1, FDR_range, n_tune, mc_size = 50){
  n <- set$n
  p <- set$p
  sigma <- 1
  
  with_seed(1,{
    if(!set$random_X){
      X <- gene_X(set$X_type, n, p, set$X_seed,
                  scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
      Sigma <- t(X) %*% X
    } else{
      Sigma <- gene_X(set$X_type, n, p, set$X_seed,
                      set$scaling_X, model_X = set$model_X,
                      signal = mu1, nonnulls = p*set$pi1, mis_model = set$X_mismodel)$Xcov.true
    }
    
    if(sel_method == "lasso"){
      X <- gene_X(set$X_type, n, p, set$X_seed,
                  scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
      beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
      y <- X %*% beta + sigma * eval(set$noise)
      # lambda_max <- max(abs(t(X) %*% (X %*% beta + sigma * rnorm(n)))) / n * ifelse(set$scaling, 1, n)
      # lambda_min <- min(lambda_max * 0.9, 1 * sigma / n * ifelse(set$scaling, 1, sqrt(n)))
      # n_lambda <- 100
      # k <- (0:(n_lambda-1)) / n_lambda
      # tune_candidates <- lambda_max * (lambda_min/lambda_max)^k
      tune_candidates <- glmnet::glmnet(X, y, intercept = T, standardize = T, family = "gaussian")$lambda
    } else if(sel_method == "FS"){
      tune_candidates <- 1:min(round(p*set$pi1/(1-max(FDR_range))), p, n)
    } else if(sel_method == "logistic"){
      X <- gene_X(set$X_type, n, p, set$X_seed,
                  scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
      beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
      y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
      tune_candidates <- glmnet::glmnet(X, y, intercept = T, standardize = T, family = "binomial")$lambda
    } else if(sel_method == "ctree"){
      X <- gene_X(set$X_type, n, p, set$X_seed,
                  scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
      beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
      y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
      res <- rpart(y ~ X, method = "class", control = rpart.control(minsplit = 5))
      tune_candidates <- unique(unname(res$variable.importance))
      # tune_max <- max(unname(res$variable.importance)*1.1)
      # tune_min <- tune_max / 100
      # tune_candidates <- exp(seq(log(tune_max), log(tune_min), length.out = 50))
    } else if(sel_method == "poisson"){
      X <- gene_X(set$X_type, n, p, set$X_seed,
                  scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
      beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
      y <- rpois(n = n, lambda = exp(X %*% beta))
      tune_candidates <- glmnet::glmnet(X, y, intercept = T, standardize = T, family = "poisson")$lambda
    } else if(sel_method == "glasso"){
      # rho_max <- max(abs((Sigma)[upper.tri(Sigma)])) * n * 50
      # X.null <- gene_X(X_type = set$X_type, n, p = NCOL(Sigma), 1,
      #                  scaling = F, model_X = T, signal = 0, nonnulls = 0, mis_model = set$X_mismodel)$X
      # Sigma.null <- t(X.null) %*% X.null
      # rho_min <- unname(quantile(abs((Sigma.null)[upper.tri(Sigma.null)]), 0.98))
      # 
      # tune_candidates <- exp(seq(log(rho_max), log(rho_min), length.out = 50)) / n
      
      X.res <- gene_X(set$X_type, n, p, set$X_seed,
                      scaling = set$scaling_X, model_X = set$model_X,
                      signal = mu1, nonnulls = p*set$pi1, cov_seed = set$cov_seed, mis_model = set$X_mismodel)
      X <- X.res$X
      rho <- mean(abs((X.res$Xcov.true)))
      
      rho_max <- rho
      rho_min <- rho
      for(i in 1:1000){
        if(sum(select_variables(X, y, rho_max, sel_method)) > 1){
          rho_max <- 2 * rho_max
        } else break
      }
      for(i in 1:1000){
        if(sum(select_variables(X, y, rho_min, sel_method)) < min(p*set$pi1/(1-max(FDR_range)), p*(p-1)/2)){
          rho_min <- rho_min / 2
        } else break
      }
      tune_candidates <- exp(seq(log(rho_max), log(rho_min), length.out = 50))
    } else{
      stop()
    }
    
    FDP <- sapply(1:mc_size, function(seed){
      X.res <- gene_X(set$X_type, n, p, seed,
                      scaling = set$scaling_X, model_X = set$model_X,
                      signal = mu1, nonnulls = p*set$pi1, cov_seed = set$cov_seed, mis_model = set$X_mismodel)
      X <- X.res$X
      if(sel_method %in% c("logistic", "ctree")){
        beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
        H0 <- beta == 0
        y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
      } else if(sel_method == "poisson"){
        beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
        H0 <- beta == 0
        y <- rpois(n = n, lambda = exp(X %*% beta))
      } else if(sel_method == "glasso"){
        invSigma <- solve(X.res$Xcov.true)
        invSigma[abs(invSigma) < 1e-8] <- 0
        H0 <- c(invSigma[upper.tri(invSigma)] == 0)
        y <- NA
      } else{
        beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
        H0 <- beta == 0
        y <- X %*% beta + eval(set$noise)
      }
      selected <- select_variables(X, y, tune_candidates, sel_method)
      FDP <- calc_sel_FDP(selected, H0)
      return(FDP)
    })
    FDR <- rowMeans(FDP)
  })
  
  
  if(sel_method != "FS"){
    tune_candidates <- log(tune_candidates)
  }
  # print(FDR)
  tune_range <- sapply(FDR_range, function(target){
    left <- max(1, max(which(FDR <= target)))
    if(sum(FDR <= target) == 0 || left == length(FDR)){
      endpoint <- tune_candidates[left]
    } else{
      if(abs(diff(FDR[c(left, left+1)])) > 1e-6){
        coefs <- lm(FDR[c(left, left+1)] ~ tune_candidates[c(left, left+1)])$coefficients
        endpoint <- (target - coefs[1]) / coefs[2]
      } else{
        endpoint <- tune_candidates[left]
      }
    }
    
    return(endpoint)
  })
  
  if(sel_method == "FS"){
    # tune_seq <- round(seq(tune_range[1], min(tune_range[2], p), length.out = n_tune))
    lower_bound <- 10
    tune_seq <- round(exp(seq(log(max(lower_bound, tune_range[1])), log(min(tune_range[2], p)), length.out = n_tune)))
    # tune_seq <- round((seq(sqrt(tune_range[1]), sqrt(min(tune_range[2], p)), length.out = n_tune))^2)
    tune_seq <- pmax(1, tune_seq)
    for(i in 2:n_tune){
      if(tune_seq[i] <= tune_seq[i-1]) tune_seq[i] <- min(tune_seq[i-1] + 1, p)
    }
  } else{
    tune_seq <- exp(seq(tune_range[1], tune_range[2], length.out = n_tune))
  }
  return(tune_seq)
}

signal_calib <- function(set, X, nreps = 50,
                         alpha = 0.2, target = 0.8,
                         n_cores = 7){
  parallel <- process_parallel(n_cores)
  
  calib_method <- if(set$calib_method == "BH"){
    BH_lm_calib
  } else if(set$calib_method == "glasso"){
    glasso_calib
  } else{
    select_calib
  }
  return(calib_method(set, X, nreps, alpha, target, parallel))
}



select_calib <- function(set, X, nreps, alpha, target, parallel){
  side <- "two"
  
  if(!set$random_X){
    n <- nrow(X)
    p <- ncol(X)
  } else{
    n <- set$n
    p <- set$p
  }
  
  if(!set$random_X){
    X_list <- list(X)
  } else{
    sample_num <- nreps
    X_list <- lapply(1:sample_num, function(i){
      gene_X(set$X_type, n, p, i,
             scaling = set$scaling, model_X = set$model_X, mis_model = set$X_mismodel)$X
    })
  }
  X_sample_num <- length(X_list)
  
  beta_list <- lapply(1:nreps, function(i){
    with_seed(i, {
      beta <- genmu(p, set$pi1, 1, set$posit_type, 1)
    })
    if (side == "right"){
      beta <- abs(beta)
    } else if (side == "left"){
      beta <- -abs(beta)
    }
    return(beta)
  })
  eps_list <- lapply(1:nreps, function(i){
    with_seed(i, eval(set$noise))
    # with_seed(i, rnorm(n))
  })
  
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  sel_power <- function(mu1){
    FDP_power_list <- forall(i = 1:nreps, .options.multicore = list(preschedule = F)) %exec% {
      H0 <- beta_list[[i]] == 0            
      beta <- beta_list[[i]] * mu1
      eps <- eps_list[[i]]
      
      X <- X_list[[(i%%X_sample_num)+1]]
      
      with_seed(i, {
        if(set$calib_method == "lasso"){
          y <- X %*% beta + eps
          tune_seq <- gene_tune_seq(set, set$calib_method, mu1,
                                    FDR_range = c(alpha/2, min(1, alpha*2)),
                                    n_tune = 50, mc_size = 3)
          # selected <- as.matrix(glmnet::glmnet(X, y, lambda = tune_seq, intercept = F)$beta != 0)
        } else if(set$calib_method == "FS"){
          y <- X %*% beta + eps
          tune_seq <- 1:p
          # res <- summary(regsubsets(x = X, y = y, nvmax = p, method = "forward"))
          # selected <- t(unname(res$which[, -1]))
        } else if(set$calib_method %in% c("logistic", "ctree")){
          y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
          tune_seq <- gene_tune_seq(set, set$calib_method, mu1,
                                    FDR_range = c(alpha/2, min(1, alpha*2)),
                                    n_tune = 50, mc_size = 3)
        } else if(set$calib_method == "poisson"){
          y <- rpois(n = n, lambda = exp(X %*% beta))
          tune_seq <- gene_tune_seq(set, set$calib_method, mu1,
                                    FDR_range = c(alpha/2, min(1, alpha*2)),
                                    n_tune = 50, mc_size = 3)
        }
        selected <- select_variables(X, y, tune_seq, method = set$calib_method)
      })
      
      FDP_power <- sapply(1:NCOL(selected), function(tune_i){
        calc_FDP_power(which(selected[, tune_i]), H0)[1:2]
      })
      
      return(FDP_power)
    }
    FDP_power <- Reduce('+', FDP_power_list) / nreps
    gap <- abs(FDP_power[1, ] - alpha)
    index <- max(which(gap == min(gap)))
    power <- FDP_power[2, index]
    mean(power) - target
  }
  
  lower <- 0
  upper <- 10
  while (upper < 100){
    with_seed(1, {
      tmp <- try(uniroot(sel_power, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper * 2
    } else {
      return(tmp)
    }
  }
  upper <- 5
  while (upper > 0.01){
    with_seed(1, {
      tmp <- try(uniroot(sel_power, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper / 2
    } else {
      return(tmp)
    }
  }
  stop("desired singal strength not found")
}

glasso_calib <- function(set, X, nreps, alpha, target, parallel){
  n <- set$n
  p <- set$p
  
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  sel_power <- function(mu1){
    FDP_power_list <- forall(i = 1:nreps, .options.multicore = list(preschedule = F)) %exec% {
      X.res <- gene_X("Inv_Sparse", n, p, i,
                      scaling = F, model_X = T,
                      signal = mu1, nonnulls = p*set$pi1, mis_model = set$X_mismodel,
                      cov_seed = ifelse(set$random_cov, NA, set$cov_seed))
      X <- X.res$X
      beta_cov <- X.res$beta_cov
      H0 <- c(beta_cov[upper.tri(beta_cov)] == 0)
      
      tune_seq <- gene_tune_seq(set, set$calib_method, mu1,
                                FDR_range = c(alpha/2, min(1, alpha*2)),
                                n_tune = 50, mc_size = 3)
      
      FDP_power <- c()
      for(rho in tune_seq){
        with_seed(i, {
          glasso.select <- select_variables(X, y = NA, rho, "glasso")
        })
        FDP_power <- cbind(FDP_power, calc_FDP_power(which(glasso.select), H0)[1:2])
        # if(min(abs(FDP_power[1, ] - alpha)) < 0.1 * alpha){
        #   break
        # }
      }
      return(FDP_power)
    }
    FDP_power <- Reduce('+', FDP_power_list) / nreps
    gap <- abs(FDP_power[1, ] - alpha)
    index <- max(which(gap == min(gap)))
    power <- FDP_power[2, index]
    mean(power) - target
  }
  
  lower <- 0
  upper <- 10
  while (upper < 100){
    with_seed(1, {
      tmp <- try(uniroot(sel_power, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper * 2
    } else {
      return(tmp)
    }
  }
  upper <- 5
  while (upper > 0.01){
    with_seed(1, {
      tmp <- try(uniroot(sel_power, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper / 2
    } else {
      return(tmp)
    }
  }
  stop("desired singal strength not found")
}

## calibrate signal strength, modified from Lihua
BH_lm_calib <- function(set, X, nreps, alpha, target, parallel){
  side <- "two"
  
  if(!set$random_X){
    n <- nrow(X)
    p <- ncol(X)
  } else{
    n <- set$n
    p <- set$p
  }
  
  if(!set$random_X){
    X_list <- list(X)
    Sigma_list <- list(solve(t(X) %*% X))
  } else{
    sample_num <- 5
    X_list <- lapply(1:sample_num, function(i){
      gene_X(set$X_type, n, p, i,
             scaling = set$scaling, model_X = set$model_X)$X
    })
    Sigma_list <- lapply(X_list, function(X){
      solve(t(X) %*% X)
    })
  }
  X_sample_num <- length(X_list)
  
  beta_list <- lapply(1:nreps, function(i){
    with_seed(1, {
      beta <- genmu(p, set$pi1, 1, set$posit_type, 1)
    })
    if (side == "right"){
      beta <- abs(beta)
    } else if (side == "left"){
      beta <- -abs(beta)
    }
    return(beta)
  })
  eps_list <- lapply(1:nreps, function(i){
    with_seed(1, eval(set$noise))
    # with_seed(i, rnorm(n))
  })
  
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  BH_power <- function(mu1){
    power <- unlist(forall(i = 1:nreps, .options.multicore = list(preschedule = F)) %exec% {
      H0 <- beta_list[[i]] == 0            
      beta <- beta_list[[i]] * mu1
      eps <- eps_list[[i]]
      
      X <- X_list[[(i%%X_sample_num)+1]]
      Sigma <- Sigma_list[[(i%%X_sample_num)+1]]
      
      y <- X %*% beta + eps
      
      rejs_BH <- BH_lm(y, X, side = "two", alpha, Sigma = Sigma)$rejs
      power_sample <- calc_FDP_power(rejs_BH, H0)[2]
      
      return(power_sample)
    })
    mean(power) - target
  }
  
  lower <- 0
  upper <- 10
  while (upper < 100){
    with_seed(1, {
      tmp <- try(uniroot(BH_power, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper * 2
    } else {
      return(tmp)
    }
  }
  stop("desired singal strength not found")
}

genmu <- function(n, pi1, mu1,
                  posit_type = c("random", "rand_block5", "equi", "head"),
                  mu_type = 1:3){
  m <- ceiling(n * pi1)
  posit_type <- posit_type[1]
  mu_type <- mu_type[1]
  if (posit_type == "random"){
    inds <- sample(n, m, replace = F)
  } else if (posit_type == "rand_block5"){
    block_size <- 5
    if(n %% block_size != 0) stop("#variables not a multiple of block size 5")
    if(m > n/block_size) stop("#non-null is greater than #blocks")
    inds <- (sample(n/block_size, m, replace = F) - 1) * block_size + 1
  } else if (posit_type == "equi"){
    inds <- seq(1, n, floor(1 / pi1))[1:m]
  } else if (posit_type == "head"){
    inds <- 1:m
  }
  mu <- rep(0, n)
  altmu <- switch(mu_type,
                  `1` = rep(1, m),
                  `2` = rnorm(m),
                  `3` = rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1))
  mu[inds] <- mu1 * altmu
  mu * (1+rexp(n))/2
}


robust_noise_calib <- function(set, X, mu1, alpha = 0.05, nreps = 50, n_cores = 1){
  
  if(!set$random_X){
    n <- nrow(X)
    p <- ncol(X)
  } else{
    n <- set$n
    p <- set$p
  }
  
  if(!set$random_X){
    X_list <- list(X)
    Sigma_list <- list(solve(t(X) %*% X))
  } else{
    sample_num <- 5
    X_list <- lapply(1:sample_num, function(i){
      gene_X(set$X_type, n, p, i,
             scaling = set$scaling, model_X = set$model_X)$X
    })
    Sigma_list <- lapply(X_list, function(X){
      solve(t(X) %*% X)
    })
  }
  X_sample_num <- length(X_list)
  
  beta_list <- lapply(1:nreps, function(i){
    with_seed(1, {
      beta <- genmu(p, set$pi1, 1, set$posit_type, 1)
    })
    return(beta)
  })
  eps_list <- lapply(1:nreps, function(i){
    with_seed(i, rnorm(n))
  })
  
  parallel <- process_parallel(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  single_test_error <- function(mult){
    null_pvals <- unlist(forall(i = 1:nreps, .options.multicore = list(preschedule = F)) %exec% {
      H0 <- beta_list[[i]] == 0            
      beta <- beta_list[[i]] * mu1
      eps <- eps_list[[i]]
      
      X <- X_list[[(i%%X_sample_num)+1]]
      Sigma <- Sigma_list[[(i%%X_sample_num)+1]]
      
      # sd_mult <- rgamma(n, shape = 1/var_of_var, scale = var_of_var)
      # sd_mult <- runif(n, 1 - sqrt(12*var_of_var)/2, 1 + sqrt(12*var_of_var)/2)
      y <- X %*% beta + eps * exp(mult * X %*% beta)
      # y <- X %*% beta + rt(n, df = var_of_var)
      
      tvals <- lm_to_t(y, X, Sigma)
      pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")
      
      return(pvals[H0])
    })
    # mean(null_pvals < alpha) - 2 * alpha
    null_pvals
  }
  
  browser()
  null_pvals <- single_test_error(1)
  hist(null_pvals, breaks = seq(0, 1, by = 0.05))
  
  lower <- 0.01
  upper <- 0.5
  while (upper <= 10){
    with_seed(1, {
      tmp <- try(uniroot(single_test_error, c(lower, upper))$root)
    })
    if (class(tmp) == "try-error"){
      upper <- upper * 2
    } else {
      return(tmp)
    }
  }
  warning("desired singal strength not found. Using var_of_sd = 10")
  return(10)
}


calc_lasso_FDP <- function(lasso.fit, H0){
  FDPs <- sapply(1:length(lasso.fit$lambda), function(lambda_i){
    selected <- which(abs(lasso.fit$beta[, lambda_i]) > 0)
    nsel <- length(selected)
    false_discovery <- length(intersect(selected, which(H0)))
    FDP <- false_discovery / max(nsel, 1)
    
    return(FDP)
  })
  
  return(FDPs)
}

calc_sel_FDP <- function(selects, H0){
  FDPs <- sapply(1:NCOL(selects), function(tune_i){
    selected <- which(selects[, tune_i])
    nsel <- length(selected)
    false_discovery <- length(intersect(selected, which(H0)))
    FDP <- false_discovery / max(nsel, 1)
    
    return(FDP)
  })
  
  return(FDPs)
}

calc_sel_power <- function(selects, H0){
  powers <- sapply(1:NCOL(selects), function(tune_i){
    selected <- which(selects[, tune_i])
    nsel <- length(selected)
    false_discovery <- length(intersect(selected, which(H0)))
    power <- (nsel - false_discovery) / max(sum(!H0), 1)
    
    return(power)
  })
  
  return(powers)
}

calc_FDP_power <- function(rejs, H0, sign_predict = NULL, sign_beta = NULL){
  nrejs <- length(rejs)
  
  false_discovery <- length(intersect(rejs, which(H0)))
  true_discovery <- nrejs - false_discovery
  
  FDP <- false_discovery / max(nrejs, 1)
  power <- true_discovery / max(sum(!H0), 1)
  
  if(!is.null(sign_predict)){
    false_dir <- sum((sign_predict * sign_beta)[rejs] <= 0)
    true_dir <- nrejs - false_dir
    
    FDP_dir <- false_dir / max(nrejs, 1)
    power_dir <- true_dir / length(sign_beta)
  } else{
    FDP_dir <- NA
    power_dir <- NA
  }
  
  return(c(FDP, power, FDP_dir, power_dir))
}

calc_sample_deviate <- function(hFDR.samples, FDP.samples){
  # deviates <- abs(hFDR.samples - mean(hFDR.samples))
  deviates <- abs(hFDR.samples - mean(FDP.samples))
  # deviates <- abs(hFDR.samples - FDP.samples)
  
  return(deviates)
}

calc_deviate <- function(hFDR.samples, FDR.samples, measure = "std"){
  # deviates <- sapply(1:NROW(hFDR.samples), function(lambda_i){
  #   quantile(calc_sample_deviate(hFDR.samples[lambda_i, ], FDP.samples[lambda_i, ]), 1-alpha)
  # })
  
  if(measure == "std"){
    std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
      sqrt(var(hFDR.samples[lambda_i, ]))
    })
    low <- -std
    up <- std
    type <- "add"
  } else if(measure == "logstd"){
    std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
      sqrt(var(log(hFDR.samples[lambda_i, ] + 1e-8)))
    })
    low <- exp(-std)
    up <- exp(std)
    type <- "mult"
  } else if(measure == "quantile"){
    low <- sapply(1:NROW(hFDR.samples), function(lambda_i){
      quantile(hFDR.samples[lambda_i, ], probs = c(0.16))
    })
    up <- sapply(1:NROW(hFDR.samples), function(lambda_i){
      quantile(hFDR.samples[lambda_i, ], probs = c(0.84))
    })
    type <- "val"
  }
  
  return(list(low = low, up = up, type = type))
}

compress_dev <- function(dev){
  digits <- 3
  paste(round(dev$low, digits), round(dev$up, digits), dev$type)
}

extract_dev <- function(strs, index){
  sapply(strs, function(str){
    as.numeric(str_split(str, " ")[[1]][index])
  })
}




# compute the p-values of t-statistics
pvals_t <- function(tvals, df, side = "two"){
  if (side == "right"){
    pvals <- pt(tvals, df = df, lower.tail = FALSE)
  } else if (side == "left"){
    pvals <- pt(tvals, df = df, lower.tail = TRUE)
  } else if (side == "two"){
    pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
  }
  
  return(pvals)
}

# compute the t-vals of a linear regression test problem
lm_to_t <- function(y, X, Sigma = NULL){
  n <- NROW(X)
  p <- NCOL(X)
  
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = T, scale = F)
  
  if(is.null(Sigma)){
    Sigma <- solve(t(X) %*% X)
  }
  Xy <- t(X) %*% y
  df <- n - p - 1
  
  zvals <- Sigma %*% Xy
  sigmahat <- as.numeric(sqrt((sum(y^2) - t(Xy) %*% zvals) / df))
  tvals <- zvals / sqrt(diag(Sigma)) / sigmahat
  
  return(list(tvals = tvals, df = df))
}





# the complement of a union of disjoint intervals.
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
# univ_left, univ_right: the universal interval
interval_complement <- function(intervals, univ_left = 0, univ_right = 1){
  left <- c(univ_left, intervals$right)
  right <- c(intervals$left, univ_right)
  
  left_unique <- NULL
  right_unique <- NULL
  for(ind in 1:length(left)){
    if(left[ind] != right[ind]){
      left_unique <- c(left_unique, left[ind])
      right_unique <- c(right_unique, right[ind])
    }
  }
  
  complement <- list(left = left_unique, right = right_unique)
  
  return(complement)
}

# the intersection of two "union of disjoint intervals".
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_intersect <- function(intervals1, intervals2){
  left <- NULL
  right <- NULL
  intv1_num <- length(intervals1$left)
  intv2_num <- length(intervals2$left)
  if(min(intv1_num, intv2_num) > 0){
    for(ind1 in 1:intv1_num){
      for(ind2 in 1:intv2_num){
        low <- max(intervals1$left[ind1], intervals2$left[ind2])
        up <- min(intervals1$right[ind1], intervals2$right[ind2])
        if(low < up){
          left <- c(left, low)
          right <- c(right, up)
        }
      }
    }
  }
  intersect <- list(left = left, right = right)
  
  return(intersect)
}

# the union of two "union of disjoint intervals".
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_union <- function(intervals1, intervals2){
  complement1 <- interval_complement(intervals1, univ_left = -Inf, univ_right = Inf)
  complement2 <- interval_complement(intervals2, univ_left = -Inf, univ_right = Inf)
  
  complements_intersect <- interval_intersect(complement1, complement2)
  
  union <- interval_complement(complements_intersect, univ_left = -Inf, univ_right = Inf)
  
  return(union)
}

# intervals1 \setminus intervals2
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_minus <- function(intervals1, intervals2){
  if(is.null(intervals1$left)){
    return(list(left = NULL, right = NULL))
  } else if(is.null(intervals2$left)){
    return(list(left = intervals1$left, right = intervals1$right))
  } else{
    intv2_comp <- interval_complement(intervals2,
                                      univ_left = min(c(intervals1$left, intervals2$left)),
                                      univ_right = max(c(intervals1$right, intervals2$right)))
    minus <- interval_intersect(intervals1, intv2_comp)
    return(minus)
  }
}

# check if x is in intervals
x_in_intervals <- function(x, intervals){
  if(is.null(intervals$left)){
    return(FALSE)
  }
  int_num <- length(intervals$left)
  for(i in 1:int_num){
    if(x >= intervals$left[i] & x <= intervals$right[i]){
      return(TRUE)
    }
  }
  return(FALSE)
}



# sequence generator for a sequential iterator,
# need such format similar to the one used in foreach
iterate_seq <- function(...){
  arg <- list(...)[1]
  argname <- names(arg)
  arg <- arg[[1]]
  
  return(list(arg = arg, argname = argname))
}

# binary operator for a sequential iterator,
# need such format similar to the one used in foreach
`%do_seq%` <- function(obj, expr){
  result <- NULL
  if(length(obj$arg) > 0){
    for(item in obj$arg){
      envir <- parent.frame()
      envir[[obj$argname]] <- item
      result <- c(result, list(eval(substitute(expr), envir = envir)))
    }
  }
  return(result)
}

# load packages and prepare snips for parallel computing
process_parallel <- function(n_cores){
  if(n_cores > 1 && !requireNamespace("doParallel", quietly=T)) {
    warning("doParallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("parallel", quietly=T)) {
    warning("parallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("foreach", quietly=T)) {
    warning("foreach is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1){
    cores_avail <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if(n_cores > cores_avail){
      warning(paste("The requested number of cores is not available. Using instead", cores_avail, "cores."), immediate. = T)
      n_cores <- cores_avail
    }
  }
  if(n_cores > 1){
    doParallel::registerDoParallel(n_cores)
    
    forall <- foreach::foreach
    `%exec%` <- foreach::`%dopar%`
  } else{
    forall <- iterate_seq
    `%exec%` <- `%do_seq%`
  }
  
  return(list(iterator = forall, connector = `%exec%`))
}


error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}










try_repeat <- function(expr, default = NULL, n_times = 100){
  success <- F
  for(iter in 1:n_times){
    tryCatch({
      eval(expr, envir = parent.frame())
      success <- T
    }, error = function(msg){}, warning = function(msg){})
    if(success) break
    Sys.sleep(0.01)
  }
  if(!success) eval(default, envir = parent.frame())
  
  return(success)
}

update_count <- function(expr_name, action){
  file_name <- here("data", "temp", paste0("progress-", expr_name, ".RData"))
  
  if(action == "start"){
    iters_done <- 0
    start_time <- Sys.time()
    save(iters_done, start_time, file = file_name)
  } else if(action == "progress"){
    success <- try_repeat(load(file = file_name))
    if(success){
      iters_done <- iters_done + 1
      try_repeat(save(iters_done, start_time, file = file_name))
    } else{
      iters_done <- "(can't access progress record)"
    }
  } else if(action == "end"){
    load(file = file_name)
    iters_done <- Sys.time() - start_time  # abuse of var name
    file.remove(file_name)
  }
  
  return(iters_done)
}

print_progress <- function(expr_name, X_title, alpha, iters_done, sample_size, action){
  file_name <- here("data", "temp", paste0("progress-", expr_name, ".txt"))
  
  try_num <- ifelse(action == "progress", 100, 1)
  try_repeat(progress_record <- readLines(file_name), 
             default = {progress_record <- NULL},
             n_times = try_num)
  
  if(action == "start"){
    progress_record <- c(progress_record,
                         paste0("-------------- ", X_title, " --------------"),
                         paste0("alpha: ", alpha, " -- start: ", Sys.time()),
                         paste0("alpha: ", alpha, " -- 0/", sample_size, " done."))
    writeLines(progress_record, con = file_name)
  } else if(action == "progress"){
    if(length(progress_record) > 1){
      progress_record <- progress_record[1:(length(progress_record)-1)]
      progress_record <- c(progress_record,
                           paste0("alpha: ", alpha, " -- ", iters_done, "/", sample_size, " done."))
      try_repeat(writeLines(progress_record, con = file_name))
    }
  } else if(action == "end"){
    runtime <- round(as.double(iters_done), digits = 2)
    time_unit <- units(iters_done)
    progress_record <- c(progress_record,
                         paste0("alpha: ", alpha, " -- end: ", Sys.time()),
                         paste("--------- runtime:", runtime, time_unit, "---------"),
                         "")
    writeLines(progress_record, con = file_name)
  }
  
  
}

update_progress <- function(expr_name, X_title, alpha, sample_size, action){
  iters_done <- update_count(expr_name, action)
  print_progress(expr_name, X_title, alpha, iters_done, sample_size, action)
}

# borrowed from knockoff
# Evaluate an expression with the given random seed, then restore the old seed
with_seed = function(seed, expr) {
  seed.old = if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed'))
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}



# convert a digit to a string
digit_to_char <- function(number){
  return(str_replace(as.character(number), "\\.", "d"))
}

# convert a string to a digit
char_to_digit <- function(str){
  gsub(paste0("([0-9]+)d([0-9]+)"),
       paste0("\\1.\\2"),
       str)
}

int_to_str <- function(number, digits){
  num_char <- nchar(number)
  if(num_char > digits) return(as.character(number))
  paste0(rep(0, digits-num_char), number)
}

parse_name <- function(str){
  str <- char_to_digit(str)
  str <- str_replace(str, "_L_", "(")
  str <- str_replace(str, "_R_", ")")
  str <- str_replace(str, "_D_", "-")
  str <- str_replace(str, "_STAR", "*")
}

