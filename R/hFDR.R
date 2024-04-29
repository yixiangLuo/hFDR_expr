library(here)
library(glasso)
library(rpart)
library(CVXR)
library(latex2exp)

source(here("R", "lasso_hFDR.R"))
source(here("R", "methods.R"))
source(here("R", "hFDR_FS.R"))
source(here("R", "utils.R"))

VALID_METHOD_NAMES <- c("lasso", "elastic", "FS", "LARS", "logistic", "ctree", "poisson", "glasso")

hFDR.fixedX <- function(X, y, tune_seq, guess_null_by,
                        method = VALID_METHOD_NAMES,
                        mc_size = 20, n_cores = 1, decompose = FALSE){
  
  method <- match.arg(method)
  p <- NCOL(X)
  n <- NROW(X)
  n_tune <- length(tune_seq)
  
  if(method == "lasso"){
    y <- scale(y, center = T, scale = F)
    X <- scale(X, center = T, scale = F)
    X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
    
    lasso.fit <- glmnet::glmnet(X, y, lambda = tune_seq,
                                intercept = T, standardize = T,
                                family = "gaussian")
    
    glmnet.pack <- pack_glmnet(X, y, lasso.fit)
    data.pack <- process_data(X, y)
    
    hFDR.res <- calc_lasso_hFDR(glmnet.pack, data.pack, guess_null_by, couple = T, n_cores = n_cores)
    
    if(decompose == T) {
      hFDR <- hFDR.res$DRj_lambda
    } else {
      hFDR <- hFDR.res$hFDR_lambda
    }
  } else if(method == "FS" && guess_null_by == "p-value"){
    hFDR <- FS_hFDR.fixedX(X, y, tune_seq, guess_null_by, decompose = decompose)
  } else{
    
    parallel <- process_parallel(n_cores)
    forall <- parallel$iterator
    `%exec%` <- parallel$connector
    
    if(guess_null_by == "p-value"){
      guess.res <- guess_null(X, y, tune_seq, Xcov.true = NA, method = method,
                              model = "fixed-X")
      
      hNulls <- guess.res$nulls
      filter_prob_unselect <- guess.res$probs
      Rs <- guess.res$Rs
      # Rs <- matrix(1, p, n_tune)
    } else if(guess_null_by == "selection"){
      filter_ind <= ceiling(2/3 * length(tune_seq))
      
      filter.res <- select_variables(X, y, tune_seq[filter_ind], method)
      hNulls <- matrix(rep(!filter.res, n_tune), ncol = n_tune)
      filter_prob_unselect <- matrix(1, p, n_tune)
      Rs <- matrix(1, p, n_tune)
    }
    
    DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
    
    
    for(j in which(rowSums(hNulls) > 0)){
      # if(j != 6) next
      y.samples <- sampler.fixedX(X, y, j, mc_size)
      # mc_n <- pmax(1, ceiling(mc_size / Rs[j, ]))
      mc_n <- mc_size
      
      DPj_tune <- forall(mc_i = 1:mc_size, .options.multicore = list(preschedule = F)) %exec% {
        X.sample <- X
        y.sample <- y.samples[, mc_i]
        
        DPj_mc <- rep(NA, n_tune)
        tune_seq_subset <- hNulls[j, ] & (mc_i <= mc_n)
        if(sum(tune_seq_subset) > 0){
          res.sample <- select_variables(X.sample, y.sample, tune_seq[tune_seq_subset], method)
          DPj_mc[tune_seq_subset] <- res.sample[j, ] / pmax(1, colSums(res.sample))
        }
        return(DPj_mc[hNulls[j, ]])
      }
      DPj_tune <- do.call(cbind, DPj_tune)
      
      # H = X[, -j]%*%solve(t(X[, -j])%*%X[, -j])%*%t(X[, -j])
      # Xjj=H %*% X[, j]
      # Xjj = Xjj / sqrt(sum(Xjj^2))
      # Ho = diag(n) - H
      # vj = Ho %*% X[, j]
      # vj = vj / sqrt(sum(vj^2))
      # # print(t(Xjj) %*% y.samples)
      # threshold = max(abs(t(X[, -j]) %*% y))
      # print(sum(abs(t(X[, j]) %*% y.samples) >= threshold)/mc_size)
      # browser()
      DRj_tune[j, hNulls[j, ]] <- rowMeans(DPj_tune, na.rm = T)
      if(guess_null_by == "selection"){
        filter_prob_unselect[j, ] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
      }
    }
    if(decompose) {
      hFDR <- DRj_tune / filter_prob_unselect
    } else {
      hFDR <- colSums(DRj_tune / filter_prob_unselect)
    }
  }
  
  return(hFDR)
}


hFDR.modelX <- function(X, y, Xcov.true, tune_seq, guess_null_by,
                        method = VALID_METHOD_NAMES, mc_size = 20, n_cores = 1){
  
  method <- match.arg(method)
  p <- NCOL(X)
  n <- NROW(X)
  
  if(method == "FS" && guess_null_by == "p-value"){
    hFDR <- FS_hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by, mc_size)
    return(hFDR)
  }
  
  parallel <- process_parallel(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  if(guess_null_by == "p-value"){
    guess.res <- guess_null(X, y, tune_seq, Xcov.true, method,
                            model = "model-X")
    
    hNulls <- guess.res$nulls
    filter_prob_unselect <- guess.res$probs
  } else if(guess_null_by == "selection"){
    filter_ind <= ceiling(2/3 * length(tune_seq))
      
    filter.res <- select_variables(X, y, tune_seq[filter_ind], method)
    hNulls <- matrix(rep(!filter.res, n_tune), ncol = n_tune)
    filter_prob_unselect <- matrix(1, p, n_tune)
  }
  
  
  DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
  
  for(j in which(rowSums(hNulls) > 0)){
    DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
    X.samples <- sampler.modelX(X, Xcov.true, j, mc_size)
    
    DPj_tune <- forall(mc_i = 1:mc_size, .options.multicore = list(preschedule = F)) %exec% {
      X.sample <- X
      X.sample[, j] <- X.samples[, mc_i]
      y.sample <- y
      res.sample <- select_variables(X.sample, y.sample, tune_seq[hNulls[j, ]], method)
      DPj_mc <- res.sample[j, ] / pmax(1, colSums(res.sample))
      return(DPj_mc)
    }
    DPj_tune <- do.call(cbind, DPj_tune)
    
    DRj_tune[j, hNulls[j, ]] <- rowMeans(DPj_tune)
    if(guess_null_by == "selection"){
      filter_prob_unselect[j, ] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
    }
    
  }
  hFDR <- colSums(DRj_tune / filter_prob_unselect)
  
  return(hFDR)
}


hFDR.glasso <- function(X, tune_seq, guess_null_by,
                        method = c("glasso"), mc_size = 20, n_cores = 1){
  
  method <- match.arg(method)
  p <- NCOL(X)
  pair_num <- p*(p-1)/2
  
  parallel <- process_parallel(n_cores)
  forall <- parallel$iterator
  `%exec%` <- parallel$connector
  
  if(guess_null_by == "p-value"){
    guess.res <- guess_null(X, y = NA, tune_seq, Xcov.true = NA, method = method,
                            model = "graphGauss")
    
    hNulls <- guess.res$nulls
    filter_prob_unselect <- guess.res$probs
  } else if(guess_null_by == "selection"){
    filter_ind <= ceiling(2/3 * length(tune_seq))
    
    filter.res <- select_variables(X, y, tune_seq[filter_ind], method)
    hNulls <- matrix(rep(!filter.res, n_tune), ncol = n_tune)
    filter_prob_unselect <- matrix(1, p, n_tune)
  }
  
  DRj_tune <- matrix(0, nrow = pair_num, ncol = length(tune_seq))
  
  for(var_i in which(rowSums(hNulls) > 0)){
    
    # var_i = (j-1)*(j-2)/2 + i
    j <- sum(((2:p)-1)*((2:p)-2)/2 < var_i) + 1
    i <- var_i - (j-1)*(j-2)/2
    
    X.samples <- sampler.glassoX(X, i, j, mc_size)
    
    DPj_tune <- forall(mc_i = 1:mc_size, .options.multicore = list(preschedule = F)) %exec% {
      X.sample <- X
      X.sample[, i] <- X.samples$Xi[, mc_i]
      X.sample[, j] <- X.samples$Xj[, mc_i]
      res.sample <- select_variables(X.sample, NA, tune_seq[hNulls[var_i, ]], method)
      DPj_mc <- res.sample[var_i, ] / pmax(1, colSums(res.sample))
      return(DPj_mc)
    }
    DPj_tune <- do.call(cbind, DPj_tune)
    
    DRj_tune[var_i, hNulls[var_i, ]] <- rowMeans(DPj_tune)
    if(guess_null_by == "selection"){
      filter_prob_unselect[var_i, ] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
    }
  }
  hFDR <- colSums(DRj_tune / filter_prob_unselect)
  
  return(hFDR)
}


guess_null <- function(X, y, tune_seq, Xcov.true = NA, method = VALID_METHOD_NAMES,
                       model = c("fixed-X", "model-X", "graphGauss")){
  
  n <- NROW(X)
  p <- NCOL(X)
  n_tune <- length(tune_seq)
  
  method <- match.arg(method)
  model <- match.arg(model)
  
  if(model == "graphGauss"){
    vars_num <- p*(p-1)/2
  } else{
    vars_num <- p
  }
  
  thresholds <- matrix(NA, nrow = vars_num, ncol = n_tune)
  nulls <- matrix(NA, nrow = vars_num, ncol = n_tune)
  probs <- matrix(NA, nrow = vars_num, ncol = n_tune)
  Rs <- matrix(NA, nrow = vars_num, ncol = n_tune)
  
  if(model == "fixed-X"){
    tvals <- lm_to_t(y, X)
    pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")
  } else if(model == "model-X"){
    pvals <- sapply(1:p, function(j){
      trans.mat <- Xcov.true[j, -j] %*% solve(Xcov.true[-j, -j])
      mean.cond <- c(trans.mat %*% t(X[, -j]))
      cov.cond <- c(Xcov.true[j, j] - trans.mat %*% Xcov.true[-j, j])
      
      Xjy_mean.cond <- sum(y * mean.cond)
      Xjy_std.cond <- sqrt(sum(y^2) * cov.cond)
      
      pval <- 2 * pnorm(-abs(sum(X[, j] * y) - Xjy_mean.cond), sd = Xjy_std.cond)
    })
  } else if(model == "graphGauss"){
    pvals <- matrix(1, p, p)
    for(j in 2:p){
      for(i in 1:(j-1)){
        tvals <- lm_to_t(X[, i], X[, -i])
        pvals[i, j] <- pvals_t(tvals$tvals[j-1], tvals$df)
        # tvals <- lm_to_t(X[, j], X[, -j])
        # pvals[i, j] <- pvals_t(tvals$tvals[i], tvals$df)
      }
    }
    pvals <- pvals[upper.tri(pvals)]
  }
  
  for(var_i in 1:vars_num){
    for(tune_i in 1:n_tune){
      # if(model == "fixed-X"){
      #   y.Sj <- sampler.fixedX(X, y, var_i, 1)
      #   X.Sj <- X
      # } else if(model == "model-X"){
      #   X.Sj <- X
      #   X.Sj[, var_i] <- sampler.modelX(X, Xcov.true, var_i, 1)
      #   y.Sj <- y
      # } else if(model == "graphGauss"){
      #   j <- sum(((2:p)-1)*((2:p)-2)/2 < var_i) + 1
      #   i <- var_i - (j-1)*(j-2)/2
      #   X.Sj <- X
      #   X.samples <- sampler.glassoX(X, i, j, 1)
      #   X.Sj[, i] <- X.samples$Xi[, 1]
      #   X.Sj[, j] <- X.samples$Xj[, 1]
      #   y.Sj <- NA
      # }
      # selected <- select_variables(X.Sj, y.Sj, tune_seq[tune_i], method)
      # 
      # Rj <- max(1, sum(selected))
      # threshold <- 1 / (vars_num / Rj + 1)
      # Rs[var_i, tune_i] <- Rj
      threshold <- 0.1
      Rs[var_i, tune_i] <- 1
      thresholds[var_i, tune_i] <- threshold
      nulls[var_i, tune_i] <- (pvals[var_i] > threshold)
      probs[var_i, tune_i] <- 1 - threshold
    }
  }
  
  return(list(nulls = nulls, probs = probs, thresholds = thresholds, Rs = Rs, pvals = pvals))
}

select_variables <- function(X, y, tune_seq, method = VALID_METHOD_NAMES){
  method <- match.arg(method)
  n_tune <- length(tune_seq)
  
  if(method == "lasso"){
    res <- glmnet::glmnet(X, y, lambda = tune_seq,
                          intercept = T, standardize = T,
                          family = "gaussian")
    res <- as.matrix(res$beta != 0)
  } else if(method == "elastic"){
    res <- glmnet::glmnet(X, y, alpha = 0.5, lambda = tune_seq,
                          intercept = T, standardize = T,
                          family = "gaussian")
    res <- as.matrix(res$beta != 0)
  } else if(method == "FS"){
    res <- forward_stepwise(X, y, tune_seq)
  } else if(method == "LARS"){
    res <- least_angle(X, y, tune_seq)$selected
  } else if(method == "logistic"){
    res <- glmnet::glmnet(X, y, lambda = tune_seq,
                          intercept = T, standardize = T,
                          family = "binomial")
    res <- as.matrix(res$beta != 0)
  } else if(method == "ctree"){
    res <- rpart(y ~ X, method = "class", control = rpart.control(minsplit = 10, xval = 0))
    importance <- rep(0, NCOL(X))
    for(var_name in names(res$variable.importance)){
      ind <- as.integer(gsub("X", "", var_name))
      importance[ind] <- res$variable.importance[[var_name]]
    }
    res <- sapply(tune_seq, function(tune_val){
      importance >= tune_val
    })
  } else if(method == "poisson"){
    res <- glmnet::glmnet(X, y, lambda = tune_seq,
                          intercept = T, standardize = T,
                          family = "poisson")
    res <- as.matrix(res$beta != 0)
  } else if(method == "glasso"){
    S <- var(X)
    res <- sapply(tune_seq, function(rho){
      inv.Sigma.est <- glasso(S, rho = rho)$wi
      inv.Sigma.est[upper.tri(inv.Sigma.est)] != 0
    })
  } else{
    stop()
  }
  res <- matrix(c(res), ncol = n_tune)
  
  return(res)
}


sampler.fixedX <- function(X, y, j, sample_size){
  sampler.core(y, X[, -j], sample_size)
}

sampler.modelX <- function(X, Sigma, j, sample_size){
  n <- NROW(X)
  
  trans.mat <- Sigma[j, -j] %*% solve(Sigma[-j, -j])
  mean.cond <- c(trans.mat %*% t(X[, -j]))
  cov.cond <- c(Sigma[j, j] - trans.mat %*% Sigma[-j, j])
  
  Xj.sample <- mean.cond + sqrt(cov.cond) * matrix(rnorm(n * sample_size), nrow = n)
}

sampler.glassoX <- function(X, i, j, sample_size){
  Xi.samples <- sampler.core(X[, i], X[, -c(i, j)], sample_size)
  Xj.samples <- sampler.core(X[, j], X[, -c(i, j)], sample_size)
  
  return(list(Xi = Xi.samples, Xj = Xj.samples))
}

sampler.core <- function(var, basis, sample_size){
  n <- length(var)
  
  var.proj <- lm(var ~ basis + 1)$fitted.values
  
  radius <- sqrt(sum(var^2) - sum(var.proj^2))
  
  var.sample <- matrix(rnorm(n * sample_size), nrow = n)
  
  var.sample.proj <-lm(var.sample ~ basis + 1)$fitted.values
  
  var.sample <- var.sample - var.sample.proj
  var.sample <- scale(var.sample, center = FALSE, scale = sqrt(colSums(var.sample^2)) / radius)
  var.sample <- var.proj + var.sample
  
  subspace_samples <- is.na(colSums(var.sample))
  if(sum(subspace_samples) > sample_size/2) stop()
  for(bad_i in which(subspace_samples)){
    var.sample[, bad_i] <- var.sample[, sample(which(!subspace_samples), 1)]
  }
  
  return(var.sample)
}


std.hFDR.fixedX <- function(X, y, tune_seq, guess_null_by,
                            method = VALID_METHOD_NAMES,
                            mc_size = 10, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)
  n_tune <- length(tune_seq)
  
  mse.est <- cv.predict(X, y, method, tune_seq, nfold = 10, rescale = T)
  for(cand in order(mse.est$cv.mean)){
    tune_star <- tune_seq[cand]
    model <- select_variables(X, y, tune_star, method)
    if(sum(model) > 0) break
  }
  # if(method == "lasso"){
  #   lasso.cv <- cv.glmnet(X, y, nfolds = 5, intercept = F, standardize = F,
  #                         standardize.response = F, family = "gaussian")
  #   beta_hat.lasso <- lasso.cv$glmnet.fit$beta[, lasso.cv$index[1]]
  #   model <- beta_hat.lasso != 0
  # } else{
  #   mse.est <- cv.predict(X, y, method, tune_seq, nfold = 10, rescale = T)
  #   tune_star <- tune_seq[which.min(mse.est$cv.mean)]
  #   model <- select_variables(X, y, tune_star, method)
  # }
  beta_hat <- rep(0, p)
  if(sum(model) > 0){
    beta_hat[model] <- lm(y ~ X[, model] + 1)$coefficients[-1]
  }
  
  RSS <- sum(lm(y ~ X + 1)$residuals^2)
  sigma_hat <- sqrt(RSS / (n-p-1))
  
  hFDR.samples <- matrix(NA, n_tune, mc_size)
  
  for(boot_i in 1:mc_size){
    X.boot <- X
    y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
    
    hFDR.boot <- hFDR.fixedX(X.boot, y.boot, tune_seq, guess_null_by,
                             method, n_cores = n_cores)
    
    hFDR.samples[, boot_i] <- hFDR.boot
  }
  
  std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })
  
  return(std)
}

std.hFDR.modelX <- function(X, y, Xcov.true, tune_seq, guess_null_by,
                            method = VALID_METHOD_NAMES,
                            mc_size = 10, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)
  n_tune <- length(tune_seq)
  
  hFDR.samples <- matrix(NA, n_tune, mc_size)
  
  # pvals <- guess_null(X, y, tune_seq, Xcov.true, method, model = "model-X")$pvals
  # nulls.guess <- which(pvals > 0.1)
  
  # nulls.guess <- (1:p)[-BH_weighted(pvals, alpha = 0.5)$rejs]]
  
  mse.est <- cv.predict(X, y, method, tune_seq, nfold = 10, rescale = F)
  for(cand in order(mse.est$cv.mean)){
    tune_star <- tune_seq[cand]
    model <- select_variables(X, y, tune_star, method)
    if(sum(model) > 0) break
  }
  nulls.guess <- which(!model)
  
  # nulls.guess <- intersect(which(!model), nulls.guess)
  
  n_nulls <- length(nulls.guess)
  
  if(n_nulls > 0){
    trans.mat <- Xcov.true[nulls.guess, -nulls.guess] %*% solve(Xcov.true[-nulls.guess, -nulls.guess])
    cov.cond <- Xcov.true[nulls.guess, nulls.guess] - trans.mat %*% Xcov.true[-nulls.guess, nulls.guess]
    R.cond <- chol(cov.cond)
  }
  
  
  for(boot_i in 1:mc_size){
    bs.samples <- sample(1:n, n, replace = T)
    X.boot <- X[bs.samples, ]
    if(n_nulls > 0 && n_nulls < p){
      mean.cond <- X.boot[, -nulls.guess] %*% t(trans.mat)
      X.boot[, nulls.guess] <- mean.cond + matrix(rnorm(n*n_nulls), n) %*% R.cond
    }
    y.boot <- y[bs.samples]
    
    hFDR.boot <- hFDR.modelX(X.boot, y.boot, Xcov.true, tune_seq, guess_null_by,
                             method, n_cores = n_cores)
    
    hFDR.samples[, boot_i] <- hFDR.boot
  }
  
  std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })
  
  return(std)
}

std.hFDR.GG <- function(X, tune_seq, guess_null_by,
                        method = c("glasso"),
                        mc_size = 10, n_cores = 1){
  p <- NCOL(X)
  n <- NROW(X)
  pair_num <- p*(p-1)/2
  n_tune <- length(tune_seq)
  
  mse.est <- cv.ggraph(X, method, tune_seq, use_mle = T)
  tune_star <- tune_seq[which.min(mse.est$cv.mean)]
  model <- select_variables(X, y = NA, tune_star, method)
  
  hPrecision <- gg_mle_constrained_precise(X, model)
  hSigma <- base::solve(hPrecision)
  R <- chol(hSigma)
  
  hFDR.samples <- matrix(NA, n_tune, mc_size)
  
  for(boot_i in 1:mc_size){
    X.boot <- matrix(rnorm(n*p), n) %*% R
    
    hFDR.boot <- hFDR.glasso(X.boot, tune_seq, guess_null_by,
                             method, n_cores = n_cores)
    
    hFDR.samples[, boot_i] <- hFDR.boot
  }
  
  std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
    sqrt(var(hFDR.samples[lambda_i, ]))
  })
  
  return(std)
}

gg_mle_constrained_precise <- function(X, model){
  p <- NCOL(X)
  n <- NROW(X)
  
  # https://cvxr.rbind.io/cvxr_examples/cvxr_sparse_inverse_covariance_estimation/
  S <- cov(X)
  hPrecision <- Variable(p, p, PSD = TRUE)
  obj <- Maximize(log_det(hPrecision) - matrix_trace(S %*% hPrecision))
  constraints <- lapply(which(!model), function(var_i){
    j <- sum(((2:p)-1)*((2:p)-2)/2 < var_i) + 1
    i <- var_i - (j-1)*(j-2)/2
    hPrecision[i, j] == 0
  })
  problem <- Problem(obj, constraints)
  result <- psolve(problem)
  hPrecision <- result$getValue(hPrecision)
  
  return(hPrecision)
}

gg_mle_constrained <- function(X, model){
  p <- NCOL(X)
  n <- NROW(X)
  
  glasso.res <- glasso(var(X), rho = 0)
  hPrecision <- glasso.res$wi
  hPrecision <- (t(hPrecision) + hPrecision)/2
  for(var_i in which(!model)){
    j <- sum(((2:p)-1)*((2:p)-2)/2 < var_i) + 1
    i <- var_i - (j-1)*(j-2)/2
    hPrecision[i, j] <- 0
    hPrecision[j, i] <- 0
  }
  eigvals <- eigen(hPrecision)$values
  min_eig <- min(eigvals)
  if(min_eig <= 0){
    hPrecision <- hPrecision + ((max(eigvals)-min_eig)*1e-5 - min_eig) * diag(p) # make condition number 1e5
  }
  
  return(hPrecision)
}


std.hFDR.bs <- function(X, y, FDR_estimator, mc_size = 20){
  n <- NROW(X)
  set.seed(1)
  
  hFDR.bs <- sapply(1:mc_size, function(mc_i){
    samples.mc <- sample(n, n, replace = T)
    X.bs <- X[samples.mc, ]
    if(!any(is.na(y))) y.bs <- y[samples.mc]
    
    FDR_estimator(X.bs, y.bs)
  })
  std.est <- sqrt(rowSums((hFDR.bs - rowMeans(hFDR.bs))^2) / (mc_size - 1))
  
  return(std.est)
}


model_predict <- function(X_new, X, y, tune_seq, method = VALID_METHOD_NAMES, use_mle = T){
  method <- match.arg(method)
  tol <- 1e-6
  
  if(method %in% c("lasso", "elastic", "FS", "LARS")){
    if(method == "lasso"){
      res <- glmnet::glmnet(X, y, lambda = tune_seq,
                            intercept = T, standardize = T,
                            family = "gaussian")
      selected <- (abs(res$beta) > tol)
      # prediction <- predict(res, newx = X_new, s = tune_seq, type = "response")
      # return(prediction)
    } else if(method == "elastic"){
      res <- glmnet::glmnet(X, y, alpha = 0.5, lambda = tune_seq,
                            intercept = T, standardize = T,
                            family = "gaussian")
      selected <- (abs(res$beta) > tol)
    } else if(method == "FS"){
      selected <- forward_stepwise(X, y, tune_seq)
    } else if(method == "LARS"){
      res <- least_angle(X, y, tune_seq)
      selected <- res$selected
      res <- res$lar.obj
    }
    
    prediction <- sapply(1:NCOL(selected), function(tune_i){
      if(use_mle | method == "FS"){
        coef <- rep(0, NROW(selected))
        a0 <- 0
        model <- selected[, tune_i]
        if(sum(model) > 0){
          lm.coefs <- lm(y ~ X[, model] + 1)$coefficients
          coef[model] <- lm.coefs[-1]
          a0 <- lm.coefs[1]
        }
        if(any(is.na(coef))){
          coef[is.na(coef)] <- 0
          warning("singular X in CV")
        }
        X_new %*% coef + a0
      } else{
        if(method != "LARS"){
          X_new %*% res$beta[, tune_i] + res$a0[tune_i]
        }
        else{
          n <- NROW(X)
          predict(res, newx = X_new, s = tune_seq[tune_i] * sqrt(n), mode = "lambda")$fit
        }
      }
    })
  }
  else if(method %in% c("logistic", "poisson")){
    family <- if(method == "logistic") "binomial" else if(method == "poisson") "poisson" else stop()
    res <- glmnet::glmnet(X, y, lambda = tune_seq,
                          intercept = T, standardize = T,
                          family = family)
    selected <- (abs(res$beta) > tol)
    if(use_mle){
      prediction <- sapply(1:NCOL(selected), function(tune_i){
        model <- selected[, tune_i]
        if(sum(model) <= 1){
          return(rep(Inf, NROW(X_new)))
        }
        mle <- glmnet::glmnet(X[, model], y, lambda = 0,
                              intercept = T, standardize = T,
                              family = family)
        c(predict(mle, newx = X_new[, model], s = 0, type = "response"))
      })
    } else{
      prediction <- predict(res, newx = X_new, s = tune_seq, type = "response")
    }
    
  }
  else{
    stop()
  }
  
  return(prediction)
}

cv.predict <- function(X, y, method, tune_seq, nfold = 10, rescale = F, use_mle = T, fit_model = F){
  n <- NROW(X)
  p <- NCOL(X)
  
  batch_size <- round(n / nfold)
  batches <- lapply(1:nfold, function(f_i){
    if(f_i < nfold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
    else ((f_i-1) * batch_size + 1) : n
  })
  cv.errors <- sapply(1:nfold, function(f_i){
    test_batch <- batches[[f_i]]
    X.tr <- X[-test_batch, ]
    y.tr <- y[-test_batch]
    
    if(rescale){  
      # X.tr <- scale(X.tr, center = T, scale = F)
      X.tr <- scale(X.tr, center = F, scale = sqrt(colSums(X.tr^2) / colSums(X^2)))
    }
    
    prediction <- model_predict(X[test_batch, ], X.tr, y.tr, tune_seq, method, use_mle)
    
    test_error <- sapply(1:NCOL(prediction), function(tune_i){
      mean((y[test_batch] - prediction[, tune_i])^2, na.rm = T)
    })
    
    return(test_error)
  })
  MSE.mean <- rowMeans(cv.errors)
  MSE.std <- sapply(1:NROW(cv.errors), function(row){
    sqrt(var(cv.errors[row, ]) / nfold)
  })
  
  ind.min <- which.min(MSE.mean)
  ind.1se <- min(which(MSE.mean[1:ind.min] <= MSE.mean[ind.min] + MSE.std[ind.min]))
  
  selected <- if(fit_model){
    select_variables(X, y, tune_seq, method)
  } else NA
  
  cv.predict <- structure(list(call = match.call(),
                               X = X,
                               y = y,
                               tune_seq = tune_seq,
                               method = method,
                               cv.mean = MSE.mean,
                               cv.std = MSE.std,
                               ind.min = ind.min,
                               ind.1se = ind.1se,
                               selected = selected,
                               use_mle = use_mle),
                          class = 'cv.obj')
  
  return(cv.predict)
}

cv.ggraph <- function(X, method, tune_seq, nfold = 10, use_mle = F, fit_model = F){
  if(method != "glasso") stop()
  
  n <- NROW(X)
  p <- NCOL(X)
  
  batch_size <- round(n / nfold)
  batches <- lapply(1:nfold, function(f_i){
    if(f_i < nfold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
    else ((f_i-1) * batch_size + 1) : n
  })
  
  cv.res <- lapply(tune_seq, function(rho){
    cv.errors <- sapply(1:nfold, function(f_i){
      test_batch <- batches[[f_i]]
      X.bs <- X[-test_batch, ]
      X.bs <- scale(X.bs, center = T, scale = F)
      X.bs <- scale(X.bs, center = F, scale = sqrt(colSums(X^2) / (NROW(X) - 1)))
      
      if(use_mle){
        model <- select_variables(X.bs, y = NA, rho, method)
        hPrecision <- gg_mle_constrained(X.bs, model)
        hSigma <- base::solve(hPrecision)
      } else{
        S <- var(X.bs)
        glasso.res <- glasso(S, rho)
        hPrecision <- glasso.res$wi
        hSigma <- glasso.res$w
      }
      
      X.predict <- sapply(1:p, function(j){
        hSigma[j, -j] %*% solve(hSigma[-j, -j]) %*% t(X[test_batch, -j])
      })
      MSE <- mean(colMeans(X[test_batch, ] - X.predict)^2)
      
      logP <- -(log(det(hPrecision)) - sum(diag(var(X[test_batch, ]) %*% hPrecision)))
      
      return(c(MSE, logP))
    })
    MSE.mean <- mean(cv.errors[1, ])
    MSE.std <- sqrt(var(cv.errors[1, ]) / nfold)
    
    logP.mean <- mean(cv.errors[2, ])
    logP.std <- sqrt(var(cv.errors[2, ]) / nfold)
    
    return(list(MSE.mean = MSE.mean, MSE.std = MSE.std,
                logP.mean = logP.mean, logP.std = logP.std))
  })
  
  cv.MSE <- list(mean = sapply(cv.res, function(res){res$MSE.mean}),
                        std = sapply(cv.res, function(res){res$MSE.std}))
  cv.logP <- list(mean = sapply(cv.res, function(res){res$logP.mean}),
                         std = sapply(cv.res, function(res){res$logP.std}))
  
  ind.min <- which.min(cv.logP$mean)
  ind.1se <- min(which(cv.logP$mean[1:ind.min] <= cv.logP$mean[ind.min] + cv.logP$std[ind.min]))
  
  selected <- if(fit_model){
    select_variables(X, y, tune_seq, method)
  } else NA
  
  cv.ggraph <- structure(list(call = match.call(),
                              X = X,
                              tune_seq = tune_seq,
                              method = method,
                              cv.mean = cv.logP$mean,
                              cv.std = cv.logP$std,
                              ind.min = ind.min,
                              ind.1se = ind.1se,
                              selected = selected,
                              cv.MSE = cv.MSE),
                         class = 'cv.obj')
  
  return(cv.ggraph)
}

## decompose only implemented for fixed-X method
hFDR_on_cv <- function(cv.obj, model = c("fixed-X", "model-X", "graphGauss"), n_tune = 20, n_cores = 1, decompose = FALSE, skip_se = FALSE){
  if(!all(class(cv.obj) == "cv.obj")){
    stop()
  }
  
  X <- cv.obj$X
  y <- cv.obj$y
  method <- cv.obj$method
  
  n_tune.cv <- length(cv.obj$tune_seq)
  dist <- cv.obj$ind.min - cv.obj$ind.1se
  ind.end <- max(cv.obj$ind.min + 3*dist, (n_tune.cv + 0*cv.obj$ind.min) %/% 1)
  ind.end <- min(ind.end, n_tune.cv)
  ind.dense_end <- min(cv.obj$ind.min + max(dist, round(cv.obj$ind.min/3)), ind.end)
  n_tune.dense <- round(ind.dense_end*2 / (ind.dense_end + ind.end) * n_tune)
  tune.ind <- c(round(seq(from = 1, to = ind.dense_end, length.out = n_tune.dense)),
                round(seq(from = ind.dense_end, to = ind.end, length.out = n_tune-n_tune.dense+1)))
  tune.ind <- unique(tune.ind)
  tune_seq <- cv.obj$tune_seq[tune.ind]
  
  model <- match.arg(model)
  guess_null_by <- "p-value"
  
  if(model == "fixed-X"){
    hFDR <- hFDR.fixedX(X, y, tune_seq, guess_null_by, method, n_cores = n_cores, decompose = decompose)
    if(skip_se == TRUE) {
      std <- rep(NA,n_tune)
    } else {
      std <- std.hFDR.fixedX(X, y, tune_seq, guess_null_by, method, n_cores = n_cores)
    }
  } else if(model == "model-X"){
    Xcov.est <- cov(X)
    hFDR <- hFDR.modelX(X, y, Xcov.est, tune_seq, guess_null_by, method, n_cores = n_cores)
    if(skip_se == TRUE) {
      std <- rep(NA,n_tune)
    } else {
      std <- std.hFDR.modelX(X, y, Xcov.est, tune_seq, guess_null_by, method, n_cores = n_cores)
    }
  } else if(model == "graphGauss"){
    hFDR <- hFDR.glasso(X, tune_seq, guess_null_by, method, n_cores = n_cores)
    if(skip_se == TRUE) {
      std <- rep(NA,n_tune)
    } else {
      std <- std.hFDR.GG(X, tune_seq, guess_null_by, method, n_cores = n_cores)
    }
  } else stop()
  
  hFDR.res <- structure(list(call = match.call(),
                             X = X,
                             y = y,
                             tune_seq = tune_seq,
                             hFDR = hFDR,
                             hFDR.std = std),
                        class = 'hFDR')
  
  return(hFDR.res)
}


plot.hFDR <- function(hFDR.obj, cv.obj, errors = NA, sign.tune = -1, xlab = NA, ylab = NA,
                      show_cv = T, show_hFDR = T, show_FPR = F, show_FDR = T, show_FDP = show_FDR,
                      log_cv = T, log_x = F, show_legend = T, show_axes = T, show_extra_axes = show_axes, 
                      FDP.leg = "FDP", leg_bty="o", main = NULL){
  if(is.na(xlab)) xlab <- if(cv.obj$method %in% c("FS")) "number of steps" else expression(-log(lambda))
  if(is.na(ylab)) ylab <- expression("False Discovery Rate")
  # if(sign.tune < 0) xlab <- paste("-", xlab, sep = "")
  
  n <- NROW(cv.obj$X)
  trans_x <- if(!(cv.obj$method %in% c("FS")) | log_x) {function(x) log(n*x)} else identity
  
  
  plot.range <- range(c(hFDR.obj$hFDR-hFDR.obj$hFDR.std, hFDR.obj$hFDR+hFDR.obj$hFDR.std))
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  plot.args <- list(x = sign.tune * trans_x(hFDR.obj$tune_seq),
                    y = hFDR.obj$hFDR,
                    xlab = "", ylab = "",
                    xlim = range(sign.tune * trans_x(cv.obj$tune_seq)),
                    ylim = plot.range,
                    type = "n")
  
  # "estimated FDR and scaled CV MSE"
  # new.args <- list(...)
  # if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  # use log scale for MSE
  ep <- 1e-8
  mse_scale <- if(log_cv) function(x) log(pmax(x,ep)) else identity
  mse.mean <- mse_scale(cv.obj$cv.mean)
  mse.low <- mse_scale(cv.obj$cv.mean-cv.obj$cv.std)
  mse.up <- mse_scale(cv.obj$cv.mean+cv.obj$cv.std)
  cv.range <- range(c(mse.low, mse.up))
  cv.transfer <- function(mse){
    (mse - min(cv.range)) / abs(diff(cv.range)) * abs(diff(plot.range)) + min(plot.range)
  }
  
  leg <- c()
  lty <- c()
  lwd <- c()
  pch <- c()
  col <- c()
  
  do.call("plot", plot.args)
  
  if(show_cv){
    lines(x = sign.tune * trans_x(cv.obj$tune_seq),
          y = cv.transfer(mse.mean),
          col = "dodgerblue3")
    error.bars(sign.tune * trans_x(cv.obj$tune_seq),
               cv.transfer(mse.up), cv.transfer(mse.low),
               width = 0.01, col = alpha("dodgerblue3", 0.1))
    
    abline(v = sign.tune * trans_x(cv.obj$tune_seq[cv.obj$ind.min]), lty = 1, col = "dodgerblue3")
    abline(v = sign.tune * trans_x(cv.obj$tune_seq[cv.obj$ind.1se]), lty = 2, col = "dodgerblue3")
    
    # leg <- c(leg, "-loglikelihood")
    leg <- c(leg, "CV MSE")
    lty <- c(lty, 1)
    lwd <- c(lwd, 1)
    pch <- c(pch, 26)
    col <- c(col, "dodgerblue3")
  }
  
  if(!any(is.na(errors))){
    if(show_FPR){
      lines(x = sign.tune * trans_x(cv.obj$tune_seq),
            y = errors$FPP,
            col = "orange3", lty = 2)
      lines(x = sign.tune * trans_x(cv.obj$tune_seq),
            y = errors$FPR,
            col = "orange3")
      
      leg <- c(leg, "Type II error", "Type II error rate")
      lty <- c(lty, 2, 1)
      lwd <- c(lwd, 1, 1)
      pch <- c(pch, 26, 26)
      col <- c(col, "orange3", "orange3")
    }
    
    if(!is.null(errors$FDP) && show_FDP){
      lines(x = sign.tune * trans_x(cv.obj$tune_seq),
            y = errors$FDP,
            col = "black", lty = 2)
      
      leg <- c(leg, FDP.leg)
      lty <- c(lty, 2)
      lwd <- c(lwd, 1)
      pch <- c(pch, 26)
      col <- c(col, "black")
    }
    if(!is.null(errors$FDR) && show_FDR){
      lines(x = sign.tune * trans_x(cv.obj$tune_seq),
            y = errors$FDR,
            col = "black")
      
      leg <- c(leg, "FDR")
      lty <- c(lty, 1)
      lwd <- c(lwd, 1)
      pch <- c(pch, 26)
      col <- c(col, "black")
    }
  }
  
  
  if(show_hFDR){
    error.bars(sign.tune * trans_x(hFDR.obj$tune_seq),
               hFDR.obj$hFDR+hFDR.obj$hFDR.std, hFDR.obj$hFDR-hFDR.obj$hFDR.std,
               width = 0.01, alpha("red", 0.3))
    points(sign.tune*trans_x(hFDR.obj$tune_seq), hFDR.obj$hFDR,
           pch = 20, col = "red")
    
    leg <- c(leg, TeX("$\\widehat{FDR}$"))
    lty <- c(lty, 1)
    lwd <- c(lwd, 0)
    pch <- c(pch, 19)
    col <- c(col, "red")
  }
  
  if(show_extra_axes) {
    axis(side = 3, at = sign.tune*trans_x(cv.obj$tune_seq),
         labels = paste(colSums(cv.obj$selected)), tick = FALSE, line = -0.5)
    mtext("# selections", side = 3, line = 2)
  
    cv.ticks <- if(log_cv){ seq(exp(cv.range[1]), exp(cv.range[2]), length.out = 6) }
    else { seq(cv.range[1], cv.range[2], length.out = 6) }
    axis(side = 4, at = cv.transfer(mse_scale(cv.ticks)),
         labels = formatC(cv.ticks, format = "g", digits = 2), tick = T, line = 0)
    mtext("CV MSE", side = 4, line = 2.5)  # "-log likelihood (cv)"
  }
  
  legend.pos <- legend("bottomright", inset = 0.05,
                       legend = leg, bg = "white",
                       lty = lty, lwd = lwd,
                       pch = pch, col = col, plot = F, y.intersp=1.2)[[1]]
  cap <- min((hFDR.obj$hFDR-hFDR.obj$hFDR.std)[sign.tune * trans_x(hFDR.obj$tune_seq) >= legend.pos$left])
  if(!any(is.na(errors))) cap <- min(cap, errors$FDP[sign.tune * trans_x(cv.obj$tune_seq) >= legend.pos$left])
  floor <- max(cv.transfer(mse.up[sign.tune * trans_x(cv.obj$tune_seq) >= legend.pos$left]))
  top <- if(cap > floor) (legend.pos$top + cap) / 2 else legend.pos$top
  if(show_legend) {
    legend(x = legend.pos$left, y = top, inset = 0.05,
         legend = leg, bg = "white",
         lty = lty, lwd = lwd,
         pch = pch, col = col, y.intersp=1.2, bty=leg_bty)
  }
  
  #par(mar = c(4,4,3,4), cex = 1, cex.axis = 1, cex.lab = 1)
  if(show_extra_axes) {
    title(xlab = xlab, ylab = ylab, line = 2.5, cex.lab = 1)
  } else {
    title(xlab = xlab, ylab = ylab, line = 2.5, cex.lab = 1)
    title(main = main, line = 1, cex.lab = 1)
  }
  invisible()
}






# std.hFDR.modelX <- function(X, y, Xcov.true, tune_seq, guess_null_by,
#                             method = c("lasso", "elastic", "FS"),
#                             mc_size = 20, n_cores = 1){
#   p <- NCOL(X)
#   n <- NROW(X)
#   n_tune <- length(tune_seq)
#   
#   hFDR.samples <- matrix(NA, n_tune, mc_size)
#   
#   for(boot_i in 1:mc_size){
#     bs.samples <- sample(1:n, n, replace = T)
#     X.boot <- X[bs.samples, ]
#     y.boot <- y[bs.samples]
#     
#     hFDR.boot <- hFDR.modelX(X.boot, y.boot, Xcov.true, tune_seq, guess_null_by,
#                              method, mc_size = 20, n_cores)
#     
#     hFDR.samples[, boot_i] <- hFDR.boot
#   }
#   
#   std <- sapply(1:NROW(hFDR.samples), function(lambda_i){
#     sqrt(var(hFDR.samples[lambda_i, ]))
#   })
#   
#   return(std)
# }

# std.hFDR.cv <- function(X, FDR_estimator, std.fold = 10){
#   n <- NROW(X)
#   
#   batch_size <- round(n / std.fold)
#   batches <- lapply(1:std.fold, function(f_i){
#     if(f_i < std.fold) ((f_i-1) * batch_size + 1) : (f_i * batch_size)
#     else ((f_i-1) * batch_size + 1) : n
#   })
#   hFDR.bs <- sapply(1:std.fold, function(f_i){
#     X.bs <- X[-batches[[f_i]], ]
#     X.bs <- scale(X.bs, center = T, scale = F)
#     X.bs <- scale(X.bs, center = F, scale = sqrt(colSums(X^2) / (NROW(X) - 1)))
#     
#     FDR_estimator(X.bs)
#   })
#   std.est <- sqrt(rowSums((hFDR.bs - rowMeans(hFDR.bs))^2) / (std.fold - 1))
#   
#   return(std.est)
# }
# 
# 
# 
# hFDR.fixedX <- function(X, y, tune_seq, guess_null_by,
# method = c("lasso", "elastic", "FS"),
# mc_size = 50, n_cores = 1){
#   mc_size <- 50
#   
#   method <- match.arg(method)
#   p <- NCOL(X)
#   n <- NROW(X)
#   
#   if(method == "lasso"){
#     lasso.fit <- glmnet::glmnet(X, y, lambda = tune_seq,
#                                 intercept = F, standardize = F,
#                                 standardize.response = F, family = "gaussian")
#     
#     glmnet.pack <- pack_glmnet(X, y, lasso.fit)
#     data.pack <- process_data(X, y)
#     
#     hFDR.res <- calc_lasso_hFDR(glmnet.pack, data.pack, guess_null_by, couple = T, n_cores = n_cores)
#     
#     hFDR <- hFDR.res$hFDR_lambda
#   } else{
#     
#     parallel <- process_parallel(n_cores)
#     forall <- parallel$iterator
#     `%exec%` <- parallel$connector
#     
#     if(guess_null_by == "p-value"){
#       data.pack <- process_data(X, y)
#       # guess.res <- guess_null(data.pack, method = guess_null_by)
#       
#       hNulls <- guess.res$nulls
#       filter_prob_unselect <- guess.res$probs
#     } else if(guess_null_by == "selection"){
#       filter_ind <= ceiling(2/3 * length(tune_seq))
#       
#       filter.res <- select_variables(X, y, tune_seq[filter_ind], method)
#       hNulls <- which(!filter.res)
#       filter_prob_unselect <- rep(1, p)
#     }
#     
#     DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
#     
#     
#     for(j in hNulls){
#       y.samples <- sampler.fixedX(X, y, j, mc_size)
#       
#       DPj_tune <- forall(mc_i = 1:mc_size, .options.multicore = list(preschedule = F)) %exec% {
#         X.sample <- X
#         y.sample <- y.samples[, mc_i]
#         res.sample <- select_variables(X.sample, y.sample, tune_seq, method)
#         DPj_mc <- res.sample[j, ] / pmax(1, colSums(res.sample))
#         return(DPj_mc)
#       }
#       DPj_tune <- do.call(cbind, DPj_tune)
#       
#       DRj_tune[j, ] <- rowMeans(DPj_tune)
#       if(guess_null_by == "selection"){
#         filter_prob_unselect[j] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
#       }
#     }
#     hFDR <- colSums(DRj_tune / filter_prob_unselect)
#   }
#   
#   return(hFDR)
# }
# 
# 
# hFDR.modelX <- function(X, y, Xcov.true, tune_seq, guess_null_by,
#                         method = c("lasso", "FS"), mc_size = 50){
#   
#   method <- match.arg(method)
#   p <- NCOL(X)
#   n <- NROW(X)
#   
#   if(guess_null_by == "p-value"){
#     # CRT p-value
#     threshold <- 0.5
#     pvals <- sapply(1:p, function(j){
#       trans.mat <- Xcov.true[j, -j] %*% solve(Xcov.true[-j, -j])
#       mean.cond <- c(trans.mat %*% t(X[, -j]))
#       cov.cond <- c(Xcov.true[j, j] - trans.mat %*% Xcov.true[-j, j])
#       
#       Xjy_mean.cond <- sum(y * mean.cond)
#       Xjy_std.cond <- sqrt(sum(y^2) * cov.cond)
#       
#       pval <- 2 * pnorm(-abs(sum(X[, j] * y) - Xjy_mean.cond), sd = Xjy_std.cond)
#     })
#     hNulls <- which(pvals > threshold)
#     filter_prob_unselect <- rep(1 - threshold, p)
#   } else if(guess_null_by == "selection"){
#     filter_ind <= ceiling(2/3 * length(tune_seq))
#     
#     filter.res <- select_variables(X, y, tune_seq[filter_ind], method)
#     hNulls <- which(!filter.res)
#     filter_prob_unselect <- rep(1, p)
#   }
#   
#   
#   DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
#   
#   for(j in hNulls){
#     DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
#     X.samples <- sampler.modelX(X, Xcov.true, j, mc_size)
#     
#     for(mc_i in 1:mc_size){
#       X.sample <- X
#       X.sample[, j] <- X.samples[, mc_i]
#       y.sample <- y
#       res.sample <- select_variables(X.sample, y.sample, tune_seq, method)
#       DPj_tune[, mc_i] <- res.sample[j, ] / pmax(1, colSums(res.sample))
#     }
#     
#     DRj_tune[j, ] <- rowMeans(DPj_tune)
#     if(guess_null_by == "selection"){
#       filter_prob_unselect[j] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
#     }
#     
#   }
#   hFDR <- colSums(DRj_tune / filter_prob_unselect)
#   
#   return(hFDR)
# }
# 
# 
# hFDR.glasso <- function(X, tune_seq, guess_null_by,
#                         method = c("glasso"), mc_size = 50){
#   
#   method <- match.arg(method)
#   p <- NCOL(X)
#   pair_num <- p*(p-1)/2
#   
#   if(guess_null_by == "p-value"){
#     pvals <- matrix(1, p, p)
#     for(j in 2:p){
#       for(i in 1:(j-1)){
#         tvals <- lm_to_t(X[, i], X[, -i])
#         pvals[i, j] <- pvals_t(tvals$tvals[j-1], tvals$df)
#       }
#     }
#     threshold <- 0.5
#     filter.res <- pvals[upper.tri(pvals)] <= threshold
#     filter_prob_unselect <- rep(1-threshold, pair_num)
#   } else if(guess_null_by == "selection"){
#     filter_ind <= ceiling(2/3 * length(tune_seq))
#     filter.res <- select_variables(X, y = NA, tune_seq[filter_ind], method)
#     filter_prob_unselect <- rep(1, pair_num)
#   }
#   
#   hNulls <- matrix(F, p, p)
#   hNulls[upper.tri(hNulls)] <- !filter.res
#   
#   DRj_tune <- matrix(0, nrow = pair_num, ncol = length(tune_seq))
#   
#   for(j in 2:p){
#     for(i in 1:(j-1)){
#       if(!hNulls[i,j]) next
#       
#       DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
#       X.samples <- sampler.glassoX(X, i, j, mc_size)
#       
#       for(mc_i in 1:mc_size){
#         X.sample <- X
#         X.sample[, i] <- X.samples$Xi[, mc_i]
#         X.sample[, j] <- X.samples$Xj[, mc_i]
#         res.sample <- select_variables(X.sample, NA, tune_seq, method)
#         DPj_tune[, mc_i] <- res.sample[(j-1)*(j-2)/2 + i, ] / pmax(1, colSums(res.sample))
#       }
#       
#       DRj_tune[(j-1)*(j-2)/2 + i, ] <- rowMeans(DPj_tune)
#       if(guess_null_by == "selection"){
#         filter_prob_unselect[(j-1)*(j-2)/2 + i] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
#       }
#     }
#   }
#   hFDR <- colSums(DRj_tune / filter_prob_unselect)
#   
#   return(hFDR)
# }
