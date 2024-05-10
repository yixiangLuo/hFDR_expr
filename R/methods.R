library(rstan)
library(leaps)
library(lars)
# library(rstanarm)


calc_hFDR <- function(X, y, tune_seq, settings){
    
    method_names <- settings$method_names
    select <- settings$sel_method
    model <- settings$model
    
    if(model == "modelX"){
      X.mean <- rep(0, p)
      X.cov <- Xcov.true
      
      modelX <- list(
        sampler.modelX = function(X, ind, sample_size, storage = NULL){
          sampler.modelX.gauss(X, ind, X.mean = X.mean, X.cov = X.cov, sample_size, storage)
        },
        pred_fit = function(X.new, X, y, lambda){
          model_predict(X.new, X, y, lambda, method = settings$sel_method, use_mle = T)
        }
      )
      psi <- function(X, y, variables){
        psi.modelX.gauss(X, y, variables, X.mean = X.mean, X.cov = X.cov, threshold = 0.1)
      }
    } else{
      modelX <- NULL
      psi <- "pval"
    }
    
    if(!(select %in% c("lasso", "fs"))){
      if(model == "gaussgraph"){
        select <- function(X, lambda){
          select_variables(X, NULL, lambda, method = settings$sel_method)
        }
      } else{
        select <- function(X, y, lambda){
          select_variables(X, y, lambda, method = settings$sel_method)
        }
      }
    }
    
    hFDR.res <- hFDR(X, y, model = model, modelX = modelX, select = select,
                     lambda = tune_seq, psi = psi, se = ("std" %in% method_names),
                     n_sample.hfdr = 20, n_sample.se = 10, n_cores = 1)
    
    hFDR_res <- list()
    hFDR_res$hFDR <- hFDR.res$hFDR
    if("std" %in% method_names){
      hFDR_res$std <- hFDR.res$hFDR.se
    }
    
    return(hFDR_res)
}


VALID_METHOD_NAMES <- c("lasso", "elastic", "fs", "LARS", "logistic", "ctree", "poisson", "glasso")

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
  } else if(method == "fs"){
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


model_predict <- function(X_new, X, y, tune_seq, method = VALID_METHOD_NAMES, use_mle = T){
  method <- match.arg(method)
  tol <- 1e-6
  
  if(method %in% c("lasso", "elastic", "fs", "LARS")){
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
    } else if(method == "fs"){
      selected <- forward_stepwise(X, y, tune_seq)
    } else if(method == "LARS"){
      res <- least_angle(X, y, tune_seq)
      selected <- res$selected
      res <- res$lar.obj
    }
    
    prediction <- sapply(1:NCOL(selected), function(tune_i){
      if(use_mle | method == "fs"){
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


forward_stepwise <- function(X, y, n_sels, use_leaps = T, force.out = NULL){
    n <- NROW(X)
    p <- NCOL(X)
    max_step <- max(n_sels)
    if(max_step > p) stop("Cannot select more than NCOL(X) variables")
    
    X <- scale(X, center = T, scale = F)
    y <- scale(y, center = T, scale = F)
    # # QR <- qr(X)
    # # pivot_back <- sort.list(QR$pivot)
    # # Q <- qr.Q(QR, complete = F)
    # # X <- qr.R(QR, complete = F)[, pivot_back]
    # # y <- as.vector(matrix(y, nrow=1) %*% Q)
    X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
    
    if(use_leaps){
        max_step <- max(n_sels)
        res <- summary(regsubsets(x = X, y = y, nvmax = max_step, method = "forward"))
        # dof <- (n-2):(n-(1+max_step))
        # F_stat <- -diff(c(var(y)*n, res$rss)) / (res$rss / dof)
        # F_pvals <- sapply(1:length(F_stat), function(step_i){
        #     pf(F_stat[step_i], 1, dof[step_i], lower.tail = F)
        # })
        sel_mat <- t(unname(res$which[n_sels, -1]))
    } else{
        # corr <- abs(as.vector(matrix(y, nrow=1) %*% X))
        # thres <- sort(corr, decreasing = T)[min(p, 2*max_step)]
        # candidates <- which(corr >= thres)
        candidates <- 1:p
        sel_mat <- matrix(F, nrow = p, ncol = max_step)
        
        for(step in 1:max_step){
            sel <- which.max(abs(as.vector(matrix(y, nrow=1) %*% X[, candidates])))
            X_sel <- X[, candidates[sel]]
            sel_mat[candidates[sel], step:max_step] <- T
            
            candidates <- candidates[-sel]
            y <- y - sum(y * X_sel) * X_sel
            for(j in candidates){
                X[, j] <- X[, j] - sum(X[, j] * X_sel) * X_sel
                X[, j] <- X[, j] / sqrt(sum(X[, j]^2))
            }
        }
        sel_mat <- sel_mat[, n_sels]
    }
    
    return(sel_mat)
}

forward_stepwise_backup <- function(X, y, n_sels){
    # https://cran.r-project.org/web/packages/StepReg/StepReg.pdf
    
    p <- NCOL(X)
    max_step <- max(n_sels)
    if(max_step > p) stop("Cannot select more than NCOL(X) variables")
    
    data <- cbind(data.frame(y = y), data.frame(unname(X)))
    intercept_only <- lm(y ~ 1, data = data)
    all <- lm(y ~ ., data = data)
    forward <- step(intercept_only, direction = "forward",
                    scope = formula(all), trace = 0, steps = max_step, k = 0)
    new_sel <- forward$anova$Step[-1]
    new_sel <- unname(sapply(substring(new_sel, 4), as.integer))
    
    sel_mat <- matrix(F, nrow = p, ncol = max_step)
    for(step in 1:max_step){
        sel_mat[new_sel[step], step:max_step] <- T
    }
    sel_mat <- sel_mat[, n_sels]
    
    return(sel_mat)
}


least_angle <- function(X, y, tune_seq){
    n <- NROW(X)
    p <- NCOL(X)
    n_tune <- length(tune_seq)
    
    tune_seq <- tune_seq * sqrt(n)
    
    res <- lars(X, y, type = "lar", intercept = T, normalize = T)
    
    entry_lambda <- sapply(1:p, function(j){
      if(res$entry[j] > 0) res$lambda[res$entry[j]] else Inf
    })
    selected <- sapply(1:n_tune, function(tune_i){
      entry_lambda >= tune_seq[tune_i]
    })
    # a0_step <- drop(res$lar.obj$beta[1, ])
    # beta_step <- t(res$lar.obj$beta[-1, ])
    
    # # be careful about the prediction function in cv.predict()
    # predict <- function(lar.obj, X_new, tune_seq){
    #   predict(lar.obj, newx = X_new, s = tune_seq * sqrt(n), mode = "lambda")$fit
    # }
    
    return(list(selected = selected, lar.obj = res))
}


# weighted BH method
BH_weighted <- function(pvals, alpha,
                        weights = rep(1, length(pvals)) / length(pvals)){
    n <- length(pvals)
    
    adjust_pvals <- sort(pvals / weights) / (1:n)
    nrejs <- max(0, which(adjust_pvals <= alpha))
    
    rejs <- which(pvals <= nrejs * alpha * weights)
    
    return(list(nrejs = nrejs, rejs = rejs))
}


# apply BH method to a linear regression test problem
BH_lm <- function(y, X, side = "two", alpha,
                  weights = rep(1, NCOL(X)) / NCOL(X),
                  Sigma = NULL){
    
    t_result <- lm_to_t(y, X, Sigma)
    pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")
    
    BH_result <-  BH_weighted(pvals, alpha, weights)
    
    return(BH_result)
}