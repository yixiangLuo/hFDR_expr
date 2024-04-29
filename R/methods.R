library(rstan)
library(leaps)
library(lars)
# library(rstanarm)


get_multi_method_list <- function(X, y, settings){
    
    method_names <- settings$method_names
    sel_method <- settings$sel_method
    model <- settings$model
    methods <- list()
    
    n <- settings$n
    p <- settings$p
    
    guess_null_by <- "p-value"
    
    if("hFDR" %in% method_names){
        methods$hFDR <- function(X, y, tune_seq, Xcov.true = NA){
            hFDR <- if(model == "fixed-X"){
                hFDR.fixedX(X, y, tune_seq, guess_null_by, method = sel_method)
            } else if(model == "model-X"){
                hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by, method = sel_method)
            } else if(model == "graphGauss" && sel_method == "glasso"){
                hFDR.glasso(X, tune_seq, guess_null_by)
            } else stop()
            
            return(hFDR)
        }
    }
    
    if("std" %in% method_names){
        methods$std <- function(X, y, tune_seq, Xcov.true = NA){
            # std <- if(model == "fixed-X" && sel_method == "lasso"){
            #     lasso.fit <- glmnet::glmnet(X, y, lambda = tune_seq,
            #                                 intercept = F, standardize = F,
            #                                 standardize.response = F, family = "gaussian")
            #     glmnet.pack <- pack_glmnet(X, y, lasso.fit)
            #     data.pack <- process_data(X, y)
            #     
            #     est_lasso_hFDR_std(glmnet.pack, data.pack,
            #                        guess_null_by, measure = "std",
            #                        couple = F)$deviate.est$up
            # } else{
            #     # std.hFDR.bs(X, y, function(X, y){
            #     #     if(sel_method == "lasso"){
            #     #         hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by, method = "lasso")
            #     #     } else if(sel_method == "FS"){
            #     #         if(n > p){
            #     #             hFDR.fixedX(X, y, tune_seq, guess_null_by, method = "FS")
            #     #         } else{
            #     #             hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by, method = "FS")
            #     #         }
            #     #     } else if(sel_method == "glasso"){
            #     #         hFDR.glasso(X, tune_seq, guess_null_by)
            #     #     }
            #     # }, mc_size = 20)
            #     if(model == "fixed-X"){
            #         std.hFDR.fixedX(X, y, tune_seq, guess_null_by,
            #                         method = sel_method, mc_size = 20, n_cores = 1)
            #     } else if(model == "model-X"){
            #         std.hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by,
            #                         method = sel_method,
            #                         mc_size = 20, n_cores = 1)
            #     } else if(model == "graphGauss" && sel_method == "glasso"){
            #         std.hFDR.GG(X, tune_seq, guess_null_by,
            #                     method = sel_method,
            #                     mc_size = 20, n_cores = 1)
            #     } else stop()
            # }
            std <- if(model == "fixed-X"){
                std.hFDR.fixedX(X, y, tune_seq, guess_null_by,
                                method = sel_method, n_cores = 1)
            } else if(model == "model-X"){
                std.hFDR.modelX(X, y, Xcov.true, tune_seq, guess_null_by,
                                method = sel_method,
                                n_cores = 1)
            } else if(model == "graphGauss" && sel_method == "glasso"){
                std.hFDR.GG(X, tune_seq, guess_null_by,
                            method = sel_method,
                            n_cores = 1)
            } else stop()
            
            return(compress_dev(list(low = -std, up = std, type = "add")))
        }
    }
    
    for(method_name in setdiff(method_names, c("hFDR", "std"))){
        methods[[method_name]] <- function(X, y, tune_seq){
            hFDR.std.est <- est_lasso_hFDR_std(glmnet.pack, data.pack, guess_null_by, method_name, couple = F)
            hFDR.std.est <- compress_dev(hFDR.std.est)
            return(hFDR.std.est)
        }
    }
    
    return(list(methods = methods, side_info = NA))
}

multi_method_color <- c("black", "#e41a1c", "#377eb8", "#377eb8", "#984ea3", "red", "#984ea3", "#006d2c", "orange1", "grey", "#006d2c")
names(multi_method_color) <- c("FDP", "hFDR", "std", "Poincare", "Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara", "Efron_Stein", "Efron_Stein.oracle", "Bootstrap.Bayes")

multi_method_shape <- c(1, 19, 15, 19, 19, 19, 19, 19, 19, 19, 19)
names(multi_method_shape) <- c("FDP", "hFDR", "std", "Poincare", "Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara", "Efron_Stein", "Efron_Stein.oracle", "Bootstrap.Bayes")



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