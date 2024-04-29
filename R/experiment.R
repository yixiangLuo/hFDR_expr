library(abind)    # abind

library(foreach)
library(doParallel)


# compute fdr and power for methods on a linear problem
get_hFDR <- function(set, sample_size, n_cores,
                     expr_name, fig_title){
    
    n <- set$n
    p <- set$p
    
    # get mu1 (and X if not random)
    X <- if(!set$random_X){
      gene_X(set$X_type, n, p, set$X_seed,
             scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
    } else { NA }
    mu1 <- signal_calib(set, X, nreps = 50,
                        alpha = set$target_at_alpha, target = set$target,
                        n_cores = n_cores)
    # print(mu1)
    # X.data <- gene_X(set$X_type, n, p, -1,
    #                  set$scaling_X, set$model_X,
    #                  signal = mu1, nonnulls = p*set$pi1, -2,
    #                  mis_model = set$X_mismodel)
    # X <- X.data$X
    # print(mean(sqrt(diag(solve(t(X)%*%X)))))
    # print(mean(mu1/sqrt(diag(solve(t(X)%*%X)))))
    # browser()
    
    tune_seq <- gene_tune_seq(set, set$sel_method, mu1, set$FDR_range, set$n_tune)
    
    # var_of_var <- robust_noise_calib(set, X, mu1, alpha = 0.05, nreps = 100, n_cores = n_cores)

    # update_progress(expr_name, fig_title, alpha, sample_size, action = "start")
    
    registerDoParallel(n_cores)
    
    results <- foreach(iter = 1:sample_size, .options.multicore = list(preschedule = F)) %dopar% {
    # results <- lapply(1:sample_size, function(iter){
      
        set.seed(10000 + iter)
      
        print(paste0(iter, "-start"))
        
        if(set$random_X){
          cov_seed <- if(set$random_cov){-iter} else{set$cov_seed}
          X.data <- gene_X(set$X_type, n, p, iter,
                           set$scaling_X, set$model_X,
                           signal = mu1, nonnulls = p*set$pi1, cov_seed,
                           mis_model = set$X_mismodel)
          X <- X.data$X
          Xcov.true <- X.data$Xcov.true
        } else{
          Xcov.true <- NA
        }

        if(set$sel_method == "glasso"){
          invSigma <- solve(Xcov.true)
          invSigma[abs(invSigma) < 1e-8] <- 0
          H0 <- c(invSigma[upper.tri(invSigma)] == 0)
          y <- NA
        } else if(set$sel_method %in% c("logistic", "ctree")){
          beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
          H0 <- beta == 0
          y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
        } else if(set$sel_method == "poisson"){
          beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
          H0 <- beta == 0
          y <- rpois(n = n, lambda = exp(X %*% beta))
        } else{
          beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
          # beta <- beta * (1 + rexp(p, 1))/2
          # beta <- beta * sample(c(-1,1), p, replace = T)
          H0 <- beta == 0
          
          eval(set$beta_permute)
          
          if(!set$y_mismodel){
            y <- X %*% beta + eval(set$noise)
          } else if(set$y_mismodel == 1){
            hetero <- exp(abs(X) %*% abs(beta))
            hetero <- hetero / sqrt(sum(hetero^2)/n)  # normalize to have similar signal-to-noise ratio
            y <- X %*% beta + eval(set$noise) * hetero
          }
          
          # save(X, y, alpha, file = "debug.RData")
          
          # signals <- which(!H0)
          # for(interacts in 1:3){
          #   pair <- sample(signals, 2)
          #   y <- y + X[, pair[1]] * X[, pair[2]] * mu1 / 2
          # }
          
        }
        

        methods_hFDR <- list()

        select.res <- select_variables(X, y, tune_seq, set$sel_method)

        oracle_FDP <- calc_sel_FDP(select.res, H0)
        methods_hFDR[["FDP"]] <- oracle_FDP

        typeII_error <- 1 - calc_sel_power(select.res, H0)
        methods_hFDR[["typeII_error"]] <- typeII_error
        
        method_list <- get_method_list(X, y, set)
        side_info <- c(list(beta = beta), set, method_list$side_info)
        method_list <- method_list$methods

        for(method_i in 1:length(method_list)){
          method_hFDR <- method_list[[method_i]](X, y, tune_seq, Xcov.true)
          methods_hFDR[[names(method_list[method_i])]] <- method_hFDR
        }

        methods_hFDR <- as.data.frame(methods_hFDR) %>%
          mutate(tune_index = 1:length(tune_seq), sample_id = iter)
        
        print(paste0(iter, "-end"))

        return(list(methods_hFDR = methods_hFDR, side_info = side_info))
    }

    # update_progress(expr_name, fig_title, alpha, sample_size, action = "end")

    hFDR_res <- lapply(results, function(result){
      result$methods_hFDR
    })
    side_info <- lapply(results, function(result){
      result$side_info
    })
    results <- do.call("rbind", hFDR_res)
    
    # if(set$sel_method == "FS") tune_seq <- 1/exp(tune_seq)

    return(list(hFDR = results, tune_seq = tune_seq, side_info = side_info))
}







# if(set$sel_method == "lasso"){
#   tune_seq <- gene_lambda_seq(Sigma, beta = genmu(p, set$pi1, mu1, set$posit_type, 1),
#                               set$nlambda, sigma = 1, n, scaling = set$scaling_X)
# } else if(set$sel_method == "FS"){
#   tune_seq <- gene_FS_tune_seq(set$X_type, n, p, set$scaling_X, set$model_X,
#                                set$pi1, mu1, set$posit_type,
#                                ntune = set$nlambda)
# } else if(set$sel_method == "glasso"){
#   tune_seq <- gene_glasso_rho_seq(Sigma, n, ntune = 10)
# }