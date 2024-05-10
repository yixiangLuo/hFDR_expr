library(abind)    # abind

library(foreach)
library(doParallel)


get_hFDR <- function(set, sample_size, n_cores,
                     expr_name, fig_title){
    
    n <- set$n
    p <- set$p
    
    X <- if(!set$random_X){
      gene_X(set$X_type, n, p, set$X_seed,
             scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
    } else { NA }
    mu1 <- signal_calib(set, X, nreps = 50,
                        alpha = set$target_at_alpha, target = set$target,
                        n_cores = n_cores)
    
    tune_seq <- gene_tune_seq(set, set$sel_method, mu1, set$FDR_range, set$n_tune)
    
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
          H0 <- beta == 0
          
          eval(set$beta_permute)
          
          if(!set$y_mismodel){
            y <- X %*% beta + eval(set$noise)
          } else if(set$y_mismodel == 1){
            hetero <- exp(abs(X) %*% abs(beta))
            hetero <- hetero / sqrt(sum(hetero^2)/n)  # normalize to have similar signal-to-noise ratio
            y <- X %*% beta + eval(set$noise) * hetero
          }
        }

        hFDR_res <- list()

        select.res <- select_variables(X, y, tune_seq, set$sel_method)

        oracle_FDP <- calc_sel_FDP(select.res, H0)
        hFDR_res[["FDP"]] <- oracle_FDP

        typeII_error <- 1 - calc_sel_power(select.res, H0)
        hFDR_res[["typeII_error"]] <- typeII_error
        
        hFDR_res <- c(hFDR_res, calc_hFDR(X, y, tune_seq, set))

        hFDR_res <- as.data.frame(hFDR_res) %>%
          mutate(tune_index = 1:length(tune_seq), sample_id = iter)
        
        side_info <- c(list(beta = beta), set)
        
        print(paste0(iter, "-end"))

        return(list(hFDR_res = hFDR_res, side_info = side_info))
    }

    # update_progress(expr_name, fig_title, alpha, sample_size, action = "end")

    hFDR_res <- lapply(results, function(result){
      result$hFDR_res
    })
    side_info <- lapply(results, function(result){
      result$side_info
    })
    results <- do.call("rbind", hFDR_res)

    return(list(hFDR = results, tune_seq = tune_seq, side_info = side_info))
}





