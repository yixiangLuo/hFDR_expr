library(here)
library(tidyverse)

source(here("R", "cluster", "cluster_utils.R"))

# read cmd args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
  experiment <- args[1]
  job_id <- as.integer(args[2])
} else {
  stop("accept two paramters: expr_name and job id.")
}

# experiment <- "test"
# job_id <- 1

# load packages and data
source(here("R", "settings", paste0(experiment, ".R")))
load(here("data", "temp", experiment, "settings.RData"))

# the indeces of experiments this job should do
job_trial_indeces <- get_data_indeces(job_id, index_data)


# do the experiments
for(trial_index in job_trial_indeces){
  # if(file.exists(here("data", "temp", experiment, paste0(trial_index$trial_id, ".RData")))){
  #   next
  # }
  
  # record the start of the program
  cat(as.integer(Sys.time()),
      file = here("data", "temp", experiment, "progress", trial_index$trial_id))
    
  iter <- trial_index$sample_index
  set <- settings[[trial_index$expr_index]]
  
  n <- set$n
  p <- set$p
  
  # get mu1 (and X if not random)
  X <- if(!set$random_X){
    gene_X(set$X_type, n, p, set$X_seed,
                scaling = set$scaling_X, model_X = set$model_X, mis_model = set$X_mismodel)$X
  } else { NA }
  mu1 <- signal_calib(set, X, nreps = 50,
                      alpha = set$target_at_alpha, target = set$target,
                      n_cores = 1)
  tune_seq <- gene_tune_seq(set, set$sel_method, mu1, set$FDR_range, set$n_tune)
  
  
  set.seed(10000 + iter)
  
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
    } else{
      hetero <- exp(abs(X) %*% abs(beta))
      hetero <- hetero / sqrt(sum(hetero^2)/n)  # normalize to have similar signal-to-noise ratio
      y <- X %*% beta + eval(set$noise) * hetero
    }
    
    # signals <- which(!H0)
    # for(interacts in 1:3){
    #   pair <- sample(signals, 2)
    #   y <- y + X[, pair[1]] * X[, pair[2]] * mu1 / 2
    # }
    
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
  
  # save results
  save(hFDR_res, tune_seq, side_info,
       file = here("data", "temp", experiment, paste0(trial_index$trial_id, ".RData")))
  
  # delete the "in-progress" record
  unlink(here("data", "temp", experiment, "progress", trial_index$trial_id))
}







