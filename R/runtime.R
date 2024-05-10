#!/usr/bin/env Rscript

library(here)
library(tictoc)
library(tidyverse)

source(here("R", "experiment.R"))
source(here("R", "hFDR.R"))
source(here("R", "methods.R"))
source(here("R", "utils.R"))

expr_names <- c("fixedX-lasso", "fixedX-lasso-extend", "modelX-lasso", "modelX-logistic", "glasso")
# expr_names <- c("test")

runtime <- sapply(expr_names, function(expr){list()})

for(expr in expr_names){
  print(expr)
  source(here("R", "settings", paste0(expr, ".R")))
  method_names <- c("hFDR")
  sample_size <- 10
  
  for(iter_i in 1:length(fig_var$value)){
    
    set <- list(sel_method = sel_method, model = model,
                n = n_seq[[iter_i]], p = p_seq[[iter_i]],
                pi1 = pi1_seq[[iter_i]],
                X_type = X_types[[iter_i]],
                model_X = model_Xs[[iter_i]],
                random_X = random_Xs[[iter_i]],
                random_cov = random_covs[[iter_i]],
                scaling_X = scaling_Xs[[iter_i]],
                X_seed = X_seed,
                cov_seed = cov_seed,
                posit_type = posit_types[[iter_i]],
                X_mismodel = X_mismodel,
                y_mismodel = y_mismodel,
                target_at_alpha = target_at_alpha,
                target = targets[[iter_i]],
                calib_method = calib_methods[[iter_i]],
                noise = noises[[iter_i]],
                method_names = method_names,
                FDR_range = FDR_range,
                n_tune = n_tune)
    
    print(paste0(fig_var$name, ": ", fig_var$value[iter_i]))
    
    n <- set$n
    p <- set$p
    
    # get mu1 (and X if not random)
    X <- if(!set$random_X){
      gene_X(set$X_type, n, p, set$X_seed,
             scaling = set$scaling_X, model_X = set$model_X)$X
    } else { NA }
    
    start.time <- Sys.time()
    mu1 <- signal_calib(set, X, nreps = 20,
                        alpha = set$target_at_alpha, target = set$target,
                        n_cores = n_cores)
    tune_seq <- gene_tune_seq(set, set$sel_method, mu1, set$FDR_range, set$n_tune)
    end.time <- Sys.time()
    print(paste0("    Preparation: ", round(as.numeric(difftime(end.time, start.time, units = "secs")), digits = 1)))
    
    time.taken.samples <- sapply(1:sample_size, function(iter){
      
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
      }
      
      start.time <- Sys.time()
      calc_hFDR(X, y, tune_seq, set)
      end.time <- Sys.time()
      time.taken <- as.numeric(difftime(end.time, start.time, units = "secs"))
      
      print(paste0("    ", iter, ": ", time.taken))
      
      return(time.taken)
    })
    
    runtime[[expr]][[fig_var$value[iter_i]]] <- list(runtime = time.taken.samples, set = set)
    save(runtime, file = here("data", "runtime.RData"))
  }
}

for(expr in names(runtime)){
  print(expr)
  for(X_type in names(runtime[[expr]])){
    print(paste0("    ", X_type, ": ", round(mean(runtime[[expr]][[X_type]]$runtime), digits = 1)))
  }
}

data <- data.frame()
for(expr in names(runtime)){
  for(X_type in names(runtime[[expr]])){
    setting <- paste0(expr, " ", X_type)
    # setting <- as.integer(X_type)
    data <- rbind(data, data.frame(setting = setting, time = runtime[[expr]][[X_type]]$runtime))
  }
}
ggplot(data, aes(x = setting, y = time)) +
  geom_boxplot() +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "", y = "running time (seconds)")
ggsave(filename = here("figs", "runtime.pdf"),
       width = 5, height = 6)

# lm(log(time)~log(setting), data)$coefficients