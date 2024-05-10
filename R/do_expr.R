#!/usr/bin/env Rscript

library(here)
library(hFDR)

source(here("R", "experiment.R"))
source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))

source(here("R", "settings", "test.R"))


args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  stop("One and only one experiment at a time.")
} else if (length(args) == 1) {
  source(here("R", "settings", paste0(args[1], ".R")))
}

file_name <- here("data", "temp", paste0("progress-", experiment, ".txt"))
if (file.exists(file_name)) {
  file.remove(file_name)
  invisible()
}


results <- lapply(1:length(fig_var$value), function(iter_i){

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
  result <- get_hFDR(set, sample_size, n_cores,
                     expr_name, fig_title = fig_titles[iter_i])

  save(result, set,
       file = here("data", "temp", paste0(experiment, "-", iter_i, ".RData")))

  return(result)
})

names(results) <- char_to_digit(fig_var$value)

save(results, fig_var, file = here("data", paste0(experiment, ".RData")))
for(iter_i in 1:length(fig_var$value)){
  file.remove(here("data", "temp", paste0(experiment, "-", iter_i, ".RData")))
}
