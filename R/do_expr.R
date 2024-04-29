#!/usr/bin/env Rscript

library(here)

source(here("R", "experiment.R"))
source(here("R", "hFDR.R"))
source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))

# source(here("R", "settings", "test.R"))
# source(here("R", "settings", "fixedX-lasso.R"))
# source(here("R", "settings", "glasso.R"))
# source(here("R", "settings", "fixedX-lasso-asymp.R"))
# source(here("R", "settings", "modelX-logistic.R"))
# source(here("R", "settings", "modelX-ctree.R"))


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
              beta_permute = beta_permutes[[iter_i]],
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


method_names <- c("hFDR")
# hFDR (mean and quantile) v.s. true FDR
draw_hFDR(experiment, fig_var, band_mass = 0.95,
          method_names, multi_method_color, multi_method_shape)

method_names <- c("hFDR", "std")
# estimated std of hFDR (mean and quantile) v.s. true std of hFDR
draw_dev(experiment, fig_var, band_mass = 0.95, method_names, method_colors)



# draw_hFDR_asymp(experiment = "fixedX-", sel_methods = c("lasso"), X_types = c("IID_Normal", "X_AR"), p_seq = p_seq_org, at_tune_index = 1)

# # hFDR - FDP(mean and quantile) v.s. FDP - FDP = 0
# draw_FDP_diff(experiment, fig_var, band_mass = 0.95,
#               method_names, multi_method_color, multi_method_shape)

# # samples of hFDR v.s. true FDR
# draw_est(experiment, fig_var, sample = 1,
#          method_names, multi_method_color, multi_method_shape)
# 
# method_names <- c("hFDR", "std")
# # estimated std of hFDR (mean and quantile) v.s. true std of hFDR
# draw_dev(experiment, fig_var, band_mass = 0.95, method_names, method_colors)

# # samples of estimated std of hFDR v.s. true std of hFDR
# draw_dev_est(experiment, fig_var, lambda.index = 5, exaggerate = 1)


# # hFDR +- std v.s. true FDP for individual problems
# for(sample_id in 1:sample_size){
#   print(sample_id)
#   draw_predict(experiment, fig_var, sample = sample_id, exaggerate = 1)
#   Sys.sleep(1)
# }
# # gif of draw_predict for each individual problem
# draw_predict_all(experiment, fig_var, sample_size = 50, exaggerate = 1)





## for developing, not used now


# draw_beta_sample(experiment, fig_var, lambda.index = 5, fig_var_index = 1)

# draw_std_samples(experiment, fig_var, lambda.index = 5)

# draw_corr(experiment, fig_var, x_var = "pred")
# draw_std_samples(experiment, fig_var, lambda.index = 10)