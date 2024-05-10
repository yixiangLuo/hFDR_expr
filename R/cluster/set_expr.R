#!/usr/bin/env Rscript

library(here)

source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "cluster", "cluster_utils.R"))

# source(here("R", "settings", "test.R"))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
  source(here("R", "settings", paste0(args[1], ".R")))
  n_jobs <- as.integer(args[2])
} else {
  stop("accept two paramters: expr_setting and maximal number of jobs.")
}

settings <- list()
results <- list()

# create data for each experiment to access
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
  
  settings[[iter_i]] <- set
  
  results[[iter_i]] <- list()
  for(sample_i in 1:sample_size){
    results[[iter_i]]$res <- list()
  }
}


# info used to get the index of the data record for each experiment
index_data <- list(n_expr = length(fig_var$value),
                   sample_len = sample_size, n_jobs = n_jobs)
trial_num <- length(fig_var$value) * sample_size

# prepare directories to store the data
unlink(here("data", "temp", experiment), recursive = TRUE) 
dir.create(here("data", "temp", experiment), showWarnings = F)
dir.create(here("data", "temp", experiment, "progress"), showWarnings = F)
save(settings, results, index_data,
     file = here("data", "temp", experiment, "settings.RData"))
save(trial_num, file = here("data", "temp", experiment, "trial_num.RData"))

unlink(here("log", experiment), recursive = TRUE) 
dir.create(here("log", experiment), showWarnings = F)


## create bash file
# read template
bash_temp <- here("jobs", "expr_template.sh")
con <- file(bash_temp, open="r")
expr_bash <- readLines(con)
close(con)

# specify file content
expr_bash[3] <- paste0("#SBATCH --job-name=", experiment)
expr_bash[4] <- paste0("#SBATCH --output=../log/", experiment, "/%a.out")
expr_bash[5] <- paste0("#SBATCH --array=1-", n_jobs)
expr_bash[14] <- paste0("EXPERIMENT=", experiment)

# write to bash file
expr_bash_file <- here("jobs", paste0(experiment, "_job.sh"))
con <- file(expr_bash_file, open="w")
writeLines(expr_bash, con)
close(con)

