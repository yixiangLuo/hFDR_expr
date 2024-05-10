library(here)
library(abind)

source(here("R", "utils.R"))
source(here("R", "cluster", "cluster_utils.R"))
# source(here("R", "plot.R"))

# read cmd args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
    experiment <- args[1]
} else {
    stop("accept one paramter: expr_name.")
}

# load packages and data
source(here("R", "settings", paste0(experiment, ".R")))
load(here("data", "temp", experiment, "settings.RData"))

# record computed result in the data structure
trial_num <- index_data$n_expr * index_data$sample_len
for(trial_id in 1:trial_num){
    trial_index <- get_expr_index(trial_id, index_data)
    
    load(here("data", "temp", experiment, paste0(trial_index$trial_id, ".RData")))
    results[[trial_index$expr_index]]$res[[trial_index$sample_index]] <- list(methods_hFDR = methods_hFDR, side_info = side_info, tune_seq = tune_seq)
}

expr_results <- lapply(1:index_data$n_expr, function(expr_iter){
    
    set <- settings[[expr_iter]]
    
    hFDR_res <- lapply(results[[expr_iter]]$res, function(result){
        result$hFDR_res
    })
    side_info <- lapply(results[[expr_iter]]$res, function(result){
        result$side_info
    })
    hFDR_res <- do.call("rbind", hFDR_res)
    
    tune_seq <- results[[expr_iter]]$res[[1]]$tune_seq
    
    result <- list(hFDR = hFDR_res, tune_seq = tune_seq, side_info = side_info)
    
    return(result)
})

results <- expr_results
names(results) <- char_to_digit(fig_var$value)

save(results, fig_var, file = here("data", paste0(experiment, ".RData")))


# method_names <- c("VR_pi0")  # "VR_pi0", "VR", "V_pi0", "V"
# draw_hFDR(experiment, fig_var, band_mass = 0.95,
#           method_names, multi_method_color, multi_method_shape)
# # draw_diff(experiment, fig_var, band_mass = 0.95,
# #           method_names, multi_method_color, multi_method_shape)
# 
# draw_dev(experiment, fig_var, band_mass = 0.95)
# draw_dev_pred(experiment, fig_var, lambda.index = 5)





