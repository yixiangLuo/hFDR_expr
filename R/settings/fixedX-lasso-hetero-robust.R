library(here)

source(here("R", "homotopy_lasso.R"))
source(here("R", "lasso_hFDR.R"))
source(here("R", "hFDR.R"))

library(glmnet)   # lasso
# library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "fixedX-lasso-hetero-robust"
sel_method <- "lasso" # lasso, FS, glasso
model <- "fixed-X"

p_seq <- 500
n_seq <- 3*p_seq
pi1_seq <- 30 / p_seq

X_seed <- 2023
cov_seed <- -2023

X_types <- c("IID_Normal") # "IID_Normal", "MCC", "MCC_Block", "X_AR"
posit_types <- rep("random", length(X_types))
model_Xs <- !(X_types %in% c("MCC"))
random_Xs <- !(X_types %in% c("MCC"))
random_covs <- c(F)
scaling_Xs <- c(F)

X_mismodel <- F
y_mismodel <- T

target_at_alpha <- 0.2
# targets <- c(0.5)
targets <- c(0.8)
calib_methods <- c("lasso")

beta_permutes <- NA
noises <- c(quote(rnorm(n)))

# fig_var <- list(name = "n / m = ", value = n_seq / p_seq)
# fig_var <- list(name = "Signal strength: BH(0.05) has power", value = targets)
# fig_var <- list(name = "sparsity: pi1 = ", value = pi1_seq)
fig_var <- list(name = "robust", value = "Heteroscedasticity")

makeup_vectors(p_seq = p_seq, n_seq = n_seq,
               X_types = X_types, posit_types = posit_types,
               model_Xs = model_Xs, random_Xs = random_Xs,
               random_covs = random_covs, scaling_Xs = scaling_Xs,
               pi1_seq = pi1_seq, targets = targets, calib_methods = calib_methods,
               beta_permutes = beta_permutes, noises = noises)

FDR_range <- c(0.01, 0.75)
n_tune <- 10
sample_size <- 200

n_cores <- 14


get_method_list <- get_multi_method_list
method_names <- c("hFDR", "std")

fig_titles <- paste0("fig: ", fig_var$value)

# method_colors <- unname(multi_method_color[method_names])
# method_shapes <- unname(multi_method_shape[method_names])










