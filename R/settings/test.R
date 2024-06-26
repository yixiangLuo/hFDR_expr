library(here)
library(hFDR)

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "test"
sel_method <- "lasso" # lasso, FS, glasso
model <- "gausslinear"

p_seq <- 100
n_seq <- 3*p_seq
pi1_seq <- 10 / p_seq

X_seed <- 2023
cov_seed <- -2023

# "PEC_0", "PEC_2", "PEC_4", "PEC_6", "PEC_8", "PEC_98"
X_types <- c("IID_Normal") # "IID_Normal", "MCC", "MCC_Block", "X_AR"
posit_types <- rep("random", length(X_types))  # equi
model_Xs <- X_types %in% c("IID_Normal", "X_AR")
random_Xs <- X_types %in% c("IID_Normal", "X_AR")
random_covs <- c(F)
scaling_Xs <- c(F)

X_mismodel <- F
y_mismodel <- F

target_at_alpha <- 0.2
# targets <- c(0.5)
targets <- c(0.8)
calib_methods <- c("lasso")

beta_permutes <- NA
noises <- c(quote(rnorm(n)))
# noises <- c(quote(rnorm(n, mean = 0, sd = (1 + rexp(n, 1))/2)))
# noises <- c(quote(rnorm(n, mean = 0, sd = runif(n, 1e-6, 2))))

# fig_var <- list(name = "n / m = ", value = n_seq / p_seq)
# fig_var <- list(name = "Signal strength: BH(0.05) has power", value = targets)
# fig_var <- list(name = "sparsity: pi1 = ", value = pi1_seq)
fig_var <- list(name = "Design Matrix", value = X_types)

makeup_vectors(p_seq = p_seq, n_seq = n_seq,
               X_types = X_types, posit_types = posit_types,
               model_Xs = model_Xs, random_Xs = random_Xs,
               random_covs = random_covs, scaling_Xs = scaling_Xs,
               pi1_seq = pi1_seq, targets = targets, calib_methods = calib_methods,
               beta_permutes = beta_permutes, noises = noises)

FDR_range <- c(0.01, 0.75)
n_tune <- 10
sample_size <- 100

n_cores <- 14


method_names <- c("hFDR", "std") # "std"
# "hFDR", "Poincare", "Bootstrap.ols", "Bootstrap.lasso_ols", "Bootstrap.lasso_ols.1se", "Bootstrap.nonpara", "Efron_Stein", "Efron_Stein.oracle", "Bootstrap.Bayes"

fig_titles <- paste0("fig: ", fig_var$value)

# method_colors <- unname(multi_method_color[method_names])
# method_shapes <- unname(multi_method_shape[method_names])










