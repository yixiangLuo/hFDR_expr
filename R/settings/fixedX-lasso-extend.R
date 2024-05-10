library(here)
library(hFDR)

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "fixedX-lasso-extend"
sel_method <- "lasso"
model <- "gausslinear"

p_seq <- c(500, 500, 60, 500)
n_seq <- c(1500, 1500, 180, 600)
pi1_seq <- c(30, 10, 30, 30) / p_seq

X_seed <- 2023
cov_seed <- -2023

X_types <- c("IID_Normal")
posit_types <- rep("random", length(X_types))
model_Xs <- !(X_types %in% c("MCC"))
random_Xs <- !(X_types %in% c("MCC"))
random_covs <- c(F)
scaling_Xs <- c(F)

X_mismodel <- F
y_mismodel <- F

target_at_alpha <- 0.2
targets <- c(0.3, 0.8, 0.8, 0.8)
calib_methods <- c("lasso")

noises <- c(quote(rnorm(n)))

fig_var <- list(name = "subfig", value = c("Weaker Signals", "Very Sparse", "Dense Signals", "Large Aspect Ratio"))

makeup_vectors(p_seq = p_seq, n_seq = n_seq,
               X_types = X_types, posit_types = posit_types,
               model_Xs = model_Xs, random_Xs = random_Xs,
               random_covs = random_covs, scaling_Xs = scaling_Xs,
               pi1_seq = pi1_seq, targets = targets, calib_methods = calib_methods,
               noises = noises)

FDR_range <- c(-0.01, 0.75)
n_tune <- 10
sample_size <- 200

n_cores <- 14


method_names <- c("hFDR", "std")

fig_titles <- paste0("fig: ", fig_var$value)











