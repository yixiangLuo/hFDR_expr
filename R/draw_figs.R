#!/usr/bin/env Rscript

library(here)

source(here("R", "plot.R"))


aggregate <- list("fixedX-lasso" = c("IID_Normal", "X_AR", "Coef_AR", "Sparse"))
# aggregate <- list("fixedX-lasso-extend" = c("Weaker Signals", "Dense Signals", "Large Aspect Ratio", "Very Sparse"))
# aggregate <- list("modelX-lasso" = c("IID_Normal", "X_AR"),
#                   "modelX-logistic" = c("IID_Normal", "X_AR"))
# aggregate <- list("glasso" = c("Inv_Sparse"))
# aggregate <- list("fixedX-lasso" = c("MCC"))
# aggregate <- list("fixedX-lasso-hetero-robust" = c("Heteroscedasticity"),
#                   "fixedX-lasso-nonGauss-robust" = c("t-noise"),
#                   "modelX-logistic-robust" = c("Exponential"),
#                   "glasso-robust" = c("X Perturbation"))
# aggregate <- list("test" = c("IID_Normal"))

draw_aggregate(aggregate, band_mass = 0.9)

