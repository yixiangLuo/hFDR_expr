#!/usr/bin/env Rscript

library(here)

source(here("R", "plot.R"))


# aggregate <- list("fixedX-lasso" = c("IID_Normal", "X_AR", "Coef_AR", "Sparse"))
# aggregate <- list("fixedX-lasso-extend" = c("Weaker Signals", "Dense Signals", "Large Aspect Ratio", "Very Sparse"))
# aggregate <- list("modelX-lasso" = c("IID_Normal", "X_AR"),
#                   "modelX-logistic" = c("IID_Normal", "X_AR"))
# aggregate <- list("glasso" = c("Inv_Sparse"))
# aggregate <- list("fixedX-lasso" = c("MCC"))
aggregate <- list("fixedX-lasso-hetero-robust" = c("Heteroscedasticity"),
                  "fixedX-lasso-nonGauss-robust" = c("t-noise"),
                  "modelX-logistic-robust" = c("Exponential"),
                  "glasso-robust" = c("X Perturbation"))

draw_aggregate(aggregate, band_mass = 0.9)


# draw_hFDR_aggre(experiment = "fixedX-", sel_methods = c("lasso", "FS"), X_types = c("IID_Normal", "X_AR"),
#                 band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_hFDR_aggre(experiment = "fixedX-", sel_methods = c("lasso", "FS"), X_types = c("MCC"),
#                 band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_hFDR_aggre(experiment = "modelX-", sel_methods = c("lasso", "logistic"), X_types = c("IID_Normal", "X_AR"),
#                 band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_hFDR_aggre(experiment = "", sel_methods = c("glasso"), X_types = c("Inv_Sparse"),
#                 band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# 
# 
# draw_dev_aggre(experiment = "fixedX-", sel_methods = c("lasso", "FS"), X_types = c("IID_Normal", "X_AR"),
#               band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_dev_aggre(experiment = "fixedX-", sel_methods = c("lasso", "FS"), X_types = c("MCC"),
#               band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_dev_aggre(experiment = "modelX-", sel_methods = c("lasso", "logistic"), X_types = c("IID_Normal", "X_AR"),
#               band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
# draw_dev_aggre(experiment = "", sel_methods = c("glasso"), X_types = c("Inv_Sparse"),
#               band_mass = 0.95, method_names, multi_method_color, multi_method_shape)
