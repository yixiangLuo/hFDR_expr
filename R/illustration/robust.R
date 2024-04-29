# load in useful packages
library(here)

library(profvis)
library(tictoc)

library(foreach)
library(doParallel)


source(here("R", "hFDR.R"))
source(here("R", "utils.R"))




# sufficient statistics

method <- "glasso"
df <- 3
show <- "pdf"

if(method == "logistic"){
  p <- 100
  n <- 3*p
  pi1 <- 10 / p
  mu1 <- 0.5627614
  X_type <- "X_AR"
  method <- "logistic"
} else if(method == "glasso"){
  p <- 15
  n <- 50*p
  pi1 <- 10 / p
  mu1 <- 6.282682
  X_type <- "Inv_Sparse"
  method <- "glasso"
} else stop()




set <- list(n = n, p = p, pi1 = pi1, X_type = X_type, random_X = T, scaling_X = F,
            X_mismodel = T, y_mismodel = F,
            model_X = T, posit_type = "random", sel_method = method, 
            calib_method = method, random_cov = T, X_seed = 2023, cov_seed = -2023)


# mu1 <- signal_calib(set, NA, nreps = 50, alpha = 0.2, target = 0.8, n_cores = 14)



gene_suff_samples <- function(deviate, sample_size){
  sapply(1:sample_size, function(iter){
    set.seed(10000 + iter)
    cov_seed <- -iter
    X.data <- gene_X(set$X_type, n, p, iter,
                     set$scaling_X, set$model_X,
                     signal = mu1, nonnulls = p*set$pi1, cov_seed = cov_seed,
                     mis_model = deviate)
    X <- X.data$X
    Xcov.true <- X.data$Xcov.true
    
    if(method == "logistic"){
      beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
      H0 <- beta == 0
      y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
      return(t(X) %*% y)
    } else if(method == "glasso"){
      Sigma <- t(X) %*% X
      lower_tri <- lower.tri(Sigma)
      return(c(Sigma[lower_tri]))
    } else stop()
  })
}

n_sample <- 100
results.auth <- gene_suff_samples(deviate = F, n_sample)
results.dev <- gene_suff_samples(deviate = T, n_sample)
variables <- 1:NROW(results.auth)

ecdf.auth <- ecdf(c(results.auth[variables, ]))
ecdf.dev <- ecdf(c(results.dev[variables, ]))
data <- c(results.auth[variables, ], results.dev[variables, ])
# x_range <- unname(quantile(data, c(0.01, 0.99)))
x_range <- c(min(data), max(data))
x <- seq(x_range[1], x_range[2], length.out = 100)

epdf.auth <- density(c(results.auth[variables, ]))
epdf.dev <- density(c(results.dev[variables, ]))

if(show == "cdf"){
  plot(x, ecdf.auth(x), xlim = x_range, main = "marginal", type = "n")
  lines(x, ecdf.auth(x), col = "blue")
  points(x, ecdf.dev(x), col = "red")
} else{
  plot(epdf.auth, xlim = x_range, main = "marginal", type = "n")
  lines(epdf.auth, col = "blue")
  points(epdf.dev, col = "red")
  
  plot(epdf.auth$x, log(epdf.auth$y), xlim = x_range, main = "marginal", type = "n")
  lines(epdf.auth$x, log(epdf.auth$y), col = "blue")
  points(epdf.dev$x, log(epdf.dev$y), col = "red")
}







# 
# variables <- 1
# n_sample <- 1
# n_mc <- 1000
# 
# 
# gene_cond_suff_samples <- function(deviate, sample_size = n_sample){
#   sapply(1:sample_size, function(iter){
#     set.seed(10000 + iter)
#     cov_seed <- -iter
#     X.data <- gene_X(set$X_type, n, p, iter,
#                      set$scaling_X, set$model_X,
#                      signal = mu1, nonnulls = p*set$pi1, cov_seed,
#                      mis_model = T)
#     X <- X.data$X
#     Xcov.true <- X.data$Xcov.true
# 
#     if(method == "logistic"){
#       beta <- genmu(p, set$pi1, mu1, set$posit_type, 1)
#       H0 <- beta == 0
#       y <- rbinom(n, 1, 1 / (1 + exp(-X %*% beta)))
# 
#       print("variables are null: ", variables %in% which(H0))
# 
#       samples <- sapply(variables, function(j){
#         if(deviate){
#           Sigma22_inv <- solve(Xcov.true[-j, -j])
#           mean.cond <- c(Xcov.true[j, -j] %*% Sigma22_inv %*% t(X[, -j]))
#           cov.cond <- c(Xcov.true[j, j] - Xcov.true[j, -j] %*% Sigma22_inv %*% Xcov.true[-j, j])
# 
#           Xj.samples <- sapply(1:n, function(row){
#             df_tilde <- df + p - 1
#             mu.cond <- mean.cond[row]
#             d2 <- c(t(X[row, -j])) %*% Sigma22_inv %*% X[row, -j]
#             c(mu.cond + rt.mult(n_mc, 1, df_tilde, (df+d2)/df_tilde * cov.cond))
#           })
#           Xj.samples <- t(Xj.samples)
#         } else{
#           Xj.samples <- sampler.modelX(X, Xcov.true, j, n_mc)
#         }
#         t(Xj.samples) %*% y
#       })
# 
#       return(c(samples))
#     } else if(method == "glasso"){
#       samples <- sapply(variables, function(var_i){
#         j <- sum(((2:p)-1)*((2:p)-2)/2 < var_i) + 1
#         i <- var_i - (j-1)*(j-2)/2
#         if(deviate){
#           stop()
#         } else{
#           stop()
#         }
#       })
# 
#       return(c(samples))
#     } else stop()
#   })
# }
# 
# results.auth <- gene_cond_suff_samples(deviate = F)
# results.dev <- gene_cond_suff_samples(deviate = T)
# 
# ecdf.auth <- ecdf(c(results.auth))
# ecdf.dev <- ecdf(c(results.dev))
# data <- c(results.auth, results.dev)
# x_range <- unname(quantile(data, c(0.05, 0.95)))
# x <- seq(x_range[1], x_range[2], length.out = 100)
# 
# if(show == "cdf"){
#   plot(x, ecdf.auth(x), xlim = x_range, main = "conditional", type = "n")
#   lines(x, ecdf.auth(x), col = "blue")
#   points(x, ecdf.dev(x), col = "red")
# } else{
#   plot(density(c(results.auth)), xlim = x_range, main = "conditional", type = "n")
#   lines(density(c(results.auth)), col = "blue")
#   points(density(c(results.dev)), col = "red")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# sampler.core.t <- function(i, j, Xcov.true, sample_size){
#   precision.null <- solve(Xcov.true)
#   precision.null[i, j] <- 0
#   precision.null[j, i] <- 0
#   Xcov.true <- solve(precision.null)
# 
#   basis <- X[, -c(i,j)]
#   proj_mat <- basis %*% solve(t(basis) %*% basis) %*% t(basis)
# 
#   i.proj <- proj_mat %*% X[, i]
#   j.proj <- proj_mat %*% X[, j]
# 
#   i.radius <- sqrt(sum(X[, i]^2) - sum(i.proj^2))
#   j.radius <- sqrt(sum(X[, j]^2) - sum(j.proj^2))
# 
#   var.sample <- matrix(rnorm(n * sample_size), nrow = n)
# 
#   var.sample.proj <-lm(var.sample ~ basis + 0)$fitted.values
# 
#   var.sample <- var.sample - var.sample.proj
#   var.sample <- scale(var.sample, center = FALSE, scale = sqrt(colSums(var.sample^2)) / radius)
#   var.sample <- var.proj + var.sample
# 
#   subspace_samples <- is.na(colSums(var.sample))
#   if(sum(subspace_samples) > sample_size/2) stop()
#   for(bad_i in which(subspace_samples)){
#     var.sample[, bad_i] <- var.sample[, sample(which(!subspace_samples), 1)]
#   }
# 
#   return(var.sample)
# }