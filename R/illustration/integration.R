# load in useful packages
library(here)
library(tidyverse)


source(here("R", "hFDR.R"))
source(here("R", "utils.R"))



# investigate the conditional distribution (fixed-X lasso)
cond_distribution <- function(X, y, lambda, j, n_samples = 100, n_cores = 14){
  p <- NCOL(X)
  n <- NROW(X)
  
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2)/n))
  
  tvals <- lm_to_t(y, X)
  pvals <- pvals_t(tvals$tvals, tvals$df, side = "two")
  print(paste0("p-value: ", pvals[j]))
  
  lasso.fit <- glmnet::glmnet(X, y, lambda = lambda,
                              intercept = T, standardize = T,
                              family = "gaussian")
  
  glmnet.pack <- pack_glmnet(X, y, lasso.fit)
  data.pack <- process_data(X, y)
  
  lambda_i <- 1
  lasso.pack <- list(X = X, Xy = glmnet.pack$Xy,
                     beta_lasso = glmnet.pack$betas[, lambda_i],
                     selected = glmnet.pack$selected.list[[lambda_i]],
                     Sigma = glmnet.pack$Sigma,
                     Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
                     lambda = glmnet.pack$lambdas[lambda_i],
                     tol = 1e-8)
  Xjy_range <- c(data.pack$Xy_bound[, j])
  
  lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
  Xjy <- lasso.homopath$Xjy_nodes
  vjy <- Xjy_to_vjy(Xjy, data.pack$trans, j)
  tval.nodes <- vjy / sqrt(data.pack$RSS_X / (n-p))
  lasso.beta.j <- lasso.homopath$beta_at_nodes[j, ]
  
  interp <- approx(tval.nodes, lasso.beta.j, n = n_samples)
  tval <- interp$x
  lasso.beta.j <- interp$y
  
  pval <- pvals_t(tval, n-p, side = "two")
  integrand <- sapply(tval, function(t){
    ind <- max(max(which(t > tval.nodes)), 1)
    selected <- unname(lasso.homopath$selected_at_nodes_right[[ind]])
    (j %in% selected) / max(1, length(selected))
  })
  results <- list(tval = tval, pval = pval, lasso.beta.j = lasso.beta.j, integrand = integrand)
  
  # # by monte carlo, sometimes noisy due to glmnet, can be used for sanity check
  # parallel <- process_parallel(n_cores)
  # forall <- parallel$iterator
  # `%exec%` <- parallel$connector
  # 
  # y.samples <- sampler.fixedX(X, y, j, n_samples)
  # results <- forall(mc_i = 1:n_samples, .options.multicore = list(preschedule = F)) %exec% {
  #   X.sample <- X
  #   y.sample <- y.samples[, mc_i]
  # 
  #   tvals <- lm_to_t(y.sample, X.sample)
  #   tval <- tvals$tvals[j]
  #   pval <- pvals_t(tvals$tvals, tvals$df, side = "two")[j]
  # 
  #   lasso.fit <- glmnet::glmnet(X.sample, y.sample, lambda = lambda,
  #                               intercept = T, standardize = T,
  #                               family = "gaussian")
  #   selected <- which(abs(lasso.fit$beta) > 1e-8)
  # 
  #   integrand <- (j %in% selected) / max(1, length(selected))
  # 
  #   lasso.beta <- c(as.matrix(lasso.fit$beta))
  #   lasso.beta.j <- lasso.beta[j]
  # 
  #   res <- c(tval, pval, lasso.beta.j, integrand)  # t_j, p_j, \hbeta_j, 1{j \in R}/R
  #   names(res) <- c("tval", "pval", "lasso.beta.j", "integrand")
  #   return(res)
  # }
  # results <- as.data.frame(t(do.call(cbind, results))) %>% arrange(tval)
}


# example usage
seed <- 1
set.seed(seed)
p <- 100
n <- 3*p
pi1 <- 10 / p
mu1 <- 0.5
X <- gene_X("X_AR_Strong", n, p, X_seed = seed, scaling = F, model_X = T)$X 
beta <- genmu(p, pi1, mu1, "random", 1) * (1 + rexp(p, 1))/2
H0 <- beta == 0

y <- X %*% beta + rnorm(n)
lambda <- 10 / n

lambda <- 30 / n
j <- 1

# n_sample is the number points in the plot
res <- cond_distribution(X, y, lambda, j, n_samples = 1000, n_cores = 14)

# due to numerical error of glmnet, some figures may be noisy

pdf(here("figs","showcalc.pdf"),width=7,height=4)
zeta <- 0.1
thresh <- qt(zeta/2,df=n-p,lower.tail=FALSE)
plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
     type = "n",main="Estimating the FDR Contribution",xlab="t-statistic",ylab="",
     xlim=c(-4,5))
#plot(res$tval, res$lasso.beta.j, type = "n",)
polygon(c(res$tval,rev(res$tval)),c(dt(res$tval,df=n-p),rep(0,length(res$tval))),
        border="blue")
abline(h=0,col="darkgray")
abline(v=0,col="darkgray")
#abline(v=c(-thresh,thresh),col="blue",lty=2)
final.integrand <- dt(res$tval,df=n-p)*res$integrand
which.plot <- which(final.integrand>0)
polygon(res$tval[c(which.plot,rev(which.plot))],
        c(final.integrand[which.plot],rep(0,length(which.plot))),
        col = "purple",border="purple")
lines(res$tval, res$lasso.beta.j, col = "black")
lines(res$tval,dt(res$tval,df=n-p),col="blue")
lines(res$tval, res$integrand, col = "red")
legend("topleft",col=c("black","red","blue","purple"),
       lty=1,lwd=c(1,1,1,10),bty="n",
       legend=c(expression(hat(beta)[j]),
                TeX("$1\\{j \\in R\\}/R$"),
                TeX("$t_{n-d} \\;density$"),
                TeX("$\\widehat{FDR}_j$")))
dev.off()


pdf(here("figs","talk","showcalc1.pdf"),width=7,height=4)
plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
     type = "n",main="Estimating the FDR Contribution",xlab="t-statistic",ylab="",
     xlim=c(-4,5))
abline(h=0,col="darkgray")
abline(v=0,col="darkgray")
lines(res$tval, res$lasso.beta.j, col = "black")
legend("topleft",col=c("black"),
       lty=1,lwd=c(1),bty="n",
       legend=c(expression(hat(beta)[j])))
dev.off()

pdf(here("figs","talk","showcalc2.pdf"),width=7,height=4)
plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
     type = "n",main="Estimating the FDR Contribution",xlab="t-statistic",ylab="",
     xlim=c(-4,5))
abline(h=0,col="darkgray")
abline(v=0,col="darkgray")
lines(res$tval, res$lasso.beta.j, col = "black")
lines(res$tval, res$integrand, col = "red")
legend("topleft",col=c("black","red"),
       lty=1,lwd=c(1,1),bty="n",
       legend=c(expression(hat(beta)[j]),
                TeX("$1\\{j \\in R\\}/R$")))
dev.off()

pdf(here("figs","talk","showcalc3.pdf"),width=7,height=4)
plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
     type = "n",main="Estimating the FDR Contribution",xlab="t-statistic",ylab="",
     xlim=c(-4,5))
abline(h=0,col="darkgray")
abline(v=0,col="darkgray")
lines(res$tval, res$lasso.beta.j, col = "black")
lines(res$tval, res$integrand, col = "red")
lines(res$tval,dt(res$tval,df=n-p),col="blue")
legend("topleft",col=c("black","red","blue"),
       lty=1,lwd=c(1,1,1),bty="n",
       legend=c(expression(hat(beta)[j]),
                TeX("$1\\{j \\in R\\}/R$"),
                TeX("$t_{n-d} \\;density$")))
dev.off()

pdf(here("figs","talk","showcalc4.pdf"),width=7,height=4)
plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
     type = "n",main="Estimating the FDR Contribution",xlab="t-statistic",ylab="",
     xlim=c(-4,5))
abline(h=0,col="darkgray")
abline(v=0,col="darkgray")
lines(res$tval, res$lasso.beta.j, col = "black")
lines(res$tval, res$integrand, col = "red")
lines(res$tval,dt(res$tval,df=n-p),col="blue")
final.integrand <- dt(res$tval,df=n-p)*res$integrand
which.plot <- which(final.integrand>0)
polygon(res$tval[c(which.plot,rev(which.plot))],
        c(final.integrand[which.plot],rep(0,length(which.plot))),
        col = "purple",border="purple")
legend("topleft",col=c("black","red","blue","purple"),
       lty=1,lwd=c(1,1,1,10),bty="n",
       legend=c(expression(hat(beta)[j]),
                TeX("$1\\{j \\in R\\}/R$"),
                TeX("$t_{n-d} \\;density$"),
                TeX("$\\widehat{FDR}_j$")))
dev.off()

for (stage in 1:4) {
  pdf(here("figs","talk",paste0("showcalc",stage,".pdf")),width=5,height=3)
  par(mar=c(4.5,2.5,0.1,0.1))
  plot(c(res$tval,res$tval), c(res$lasso.beta.j,dt(res$tval,df=n-p)), 
       type = "n",#main="Estimating the FDR Contribution",
       xlab="t-statistic",ylab="",
       xlim=c(-4,5))
  abline(h=0,col="darkgray")
  abline(v=0,col="darkgray")
  lines(res$tval, res$lasso.beta.j, col = "black")
  if(stage >= 2) {
    lines(res$tval, res$integrand, col = "red")
  }
  if(stage >= 3) {
    lines(res$tval,dt(res$tval,df=n-p),col="blue")
  }
  if(stage >= 4) {
    final.integrand <- dt(res$tval,df=n-p)*res$integrand
    which.plot <- which(final.integrand>0)
    polygon(res$tval[c(which.plot,rev(which.plot))],
            c(final.integrand[which.plot],rep(0,length(which.plot))),
            col = "purple",border="purple")
  }
  legend("topleft",lty=1,bty="n",
         col=c("black","red","blue","purple")[1:stage],
         lwd=c(1,1,1,10)[1:stage],
         legend=c(expression(hat(theta)[j]),
                  TeX("$1\\{j \\in R\\}/R$"),
                  TeX("$t_{n-d} \\;density$"),
                  TeX("$\\widehat{FDR}_j$"))[1:stage])
  dev.off()
}
