library(Matrix)
library(here)
source(here("R", "utils.R"))



lasso_homotopy_step <- function(lasso.pack, j, Xjy_range, direction){
  
  Xy <- lasso.pack$Xy
  beta_lasso <- lasso.pack$beta_lasso
  selected <- lasso.pack$selected
  Sigma <- lasso.pack$Sigma
  Sigma.selected.inv <- lasso.pack$Sigma.selected.inv
  lambda <- lasso.pack$lambda
  tol <- lasso.pack$tol
  
  beta_lasso <- as.vector(beta_lasso)
  if(!is.null(Sigma.selected.inv)){
    Sigma.selected.inv <- matrix(Sigma.selected.inv, nrow = NROW(Sigma.selected.inv))
  }
  direction <- sign(direction)
  p <- NCOL(Sigma)
  n_sel <- length(selected)
  not_selected <- which(!(1:p %in% selected))
  
  if(n_sel > 0){
    j_index <- which(selected == j)
    slope.beta <- if(length(j_index) > 0){
      c(Sigma.selected.inv[, j_index])
    } else{
      rep(0, n_sel)
    }
    intercept.beta <- beta_lasso[selected]
    PoC.beta <- (0 - intercept.beta) / slope.beta
    PoC.beta[is.na(PoC.beta)] <- 0
    # PoC.beta[abs(PoC.beta) < tol] <- direction * tol
    PoC.beta <- dir_trunc(PoC.beta, direction)
  } else{
    PoC.beta <- NULL
  }
  
  if(n_sel < p){
    slope.lambda <- (not_selected == j) - 
      ifelse(n_sel > 0, matrix(Sigma[not_selected, selected], ncol = n_sel) 
             %*% slope.beta, 0)
    intercept.lambda <- Xy[not_selected] - Sigma[not_selected, ] %*% beta_lasso
    intercept.lambda <- sign(intercept.lambda) * pmin(abs(intercept.lambda), lambda*(1-tol))
    PoC.lambda <- (sign(slope.lambda) * direction * lambda - intercept.lambda) / slope.lambda
    PoC.lambda[is.na(PoC.lambda)] <- 0
    # PoC.lambda[abs(PoC.lambda) < tol] <- direction * tol
    PoC.lambda <- dir_trunc(PoC.lambda, direction)
  } else{
    PoC.lambda <- NULL
  }
  
  PoC <- rep(NA, p)
  PoC[c(selected, not_selected)] <- c(PoC.beta, PoC.lambda)
  PoC_index <- which.min(abs(PoC))
  PoC <- PoC[PoC_index]
  
  if(abs(PoC) < Inf){
    exclude <- which(selected %in% PoC_index)
    next.selected <- if(length(exclude) > 0){
      # sort(setdiff(selected, PoC_index))
      selected[-exclude]
    } else{
      insert_index <- sum(selected < PoC_index)
      # sort(union(selected, PoC_index))
      append(selected, PoC_index, after = insert_index)
    }
    next.Xjy <- Xy[j] + PoC
  } else{
    next.selected <- selected
    next.Xjy <- ifelse(PoC > 0, max(Xjy_range), min(Xjy_range))
    PoC <- next.Xjy - Xy[j]
  }
  next.beta_lasso <- beta_lasso
  if(n_sel > 0){
    next.beta_lasso[selected] <- beta_lasso[selected] + slope.beta * PoC
    next.beta_lasso[abs(next.beta_lasso) <= tol] <- 0
  }
  if(any(is.na(next.beta_lasso))) browser()
  # if(abs(PoC) < 1e-7) browser()
  
  stop <- if(direction > 0){
    next.Xjy >= max(Xjy_range)
  } else{
    next.Xjy <= min(Xjy_range)
  }
  
  X <- lasso.pack$X
  sen_coef <- if(n_sel == 0 || j %in% selected){
    rep(0, NROW(X))
  } else{
    c(t(X[, j]) %*% X[, selected] %*% Sigma.selected.inv %*% t(X[, selected]))
  }
  
  return(list(Xjy = next.Xjy, beta_lasso = next.beta_lasso,
              selected = next.selected, stop = stop, sen_coef = sen_coef))
}

update_selectedInv <- function(Sigma, Sigma.selected.inv, selected, next.selected){
  if(length(next.selected) == 0){
    return(NULL)
  }
  if(length(selected) == 0){
    return(solve(Sigma[next.selected, next.selected]))
  }
  
  include <- next.selected[!(next.selected %in% selected)]
    # setdiff(next.selected, selected)
  exclude <- selected[!(selected %in% next.selected)]
    # setdiff(selected, next.selected)
  
  if(length(include) + length(exclude) != 1){
    stop()
  }
  if(length(include) > 0){
    l <- length(selected)
    update.index <- order(c(selected, include))
    ## inverse of [A, b; b^T, d] by Woodbury formula
    # original = [A, 0; 0, d]
    # org.inv <- bdiag(Sigma.selected.inv, 1/Sigma[include, include])
    org.inv <- cbind(rbind(Sigma.selected.inv, rep(0, l)), c(rep(0, l), 1/Sigma[include, include]))
    # U = [0, b; 1, 0], V = [0, b^T; 1, 0]
    # U <- sparseMatrix(i = 1:(1+l), j = c(rep(2, l), 1),
    #                   x = c(Sigma[selected, include], 1))
    # V <- sparseMatrix(i = c(rep(1, l), 2), j = 1:(1+l),
    #                   x = c(Sigma[selected, include], 1))
    U <- matrix(c(rep(0, l), 1, Sigma[selected, include], 0), ncol = 2)
    # V <- t(matrix(c(Sigma[selected, include], 0, rep(0, l), 1), ncol = 2))
    V <- rbind(c(Sigma[selected, include], 0), c(rep(0, l), 1))
    org.inv_U <- org.inv %*% U
    update.inv <- org.inv - org.inv_U %*% solve_mat22(diag(2) + V %*% org.inv_U) %*% (V %*% org.inv)
    update.inv <- update.inv[update.index, update.index]
  } else if(length(exclude) > 0){
    l <- length(next.selected)
    update.index <- rank(c(next.selected, exclude))
    ## inverse of [A, 0; 0, d] by Woodbury formula
    # original = [A, b; b^T, d]
    org.inv <- Sigma.selected.inv[update.index, update.index]
    # U = [0, b; 1, 0], V = [0, b^T; 1, 0]
    # U <- sparseMatrix(i = 1:(1+l), j = c(rep(2, l), 1),
    #                   x = c(-Sigma[next.selected, exclude], 1))
    # V <- sparseMatrix(i = c(rep(1, l), 2), j = 1:(1+l),
    #                   x = c(-Sigma[next.selected, exclude], 1))
    U <- matrix(c(rep(0, l), 1, -Sigma[next.selected, exclude], 0), ncol = 2)
    # V <- t(matrix(c(-Sigma[next.selected, exclude], 0, rep(0, l), 1), ncol = 2))
    V <- rbind(c(-Sigma[next.selected, exclude], 0), c(rep(0, l), 1))
    org.inv_U <- org.inv %*% U
    update.inv <- org.inv - org.inv_U %*% solve_mat22(diag(2) + V %*% org.inv_U) %*% (V %*% org.inv)
    update.inv <- update.inv[1:l, 1:l]
  }
  
  
  # org <- solve(diag(2) + V %*% org.inv_U)
  # mine <- solve_mat22(diag(2) + V %*% org.inv_U)
  # if(max(abs(org - mine)) > 1e-2) browser()
  
  # next.Sigma.selected.inv <- solve(Sigma[next.selected, next.selected])
  
  return(update.inv)
}

lasso_homotopy <- function(lasso.pack, j,
                           Xjy_range = c(-max(abs(Xy)), max(abs(Xy)))){
  
  Sigma <- lasso.pack$Sigma
  
  p <- length(lasso.pack$Xy)
  lasso.pack$beta_lasso <- as.vector(lasso.pack$beta_lasso)
  
  selected <- lasso.pack$selected
  not_selected <- setdiff(1:p, selected)
  if(length(selected) == 0){
    lasso.pack$Sigma.selected.inv <- NULL
  } else if(is.null(lasso.pack$Sigma.selected.inv)){
    lasso.pack$Sigma.selected.inv <- solve(Sigma[selected, selected])
  }
    
  Xjy_nodes <- c()
  beta_at_nodes <- c()
  selected.list <- list()
  Sigma.selected.inv.list <- list()
  sen_coef <- 0
  
  for(direction in c(1, -1)){
    stop <- F
    lasso.pack.cur <- lasso.pack
    
    while (!stop) {
      if(direction < 0){
        selected.list <- c(selected.list, list(lasso.pack.cur$selected))
        Sigma.selected.inv.list <- c(Sigma.selected.inv.list, list(lasso.pack.cur$Sigma.selected.inv))
      }
      
      res_at_node <- lasso_homotopy_step(lasso.pack.cur, j, Xjy_range, direction)
      Xjy_nodes <- c(Xjy_nodes, res_at_node$Xjy)
      beta_at_nodes <- cbind(beta_at_nodes, res_at_node$beta_lasso)
      if(sum(abs(res_at_node$sen_coef)) > 0)
        sen_coef <- res_at_node$sen_coef
      
      stop <- res_at_node$stop
      if(!stop){
        lasso.pack.cur$X <- lasso.pack$X
        lasso.pack.cur$Sigma.selected.inv <- update_selectedInv(Sigma,
                                                                lasso.pack.cur$Sigma.selected.inv,
                                                                lasso.pack.cur$selected,
                                                                res_at_node$selected)
        lasso.pack.cur$Xy[j] <- res_at_node$Xjy
        lasso.pack.cur$beta_lasso <- res_at_node$beta_lasso
        lasso.pack.cur$selected <- res_at_node$selected
      }
      
      if(direction > 0){
        selected.list <- c(selected.list, list(lasso.pack.cur$selected))
        Sigma.selected.inv.list <- c(Sigma.selected.inv.list, list(lasso.pack.cur$Sigma.selected.inv))
      }
    }
  }
  
  order_index <- order(Xjy_nodes)
  Xjy_nodes <- Xjy_nodes[order_index]
  beta_at_nodes <- beta_at_nodes[, order_index]
  selected.list <- selected.list[order_index]
  Sigma.selected.inv.list <- Sigma.selected.inv.list[order_index]
  
  lasso.homopath <- structure(list(Xjy_nodes = Xjy_nodes,
                                   beta_at_nodes = beta_at_nodes,
                                   selected_at_nodes_right = selected.list,
                                   matrix_at_nodes_right = Sigma.selected.inv.list,
                                   sen_coef = sen_coef),
                              class = 'lasso.homopath')
  
  return(lasso.homopath)
}

homotopy_fit <- function(Xjy, lasso.homopath){
  Xjy_nodes <- lasso.homopath$Xjy_nodes
  beta_at_nodes <- lasso.homopath$beta_at_nodes
  if(Xjy > max(Xjy_nodes) || Xjy < min(Xjy_nodes)){
    stop("provided homotopy path doesn't cover the desired Xjy.")
  }
  left <- min(max(which(Xjy >= Xjy_nodes)), length(Xjy_nodes)-1)
  right <- left + 1
  beta_lasso <- ((Xjy - Xjy_nodes[left]) * beta_at_nodes[, right] + 
    (Xjy_nodes[right] - Xjy) * beta_at_nodes[, left]) / (Xjy_nodes[right] - Xjy_nodes[left])
  
  return(list(beta_lasso = beta_lasso,
              selected = lasso.homopath$selected_at_nodes_right[[left]],
              Sigma.selected.inv = lasso.homopath$matrix_at_nodes_right[[left]]))
}