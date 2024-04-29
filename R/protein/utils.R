library(dagitty)


protein_names <- c("Raf", "Mek", "Plcg", "PIP2", "PIP3", "Erk", "Akt", "PKA",
                   "PKC", "P38", "Jnk")

source_dag <- "Sachs" # "Friedman"
if(source_dag == "Sachs"){
  # read from Figure 3(A) in https://www.esalq.usp.br/lepse/imgs/conteudo_thumb/Causal-Protein-Signaling-Networks-Derived-from-Multiparameter-Single-Cell-Data.pdf
  DAG <- dagitty("dag {
    Raf
    Mek
    Plcg
    PIP2
    PIP3
    Erk
    Akt
    PKA
    PKC
    P38
    Jnk
    Raf -> Mek
    PKA -> Raf
    PKC -> Raf
    Mek -> Erk
    PKA -> Mek
    PKC -> Mek
    Plcg -> PIP2
    Plcg -> PIP3
    Plcg -> PKC
    PIP3 -> PIP2
    PIP2 -> PKC
    PIP3 -> Akt
    Erk -> Akt
    PKA -> Erk
    PKA -> Akt
    PKA -> P38
    PKA -> Jnk
    PKC -> P38
    PKC -> Jnk
    PKC -> PKA
  }")
  # DAG edges
  pairs.DAG <- list(c("Raf", "Mek"), c("PKA", "Raf"), c("PKC", "Raf"),
                    c("Mek", "Erk"), c("PKA", "Mek"), c("PKC", "Mek"),
                    c("Plcg", "PIP2"), c("Plcg", "PIP3"), c("Plcg", "PKC"),
                    c("PIP3", "PIP2"), c("PIP2", "PKC"),
                    c("PIP3", "Akt"),
                    c("Erk", "Akt"), c("PKA", "Erk"),
                    c("PKA", "Akt"),
                    c("PKA", "P38"), c("PKA", "Jnk"),
                    c("PKC", "P38"), c("PKC", "Jnk"), c("PKC", "PKA"))
} else{
  # read from Figure 2 in https://academic.oup.com/biostatistics/article-pdf/9/3/432/17742149/kxm045.pdf
  # miss the PKC -> PKA edge in Sachs
  DAG <- dagitty("dag {
    Raf
    Mek
    Plcg
    PIP2
    PIP3
    Erk
    Akt
    PKA
    PKC
    P38
    Jnk
    Raf -> Mek
    PKA -> Raf
    PKC -> Raf
    Mek -> Erk
    PKA -> Mek
    PKC -> Mek
    Plcg -> PIP2
    Plcg -> PIP3
    Plcg -> PKC
    PIP3 -> PIP2
    PIP2 -> PKC
    PIP3 -> Akt
    Erk -> Akt
    PKA -> Erk
    PKA -> Akt
    PKA -> P38
    PKA -> Jnk
    PKC -> P38
    PKC -> Jnk
  }")
  # DAG edges
  pairs.DAG <- list(c("Raf", "Mek"), c("PKA", "Raf"), c("PKC", "Raf"),
                    c("Mek", "Erk"), c("PKA", "Mek"), c("PKC", "Mek"),
                    c("Plcg", "PIP2"), c("Plcg", "PIP3"), c("Plcg", "PKC"),
                    c("PIP3", "PIP2"), c("PIP2", "PKC"),
                    c("PIP3", "Akt"),
                    c("Erk", "Akt"), c("PKA", "Erk"),
                    c("PKA", "Akt"),
                    c("PKA", "P38"), c("PKA", "Jnk"),
                    c("PKC", "P38"), c("PKC", "Jnk"))
}

# Moralized graph: graph for conditional independence: 2(Sachs) / 3(Friedman) more edges, same moralized undirected graph
moral.DAG <- moralize(DAG)
# print moral.DAG and hard code its edges below
pairs.true <- list(c("Akt", "Erk"), c("Akt", "PIP3"), c("Akt", "PKA"),
                   c("Erk", "Mek"), c("Erk", "PIP3"), c("Erk", "PKA"),
                   c("Jnk", "PKA"), c("Jnk", "PKC"),
                   c("Mek", "PKA"), c("Mek", "PKC"), c("Mek", "Raf"),
                   c("P38", "PKA"), c("P38", "PKC"),
                   c("PIP2", "PIP3"), c("PIP2", "PKC"), c("PIP2", "Plcg"),
                   c("PIP3", "PKA"), c("PIP3", "Plcg"),
                   c("PKA", "PKC"), c("PKA", "Raf"),
                   c("PKC", "Plcg"), c("PKC", "Raf"))


adjmat.DAG <- matrix(F, 11, 11)
for(pair in pairs.DAG){
  i <- which(protein_names == pair[1])
  j <- which(protein_names == pair[2])
  adjmat.DAG[i, j] <- T
}

adjmat.true <- matrix(F, 11, 11)
for(pair in pairs.true){
  i <- which(protein_names == pair[1])
  j <- which(protein_names == pair[2])
  adjmat.true[i, j] <- T
  adjmat.true[j, i] <- T
}
null_pairs <- which(adjmat.true[upper.tri(adjmat.true)] == 0)

read_data <- function(exprs = 1:9){
  files <- list.files(path = here("data", "protein"), pattern = "*.xls")
  
  data <- lapply(files, function(file){
    data <- read_excel(here("data", "protein", file))
    names(data) <- protein_names
    return(data)
  })
  
  X <- do.call(rbind, data[exprs])
  X <- scale(X, center = T, scale = F)
  X <- scale(X, center = F, scale = sqrt(colSums(X^2) / (NROW(X) - 1)))
  
  n <- NROW(X)
  X <- X[sample(n, n), ]
  
  return(X)
}

gene_rho <- function(X, n_tune = 100){
  Xcov <- cov(X)
  rho <- mean(abs(Xcov))
  p <- NCOL(X)
  
  rho_max <- rho
  rho_min <- rho
  y <- NA
  sel_method <- "glasso"
  for(i in 1:1000){
    if(sum(select_variables(X, y, rho_max, sel_method)) > 1){
      rho_max <- 1.1 * rho_max
    } else break
  }
  for(i in 1:1000){
    if(sum(select_variables(X, y, rho_min, sel_method)) < p*(p-1)/2){
      rho_min <- rho_min / 1.1
    } else break
  }
  tune_seq <- exp(seq(log(rho_max), log(rho_min), length.out = n_tune))
  
  return(tune_seq)
}


plot.graph <- function(Theta, mode = "undirected"){
  adj.graph <- graph_from_adjacency_matrix(Theta != 0, mode = mode, diag = F)
  plot.igraph(adj.graph, size = 25, vertex.size = 15, 
              vertex.label = protein_names, vertex.label.dist = 6, vertex.label.cex = 2, 
              vertex.label.degree = -pi/2 - seq(0, 2*pi, length.out = 11),
              layout = function(graph){
                layout_in_circle(graph) %*% matrix(c(0, -1, 1, 0), 2)
              })
}


# compute FDP
FDP <- function(selected){
  length(intersect(selected, null_pairs)) / max(1, length(selected))
}




plot.glasso.hFDR <- function(hFDR, FDP, cv.glasso, num.selected, lambda, sign.lambda = -1,...){
  xlab <- expression(Log(lambda))
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  
  nlambda <- length(lambda)
  
  plot.range <- range(hFDR$mean - hFDR$std, hFDR$mean + hFDR$std, FDP)
  plot.range[1] <- max(plot.range[1], 0)
  plot.range[2] <- min(plot.range[2], 1)
  plot.args <- list(x = sign.lambda * log(lambda),
                    y = hFDR$mean,
                    xlim = range(sign.lambda*log(lambda)),
                    ylim = plot.range,
                    xlab = xlab, ylab = "estimated FDR and scaled CV Error", type = "n")
  new.args <- list(...)
  if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  cv.range <- abs(diff(range(c(cv.glasso$mean - cv.glasso$std, cv.glasso$mean + cv.glasso$std))))
  do.call("plot", plot.args)
  lines(x = sign.lambda * log(lambda),
        y = (cv.glasso$mean - min(cv.glasso$mean - cv.glasso$std)) / cv.range * abs(diff(plot.range)) + min(plot.range),
        col = "dodgerblue3")
  error.bars(sign.lambda * log(lambda),
             (cv.glasso$mean + cv.glasso$std - min(cv.glasso$mean - cv.glasso$std)) / cv.range * abs(diff(plot.range)) + min(plot.range),
             (cv.glasso$mean - cv.glasso$std - min(cv.glasso$mean - cv.glasso$std)) / cv.range * abs(diff(plot.range)) + min(plot.range),
             width = 0.01, col = alpha("dodgerblue3", 0.1))
  lines(x = sign.lambda * log(lambda),
        y = FDP,
        col = "black", lty = 2)
  error.bars(sign.lambda * log(lambda),
             hFDR$mean + hFDR$std, hFDR$mean - hFDR$std,
             width = 0.01, alpha("red", 0.3))
  points(sign.lambda*log(lambda), hFDR$mean,
         pch = 20, col = "red")
  axis(side = 3, at = sign.lambda*log(lambda),
       labels = paste(num.selected), tick = FALSE, line = 0)
  # abline(v = sign.lambda * log(lambda.min), lty = 3)
  # abline(v = sign.lambda * log(lambda.1se), lty = 3)
  legend("bottomright", inset = 0.05,
         legend = c("CV MSE", "FDP", "hFDR"),
         lty = c(1, 2, 1), lwd = c(1, 1, 0),
         pch = c(26, 26, 19), col = c("dodgerblue3", "black", "red"))
  invisible()
}

error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}


# compare.data <- data.frame(log_rho = log(tune_seq), value = FDR.est_seq, type = "hFDR")
# compare.data <- rbind(compare.data, data.frame(log_rho = log(tune_seq), value = FDP.true_seq, type = "'True' FDP"))
# compare.data <- rbind(compare.data, data.frame(log_rho = log(tune_seq), value = num.selected / 20, type = "# selected"))
# ggplot(compare.data, aes(x = log_rho, y = value, color = type)) +
#   geom_line() +
#   geom_point() +
#   scale_x_reverse() +
#   labs(x = "log(rho)", y = "Estiamted FDR / 'True' FDP")