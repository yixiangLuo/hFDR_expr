library(here)
library(glasso)
library(dagitty)

source(here("R", "methods.R"))


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


