library(tidyverse)

load_gene_info = function() read_csv('data/gene-info.csv')

load_correlations = function () read_csv('data/corrs.csv')
load_corr_matrix = function () load_correlations() |> as.matrix()

# --- Network construction ---
connectivity = function(m) rowSums(m, na.rm=TRUE)

adjacency = function(corr_matrix, beta) {
  a = abs(corr_matrix)
  0 -> a[is.na(a)]
  0 -> diag(a)

  a^beta
}

tom = function(corr_matrix, beta) {
  a = adjacency(corr_matrix, beta)

  k = connectivity(a)
  k_min = outer(k, k, pmin)

  (a %*% a + a) / (k_min + 1 - a)
}

# --- Scale-free topology ---
scale_free_r2 = function(m) {
  k = connectivity(m)
  binned_k = hist(k, plot=FALSE)
  p_k = binned_k$counts / sum(binned_k$counts)

  nonzero = p_k > 0
  p_k[nonzero] -> p_k
  use_k = binned_k$mids[nonzero]

  cor(log10(p_k), log10(use_k))^2
}

tabulate_betas = function(corr_matrix, betas = 1:20) {
  map_dfr(betas, function(beta) {
    a = adjacency(corr_matrix, beta)
    tibble(beta = beta, r2 = scale_free_r2(a), mean_k = mean(connectivity(a)))
  })
}

# --- Clustering ---
cluster_genes = function(tom) {
  distances = as.dist(1 - tom)
  hclust(distances, method="average")
}