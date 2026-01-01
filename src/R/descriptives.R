library(ggplot2)

source('src/R/create_network.R')

create_dataset = function(gene_info, correlations) {
  select(gene_info, transcript, is_canonical) |>
    mutate(is_canonical = if_else(is_canonical == "canonical", 1, 0)) -> 
    canonical_lookup

  mutate(correlations, gene1 = colnames(correlations)) |>
    pivot_longer(-gene1, names_to = 'gene2', values_to="correlation") |>
    filter(gene1 < gene2) |>
    filter(!is.na(correlation)) |>
    left_join(canonical_lookup, by = c('gene1' = 'transcript')) |>
    rename(is_canonical1 = is_canonical) |>
    left_join(canonical_lookup, by = c('gene2' = 'transcript')) |>
    rename(is_canonical2 = is_canonical) |>
    mutate(corr_type = case_when(
      is_canonical1 == 0 & is_canonical2 == 0 ~ "both_no",
      is_canonical1 == 1 & is_canonical2 == 1 ~ "both_yes",
      is_canonical1 != is_canonical2 ~ "unsame"))
}

plot_histograms <- function(x) {
  ggplot(x, aes(x = correlation, fill=factor(corr_type))) +
    geom_histogram(show.legend = FALSE) +
    geom_vline(xintercept = 0.5) +
    facet_wrap(~corr_type, scales='free_y')
}
