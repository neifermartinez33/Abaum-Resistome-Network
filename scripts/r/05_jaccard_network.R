library(data.table)
library(igraph)

cat("=== Phase 8: Jaccard Filtered ===\n")
cat(format(Sys.time()), "\n\n")

# Load matrix
mat <- fread("~/matrix_output_v2/matrix1_gene_families.csv")
genomes <- mat[[1]]
mat[, 1 := NULL]
M <- as.matrix(mat)
rownames(M) <- genomes

cat("Dimensions:", nrow(M), "x", ncol(M), "\n")
cat("Families:", colnames(M), "\n\n")

# Filter 1: prevalence >= 5% (already met by design, but we verify)
prev <- colMeans(M)
cat("Prevalences min/max:", round(min(prev),3), "/", round(max(prev),3), "\n")

# Calculate pairwise Jaccard
n_genes <- ncol(M)
genes <- colnames(M)
results <- list()

for (i in 1:(n_genes-1)) {
  for (j in (i+1):n_genes) {
    a <- M[, i]
    b <- M[, j]
    inter <- sum(a == 1 & b == 1)
    union <- sum(a == 1 | b == 1)
    jacc <- if (union == 0) 0 else inter / union
    cooc <- inter / nrow(M)  # co-occurrence
    results[[length(results)+1]] <- data.frame(
      gene_a = genes[i], gene_b = genes[j],
      jaccard = round(jacc, 4),
      cooccurrence = round(cooc, 4),
      n_shared = inter
    )
  }
}

df <- rbindlist(results)
cat("\nTotal pairs calculated:", nrow(df), "\n")

# Filters: co-occurrence >= 10% AND Jaccard > 0.3
df_filtered <- df[cooccurrence >= 0.10 & jaccard > 0.30]
cat("Pairs after filter (cooc>=10%, Jaccard>0.3):", nrow(df_filtered), "\n")

# Build igraph network
g <- graph_from_data_frame(df_filtered[, .(gene_a, gene_b, jaccard)],
                            directed = FALSE,
                            vertices = genes)
E(g)$weight <- df_filtered$jaccard

cat("\n=== Jaccard Network Topology ===\n")
cat("Nodes:", vcount(g), "\n")
cat("Edges:", ecount(g), "\n")
cat("Density:", round(edge_density(g), 4), "\n")
cat("Components:", components(g)$no, "\n")
cat("Isolated nodes:", sum(degree(g) == 0), "\n")

# Save outputs
fwrite(df, "~/matrix_output_v2/jaccard_all_pairs.csv")
fwrite(df_filtered, "~/matrix_output_v2/jaccard_filtered_edges.csv")
saveRDS(g, "~/matrix_output_v2/jaccard_network.rds")

cat("\nFiles saved:\n")
cat("  jaccard_all_pairs.csv\n")
cat("  jaccard_filtered_edges.csv\n")
cat("  jaccard_network.rds\n")
cat("\nDone:", format(Sys.time()), "\n")
