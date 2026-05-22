library(data.table)
library(bnlearn)
library(igraph)
library(parallel)

cat("=== Phase 8: Bayesian Network (bnlearn) - Strict Threshold ===\n")
cat(format(Sys.time()), "\n\n")

# Load matrix
mat <- fread("~/matrix_output_v2/matrix1_gene_families.csv")
genomes <- mat[[1]]
mat[, 1 := NULL]
M <- as.matrix(mat)
rownames(M) <- genomes

cat("Dimensions:", nrow(M), "x", ncol(M), "\n")

# Convert to factor data.frame
df <- as.data.frame(M)
for (col in colnames(df)) {
  df[[col]] <- factor(df[[col]], levels = c(0, 1))
}
cat("Data converted to factors\n")

# Filter variables with prevalence < 1% or > 99%
prev <- colMeans(M)
keep <- prev > 0.01 & prev < 0.99
cat("Variables retained for bnlearn:", sum(keep), "/", ncol(M), "\n")
df_bn <- df[, keep]

# Create cluster
cl <- makeCluster(32)
cat("\nCluster created with 32 nodes\n")

# Bootstrap Hill Climbing, 500 reps, BDe score
cat("Starting bootstrap Hill Climbing (500 reps, 32 cores)...\n")
cat(format(Sys.time()), "\n\n")

set.seed(42)
boot_strength <- tryCatch({
  boot.strength(
    df_bn,
    R              = 500,
    algorithm      = "hc",
    algorithm.args = list(score = "bde", iss = 10),
    cpdag          = TRUE,
    cluster        = cl
  )
}, error = function(e) {
  stopCluster(cl)
  stop(e)
})

stopCluster(cl)
cat("Bootstrap completed:", format(Sys.time()), "\n\n")

# Save complete boot_strength for reference
fwrite(as.data.table(boot_strength), "~/matrix_output_v2/bnlearn_boot_strength.csv")
saveRDS(boot_strength, "~/matrix_output_v2/bnlearn_boot_strength.rds")

# Score distribution
cat("=== Strength Score Distribution ===\n")
breaks <- c(0, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0)
for (i in 1:(length(breaks)-1)) {
  n <- sum(boot_strength$strength >= breaks[i] & boot_strength$strength < breaks[i+1])
  cat(sprintf("  [%.2f - %.2f): %d edges\n", breaks[i], breaks[i+1], n))
}
cat("  [1.0]:", sum(boot_strength$strength == 1.0), "edges\n\n")

# Strict threshold: strength >= 0.85 AND direction >= 0.5
THRESHOLD <- 0.85
cat("Applying threshold: strength >=", THRESHOLD, "and direction >= 0.5\n")

avg_net <- averaged.network(boot_strength, threshold = THRESHOLD)
edges_bn <- as.data.frame(arcs(avg_net))
colnames(edges_bn) <- c("from", "to")

# Add strength and direction
edges_bn$strength <- boot_strength$strength[
  match(paste(edges_bn$from, edges_bn$to),
        paste(boot_strength$from, boot_strength$to))
]
edges_bn$direction <- boot_strength$direction[
  match(paste(edges_bn$from, edges_bn$to),
        paste(boot_strength$from, boot_strength$to))
]

# Filter direction >= 0.5
edges_bn <- edges_bn[edges_bn$direction >= 0.5, ]
cat("Final edges (strength>=0.85, direction>=0.5):", nrow(edges_bn), "\n")

# Build undirected igraph network
all_nodes <- colnames(df_bn)
g_bn <- graph_from_data_frame(edges_bn, directed = FALSE, vertices = all_nodes)

cat("\n=== bnlearn Topology (threshold=0.85) ===\n")
cat("Nodes:", vcount(g_bn), "\n")
cat("Edges:", ecount(g_bn), "\n")
cat("Density:", round(edge_density(g_bn), 4), "\n")
cat("Components:", components(g_bn)$no, "\n")
cat("Isolated nodes:", sum(degree(g_bn) == 0), "\n")

# Overlap with Jaccard
jacc <- fread("~/matrix_output_v2/jaccard_filtered_edges.csv")
jacc_pairs <- paste(pmin(jacc$gene_a, jacc$gene_b),
                    pmax(jacc$gene_a, jacc$gene_b))
bn_pairs <- paste(pmin(edges_bn$from, edges_bn$to),
                  pmax(edges_bn$from, edges_bn$to))
overlap <- sum(bn_pairs %in% jacc_pairs)
cat("\nOverlap Jaccard-bnlearn:", overlap, "/", nrow(edges_bn), "edges\n")

cat("\n=== Top 15 edges by strength ===\n")
top15 <- edges_bn[order(-edges_bn$strength), ][1:min(15, nrow(edges_bn)), ]
print(top15)

# Save
saveRDS(avg_net, "~/matrix_output_v2/bnlearn_avg_network.rds")
saveRDS(g_bn,    "~/matrix_output_v2/bnlearn_network.rds")
fwrite(as.data.table(edges_bn), "~/matrix_output_v2/bnlearn_edges.csv")

cat("\nFiles saved:\n")
cat("  bnlearn_boot_strength.rds/.csv\n")
cat("  bnlearn_avg_network.rds\n")
cat("  bnlearn_network.rds\n")
cat("  bnlearn_edges.csv\n")
cat("\nDone:", format(Sys.time()), "\n")
