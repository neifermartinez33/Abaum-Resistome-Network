library(data.table)
library(igraph)

cat("=== Phase 9: Topological Analysis ===\n")
cat(format(Sys.time()), "\n\n")

g <- readRDS("~/matrix_output_v2/network_consensus.rds")
cat("Network loaded:", vcount(g), "nodes,", ecount(g), "edges\n\n")

# ── 1. CENTRALITY MEASURES ───────────────────────────────────
cat("=== 1. Centrality Measures ===\n")
deg       <- degree(g)
btw       <- betweenness(g, normalized = TRUE)
clo       <- closeness(g, normalized = TRUE)
eigen_c   <- eigen_centrality(g)$vector
trans_loc <- transitivity(g, type = "local")
trans_loc[is.nan(trans_loc)] <- 0

centrality <- data.table(
  gene        = V(g)$name,
  degree      = deg,
  betweenness = round(btw, 6),
  closeness   = round(clo, 6),
  eigenvector = round(eigen_c, 6),
  clustering  = round(trans_loc, 6)
)
centrality <- centrality[order(-degree)]
cat("\nTop 15 nodes by degree:\n")
print(centrality[1:15])

# ── 2. HUB IDENTIFICATION ────────────────────────────────────
cat("\n=== 2. Hub Identification (top 20% degree OR top 20% betweenness) ===\n")
deg_thresh <- quantile(deg[deg > 0], 0.80)
btw_thresh <- quantile(btw[btw > 0], 0.80)
cat("Degree threshold (80th percentile):", deg_thresh, "\n")
cat("Betweenness threshold (80th percentile):", round(btw_thresh, 6), "\n")
hubs <- centrality[(degree >= deg_thresh) | (betweenness >= btw_thresh)]
hubs <- hubs[order(-degree)]
cat("Hubs identified:", nrow(hubs), "\n")
print(hubs)

# ── 3. GLOBAL NETWORK STATISTICS ────────────────────────────
cat("\n=== 3. Global Network Statistics ===\n")
g_conn <- induced_subgraph(g, V(g)[deg > 0])
cat("Global clustering (transitivity):", round(transitivity(g, type="global"), 4), "\n")
cat("Mean local clustering:", round(mean(trans_loc[deg > 0]), 4), "\n")
cat("Diameter (connected subgraph):", diameter(g_conn), "\n")
cat("Mean path length:", round(mean_distance(g_conn), 4), "\n")
cat("Degree assortativity:", round(assortativity_degree(g, directed=FALSE), 4), "\n")

# ── 4. DEGREE DISTRIBUTION ──────────────────────────────────
cat("\n=== 4. Degree Distribution / Power Law Fit ===\n")
deg_nonzero <- deg[deg > 0]
fit <- fit_power_law(deg_nonzero)
cat("Alpha exponent:", round(fit$alpha, 4), "\n")
cat("KS statistic:", round(fit$KS.stat, 4), "\n")
ks_p <- fit$KS.p
if (is.null(ks_p) || is.na(ks_p)) {
  cat("KS p-value: Not calculable (small sample n=", length(deg_nonzero), ")\n")
} else {
  cat("KS p-value:", round(ks_p, 4), "\n")
}
cat("xmin:", fit$xmin, "\n")
cat("Interpretation: alpha =", round(fit$alpha, 2),
    ">> 3.0 indicates distribution is NOT classical scale-free\n")

# ── 5. COMMUNITY DETECTION ──────────────────────────────────
cat("\n=== 5. Community Detection ===\n")

set.seed(42)
louvain  <- cluster_louvain(g, weights = E(g)$weight)
walktrap <- cluster_walktrap(g, weights = E(g)$weight)
infomap  <- cluster_infomap(g, e.weights = E(g)$weight)

cat("Louvain  - Modularity:", round(modularity(louvain), 4),
    "| Communities:", length(louvain), "\n")
cat("Walktrap - Modularity:", round(modularity(walktrap), 4),
    "| Communities:", length(walktrap), "\n")
cat("Infomap  - Modularity:", round(modularity(infomap), 4),
    "| Communities:", length(infomap), "\n")

nmi_lw <- compare(louvain, walktrap, method = "nmi")
nmi_li <- compare(louvain, infomap,  method = "nmi")
cat("NMI Louvain-Walktrap:", round(nmi_lw, 4), "\n")
cat("NMI Louvain-Infomap: ", round(nmi_li, 4), "\n")

cat("\nLouvain communities with >= 3 nodes:\n")
louvain_sizes <- table(membership(louvain))
big <- as.integer(names(louvain_sizes[louvain_sizes >= 3]))
for (cid in big) {
  members <- V(g)$name[membership(louvain) == cid]
  cat(sprintf("  Community %d (%d nodes): %s\n",
              cid, length(members), paste(members, collapse=", ")))
}

# ── 6. SAVE OUTPUTS ─────────────────────────────────────────
V(g)$degree             <- deg
V(g)$betweenness        <- btw
V(g)$closeness          <- clo
V(g)$eigenvector        <- eigen_c
V(g)$clustering         <- trans_loc
V(g)$community_louvain  <- membership(louvain)
V(g)$community_walktrap <- membership(walktrap)
V(g)$is_hub             <- as.integer(V(g)$name %in% hubs$gene)

saveRDS(g, "~/matrix_output_v2/network_topology.rds")
fwrite(centrality, "~/matrix_output_v2/topology_centrality.csv")

community_table <- data.table(
  gene                = V(g)$name,
  community_louvain   = membership(louvain),
  community_walktrap  = membership(walktrap),
  community_infomap   = membership(infomap),
  is_hub              = V(g)$is_hub
)
fwrite(community_table, "~/matrix_output_v2/topology_communities.csv")

cat("\nFiles saved:\n")
cat("  network_topology.rds\n")
cat("  topology_centrality.csv\n")
cat("  topology_communities.csv\n")
cat("\nDone:", format(Sys.time()), "\n")
