library(igraph)
library(data.table)

base <- "/home/miguel/Abaum_Resistome_Network"

# Load consensus network
net <- fread(file.path(base, "results/supplementary/tables_submission/TableS5_network_consensus_edges.csv"))

# Build igraph object with all 41 genes
all_genes <- unique(c(net$Gene_A, net$Gene_B))

# Add isolated nodes from topology
cent <- fread(file.path(base, "results/tables/topology_centrality.csv"))
all_genes <- unique(c(all_genes, cent$gene))

# Create edges (only non-NA Jaccard = Jaccard+bnlearn edges)
edges_connected <- net[Validation_source == "Jaccard+bnlearn" | 
                       Validation_source == "bnlearn_only"]

g <- graph_from_data_frame(
  edges_connected[, .(Gene_A, Gene_B)],
  directed = FALSE,
  vertices = all_genes
)

cat("=== NETWORK STATS ===\n")
cat("Nodes:", vcount(g), "\n")
cat("Edges:", ecount(g), "\n")
cat("Components:", components(g)$no, "\n")

set.seed(42)

# LOUVAIN
comm_louvain <- cluster_louvain(g)
q_louvain    <- modularity(comm_louvain)
cat("\n=== LOUVAIN ===\n")
cat("Q:", round(q_louvain, 4), "\n")
cat("Communities:", length(comm_louvain), "\n")
cat("Sizes:", table(sizes(comm_louvain)), "\n")

# WALKTRAP
comm_walktrap <- cluster_walktrap(g)
q_walktrap    <- modularity(comm_walktrap)
nmi_lw        <- compare(comm_louvain, comm_walktrap, method="nmi")
cat("\n=== WALKTRAP ===\n")
cat("Q:", round(q_walktrap, 4), "\n")
cat("Communities:", length(comm_walktrap), "\n")
cat("NMI vs Louvain:", round(nmi_lw, 4), "\n")

# INFOMAP
comm_infomap <- cluster_infomap(g)
q_infomap    <- modularity(comm_infomap)
nmi_li       <- compare(comm_louvain, comm_infomap, method="nmi")
cat("\n=== INFOMAP ===\n")
cat("Q:", round(q_infomap, 4), "\n")
cat("Communities:", length(comm_infomap), "\n")
cat("NMI vs Louvain:", round(nmi_li, 4), "\n")

# SUMMARY TABLE
summary_df <- data.frame(
  Algorithm     = c("Louvain", "Walktrap", "Infomap"),
  Modularity_Q  = round(c(q_louvain, q_walktrap, q_infomap), 4),
  N_communities = c(length(comm_louvain),
                    length(comm_walktrap),
                    length(comm_infomap)),
  NMI_vs_Louvain = round(c(1.000, nmi_lw, nmi_li), 4),
  Seed = 42
)

cat("\n=== SUMMARY ===\n")
print(summary_df)

# Save
out <- file.path(base,
  "results/supplementary/tables_submission/TableS12_community_validation.csv")
write.csv(summary_df, out, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out))
