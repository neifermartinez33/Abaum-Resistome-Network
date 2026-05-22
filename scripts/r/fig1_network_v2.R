suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(graphlayouts)
  library(ggforce)
  library(scales)
})

# ── PATHS ────────────────────────────────────────────────────
base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── DATA ────────────────────────────────────────────────────
edges_dt <- fread(file.path(tables, "network_consensus.csv"))
central  <- fread(file.path(tables, "topology_centrality.csv"))
comm_dt  <- fread(file.path(tables, "topology_communities.csv"))

nodes <- merge(central, comm_dt[, .(gene, community_louvain, is_hub)],
               by = "gene", all = TRUE)
nodes[is.na(degree), degree := 0]
nodes[, comm_lab := fcase(
  community_louvain == 1, "C1",
  community_louvain == 2, "C2",
  community_louvain == 3, "C3",
  community_louvain == 4, "C4",
  community_louvain == 5, "C5",
  default = "Isolated"
)]

# ── SEPARATE CONNECTED / ISOLATED ────────────────────────────
connected_v <- nodes[degree > 0, gene]
isolated_v  <- sort(nodes[degree == 0, gene])   # sorted alphabetically
nodes_conn  <- nodes[gene %in% connected_v]
edges_conn  <- edges_dt[gene_a %in% connected_v & gene_b %in% connected_v]

# ── GRAPH (connected nodes only) ─────────────────────────────
g <- graph_from_data_frame(
  d        = edges_conn[, .(from = gene_a, to = gene_b,
                             weight = jaccard, edge_type = source)],
  directed = FALSE,
  vertices = nodes_conn[, .(name = gene, degree, betweenness,
                             comm_lab, is_hub)]
)

# ── CONNECTED COMPONENT LAYOUT ───────────────────────────────
set.seed(42)
w <- E(g)$weight
w[is.nan(w) | is.na(w)] <- median(w, na.rm = TRUE)
E(g)$weight <- w
pos <- layout_with_stress(g, weights = sqrt(E(g)$weight))
# Normalize to [-0.72, 0.72] — leaves outer ring available
pos <- norm_coords(pos, xmin = -0.72, xmax = 0.72, ymin = -0.72, ymax = 0.72)

# Manually reposition blandm so its 3 edges are visible
# Its partners (arr-2, aph(3'), blatem) are in the lower core
# Move it down-left, away from the dense cluster
blandm_idx <- which(V(g)$name == "blandm")
if (length(blandm_idx) > 0) {
  pos[blandm_idx, 1] <- 0.35   # x: further right
  pos[blandm_idx, 2] <- -0.80  # y: below the core
}

lay <- create_layout(g, layout = "manual",
                     x = pos[, 1], y = pos[, 2])

# ── OUTER RING POSITIONS (isolated nodes) ─────────────────────
# Radius 1.18: sufficiently distant from the component [-0.72, 0.72]
# Nodes are distributed at equidistant angles
n_iso  <- length(isolated_v)
angles <- seq(pi/2, pi/2 + 2*pi, length.out = n_iso + 1)[-(n_iso + 1)]
r_iso  <- 1.05
iso_x  <- cos(angles) * r_iso
iso_y  <- sin(angles) * r_iso

# ── PALETTE ───────────────────────────────────────────────────
comm_pal <- c(
  "C1" = "#D55E00", "C2" = "#E69F00", "C3" = "#009E73",
  "C4" = "#0072B2", "C5" = "#CC79A7", "Isolated" = "#AAAAAA"
)
comm_full <- c(
  "C1" = "C1: Beta-lactam / Aminoglycoside",
  "C2" = "C2: Aminoglycoside / PBP3",
  "C3" = "C3: MDR Core",
  "C4" = "C4: Quinolone / OXA",
  "C5" = "C5: Macrolide"
)

# ── ELLIPSE DATA ──────────────────────────────────────────────
ell_data <- data.frame(x = lay$x, y = lay$y, comm_lab = lay$comm_lab)
ell_data <- ell_data[ell_data$comm_lab != "Isolated", ]
comm_counts <- table(ell_data$comm_lab)
ell_data    <- ell_data[ell_data$comm_lab %in%
                          names(comm_counts[comm_counts >= 3]), ]

# ── LEGEND PARAMETERS ─────────────────────────────────────────
# Shifted further left: lx0 = -1.55 (canvas expanded to -1.65)
lx0  <- -1.58
ly0  <- 1.42
lsep <- 0.100   # more compact

# ── FIGURE ───────────────────────────────────────────────────
p <- ggraph(lay) +

  # Dotted border ellipses
  geom_mark_ellipse(
    data = ell_data,
    aes(x = x, y = y, group = comm_lab, color = comm_lab),
    fill = NA, linewidth = 0.60, linetype = "dashed",
    expand = unit(6, "mm"),
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = comm_pal, guide = "none") +

  # bnlearn-only edges — fixed alpha: NaN in jaccard made them invisible
  geom_edge_link(
    aes(filter = (edge_type != "Jaccard+bnlearn")),
    edge_alpha  = 0.35,
    color       = "#C5C5C5",
    width       = 0.18,
    show.legend = FALSE
  ) +

  # Consensus edges
  geom_edge_link(
    aes(filter = (edge_type == "Jaccard+bnlearn"),
        edge_colour = weight, edge_alpha = weight),
    width = 0.85, show.legend = FALSE
  ) +
  scale_edge_colour_gradientn(
    colours = c("#C9E8F5","#6BAED6","#2171B5","#084594","#03224C"),
    limits = c(0.3, 1.0), guide = "none"
  ) +
  scale_edge_alpha(range = c(0.18, 0.88), guide = "none") +

  # Connected non-hub nodes
  geom_node_point(
    aes(filter = (is_hub == 0 & degree > 0),
        size   = 2.0 + degree * 0.35, fill = comm_lab),
    shape = 21, color = "white", stroke = 0.55, alpha = 0.95,
    show.legend = FALSE
  ) +

  # Hubs
  geom_node_point(
    aes(filter = (is_hub == 1),
        size   = 4.0 + degree * 0.50, fill = comm_lab),
    shape = 23, color = "white", stroke = 0.5, alpha = 1.0,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = comm_pal, guide = "none") +
  scale_size_identity() +

  # Non-hub gene labels
  geom_node_text(
    aes(filter = (degree > 0 & is_hub == 0), label = name, color = comm_lab),
    size = 2.2, repel = TRUE, bg.color = "white", bg.r = 0.10,
    max.overlaps = 15, seed = 42, fontface = "plain", show.legend = FALSE,
    point.padding = unit(0.18, "lines"), box.padding = unit(0.08, "lines"),
    force = 0.25, force_pull = 1.5
  ) +

  # Hub labels
  geom_node_text(
    aes(filter = (is_hub == 1), label = name, color = comm_lab),
    size = 2.8, fontface = "bold", repel = TRUE,
    bg.color = "white", bg.r = 0.12,
    max.overlaps = 15, seed = 42, show.legend = FALSE,
    point.padding = unit(0.32, "lines"), box.padding = unit(0.10, "lines"),
    force = 0.35, force_pull = 1.5
  ) +

  scale_color_manual(values = comm_pal, guide = "none") +

  # ── ISOLATED NODES: outer ring ──────────────────────────────
  annotate("point",
           x = iso_x, y = iso_y,
           shape = 21, size = 2.2,
           fill = "#AAAAAA", color = "white", stroke = 0.45) +
  # Isolated node labels: oriented radially outward
  annotate("text",
           x        = iso_x * 1.10,
           y        = iso_y * 1.10,
           label    = isolated_v,
           size     = 2.5,
           fontface = "bold",
           color    = "#666666",
           hjust    = ifelse(iso_x >= 0, 0, 1)) +

  # ── LEGEND — upper left corner ───────────────────────────────
  # Node type
  annotate("text", x = lx0, y = ly0,
           label = "Node type", hjust = 0, fontface = "bold",
           size = 3.5, color = "#222222") +
  annotate("point", x = lx0 + 0.06, y = ly0 - lsep,
           shape = 21, size = 4.5, fill = "#777777",
           color = "white", stroke = 0.6) +
  annotate("text", x = lx0 + 0.16, y = ly0 - lsep,
           label = "Gene", hjust = 0, size = 3.1, color = "#333333") +
  annotate("point", x = lx0 + 0.06, y = ly0 - 2*lsep,
           shape = 23, size = 5.2, fill = "#777777",
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 2*lsep,
           label = "Hub", hjust = 0, size = 3.1, color = "#333333") +

  # Community
  annotate("text", x = lx0, y = ly0 - 3.1*lsep,
           label = "Community", hjust = 0, fontface = "bold",
           size = 3.5, color = "#222222") +
  annotate("point", x = lx0 + 0.06, y = ly0 - 4.1*lsep,
           shape = 21, size = 4.5, fill = comm_pal["C1"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 4.1*lsep,
           label = comm_full["C1"], hjust = 0, size = 2.9,
           color = comm_pal["C1"]) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 5.1*lsep,
           shape = 21, size = 4.5, fill = comm_pal["C2"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 5.1*lsep,
           label = comm_full["C2"], hjust = 0, size = 2.9,
           color = comm_pal["C2"]) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 6.1*lsep,
           shape = 21, size = 4.5, fill = comm_pal["C3"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 6.1*lsep,
           label = comm_full["C3"], hjust = 0, size = 2.9,
           color = comm_pal["C3"]) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 7.1*lsep,
           shape = 21, size = 4.5, fill = comm_pal["C4"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 7.1*lsep,
           label = comm_full["C4"], hjust = 0, size = 2.9,
           color = comm_pal["C4"]) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 8.1*lsep,
           shape = 21, size = 4.5, fill = comm_pal["C5"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.16, y = ly0 - 8.1*lsep,
           label = comm_full["C5"], hjust = 0, size = 2.9,
           color = comm_pal["C5"]) +

  # ── JACCARD — lower left, smaller ────────────────────────────
  annotate("text", x = lx0, y = -1.22,
           label = "Jaccard index", hjust = 0, fontface = "bold",
           size = 2.9, color = "#222222") +
  annotate("rect",
           xmin = lx0 + seq(0, 0.56, length.out = 20),
           xmax = lx0 + seq(0, 0.56, length.out = 20) + 0.031,
           ymin = -1.33, ymax = -1.27,
           fill = colorRampPalette(
             c("#C9E8F5","#6BAED6","#2171B5","#084594","#03224C"))(20),
           color = NA) +
  annotate("rect",
           xmin = lx0, xmax = lx0 + 0.59,
           ymin = -1.33, ymax = -1.27,
           fill = NA, color = "#AAAAAA", linewidth = 0.3) +
  annotate("text", x = lx0,        y = -1.38, label = "0.3",
           hjust = 0.5, size = 2.2, color = "#555555") +
  annotate("text", x = lx0 + 0.30, y = -1.38, label = "0.65",
           hjust = 0.5, size = 2.2, color = "#555555") +
  annotate("text", x = lx0 + 0.59, y = -1.38, label = "1.0",
           hjust = 0.5, size = 2.2, color = "#555555") +

  coord_fixed(xlim = c(-1.65, 1.45), ylim = c(-1.50, 1.50)) +
  labs(title = NULL, subtitle = NULL) +
  theme_graph(base_family = "sans") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "none",
    plot.margin      = margin(8, 8, 8, 8)
  )

# ── EXPORT ───────────────────────────────────────────────────
out_path <- file.path(out_dir, "Fig1_consensus_network_communities.png")
ggsave(out_path, p, width = 14, height = 11, dpi = 300, bg = "white")
cat("OK: Fig1 saved at:", out_path, "\n")
