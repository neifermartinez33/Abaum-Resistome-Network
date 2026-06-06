suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(graphlayouts)
  library(ggforce)
  library(scales)
  library(showtext)
})

font_add("TeX Gyre Pagella",
         regular    = "/home/miguel/.fonts/texgyrepagella-regular.otf",
         bold       = "/home/miguel/.fonts/texgyrepagella-bold.otf",
         italic     = "/home/miguel/.fonts/texgyrepagella-italic.otf",
         bolditalic = "/home/miguel/.fonts/texgyrepagella-bolditalic.otf")
showtext_auto()
showtext_opts(dpi = 600)
FONT <- "TeX Gyre Pagella"

base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

# WHO whitelist — ONLY truly isolated genes that appear as named grey nodes
# in the reference image. Connected genes all receive numeric labels (1-26).
# Critically: blaoxa_235like, blaoxa_241like, blaoxa_58like, blandm are
# numbered 23-26 in the reference → they are NOT in this list.
# blakpc is an isolated grey node → in this list (gets no number, grey label).
who_genes <- c(
  "blakpc",          # isolated, grey, left-middle
  "blaoxa_143like",  # isolated, grey, left-middle
  "blavim",          # isolated, grey, left-low
  "mcr-4.3",         # isolated, grey, bottom-left
  "mcr-4.7",         # isolated, grey, bottom-center
  "tet(x3)",         # isolated, grey, right-top
  "tet(x5)",         # isolated, grey, top-right
  "pmrb",            # isolated, grey, right-mid
  "rpob",            # isolated, grey, right-mid
  "merr",            # isolated, grey, bottom-right
  "mert",            # isolated, grey, right-low
  "arr",             # isolated, grey, top-center
  "arr-3"            # isolated, grey, top-left-center
)

# Number assignment: community order C1->C5, degree desc within community
nodes_conn <- nodes[degree > 0]
comm_order <- c("C1","C2","C3","C4","C5")
nodes_conn[, comm_ord := match(comm_lab, comm_order)]
setorder(nodes_conn, comm_ord, -degree)
nodes_conn[, node_num := NA_integer_]
counter <- 1L
for (i in seq_len(nrow(nodes_conn))) {
  if (!nodes_conn$gene[i] %in% who_genes) {
    nodes_conn$node_num[i] <- counter
    counter <- counter + 1L
  }
}

isolated_v  <- sort(nodes[degree == 0, gene])
connected_v <- nodes_conn$gene
edges_conn  <- edges_dt[gene_a %in% connected_v & gene_b %in% connected_v]

g <- graph_from_data_frame(
  d        = edges_conn[, .(from = gene_a, to = gene_b,
                             weight = jaccard, edge_type = source)],
  directed = FALSE,
  vertices = nodes_conn[, .(name = gene, degree, betweenness,
                             comm_lab, is_hub, node_num)]
)

set.seed(42)
w <- E(g)$weight
w[is.nan(w) | is.na(w)] <- median(w, na.rm = TRUE)
E(g)$weight <- w
pos <- layout_with_stress(g, weights = sqrt(E(g)$weight))
pos <- norm_coords(pos, xmin = -0.72, xmax = 0.72,
                   ymin = -0.72, ymax = 0.72)

blandm_idx <- which(V(g)$name == "blandm")
if (length(blandm_idx) > 0) {
  pos[blandm_idx, 1] <- 0.35
  pos[blandm_idx, 2] <- -0.80
}

lay <- create_layout(g, layout = "manual",
                     x = pos[, 1], y = pos[, 2])

# Isolated nodes: fixed positions matching reference image
# Reference layout (clockwise from top):
# arr (top-center), tet(x5) (top-right), tet(x3) (right-top),
# rpob (right-mid-top), pmrb (right-mid), mert (right-mid-low),
# merr (bottom-right), mcr-4.7 (bottom-center), mcr-4.3 (bottom-left),
# blavim (left-low), blaoxa_143like (left-mid), blakpc (left-mid-top),
# arr-3 (top-left-center)
iso_positions <- list(
  "arr"            = c( 0.10,  1.08),
  "arr-3"          = c(-0.38,  1.00),
  "tet(x5)"        = c( 0.72,  1.05),
  "tet(x3)"        = c( 1.10,  0.60),
  "rpob"           = c( 1.18,  0.18),
  "pmrb"           = c( 1.18, -0.20),
  "mert"           = c( 1.10, -0.52),
  "merr"           = c( 0.78, -0.88),
  "mcr-4.7"        = c( 0.18, -1.05),
  "mcr-4.3"        = c(-0.42, -0.95),
  "blavim"         = c(-0.90, -0.62),
  "blaoxa_143like" = c(-1.05, -0.10),
  "blakpc"         = c(-1.00,  0.28)
)

# Build iso data frame in the order of isolated_v
iso_df <- data.frame(
  gene = isolated_v,
  stringsAsFactors = FALSE
)
iso_df$x <- sapply(iso_df$gene, function(g) {
  if (g %in% names(iso_positions)) iso_positions[[g]][1] else NA_real_
})
iso_df$y <- sapply(iso_df$gene, function(g) {
  if (g %in% names(iso_positions)) iso_positions[[g]][2] else NA_real_
})
# Drop any with no position defined
iso_df <- iso_df[!is.na(iso_df$x), ]

# Label hjust: right-side nodes get hjust=0, left-side get hjust=1
iso_df$hjust <- ifelse(iso_df$x >= 0, 0, 1)
# Label offset direction
iso_df$lx <- iso_df$x + ifelse(iso_df$x >= 0,  0.10, -0.10)
iso_df$ly <- iso_df$y + 0.04

comm_pal <- c(
  "C1" = "#D55E00", "C2" = "#E69F00", "C3" = "#009E73",
  "C4" = "#0072B2", "C5" = "#CC79A7", "Isolated" = "#AAAAAA"
)
comm_full <- c(
  "C1" = "C1: Beta-lactam / Aminoglycoside",
  "C2" = "C2: Aminoglycoside / PBP3",
  "C3" = "C3: MDR core",
  "C4" = "C4: Quinolone / OXA",
  "C5" = "C5: Macrolide"
)

ell_data <- data.frame(x = lay$x, y = lay$y, comm_lab = lay$comm_lab)
ell_data <- ell_data[ell_data$comm_lab != "Isolated", ]
comm_counts <- table(ell_data$comm_lab)
ell_data <- ell_data[ell_data$comm_lab %in%
                       names(comm_counts[comm_counts >= 3]), ]

# Gene index — 5 columns
num_genes <- nodes_conn[!is.na(node_num)][order(node_num)]
n_total   <- nrow(num_genes)
n_cols    <- 5L
n_rows    <- ceiling(n_total / n_cols)

# Legend position
lx0  <- -1.78
ly0  <-  1.50
lsep <-  0.115

jacc_shift <- 10.1

# Gene index position — below network, with enough column spacing
ix_y0    <- -1.50
ix_dy    <- -0.095
# Fixed column x-positions (absolute, not relative) — wider spacing
ix_col_x <- c(-1.78, -0.96, -0.14,  0.68,  1.22)

p <- ggraph(lay) +

  geom_mark_ellipse(
    data = ell_data,
    aes(x = x, y = y, group = comm_lab, color = comm_lab),
    fill = NA, linewidth = 0.60, linetype = "dashed",
    expand = unit(7, "mm"),
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = comm_pal, guide = "none") +

  # bnlearn-only edges (light grey)
  geom_edge_link(
    aes(filter = (edge_type != "Jaccard+bnlearn")),
    edge_alpha = 0.35, color = "#C5C5C5",
    width = 0.16, show.legend = FALSE
  ) +
  # Consensus edges (blue gradient)
  geom_edge_link(
    aes(filter      = (edge_type == "Jaccard+bnlearn"),
        edge_colour = weight,
        edge_alpha  = weight),
    width = 0.85, show.legend = FALSE
  ) +
  scale_edge_colour_gradientn(
    colours = c("#C9E8F5","#6BAED6","#2171B5","#084594","#03224C"),
    limits  = c(0.3, 1.0), guide = "none"
  ) +
  scale_edge_alpha(range = c(0.18, 0.88), guide = "none") +

  # Non-hub nodes (circles)
  geom_node_point(
    aes(filter = (is_hub == 0 & degree > 0),
        size   = 2.5 + degree * 0.40,
        fill   = comm_lab),
    shape = 21, color = "white", stroke = 0.55, alpha = 0.95
  ) +
  # Hub nodes (diamonds)
  geom_node_point(
    aes(filter = (is_hub == 1),
        size   = 4.8 + degree * 0.55,
        fill   = comm_lab),
    shape = 23, color = "white", stroke = 0.55, alpha = 1.0
  ) +
  scale_fill_manual(values = comm_pal, guide = "none") +
  scale_size_identity() +

  # Numbers beside non-hub nodes — BLACK, upper-right nudge
  geom_node_text(
    aes(filter = (!is.na(node_num) & is_hub == 0),
        label  = node_num),
    size        = 5.0,
    fontface    = "bold",
    family      = FONT,
    color       = "#111111",
    nudge_x     = 0.028,
    nudge_y     = 0.028,
    show.legend = FALSE
  ) +
  # Numbers inside hub nodes — WHITE, centred
  geom_node_text(
    aes(filter = (!is.na(node_num) & is_hub == 1),
        label  = node_num),
    size        = 5.0,
    fontface    = "bold",
    family      = FONT,
    color       = "#111111",
    show.legend = FALSE
  ) +

  # Isolated nodes — grey circles at fixed positions
  annotate("point",
           x     = iso_df$x,
           y     = iso_df$y,
           shape = 21, size = 4.5,
           fill  = "#AAAAAA", color = "white", stroke = 0.45) +
  # Isolated node labels — grey text
  annotate("text",
           x        = iso_df$lx,
           y        = iso_df$ly,
           label    = iso_df$gene,
           size     = 5.0,
           fontface = "plain",
           family   = FONT,
           color    = "#555555",
           hjust    = iso_df$hjust) +

  # TOP-LEFT LEGEND
  annotate("text", x = lx0, y = ly0,
           label = "Node type", hjust = 0, fontface = "bold",
           size = 5.0, color = "#222222", family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - lsep,
           shape = 21, size = 4.8, fill = "#777777",
           color = "white", stroke = 0.6) +
  annotate("text", x = lx0 + 0.18, y = ly0 - lsep,
           label = "Gene", hjust = 0, size = 5.0,
           color = "#333333", family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 2*lsep,
           shape = 23, size = 5.5, fill = "#777777",
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 2*lsep,
           label = "Hub", hjust = 0, size = 5.0,
           color = "#333333", family = FONT) +

  annotate("text", x = lx0, y = ly0 - 3.2*lsep,
           label = "Community", hjust = 0, fontface = "bold",
           size = 5.0, color = "#222222", family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 4.2*lsep,
           shape = 21, size = 4.8, fill = comm_pal["C1"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 4.2*lsep,
           label = comm_full["C1"], hjust = 0, size = 5.0,
           color = comm_pal["C1"], family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 5.2*lsep,
           shape = 21, size = 4.8, fill = comm_pal["C2"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 5.2*lsep,
           label = comm_full["C2"], hjust = 0, size = 5.0,
           color = comm_pal["C2"], family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 6.2*lsep,
           shape = 21, size = 4.8, fill = comm_pal["C3"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 6.2*lsep,
           label = comm_full["C3"], hjust = 0, size = 5.0,
           color = comm_pal["C3"], family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 7.2*lsep,
           shape = 21, size = 4.8, fill = comm_pal["C4"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 7.2*lsep,
           label = comm_full["C4"], hjust = 0, size = 5.0,
           color = comm_pal["C4"], family = FONT) +
  annotate("point", x = lx0 + 0.06, y = ly0 - 8.2*lsep,
           shape = 21, size = 4.8, fill = comm_pal["C5"],
           color = "white", stroke = 0.5) +
  annotate("text", x = lx0 + 0.18, y = ly0 - 8.2*lsep,
           label = comm_full["C5"], hjust = 0, size = 5.0,
           color = comm_pal["C5"], family = FONT) +

  # Jaccard gradient bar
  annotate("text", x = lx0, y = ly0 - (9.8 + jacc_shift)*lsep,
           label = "Jaccard index", hjust = 0, fontface = "bold",
           size = 5.0, color = "#222222", family = FONT) +
  annotate("rect",
           xmin = lx0 + seq(0, 0.56, length.out = 20),
           xmax = lx0 + seq(0, 0.56, length.out = 20) + 0.031,
           ymin = ly0 - (11.2 + jacc_shift)*lsep,
           ymax = ly0 - (10.5 + jacc_shift)*lsep,
           fill = colorRampPalette(
             c("#C9E8F5","#6BAED6","#2171B5","#084594","#03224C"))(20),
           color = NA) +
  annotate("rect",
           xmin = lx0, xmax = lx0 + 0.59,
           ymin = ly0 - (11.2 + jacc_shift)*lsep,
           ymax = ly0 - (10.5 + jacc_shift)*lsep,
           fill = NA, color = "#AAAAAA", linewidth = 0.3) +
  annotate("text", x = lx0,        y = ly0 - (11.7 + jacc_shift)*lsep,
           label = "0.3",  hjust = 0.5, size = 5.0,
           color = "#555555", family = FONT) +
  annotate("text", x = lx0 + 0.30, y = ly0 - (11.7 + jacc_shift)*lsep,
           label = "0.65", hjust = 0.5, size = 5.0,
           color = "#555555", family = FONT) +
  annotate("text", x = lx0 + 0.59, y = ly0 - (11.7 + jacc_shift)*lsep,
           label = "1.0",  hjust = 0.5, size = 5.0,
           color = "#555555", family = FONT) +

  # GENE INDEX — 5 columns, fixed x positions, coloured by community
  annotate("text", x = -1.78, y = ix_y0 + 0.08,
           label = "Gene index", hjust = 0, fontface = "bold",
           size = 5.0, color = "#222222", family = FONT) +
  {
    col_idx <- ceiling(seq_len(n_total) / n_rows)   # 1..5
    row_idx <- (seq_len(n_total) - 1) %% n_rows     # 0-based row
    x_pos   <- ix_col_x[col_idx]
    y_pos   <- ix_y0 + ix_dy * row_idx
    annotate("text",
             x        = x_pos,
             y        = y_pos,
             label    = paste0(num_genes$node_num, ". ", num_genes$gene),
             hjust    = 0,
             size     = 5.0,
             family   = FONT,
             color    = comm_pal[num_genes$comm_lab])
  } +

  coord_fixed(expand = FALSE, xlim = c(-1.82, 1.50),
              ylim = c(ix_y0 + ix_dy * n_rows - 0.10, 1.60)) +
  labs(title = NULL, subtitle = NULL) +
  theme_graph(base_family = FONT) +
  theme(
    text             = element_text(family = FONT),
    plot.background  = element_rect(fill = "white", color = NA, linewidth = 0),
    panel.background = element_rect(fill = "white", color = NA, linewidth = 0),
    panel.border     = element_blank(),
    plot.margin      = margin(0, 0, 0, 0, unit = "pt")
  )

png_path  <- file.path(out_dir, "Fig1_consensus_network_communities.png")
tiff_path <- file.path(out_dir, "Fig1_consensus_network_communities.tiff")

# Save with zero padding: use ragg for clean transparent-border output
ggsave(png_path,  p, width = 12, height = 13, dpi = 600, bg = "white")
ggsave(tiff_path, p, width = 12, height = 13, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_path,  "\n")
cat("\u2713 TIFF:", tiff_path, "\n")
