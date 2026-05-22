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
out_dir <- file.path(base, "results/figures/supplementary")
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

# ── RESISTANCE CLASS ──────────────────────────────────────────
res_class <- c(
  "aac(3)"        = "Aminoglycoside",
  "aac(6')"       = "Aminoglycoside",
  "aada"          = "Aminoglycoside",
  "ant(2'')"      = "Aminoglycoside",
  "aph(3'')"      = "Aminoglycoside",
  "aph(3')"       = "Aminoglycoside",
  "aph(6)"        = "Aminoglycoside",
  "arma"          = "Aminoglycoside",
  "abaf"          = "Fosfomycin",
  "adec"          = "Efflux pump",
  "blaimp"        = "Carbapenem",
  "blakpc"        = "Carbapenem",
  "blandm"        = "Carbapenem",
  "blavim"        = "Carbapenem",
  "blaoxa_23like" = "Carbapenem",
  "blaoxa_24like" = "Carbapenem",
  "blaoxa_58like" = "Carbapenem",
  "blaoxa_143like"= "Carbapenem",
  "blaoxa_235like"= "Carbapenem",
  "blatem"        = "Beta-lactam",
  "catb"          = "Phenicol",
  "dfra"          = "Trimethoprim",
  "ftsi"          = "Beta-lactam",
  "gyra"          = "Quinolone",
  "parc"          = "Quinolone",
  "mcr-4.3"       = "Colistin",
  "mcr-4.7"       = "Colistin",
  "merr"          = "Metal",
  "mert"          = "Metal",
  "mph(e)"        = "Macrolide",
  "msr(e)"        = "Macrolide",
  "pmrb"          = "Colistin",
  "qacedelta"     = "Disinfectant",
  "rpob"          = "Rifampicin",
  "sul"           = "Sulfonamide",
  "tet(b)"        = "Tetracycline",
  "tet(x3)"       = "Tetracycline",
  "tet(x5)"       = "Tetracycline",
  "arr"           = "Rifampicin",
  "arr-2"         = "Rifampicin",
  "arr-3"         = "Rifampicin"
)
nodes[, res_class := res_class[gene]]
nodes[is.na(res_class), res_class := "Other"]

# ── CLASS PALETTE — Paul Tol Bright ───────────────────────────
class_pal <- c(
  "Aminoglycoside" = "#E69F00",
  "Beta-lactam"    = "#56B4E9",
  "Carbapenem"     = "#7B2FBE",
  "Efflux pump"    = "#009E73",
  "Quinolone"      = "#CC6600",
  "Macrolide"      = "#CC79A7",
  "Tetracycline"   = "#8B6914",
  "Colistin"       = "#882255",
  "Fosfomycin"     = "#44AA99",
  "Trimethoprim"   = "#4477AA",
  "Sulfonamide"    = "#117733",
  "Disinfectant"   = "#BBAA00",
  "Rifampicin"     = "#332288",
  "Phenicol"       = "#AA3377",
  "Metal"          = "#888888",
  "Other"          = "#BBBBBB"
)

# ── COMMUNITY PALETTE (ellipses only) ─────────────────────────
comm_pal <- c(
  "C1" = "#D55E00", "C2" = "#E69F00", "C3" = "#009E73",
  "C4" = "#0072B2", "C5" = "#CC79A7", "Isolated" = "#BBBBBB"
)

# ── GRAPH: connected nodes only ───────────────────────────────
connected_v <- nodes[degree > 0, gene]
isolated_v  <- sort(nodes[degree == 0, gene])
nodes_conn  <- nodes[gene %in% connected_v]
edges_conn  <- edges_dt[gene_a %in% connected_v & gene_b %in% connected_v]

g <- graph_from_data_frame(
  d        = edges_conn[, .(from = gene_a, to = gene_b,
                             weight = jaccard, edge_type = source)],
  directed = FALSE,
  vertices = nodes_conn[, .(name = gene, degree, betweenness,
                             res_class, comm_lab, is_hub)]
)

# ── LAYOUT identical to Fig1 ──────────────────────────────────
set.seed(42)
w <- E(g)$weight
w[is.nan(w) | is.na(w)] <- median(w, na.rm = TRUE)
E(g)$weight <- w
pos <- layout_with_stress(g, weights = sqrt(E(g)$weight))
pos <- norm_coords(pos, xmin = -0.72, xmax = 0.72,
                   ymin = -0.72, ymax = 0.72)
lay <- create_layout(g, layout = "manual",
                     x = pos[, 1], y = pos[, 2])

# ── Community ELLIPSES (same criteria as Fig1) ────────────────
ell_data <- data.frame(x = lay$x, y = lay$y,
                       comm_lab = lay$comm_lab)
ell_data <- ell_data[ell_data$comm_lab != "Isolated", ]
comm_counts <- table(ell_data$comm_lab)
ell_data    <- ell_data[ell_data$comm_lab %in%
                          names(comm_counts[comm_counts >= 3]), ]

# ── ISOLATED RING identical to Fig1 ──────────────────────────
n_iso  <- length(isolated_v)
angles <- seq(pi/2, pi/2 + 2*pi, length.out = n_iso + 1)[-(n_iso + 1)]
r_iso  <- 1.05
iso_x  <- cos(angles) * r_iso
iso_y  <- sin(angles) * r_iso

# ── LEGEND: same parameters as Fig1 ──────────────────────────
lx0  <- -1.58
ly0  <- 1.42
lsep <- 0.100

# ── FIGURE ───────────────────────────────────────────────────
p <- ggraph(lay) +

  # Community ellipses: dotted border in community color
  # Same as Fig1 — reader can compare directly
  geom_mark_ellipse(
    data      = ell_data,
    aes(x = x, y = y, group = comm_lab, color = comm_lab),
    fill      = NA, linewidth = 0.60, linetype = "dashed",
    expand    = unit(7, "mm"),
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = comm_pal, guide = "none") +

  # Community labels at the periphery (C1-C5)
  # Allows reader to associate class color with community
  {
    centroids <- aggregate(cbind(x, y) ~ comm_lab, data = ell_data,
                           FUN = mean)
    gx <- mean(lay$x); gy <- mean(lay$y)
    centroids$dx  <- centroids$x - gx
    centroids$dy  <- centroids$y - gy
    centroids$len <- sqrt(centroids$dx^2 + centroids$dy^2)
    centroids$dx  <- centroids$dx / centroids$len
    centroids$dy  <- centroids$dy / centroids$len
    max_r <- sapply(centroids$comm_lab, function(cl) {
      sub <- ell_data[ell_data$comm_lab == cl, ]
      cx  <- centroids$x[centroids$comm_lab == cl]
      cy  <- centroids$y[centroids$comm_lab == cl]
      max(sqrt((sub$x - cx)^2 + (sub$y - cy)^2))
    })
    centroids$lx <- pmin(pmax(
      centroids$x + centroids$dx * (max_r + 0.18), -1.38), 1.38)
    centroids$ly <- pmin(pmax(
      centroids$y + centroids$dy * (max_r + 0.18), -1.38), 1.38)
    annotate("text",
             x = centroids$lx, y = centroids$ly,
             label    = centroids$comm_lab,
             color    = comm_pal[centroids$comm_lab],
             fontface = "bold", size = 3.2)
  } +

  # bnlearn-only edges
  geom_edge_link(
    aes(filter     = (edge_type != "Jaccard+bnlearn"),
        edge_alpha = weight * 0.35),
    color = "#C5C5C5", width = 0.16, show.legend = FALSE
  ) +

  # Consensus edges
  geom_edge_link(
    aes(filter      = (edge_type == "Jaccard+bnlearn"),
        edge_colour = weight, edge_alpha = weight),
    width = 0.85, show.legend = FALSE
  ) +
  scale_edge_colour_gradientn(
    colours = c("#C9E8F5","#6BAED6","#2171B5","#084594","#03224C"),
    limits  = c(0.3, 1.0), guide = "none"
  ) +
  scale_edge_alpha(range = c(0.18, 0.88), guide = "none") +

  # ALL nodes as circles (no diamond) — color by class
  geom_node_point(
    aes(filter = (degree > 0),
        size   = 2.0 + degree * 0.38,
        fill   = res_class),
    shape = 21, color = "white", stroke = 0.55, alpha = 0.95,
    show.legend = TRUE
  ) +

  scale_fill_manual(
    values = class_pal,
    name   = "Resistance class",
    guide  = guide_legend(
      override.aes = list(shape = 21, size = 4.5,
                          stroke = 0.5, alpha = 1),
      ncol         = 1,
      keyheight    = unit(0.40, "cm"),
      keywidth     = unit(0.40, "cm"),
      label.theme  = element_text(size = 8.0, face = "italic",
                                  color = "#222222")
    )
  ) +
  scale_size_identity() +

  # Non-hub gene labels — size 2.3 ~ 8.3 pt
  geom_node_text(
    aes(filter   = (degree > 0 & is_hub == 0),
        label    = name, color = res_class),
    size = 2.3, repel = TRUE, bg.color = "white", bg.r = 0.10,
    max.overlaps = 15, seed = 42, fontface = "plain",
    show.legend = FALSE,
    point.padding = unit(0.18, "lines"),
    box.padding   = unit(0.08, "lines"),
    force = 0.25, force_pull = 1.5
  ) +
  geom_node_text(
    aes(filter = (is_hub == 1), label = name, color = res_class),
    size = 2.8, fontface = "bold", repel = TRUE,
    bg.color = "white", bg.r = 0.12,
    max.overlaps = 15, seed = 42, show.legend = FALSE,
    point.padding = unit(0.32, "lines"),
    box.padding   = unit(0.10, "lines"),
    force = 0.35, force_pull = 1.5
  ) +
  scale_color_manual(values = class_pal, guide = "none") +

  # ── ISOLATED NODES ────────────────────────────────────────────
  annotate("point", x = iso_x, y = iso_y,
           shape = 21, size = 2.2,
           fill = "#AAAAAA", color = "white", stroke = 0.45) +
  annotate("text",
           x = iso_x * 1.10, y = iso_y * 1.10,
           label    = isolated_v,
           size     = 2.5, fontface = "bold", color = "#666666",
           hjust    = ifelse(iso_x >= 0, 0, 1)) +

  # ── MORPHOLOGY LEGEND — upper left corner ─────────────────────
  annotate("text", x = lx0, y = ly0,
           label = "Node type", hjust = 0, fontface = "bold",
           size = 2.9, color = "#222222") +
  annotate("point", x = lx0 + 0.06, y = ly0 - lsep,
           shape = 21, size = 3.8, fill = "#777777",
           color = "white", stroke = 0.6) +
  annotate("text", x = lx0 + 0.16, y = ly0 - lsep,
           label = "Gene", hjust = 0, size = 2.6, color = "#333333") +
  annotate("segment", x = lx0 + 0.04, xend = lx0 + 0.18,
           y = ly0 - 2*lsep, yend = ly0 - 2*lsep,
           color = "#2171B5", linewidth = 0.8) +
  annotate("text", x = lx0 + 0.22, y = ly0 - 2*lsep,
           label = "Jaccard ∩ bnlearn", hjust = 0,
           size = 2.6, color = "#444444") +
  annotate("segment", x = lx0 + 0.04, xend = lx0 + 0.18,
           y = ly0 - 3*lsep, yend = ly0 - 3*lsep,
           color = "#CCCCCC", linewidth = 0.4) +
  annotate("text", x = lx0 + 0.22, y = ly0 - 3*lsep,
           label = "bnlearn only", hjust = 0,
           size = 2.6, color = "#999999") +

  # ── JACCARD — lower left ──────────────────────────────────────
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
  annotate("rect", xmin = lx0, xmax = lx0 + 0.59,
           ymin = -1.33, ymax = -1.27,
           fill = NA, color = "#AAAAAA", linewidth = 0.3) +
  annotate("text", x = lx0,        y = -1.38, label = "0.3",
           hjust = 0.5, size = 2.3, color = "#555555") +
  annotate("text", x = lx0 + 0.30, y = -1.38, label = "0.65",
           hjust = 0.5, size = 2.3, color = "#555555") +
  annotate("text", x = lx0 + 0.59, y = -1.38, label = "1.0",
           hjust = 0.5, size = 2.3, color = "#555555") +

  coord_fixed(xlim = c(-1.65, 1.45), ylim = c(-1.50, 1.50)) +
  labs(title = NULL, subtitle = NULL) +
  theme_graph(base_family = "sans") +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "right",
    legend.background = element_rect(fill      = "#FAFAFA",
                                     color     = "#E0E0E0",
                                     linewidth = 0.3),
    legend.margin     = margin(6, 8, 6, 8),
    legend.title      = element_text(size = 9.0, face = "bold",
                                     color = "#222222"),
    legend.key.size   = unit(0.42, "cm"),
    plot.margin       = margin(8, 8, 8, 8)
  )

# ── EXPORT PNG + TIFF ─────────────────────────────────────────
png_path  <- file.path(out_dir, "FigS1_network_resistance_class.png")
tiff_path <- file.path(out_dir, "FigS1_network_resistance_class.tiff")

ggsave(png_path,  p, width = 14, height = 11, dpi = 300, bg = "white")
ggsave(tiff_path, p, width = 14, height = 11, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("OK: PNG :", png_path,  "\n")
cat("OK: TIFF:", tiff_path, "\n")
