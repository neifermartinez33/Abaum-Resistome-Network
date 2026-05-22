cat > /home/miguel/Abaum_Resistome_Network/scripts/r/figS2_topology_v2.R << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

topo <- fread(file.path(tables, "temporal_topology.csv"))
topo <- topo[year >= 2003 & year <= 2024]

theme_s2 <- function(show_x = FALSE) {
  t <- theme_minimal(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 9, color = "#333333"),
      axis.title.y     = element_text(size = 10, color = "#222222", margin = margin(r = 6)),
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#FAFAFA", color = NA),
      plot.margin      = margin(4, 14, 2, 14),
      legend.position  = "none"
    )
  if (show_x) {
    t <- t + theme(
      axis.text.x  = element_text(size = 9, color = "#333333", angle = 45, hjust = 1),
      axis.title.x = element_text(size = 10, color = "#222222", margin = margin(t = 6))
    )
  } else {
    t <- t + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
  }
  t
}

col_density    <- "#0072B2"
col_clustering <- "#009E73"
col_degree     <- "#CC6600"
col_genomes    <- "#999999"
xpan_min <- 2019.8
xpan_max <- 2021.5
xrng     <- c(2002.5, 2024.5)
base_clustering <- topo[year == 2003, clustering]
base_degree     <- topo[year == 2003, mean_degree]

p1 <- ggplot(topo, aes(x = year, y = density)) +
  annotate("rect", xmin = xpan_min, xmax = xpan_max, ymin = -Inf, ymax = Inf, fill = "#FFF3CD", color = NA, alpha = 0.7) +
  geom_area(fill = alpha(col_density, 0.10), color = NA) +
  geom_smooth(method = "loess", se = TRUE, span = 0.55, color = alpha(col_density, 0.50), fill = alpha(col_density, 0.12), linewidth = 0.65) +
  geom_line(color = col_density, linewidth = 0.90) +
  geom_point(fill = "white", color = col_density, shape = 21, size = 2.5, stroke = 1.1) +
  annotate("text", x = xpan_min + 0.1, y = 0.060, label = "Pandemic", size = 2.3, hjust = 0, color = "#AA8800", fontface = "italic") +
  scale_x_continuous(breaks = seq(2003, 2024, 3), limits = xrng) +
  scale_y_continuous(limits = c(0.030, 0.063), labels = number_format(accuracy = 0.001)) +
  theme_s2(show_x = FALSE) +
  labs(y = "Network density")

p2 <- ggplot(topo, aes(x = year, y = clustering)) +
  annotate("rect", xmin = xpan_min, xmax = xpan_max, ymin = -Inf, ymax = Inf, fill = "#FFF3CD", color = NA, alpha = 0.7) +
  geom_area(fill = alpha(col_clustering, 0.10), color = NA) +
  geom_smooth(method = "loess", se = TRUE, span = 0.55, color = alpha(col_clustering, 0.50), fill = alpha(col_clustering, 0.12), linewidth = 0.65) +
  geom_line(color = col_clustering, linewidth = 0.90) +
  geom_point(fill = "white", color = col_clustering, shape = 21, size = 2.5, stroke = 1.1) +
  geom_hline(yintercept = base_clustering, linetype = "dotted", color = "#AAAAAA", linewidth = 0.5) +
  annotate("text", x = 2002.6, y = base_clustering + 0.012, label = "2003 baseline", size = 2.1, hjust = 0, color = "#AAAAAA", fontface = "italic") +
  scale_x_continuous(breaks = seq(2003, 2024, 3), limits = xrng) +
  scale_y_continuous(limits = c(0.18, 0.54), labels = number_format(accuracy = 0.01)) +
  theme_s2(show_x = FALSE) +
  labs(y = "Clustering\ncoefficient")

p3 <- ggplot(topo, aes(x = year, y = mean_degree)) +
  annotate("rect", xmin = xpan_min, xmax = xpan_max, ymin = -Inf, ymax = Inf, fill = "#FFF3CD", color = NA, alpha = 0.7) +
  geom_area(fill = alpha(col_degree, 0.10), color = NA) +
  geom_smooth(method = "loess", se = TRUE, span = 0.55, color = alpha(col_degree, 0.50), fill = alpha(col_degree, 0.12), linewidth = 0.65) +
  geom_line(color = col_degree, linewidth = 0.90) +
  geom_point(fill = "white", color = col_degree, shape = 21, size = 2.5, stroke = 1.1) +
  geom_hline(yintercept = base_degree, linetype = "dotted", color = "#AAAAAA", linewidth = 0.5) +
  annotate("text", x = 2002.6, y = base_degree + 0.065, label = "2003 baseline", size = 2.1, hjust = 0, color = "#AAAAAA", fontface = "italic") +
  scale_x_continuous(breaks = seq(2003, 2024, 3), limits = xrng) +
  scale_y_continuous(limits = c(1.20, 2.45), labels = number_format(accuracy = 0.01)) +
  theme_s2(show_x = FALSE) +
  labs(y = "Mean degree")

p4 <- ggplot(topo, aes(x = year, y = n_genomes)) +
  annotate("rect", xmin = xpan_min, xmax = xpan_max, ymin = -Inf, ymax = Inf, fill = "#FFF3CD", color = NA, alpha = 0.7) +
  geom_col(fill = col_genomes, alpha = 0.80, width = 0.75) +
  geom_text(aes(label = format(n_genomes, big.mark = ",")), vjust = -0.35, size = 1.85, color = "#555555", angle = 90, hjust = 0) +
  scale_x_continuous(breaks = seq(2003, 2024, 3), limits = xrng) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.32))) +
  theme_s2(show_x = TRUE) +
  theme(plot.margin = margin(2, 14, 10, 14)) +
  labs(x = "Collection year", y = "Genomes (n)")

final <- (p1 / p2 / p3 / p4) +
  plot_layout(heights = c(1, 1, 1, 0.65)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

out_path <- file.path(out_dir, "FigS2_temporal_network_topology.png")
ggsave(out_path, final, width = 11, height = 13, dpi = 300, bg = "white")
cat("OK: FigS2 saved at:", out_path, "\n")
EOF
