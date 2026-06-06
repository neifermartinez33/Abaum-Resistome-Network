suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
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
out_dir <- file.path(base, "results/figures/supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

topo <- fread(file.path(tables, "temporal_topology.csv"))
topo <- topo[year >= 2003 & year <= 2024]

# RELATIVE CHANGE vs 2003 BASELINE
base_density    <- topo[year == 2003, density]
base_clustering <- topo[year == 2003, clustering]
base_degree     <- topo[year == 2003, mean_degree]

topo[, delta_density    := (density     - base_density)    / base_density    * 100]
topo[, delta_clustering := (clustering  - base_clustering) / base_clustering * 100]
topo[, delta_degree     := (mean_degree - base_degree)     / base_degree     * 100]

# LONG FORMAT for heatmap
metrics <- c("Network density", "Clustering coefficient", "Mean degree")

long <- rbind(
  topo[, .(year, metric = "Network density",       delta = delta_density,
           raw = density,    raw_lab = sprintf("%.3f", density))],
  topo[, .(year, metric = "Clustering coefficient",delta = delta_clustering,
           raw = clustering, raw_lab = sprintf("%.2f", clustering))],
  topo[, .(year, metric = "Mean degree",           delta = delta_degree,
           raw = mean_degree,raw_lab = sprintf("%.2f", mean_degree))]
)
long[, metric := factor(metric, levels = rev(metrics))]

# Clamp delta for colour scale (-40 to +40 %)
long[, delta_clamp := pmax(pmin(delta, 40), -40)]

# Text colour: white on extreme cells, dark on light
long[, txt_col := ifelse(abs(delta_clamp) > 22, "white", "#222222")]

# ALTERNATING BANDS
years <- sort(unique(long$year))
n_yr  <- length(years)
band_df <- data.frame(
  xmin = years[seq(1, n_yr, 2)] - 0.5,
  xmax = years[seq(1, n_yr, 2)] + 0.5
)

# PANEL A: HEATMAP
p_heat <- ggplot(long, aes(x = year, y = metric)) +

  # Alternating column bands
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "#F7F7F7", color = NA, inherit.aes = FALSE
  ) +

  # Heatmap tiles — divergent palette centred at 0
  geom_tile(
    aes(fill = delta_clamp),
    color = "white", linewidth = 0.8,
    width = 0.92, height = 0.88
  ) +
  scale_fill_gradientn(
    colours  = c("#2166AC","#4393C3","#92C5DE","#D1E5F0",
                 "#F7F7F7",
                 "#FDDBC7","#F4A582","#D6604D","#B2182B"),
    limits   = c(-40, 40),
    breaks   = c(-40, -20, 0, 20, 40),
    labels   = c("-40%", "-20%", "0%", "+20%", "+40%"),
    name     = "Change vs.\n2003 baseline",
    guide    = guide_colourbar(
      barwidth       = unit(0.5, "cm"),
      barheight      = unit(4.0, "cm"),
      title.theme    = element_text(size = 18.0, family = FONT,
                                    face = "bold", lineheight = 1.2,
                                    color = "#222222"),
      label.theme    = element_text(size = 16.0, family = FONT,
                                    color = "#444444"),
      ticks.linewidth = 0.5,
      frame.colour   = "#AAAAAA",
      frame.linewidth = 0.4
    )
  ) +

  # Raw value labels inside tiles — rotated 90° to prevent overlap across 22 years
  geom_text(
    aes(label = raw_lab, color = txt_col),
    size     = 5.8, fontface = "bold", family = FONT,
    angle    = 90
  ) +
  scale_color_identity() +

  # Baseline 2003 marker
  annotate("rect",
           xmin = 2002.54, xmax = 2003.46,
           ymin = 0.52, ymax = 3.48,
           fill = NA, color = "#555555",
           linewidth = 1.2) +
  annotate("text",
           x = 2003, y = 3.62,
           label = "Baseline",
           size = 5.8, hjust = 0.5, family = FONT,
           color = "#555555", fontface = "italic") +

  scale_x_continuous(
    breaks = seq(2003, 2024, 3),
    expand = expansion(add = c(0.55, 0.55))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.55, 0.75))) +

  labs(x = NULL, y = NULL) +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 19.0, family = FONT,
                                     color = "#111111", face = "bold",
                                     hjust = 1),
    legend.position   = "right",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.margin       = margin(t = 16, r = 8, b = 4, l = 8)
  )

# PANEL B: N GENOMES BAR
p_bar <- ggplot(topo, aes(x = year, y = n_genomes)) +

  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "#F7F7F7", color = NA, inherit.aes = FALSE
  ) +

  geom_col(
    fill = "#888888", alpha = 0.75, width = 0.82
  ) +

  # n labels above bars — rotated 90° to prevent crowding across 22 years
  geom_text(
    aes(label = format(n_genomes, big.mark = ",")),
    vjust  = 0.4,
    hjust  = -0.1,
    angle  = 90,
    family = FONT, color = "#555555", fontface = "plain", size = 5.8
  ) +

  scale_x_continuous(
    breaks = seq(2003, 2024, 3),
    expand = expansion(add = c(0.55, 0.55))
  ) +
  scale_y_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.38))
  ) +

  labs(x = "Collection year", y = "Genomes (n)") +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 18.0, family = FONT,
                                      color = "#333333",
                                      angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 17.0, family = FONT,
                                      color = "#333333"),
    axis.title.x       = element_text(size = 19.0, family = FONT,
                                      color = "#222222",
                                      margin = margin(t = 8)),
    axis.title.y       = element_text(size = 18.0, family = FONT,
                                      color = "#222222",
                                      margin = margin(r = 6)),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 4, r = 8, b = 12, l = 8)
  )

# COMPOSITION
final <- (p_heat / p_bar) +
  plot_layout(heights = c(2.2, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# EXPORT
png_path  <- file.path(out_dir, "FigS1_temporal_network_topology.png")
tiff_path <- file.path(out_dir, "FigS1_temporal_network_topology.tiff")
ggsave(png_path,  final, width = 14, height = 8, dpi = 600, bg = "white")
ggsave(tiff_path, final, width = 14, height = 8, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_path,  "\n")
cat("\u2713 TIFF:", tiff_path, "\n")
