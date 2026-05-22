suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ── PATHS ────────────────────────────────────────────────────
base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

geo <- fread(file.path(tables, "geo_regional_burden.csv"))
geo <- geo[order(mean_args)]
geo[, region_f := factor(region, levels = region)]

# N label with correct grammar
geo[, n_lab := paste0(
  "n = ", format(n_genomes, big.mark = ","),
  "  │  ", n_countries,
  ifelse(n_countries == 1, " country", " countries")
)]

# ── GLOBAL MEAN ───────────────────────────────────────────────
mean_global <- mean(geo$mean_args)

# ── ALTERNATING BANDS ─────────────────────────────────────────
n_reg   <- nrow(geo)
band_df <- data.frame(
  ymin = seq(0.5, n_reg - 0.5, by = 2),
  ymax = seq(1.5, n_reg + 0.5, by = 2)
)

# ── REGION PALETTE: sequential burden scale ───────────────────
# Color encodes burden (mean_args) — from lowest to highest
# Blue-teal-orange palette: perceptually uniform, without saturation
region_cols <- colorRampPalette(
  c("#4575B4", "#74ADD1", "#ABD9E9", "#FEE090", "#FDAE61", "#F46D43", "#D73027")
)(n_reg)
names(region_cols) <- levels(geo$region_f)

# ── GENE PALETTE ──────────────────────────────────────────────
gene_pal <- c(
  "blaOXA-23-like" = "#1B7837",   # dark green — OXA carbapenemase
  "blaNDM"         = "#762A83",   # purple — metallo-carbapenemase
  "gyrA"           = "#E08214"    # orange — quinolone QRDR
)

# ── PANEL A: BURDEN by region ─────────────────────────────────
p_left <- ggplot(geo, aes(x = mean_args, y = region_f)) +

  # Alternating background bands
  geom_rect(
    data        = band_df,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = Inf),
    fill        = "#F5F5F5", color = NA,
    inherit.aes = FALSE
  ) +

  geom_col(aes(fill = region_f), width = 0.72, alpha = 0.92) +
  scale_fill_manual(values = region_cols, guide = "none") +

  # Numeric value outside the bar
  geom_text(
    aes(label = sprintf("%.2f", mean_args)),
    hjust    = -0.12,
    size     = 3.1,
    fontface = "bold",
    color    = "#222222"
  ) +

  # N and countries — inside the bar, near the right of origin
  geom_text(
    aes(x = 0.25, label = n_lab),
    hjust    = 0,
    size     = 2.85,
    color    = "white",
    fontface = "italic"
  ) +

  # Global mean line
  geom_vline(
    xintercept = mean_global,
    linetype   = "dashed",
    color      = "#555555",
    linewidth  = 0.65
  ) +
  annotate(
    "text",
    x = mean_global + 0.18, y = 0.55,
    label    = sprintf("Global mean\n%.2f ARGs", mean_global),
    size     = 2.85, hjust = 0,
    color    = "#555555", fontface = "italic"
  ) +

  scale_x_continuous(
    limits = c(0, 16),
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +

  labs(x = "Mean ARG families per genome", y = NULL) +

  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10, face = "bold",
                                      color = "#111111"),
    axis.text.x        = element_text(size = 9,  color = "#333333"),
    axis.title.x       = element_text(size = 9.5, color = "#333333",
                                      margin = margin(t = 6)),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 4, 10, 12)
  )

# ── PANEL B: key gene prevalence ─────────────────────────────
geo_long <- melt(
  geo[, .(region_f, blaoxa_23like, blandm, gyra)],
  id.vars       = "region_f",
  variable.name = "gene",
  value.name    = "prev"
)
geo_long[, gene := factor(
  gene,
  levels = c("gyra", "blaoxa_23like", "blandm"),
  labels = c("gyrA", "blaOXA-23-like", "blaNDM")
)]
geo_long[, lab := ifelse(prev >= 0.02,
                          paste0(round(prev * 100, 0), "%"), "")]

p_right <- ggplot(geo_long, aes(x = prev, y = region_f, fill = gene)) +

  # Alternating bands — aligned with left panel
  geom_rect(
    data        = band_df,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = Inf),
    fill        = "#F5F5F5", color = NA,
    inherit.aes = FALSE
  ) +

  geom_col(
    position = position_dodge(0.78),
    width    = 0.70,
    alpha    = 0.92
  ) +

  geom_text(
    aes(label = lab),
    position  = position_dodge(0.78),
    hjust     = -0.10,
    size      = 2.45,
    fontface  = "bold",
    color     = "#222222"
  ) +

  # 50% reference line
  geom_vline(
    xintercept = 0.50,
    linetype   = "dashed",
    color      = "#BBBBBB",
    linewidth  = 0.50
  ) +
  annotate("text", x = 0.515, y = 0.55,
           label = "50%", size = 2.85, hjust = 0,
           color = "#BBBBBB", fontface = "italic") +

  scale_fill_manual(
    values = gene_pal,
    name   = "Resistance marker",
    guide  = guide_legend(
      keywidth  = unit(0.5, "cm"),
      keyheight = unit(0.4, "cm"),
      label.theme = element_text(size = 9, face = "italic")
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.22),
    breaks = seq(0, 1.0, 0.25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +

  labs(x = "Prevalence", y = NULL) +

  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.justification = "left",
    legend.title       = element_text(size = 9.5, face = "bold",
                                      color = "#222222"),
    legend.text        = element_text(size = 9.0),
    legend.key.spacing = unit(0.5, "cm"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 9,  color = "#333333"),
    axis.title.x       = element_text(size = 9.5, color = "#333333",
                                      margin = margin(t = 6)),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 14, 10, 2)
  )

# ── COMPOSITION ───────────────────────────────────────────────
final <- (p_left | p_right) +
  plot_layout(widths = c(1.25, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# ── EXPORT PNG + TIFF ─────────────────────────────────────────
png_path  <- file.path(out_dir, "Fig4_geo_regional_burden.png")
tiff_path <- file.path(out_dir, "Fig4_geo_regional_burden.tiff")
ggsave(png_path,  final, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(tiff_path, final, width = 14, height = 7, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("OK: PNG :", png_path,  "\n")
cat("OK: TIFF:", tiff_path, "\n")
