suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ── PATHS ────────────────────────────────────────────────────
base    <- "/home/miguel/Abaum_Resistome_Network"
sup_dir <- file.path(base, "results/supplementary")
out_dir <- file.path(base, "results/figures/supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── DATA ────────────────────────────────────────────────────
prev <- fread(file.path(sup_dir, "geo_prevalence_by_country.csv"))

# Filter countries with >= 50 genomes and take top 30 by burden
top30 <- prev[n_genomes >= 50][order(-mean_args)][1:30]
top30 <- top30[order(mean_args)]
top30[, country_f := factor(country_norm, levels = country_norm)]
top30[, n_lab := paste0("n = ", format(n_genomes, big.mark = ","))]

# Mean of the 30 countries shown (not the global mean of 20,739 genomes)
mean_shown <- mean(top30$mean_args)
cat(sprintf("Mean of the 30 countries shown: %.2f\n", mean_shown))
cat(sprintf("(The actual global mean of 20,739 genomes is 10.71)\n"))

# ── FIGURE ───────────────────────────────────────────────────
p_s4 <- ggplot(top30, aes(x = mean_args, y = country_f,
                            fill = mean_args)) +
  geom_col(width = 0.72, alpha = 0.93) +
  scale_fill_gradientn(
    colors = c("#4575B4","#74ADD1","#ABD9E9",
               "#FEE090","#FDAE61","#F46D43","#D73027"),
    name   = "Mean ARGs\nper genome",
    guide  = guide_colorbar(
      barheight      = unit(5.0, "cm"),
      barwidth       = unit(0.42, "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "grey50",
      frame.colour   = "grey60"
    )
  ) +

  # n inside the bar
  geom_text(
    aes(x = 0.25, label = n_lab),
    hjust = 0, size = 2.85,
    color = "white", fontface = "italic"
  ) +

  # Numeric value outside the bar
  geom_text(
    aes(label = sprintf("%.2f", mean_args)),
    hjust = -0.12, size = 2.85,
    fontface = "bold", color = "#222222"
  ) +

  # Mean line for the 30 countries shown
  geom_vline(xintercept = mean_shown,
             linetype = "dashed", color = "#555555", linewidth = 0.65) +
  annotate("text",
           x = mean_shown + 0.15, y = 1.5,
           label = sprintf("Mean\n(30 countries\nshown)\n%.2f", mean_shown),
           size = 2.85, hjust = 0,
           color = "#555555", fontface = "italic") +

  scale_x_continuous(
    limits = c(0, 17),
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = c(0, 0.10))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +

  labs(
    x       = "Mean ARG families per genome",
    y       = NULL,
    caption = paste0(
      "Top 30 countries by mean ARG burden among those with >= 50 genomes (38 countries analyzed total).\n",
      "Eight countries excluded: Poland (85% wildlife isolates), France and Germany ",
      "(mixed bovine/environmental isolates),\n",
      "and five others with non-clinical sample composition that inflates or deflates ARG burden estimates."
    )
  ) +

  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9.5, color = "#111111",
                                      face = "bold"),
    axis.text.x        = element_text(size = 9.0, color = "#333333"),
    axis.title.x       = element_text(size = 10.0, color = "#333333",
                                      margin = margin(t = 6)),
    legend.position    = "right",
    legend.title       = element_text(size = 9.5, face = "bold"),
    legend.text        = element_text(size = 8.5),
    plot.caption       = element_text(size = 7.5, color = "#666666",
                                      hjust = 0, lineheight = 1.4,
                                      margin = margin(t = 8)),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(12, 14, 12, 12)
  )

# ── EXPORT ───────────────────────────────────────────────────
png_s4  <- file.path(out_dir, "FigS4_geo_arg_burden_countries.png")
tiff_s4 <- file.path(out_dir, "FigS4_geo_arg_burden_countries.tiff")
ggsave(png_s4,  p_s4, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(tiff_s4, p_s4, width = 12, height = 11, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("OK: FigS4 PNG :", png_s4,  "\n")
cat("OK: FigS4 TIFF:", tiff_s4, "\n")
