suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ── PATHS ────────────────────────────────────────────────────
base    <- "/home/miguel/Abaum_Resistome_Network"
tables  <- file.path(base, "results/tables")
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

conv <- fread(file.path(tables, "one_health_convergence.csv"))

# ── TOP 20 by convergence score ───────────────────────────────
top20 <- conv[convergence_score > 0][order(-convergence_score)][1:20]
top20[, gene_f := factor(gene, levels = rev(gene))]

# ── LONG FORMAT ───────────────────────────────────────────────
cl <- melt(
  top20[, .(gene_f, human_prev, animal_prev, env_prev)],
  id.vars       = "gene_f",
  variable.name = "comp",
  value.name    = "prev"
)
cl[, comp := factor(
  comp,
  levels = c("human_prev", "animal_prev", "env_prev"),
  labels = c("Human", "Animal", "Environment")
)]

# Labels: all bars >= 3%
cl[, lab := ifelse(prev > 0, paste0(round(prev * 100, 0), "%"), "")]

# ── One Health PALETTE ────────────────────────────────────────
comp_pal <- c(
  "Human"       = "#2166AC",
  "Animal"      = "#D6604D",
  "Environment" = "#4DAC26"
)

# ── ALTERNATING BACKGROUND BANDS for row readability ──────────
n_genes  <- nlevels(cl$gene_f)
band_df  <- data.frame(
  ymin = seq(0.5, n_genes - 0.5, by = 2),
  ymax = seq(1.5, n_genes + 0.5, by = 2)
)

# ── FIGURE ───────────────────────────────────────────────────
p <- ggplot(cl, aes(x = prev, y = gene_f, fill = comp)) +

  # Soft alternating bands — improves readability without distraction
  geom_rect(
    data        = band_df,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = Inf),
    fill        = "#F7F7F7",
    color       = NA,
    inherit.aes = FALSE
  ) +

  # Bars
  geom_col(
    position = position_dodge(0.78),
    width    = 0.72,
    alpha    = 0.93
  ) +

  # Labels always outside the bar — avoids clipping in short bars
  geom_text(
    aes(label = lab,
        color = "#222222"),
    position    = position_dodge(0.78),
    hjust       = -0.12,
    size        = 2.85,
    fontface    = "bold",
    show.legend = FALSE
  ) +
  scale_color_identity() +

  # 10% reference line — minimum epidemiological relevance threshold
  geom_vline(
    xintercept = 0.10,
    linetype   = "dashed",
    color      = "#999999",
    linewidth  = 0.50
  ) +
  annotate(
    "text", x = 0.103, y = 0.55,
    label    = "10%",
    size     = 2.85, hjust = 0,
    color    = "#999999", fontface = "italic"
  ) +

  scale_fill_manual(
    values = comp_pal,
    name   = "Compartment",
    guide  = guide_legend(
      keywidth     = unit(0.55, "cm"),
      keyheight    = unit(0.42, "cm"),
      override.aes = list(alpha = 1),
      label.theme  = element_text(size = 9.5)
    )
  ) +

  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.18),
    breaks = seq(0, 1.0, 0.25),
    expand = expansion(mult = c(0, 0.01))
  ) +

  scale_y_discrete(expand = expansion(add = c(0.5, 0.7))) +

  labs(
    x = "Mean prevalence\n(Human: n = 6,446   │   Animal: n = 284   │   Environment: n = 397)",
    y = NULL
  ) +

  theme_minimal(base_size = 11) +
  theme(
    # Legend on top, centered, compact
    legend.position    = "top",
    legend.justification = "center",
    legend.title       = element_text(size = 10, face = "bold",
                                      color = "#222222"),
    legend.text        = element_text(size = 9.5, color = "#333333"),
    legend.key.spacing = unit(0.6, "cm"),
    legend.box.spacing = unit(0.15, "cm"),
    legend.margin      = margin(b = 4),

    # Grid
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor.x = element_blank(),

    # Axes
    axis.text.y  = element_text(size = 9.5, face = "italic",
                                color = "#111111"),
    axis.text.x  = element_text(size = 9.5, color = "#333333"),
    axis.title.x = element_text(size = 9.0, color = "#555555",
                                margin    = margin(t = 8),
                                lineheight = 1.5),

    # Backgrounds
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Margins
    plot.margin = margin(t = 12, r = 20, b = 12, l = 12)
  )

# ── EXPORT PNG + TIFF ─────────────────────────────────────────
png_path  <- file.path(out_dir, "Fig3_one_health_convergence.png")
tiff_path <- file.path(out_dir, "Fig3_one_health_convergence.tiff")
ggsave(png_path,  p, width = 13, height = 10, dpi = 300, bg = "white")
ggsave(tiff_path, p, width = 13, height = 10, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("OK: PNG :", png_path,  "\n")
cat("OK: TIFF:", tiff_path, "\n")
