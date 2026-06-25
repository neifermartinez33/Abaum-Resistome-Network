suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
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

base <- "/home/miguel/Abaum_Resistome_Network"
prev <- fread(file.path(base, "results/tables/temporal_prevalence_by_year.csv"))
topo <- fread(file.path(base, "results/tables/temporal_topology.csv"))

prev <- merge(prev, topo[, .(year, n_genomes)], on = "year")
prev <- prev[n_genomes >= 100]

genes_up <- c("ftsi", "arma", "blandm", "blaoxa_23like", "msr(e)", "mph(e)")

labels <- c(
  "ftsi"          = "ftsI (R²=0.830)",
  "arma"          = "armA (R²=0.633)",
  "blandm"        = "blaNDM (R²=0.762)",
  "blaoxa_23like" = "blaOXA-23-like (R²=0.432)",
  "msr(e)"        = "msr(E) (R²=0.576)",
  "mph(e)"        = "mph(E) (R²=0.553)"
)

plots <- list()
for (g in genes_up) {
  df <- data.frame(year = prev$year, prev = prev[[g]])
  mod <- lm(prev ~ year, data = df)
  df$fitted  <- fitted(mod)
  df$std_res <- rstandard(mod)
  pv   <- summary(mod)$coefficients[2, 4]
  plab <- ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv))

  p <- ggplot(df, aes(x = fitted, y = std_res)) +

    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted",
               color = "#CC0000", linewidth = 0.5) +

    geom_smooth(method = "loess", se = TRUE,
                color = "#E84855", fill = "#E84855",
                alpha = 0.15, linewidth = 1.0) +

    geom_point(color = "#2E4057", size = 3.2, alpha = 0.85) +

    labs(
      title    = labels[g],
      subtitle = plab,
      x        = "Fitted values",
      y        = "Standardized residuals"
    ) +

    theme_bw(base_size = 22, base_family = FONT) +
    theme(
      plot.title       = element_text(face = "bold.italic", size = 19,
                                      family = FONT, color = "#111111"),
      plot.subtitle    = element_text(size = 17, family = FONT,
                                      color = "gray40",
                                      margin = margin(b = 4)),
      axis.title       = element_text(size = 18.0, family = FONT,
                                      color = "#333333"),
      axis.text        = element_text(size = 17.0, family = FONT,
                                      color = "#444444"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.4),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(t = 10, r = 14, b = 8, l = 8)
    )

  plots[[g]] <- p
}

fig <- wrap_plots(plots, ncol = 3) &
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png  <- file.path(base,
  "results/figures/supplementary/FigS4_regression_residuals.png")
out_tiff <- file.path(base,
  "results/figures/supplementary/FigS4_regression_residuals.tiff")

ggsave(out_png,  fig, width = 14, height = 9, dpi = 600, bg = "white")
ggsave(out_tiff, fig, width = 14, height = 9, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", out_png,  "\n")
cat("\u2713 TIFF:", out_tiff, "\n")
