library(ggplot2)
library(data.table)
library(patchwork)

base <- "/home/miguel/Abaum_Resistome_Network"
prev <- fread(file.path(base, "results/tables/temporal_prevalence_by_year.csv"))
topo <- fread(file.path(base, "results/tables/temporal_topology.csv"))

prev <- merge(prev, topo[, .(year, n_genomes)], on="year")
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
for(g in genes_up) {
  df <- data.frame(year = prev$year, prev = prev[[g]])
  mod <- lm(prev ~ year, data = df)
  df$fitted   <- fitted(mod)
  df$residual <- residuals(mod)
  df$std_res  <- rstandard(mod)
  pv   <- summary(mod)$coefficients[2,4]
  plab <- ifelse(pv < 0.001, "p<0.001", sprintf("p=%.3f", pv))

  p <- ggplot(df, aes(x=fitted, y=std_res)) +
    geom_hline(yintercept=0, linetype="dashed",
               color="gray50", linewidth=0.5) +
    geom_hline(yintercept=c(-2,2), linetype="dotted",
               color="red", linewidth=0.4) +
    geom_point(color="#2E4057", size=2.5, alpha=0.8) +
    geom_smooth(method="loess", se=TRUE, color="#E84855",
                fill="#E84855", alpha=0.15, linewidth=0.8) +
    labs(
      title    = labels[g],
      subtitle = plab,
      x        = "Fitted values",
      y        = "Standardized residuals"
    ) +
    theme_bw(base_size=10) +
    theme(
      plot.title       = element_text(face="bold.italic", size=9),
      plot.subtitle    = element_text(size=8, color="gray40"),
      panel.grid.minor = element_blank()
    )
  plots[[g]] <- p
}

fig <- wrap_plots(plots, ncol=3)

out_png  <- file.path(base,
  "results/figures/supplementary/FigS5_regression_residuals.png")
out_tiff <- file.path(base,
  "results/figures/supplementary/FigS5_regression_residuals.tiff")

ggsave(out_png,  fig, width=12, height=8, dpi=300)
ggsave(out_tiff, fig, width=12, height=8, dpi=300, compression="lzw")

cat("Saved:\n", out_png, "\n", out_tiff, "\n")
