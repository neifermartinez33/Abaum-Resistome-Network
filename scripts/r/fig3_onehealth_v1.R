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
out_dir <- file.path(base, "results/figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

conv <- fread(file.path(tables, "one_health_convergence.csv"))

# TOP 20 ordenados por human_prev (eje Y más informativo)
top20 <- conv[convergence_score > 0][order(-convergence_score)][1:20]
top20 <- top20[order(human_prev)]
top20[, gene_f := factor(gene, levels = gene)]

# CLASE DE RESISTENCIA (para color del convergence score)
res_class <- c(
  "msr(E)"        = "Macrolide",
  "armA"          = "Aminoglycoside",
  "mph(E)"        = "Macrolide",
  "ftsI"          = "Beta-lactam",
  "blaOXA-23-like" = "Carbapenem",
  "tet(B)"        = "Tetracycline",
  "blaNDM"        = "Carbapenem",
  "blaOXA-24-like" = "Carbapenem",
  "arr-2"         = "Rifampicin",
  "pmrB"          = "Colistin",
  "aac(3)"        = "Aminoglycoside",
  "ant(2'')"      = "Aminoglycoside",
  "qacEΔ1"        = "Disinfectant",
  "sul"           = "Sulfonamide",
  "aph(3')"       = "Aminoglycoside",
  "aadA"          = "Aminoglycoside",
  "adeC"          = "Efflux",
  "blaTEM"        = "Beta-lactam",
  "merT"          = "Metal",
  "merR"          = "Metal",
  "blaOXA-58-like" = "Carbapenem",
  "abaF"          = "Efflux",
  "gyrA"          = "Fluoroquinolone",
  "parC"          = "Fluoroquinolone",
  "aph(6)"        = "Aminoglycoside",
  "aph(3'')"      = "Aminoglycoside"
)

top20[, res_class := res_class[gene]]
top20[is.na(res_class), res_class := "Other"]

class_pal <- c(
  "Aminoglycoside"  = "#E69F00",
  "Beta-lactam"     = "#D55E00",
  "Carbapenem"      = "#8B0000",
  "Efflux"          = "#0072B2",
  "Fluoroquinolone" = "#009E73",
  "Macrolide"       = "#CC79A7",
  "Metal"           = "#888888",
  "Sulfonamide"     = "#56B4E9",
  "Tetracycline"    = "#F0E442",
  "Disinfectant"    = "#999933",
  "Rifampicin"      = "#AA4499",
  "Colistin"        = "#44AA99",
  "Other"           = "#AAAAAA"
)

# PALETA COMPARTIMENTOS
comp_pal <- c(
  "Human"       = "#2166AC",
  "Animal"      = "#D6604D",
  "Environment" = "#4DAC26"
)

# BANDAS ALTERNADAS
n_genes <- nlevels(top20$gene_f)
band_df <- data.frame(
  ymin = seq(0.5, n_genes - 0.5, by = 2),
  ymax = seq(1.5, n_genes + 0.5, by = 2)
)

# DUMBBELL: segmento entre min(Animal,Env) y Human
top20[, nonhuman_max := pmax(animal_prev, env_prev)]
top20[, nonhuman_min := pmin(animal_prev, env_prev)]

# Shade between Human and nonhuman_max (the gap)
top20[, gap_lo := pmin(human_prev, nonhuman_max)]
top20[, gap_hi := pmax(human_prev, nonhuman_max)]
top20[, human_dominates := human_prev >= nonhuman_max]

# PANEL A: DUMBBELL
p_left <- ggplot(top20, aes(y = gene_f)) +

  # Bandas alternadas
  geom_rect(
    data        = band_df,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
    fill        = "#F4F6F8", color = NA, inherit.aes = FALSE
  ) +

  # Segmento gris completo (rango total min→max)
  geom_segment(
    aes(x    = pmin(human_prev, nonhuman_min),
        xend = pmax(human_prev, nonhuman_max),
        yend = gene_f),
    color = "#DDDDDD", linewidth = 2.0, lineend = "round"
  ) +

  # Franja sombreada entre Human y nonhuman_max (la "brecha")
  geom_segment(
    aes(x    = gap_lo, xend = gap_hi, yend = gene_f,
        color = human_dominates),
    linewidth = 2.0, lineend = "round", alpha = 0.35
  ) +
  scale_color_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
    guide  = "none"
  ) +

  # Línea de referencia 10%
  geom_vline(xintercept = 0.10, linetype = "dashed",
             color = "#BBBBBB", linewidth = 0.45) +
  # Línea de referencia 50%
  geom_vline(xintercept = 0.50, linetype = "dotted",
             color = "#CCCCCC", linewidth = 0.40) +

  # Punto Animal
  geom_point(
    aes(x = animal_prev),
    shape = 22, size = 3.8, fill = "#D6604D",
    color = "white", stroke = 0.7
  ) +
  # Punto Environment
  geom_point(
    aes(x = env_prev),
    shape = 23, size = 3.8, fill = "#4DAC26",
    color = "white", stroke = 0.7
  ) +
  # Punto Human (encima — más prominente)
  geom_point(
    aes(x = human_prev),
    shape = 21, size = 4.5, fill = "#2166AC",
    color = "white", stroke = 0.8
  ) +

  # Etiqueta Human
  geom_text(
    aes(x = human_prev,
        label = paste0(round(human_prev * 100, 0), "%")),
    nudge_y  = 0.40,
    nudge_x  = 0.010,
    size     = 6.5, family = FONT,
    fontface = "bold", color = "#2166AC", hjust = 0
  ) +

  # Anotaciones de referencia
  annotate("text", x = 0.103, y = 0.25, label = "10%",
           size = 6.5, hjust = 0, family = FONT,
           color = "#BBBBBB", fontface = "italic") +
  annotate("text", x = 0.503, y = 0.25, label = "50%",
           size = 6.5, hjust = 0, family = FONT,
           color = "#CCCCCC", fontface = "italic") +

  # Leyenda manual (puntos de los 3 compartimentos)


  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.08),
    breaks = seq(0, 1.0, 0.25),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_discrete(expand = expansion(add = c(1.0, 2.2))) +

  coord_cartesian(clip = "off") +

  labs(
    x = "Mean prevalence",
    y = NULL
  ) +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    legend.position      = "none",
    panel.grid.major.y   = element_blank(),
    panel.grid.major.x   = element_line(color = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor.x   = element_blank(),
    axis.text.y  = element_text(size = 19.0, face = "italic",
                                family = FONT, color = "#111111"),
    axis.text.x  = element_text(size = 18.0, family = FONT, color = "#333333"),
    axis.title.x = element_text(size = 17.5, family = FONT, color = "#555555",
                                margin = margin(t = 8), lineheight = 1.5),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 12, r = 6, b = 12, l = 12)
  )

# PANEL B: CONVERGENCE SCORE por clase de resistencia
# Normalizar convergence_score para mostrar como barra 0-1
top20[, cs_scaled := convergence_score / max(convergence_score)]

p_right <- ggplot(top20, aes(x = cs_scaled, y = gene_f, fill = res_class)) +

  # Bandas alternadas
  geom_rect(
    data        = band_df,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
    fill        = "#F4F6F8", color = NA, inherit.aes = FALSE
  ) +

  # Barras de convergence score
  geom_col(width = 0.55, alpha = 0.88) +

  # Etiqueta del score (valor original)
  geom_text(
    aes(label = sprintf("%.2f", convergence_score)),
    hjust    = -0.12,
    size     = 6.5,
    family   = FONT,
    fontface = "bold",
    color    = "#333333"
  ) +

  # Línea vertical de referencia media
  geom_vline(
    xintercept = mean(top20$cs_scaled),
    linetype   = "dashed",
    color      = "#AAAAAA",
    linewidth  = 0.45
  ) +
  annotate("text",
           x = mean(top20$cs_scaled) + 0.02,
           y = 0.25,
           label = "mean", size = 6.5,
           family = FONT, color = "#AAAAAA", fontface = "italic", hjust = 0) +

  scale_fill_manual(
    values = class_pal,
    name   = "Resistance\nclass"
  ) +
  scale_x_continuous(
    limits = c(0, 1.35),
    breaks = seq(0, 1.0, 0.5),
    labels = c("0", "0.5", "1.0"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = c(1.0, 2.2))) +

  labs(x = "Convergence score", y = NULL) +

  guides(
    fill = guide_legend(
      keywidth  = unit(0.75, "cm"),
      keyheight = unit(0.70, "cm"),
      ncol      = 1
    )
  ) +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 19.0, face = "bold",
                                    family = FONT, color = "#222222",
                                    lineheight = 1.3),
    legend.text      = element_text(size = 17.5, family = FONT),
    legend.key.spacing = unit(0.2, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor.x = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 18.0, family = FONT, color = "#333333"),
    axis.title.x       = element_text(size = 17.5, family = FONT, color = "#555555",
                                      margin = margin(t = 8), lineheight = 1.5),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 12, r = 14, b = 12, l = 4)
  )

# COMPOSICIÓN FINAL
final <- (p_left | p_right) +
  plot_layout(widths = c(2.2, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# EXPORT
png_path  <- file.path(out_dir, "Fig3_one_health_convergence.png")
tiff_path <- file.path(out_dir, "Fig3_one_health_convergence.tiff")
ggsave(png_path,  final, width = 16, height = 10, dpi = 600, bg = "white")
ggsave(tiff_path, final, width = 16, height = 10, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_path,  "\n")
cat("\u2713 TIFF:", tiff_path, "\n")
