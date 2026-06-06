suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

geo <- fread(file.path(tables, "geo_regional_burden.csv"))
geo <- geo[order(-mean_args)]
geo[, region_f := factor(region, levels = rev(region))]
geo[, rank_n   := seq_len(.N)]

mean_global <- mean(geo$mean_args)
n_reg       <- nrow(geo)

# NORMALISE each variable independently 0-1
geo[, args_norm := (mean_args    - min(mean_args))    / (max(mean_args)    - min(mean_args))]
geo[, gyra_norm := (gyra         - min(gyra))         / (max(gyra)         - min(gyra))]
geo[, oxa_norm  := (blaoxa_23like- min(blaoxa_23like))/ (max(blaoxa_23like)- min(blaoxa_23like))]
geo[, ndm_norm  := (blandm       - min(blandm))       / (max(blandm)       - min(blandm))]

# LONG FORMAT — one row per cell
cells <- rbind(
  geo[, .(region_f, col_id = 1L, col_name = "Mean ARG families",
          fill_val = args_norm,
          raw_val  = mean_args,
          lab      = sprintf("%.2f", mean_args),
          dark     = args_norm  > 0.52)],
  geo[, .(region_f, col_id = 2L, col_name = "gyrA",
          fill_val = gyra_norm,
          raw_val  = gyra,
          lab      = paste0(round(gyra * 100), "%"),
          dark     = gyra_norm  > 0.52)],
  geo[, .(region_f, col_id = 3L, col_name = "blaOXA-23-like",
          fill_val = oxa_norm,
          raw_val  = blaoxa_23like,
          lab      = paste0(round(blaoxa_23like * 100), "%"),
          dark     = oxa_norm   > 0.52)],
  geo[, .(region_f, col_id = 4L, col_name = "blaNDM",
          fill_val = ndm_norm,
          raw_val  = blandm,
          lab      = paste0(round(blandm * 100), "%"),
          dark     = ndm_norm   > 0.52)]
)
cells[, col_id := as.numeric(col_id)]

# COLOUR PALETTES per column
cells[, fill_enc := fcase(
  col_id == 1L, fill_val * 0.235,
  col_id == 2L, 0.255 + fill_val * 0.235,
  col_id == 3L, 0.510 + fill_val * 0.235,
  col_id == 4L, 0.765 + fill_val * 0.235
)]

pal_master <- c(
  colorRampPalette(c("#FFF5F0","#FC8D59","#990000"))(60),
  colorRampPalette(c("#FFFDE7","#FDB863","#B35806"))(60),
  colorRampPalette(c("#F7FCF5","#74C476","#005A32"))(60),
  colorRampPalette(c("#FCFBFD","#B87FD4","#3F007D"))(60)
)

# N LABEL
geo[, n_lab := paste0(format(n_genomes, big.mark = ",", trim = TRUE),
                      " (", n_countries,
                      ifelse(n_countries == 1, " country)", " countries)"))]

# FIGURE
p <- ggplot(cells, aes(x = col_id, y = region_f)) +

  # Tiles
  geom_tile(
    aes(fill = fill_enc),
    color = "white", linewidth = 1.0,
    width = 0.96, height = 0.92
  ) +
  scale_fill_gradientn(
    colours = pal_master,
    limits  = c(0, 1),
    guide   = "none"
  ) +

  # Cell text — white on dark, dark on light
  geom_text(
    aes(label = lab,
        color = dark),
    size     = 6.2, fontface = "bold", family = FONT
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "#2C2C2C"),
    guide  = "none"
  ) +

  # Vertical separator between burden and markers
  # Column header background rect
  annotate("rect",
           xmin = 0.52, xmax = 4.48,
           ymin = n_reg + 0.58, ymax = n_reg + 1.40,
           fill = "#F2F2F2", color = NA) +
  annotate("segment",
           x = 0.52, xend = 4.48,
           y = n_reg + 0.58, yend = n_reg + 0.58,
           color = "#CCCCCC", linewidth = 0.6) +

  # Column headers
  annotate("text", x = 1, y = n_reg + 1.00,
           label = "Mean ARG\nfamilies",
           size = 6.0, hjust = 0.5, family = FONT,
           fontface = "bold", color = "#333333", lineheight = 1.1) +
  annotate("text", x = 2, y = n_reg + 1.00,
           label = "gyrA",
           size = 6.0, hjust = 0.5, family = FONT,
           fontface = "bold.italic", color = "#B35806") +
  annotate("text", x = 3, y = n_reg + 1.00,
           label = "blaOXA-23-like",
           size = 6.0, hjust = 0.5, family = FONT,
           fontface = "bold.italic", color = "#005A32") +
  annotate("text", x = 4, y = n_reg + 1.00,
           label = "blaNDM",
           size = 6.0, hjust = 0.5, family = FONT,
           fontface = "bold.italic", color = "#3F007D") +



  # "Genomes" column header
  annotate("text",
           x = 4.65, y = n_reg + 1.00,
           label = "Genomes",
           size = 6.0, hjust = 0.0, family = FONT,
           fontface = "bold", color = "#333333") +

  # n= labels (right)
  annotate("text",
           x     = rep(4.65, n_reg),
           y     = geo$region_f,
           label = geo$n_lab,
           size  = 5.8, hjust = 0.0, family = FONT,
           color = "#666666", fontface = "italic") +

  # Global mean label below column 1
  annotate("text",
           x = 1, y = 0.30,
           label = sprintf("Global mean: %.2f", mean_global),
           size = 5.8, hjust = 0.5, family = FONT,
           color = "#990000", fontface = "italic") +

  # LEYENDAS DEGRADADAS INFERIORES (TEXTO DESPLAZADO HACIA ARRIBA)
  
  # Leyenda: gyrA
  annotate("rect",
           xmin  = seq(1.62, 2.42, length.out = 31)[-(31)],
           xmax  = seq(1.62, 2.42, length.out = 31)[-1],
           ymin  = -0.70, ymax = -0.46,
           fill  = colorRampPalette(c("#FFFDE7","#FDB863","#B35806"))(30),
           color = NA) +
  annotate("rect", xmin=1.62, xmax=2.42, ymin=-0.70, ymax=-0.46,
           fill=NA, color="#999999", linewidth=0.3) +
  annotate("text", x=2.02, y=-0.30,  # Modificado de -0.40 a -0.30
           label=paste0(round(min(geo$gyra)*100), "%  –  ", round(max(geo$gyra)*100), "%"),
           size=6.0, hjust=0.5, family=FONT, color="#555555") +

  # Leyenda: blaOXA-23-like
  annotate("rect",
           xmin  = seq(2.62, 3.42, length.out = 31)[-(31)],
           xmax  = seq(2.62, 3.42, length.out = 31)[-1],
           ymin  = -0.70, ymax = -0.46,
           fill  = colorRampPalette(c("#F7FCF5","#74C476","#005A32"))(30),
           color = NA) +
  annotate("rect", xmin=2.62, xmax=3.42, ymin=-0.70, ymax=-0.46,
           fill=NA, color="#999999", linewidth=0.3) +
  annotate("text", x=3.02, y=-0.30,  # Modificado de -0.40 a -0.30
           label=paste0(round(min(geo$blaoxa_23like)*100), "%  –  ", round(max(geo$blaoxa_23like)*100), "%"),
           size=6.0, hjust=0.5, family=FONT, color="#555555") +

  # Leyenda: blaNDM
  annotate("rect",
           xmin  = seq(3.62, 4.42, length.out = 31)[-(31)],
           xmax  = seq(3.62, 4.42, length.out = 31)[-1],
           ymin  = -0.70, ymax = -0.46,
           fill  = colorRampPalette(c("#FCFBFD","#B87FD4","#3F007D"))(30),
           color = NA) +
  annotate("rect", xmin=3.62, xmax=4.42, ymin=-0.70, ymax=-0.46,
           fill=NA, color="#999999", linewidth=0.3) +
  annotate("text", x=4.02, y=-0.30,  # Modificado de -0.40 a -0.30
           label=paste0(round(min(geo$blandm)*100), "%  –  ", round(max(geo$blandm)*100), "%"),
           size=6.0, hjust=0.5, family=FONT, color="#555555") +

  # Configuración de escalas y temas
  scale_x_continuous(
    limits = c(-0.2, 7.5),
    expand = expansion(mult = c(0, 0)),
    breaks = NULL
  ) +
  scale_y_discrete(expand = expansion(add = c(1.4, 1.8))) +

  labs(x = NULL, y = NULL,
       caption = "") +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 19.0, face = "bold",
                                     family = FONT, color = "#111111",
                                     hjust = 1),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.caption      = element_text(size = 14.0, family = FONT, color = "#888888",
                                     hjust = 0, margin = margin(t = 8)),
    plot.margin       = margin(t = 20, r = 20, b = 12, l = 10)
  )

# EXPORT
png_path  <- file.path(out_dir, "Fig4_geo_regional_burden.png")
tiff_path <- file.path(out_dir, "Fig4_geo_regional_burden.tiff")
ggsave(png_path,  p, width = 14, height = 8.5, dpi = 600, bg = "white")
ggsave(tiff_path, p, width = 14, height = 8.5, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_path,  "\n")
cat("\u2713 TIFF:", tiff_path, "\n")
