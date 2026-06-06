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
sup_dir <- file.path(base, "results/supplementary")
out_dir <- file.path(base, "results/figures/supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# DATA
wl   <- fread(file.path(base, "results/supplementary/tables_submission/TableS7_whitelist_genes_by_country.csv"))
prev <- fread(file.path(base, "results/supplementary/tables_submission/TableS6_geo_ARG_burden_by_country.csv"))

wl <- wl[!country_norm %in% c("not applicable", "missing", "NA")]
wl <- wl[!is.na(country_norm)]
wl <- wl[total_wl > 0]
wl <- merge(wl, prev[, .(country_norm, n_genomes)],
            by = "country_norm", all.x = TRUE)
wl <- wl[!is.na(n_genomes) & n_genomes >= 30]

# GENES
wl_genes <- c("blandm","blaoxa_23like","blaoxa_24like","blaoxa_58like",
              "blaoxa_143like","blaoxa_235like","blaimp","blakpc",
              "blavim","mcr-4.3","mcr-4.7","pmrb","tet(x3)","tet(x5)")

wl_prev <- copy(wl)
setnames(wl_prev, "tet(x3)", "tet_x3")
setnames(wl_prev, "tet(x5)", "tet_x5")
wl_genes <- gsub("tet\\(x3\\)", "tet_x3", wl_genes)
wl_genes <- gsub("tet\\(x5\\)", "tet_x5", wl_genes)

gene_labels <- c(
  "blandm"        = "blaNDM",
  "blaoxa_23like" = "blaOXA-23",
  "blaoxa_24like" = "blaOXA-24",
  "blaoxa_58like" = "blaOXA-58",
  "blaoxa_143like"= "blaOXA-143",
  "blaoxa_235like"= "blaOXA-235",
  "blaimp"        = "blaIMP",
  "blakpc"        = "blaKPC",
  "blavim"        = "blaVIM",
  "mcr-4.3"       = "mcr-4.3",
  "mcr-4.7"       = "mcr-4.7",
  "tet_x3"        = "tet(X3)",
  "tet_x5"        = "tet(X5)",
  "pmrb"          = "pmrB"
)

# Gene class for header colour bands
gene_class <- c(
  "blandm"        = "Carbapenem",
  "blaoxa_23like" = "Carbapenem",
  "blaoxa_24like" = "Carbapenem",
  "blaoxa_58like" = "Carbapenem",
  "blaoxa_143like"= "Carbapenem",
  "blaoxa_235like"= "Carbapenem",
  "blaimp"        = "Carbapenem",
  "blakpc"        = "Carbapenem",
  "blavim"        = "Carbapenem",
  "mcr-4.3"       = "Colistin",
  "mcr-4.7"       = "Colistin",
  "pmrb"          = "Colistin",
  "tet_x3"        = "Tetracycline",
  "tet_x5"        = "Tetracycline"
)

class_pal <- c(
  "Carbapenem"  = "#7B2FBE",
  "Colistin"    = "#882255",
  "Tetracycline"= "#8B6914"
)

# Calculate prevalence
for (g in wl_genes) {
  if (g %in% colnames(wl_prev))
    set(wl_prev, j = g, value = wl_prev[[g]] / wl_prev$n_genomes)
}

# Sort countries by total WHO gene prevalence (descending)
wl_prev[, total_prev := 0]
for (g in wl_genes) {
  if (g %in% colnames(wl_prev)) {
    vals <- wl_prev[[g]]; vals[is.na(vals)] <- 0
    wl_prev[, total_prev := total_prev + vals]
  }
}
wl_prev <- wl_prev[order(-total_prev)]
wl_prev[, country_f := factor(country_norm, levels = rev(country_norm))]

# Long format
wl_long <- rbindlist(lapply(wl_genes, function(g) {
  if (!g %in% colnames(wl_prev)) return(NULL)
  data.table(country_f = wl_prev$country_f,
             gene = g, prev = wl_prev[[g]])
}))
wl_long[is.na(prev) | prev == 0, prev := NA]

genes_present <- wl_genes[wl_genes %in% unique(wl_long$gene)]
wl_long[, gene_f := factor(gene,
                            levels = genes_present,
                            labels = gene_labels[genes_present])]

# n_genomes for right-side bar
n_df <- wl_prev[, .(country_f, n_genomes, total_prev)]

# Gene detection frequency for top bar
gene_freq <- wl_long[!is.na(prev), .N, by = gene_f]
gene_freq[, pct := N / nrow(wl_prev)]

# Number of countries
n_countries <- nlevels(wl_long$country_f)
n_genes     <- length(genes_present)

# FIGURE
p <- ggplot(wl_long, aes(x = gene_f, y = country_f, fill = prev)) +

  # Heatmap tiles
  geom_tile(color = "white", linewidth = 0.30) +

  # Prevalence labels — only ≥5% to reduce crowding in dense rows
  geom_text(
    aes(label = ifelse(!is.na(prev) & prev >= 0.05,
                       paste0(round(prev * 100), "%"), "")),
    size     = 6.2, family = FONT,
    fontface = "bold", color = "#222222",
    show.legend = FALSE
  ) +

  # Colour scale
  scale_fill_gradientn(
    colors   = c("#FFFEF0","#FEC44F","#FE9929","#EC7014","#CC4C02","#8C2D04"),
    values   = rescale(c(0, 0.10, 0.30, 0.55, 0.75, 1.0)),
    limits   = c(0, 1),
    na.value = "#F0F0F0",
    name     = "Prevalence",
    labels   = percent_format(accuracy = 1),
    breaks   = c(0, 0.25, 0.50, 0.75, 1.0),
    guide    = guide_colorbar(
      barheight      = unit(5.5, "cm"),
      barwidth       = unit(0.50, "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "grey50",
      frame.colour   = "grey60",
      title.theme    = element_text(size = 19.0, face = "bold",
                                    family = FONT, color = "#222222"),
      label.theme    = element_text(size = 17.5, family = FONT,
                                    color = "#444444")
    )
  ) +

  # Gene class colour band at top (as annotation above plot)
  annotate("rect",
           xmin = seq_along(genes_present) - 0.46,
           xmax = seq_along(genes_present) + 0.46,
           ymin = n_countries + 0.55,
           ymax = n_countries + 0.95,
           fill = class_pal[gene_class[genes_present]],
           color = NA) +

  # n_genomes bar at right side (as annotation)
  annotate("rect",
           xmin = n_genes + 1.20,
           xmax = n_genes + 1.20 + (n_df$n_genomes / max(n_df$n_genomes)) * 1.2,
           ymin = as.numeric(n_df$country_f) - 0.40,
           ymax = as.numeric(n_df$country_f) + 0.40,
           fill = "#AAAAAA", color = NA, alpha = 0.7) +

  # n_genomes labels
  annotate("text",
           x     = n_genes + 1.20 + (n_df$n_genomes / max(n_df$n_genomes)) * 1.2 + 0.08,
           y     = as.numeric(n_df$country_f),
           label = format(n_df$n_genomes, big.mark = ","),
           hjust = 0, size = 6.5, family = FONT,
           color = "#666666", fontface = "italic") +

  # "Genomes" header above bar
  annotate("text",
           x = n_genes + 2.20, y = n_countries + 1.35,
           label = "Genomes", hjust = 0.5,
           size = 7.0, family = FONT,
           fontface = "bold", color = "#444444") +

  # Gene class labels — all at same height, no stagger needed
  annotate("text",
           x     = c(mean(which(gene_class[genes_present] == "Carbapenem")),
                     mean(which(gene_class[genes_present] == "Colistin")),
                     mean(which(gene_class[genes_present] == "Tetracycline"))),
           y     = n_countries + 1.35,
           label = c("Carbapenem", "Colistin", "Tetracycline"),
           hjust = 0.5, size = 7.0, family = FONT,
           fontface = "bold",
           color = class_pal[c("Carbapenem","Colistin","Tetracycline")]) +



  scale_x_discrete(expand = expansion(add = c(0.5, 2.8))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.8))) +

  coord_cartesian(clip = "off") +

  labs(x = NULL, y = NULL) +

  theme_minimal(base_size = 22, base_family = FONT) +
  theme(
    axis.text.x      = element_text(size = 19.0, angle = 45, hjust = 1,
                                    face = "italic", family = FONT,
                                    color = "#222222"),
    axis.text.y      = element_text(size = 19.0, family = FONT,
                                    color = "#222222"),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.margin    = margin(l = 90),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 28, r = 60, b = 28, l = 12)
  )

# EXPORT
png_s3  <- file.path(out_dir, "FigS2_geo_whitelist_heatmap.png")
tiff_s3 <- file.path(out_dir, "FigS2_geo_whitelist_heatmap.tiff")
ggsave(png_s3,  p, width = 14, height = 15, dpi = 600, bg = "white")
ggsave(tiff_s3, p, width = 14, height = 15, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("\u2713 PNG :", png_s3,  "\n")
cat("\u2713 TIFF:", tiff_s3, "\n")
