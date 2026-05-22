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
wl   <- fread(file.path(sup_dir, "geo_whitelist_by_country.csv"))
prev <- fread(file.path(sup_dir, "geo_prevalence_by_country.csv"))

# Filter invalid entries
wl <- wl[!country_norm %in% c("not applicable", "missing", "NA")]
wl <- wl[total_wl > 0]

# Add n_genomes and filter countries with >= 30 genomes
wl <- merge(wl, prev[, .(country_norm, n_genomes)],
            by = "country_norm", all.x = TRUE)
wl <- wl[!is.na(n_genomes) & n_genomes >= 30]

# ── GENES IN EXACT FILE ORDER ─────────────────────────────────
# Manually verified: this is the order in geo_whitelist_by_country.csv
wl_genes <- c("blandm","blaoxa_23like","blaoxa_24like","blaoxa_58like",
              "blaoxa_143like","blaoxa_235like","blaimp","blakpc",
              "blavim","mcr-4.3","mcr-4.7","tet(x3)","tet(x5)","pmrb")

# Calculate prevalence by dividing by n_genomes
# Use explicit copy to avoid altering the original data.table
wl_prev <- copy(wl)

# ── RENAME columns with parentheses — avoids unstable behavior ──
# data.table::set() and [[ indexing have bugs with names like tet(x3)
setnames(wl_prev, "tet(x3)", "tet_x3")
setnames(wl_prev, "tet(x5)", "tet_x5")

# Update gene vector with safe names
wl_genes <- gsub("tet\\(x3\\)", "tet_x3", wl_genes)
wl_genes <- gsub("tet\\(x5\\)", "tet_x5", wl_genes)

# Labels updated with safe names as keys
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
for (g in wl_genes) {
  if (g %in% colnames(wl_prev)) {
    set(wl_prev, j = g, value = wl_prev[[g]] / wl_prev$n_genomes)
  }
}

# Sort countries by sum of whitelist gene prevalences
# total_prev — explicit loop to avoid reordering with sapply
# sapply over names with parentheses like tet(x3) generates unstable matrices
wl_prev[, total_prev := 0]
for (g in wl_genes) {
  if (g %in% colnames(wl_prev)) {
    vals <- wl_prev[[g]]
    vals[is.na(vals)] <- 0
    wl_prev[, total_prev := total_prev + vals]
  }
}
wl_prev <- wl_prev[order(-total_prev)]
wl_prev[, country_f := factor(country_norm, levels = rev(country_norm))]

# ── LONG FORMAT — WITHOUT melt(), access by name ──────────────
# Manual gene-by-gene construction to avoid data.table reordering
wl_long <- rbindlist(lapply(wl_genes, function(g) {
  if (!g %in% colnames(wl_prev)) return(NULL)
  data.table(
    country_f = wl_prev$country_f,
    gene      = g,
    prev      = wl_prev[[g]]   # name access — immune to reordering
  )
}))
wl_long[is.na(prev) | prev == 0, prev := NA]

# gene_labels already defined above with tet_x3/tet_x5 — do not redefine here

genes_present <- wl_genes[wl_genes %in% unique(wl_long$gene)]
wl_long[, gene := factor(gene,
                          levels = genes_present,
                          labels = gene_labels[genes_present])]

# ── CONSOLE VERIFICATION ──────────────────────────────────────
cat("\n=== Taiwan tet(X3) and tet(X5) Verification ===\n")
tw <- wl_long[as.character(country_f) == "Taiwan" & gene %in% c("tet(X3)","tet(X5)"),
              .(gene = as.character(gene), prev)]
print(tw)
vim_check <- wl_long[gene == "blaVIM" & !is.na(prev),
                     .(country = as.character(country_f), prev)]
print(vim_check)

cat("\n=== mcr-4.3 Verification (expected: China, Brazil, South Korea, Czech Rep) ===\n")
mcr_check <- wl_long[gene == "mcr-4.3" & !is.na(prev),
                     .(country = as.character(country_f), prev)]
print(mcr_check)

# ── FIGURE ───────────────────────────────────────────────────
p_s3 <- ggplot(wl_long, aes(x = gene, y = country_f, fill = prev)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(!is.na(prev) & prev >= 0.05,
                       paste0(round(prev * 100, 0), "%"), "")),
    size = 2.85, color = "#222222", show.legend = FALSE
  ) +
  scale_fill_gradientn(
    colors   = c("#FFFEF0","#FEC44F","#FE9929","#EC7014","#CC4C02","#8C2D04"),
    values   = rescale(c(0, 0.10, 0.30, 0.55, 0.75, 1.0)),
    limits   = c(0, 1),
    na.value = "#F0F0F0",
    name     = "Prevalence",
    labels   = percent_format(accuracy = 1),
    breaks   = c(0, 0.25, 0.50, 0.75, 1.0),
    guide    = guide_colorbar(
      barheight      = unit(5.0, "cm"),
      barwidth       = unit(0.42, "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "grey50",
      frame.colour   = "grey60"
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  annotate("point", x = 0.55, y = -0.9,
           shape = 22, size = 3.5, fill = "#F0F0F0", color = "grey60") +
  annotate("text", x = 0.85, y = -0.9,
           label = "Not detected (or < 30 genomes)",
           hjust = 0, size = 2.85, color = "#666666") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10.5) +
  theme(
    axis.text.x      = element_text(size = 9.0, angle = 45, hjust = 1,
                                    face = "italic", color = "#222222"),
    axis.text.y      = element_text(size = 8.5, color = "#222222"),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 9.5, face = "bold"),
    legend.text      = element_text(size = 8.5),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(12, 12, 24, 12)
  )

# ── EXPORT ───────────────────────────────────────────────────
png_s3  <- file.path(out_dir, "FigS3_geo_whitelist_heatmap.png")
tiff_s3 <- file.path(out_dir, "FigS3_geo_whitelist_heatmap.tiff")
ggsave(png_s3,  p_s3, width = 13, height = 14, dpi = 300, bg = "white")
ggsave(tiff_s3, p_s3, width = 13, height = 14, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("OK: FigS3 PNG :", png_s3,  "\n")
cat("OK: FigS3 TIFF:", tiff_s3, "\n")
