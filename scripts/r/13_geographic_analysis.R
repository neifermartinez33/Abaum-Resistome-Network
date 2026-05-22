library(data.table)
library(igraph)

cat("=== Phase 13: Geographic Analysis ===\n")
cat(format(Sys.time()), "\n\n")

# ── LOAD DATA ────────────────────────────────────────────────
meta <- fread("~/genomes_pass_qc_validated.csv")
mat  <- fread("~/matrix_output_v2/matrix1_gene_families.csv")
setnames(meta, "Name", "genome_id")
setnames(mat,  colnames(mat)[1], "genome_id")
gene_cols <- colnames(mat)[-1]

# ── 1. NORMALIZE COUNTRY NAMES ───────────────────────────────
cat("=== 1. Country Normalization ===\n")

normalize_country <- function(x) {
  x <- trimws(x)
  # Extract base country before ":"
  x <- sub(":.*$", "", x)
  x <- trimws(x)
  # Standardize variants
  x <- gsub("^UK$", "United Kingdom", x)
  x <- gsub("^USA$", "United States", x)
  x <- gsub("^South Korea.*", "South Korea", x)
  x <- gsub("^China.*", "China", x)
  x <- gsub("^France.*", "France", x)
  x <- gsub("^Germany.*", "Germany", x)
  x <- gsub("^Thailand.*", "Thailand", x)
  x <- gsub("^Poland.*", "Poland", x)
  x <- gsub("^Serbia.*", "Serbia", x)
  x <- gsub("^South Africa.*", "South Africa", x)
  x <- gsub("^Unknown$|^unknown$|^NA$", NA_character_, x)
  return(x)
}

meta[, country_norm := normalize_country(country)]
country_n <- meta[!is.na(country_norm), .N, by=country_norm][order(-N)]
cat("Top 30 countries after normalization:\n")
print(country_n[1:30])

# ── MERGE ────────────────────────────────────────────────────
merged <- merge(mat, meta[, .(genome_id, country_norm)], by="genome_id")
merged_geo <- merged[!is.na(country_norm)]
cat("\nGenomes with valid country:", nrow(merged_geo), "\n")

# Focus on countries with >= 50 genomes
country_keep <- country_n[N >= 50, country_norm]
merged_top <- merged_geo[country_norm %in% country_keep]
cat("Genomes in countries >= 50:", nrow(merged_top), "\n")
cat("Countries analyzed:", length(country_keep), "\n\n")

# ── 2. ARG PREVALENCE BY COUNTRY ─────────────────────────────
cat("=== 2. ARG Prevalence by Country ===\n")

prev_country <- merged_top[, lapply(.SD, mean), 
                            by=country_norm, .SDcols=gene_cols]
prev_country[, n_genomes := merged_top[, .N, by=country_norm]$N]

# Mean ARGs per genome per country
prev_country[, mean_args := rowSums(.SD), .SDcols=gene_cols]
prev_country <- prev_country[order(-mean_args)]

cat("Top 15 countries by mean ARG count per genome:\n")
print(prev_country[1:15, .(country_norm, n_genomes, mean_args=round(mean_args,3))])

# ── 3. WHITELIST GENES BY COUNTRY ────────────────────────────
cat("\n=== 3. WHO Whitelist Genes by Country ===\n")

whitelist <- c("blandm","blaoxa_23like","blaoxa_24like","blaoxa_58like",
               "blaoxa_143like","blaoxa_235like","blaimp","blakpc","blavim",
               "mcr-4.3","mcr-4.7","tet(x3)","tet(x5)","pmrb")
wl_avail <- whitelist[whitelist %in% gene_cols]

wl_country <- merged_geo[, lapply(.SD, function(x) sum(x>0)),
                          by=country_norm, .SDcols=wl_avail]
wl_country[, total_wl := rowSums(.SD), .SDcols=wl_avail]
wl_country <- wl_country[total_wl > 0][order(-total_wl)]

cat("Countries with >= 1 WHO whitelist gene:\n")
print(wl_country[, .(country_norm, total_wl, 
                     blandm, blaoxa_23like, blaimp, blakpc, blavim,
                     `mcr-4.3`, `mcr-4.7`, `tet(x3)`, `tet(x5)`)])

# ── 4. GEOGRAPHIC NETWORK ANALYSIS ───────────────────────────
cat("\n=== 4. Network Metrics by Country ===\n")

g_ref <- readRDS("~/matrix_output_v2/network_topology.rds")
ref_edges <- as_data_frame(g_ref, what="edges")[, c("from","to")]

geo_topology <- list()
for (ctry in country_keep) {
  ctry_mat <- merged_top[country_norm == ctry, ..gene_cols]
  if (nrow(ctry_mat) < 30) next
  M <- as.matrix(ctry_mat)
  
  ew <- numeric(nrow(ref_edges))
  for (k in seq_len(nrow(ref_edges))) {
    a <- M[, ref_edges$from[k]]
    b <- M[, ref_edges$to[k]]
    inter <- sum(a==1 & b==1)
    union <- sum(a==1 | b==1)
    ew[k] <- if (union==0) 0 else inter/union
  }
  
  ke <- ref_edges[ew > 0.30, ]
  g_c <- graph_from_data_frame(ke, directed=FALSE, vertices=gene_cols)
  
  geo_topology[[ctry]] <- data.table(
    country     = ctry,
    n_genomes   = nrow(ctry_mat),
    n_edges     = ecount(g_c),
    density     = round(edge_density(g_c), 4),
    clustering  = round(transitivity(g_c, type="global"), 4),
    mean_args   = round(mean(rowSums(M)), 3)
  )
}

geo_dt <- rbindlist(geo_topology)[order(-mean_args)]
cat("\nTop 20 countries by mean ARGs/genome:\n")
print(geo_dt[1:20])

# ── 5. REGIONAL GROUPING ─────────────────────────────────────
cat("\n=== 5. ARG Burden by World Region ===\n")

region_map <- list(
  "East Asia"      = c("China","Japan","South Korea","Taiwan"),
  "Southeast Asia" = c("Thailand","Vietnam","Cambodia","Singapore","Indonesia"),
  "South Asia"     = c("India","Pakistan","Nepal","Bangladesh","Sri Lanka"),
  "Europe"         = c("Germany","France","United Kingdom","Belgium","Greece",
                       "Italy","Spain","Romania","Hungary","Poland","Serbia","Russia"),
  "North America"  = c("United States","Canada","Mexico"),
  "South America"  = c("Brazil","Colombia","Argentina","Peru"),
  "Middle East"    = c("Saudi Arabia","Israel","Iran","Turkey"),
  "Oceania"        = c("Australia","New Zealand"),
  "Africa"         = c("South Africa","Kenya","Tunisia","Nigeria","Ghana","Egypt"),
  "Central Asia"   = c("Afghanistan")
)

region_dt <- rbindlist(lapply(names(region_map), function(reg) {
  countries_in_reg <- region_map[[reg]]
  reg_data <- merged_geo[country_norm %in% countries_in_reg]
  if (nrow(reg_data) == 0) return(NULL)
  prev <- colMeans(reg_data[, ..gene_cols])
  data.table(
    region     = reg,
    n_genomes  = nrow(reg_data),
    n_countries = length(intersect(countries_in_reg,
                                   unique(merged_geo$country_norm))),
    mean_args  = round(sum(prev), 3),
    blaoxa_23like = round(prev["blaoxa_23like"], 3),
    blandm     = round(prev["blandm"], 3),
    gyra       = round(prev["gyra"], 3)
  )
}))[order(-mean_args)]

cat("\nARG burden by world region:\n")
print(region_dt)

# ── 6. SAVE OUTPUTS ─────────────────────────────────────────
fwrite(prev_country, "~/matrix_output_v2/geo_prevalence_by_country.csv")
fwrite(wl_country,   "~/matrix_output_v2/geo_whitelist_by_country.csv")
fwrite(geo_dt,       "~/matrix_output_v2/geo_network_metrics.csv")
fwrite(region_dt,    "~/matrix_output_v2/geo_regional_burden.csv")

cat("\nFiles saved:\n")
cat("  geo_prevalence_by_country.csv\n")
cat("  geo_whitelist_by_country.csv\n")
cat("  geo_network_metrics.csv\n")
cat("  geo_regional_burden.csv\n")
cat("\nDone:", format(Sys.time()), "\n")
