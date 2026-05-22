library(data.table)
library(igraph)
library(vegan)

cat("=== Phase 11: One Health Analysis ===\n")
cat(format(Sys.time()), "\n\n")

# ── LOAD DATA ────────────────────────────────────────────────
meta <- fread("~/genomes_pass_qc_validated.csv")
mat  <- fread("~/matrix_output_v2/matrix1_gene_families.csv")
g_ref <- readRDS("~/matrix_output_v2/network_topology.rds")

setnames(meta, "Name", "genome_id")
setnames(mat,  colnames(mat)[1], "genome_id")

gene_cols <- colnames(mat)[-1]

# ── 1. DEFINE ONE HEALTH COMPARTMENTS ────────────────────────
cat("=== 1. One Health Compartment Assignment ===\n")

human_terms <- c("hospital", "sputum", "blood", "urine", "wound",
                 "respiratory", "stool", "rectal swab", "skin", "lung",
                 "icu", "tracheal secretion", "balf", "bal", "csf",
                 "oral swab", "abdominal fluid", "pus", "clinical isolate",
                 "wound/abscess", "upper respiratory", "lower respiratory",
                 "sputamentum", "sputums", "blood culture",
                 "cerebrospinal fluid", "cerebrospinal", "balf")

animal_terms <- c("animal", "cow", "bovine", "swine", "pig", "poultry",
                  "chicken", "horse", "dog", "cat", "sheep", "goat",
                  "fish", "bird", "veterinary", "livestock")

env_terms    <- c("soil", "environment", "sewage", "water", "river",
                  "lake", "ocean", "sediment", "effluent", "wastewater",
                  "compost", "manure", "air")

assign_compartment <- function(src) {
  s <- tolower(trimws(src))
  if (s == "" || s %in% c("not applicable", "not available", "missing",
                           "unknown", "not collected", "other", "na")) return("unknown")
  if (any(sapply(human_terms, function(t) grepl(t, s, fixed=TRUE)))) return("human")
  if (any(sapply(animal_terms, function(t) grepl(t, s, fixed=TRUE)))) return("animal")
  if (any(sapply(env_terms,    function(t) grepl(t, s, fixed=TRUE)))) return("environment")
  return("unknown")
}

meta[, compartment := sapply(isolation_source, assign_compartment)]

cat("Compartment distribution:\n")
print(table(meta$compartment))

# ── MERGE ────────────────────────────────────────────────────
merged <- merge(mat, meta[, .(genome_id, compartment, country, isolation_source)],
                by = "genome_id")
cat("\nGenomes after merge:", nrow(merged), "\n")

# Focus on known compartments (exclude unknown)
merged_known <- merged[compartment != "unknown"]
cat("Genomes with known compartment:", nrow(merged_known), "\n")
print(table(merged_known$compartment))

# ── 2. ARG PREVALENCE BY COMPARTMENT ─────────────────────────
cat("\n=== 2. ARG Prevalence by Compartment ===\n")

prev_by_comp <- merged_known[, lapply(.SD, mean), 
                              by = compartment, .SDcols = gene_cols]
prev_by_comp_long <- melt(prev_by_comp, id.vars = "compartment",
                           variable.name = "gene", value.name = "prevalence")

# Top differences between human and animal/environment
human_prev <- prev_by_comp[compartment == "human", ..gene_cols]
animal_prev <- prev_by_comp[compartment == "animal", ..gene_cols]
env_prev    <- prev_by_comp[compartment == "environment", ..gene_cols]

if (nrow(human_prev) > 0 && nrow(animal_prev) > 0) {
  diff_ha <- as.numeric(human_prev) - as.numeric(animal_prev)
  names(diff_ha) <- gene_cols
  cat("\nTop genes enriched in human vs animal (difference in prevalence):\n")
  print(sort(diff_ha, decreasing=TRUE)[1:10])
  cat("\nTop genes enriched in animal vs human:\n")
  print(sort(diff_ha)[1:10])
}

if (nrow(human_prev) > 0 && nrow(env_prev) > 0) {
  diff_he <- as.numeric(human_prev) - as.numeric(env_prev)
  names(diff_he) <- gene_cols
  cat("\nTop genes enriched in human vs environment:\n")
  print(sort(diff_he, decreasing=TRUE)[1:10])
}

# ── 3. STATISTICAL TESTS (Kruskal-Wallis per gene) ───────────
cat("\n=== 3. Statistical Differences Across Compartments ===\n")

kw_results <- list()
for (gene in gene_cols) {
  groups <- merged_known[, .(val = get(gene), compartment)]
  if (length(unique(groups$val)) > 1) {
    kt <- kruskal.test(val ~ compartment, data = groups)
    kw_results[[gene]] <- data.table(
      gene    = gene,
      H_stat  = round(kt$statistic, 4),
      pvalue  = round(kt$p.value, 6),
      sig     = ifelse(kt$p.value < 0.001, "***",
                ifelse(kt$p.value < 0.01,  "**",
                ifelse(kt$p.value < 0.05,  "*", "ns")))
    )
  }
}

kw_dt <- rbindlist(kw_results)
kw_dt <- kw_dt[order(pvalue)]
cat("Kruskal-Wallis results (top 20):\n")
print(kw_dt[1:20])
cat("\nSignificant (p<0.05):", sum(kw_dt$pvalue < 0.05), "/", nrow(kw_dt), "\n")

# ── 4. COMPARTMENT-SPECIFIC SUBNETWORKS ──────────────────────
cat("\n=== 4. Compartment-Specific Subnetworks ===\n")

ref_edges <- as_data_frame(g_ref, what = "edges")[, c("from", "to")]

comp_topology <- list()
for (comp in c("human", "animal", "environment")) {
  comp_mat <- merged_known[compartment == comp, ..gene_cols]
  if (nrow(comp_mat) < 30) next
  M_comp <- as.matrix(comp_mat)

  edge_weights <- numeric(nrow(ref_edges))
  for (k in seq_len(nrow(ref_edges))) {
    a <- M_comp[, ref_edges$from[k]]
    b <- M_comp[, ref_edges$to[k]]
    inter <- sum(a == 1 & b == 1)
    union <- sum(a == 1 | b == 1)
    edge_weights[k] <- if (union == 0) 0 else inter / union
  }

  keep_edges <- ref_edges[edge_weights > 0.30, ]
  g_comp <- graph_from_data_frame(keep_edges, directed = FALSE,
                                   vertices = gene_cols)

  comp_topology[[comp]] <- data.table(
    compartment = comp,
    n_genomes   = nrow(comp_mat),
    n_edges     = ecount(g_comp),
    density     = round(edge_density(g_comp), 4),
    clustering  = round(transitivity(g_comp, type = "global"), 4),
    components  = components(g_comp)$no
  )

  cat(sprintf("\n%s subnetwork (%d genomes):\n", comp, nrow(comp_mat)))
  cat("  Edges:", ecount(g_comp), "| Density:", round(edge_density(g_comp), 4),
      "| Clustering:", round(transitivity(g_comp, type="global"), 4), "\n")

  # Unique edges vs consensus
  comp_edges <- as_data_frame(g_comp, what="edges")
  comp_pairs <- paste(pmin(comp_edges$from, comp_edges$to),
                      pmax(comp_edges$from, comp_edges$to))
  ref_pairs  <- paste(pmin(ref_edges$from, ref_edges$to),
                      pmax(ref_edges$from, ref_edges$to))
  cat("  Edges shared with consensus:", sum(comp_pairs %in% ref_pairs), "\n")
  saveRDS(g_comp, sprintf("~/matrix_output_v2/network_%s.rds", comp))
}

comp_topo_dt <- rbindlist(comp_topology)
cat("\nCompartment topology summary:\n")
print(comp_topo_dt)

# ── 5. BETA DIVERSITY BETWEEN COMPARTMENTS ───────────────────
cat("\n=== 5. Resistome Beta-Diversity (Bray-Curtis) ===\n")

# Sample up to 500 per compartment for efficiency
set.seed(42)
sampled <- merged_known[compartment != "unknown",
           .SD[sample(.N, min(.N, 500))], by = compartment]

M_beta <- as.matrix(sampled[, ..gene_cols])
groups <- sampled$compartment

bc_dist <- vegdist(M_beta, method = "bray")
cat("Bray-Curtis distance matrix computed\n")

# PERMANOVA
perm <- adonis2(bc_dist ~ groups, permutations = 999)
cat("\nPERMANOVA results:\n")
print(perm)

# Betadisper (homogeneity of dispersion)
bd <- betadisper(bc_dist, groups)
cat("\nBeta-dispersion by compartment:\n")
print(bd$group.distances)
bd_test <- permutest(bd, permutations = 999)
cat("Betadisper permutation test p-value:", round(bd_test$tab$`Pr(>F)`[1], 4), "\n")

# ── 6. SHARED ARG PROFILES ───────────────────────────────────
cat("\n=== 6. Shared ARG Profiles Across Compartments ===\n")

# Genes present in ALL 3 compartments (>10% prevalence in each)
comp_prev <- merged_known[compartment %in% c("human","animal","environment"),
             lapply(.SD, mean), by = compartment, .SDcols = gene_cols]

if (nrow(comp_prev) == 3) {
  shared_genes <- gene_cols[sapply(gene_cols, function(g) {
    all(comp_prev[[g]] > 0.10)
  })]
  cat("Genes with >10% prevalence in ALL 3 compartments (", length(shared_genes), "):\n")
  cat(paste(shared_genes, collapse=", "), "\n")

  # One Health convergence score per gene
  conv_score <- sapply(gene_cols, function(g) {
    vals <- as.numeric(comp_prev[[g]])
    min(vals) / max(vals)  # ratio min/max: 1 = perfect convergence
  })
  conv_dt <- data.table(gene = gene_cols, convergence_score = round(conv_score, 4))
  conv_dt <- conv_dt[order(-convergence_score)]
  cat("\nTop 10 genes by One Health convergence score (min/max prevalence ratio):\n")
  print(conv_dt[1:10])
  fwrite(conv_dt, "~/matrix_output_v2/one_health_convergence.csv")
}

# ── 7. SAVE OUTPUTS ─────────────────────────────────────────
fwrite(kw_dt,       "~/matrix_output_v2/one_health_kruskal.csv")
fwrite(comp_topo_dt,"~/matrix_output_v2/one_health_topology.csv")
fwrite(prev_by_comp,"~/matrix_output_v2/one_health_prevalence.csv")

cat("\nFiles saved:\n")
cat("  one_health_kruskal.csv\n")
cat("  one_health_topology.csv\n")
cat("  one_health_prevalence.csv\n")
cat("  one_health_convergence.csv\n")
cat("  network_human.rds\n")
cat("  network_animal.rds\n")
cat("  network_environment.rds\n")
cat("\nDone:", format(Sys.time()), "\n")
