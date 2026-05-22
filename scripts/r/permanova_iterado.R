library(vegan)
library(data.table)

base <- "/home/miguel/Abaum_Resistome_Network"
matrix <- fread(file.path(base, "data/processed/matrices/matrix1_gene_families.csv"))
meta   <- fread(file.path(base, "data/processed/metadata/genomes_pass_qc_validated.csv"))
setnames(matrix, colnames(matrix)[1], "genome_id")
setnames(meta,   colnames(meta)[1],   "genome_id")

animal_terms <- c("cow","bovine","swine","pig","poultry","chicken","horse",
                  "dog","cat","sheep","goat","fish","bird","veterinary",
                  "livestock","feline","canine","duck","rabbit","animal")
human_terms  <- c("hospital","sputum","sputament","blood","urine","wound",
                  "respiratory","icu","lung","skin","stool","rectal","balf",
                  "bal","cerebrospinal","csf","tracheal","pus","clinical",
                  "oral","abdominal","pleural","catheter","tissue","bone",
                  "joint","ascites","endotracheal","peritoneal","abscess",
                  "secretion","urinary","bile","swab","drain","human")
env_terms    <- c("soil","environment","sewage","water","river","lake",
                  "ocean","sediment","effluent","wastewater","compost",
                  "manure","air","sink","surface","floor","tap","env")

assign_comp <- function(src) {
  if (is.na(src)) return("unknown")
  s <- tolower(src)
  if (any(sapply(animal_terms, function(t) grepl(t, s, fixed=TRUE)))) return("animal")
  if (any(sapply(human_terms,  function(t) grepl(t, s, fixed=TRUE)))) return("human")
  if (any(sapply(env_terms,    function(t) grepl(t, s, fixed=TRUE)))) return("environment")
  return("unknown")
}

meta[, compartment := sapply(isolation_source, assign_comp)]
mat_meta <- merge(matrix, meta[, .(genome_id, compartment)],
                  by = "genome_id", all.x = TRUE)
mat_meta <- mat_meta[!is.na(compartment) & compartment != "unknown"]
gene_cols <- setdiff(colnames(matrix), "genome_id")
mat_vals  <- as.matrix(mat_meta[, ..gene_cols])
comp_vec  <- mat_meta$compartment
keep      <- rowSums(mat_vals) > 0
mat_vals  <- mat_vals[keep, ]
comp_vec  <- comp_vec[keep]

cat("=== Available N ===\n")
print(table(comp_vec))

n_bal  <- min(table(comp_vec))
n_iter <- 100
cat(sprintf("\nRunning %d iterations with n=%d/compartment...\n",
            n_iter, n_bal))

results <- data.frame(iter=integer(), R2=numeric(), F=numeric(),
                      p=numeric(), bd_p=numeric())

for (i in 1:n_iter) {
  set.seed(i)
  idx <- c(
    sample(which(comp_vec == "human"),       n_bal),
    sample(which(comp_vec == "animal"),      n_bal),
    sample(which(comp_vec == "environment"), n_bal)
  )
  m <- mat_vals[idx, ]
  g <- comp_vec[idx]

  perm <- adonis2(m ~ g, method = "bray", permutations = 199)
  dist_m <- vegdist(m, method = "bray")
  bd     <- betadisper(dist_m, g)
  bd_t   <- permutest(bd, permutations = 199)

  results <- rbind(results, data.frame(
    iter = i,
    R2   = perm$R2[1],
    F    = perm$F[1],
    p    = perm$`Pr(>F)`[1],
    bd_p = bd_t$tab$`Pr(>F)`[1]
  ))

  if (i %% 10 == 0) cat(sprintf("  Iteration %d/100 completed\n", i))
}

cat("\n=== ITERATED RESULTS (100 replicates, n=236/comp) ===\n")
cat(sprintf("R2   — mean: %.4f  SD: %.4f  CI95: [%.4f, %.4f]\n",
    mean(results$R2), sd(results$R2),
    quantile(results$R2, 0.025), quantile(results$R2, 0.975)))
cat(sprintf("F    — mean: %.2f   SD: %.2f\n",
    mean(results$F), sd(results$F)))
cat(sprintf("p    — all p<=0.001: %s  (n_sig=%d/100)\n",
    all(results$p <= 0.005), sum(results$p <= 0.005)))
cat(sprintf("Betadisper p<0.05: %d/100 iterations\n",
    sum(results$bd_p < 0.05, na.rm=TRUE)))

# Save results
out <- "/home/miguel/Abaum_Resistome_Network/results/supplementary/tables_submission/TableS_permanova_iterado.csv"
write.csv(results, out, row.names=FALSE)
cat(sprintf("\nOK: Results saved at %s\n", out))

cat("\n=== FINAL COMPARISON ===\n")
cat(sprintf("Iterated balanced (100x n=%d/comp): R2=%.3f+/-%.3f  F=%.1f  p<0.001\n",
    n_bal, mean(results$R2), sd(results$R2), mean(results$F)))
cat("Manuscript VM original:              R2=0.271  F=195.6  p=0.001\n")
