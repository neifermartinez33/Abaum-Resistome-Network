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

cat("=== Available N by compartment ===\n")
print(table(comp_vec))

# STRICT BALANCED: n = minimum of the three groups
n_bal <- min(table(comp_vec))
cat(sprintf("\nStrict balanced N: %d per compartment\n", n_bal))

set.seed(42)
idx_bal <- c(
  sample(which(comp_vec == "human"),       n_bal),
  sample(which(comp_vec == "animal"),      n_bal),
  sample(which(comp_vec == "environment"), n_bal)
)
mat_bal  <- mat_vals[idx_bal, ]
comp_bal <- comp_vec[idx_bal]

cat("\n=== Final balanced N ===\n")
print(table(comp_bal))
cat(sprintf("Total: %d genomes\n", nrow(mat_bal)))

# BALANCED PERMANOVA
set.seed(42)
cat("\nRunning balanced PERMANOVA...\n")
perm_bal <- adonis2(mat_bal ~ comp_bal,
                    method = "bray", permutations = 999)
cat("\n=== BALANCED PERMANOVA (n=236/compartment) ===\n")
print(perm_bal)

# BETADISPER
dist_bal <- vegdist(mat_bal, method = "bray")
bd_bal   <- betadisper(dist_bal, comp_bal)
bd_test  <- permutest(bd_bal, permutations = 999)
cat("\n=== BETADISPER BALANCED ===\n")
print(bd_test)

# FULL PERMANOVA for comparison
cat("\nRunning full PERMANOVA (all data)...\n")
set.seed(42)
perm_full <- adonis2(mat_vals ~ comp_vec,
                     method = "bray", permutations = 999)
cat("\n=== FULL PERMANOVA (unbalanced) ===\n")
print(perm_full)

dist_full <- vegdist(mat_vals, method = "bray")
bd_full   <- betadisper(dist_full, comp_vec)
bd_full_t <- permutest(bd_full, permutations = 999)

cat("\n=== FINAL SUMMARY ===\n")
cat(sprintf("Balanced (n=%d/comp):  R2=%.3f  F=%.2f  p=%.3f  betadisper_p=%.3f\n",
    n_bal,
    perm_bal$R2[1], perm_bal$F[1],
    perm_bal$`Pr(>F)`[1], bd_test$tab$`Pr(>F)`[1]))
cat(sprintf("Full (n=%d):        R2=%.3f  F=%.2f  p=%.3f  betadisper_p=%.3f\n",
    nrow(mat_vals),
    perm_full$R2[1], perm_full$F[1],
    perm_full$`Pr(>F)`[1], bd_full_t$tab$`Pr(>F)`[1]))
cat("Manuscript VM original: R2=0.271  F=195.6  p=0.001  betadisper_p=0.539\n")
