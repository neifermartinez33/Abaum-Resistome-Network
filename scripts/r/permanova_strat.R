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

# Stratified max 500/comp (replicates VM approach)
max_n <- 500
set.seed(42)
idx_strat <- c(
  sample(which(comp_vec == "human"),       min(sum(comp_vec=="human"), max_n)),
  sample(which(comp_vec == "animal"),      min(sum(comp_vec=="animal"), max_n)),
  sample(which(comp_vec == "environment"), min(sum(comp_vec=="environment"), max_n))
)
mat_strat  <- mat_vals[idx_strat, ]
comp_strat <- comp_vec[idx_strat]

cat("\n=== Stratified N (max 500/comp) ===\n")
print(table(comp_strat))

set.seed(42)
cat("\nRunning stratified PERMANOVA...\n")
perm_strat <- adonis2(mat_strat ~ comp_strat,
                      method = "bray", permutations = 999)
cat("\n=== STRATIFIED PERMANOVA ===\n")
print(perm_strat)

dist_strat <- vegdist(mat_strat, method = "bray")
bd_strat   <- betadisper(dist_strat, comp_strat)
bd_test    <- permutest(bd_strat, permutations = 999)
cat("\n=== STRATIFIED BETADISPER ===\n")
print(bd_test)

cat("\n=== COMPARATIVE SUMMARY ===\n")
cat(sprintf("Stratified max500: R2=%.3f  F=%.1f  p=%.3f  betadisper_p=%.3f\n",
    perm_strat$R2[1], perm_strat$F[1],
    perm_strat$`Pr(>F)`[1], bd_test$tab$`Pr(>F)`[1]))
cat("Manuscript original:  R2=0.271  F=195.6  p=0.001  betadisper_p=0.539\n")
