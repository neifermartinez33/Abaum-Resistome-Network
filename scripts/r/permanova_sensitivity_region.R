library(vegan)
library(data.table)

base <- "/home/miguel/Abaum_Resistome_Network"
matrix <- fread(file.path(base, "data/processed/matrices/matrix1_gene_families.csv"))
meta   <- fread(file.path(base, "data/processed/metadata/genomes_pass_qc_validated.csv"))
setnames(matrix, colnames(matrix)[1], "genome_id")
setnames(meta,   colnames(meta)[1],   "genome_id")

# Compartment assignment
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

# Region assignment (UN Geoscheme)
assign_region <- function(country) {
  if (is.na(country)) return("Unknown")
  c <- trimws(country)
  east_asia   <- c("China","Japan","South Korea","Taiwan","Mongolia","Hong Kong","North Korea")
  south_asia  <- c("India","Pakistan","Nepal","Bangladesh","Sri Lanka","Afghanistan","Bhutan","Maldives")
  southeast   <- c("Thailand","Vietnam","Cambodia","Malaysia","Singapore","Indonesia","Philippines",
                   "Myanmar","Laos","Timor-Leste","Brunei","Viet Nam")
  central_asia<- c("Kazakhstan","Uzbekistan","Kyrgyzstan","Tajikistan","Turkmenistan")
  middle_east <- c("Saudi Arabia","Iran","Israel","Iraq","Jordan","Lebanon","Syria","Kuwait",
                   "Qatar","UAE","Oman","Yemen","Turkey","Bahrain","United Arab Emirates")
  europe      <- c("Germany","France","United Kingdom","Italy","Spain","Netherlands","Belgium",
                   "Switzerland","Austria","Poland","Czech Republic","Romania","Hungary","Greece",
                   "Serbia","Croatia","Bulgaria","Denmark","Sweden","Norway","Finland","Portugal",
                   "Russia","Ukraine","Bosnia and Herzegovina","Montenegro","Kosovo","Albania",
                   "Slovenia","Slovakia","Latvia","Lithuania","Estonia","Belarus","Moldova",
                   "North Macedonia","Luxembourg","Ireland","Iceland")
  north_am    <- c("United States","Canada","Mexico","Puerto Rico")
  south_am    <- c("Brazil","Argentina","Colombia","Chile","Peru","Venezuela","Ecuador",
                   "Bolivia","Paraguay","Uruguay","Guyana","Suriname","Guatemala",
                   "Honduras","Nicaragua","Costa Rica","Panama","El Salvador")
  africa      <- c("South Africa","Nigeria","Egypt","Kenya","Ghana","Tanzania","Ethiopia",
                   "Senegal","Togo","Benin","Sudan","Libya","Tunisia","Morocco","Algeria",
                   "Cameroon","Cote d'Ivoire","Uganda","Rwanda","Djibouti","Mozambique",
                   "Zambia","Zimbabwe","Madagascar","Angola","Congo","Somalia")
  oceania     <- c("Australia","New Zealand","Papua New Guinea","Fiji")

  if (c %in% east_asia)    return("East Asia")
  if (c %in% south_asia)   return("South Asia")
  if (c %in% southeast)    return("Southeast Asia")
  if (c %in% central_asia) return("Central Asia")
  if (c %in% middle_east)  return("Middle East")
  if (c %in% europe)       return("Europe")
  if (c %in% north_am)     return("North America")
  if (c %in% south_am)     return("South America")
  if (c %in% africa)       return("Africa")
  if (c %in% oceania)      return("Oceania")
  return("Other")
}

meta[, compartment := sapply(isolation_source, assign_comp)]
meta[, region      := sapply(country,          assign_region)]

mat_meta <- merge(matrix, meta[, .(genome_id, compartment, region)],
                  by = "genome_id", all.x = TRUE)
mat_meta <- mat_meta[!is.na(compartment) & compartment != "unknown"]
mat_meta <- mat_meta[region != "Unknown" & region != "Other"]

gene_cols <- setdiff(colnames(matrix), "genome_id")
mat_vals  <- as.matrix(mat_meta[, ..gene_cols])
comp_vec  <- mat_meta$compartment
reg_vec   <- mat_meta$region

# Remove zero-sum rows
keep     <- rowSums(mat_vals) > 0
mat_vals <- mat_vals[keep, ]
comp_vec <- comp_vec[keep]
reg_vec  <- reg_vec[keep]

cat("=== N by compartment ===\n")
print(table(comp_vec))
cat("\n=== N by region ===\n")
print(table(reg_vec))
cat("\n=== Cross-table compartment x region ===\n")
print(table(comp_vec, reg_vec))

# PERMANOVA 1: compartment only
set.seed(42)
cat("\nRunning PERMANOVA: compartment only...\n")
perm_comp <- adonis2(mat_vals ~ comp_vec,
                     method = "bray", permutations = 999)
cat("\n=== PERMANOVA: compartment only ===\n")
print(perm_comp)

# PERMANOVA 2: region only
set.seed(42)
cat("\nRunning PERMANOVA: region only...\n")
perm_reg <- adonis2(mat_vals ~ reg_vec,
                    method = "bray", permutations = 999)
cat("\n=== PERMANOVA: region only ===\n")
print(perm_reg)

# PERMANOVA 3: compartment + region (region partialled out)
set.seed(42)
cat("\nRunning PERMANOVA: region + compartment (partial)...\n")
perm_partial <- adonis2(mat_vals ~ reg_vec + comp_vec,
                        method = "bray", permutations = 999,
                        by = "margin")
cat("\n=== PERMANOVA: compartment | region (marginal) ===\n")
print(perm_partial)

# PERMANOVA 4: compartment + region (compartment partialled out)
set.seed(42)
perm_partial2 <- adonis2(mat_vals ~ comp_vec + reg_vec,
                         method = "bray", permutations = 999,
                         by = "margin")
cat("\n=== PERMANOVA: region | compartment (marginal) ===\n")
print(perm_partial2)

cat("\n=== SUMMARY: variance partitioning ===\n")
cat(sprintf("Compartment alone:              R2=%.3f  p=%.3f\n",
    perm_comp$R2[1], perm_comp$`Pr(>F)`[1]))
cat(sprintf("Region alone:                   R2=%.3f  p=%.3f\n",
    perm_reg$R2[1], perm_reg$`Pr(>F)`[1]))
cat(sprintf("Compartment | region:           R2=%.3f  p=%.3f\n",
    perm_partial$R2[which(rownames(perm_partial)=="comp_vec")],
    perm_partial$`Pr(>F)`[which(rownames(perm_partial)=="comp_vec")]))
cat(sprintf("Region | compartment:           R2=%.3f  p=%.3f\n",
    perm_partial2$R2[which(rownames(perm_partial2)=="reg_vec")],
    perm_partial2$`Pr(>F)`[which(rownames(perm_partial2)=="reg_vec")]))

# Save results
results <- data.frame(
  model = c("Compartment only","Region only",
            "Compartment | Region","Region | Compartment"),
  R2    = c(perm_comp$R2[1], perm_reg$R2[1],
            perm_partial$R2[which(rownames(perm_partial)=="comp_vec")],
            perm_partial2$R2[which(rownames(perm_partial2)=="reg_vec")]),
  F_stat = c(perm_comp$F[1], perm_reg$F[1],
             perm_partial$F[which(rownames(perm_partial)=="comp_vec")],
             perm_partial2$F[which(rownames(perm_partial2)=="reg_vec")]),
  pvalue = c(perm_comp$`Pr(>F)`[1], perm_reg$`Pr(>F)`[1],
             perm_partial$`Pr(>F)`[which(rownames(perm_partial)=="comp_vec")],
             perm_partial2$`Pr(>F)`[which(rownames(perm_partial2)=="reg_vec")])
)

out <- file.path(base, "results/supplementary/tables_submission/TableS_permanova_region.csv")
write.csv(results, out, row.names=FALSE)
cat(sprintf("\nOK: Results saved to %s\n", out))
