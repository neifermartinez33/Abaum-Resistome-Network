"""
Phase 7: Binary Presence/Absence Matrix Construction — Version 2 Final
Pipeline: Convergent evolution of the A. baumannii resistome — 20,739 genomes
Author: Neifer Miguel Martinez Durango
Date: April 2026

Design (Anti-Bias Decalogue — Final Version):
1. AMRFinderPlus backbone (--organism Acinetobacter_baumannii --plus)
   - Scope: core (ARGs) + plus (HMGs, biocides) — enables ARG-HMG co-occurrence analysis
   - Filters: identity >=90%, coverage >=60%
2. Pleiotropy fix: explode multi-class genes to long format before merge
   (e.g. msr(E) -> macrolide + streptogramin as separate rows)
3. Drug class harmonization: RGI subclasses collapsed to macro classes
   (carbapenem, cephalosporin -> beta-lactam antibiotic)
4. Asymmetric validation: AMRFinder gene confirmed if ANY of its drug classes
   detected by CARD-RGI OR ResFinder in same genome
   HMGs: single-source validation (AMRFinder only, documented as limitation)
5. OXA classifier: separates intrinsic OXA-51-like from acquired carbapenemases
   Word boundary anchors (\b) prevent prefix false positives
6. Collapse to gene family with OXA-aware classification
7. Two matrices: gene families (network, 5%-98%) + individual genes (secondary, 1%-98%)
8. Whitelist: WHO critical genes bypass prevalence filter
9. Core resistome reported separately (>98%)
10. Whitelist report with geo-temporal distribution

Methodology notes:
- HMGs included via Scope=='plus' (AMRFinder --plus flag)
  Single-source validation (no CARD-RGI/ResFinder equivalent for HMGs)
  Documented explicitly as limitation
- Efflux pumps (AdeABC, AdeIJK): excluded, AMRFinder does not catalog them
  for A. baumannii; CARD-RGI detects them only as Loose hits
"""

import logging
import re
from pathlib import Path

import pandas as pd
from argnorm.normalizers import AMRFinderPlusNormalizer, ResFinderNormalizer

# ─────────────────────────────────────────────
# 0. CONFIGURATION
# ─────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("/home/abaum-gcloud/04_binary_matrix_v2.log"),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)

# Paths
BASE      = Path("/home/abaum-gcloud/annotation")
AMR_DIR   = BASE / "amrfinderplus"
RGI_DIR   = BASE / "card_rgi"
RES_DIR   = BASE / "resfinder"
META_FILE = Path("/home/abaum-gcloud/genomes_pass_qc_validated.csv")
OUT_DIR   = Path("/home/abaum-gcloud/matrix_output_v2")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Quality thresholds
AMR_IDENTITY_THRESH = 90.0
AMR_COVERAGE_THRESH = 60.0
RES_IDENTITY_THRESH = 90.0
RES_COVERAGE_THRESH = 60.0

# Prevalence filters
PREV_MIN_MATRIX1 = 0.05   # 5%  — network matrix (gene families)
PREV_MIN_MATRIX2 = 0.01   # 1%  — secondary matrix (individual genes)
PREV_MAX         = 0.98   # 98% — above = core resistome

# Scope filter: core=ARGs, plus=HMGs/biocides/virulence
# Include both to enable ARG-HMG co-occurrence analysis (Fig. S2)
AMR_SCOPES = ["core", "plus"]

# ─────────────────────────────────────────────
# DRUG CLASS HARMONIZATION MAP
# Collapses RGI subclasses to AMRFinder macro classes
# so that drug-class-level merge works correctly
# ─────────────────────────────────────────────
DRUG_CLASS_MAP = {
    # Beta-lactam subclasses -> macro class
    "carbapenem":                          "beta-lactam antibiotic",
    "cephalosporin":                       "beta-lactam antibiotic",
    "penicillin beta-lactam":              "beta-lactam antibiotic",
    "monobactam":                          "beta-lactam antibiotic",
    "cephamycin":                          "beta-lactam antibiotic",
    "carbapenem; cephalosporin":           "beta-lactam antibiotic",
    "cephalosporin; penicillin beta-lactam": "beta-lactam antibiotic",
    # Aminoglycoside subclasses
    "aminoglycoside antibiotic":           "aminoglycoside antibiotic",
    # Fluoroquinolone
    "fluoroquinolone antibiotic":          "fluoroquinolone antibiotic",
    # Macrolide
    "macrolide antibiotic":                "macrolide antibiotic",
    # Tetracycline
    "tetracycline antibiotic":             "tetracycline antibiotic",
    "glycylcycline":                       "tetracycline antibiotic",
    # Sulfonamide
    "sulfonamide antibiotic":              "sulfonamide antibiotic",
    "diaminopyrimidine antibiotic":        "sulfonamide antibiotic",
    # Rifamycin
    "rifamycin antibiotic":                "rifamycin antibiotic",
    # Phenicol
    "phenicol antibiotic":                 "phenicol antibiotic",
    # Streptogramin
    "streptogramin antibiotic":            "streptogramin antibiotic",
    # Lincosamide
    "lincosamide antibiotic":              "lincosamide antibiotic",
    # Colistin/polymyxin
    "peptide antibiotic":                  "peptide antibiotic",
    # Fosfomycin
    "phosphonic acid antibiotic":          "phosphonic acid antibiotic",
    # Oxazolidinone
    "oxazolidinone antibiotic":            "oxazolidinone antibiotic",
    # Disinfectants/biocides
    "disinfecting agents and antiseptics": "disinfecting agents and antiseptics",
    # Metal resistance (HMGs — from AMRFinder plus scope)
    "mercury":           "mercury resistance",
    "arsenic":           "arsenic resistance",
    "copper":            "copper resistance",
    "silver":            "silver resistance",
    "nickel":            "nickel resistance",
    "tellurite":         "tellurite resistance",
    "chromium":          "chromium resistance",
    "cadmium":           "cadmium resistance",
    "zinc":              "zinc resistance",
    "cobalt":            "cobalt resistance",
    "lead":              "lead resistance",
}


# WHO critical genes whitelist — bypass prevalence filter
# Word boundary \b prevents prefix matches (blaOXA-23\b does NOT match blaOXA-234)
WHITELIST_PATTERNS = [
    r'^mcr',              # colistin resistance (plasmid-mediated, last resort)
    r'^tet\(X',           # tigecycline resistance
    r'^blaIMP',           # metallo-beta-lactamase
    r'^blaNDM',           # metallo-beta-lactamase (pandemic)
    r'^blaVIM',           # metallo-beta-lactamase
    r'^arr',              # rifamycin resistance
    r'^rpoB',             # rifamycin resistance (point mutation)
    r'^blaKPC',           # serine carbapenemase
    r'^blaOXA-23\b',      # OXA-23-like acquired carbapenemase
    r'^blaOXA-24\b',      # OXA-24/40-like acquired carbapenemase
    r'^blaOXA-40\b',      # OXA-24-like variant
    r'^blaOXA-58\b',      # OXA-58-like acquired carbapenemase
    r'^blaOXA-96\b',      # OXA-58-like variant
    r'^blaOXA-97\b',      # OXA-58-like variant
    r'^blaOXA-143\b',     # OXA-143-like (emerging in Brazil/Italy)
    r'^blaOXA-231\b',     # OXA-143-like variant
    r'^blaOXA-253\b',     # OXA-143-like (critical in South America)
    r'^blaOXA-235\b',     # OXA-235-like acquired carbapenemase
    r'^blaOXA-236\b',     # OXA-235-like variant
    r'^blaOXA-237\b',     # OXA-235-like variant
    r'^pmrB',             # colistin resistance (chromosomal mutation)
]


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────
def normalize_class(val: str) -> set:
    """Split and normalize drug class strings (comma or semicolon separated)."""
    if pd.isna(val):
        return set()
    classes = set()
    for sep in [',', ';']:
        if sep in str(val):
            for c in str(val).split(sep):
                c = c.strip().lower()
                if c:
                    classes.add(c)
            return classes
    return {str(val).strip().lower()}


def harmonize_drug_class(raw_class: str) -> str:
    """
    Collapse RGI subclasses to AMRFinder macro classes.
    Falls back to raw_class if not in map.
    """
    c = raw_class.strip().lower()
    return DRUG_CLASS_MAP.get(c, c)


def classify_oxa(symbol: str) -> str:
    """
    Classify blaOXA variants into biologically meaningful families.
    Separates intrinsic OXA-51-like (core genome, ~99%) from
    acquired carbapenemases (plasmid-mediated, variable).
    References:
    - Evans & Amyes (2014) FEMS Microbiol Rev
    - Poirel & Nordmann (2006) Clin Microbiol Infect
    - Higgins et al. (2010) — OXA-143-like
    - Dortet et al. (2015) — OXA-235-like
    """
    m = re.match(r'^blaOXA-?(\d+)', symbol, re.IGNORECASE)
    if not m:
        return "blaoxa_intrinsic"  # generic blaOXA -> intrinsic
    n = int(m.group(1))
    if n in {23, 225}:          return "blaoxa_23like"
    if n in {24, 40, 72, 160}:  return "blaoxa_24like"
    if n in {58, 96, 97}:       return "blaoxa_58like"
    if n in {143, 231, 253}:    return "blaoxa_143like"
    if n in {235, 236, 237}:    return "blaoxa_235like"
    if n in {2, 10, 21}:        return "blaoxa_narrow_acquired"
    return "blaoxa_intrinsic"


def extract_gene_family(symbol: str) -> str:
    """
    Derive gene family from AMRFinderPlus Element symbol.
    blaOXA variants use classify_oxa() for ontology-aware classification.
    """
    s = str(symbol).strip()
    if re.match(r'^blaOXA', s, re.IGNORECASE):
        return classify_oxa(s)
    if '_' in s and re.search(r'_[A-Z]\d+[A-Z]$', s):
        return s.split('_')[0].lower()
    m = re.match(r'^(bla[A-Za-z]+)-?\d*', s)
    if m:
        return m.group(1).lower()
    m = re.match(r'^([a-z]+\([^)]+\))', s, re.IGNORECASE)
    if m:
        return m.group(1).lower()
    m = re.match(r'^([a-zA-Z]+)\d+$', s)
    if m:
        return m.group(1).lower()
    return s.lower()


def is_whitelisted(symbol: str) -> bool:
    for pattern in WHITELIST_PATTERNS:
        if re.match(pattern, symbol, re.IGNORECASE):
            return True
    return False


# ─────────────────────────────────────────────
# 1. LOAD METADATA
# ─────────────────────────────────────────────
def load_metadata(meta_file: Path) -> pd.DataFrame:
    log.info("Loading genome metadata...")
    meta = pd.read_csv(meta_file, usecols=[
        "Name", "country", "collection_date",
        "isolation_source", "host_name", "source_db"
    ])
    meta = meta.rename(columns={"Name": "Genome_ID"})
    log.info(f"Metadata loaded: {len(meta)} genomes")
    return meta


# ─────────────────────────────────────────────
# 2. PARSE AMRFINDERPLUS (backbone)
# ─────────────────────────────────────────────
def parse_amrfinder(amr_dir: Path) -> pd.DataFrame:
    """
    Read all AMRFinderPlus .tsv files via argNorm.
    Include Scope==core (ARGs) AND Scope==plus (HMGs, biocides).
    FIX: explode multi-class genes (pleiotropy) to long format before merge.
    Each gene gets one row per drug class — ensures correct merge with validators.
    Apply quality filters: identity >=90%, coverage >=60%.
    """
    log.info("Parsing AMRFinderPlus (backbone — core + plus scopes)...")
    normalizer = AMRFinderPlusNormalizer()
    records = []

    for tsv in amr_dir.glob("*.tsv"):
        genome_id = tsv.stem
        try:
            df = normalizer.run(str(tsv))
        except Exception as e:
            log.warning(f"AMRFinder error {genome_id}: {e}")
            continue

        # Include core (ARGs) and plus (HMGs, biocides)
        df = df[df["Scope"].isin(AMR_SCOPES)]
        # Apply quality thresholds
        df = df[
            (df["% Identity to reference"] >= AMR_IDENTITY_THRESH) &
            (df["% Coverage of reference"] >= AMR_COVERAGE_THRESH)
        ]
        df = df.dropna(subset=["resistance_to_drug_classes_names"])
        if df.empty:
            continue

        df = df[["Element symbol", "Scope", "resistance_to_drug_classes_names"]].copy()
        df["Genome_ID"] = genome_id
        records.append(df)

    amr = pd.concat(records, ignore_index=True)
    amr = amr.rename(columns={
        "Element symbol": "gene_symbol",
        "resistance_to_drug_classes_names": "drug_class_raw"
    })

    # FIX: explode multi-class genes to long format (pleiotropy fix)
    # e.g. msr(E) with 'macrolide antibiotic,streptogramin antibiotic'
    # becomes two rows: one for macrolide, one for streptogramin
    amr["drug_class_list"] = amr["drug_class_raw"].apply(
        lambda x: list(normalize_class(x)) if normalize_class(x) else []
    )
    amr = amr.explode("drug_class_list")
    amr = amr.rename(columns={"drug_class_list": "drug_class"})
    amr = amr.dropna(subset=["drug_class"])
    amr["drug_class"] = amr["drug_class"].str.strip().str.lower()
    # Harmonize AMRFinder classes for symmetry with RGI and ResFinder
    amr["drug_class"] = amr["drug_class"].apply(harmonize_drug_class)

    # Add gene family and whitelist flag
    amr["gene_family"]  = amr["gene_symbol"].apply(extract_gene_family)
    amr["whitelisted"]  = amr["gene_symbol"].apply(is_whitelisted)
    amr["is_hmg"]       = amr["Scope"] == "plus"

    amr = amr[["Genome_ID", "gene_symbol", "gene_family",
               "drug_class", "whitelisted", "is_hmg"]].drop_duplicates()

    n_args = amr[~amr["is_hmg"]]["gene_symbol"].nunique()
    n_hmgs = amr[amr["is_hmg"]]["gene_symbol"].nunique()
    log.info(f"AMRFinderPlus: {amr['Genome_ID'].nunique()} genomes, "
             f"{n_args} unique ARGs, {n_hmgs} unique HMGs")
    wl = amr[amr['whitelisted']]['gene_symbol'].unique().tolist()
    log.info(f"Whitelisted genes detected: {wl if wl else 'none'}")
    return amr


# ─────────────────────────────────────────────
# 3. PARSE CARD-RGI (drug class validator)
# ─────────────────────────────────────────────
def parse_rgi_classes(rgi_dir: Path) -> pd.DataFrame:
    """
    Extract drug classes per genome from CARD-RGI (all Cut_Off).
    Apply harmonization to collapse subclasses to macro classes
    (carbapenem -> beta-lactam antibiotic) for correct merge.
    """
    log.info("Parsing CARD-RGI drug classes (validator)...")
    records = []

    for txt in rgi_dir.glob("*.txt"):
        genome_id = txt.stem
        try:
            df = pd.read_csv(txt, sep="\t", usecols=["Drug Class"], low_memory=False)
        except Exception as e:
            log.warning(f"RGI error {genome_id}: {e}")
            continue

        df = df.dropna(subset=["Drug Class"])
        if df.empty:
            continue

        classes = set()
        for val in df["Drug Class"]:
            for c in normalize_class(val):
                classes.add(harmonize_drug_class(c))

        for c in classes:
            records.append({"Genome_ID": genome_id, "drug_class": c})

    rgi_classes = pd.DataFrame(records).drop_duplicates()
    rgi_classes["in_rgi"] = True
    log.info(f"RGI: {rgi_classes['Genome_ID'].nunique()} genomes, "
             f"{rgi_classes['drug_class'].nunique()} unique drug classes")
    return rgi_classes


# ─────────────────────────────────────────────
# 4. PARSE RESFINDER (drug class validator)
# ─────────────────────────────────────────────
def parse_resfinder_classes(res_dir: Path) -> pd.DataFrame:
    """
    Extract drug classes per genome from ResFinder via argNorm.
    Apply harmonization for consistent merge with AMRFinder.
    """
    log.info("Parsing ResFinder drug classes (validator)...")
    normalizer = ResFinderNormalizer()
    records = []

    for genome_dir in res_dir.iterdir():
        if not genome_dir.is_dir():
            continue
        genome_id = genome_dir.name
        tab_file  = genome_dir / "ResFinder_results_tab.txt"
        if not tab_file.exists():
            continue

        try:
            df = normalizer.run(str(tab_file))
        except Exception as e:
            log.warning(f"ResFinder error {genome_id}: {e}")
            continue

        df = df[
            (df["Identity"] >= RES_IDENTITY_THRESH) &
            (df["Coverage"] >= RES_COVERAGE_THRESH)
        ]
        df = df.dropna(subset=["resistance_to_drug_classes_names"])
        if df.empty:
            continue

        classes = set()
        for val in df["resistance_to_drug_classes_names"]:
            for c in normalize_class(val):
                classes.add(harmonize_drug_class(c))

        for c in classes:
            records.append({"Genome_ID": genome_id, "drug_class": c})

    res_classes = pd.DataFrame(records).drop_duplicates()
    res_classes["in_res"] = True
    log.info(f"ResFinder: {res_classes['Genome_ID'].nunique()} genomes, "
             f"{res_classes['drug_class'].nunique()} unique drug classes")
    return res_classes


# ─────────────────────────────────────────────
# 5. ASYMMETRIC VALIDATION
# ─────────────────────────────────────────────
def validate_genes(amr: pd.DataFrame,
                   rgi_classes: pd.DataFrame,
                   res_classes: pd.DataFrame) -> pd.DataFrame:
    """
    Asymmetric validation at drug class level:
    ARG confirmed PRESENT if its drug_class detected by RGI OR ResFinder
    in same genome. Whitelisted genes bypass validation.
    HMGs (is_hmg=True): single-source validation (AMRFinder only).
    amr is in long format (exploded) — one row per gene-class pair.
    Merge is 1-to-1 on Genome_ID + drug_class.
    """
    log.info("Applying asymmetric validation...")

    # Split ARGs and HMGs
    args = amr[~amr["is_hmg"]].copy()
    hmgs = amr[amr["is_hmg"]].copy()

    # ── Validate ARGs via RGI/ResFinder drug class merge ──
    rgi_classes = rgi_classes.copy()
    res_classes = res_classes.copy()

    merged = args.merge(rgi_classes, on=["Genome_ID", "drug_class"], how="left")
    merged = merged.merge(res_classes, on=["Genome_ID", "drug_class"], how="left")
    merged["in_rgi"] = merged["in_rgi"].fillna(False)
    merged["in_res"] = merged["in_res"].fillna(False)
    merged["confirmed"] = merged["in_rgi"] | merged["in_res"] | merged["whitelisted"]

    confirmed_args = merged[merged["confirmed"]].copy()
    confirmed_args = confirmed_args[["Genome_ID", "gene_symbol", "gene_family",
                                     "drug_class", "whitelisted", "is_hmg"]].drop_duplicates()

    # ── HMGs: single-source validation ──
    # Documented limitation: no equivalent HMG database in CARD-RGI/ResFinder
    hmgs["confirmed"] = True  # AMRFinder with --plus is sufficient for HMGs
    confirmed_hmgs = hmgs[["Genome_ID", "gene_symbol", "gene_family",
                            "drug_class", "whitelisted", "is_hmg"]].drop_duplicates()

    confirmed = pd.concat([confirmed_args, confirmed_hmgs], ignore_index=True)
    confirmed = confirmed.drop_duplicates(subset=["Genome_ID", "gene_symbol"])

    n_total     = amr[~amr["is_hmg"]]["gene_symbol"].nunique()
    n_confirmed = confirmed_args["gene_symbol"].nunique()
    n_wl        = confirmed[confirmed["whitelisted"]]["gene_symbol"].nunique()
    n_hmg       = confirmed_hmgs["gene_symbol"].nunique()

    log.info(f"ARGs validated: {n_confirmed}/{n_total} unique gene types")
    log.info(f"HMGs included (single-source): {n_hmg} unique gene types")
    log.info(f"Whitelisted genes: {n_wl}")
    return confirmed


# ─────────────────────────────────────────────
# 6. WHITELIST DETAILED REPORT
# ─────────────────────────────────────────────
def generate_whitelist_report(confirmed: pd.DataFrame,
                               meta: pd.DataFrame,
                               n_genomes: int) -> pd.DataFrame:
    """
    Generate geo-temporal distribution report for WHO critical genes.
    """
    log.info("Generating whitelist genes detailed report...")

    wl = confirmed[confirmed["whitelisted"]].copy()
    if wl.empty:
        log.info("No whitelisted genes detected.")
        return pd.DataFrame()

    wl_meta = wl.merge(meta, on="Genome_ID", how="left")
    wl_meta["year"] = pd.to_datetime(
        wl_meta["collection_date"], errors="coerce"
    ).dt.year

    summary_rows = []
    for gene, grp in wl_meta.groupby("gene_symbol"):
        n    = grp["Genome_ID"].nunique()
        prev = n / n_genomes
        years = grp["year"].dropna()
        year_range = f"{int(years.min())}-{int(years.max())}" if not years.empty else "unknown"
        top_countries  = (grp.groupby("country")["Genome_ID"].nunique()
                          .sort_values(ascending=False).head(5).to_dict())
        iso_sources    = (grp.groupby("isolation_source")["Genome_ID"].nunique()
                          .sort_values(ascending=False).head(5).to_dict())
        source_dbs     = grp.groupby("source_db")["Genome_ID"].nunique().to_dict()
        summary_rows.append({
            "gene_symbol":       gene,
            "gene_family":       grp["gene_family"].iloc[0],
            "n_genomes":         n,
            "prevalence":        round(prev, 6),
            "year_range":        year_range,
            "top_countries":     str(top_countries),
            "isolation_sources": str(iso_sources),
            "source_dbs":        str(source_dbs)
        })

    report = pd.DataFrame(summary_rows).sort_values("n_genomes", ascending=False)
    report.to_csv(OUT_DIR / "whitelist_genes_report.csv", index=False)
    log.info(f"Whitelisted genes detected: {len(report)}")
    for _, row in report.iterrows():
        log.info(f"  {row['gene_symbol']}: {row['n_genomes']} genomes "
                 f"({row['prevalence']*100:.3f}%) | {row['year_range']} | "
                 f"{row['top_countries']}")
    return report


# ─────────────────────────────────────────────
# 7. PREVALENCE FILTER + WHITELIST
# ─────────────────────────────────────────────
def apply_prevalence_filter(df: pd.DataFrame,
                             id_col: str,
                             n_genomes: int,
                             prev_min: float,
                             label: str) -> tuple:
    """
    Split: retained | core_resistome (>98%) | whitelisted_rescue | excluded_rare.
    HMGs bypass validation but still subject to prevalence filter.
    """
    log.info(f"Applying prevalence filter for {label} ({prev_min*100:.0f}%-{PREV_MAX*100:.0f}%)...")

    prevalence    = df.groupby(id_col)["Genome_ID"].nunique() / n_genomes
    whitelist_ids = set(df[df["whitelisted"]][id_col].unique())
    core_ids      = set(prevalence[prevalence > PREV_MAX].index)
    rare_ids      = set(prevalence[prevalence < prev_min].index)
    retained_ids  = set(prevalence[(prevalence >= prev_min) & (prevalence <= PREV_MAX)].index)
    rescued_ids   = rare_ids & whitelist_ids
    retained_ids  = retained_ids | rescued_ids
    excluded_ids  = rare_ids - rescued_ids

    retained = df[df[id_col].isin(retained_ids)].copy()
    core_df  = df[df[id_col].isin(core_ids)].copy()

    prev_df = prevalence.reset_index()
    prev_df.columns = [id_col, "prevalence"]
    prev_df["n_genomes"] = (prev_df["prevalence"] * n_genomes).round().astype(int)
    prev_df["status"] = prev_df[id_col].apply(
        lambda x: "core_resistome"    if x in core_ids
        else ("retained"              if x in retained_ids and x not in rescued_ids
        else ("whitelisted_rescue"    if x in rescued_ids
        else  "excluded_rare"))
    )
    prev_df["whitelisted"] = prev_df[id_col].isin(whitelist_ids)

    log.info(f"  Retained: {len(retained_ids)} | Core (>98%): {len(core_ids)} | "
             f"Excluded rare: {len(excluded_ids)} | Whitelist rescued: {len(rescued_ids)}")
    return retained, core_df, prev_df


# ─────────────────────────────────────────────
# 8. BUILD BINARY MATRIX
# ─────────────────────────────────────────────
def build_matrix(df: pd.DataFrame,
                 id_col: str,
                 all_genomes: list,
                 label: str) -> pd.DataFrame:
    """Pivot to wide format (rows=genomes, cols=genes/families). int8 dtype."""
    log.info(f"Building binary matrix: {label}...")
    df = df.copy()
    df["present"] = 1
    matrix = df.pivot_table(
        index="Genome_ID", columns=id_col,
        values="present", aggfunc="max", fill_value=0
    ).astype("int8")
    missing = set(all_genomes) - set(matrix.index)
    if missing:
        empty = pd.DataFrame(0, index=list(missing),
                             columns=matrix.columns, dtype="int8")
        matrix = pd.concat([matrix, empty]).sort_index()
    matrix.index.name  = "Genome_ID"
    matrix.columns.name = None
    log.info(f"{label}: {matrix.shape[0]} x {matrix.shape[1]} | "
             f"Memory: {matrix.memory_usage(deep=True).sum()/1e6:.1f} MB")
    return matrix


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
if __name__ == "__main__":

    log.info("=" * 60)
    log.info("PHASE 7 v2 FINAL: BINARY MATRIX CONSTRUCTION")
    log.info("AMRFinder backbone (core+plus) + drug-class validation")
    log.info("Pleiotropy fix: explode multi-class genes")
    log.info("Drug class harmonization: RGI subclasses -> macro classes")
    log.info("HMGs included (single-source) for ARG-HMG co-occurrence")
    log.info("=" * 60)

    meta        = load_metadata(META_FILE)
    amr         = parse_amrfinder(AMR_DIR)
    rgi_classes = parse_rgi_classes(RGI_DIR)
    res_classes = parse_resfinder_classes(RES_DIR)

    all_genomes = list(set(
        amr["Genome_ID"].unique().tolist() +
        rgi_classes["Genome_ID"].unique().tolist() +
        res_classes["Genome_ID"].unique().tolist()
    ))
    n_genomes = len(all_genomes)
    log.info(f"Total genomes in analysis: {n_genomes}")

    confirmed   = validate_genes(amr, rgi_classes, res_classes)
    wl_report   = generate_whitelist_report(confirmed, meta, n_genomes)

    # ── Matrix 1: Gene families (network inference) ──
    log.info("-" * 40)
    fam = confirmed[["Genome_ID", "gene_family", "whitelisted"]].drop_duplicates()
    fam_retained, fam_core, fam_prev = apply_prevalence_filter(
        fam, "gene_family", n_genomes, PREV_MIN_MATRIX1, "gene_families"
    )
    matrix1 = build_matrix(fam_retained, "gene_family", all_genomes, "Matrix1_families")
    matrix1.to_csv(OUT_DIR / "matrix1_gene_families.csv")
    fam_prev.to_csv(OUT_DIR / "matrix1_prevalence.csv", index=False)

    # ── Matrix 2: Individual genes (secondary analysis) ──
    log.info("-" * 40)
    gene = confirmed[["Genome_ID", "gene_symbol", "whitelisted"]].drop_duplicates()
    gene_retained, gene_core, gene_prev = apply_prevalence_filter(
        gene, "gene_symbol", n_genomes, PREV_MIN_MATRIX2, "individual_genes"
    )
    matrix2 = build_matrix(gene_retained, "gene_symbol", all_genomes, "Matrix2_genes")
    matrix2.to_csv(OUT_DIR / "matrix2_individual_genes.csv")
    gene_prev.to_csv(OUT_DIR / "matrix2_prevalence.csv", index=False)

    # ── Core resistome report ──
    if not fam_core.empty:
        core_report = (fam_core.groupby("gene_family")["Genome_ID"]
                       .nunique().reset_index())
        core_report.columns = ["gene_family", "n_genomes"]
        core_report["prevalence"] = core_report["n_genomes"] / n_genomes
        core_report = core_report.sort_values("prevalence", ascending=False)
        core_report.to_csv(OUT_DIR / "core_resistome_report.csv", index=False)
        log.info(f"Core resistome (>98%): {len(core_report)} families")

    log.info("=" * 60)
    log.info("PHASE 7 v2 FINAL COMPLETE")
    log.info(f"Matrix 1 (network):   {matrix1.shape[0]} x {matrix1.shape[1]} gene families")
    log.info(f"Matrix 2 (secondary): {matrix2.shape[0]} x {matrix2.shape[1]} individual genes")
    log.info(f"Outputs: {OUT_DIR}")
    log.info("=" * 60)
