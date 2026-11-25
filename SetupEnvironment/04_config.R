# =========================================================
# Analysis configuration (parameters only, no functions)
# =========================================================
# Dataset-specific values are grouped first; generic analysis
# parameters follow. Every threshold used by the scripts lives here.

cfg <- list(

  # ---- Dataset / contrast (CHANGE PER PROJECT) ------------
  gse_id     = "GSE166409",
  ref_level  = "HC",           # reference group
  case_level = "THR",          # compared group
  counts_file = file.path(paths$data, "GSE166409_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
  case_pattern = "thrombocytopenia",   # keyword that identifies case samples in GEO metadata
  ref_pattern  = "healthy",            # keyword that identifies reference samples

  # ---- General -------------------------------------------
  seed = 1234,

  # ---- Pre-processing ------------------------------------
  min_count = 10,              # reads required in >= (smallest group size) samples
  exclude_sex_chr = TRUE,      # drop X/Y genes (justified only if sex is confounded)
  top_n_heatmap   = 150,       # genes shown in the DEG heatmap

  # ---- Differential expression ---------------------------
  padj_cutoff     = 0.05,
  lfc_cutoff      = 1,         # applied to the SHRUNKEN log2FC
  lfc_shrink_type = "apeglm",

  # ---- QC -------------------------------------------------
  min_n_outlier_rules = 8,     # automatic outlier flags only with >= this many samples

  # ---- GSEA -----------------------------------------------
  gsea = list(min_size = 15, max_size = 500, top_n = 10, min_leading = 5),

  # ---- ORA ------------------------------------------------
  ora  = list(min_size = 10, max_size = 500, min_genes = 5, top_n = 10,
              simplify_cutoff = 0.7),

  # ---- Advanced analysis (script 04, CHANGE PER PROJECT) ---
  advanced = list(
    genes_of_interest = c("ANKRD26", "WAC"),
    locus_genes       = c("ANKRD26", "WAC", "YME1L1", "MASTL", "ACBD5",
                          "RAB18", "MKX", "ARMC4", "MPP7", "BAMBI"),
    ridge_top_n       = 25,    # GO terms per direction in the ridge plot
    paper_terms       = paste0("myeloid leukocyte activation|leukocyte differentiation|",
                               "platelet (activation|aggregation|degranulation)|",
                               "response to chemokine")
  )
)

# Labels derived from the contrast (used for significance calls / plots)
cfg$label_up   <- paste("Up in",   cfg$case_level)
cfg$label_down <- paste("Down in", cfg$case_level)
cfg$label_ns   <- "Not significant"

# Palettes
group_palette <- stats::setNames(c("#CE2915", "#0096FF"),
                                 c(cfg$ref_level, cfg$case_level))

# Sanity checks: fail now, not halfway through an analysis
validate_cfg <- function(cfg) {
  stopifnot(
    is.character(cfg$gse_id), length(cfg$gse_id) == 1,
    cfg$ref_level != cfg$case_level,
    cfg$padj_cutoff > 0, cfg$padj_cutoff < 1,
    cfg$lfc_cutoff >= 0,
    cfg$min_count >= 1,
    cfg$gsea$min_size < cfg$gsea$max_size,
    cfg$ora$min_size  < cfg$ora$max_size
  )
  invisible(cfg)
}

validate_cfg(cfg)
