# =========================================================
# Project      : Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE166409 (Illumina NextSeq 500)
# Samples      : 3 patients + 3 healthy controls, PBMC
# Condition    : WAC-ANKRD26 gene fusion-associated thrombocytopenia
# Script       : 03 - Gene Set Enrichment Analysis (GSEA)
# Input        : DEGs_<GSE>.csv (script 01)
# Ranking      : DESeq2 Wald statistic (all tested genes, unshrunken)
# Collections  : KEGG, Reactome, GO (BP, MF, CC; redundancy-reduced),
#                MSigDB Hallmark / C3 (TFT:GTRD) / C7 / C8
# Outputs      : GSEA tables (CSV), lollipop plots, Hallmark enrichment plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

# Packages, paths, cfg and helper functions (save_plot, save_csv, ...)
source(here::here("SetupEnvironment", "setup.R"))

set_seed(cfg$seed)

message("=== Starting ", cfg$gse_id, " GSEA analysis ===")


# =========================================================
#                  Functions
# =========================================================

# Convert a clusterProfiler gseaResult or an fgsea table into one common
# tibble: ID, pathway, NES, pval, padj, set_size, leading_genes (list column)
standardize_gsea <- function(obj) {

  if (inherits(obj, "gseaResult")) {
    tibble::as_tibble(obj@result) |>
      dplyr::transmute(
        ID            = ID,
        pathway       = Description,
        NES           = NES,
        pval          = pvalue,
        padj          = p.adjust,
        set_size      = setSize,
        leading_genes = strsplit(core_enrichment, "/", fixed = TRUE)
      )
  } else {                                   # fgsea output
    tibble::as_tibble(as.data.frame(obj)) |>
      dplyr::transmute(
        ID            = pathway,
        pathway       = pathway,
        NES           = NES,
        pval          = pval,
        padj          = padj,
        set_size      = size,
        leading_genes = leadingEdge
      )
  }
}

# Collapse the leading-edge list column into a "/"-separated string
# and sort by padj (table ready to be written to CSV)
export_gsea_table <- function(std) {
  std |>
    dplyr::mutate(leading_edge = vapply(leading_genes, paste, character(1), collapse = "/")) |>
    dplyr::select(-leading_genes) |>
    dplyr::arrange(padj)
}

# fgsea on a list of gene sets (same size limits as clusterProfiler);
# leading-edge genes converted from Entrez IDs to symbols
run_fgsea <- function(sets, stats, min_size, max_size, seed) {

  set.seed(seed)
  res <- fgsea::fgsea(pathways = sets, stats = stats,
                      minSize = min_size, maxSize = max_size, eps = 0)

  all_ids <- unique(unlist(res$leadingEdge))
  sym_map <- stats::setNames(entrez_to_symbol(all_ids), all_ids)
  res$leadingEdge <- lapply(res$leadingEdge, function(v) unname(sym_map[v]))

  res[order(res$pval), ]
}

# Hide only the GOSemSim "semData is not provided" notice
quiet_semdata <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("semData is not provided", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (grepl("semData is not provided", conditionMessage(m), fixed = TRUE))
        invokeRestart("muffleMessage")
    }
  )
}

# GSEA on one GO ontology (BP, MF or CC). Returns a list with:
#  - full:       all tested terms
#  - simplified: significant terms with redundancy removed separately for
#                up (NES > 0) and down (NES < 0) terms
# Returns NULL if the run fails or gives no result.
run_gsego <- function(gene_ranks, ont,
                      fdr = cfg$padj_cutoff,
                      simplify_cutoff = cfg$ora$simplify_cutoff) {

  set.seed(cfg$seed)

  full <- tryCatch(
    clusterProfiler::gseGO(
      geneList      = gene_ranks,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      minGSSize     = cfg$gsea$min_size,
      maxGSSize     = cfg$gsea$max_size,
      eps           = 0,
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      seed          = TRUE,
      verbose       = FALSE
    ),
    error = function(e) {
      warning("[GSEA] gseGO (", ont, ") failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(full) || nrow(full@result) == 0) {
    message("[GSEA] GO ", ont, ": no result.")
    return(NULL)
  }

  full <- make_readable(full)

  # Significant terms only (simplify() on all terms would be slow)
  sig <- full
  sig@result <- sig@result[!is.na(sig@result$p.adjust) & sig@result$p.adjust < fdr, , drop = FALSE]

  # Remove redundant terms from a gseaResult
  simp_one <- function(obj) {
    if (nrow(obj@result) < 2) return(obj)
    quiet_semdata(clusterProfiler::simplify(obj, cutoff = simplify_cutoff,
                                            by = "p.adjust", select_fun = min))
  }

  # Up and down terms are simplified separately
  up <- sig; up@result <- sig@result[sig@result$NES > 0, , drop = FALSE]
  dn <- sig; dn@result <- sig@result[sig@result$NES < 0, , drop = FALSE]

  simp_up <- simp_one(up)
  simp_dn <- simp_one(dn)

  message(sprintf("[GSEA] GO %s: %d significant terms -> %d after simplify (cutoff %.1f).",
                  ont, nrow(sig@result),
                  nrow(simp_up@result) + nrow(simp_dn@result), simplify_cutoff))

  list(
    full       = standardize_gsea(full),
    simplified = dplyr::bind_rows(standardize_gsea(simp_up), standardize_gsea(simp_dn))
  )
}

# Lollipop plot of the top_n up (NES > 0) and top_n down (NES < 0) significant
# pathways. Returns list(plot, summary), or NULL if nothing passes the filters.
gsea_lollipop <- function(std, color, title,
                          top_n = cfg$gsea$top_n,
                          fdr_cutoff = cfg$padj_cutoff,
                          min_leading = cfg$gsea$min_leading,
                          wrap_width = 45) {

  # Significant pathways with a large enough leading edge
  df <- std |>
    dplyr::mutate(leading_n = lengths(leading_genes)) |>
    dplyr::filter(!is.na(padj), padj < fdr_cutoff, leading_n >= min_leading)

  if (nrow(df) == 0) {
    message("[GSEA] No pathway passes FDR < ", fdr_cutoff,
            " and leading edge >= ", min_leading, ": plot skipped.")
    return(NULL)
  }

  # Top pathways of one direction, ordered by padj
  pick <- function(x) {
    x |>
      dplyr::arrange(padj, dplyr::desc(abs(NES))) |>
      dplyr::slice_head(n = top_n)
  }

  df_top <- dplyr::bind_rows(pick(dplyr::filter(df, NES > 0)),
                             pick(dplyr::filter(df, NES < 0))) |>
    dplyr::mutate(
      FDR_log10 = -log10(padj),
      pathway   = stringr::str_wrap(pathway, width = wrap_width),
      pathway   = stats::reorder(pathway, NES)
    )

  p <- ggplot(df_top, aes(x = NES, y = pathway)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), colour = "grey70") +
    geom_point(aes(size = leading_n, colour = FDR_log10)) +
    scale_colour_gradient(low = "grey60", high = color, name = "-log10(FDR)") +
    scale_size_continuous(name = "Leading edge size") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          panel.grid.major.y = element_blank())

  list(plot = p, summary = df_top)
}


# =========================================================
#                  Ranked gene list
# =========================================================

deg_table <- read_deg_table(file.path(paths$tables, sprintf("DEGs_%s.csv", cfg$gse_id)))

# Genes ranked by Wald statistic (one row per Entrez ID)
ranks <- deg_table |>
  dplyr::filter(!is.na(ENTREZID), is.finite(stat)) |>
  dplyr::distinct(ENTREZID, .keep_all = TRUE) |>
  dplyr::arrange(dplyr::desc(stat))

gene_ranks <- stats::setNames(ranks$stat, ranks$ENTREZID)

stopifnot(!anyDuplicated(names(gene_ranks)), length(gene_ranks) > 1000)

message(sprintf("[GSEA] %d ranked genes (Wald statistic).", length(gene_ranks)))
print(summary(gene_ranks))


# =========================================================
#                  GSEA - KEGG (clusterProfiler)
# =========================================================

options(timeout = 300)

# pvalueCutoff = 1 keeps the complete table; significance is applied later
kegg <- clusterProfiler::gseKEGG(
  geneList      = gene_ranks,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  minGSSize     = cfg$gsea$min_size,
  maxGSSize     = cfg$gsea$max_size,
  eps           = 0,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  seed          = TRUE,
  verbose       = FALSE
)
kegg <- make_readable(kegg)


# =========================================================
#                  GSEA - Reactome
# =========================================================

set.seed(cfg$seed)
reactome <- ReactomePA::gsePathway(
  geneList      = gene_ranks,
  organism      = "human",
  minGSSize     = cfg$gsea$min_size,
  maxGSSize     = cfg$gsea$max_size,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)
reactome <- make_readable(reactome)


# =========================================================
#      GSEA - GO (BP, MF, CC) with redundancy reduction
# =========================================================

go_res <- lapply(c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC"),
                 function(o) run_gsego(gene_ranks, o))
go_res <- Filter(Negate(is.null), go_res)

go_std_full <- lapply(go_res, `[[`, "full")
go_std_simp <- lapply(go_res, `[[`, "simplified")


# =========================================================
#                  GSEA - MSigDB (msigdbr >= 10, fgsea)
# =========================================================

# Collections to run: Hallmark, C3 (TF targets), C7 (immunologic), C8 (cell types)
msigdb_spec <- list(
  Hallmark    = list(collection = "H"),
  C3_TFT_GTRD = list(collection = "C3", subcollection = "TFT:GTRD"),
  C7_Immune   = list(collection = "C7"),
  C8_CellType = list(collection = "C8")
)

msigdb_sets <- lapply(msigdb_spec, function(s) do.call(get_msigdb_sets, s))

msigdb_res <- lapply(msigdb_sets, run_fgsea,
                     stats    = gene_ranks,
                     min_size = cfg$gsea$min_size,
                     max_size = cfg$gsea$max_size,
                     seed     = cfg$seed)


# =========================================================
#                  Standardise + lollipop plots
# =========================================================

# All collections in one list (GO = simplified terms)
gsea_std <- c(
  list(KEGG = standardize_gsea(kegg), Reactome = standardize_gsea(reactome)),
  go_std_simp,
  lapply(msigdb_res, standardize_gsea)
)

# Readable names for Hallmark only (it would alter TF symbols and GSE IDs elsewhere)
gsea_std$Hallmark$pathway <- clean_hallmark(gsea_std$Hallmark$pathway)

# Colour, plot width and title for each lollipop plot
plot_spec <- list(
  KEGG        = list(color = "#762A83",     width = 11, title = "KEGG pathways"),
  Reactome    = list(color = "darkorange3", width = 11, title = "Reactome pathways"),
  GO_BP       = list(color = "#1B7837",     width = 12, title = "GO Biological Process (simplified)"),
  GO_MF       = list(color = "#B35806",     width = 12, title = "GO Molecular Function (simplified)"),
  GO_CC       = list(color = "#542788",     width = 12, title = "GO Cellular Component (simplified)"),
  Hallmark    = list(color = "firebrick",   width = 11, title = "MSigDB Hallmark"),
  C3_TFT_GTRD = list(color = "steelblue",   width = 11, title = "MSigDB C3 TF targets (GTRD)"),
  C7_Immune   = list(color = "darkgreen",   width = 15, title = "MSigDB C7 immunologic signatures"),
  C8_CellType = list(color = "indianred4",  width = 11, title = "MSigDB C8 cell-type signatures")
)

# Keep only the collections that gave a result
plot_spec <- plot_spec[names(plot_spec) %in% names(gsea_std)]

gsea_plots <- Map(
  function(std, spec) gsea_lollipop(
    std, color = spec$color,
    title = paste0("GSEA - ", spec$title, " (", cfg$case_level, " vs ", cfg$ref_level, ")")),
  gsea_std[names(plot_spec)], plot_spec
)

# Number of significant sets per collection
for (nm in names(gsea_std)) {
  n_sig <- sum(gsea_std[[nm]]$padj < cfg$padj_cutoff, na.rm = TRUE)
  message(sprintf("[GSEA] %-12s %5d tested sets, %4d with FDR < %.2f",
                  nm, nrow(gsea_std[[nm]]), n_sig, cfg$padj_cutoff))
}


# =========================================================
#                  Enrichment-score plots (Hallmark, top 3)
# =========================================================

# Three most significant Hallmark sets
top_hallmark <- gsea_std$Hallmark |>
  dplyr::filter(!is.na(padj), padj < cfg$padj_cutoff) |>
  dplyr::arrange(padj) |>
  dplyr::slice_head(n = 3)

for (i in seq_len(nrow(top_hallmark))) {
  id <- top_hallmark$ID[i]

  p <- fgsea::plotEnrichment(msigdb_sets$Hallmark[[id]], gene_ranks) +
    labs(title = clean_hallmark(id),
         subtitle = sprintf("NES = %.2f, FDR = %.2g", top_hallmark$NES[i], top_hallmark$padj[i]))

  save_plot_safe(p, filename = paste0("GSEA_ES_Hallmark_", gsub("[^A-Za-z0-9]+", "_", id), ".png"),
                 dir = paths$plots_gsea, width = 7, height = 5)
}


# =========================================================
#                  Save tables and plots
# =========================================================

# One table and one lollipop plot per collection
for (nm in names(gsea_std)) {
  save_csv_safe(export_gsea_table(gsea_std[[nm]]),
                filename = sprintf("results_%s_GSEA_%s.csv", cfg$gse_id, nm),
                dir = paths$tables_gsea)

  save_plot_safe(gsea_plots[[nm]]$plot, filename = sprintf("GSEA_%s.png", nm),
                 dir = paths$plots_gsea, width = plot_spec[[nm]]$width, height = 9)
}

# Complete GO tables (before redundancy removal)
for (nm in names(go_std_full)) {
  save_csv_safe(export_gsea_table(go_std_full[[nm]]),
                filename = sprintf("results_%s_GSEA_%s_full.csv", cfg$gse_id, nm),
                dir = paths$tables_gsea)
}

save_session_info(filename = sprintf("sessionInfo_GSEA_%s.txt", cfg$gse_id), dir = paths$logs,
                  label = paste0("Gene Set Enrichment Analysis (GSEA) - ", cfg$gse_id))

message("=================================================")
message("[PIPELINE] GSEA completed for ", cfg$gse_id,
        " (", cfg$case_level, " vs ", cfg$ref_level, ").")
message("=================================================")
