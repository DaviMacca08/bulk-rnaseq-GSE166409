# =========================================================
# Project      : Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE166409 (Illumina NextSeq 500)
# Samples      : 3 patients + 3 healthy controls, PBMC
# Condition    : WAC-ANKRD26 gene fusion-associated thrombocytopenia
# Script       : 02 - Over-Representation Analysis (ORA)
# Input        : DEGs_<GSE>.csv and universe_tested_<GSE>.csv (script 01)
# Gene lists   : all / up / down DEGs, as defined in script 01
# Databases    : GO (BP, MF, CC; redundancy-reduced), KEGG, Reactome, Hallmark
# Outputs      : ORA tables (CSV) and dot plots per database and direction
#
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

# Packages, paths, cfg and helper functions (save_plot, save_csv, ...)
source(here::here("SetupEnvironment", "setup.R"))

set_seed(cfg$seed)

message("=== Starting ", cfg$gse_id, " ORA analysis ===")


# =========================================================
#                  Functions
# =========================================================

# Run GO (BP, MF, CC), KEGG, Reactome and Hallmark enrichment for one gene
# list. Returns a named list of enrichment results, or NULL if the gene list
# is too short.
run_ora <- function(genes, universe, term2gene_h) {

  if (length(genes) < cfg$ora$min_genes) {
    message("[ORA] Fewer than ", cfg$ora$min_genes, " genes: skipped.")
    return(NULL)
  }

  # Arguments shared by all enrichment calls
  common <- list(
    universe      = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = cfg$padj_cutoff,   # applied to the adjusted p-value
    qvalueCutoff  = 1,                 # q-value filter disabled
    minGSSize     = cfg$ora$min_size,
    maxGSSize     = cfg$ora$max_size
  )

  out <- list()

  # GO: one ontology at a time, then remove redundant terms
  for (ont in c("BP", "MF", "CC")) {
    r <- do.call(clusterProfiler::enrichGO,
                 c(list(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                        ont = ont, readable = TRUE), common))
    if (!is.null(r) && nrow(as.data.frame(r)) > 1) {
      r <- clusterProfiler::simplify(r, cutoff = cfg$ora$simplify_cutoff,
                                     by = "p.adjust", select_fun = min)
    }
    out[[paste0("GO_", ont)]] <- r
  }

  # KEGG (needs internet; a failure is reported as a warning)
  out$KEGG <- tryCatch(
    make_readable(do.call(clusterProfiler::enrichKEGG,
                          c(list(gene = genes, organism = "hsa",
                                 keyType = "ncbi-geneid"), common))),
    error = function(e) { warning("[ORA] KEGG failed: ", conditionMessage(e)); NULL }
  )

  # Reactome
  out$Reactome <- do.call(ReactomePA::enrichPathway,
                          c(list(gene = genes, organism = "human",
                                 readable = TRUE), common))

  # Hallmark (custom gene sets), with readable term names
  h <- do.call(clusterProfiler::enricher,
               c(list(gene = genes, TERM2GENE = term2gene_h), common))
  h <- make_readable(h)
  if (!is.null(h) && nrow(as.data.frame(h)) > 0) {
    h@result$Description <- clean_hallmark(h@result$ID)
  }
  out$Hallmark <- h

  out
}

# Dot plot of the top enriched terms: gene ratio on x, fold enrichment as
# size, FDR as colour. Returns NULL if no term passes the FDR cutoff.
plot_enrichment <- function(res, title, color,
                            top_n = cfg$ora$top_n, wrap_width = 40,
                            fdr_cutoff = cfg$padj_cutoff) {

  if (is.null(res)) return(NULL)

  df <- as.data.frame(res)
  df <- df[!is.na(df$p.adjust) & df$p.adjust < fdr_cutoff, , drop = FALSE]

  if (nrow(df) == 0) {
    message("[PLOT] No term with FDR < ", fdr_cutoff, ": ", title)
    return(NULL)
  }

  # "k/n" strings (GeneRatio, BgRatio) -> numeric
  ratio <- function(x) {
    vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
  }

  df$GeneRatioNum <- ratio(df$GeneRatio)
  if (!"FoldEnrichment" %in% colnames(df)) {
    df$FoldEnrichment <- df$GeneRatioNum / ratio(df$BgRatio)
  }

  # Top terms by FDR; the most significant one is drawn at the top
  top <- df |>
    dplyr::arrange(p.adjust, dplyr::desc(Count)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = wrap_width),
      Description = factor(Description, levels = rev(unique(Description)))
    )

  ggplot(top, aes(x = GeneRatioNum, y = Description)) +
    geom_point(aes(size = FoldEnrichment, colour = p.adjust)) +
    scale_colour_gradient(low = color, high = "#2C7BB6", labels = scales::label_scientific()) +
    scale_size_continuous(range = c(2, 6)) +
    labs(title = title, x = "Gene ratio", y = NULL,
         size = "Fold enrichment", colour = "FDR") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))
}


# =========================================================
#                  Load data
# =========================================================

# DEG table and universe of tested genes (script 01)
deg_table <- read_deg_table(file.path(paths$tables, sprintf("DEGs_%s.csv", cfg$gse_id)))
universe  <- read_universe(file.path(paths$tables, sprintf("universe_tested_%s.csv", cfg$gse_id)))

# The up / down call is made once, in script 01
stopifnot("significance" %in% colnames(deg_table))

# DEG lists (All / Up / Down) restricted to the universe
gene_lists <- list(
  All  = deg_table$ENTREZID[deg_table$significance != cfg$label_ns],
  Up   = deg_table$ENTREZID[deg_table$significance == cfg$label_up],
  Down = deg_table$ENTREZID[deg_table$significance == cfg$label_down]
) |>
  lapply(function(x) unique(intersect(x[!is.na(x)], universe)))

message(sprintf("[ORA] Universe: %d tested genes. DEGs: %d all, %d up, %d down.",
                length(universe), length(gene_lists$All),
                length(gene_lists$Up), length(gene_lists$Down)))


# =========================================================
#                  Hallmark gene sets (msigdbr >= 10)
# =========================================================

# Gene sets as a TERM2GENE table (term, Entrez ID) for enricher()
hallmark_sets <- get_msigdb_sets("H")

term2gene_h <- data.frame(
  gs_name     = rep(names(hallmark_sets), lengths(hallmark_sets)),
  entrez_gene = unlist(hallmark_sets, use.names = FALSE),
  stringsAsFactors = FALSE
)


# =========================================================
#                  Run ORA
# =========================================================

# One enrichment run per gene list (All, Up, Down)
ora <- lapply(gene_lists, run_ora, universe = universe, term2gene_h = term2gene_h)


# =========================================================
#                  Save tables and plots
# =========================================================

# Dot colour per gene list
direction_color <- c(All = "#C44E52", Up = "#009E73", Down = "#CC79A7")

# Database names used in plot titles
db_label <- c(GO_BP = "GO Biological Process", GO_MF = "GO Molecular Function",
              GO_CC = "GO Cellular Component", KEGG = "KEGG",
              Reactome = "Reactome", Hallmark = "MSigDB Hallmark")

for (dir_name in names(ora)) {
  if (is.null(ora[[dir_name]])) next

  for (db in names(ora[[dir_name]])) {
    res <- ora[[dir_name]][[db]]

    # Results table
    save_csv_safe(as.data.frame(res),
                  filename = sprintf("ORA_%s_%s_%s.csv", cfg$gse_id, db, dir_name),
                  dir = paths$tables_ora)

    # Dot plot
    p <- plot_enrichment(
      res, color = direction_color[[dir_name]],
      title = sprintf("%s - %s vs %s (%s DEGs)", db_label[[db]],
                      cfg$case_level, cfg$ref_level, dir_name)
    )

    save_plot_safe(p, filename = sprintf("ORA_%s_%s.png", db, dir_name),
                   dir = paths$plots_ora, width = 10, height = 8)
  }
}

save_session_info(filename = sprintf("sessionInfo_ORA_%s.txt", cfg$gse_id), dir = paths$logs,
                  label = paste0("Over-Representation Analysis (ORA) - ", cfg$gse_id))

message("=================================================")
message("[PIPELINE] ORA completed for ", cfg$gse_id,
        " (", cfg$case_level, " vs ", cfg$ref_level, ").")
message("=================================================")
