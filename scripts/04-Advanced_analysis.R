# =========================================================
# Project      : Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE166409 (Illumina NextSeq 500)
# Samples      : 3 patients + 3 healthy controls, PBMC
# Condition    : WAC-ANKRD26 gene fusion-associated thrombocytopenia
# Script       : 04 - Advanced analysis (comparison with the associated paper)
# Input        : outputs of scripts 01 and 03 (DEG table, DESeq2 object,
#                GO-BP GSEA tables)
# Analyses     : 1) expression of the genes of interest (ANKRD26, WAC)
#                2) expression of the other genes in the SV locus
#                3) ridge plot of the GO-BP GSEA leading-edge genes
#                4) check of the paper's GO terms (e.g. myeloid signature)
# Outputs      : gene tables (CSV), expression plot, ridge plot, term table
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

# Packages, paths, cfg and helper functions (save_plot, save_csv, ...)
source(here::here("SetupEnvironment", "setup.R"))

set_seed(cfg$seed)

message("=== Starting ", cfg$gse_id, " advanced analysis ===")


# =========================================================
#                  Functions
# =========================================================

# Ridge plot of the shrunken log2FC of the leading-edge genes in the top
# GO terms of each direction. The label shows the number of genes per term.
# Returns NULL if no term passes the filters.
gsea_ridgeplot <- function(std, lfc_table,
                           top_n = cfg$advanced$ridge_top_n,
                           fdr_cutoff = cfg$padj_cutoff,
                           min_genes = cfg$gsea$min_leading,
                           title = "GSEA - GO Biological Process") {

  # Top terms in each direction (up: NES > 0, down: NES < 0)
  sel <- std |>
    dplyr::filter(!is.na(padj), padj < fdr_cutoff) |>
    dplyr::mutate(direction = ifelse(NES > 0, "up", "down")) |>
    dplyr::group_by(direction) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  # One shrunken log2FC per gene symbol
  lfc_map <- lfc_table |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
    dplyr::select(SYMBOL, log2FC_shrunk)

  # One row per (term, leading-edge gene) with its log2FC; terms with fewer
  # than 'min_genes' genes are dropped
  core <- sel |>
    dplyr::select(ID, pathway, NES, padj, leading_genes) |>
    tidyr::unnest(leading_genes) |>
    dplyr::rename(SYMBOL = leading_genes) |>
    dplyr::inner_join(lfc_map, by = "SYMBOL") |>
    dplyr::add_count(pathway, name = "n_genes") |>
    dplyr::filter(n_genes >= min_genes) |>
    dplyr::mutate(label = stringr::str_wrap(paste0(pathway, " (n=", n_genes, ")"), 45))

  if (nrow(core) == 0) return(NULL)

  ggplot(core, aes(x = log2FC_shrunk, y = stats::reorder(label, NES), fill = padj)) +
    ggridges::geom_density_ridges(scale = 1.1, alpha = 0.8, colour = "white",
                                  quantile_lines = TRUE, quantiles = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_gradient(low = "#B2182B", high = "#2166AC", name = "Adjusted P-value") +
    labs(title = title,
         subtitle = sprintf("Shrunken log2FC of leading-edge genes; top %d terms per direction (n = genes shown)", top_n),
         x = sprintf("Shrunken log2 fold change (%s vs %s)", cfg$case_level, cfg$ref_level),
         y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 9),
          plot.title  = element_text(face = "bold"))
}


# =========================================================
#                  Load data
# =========================================================

# Input files written by scripts 01 and 03
deg_file  <- file.path(paths$tables, sprintf("DEGs_%s.csv", cfg$gse_id))
dds_file  <- file.path(paths$results, "DESeq2_object.rds")
simp_file <- file.path(paths$tables_gsea, sprintf("results_%s_GSEA_GO_BP.csv", cfg$gse_id))
full_file <- file.path(paths$tables_gsea, sprintf("results_%s_GSEA_GO_BP_full.csv", cfg$gse_id))

input_files <- c(deg_file, dds_file, simp_file, full_file)
missing     <- input_files[!file.exists(input_files)]

if (length(missing) > 0) {
  stop("Input file(s) not found:\n", paste(missing, collapse = "\n"),
       "\nRun scripts 01-DEG and 03-GSEA first.", call. = FALSE)
}

# DEG table
deg_table <- read_deg_table(deg_file)

# DESeq2 object: normalised counts and sample table
dds         <- readRDS(dds_file)
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
coldata     <- as.data.frame(SummarizedExperiment::colData(dds))

# GO-BP GSEA, simplified terms: rebuild the leading-edge gene list per term
gobp_simp <- utils::read.csv(simp_file, stringsAsFactors = FALSE)
gobp_simp$leading_genes <- strsplit(gobp_simp$leading_edge, "/", fixed = TRUE)

# GO-BP GSEA, complete table (before redundancy removal)
gobp_full <- utils::read.csv(full_file, stringsAsFactors = FALSE)


# =========================================================
#                  Gene ranking table
# =========================================================

# Rank of every gene by padj and by |Wald statistic|, and unshrunken fold change
deg_ranked <- deg_table |>
  dplyr::mutate(
    rank_padj       = rank(padj, ties.method = "min", na.last = "keep"),
    rank_abs_stat   = rank(-abs(stat), ties.method = "min", na.last = "keep"),
    fold_change_raw = 2^log2FC_raw
  ) |>
  dplyr::select(SYMBOL, baseMean, log2FC_raw, fold_change_raw, log2FC_shrunk,
                padj, rank_padj, rank_abs_stat)


# =========================================================
#                  1) Genes of interest
# =========================================================

genes_of_interest <- cfg$advanced$genes_of_interest

# Effect size, significance and rank of each gene
goi_table <- deg_ranked |>
  dplyr::filter(SYMBOL %in% genes_of_interest) |>
  dplyr::arrange(padj)

stopifnot(all(genes_of_interest %in% goi_table$SYMBOL))
print(goi_table)

save_csv(goi_table, filename = sprintf("genes_of_interest_%s.csv", cfg$gse_id),
         dir = paths$tables)

# Normalised counts of the genes of interest, one row per sample and gene
goi_ids <- deg_table$ENTREZID[match(genes_of_interest, deg_table$SYMBOL)]

goi_counts <- norm_counts[goi_ids, , drop = FALSE]
rownames(goi_counts) <- genes_of_interest

goi_df <- tibble::as_tibble(t(goi_counts), rownames = "sample") |>
  dplyr::mutate(group = coldata[sample, "condition"]) |>
  tidyr::pivot_longer(dplyr::all_of(genes_of_interest),
                      names_to = "gene", values_to = "norm_count")

# Expression per sample (one panel per gene, log10 scale)
goi_plot <- ggplot(goi_df, aes(group, norm_count + 1, colour = group)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.12, size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_y_log10() +
  scale_colour_manual(values = group_palette) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = paste0(paste(genes_of_interest, collapse = " and "), " expression per sample"),
       x = NULL, y = "Normalised counts + 1 (log10)", colour = NULL) +
  theme_bw()

save_plot(goi_plot,
          filename = sprintf("check_%s.png", paste(genes_of_interest, collapse = "_")),
          dir = paths$plots_deg, width = 9, height = 5)


# =========================================================
#                  2) Genes in the SV locus
# =========================================================

# Same table for all genes in the locus
locus_table <- deg_ranked |>
  dplyr::filter(SYMBOL %in% cfg$advanced$locus_genes) |>
  dplyr::arrange(padj)

# Locus genes missing from the DEG table (filtered out or unannotated)
locus_missing <- setdiff(cfg$advanced$locus_genes, locus_table$SYMBOL)
if (length(locus_missing) > 0) {
  message("[LOCUS] Not in the DEG table: ", paste(locus_missing, collapse = ", "))
}

print(locus_table)

save_csv(locus_table, filename = sprintf("locus_genes_%s.csv", cfg$gse_id),
         dir = paths$tables)


# =========================================================
#                  3) Ridge plot (GO-BP GSEA)
# =========================================================

ridge_plot <- gsea_ridgeplot(gobp_simp, lfc_table = deg_table)

save_plot_safe(ridge_plot, filename = "ridgeplot_GO-BP_GSEA.png",
               dir = paths$plots_gsea, width = 12, height = 9)


# =========================================================
#                  4) Paper GO terms (e.g. myeloid signature)
# =========================================================

# NES, padj and size of the GO-BP terms reported in the paper,
# from the complete table (before redundancy removal)
paper_terms <- gobp_full |>
  dplyr::filter(grepl(cfg$advanced$paper_terms, pathway, ignore.case = TRUE)) |>
  dplyr::select(pathway, NES, padj, set_size) |>
  dplyr::arrange(padj)

print(paper_terms)

save_csv(paper_terms, filename = sprintf("paper_terms_GO_BP_%s.csv", cfg$gse_id),
         dir = paths$tables_gsea)

save_session_info(filename = sprintf("sessionInfo_advanced_%s.txt", cfg$gse_id), dir = paths$logs,
                  label = paste0("Advanced analysis - ", cfg$gse_id))

message("=================================================")
message("[PIPELINE] Advanced analysis completed for ", cfg$gse_id,
        " (", cfg$case_level, " vs ", cfg$ref_level, ").")
message("=================================================")
