# =========================================================
# Project      : Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE166409 (Illumina NextSeq 500)
# Samples      : 3 patients + 3 healthy controls, PBMC
# Condition    : WAC-ANKRD26 gene fusion-associated thrombocytopenia
# Organism     : Homo sapiens (GRCh38.p13, NCBI raw counts)
# Script       : 01 - Metadata, QC, EDA and Differential Expression (DESeq2)
# Input        : GEO metadata + raw count table (data/)
# Outputs      : DEGs_<GSE>.csv (all genes, raw + shrunken log2FC),
#                DEGs_significant_<GSE>.csv, universe_*.csv, QC table,
#                normalized counts, DESeq2 object, QC / EDA / DEG plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

# Packages, paths, cfg and helper functions (save_plot, save_csv, ...)
source(here::here("SetupEnvironment", "setup.R"))

set_seed(cfg$seed)

message("=== Starting ", cfg$gse_id, " RNA-seq analysis ===")


# =========================================================
#                  Functions
# =========================================================

# Robust z-score based on median and MAD (NA if MAD is zero)
robust_z <- function(x) {
  s <- stats::mad(x)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - stats::median(x)) / s
}


# =========================================================
#                  Metadata
# =========================================================

# Download the GEO series and keep the sample table
gse_list <- GEOquery::getGEO(cfg$gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
gse      <- if (is.list(gse_list)) gse_list[[1]] else gse_list
meta     <- Biobase::pData(gse)

# Collapse all 'characteristics_ch1*' columns into one string per sample
char_cols <- grep("^characteristics_ch1", colnames(meta), value = TRUE)
stopifnot(length(char_cols) > 0)

char_text <- apply(meta[, char_cols, drop = FALSE], 1, paste, collapse = " | ")

# Assign each sample to the case or reference group by keyword
group_chr <- dplyr::case_when(
  grepl(cfg$case_pattern, char_text, ignore.case = TRUE) ~ cfg$case_level,
  grepl(cfg$ref_pattern,  char_text, ignore.case = TRUE) ~ cfg$ref_level,
  TRUE ~ NA_character_
)

# Print the mapping to check it and to look for confounders (sex, age, batch)
message("[META] Raw characteristics -> assigned group:")
print(data.frame(sample = meta$geo_accession, group = group_chr,
                 characteristics = char_text, row.names = NULL))

if (anyNA(group_chr)) {
  stop("Some samples could not be classified into ",
       cfg$case_level, "/", cfg$ref_level, ". Check 'char_text'.")
}

# Sample table for DESeq2 (reference level first)
coldata <- data.frame(
  row.names = meta$geo_accession,
  condition = factor(group_chr, levels = c(cfg$ref_level, cfg$case_level))
)

message("[META] Sample distribution:")
print(table(coldata$condition))
stopifnot(all(table(coldata$condition) >= 2))


# =========================================================
#                  Count matrix
# =========================================================

counts_raw <- utils::read.delim(
  cfg$counts_file,
  header = TRUE, sep = "\t", row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE
)

counts <- as.matrix(counts_raw)
rm(counts_raw)

# Order the count columns as in the sample table (matched by name)
stopifnot(setequal(colnames(counts), rownames(coldata)))
counts <- counts[, rownames(coldata)]

# Basic checks on the count values
stopifnot(
  nrow(counts) > 1000,
  !anyDuplicated(rownames(counts)),
  !anyNA(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts))
)
storage.mode(counts) <- "integer"

rownames(counts) <- as.character(rownames(counts))   # Entrez IDs as text


# =========================================================
#                  DESeq2
# =========================================================

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition            # reference level = cfg$ref_level
)

# Expression filter: >= min_count reads in >= (smallest group size) samples
min_samples <- min(table(coldata$condition))
keep <- rowSums(DESeq2::counts(dds) >= cfg$min_count) >= min_samples

message(sprintf("[FILTER] Kept %d / %d genes (>= %d reads in >= %d samples).",
                sum(keep), nrow(dds), cfg$min_count, min_samples))
dds <- dds[keep, ]

# Remove genes on chromosomes X and Y (cfg$exclude_sex_chr)
if (cfg$exclude_sex_chr) {
  chr <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = rownames(dds), keytype = "ENTREZID",
    column = "CHR", multiVals = "first"
  ))
  is_sex_chr <- chr %in% c("X", "Y")

  message(sprintf("[SEX] %d X/Y-linked genes removed before DESeq2.", sum(is_sex_chr)))
  dds <- dds[!is_sex_chr, ]
}

# Fit the model
dds <- DESeq2::DESeq(dds)

coef_name <- paste0("condition_", cfg$case_level, "_vs_", cfg$ref_level)
stopifnot(coef_name %in% DESeq2::resultsNames(dds))

# Unshrunken results: Wald stat, p-values and FDR (used for the GSEA ranking)
res_raw <- DESeq2::results(
  dds, contrast = c("condition", cfg$case_level, cfg$ref_level),
  alpha = cfg$padj_cutoff
)

# Shrunken log2FC (apeglm): used for thresholds, volcano, MA plot and heatmap
stopifnot(requireNamespace("apeglm", quietly = TRUE))
res_shr <- DESeq2::lfcShrink(dds, coef = coef_name,
                             type = cfg$lfc_shrink_type, quiet = TRUE)

stopifnot(identical(rownames(res_raw), rownames(res_shr)))

message("[DESeq2] Genes with padj < ", cfg$padj_cutoff, " (no log2FC filter): ",
        sum(res_raw$padj < cfg$padj_cutoff, na.rm = TRUE))


# =========================================================
#                  Results table + annotation
# =========================================================

# Merge unshrunken and shrunken results in one table
shr_df <- as.data.frame(res_shr)

deg_table <- as.data.frame(res_raw) |>
  tibble::rownames_to_column("ENTREZID") |>
  dplyr::rename(log2FC_raw = log2FoldChange, lfcSE_raw = lfcSE)

deg_table$log2FC_shrunk <- shr_df[deg_table$ENTREZID, "log2FoldChange"]
deg_table$lfcSE_shrunk  <- shr_df[deg_table$ENTREZID, "lfcSE"]

# Gene symbol and name from the Entrez ID
deg_table$SYMBOL <- unname(suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = deg_table$ENTREZID, keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first")))

deg_table$GENENAME <- unname(suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = deg_table$ENTREZID, keytype = "ENTREZID",
  column = "GENENAME", multiVals = "first")))

# Up / down / not-significant call (padj + shrunken log2FC thresholds)
deg_table <- deg_table |>
  dplyr::mutate(
    significance = dplyr::case_when(
      !is.na(padj) & padj < cfg$padj_cutoff & log2FC_shrunk >  cfg$lfc_cutoff ~ cfg$label_up,
      !is.na(padj) & padj < cfg$padj_cutoff & log2FC_shrunk < -cfg$lfc_cutoff ~ cfg$label_down,
      TRUE ~ cfg$label_ns
    )
  ) |>
  dplyr::arrange(padj)

sig_deg <- dplyr::filter(deg_table, significance != cfg$label_ns)

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |shrunken log2FC| > %.1f): %d up, %d down (%d total); %d lack a gene symbol.",
  cfg$padj_cutoff, cfg$lfc_cutoff,
  sum(sig_deg$significance == cfg$label_up),
  sum(sig_deg$significance == cfg$label_down),
  nrow(sig_deg), sum(is.na(sig_deg$SYMBOL))
))

# Gene universes for the enrichment scripts:
#  - expressed: genes that passed the expression filter
#  - tested:    genes with a non-NA adjusted p-value (used by ORA)
universe_expressed <- data.frame(ENTREZID = rownames(dds))
universe_tested    <- data.frame(ENTREZID = deg_table$ENTREZID[!is.na(deg_table$padj)])


# =========================================================
#                  Transformations
# =========================================================

# VST without design (blind = TRUE): QC and EDA
# VST with design (blind = FALSE): DEG heatmap
vsd_qc <- DESeq2::vst(dds, blind = TRUE)
vsd_de <- DESeq2::vst(dds, blind = FALSE)

mat_qc <- SummarizedExperiment::assay(vsd_qc)
mat_de <- SummarizedExperiment::assay(vsd_de)

sample_colors <- group_palette[as.character(coldata$condition)]


# =========================================================
#                  Quality Control
# =========================================================

# --- Per-sample QC table (raw counts)

# Mean correlation of each sample with the others (blind VST)
cor_mat  <- stats::cor(mat_qc, method = "pearson")
mean_cor <- (rowSums(cor_mat) - 1) / (ncol(cor_mat) - 1)

qc_metrics <- data.frame(
  sample         = colnames(counts),
  group          = coldata$condition,
  lib_size       = colSums(counts),
  detected_genes = colSums(counts > 0),
  median_vst     = matrixStats::colMedians(mat_qc),
  mean_cor       = unname(mean_cor[colnames(counts)]),
  row.names      = NULL
)

# Outlier flags (|robust z| > 3.5) only if there are enough samples
if (nrow(qc_metrics) >= cfg$min_n_outlier_rules) {
  qc_metrics$z_lib_size <- robust_z(log10(qc_metrics$lib_size))
  qc_metrics$z_detected <- robust_z(qc_metrics$detected_genes)
  qc_metrics$z_mean_cor <- robust_z(qc_metrics$mean_cor)

  z <- cbind(qc_metrics$z_lib_size, qc_metrics$z_detected, qc_metrics$z_mean_cor)
  qc_metrics$flag_outlier <- rowSums(abs(z) > 3.5, na.rm = TRUE) > 0
} else {
  message("[QC] n = ", nrow(qc_metrics), " < ", cfg$min_n_outlier_rules,
          ": automatic outlier flags skipped. No sample removed.")
  qc_metrics$flag_outlier <- NA
}

# --- VST distributions

# One boxplot per sample
save_base_plot("boxplot_vst_distributions.png", paths$plots_qc, {
  par(mar = c(8, 4, 4, 2))
  boxplot(mat_qc, col = sample_colors, las = 2, cex.axis = 0.7, outline = FALSE,
          main = paste0("VST-transformed expression per sample (", cfg$gse_id, ")"),
          ylab = "VST expression")
  legend("topright", legend = names(group_palette), fill = group_palette, bty = "n")
}, width = 13, height = 7.5)

# One density curve per sample
dens <- apply(mat_qc, 2, stats::density)
ymax <- max(vapply(dens, function(d) max(d$y), numeric(1)))
xr   <- range(vapply(dens, function(d) range(d$x), numeric(2)))

save_base_plot("density_vst_expression.png", paths$plots_qc, {
  plot(NULL, xlim = xr, ylim = c(0, ymax * 1.05),
       xlab = "VST expression", ylab = "Density",
       main = paste0("Density of VST-transformed expression (", cfg$gse_id, ")"))
  for (i in seq_along(dens)) lines(dens[[i]], col = sample_colors[i], lwd = 1.5)
  legend("topright", legend = names(group_palette), col = group_palette,
         lwd = 2, bty = "n")
}, width = 10, height = 6.5)

# --- Library size per sample

libsize_plot <- ggplot(qc_metrics,
                       aes(x = reorder(sample, lib_size), y = lib_size, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = group_palette) +
  labs(title = paste0("Library size per sample (", cfg$gse_id, ")"),
       x = NULL, y = "Total counts assigned to genes", fill = NULL) +
  theme_bw()

save_plot(libsize_plot, filename = "library_size.png", dir = paths$plots_qc,
          width = 9, height = 6)

# --- DESeq2 diagnostics

# Dispersion estimates and fitted trend
save_base_plot("dispersion_estimates.png", paths$plots_qc, {
  DESeq2::plotDispEsts(dds, main = "DESeq2 dispersion estimates")
}, width = 8, height = 6.5)

# Raw p-value histogram
pval_hist <- ggplot(deg_table[!is.na(deg_table$pvalue), ], aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey70", colour = "white") +
  labs(title = "Raw p-value distribution", x = "p-value", y = "Genes") +
  theme_bw()

save_plot(pval_hist, filename = "pvalue_histogram.png", dir = paths$plots_qc,
          width = 7, height = 5)


# =========================================================
#           Exploratory Data Analysis (EDA)
# =========================================================

# --- PCA on the 500 most variable genes (blind VST)

top_var <- order(matrixStats::rowVars(mat_qc), decreasing = TRUE)
top_var <- top_var[seq_len(min(500, length(top_var)))]

pca_res <- stats::prcomp(t(mat_qc[top_var, ]), center = TRUE, scale. = FALSE)
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pca_var <- pca_var[pca_var > 1e-8]          # drop PCs with ~zero variance

pca_df <- data.frame(
  sample = colnames(mat_qc),
  PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
  group = coldata$condition
)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  xlab(sprintf("PC1 (%.1f%%)", pca_var[1])) +
  ylab(sprintf("PC2 (%.1f%%)", pca_var[2])) +
  scale_colour_manual(values = group_palette) +
  labs(title = paste0("PCA - top 500 variable genes (", cfg$gse_id, ")"),
       colour = NULL) +
  theme_bw()

save_plot(pca_plot, filename = "pca_plot.png", dir = paths$plots_eda)

# --- Variance explained per PC (bars) and cumulative variance (line)

scree_df <- data.frame(PC = factor(paste0("PC", seq_along(pca_var)),
                                   levels = paste0("PC", seq_along(pca_var))),
                       variance = pca_var, cumulative = cumsum(pca_var))

scree_plot <- ggplot(scree_df, aes(PC)) +
  geom_col(aes(y = variance), fill = "grey70") +
  geom_line(aes(y = cumulative, group = 1), colour = "red") +
  geom_point(aes(y = cumulative), colour = "red") +
  labs(title = "PCA scree plot", x = NULL,
       y = "% variance (bars) / cumulative % (line)") +
  theme_bw()

save_plot(scree_plot, filename = "pca_scree_plot.png", dir = paths$plots_eda,
          width = 8, height = 5)

# --- Hierarchical clustering (Euclidean distance, average linkage)

hc <- stats::hclust(stats::dist(t(mat_qc), method = "euclidean"), method = "average")

save_base_plot("hierarchical_clustering.png", paths$plots_eda, {
  par(mar = c(5, 4, 4, 2))
  plot(hc, labels = paste0(colnames(mat_qc), " (", coldata$condition, ")"),
       main = "Hierarchical clustering (Euclidean, average linkage, blind VST)",
       xlab = "", sub = "")
}, width = 11, height = 7.5)

# --- Sample correlation heatmap

annotation_col <- data.frame(Group = coldata$condition,
                             row.names = colnames(mat_qc))

save_base_plot("sample_correlation_heatmap.png", paths$plots_eda, {
  pheatmap::pheatmap(
    cor_mat,
    annotation_col    = annotation_col,
    annotation_row    = annotation_col,
    annotation_colors = list(Group = group_palette),
    show_rownames = TRUE, show_colnames = TRUE,
    main = "Sample-sample Pearson correlation (blind VST)"
  )
}, width = 11, height = 9.5)


# =========================================================
#           Plotting DEGs
# =========================================================

# --- Volcano plot (shrunken log2FC)

volcano_df <- dplyr::filter(deg_table, !is.na(padj))

# Top 15 significant genes by padj, labelled on the plot
label_df <- volcano_df |>
  dplyr::filter(significance != cfg$label_ns, !is.na(SYMBOL)) |>
  dplyr::arrange(padj, dplyr::desc(abs(log2FC_shrunk))) |>
  dplyr::slice_head(n = 15)

volcano_colors <- stats::setNames(
  c("#CE2915", "#0096FF", "grey75"),
  c(cfg$label_up, cfg$label_down, cfg$label_ns)
)

volcano_plot <- ggplot(volcano_df,
                       aes(x = log2FC_shrunk, y = -log10(padj), colour = significance)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = c(-cfg$lfc_cutoff, cfg$lfc_cutoff),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(cfg$padj_cutoff),
             linetype = "dashed", colour = "grey40") +
  ggrepel::geom_text_repel(data = label_df, aes(label = SYMBOL),
                           size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = volcano_colors) +
  labs(title = paste0("Differential expression ", cfg$case_level, " vs ",
                      cfg$ref_level, " (", cfg$gse_id, ")"),
       x = expression(shrunken ~ log[2] ~ fold ~ change),
       y = expression(-log[10] ~ adjusted ~ italic(P)),
       colour = NULL) +
  theme_bw()

save_plot(volcano_plot, filename = "volcano_plot.png", dir = paths$plots_deg,
          width = 8, height = 6)

# --- MA plot (shrunken log2FC)

save_base_plot("MA_plot.png", paths$plots_deg, {
  DESeq2::plotMA(res_shr, ylim = c(-4, 4),
                 main = paste0("MA plot (shrunken LFC): ", cfg$case_level, " vs ", cfg$ref_level))
  abline(h = c(-cfg$lfc_cutoff, cfg$lfc_cutoff), col = "grey50", lty = 2)
}, width = 8, height = 6.5)

# --- Heatmap of the top DEGs (VST with design, row Z-score)

# Top significant annotated genes, one row per gene symbol
top_tbl <- sig_deg |>
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
  dplyr::arrange(padj, dplyr::desc(abs(log2FC_shrunk))) |>
  dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
  dplyr::slice_head(n = cfg$top_n_heatmap)

message(sprintf("[DEG] %d annotated significant genes used for the heatmap.", nrow(top_tbl)))

if (nrow(top_tbl) >= 2) {

  # Row Z-score, dropping genes with zero variance; rows named by symbol
  hm <- mat_de[top_tbl$ENTREZID, , drop = FALSE]
  hm <- t(scale(t(hm)))
  hm <- hm[stats::complete.cases(hm), , drop = FALSE]
  rownames(hm) <- top_tbl$SYMBOL[match(rownames(hm), top_tbl$ENTREZID)]

  ann_col <- ComplexHeatmap::HeatmapAnnotation(
    Group = coldata$condition, col = list(Group = group_palette)
  )

  ht <- ComplexHeatmap::Heatmap(
    hm, name = "Row Z-score", top_annotation = ann_col,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_title = sprintf("Top %d DEGs (%s vs %s)", nrow(hm), cfg$case_level, cfg$ref_level),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
  )

  save_base_plot("Heatmap_TopDEGs.png", paths$plots_deg, {
    ComplexHeatmap::draw(ht)
  }, width = 11, height = 15)

} else {
  message("[DEG] Fewer than 2 annotated DEGs: heatmap skipped.")
}


# =========================================================
#                  Save outputs
# =========================================================

message("[OUTPUT] Saving tables and objects...")

save_rds(dds, filename = "DESeq2_object.rds", dir = paths$results)

save_csv(meta, filename = sprintf("metadata_%s.csv", cfg$gse_id), dir = paths$tables)

save_csv(deg_table, filename = sprintf("DEGs_%s.csv", cfg$gse_id), dir = paths$tables)
save_csv(sig_deg,   filename = sprintf("DEGs_significant_%s.csv", cfg$gse_id), dir = paths$tables)

save_csv(universe_expressed, filename = sprintf("universe_expressed_%s.csv", cfg$gse_id), dir = paths$tables)
save_csv(universe_tested,    filename = sprintf("universe_tested_%s.csv", cfg$gse_id),    dir = paths$tables)

save_csv(qc_metrics, filename = sprintf("QC_sample_metrics_%s.csv", cfg$gse_id), dir = paths$tables)

norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE)) |>
  tibble::rownames_to_column("ENTREZID")
save_csv(norm_counts, filename = sprintf("normalized_counts_%s.csv", cfg$gse_id), dir = paths$tables)

save_session_info(filename = sprintf("sessionInfo_DEG_%s.txt", cfg$gse_id), dir = paths$logs,
                  label = paste0("QC, EDA and differential expression - ", cfg$gse_id))

message("=================================================")
message("[PIPELINE] Differential expression completed for ", cfg$gse_id,
        " (", cfg$case_level, " vs ", cfg$ref_level, ").")
message("=================================================")
