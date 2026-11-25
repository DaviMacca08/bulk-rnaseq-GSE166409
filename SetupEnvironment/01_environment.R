# =========================================================
# Environment & packages
# =========================================================
# Principle: attach only what the scripts use unqualified; everything
# else is called as pkg::function() and only checked here. Fewer attached
# packages = fewer masked functions and a faster start.
# Record exact versions with renv::snapshot() after a successful run.

pkgs_attach <- c(
  "tidyverse",        # dplyr, tidyr, ggplot2, tibble, stringr, readr, purrr
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggridges"
)

pkgs_namespace <- c(
  "here", "matrixStats", "SummarizedExperiment", "AnnotationDbi",
  "apeglm", "GEOquery",
  "ggrepel", "pheatmap", "ComplexHeatmap", "scales",
  "DOSE", "ReactomePA", "msigdbr", "fgsea"
)

# Fail early, with one clear message, if anything is missing
pkgs_missing <- setdiff(c(pkgs_attach, pkgs_namespace),
                        rownames(utils::installed.packages()))

if (length(pkgs_missing) > 0) {
  stop("Missing packages: ", paste(pkgs_missing, collapse = ", "),
       "\nInstall them (BiocManager::install() for Bioconductor ones) and retry.",
       call. = FALSE)
}

suppressPackageStartupMessages(
  invisible(lapply(pkgs_attach, library, character.only = TRUE))
)

# ---------------------------------------------------------
# Resolve function masking explicitly
# ---------------------------------------------------------
# Bioconductor packages attached after tidyverse mask dplyr verbs
# (AnnotationDbi::select, S4Vectors::rename, IRanges::slice/desc, ...).
# Pin the dplyr versions so the result does not depend on load order.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice  <- dplyr::slice
desc   <- dplyr::desc
first  <- dplyr::first

options(warn = 1)   # print warnings as they happen
