# =========================================================
# Helper functions: annotation, gene sets, enrichment results
# =========================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------
# Reading analysis tables (ENTREZID always as character)
# ---------------------------------------------------------
read_deg_table <- function(path) {
  stopifnot(file.exists(path))
  utils::read.csv(path, colClasses = c(ENTREZID = "character"),
                  stringsAsFactors = FALSE)
}

read_universe <- function(path) {
  stopifnot(file.exists(path))
  x <- utils::read.csv(path, colClasses = c(ENTREZID = "character"),
                       stringsAsFactors = FALSE)
  unique(x$ENTREZID[!is.na(x$ENTREZID)])
}

# ---------------------------------------------------------
# Gene identifiers
# ---------------------------------------------------------

# Entrez -> gene symbol (falls back to the Entrez ID if unmapped)
entrez_to_symbol <- function(ids) {
  ids <- as.character(ids)
  map <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids),
                          keytype = "ENTREZID", column = "SYMBOL",
                          multiVals = "first")
  )
  out <- unname(map[ids])
  ifelse(is.na(out), ids, out)
}

# "1/2/3" (Entrez, "/"-separated) -> "SYM1/SYM2/SYM3"
entrez_string_to_symbol <- function(x, sep = "/") {
  ids <- strsplit(as.character(x), sep, fixed = TRUE)
  ids <- lapply(ids, function(v) v[!is.na(v) & nzchar(v)])
  all_ids <- unique(unlist(ids))
  if (length(all_ids) == 0) return(rep("", length(ids)))
  map <- stats::setNames(entrez_to_symbol(all_ids), all_ids)
  vapply(ids, function(v) paste(map[v], collapse = sep), character(1))
}

# Version-independent replacement for setReadable() (it moved between
# DOSE / clusterProfiler across Bioconductor releases).
# Works on enrichResult (geneID) and gseaResult (core_enrichment).
make_readable <- function(x) {
  if (is.null(x) || nrow(x@result) == 0) return(x)
  col <- if (inherits(x, "gseaResult")) "core_enrichment" else "geneID"
  x@result[[col]] <- entrez_string_to_symbol(x@result[[col]])
  x@readable <- TRUE
  x
}

# ---------------------------------------------------------
# Gene sets
# ---------------------------------------------------------

# MSigDB (msigdbr >= 10) -> named list of Entrez character vectors
get_msigdb_sets <- function(collection, subcollection = NULL) {
  args <- list(species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection
  df <- do.call(msigdbr::msigdbr, args)
  df <- df[!is.na(df$ncbi_gene), c("gs_name", "ncbi_gene")]
  lapply(split(as.character(df$ncbi_gene), df$gs_name), unique)
}

# "HALLMARK_TNFA_SIGNALING_VIA_NFKB" -> "Tnfa Signaling via Nfkb"
# Use ONLY for Hallmark: it would corrupt TF symbols / GSE identifiers.
clean_hallmark <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# ---------------------------------------------------------
# Backward-compatibility shim (delete once the scripts use paths$tables_ora)
# ---------------------------------------------------------
tables_ora_dir <- function() paths$tables_ora
