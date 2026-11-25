# =========================================================
# Project paths (single source of truth)
# =========================================================

project_root <- Sys.getenv("/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects/Bulk_RNASeq_GSE166409", unset = NA_character_)

if (is.na(project_root) || !nzchar(project_root)) {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required (or set the PROJECT_ROOT variable).", call. = FALSE)
  }
  project_root <- here::here()
}

project_root <- normalizePath(project_root, mustWork = TRUE)
results_dir  <- file.path(project_root, "results")

paths <- list(

  # project
  base  = project_root,
  data  = file.path(project_root, "data"),
  setup = file.path(project_root, "SetupEnvironment"),

  # results root
  results = results_dir,

  # plots
  plots_qc   = file.path(results_dir, "plots/QC"),
  plots_eda  = file.path(results_dir, "plots/EDA"),
  plots_deg  = file.path(results_dir, "plots/DEG"),
  plots_ora  = file.path(results_dir, "plots/ORA"),
  plots_gsea = file.path(results_dir, "plots/GSEA"),

  # tables
  tables      = file.path(results_dir, "tables"),
  tables_ora  = file.path(results_dir, "tables", "ORA_results"),
  tables_gsea = file.path(results_dir, "tables", "GSEA_results"),

  # logs
  logs = file.path(results_dir, "logs")
)

# Create output directories (never the raw data / setup folders)
output_dirs <- paths[setdiff(names(paths), c("base", "data", "setup"))]
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!dir.exists(paths$data)) {
  warning("Input data folder not found: ", paths$data, call. = FALSE)
}
