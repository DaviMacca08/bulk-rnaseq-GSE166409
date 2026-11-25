# =========================================================
# Single entry point for every analysis script:
#   source(here::here("SetupEnvironment", "setup.R"))
# =========================================================

local({
  root <- Sys.getenv("/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects/Bulk_RNASeq_Susanna", unset = NA_character_)
  if (is.na(root) || !nzchar(root)) root <- here::here()

  setup_files <- c(
    "00_paths.R",
    "01_environment.R",
    "02_seed.R",
    "03_helper_functions.R",
    "04_config.R",
    "05_helpers_enrichment.R"
  )

  for (f in setup_files) {
    source(file.path(root, "SetupEnvironment", f), local = FALSE)
  }
})
