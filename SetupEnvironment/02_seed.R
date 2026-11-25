# =========================================================
# Reproducibility
# =========================================================

#' Set the random seed with explicit RNG settings, so that results do not
#' depend on the R version's defaults.
set_seed <- function(seed = 1234) {

  set.seed(seed,
           kind        = "Mersenne-Twister",
           normal.kind = "Inversion",
           sample.kind = "Rejection")

  message("Seed set to: ", seed)
  invisible(seed)
}
