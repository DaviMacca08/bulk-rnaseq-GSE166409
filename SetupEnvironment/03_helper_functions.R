# =========================================================
# Helper functions: file output (plots, tables, objects, logs)
# =========================================================

# Internal: stop early, with a clear message, if the folder is missing
.assert_dir <- function(dir) {
  if (!dir.exists(dir)) stop("Directory does not exist: ", dir, call. = FALSE)
  invisible(dir)
}

# ---------------------------------------------------------
# ggplot / patchwork objects
# ---------------------------------------------------------
save_plot <- function(plot, filename, dir,
                      width = 12, height = 8, dpi = 300, bg = "white") {

  if (is.null(plot)) {
    message("No plot to save: ", filename)
    return(invisible(NULL))
  }

  .assert_dir(dir)
  path <- file.path(dir, filename)

  ggplot2::ggsave(filename = path, plot = plot,
                  width = width, height = height, dpi = dpi, bg = bg)

  message("Saved plot: ", filename)
  invisible(path)
}

# Same, but silently skips NULL (e.g. nothing significant to plot)
save_plot_safe <- function(p, filename, dir, ...) {
  if (is.null(p)) {
    message("[SAVE] Skipped (nothing significant to plot): ", filename)
    return(invisible(NULL))
  }
  save_plot(p, filename = filename, dir = dir, ...)
}

# ---------------------------------------------------------
# Base R graphics (pheatmap, ComplexHeatmap::draw, plot, ...)
# ---------------------------------------------------------
# Preferred: the device is ALWAYS closed, even if the code errors.
#   save_base_plot("x.png", paths$plots_qc, {
#     plot(1:10)
#   })
# Size is in inches (same convention as save_plot).
save_base_plot <- function(filename, dir, expr,
                           width = 10, height = 8, dpi = 300) {

  .assert_dir(dir)
  path <- file.path(dir, filename)

  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  force(expr)   # lazily evaluated here, i.e. with the device already open

  message("Saved plot: ", filename)
  invisible(path)
}

# Legacy pair, kept so existing scripts keep working. Size is in PIXELS.
# Prefer save_base_plot() in new code.
open_png <- function(filename, dir, width = 1200, height = 1000, res = 120) {
  .assert_dir(dir)
  grDevices::png(file.path(dir, filename), width = width, height = height, res = res)
  message("Saving plot: ", filename)
  invisible(NULL)
}

close_png <- function() {
  invisible(grDevices::dev.off())
}

# ---------------------------------------------------------
# R objects
# ---------------------------------------------------------
save_rds <- function(object, filename, dir) {

  .assert_dir(dir)
  path <- file.path(dir, filename)

  saveRDS(object, path)

  message("Saved object: ", filename)
  invisible(path)
}

# ---------------------------------------------------------
# Tables (CSV)
# ---------------------------------------------------------
save_csv <- function(object, filename, dir) {

  .assert_dir(dir)
  object <- as.data.frame(object)

  list_cols <- names(object)[vapply(object, is.list, logical(1))]
  if (length(list_cols) > 0) {
    stop("Cannot write list-columns to CSV: ", paste(list_cols, collapse = ", "),
         ". Collapse them to strings first.", call. = FALSE)
  }

  path <- file.path(dir, filename)
  utils::write.csv(object, file = path, row.names = FALSE)

  message("Saved table: ", filename)
  invisible(path)
}

# Same, but skips NULL / empty tables
save_csv_safe <- function(x, filename, dir) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    message("[SAVE] Skipped (empty table): ", filename)
    return(invisible(NULL))
  }
  save_csv(x, filename = filename, dir = dir)
}

# ---------------------------------------------------------
# Session information
# ---------------------------------------------------------
save_session_info <- function(filename = "sessionInfo.txt", dir, label = NULL) {

  .assert_dir(dir)
  path <- file.path(dir, filename)

  header <- if (is.null(label)) "Session information"
            else paste0("Session information (", label, ")")

  writeLines(
    c(header, as.character(Sys.time()), "", utils::capture.output(utils::sessionInfo())),
    path
  )

  message("Session information saved at: ", path)
  invisible(path)
}
