#!/usr/bin/env Rscript
# Run from the package root: Rscript scripts/run_all.R full
# quick = validation without model bootstrap; full = 1000 model / 5000 trend resamples.
# Saved reference results are not overwritten.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript scripts/run_all.R [full|quick]")
mode <- if (length(args)) args[1] else "full"
if (!mode %in% c("full", "quick")) stop("Mode must be full or quick")
if (!requireNamespace("ordinal", quietly = TRUE))
  stop("Install ordinal first: install.packages('ordinal')")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this driver using Rscript")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- dirname(dirname(script))
output <- file.path(root, "recomputed", mode)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
input <- file.path(root, "data", "phase4_current_results_p1204.csv")
folds <- file.path(root, "data", "folds.csv")
B <- if (mode == "full") 1000L else 0L
trend_B <- if (mode == "full") 5000L else 100L
rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
run <- function(name, ...) {
  a <- c(file.path(root, "scripts", name), as.character(c(...)))
  cat("\nRunning ", name, "\n", sep = "")
  status <- system2(rscript, args = vapply(a, shQuote, character(1)))
  if (status != 0L) stop("Step failed: ", name, " (exit ", status, ")")
}
main <- file.path(output, "main")
clean <- file.path(output, "screened_out")
run("model_revision.R", input, main, B, "all", folds)
run("model_revision.R", input, clean, 0, "screened_out", folds)
run("screening_compare.R", file.path(main, "oof_predictions.csv"),
    file.path(clean, "oof_predictions.csv"), file.path(output, "screening"))
run("table3_paired.R", input, file.path(output, "paired_tests"))
run("transition_audit.R", input, file.path(main, "oof_predictions.csv"),
    file.path(output, "transitions"))
run("longitudinal_missingness.R", input, file.path(output, "longitudinal"), trend_B)
run("plot_cv_weekly_mos.R", file.path(main, "cv_week_mos.csv"),
    file.path(output, "figures"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
writeLines(c(paste("mode:", mode), paste("model_bootstrap:", B),
             paste("trend_bootstrap:", trend_B), paste("fixed_folds:", folds)),
           file.path(output, "run_settings.txt"))
cat("\nCompleted. Recomputed outputs: ", output, "\n", sep = "")
