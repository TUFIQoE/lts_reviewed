#!/usr/bin/env Rscript
# Compare a full rerun with the supplied reference tables.
# Usage: Rscript scripts/check_results.R recomputed/full [ABSOLUTE_TOLERANCE]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L)
  stop("Usage: Rscript scripts/check_results.R OUTPUT_DIR [ABSOLUTE_TOLERANCE]")
tol <- if (length(args) == 2L) as.numeric(args[2]) else 1e-6
if (!is.finite(tol) || tol < 0) stop("Tolerance must be finite and nonnegative")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run using Rscript")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
output <- normalizePath(args[1], mustWork = TRUE)
specs <- list(
  list("main/cv_metrics.csv", c("dataset", "model"), c("N", "RMSE", "MAE", "PLCC", "SRCC")),
  list("main/cv_week_metrics.csv", c("dataset", "model"), c("N", "RMSE", "MAE", "PLCC", "SRCC")),
  list("main/cv_week_mos.csv", c("dataset", "model", "week"), c("n_votes", "actual_mos", "predicted_mos")),
  list("main/full_sample_fits.csv", c("dataset", "model"), c("N", "AIC", "a", "b")),
  list("main/paired_metric_differences.csv", c("dataset", "comparator", "metric"), c("estimate", "lower", "upper", "B")),
  list("main/parameter_weights_ci_all_weeks.csv", c("dataset", "parameter"), c("estimate", "lower", "upper", "successful", "attempted")),
  list("main/parameter_weights_ci_complete.csv", c("dataset", "parameter"), c("estimate", "lower", "upper", "successful", "attempted")),
  list("screened_out/cv_metrics.csv", c("dataset", "model"), c("N", "RMSE", "MAE", "PLCC", "SRCC")),
  list("screening/screening_matched_cv.csv", c("dataset", "model"),
       c("N", "n_participants", "train_all_RMSE", "train_clean_RMSE", "delta_RMSE",
         "train_all_MAE", "train_clean_MAE", "delta_MAE", "train_all_PLCC", "train_clean_PLCC",
         "train_all_SRCC", "train_clean_SRCC")),
  list("paired_tests/table3_paired_results.csv", "row",
       c("n_participants", "mean_difference", "t", "df", "p", "ci_lower", "ci_upper", "cohen_dz", "p_holm_12")),
  list("transitions/week_transition_and_cv_errors.csv", "week",
       c("big_between_pairs", "big_within_pairs", "odd_squared_error", "even_squared_error", "n_votes")),
  list("longitudinal/exploratory_within_participant_trends.csv", c("model", "predictor"),
       c("N", "n_participants", "estimate", "lower", "upper", "successful", "attempted"))
)
compare_table <- function(spec) {
  path <- spec[[1]]; keys <- spec[[2]]; values <- spec[[3]]
  ref <- read.csv(file.path(root, "results", path), stringsAsFactors = FALSE)
  got <- read.csv(file.path(output, path), stringsAsFactors = FALSE)
  needed <- c(keys, values)
  if (!all(needed %in% names(ref)) || !all(needed %in% names(got)))
    stop("Missing columns in ", path)
  key <- function(x) do.call(paste, c(x[keys], sep = "|"))
  kr <- key(ref); kg <- key(got)
  if (anyDuplicated(kr) || anyDuplicated(kg) || !setequal(kr, kg))
    stop("Missing, extra, or duplicated row keys in ", path)
  got <- got[match(kr, kg), , drop = FALSE]
  for (col in values) {
    a <- ref[[col]]; b <- got[[col]]
    if (!is.numeric(a) || !is.numeric(b)) stop("Nonnumeric column: ", path, "/", col)
    if (any(is.na(a) != is.na(b))) stop("NA pattern differs in ", path, "/", col)
    valid <- !is.na(a)
    if (any(!is.finite(a[valid])) || any(!is.finite(b[valid]))) stop("Nonfinite values in ", path)
    delta <- if (any(valid)) max(abs(a[valid] - b[valid])) else 0
    if (delta > tol) stop("Difference in ", path, "/", col,
                          ": maximum absolute difference = ", format(delta, digits = 8))
  }
  cat("PASS: ", path, "\n", sep = "")
}
failures <- character()
for (spec in specs) tryCatch(compare_table(spec), error = function(e) {
  failures <<- c(failures, conditionMessage(e))
  cat("FAIL: ", conditionMessage(e), "\n", sep = "")
})
if (length(failures)) stop(length(failures), " comparison(s) failed. Inspect the reported tables and sessionInfo.txt.")
cat("All reference comparisons passed at absolute tolerance ", tol, ".\n", sep = "")
