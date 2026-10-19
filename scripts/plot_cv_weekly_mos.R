# Weekly MOS and out-of-fold predictions for manuscript Figure 8.
# Base R only. Uses saved predictions without refitting models.
# From the repository root: Rscript scripts/plot_cv_weekly_mos.R results/main/cv_week_mos.csv figures
# In RStudio, set the working directory to the repository root and click Source.

# --- EDITABLE SETTINGS ---
INPUT_FILE <- "results/main/cv_week_mos.csv"
OUTPUT_DIR <- "figures"
FILE_STEM <- "cv_weekly_mos"
X_LABEL <- "Study week"
Y_LABEL <- "Weekly MOS (1-5)"
# %d is replaced by the number of individual weekly ratings.
PANEL_TITLES <- c(all_weeks = "All valid rated weeks (N = %d)",
                  complete = "Six daily videos observed (N = %d)")
SERIES_LABELS <- c(observed = "Observed MOS", opt_k = "Kumaraswamy",
                   mean = "Unweighted mean", odd = "Odd days")
MODELS <- list(all_weeks = c("opt_k", "mean"),
               complete = c("odd", "opt_k"))
# Styles follow the series position in each panel: observed, model 1, model 2.
COLORS <- c("black", "#246D96", "#D4573B")
LINE_TYPES <- c(1, 2, 3)
POINT_TYPES <- c(16, 15, 15)
LINE_WIDTH <- 2
POINT_SIZE <- 0.65
TEXT_SIZE <- 1.05
LEGEND_SIZE <- 0.85
FONT_FAMILY <- "sans"
Y_LIMITS <- c(1.7, 4.65)
Y_TICKS <- seq(2, 4.5, 0.5)
X_TICKS <- c(1, 5, 9, 13, 17, 21, 25)
WIDTH_IN <- 10
HEIGHT_IN <- 7.6
PNG_DPI <- 300

# --- LOAD AND VALIDATE DATA ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2) stop("Usage: Rscript plot_cv_weekly_mos.R INPUT.csv OUTPUT_DIR")
if (length(args) >= 1) INPUT_FILE <- args[1]
if (length(args) >= 2) OUTPUT_DIR <- args[2]
if (!file.exists(INPUT_FILE)) stop("File not found: ", INPUT_FILE)
d <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
required <- c("dataset", "model", "week", "n_votes", "actual_mos", "predicted_mos")
if (!all(required %in% names(d)))
  stop("Missing columns: ", paste(setdiff(required, names(d)), collapse = ", "))
panels <- lapply(names(MODELS), function(ds) {
  models <- MODELS[[ds]]
  if (length(models) != 2 || anyDuplicated(models))
    stop("Each panel must specify two distinct models.")
  z <- d[d$dataset == ds & d$model %in% models, required, drop = FALSE]
  if (!nrow(z) || anyNA(z)) stop("Missing data or NA values in panel: ", ds)
  nums <- c("week", "n_votes", "actual_mos", "predicted_mos")
  if (!all(vapply(z[nums], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(z[nums])))) stop("Invalid numeric values: ", ds)
  if (any(z$n_votes < 1) || any(z$n_votes != floor(z$n_votes)))
    stop("Invalid rating counts: ", ds)
  if (anyDuplicated(z[c("model", "week")])) stop("Duplicate model/week rows: ", ds)
  a <- z[z$model == models[1], ]; a <- a[order(a$week), ]
  b <- z[z$model == models[2], ]; b <- b[order(b$week), ]
  if (!identical(a$week, 1:25) || !identical(b$week, 1:25))
    stop("Expected weeks 1-25 for both models: ", ds)
  if (!identical(a$n_votes, b$n_votes) ||
      any(abs(a$actual_mos - b$actual_mos) > 1e-8))
    stop("Models are not evaluated on the same weekly outcomes: ", ds)
  values <- c(a$actual_mos, a$predicted_mos, b$predicted_mos)
  if (any(values < 1 | values > 5)) stop("MOS outside the 1-5 scale: ", ds)
  if (any(values < Y_LIMITS[1] | values > Y_LIMITS[2]))
    stop("Y_LIMITS clips the data; expand the y-axis range.")
  list(dataset = ds, week = a$week, N = sum(a$n_votes),
       values = cbind(a$actual_mos, a$predicted_mos, b$predicted_mos),
       series = c("observed", models))
})

# --- DRAW THE TWO PANELS ---
draw_plot <- function() {
  par(mfrow = c(2, 1), mar = c(3.4, 4.7, 2.3, 0.8),
      family = FONT_FAMILY, cex = TEXT_SIZE, las = 1)
  for (i in seq_along(panels)) {
    p <- panels[[i]]
    plot(p$week, p$values[, 1], type = "n", xlim = c(1, 25), ylim = Y_LIMITS,
         xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    abline(h = Y_TICKS, col = "#DDDDDD", lwd = 0.8)
    axis(1, at = X_TICKS, labels = if (i == length(panels)) X_TICKS else FALSE)
    axis(2, at = Y_TICKS)
    box()
    title(main = sprintf(PANEL_TITLES[[p$dataset]], p$N), adj = 0, cex.main = 1)
    mtext(Y_LABEL, side = 2, line = 3.1, las = 0)
    if (i == length(panels)) mtext(X_LABEL, side = 1, line = 2.1)
    for (j in 1:3)
      lines(p$week, p$values[, j], type = "o", col = COLORS[j],
            lty = LINE_TYPES[j], pch = POINT_TYPES[j],
            lwd = LINE_WIDTH, cex = POINT_SIZE)
    legend("bottomleft", legend = unname(SERIES_LABELS[p$series]),
           col = COLORS, lty = LINE_TYPES, pch = POINT_TYPES,
           lwd = LINE_WIDTH, pt.cex = POINT_SIZE,
           cex = LEGEND_SIZE, horiz = TRUE, bty = "n")
  }
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(kind) {
  path <- file.path(OUTPUT_DIR, paste0(FILE_STEM, ".", kind))
  if (kind == "pdf") {
    pdf(path, width = WIDTH_IN, height = HEIGHT_IN, useDingbats = FALSE)
  } else {
    png(path, width = WIDTH_IN, height = HEIGHT_IN, units = "in", res = PNG_DPI)
  }
  on.exit(dev.off())
  draw_plot()
  message("Saved: ", normalizePath(path, mustWork = FALSE))
}
save_plot("pdf")
save_plot("png")
