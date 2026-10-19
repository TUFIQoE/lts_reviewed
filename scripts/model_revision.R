#!/usr/bin/env Rscript
# Reanalysis of the 25-week QoE study. Run from any directory:
# Rscript model_revision.R path/to/phase4_current_results_p1204.csv output_dir
# Dependencies: ordinal (CRAN). R >= 4.1.

if (!requireNamespace("ordinal", quietly = TRUE)) {
  stop("Install the 'ordinal' package: install.packages('ordinal')")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 5L) stop("Usage: Rscript model_revision.R INPUT_CSV OUTPUT_DIR [BOOT_REPS] [all|screened_out] [FULL_FOLDS_CSV]")
input <- args[1]; out <- args[2]
B <- if (length(args) >= 3L) as.integer(args[3]) else 300L
if (is.na(B) || B < 0L) stop("BOOT_REPS must be a nonnegative integer")
cohort <- if (length(args) >= 4L) args[4] else "all"
if (!cohort %in% c("all", "screened_out")) stop("Cohort must be all or screened_out")
fold_input <- if (length(args) >= 5L) args[5] else NULL
if (cohort == "screened_out" && is.null(fold_input))
  stop("For screened_out, supply folds.csv from the corrected all-participant run")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
days <- c("mo", "tu", "we", "th", "fr", "sa")
needed <- c("external_id", "week", "n", "q2", days)
if (!all(needed %in% names(d))) stop("Missing columns: ", paste(setdiff(needed, names(d)), collapse = ", "))
if (anyDuplicated(d[c("external_id", "week")])) stop("Duplicate participant-week rows")
if (any(!is.na(d$q2) & !(d$q2 %in% 1:5))) stop("q2 must be in 1:5")
d <- d[!is.na(d$q2) & !is.na(d$external_id), , drop = FALSE]
# Source-log audit: participant 343ab, week 14 has no dedicated weekly rating
# video. The q2=3 value belongs to daily video 099, not to the weekly outcome.
invalid_outcome <- d[d$external_id == "343ab" & d$week == 14L, , drop = FALSE]
if (nrow(invalid_outcome) != 1L || invalid_outcome$q2 != 3L || invalid_outcome$n != 6L) {
  stop("Expected anomalous source record 343ab/week 14 has changed; re-audit before analysis")
}
write.csv(invalid_outcome, file.path(out, "excluded_invalid_outcome.csv"), row.names=FALSE)
d <- d[!(d$external_id == "343ab" & d$week == 14L), , drop = FALSE]
d$os <- ordered(d$q2, levels = 1:5)
d$record_id <- seq_len(nrow(d))
screened_ids <- c("f1308", "bb8d0", "343ab", "9d7d1", "3fe63", "57ea9")
if (!all(screened_ids %in% d$external_id)) stop("Screened participant IDs changed; re-audit screening")
screened_counts <- data.frame(external_id=screened_ids,
  n_weekly_ratings=vapply(screened_ids, function(id) sum(d$external_id == id), integer(1)),
  n_complete_weeks=vapply(screened_ids,
    function(id) sum(d$external_id == id & d$n == 7L), integer(1)))
write.csv(screened_counts, file.path(out, "screened_participant_counts.csv"), row.names=FALSE)
if (cohort == "screened_out") {
  d <- d[!(d$external_id %in% screened_ids), , drop=FALSE]
  if (any(d$external_id %in% screened_ids) ||
      length(unique(d$external_id)) != 32L || nrow(d) != 709L)
    stop("Screening exclusion failed or cohort counts changed; inspect input")
}
if (any(rowSums(!is.na(d[days])) == 0L)) stop("At least one row has no daily quality value")

# The paper's n=7 means a complete week: six daily videos plus final answer.
# Report any remaining disagreement between n and daily fields.
discrepant <- d[(d$n == 7L) != (rowSums(!is.na(d[days])) == 6L),
                c("external_id", "week", "n", days), drop = FALSE]
if (nrow(discrepant)) write.csv(discrepant, file.path(out, "n_daily_discrepancies.csv"), row.names=FALSE)

add_predictors <- function(x) {
  m <- as.matrix(x[days]); storage.mode(m) <- "double"
  x$mean_q <- rowMeans(m, na.rm = TRUE)
  x$min_q <- apply(m, 1L, function(z) min(z, na.rm = TRUE))
  x$m12 <- rowMeans(m[, 1:2, drop = FALSE], na.rm = TRUE)
  x$m34 <- rowMeans(m[, 3:4, drop = FALSE], na.rm = TRUE)
  x$m56 <- rowMeans(m[, 5:6, drop = FALSE], na.rm = TRUE)
  x
}
d <- add_predictors(d)

steps <- (0:6) / 6
kum_weights <- function(a, b) diff(1 - (1 - steps^a)^b)
kum_score <- function(x, a, b) {
  m <- as.matrix(x[days]); storage.mode(m) <- "double"
  w <- matrix(kum_weights(a, b), nrow(m), 6L, byrow = TRUE)
  w[is.na(m)] <- 0
  den <- rowSums(w)
  if (any(den <= 0)) stop("Kumaraswamy: zero observed weight")
  rowSums(m * w, na.rm = TRUE) / den
}

forms <- list(
  all = os ~ mo + tu + we + th + fr + sa,
  agg = os ~ m12 + m34 + m56,
  odd = os ~ mo + we + fr,
  even = os ~ tu + th + sa,
  mean = os ~ mean_q,
  min = os ~ min_q,
  opt_k = os ~ kum_q
)
starts <- rbind(c(1,1), c(.6,1), c(.9,1), c(1,.5), c(5,.9),
                c(.5,.5), c(.9,.9), c(5,5), c(3,6), c(2,6))

fit_clm <- function(formula, x) {
  fit <- suppressWarnings(ordinal::clm(formula, data = x, link = "logit", Hess = FALSE))
  if (!is.finite(as.numeric(logLik(fit)))) stop("Nonfinite log-likelihood")
  fit
}
fit_model <- function(name, train) {
  if (name != "opt_k") return(list(fit = fit_clm(forms[[name]], train), a = NA_real_, b = NA_real_))
  objective <- function(p) {
    z <- train
    z$kum_q <- tryCatch(kum_score(z, p[1], p[2]), error = function(e) rep(NA_real_, nrow(z)))
    if (any(!is.finite(z$kum_q))) return(1e12)
    tryCatch(as.numeric(AIC(fit_clm(forms$opt_k, z))), error = function(e) 1e12)
  }
  candidates <- lapply(seq_len(nrow(starts)), function(i) {
    tryCatch(optim(starts[i, ], objective, method = "L-BFGS-B",
                   lower = c(1e-5, 1e-5), upper = c(100, 100),
                   control = list(maxit = 1000, factr = 1e7)), error = function(e) NULL)
  })
  ok <- Filter(function(z) !is.null(z) && z$convergence == 0L &&
                 is.finite(z$value) && z$value < 1e11, candidates)
  if (!length(ok)) stop("Kumaraswamy optimization failed for all starts")
  best <- ok[[which.min(vapply(ok, `[[`, numeric(1), "value"))]]
  train$kum_q <- kum_score(train, best$par[1], best$par[2])
  list(fit = fit_clm(forms$opt_k, train), a = best$par[1], b = best$par[2],
       objective = best$value, successful_starts = length(ok), convergence = best$convergence)
}

prob_matrix <- function(fit, newdata) {
  # ordinal::predict.clm returns only P(observed category) when the response
  # column is present. Omit it to request the full 1..5 probability matrix.
  newdata$os <- NULL
  z <- predict(fit, newdata = newdata, type = "prob")
  if (is.list(z) && !is.null(z$fit)) z <- z$fit
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || nrow(z) != nrow(newdata)) stop("Unexpected clm probability output")
  p <- matrix(0, nrow(z), 5L)
  if (!is.null(colnames(z))) {
    nm <- gsub("[^1-5]", "", colnames(z))
    if (length(unique(nm)) != ncol(z) || any(!nm %in% as.character(1:5)))
      stop("Unexpected clm category labels")
    p[, as.integer(nm)] <- z
  } else if (ncol(z) == 5L) p <- z else stop("Unexpected number of clm categories")
  if (any(!is.finite(p)) || any(abs(rowSums(p) - 1) > 1e-5)) stop("Invalid predicted probabilities")
  p
}
predict_model <- function(object, name, test) {
  if (name == "opt_k") test$kum_q <- kum_score(test, object$a, object$b)
  as.vector(prob_matrix(object$fit, test) %*% (1:5))
}

metric <- function(y, yhat) {
  c(N = length(y), RMSE = sqrt(mean((y-yhat)^2)), MAE = mean(abs(y-yhat)),
    PLCC = if (sd(y) > 0 && sd(yhat) > 0) cor(y,yhat) else NA_real_,
    SRCC = if (length(unique(y)) > 1 && length(unique(yhat)) > 1)
      suppressWarnings(cor(y,yhat,method="spearman")) else NA_real_)
}
write <- function(x, filename) write.csv(x, file.path(out, filename), row.names = FALSE, na = "")
set.seed(20260925)
if (is.null(fold_input)) {
  ids <- sample(unique(d$external_id))
  fold_by_id <- setNames(rep(1:5, length.out = length(ids)), ids)
} else {
  prior <- read.csv(fold_input, stringsAsFactors=FALSE)
  if (!all(c("external_id", "fold") %in% names(prior)) ||
      anyDuplicated(prior$external_id) || anyNA(prior[c("external_id", "fold")]) ||
      any(!prior$fold %in% 1:5) ||
      !all(unique(d$external_id) %in% prior$external_id))
    stop("Invalid folds CSV: expected unique participant IDs and folds 1..5 covering the cohort")
  fold_by_id <- setNames(prior$fold, prior$external_id)
  fold_by_id <- fold_by_id[names(fold_by_id) %in% unique(d$external_id)]
}
write(data.frame(external_id = names(fold_by_id), fold = unname(fold_by_id)), "folds.csv")

sets <- list(complete = list(data = d[d$n == 7L, , drop = FALSE], names = names(forms)),
             all_weeks = list(data = d, names = c("mean", "min", "opt_k")))
all_preds <- list(); fits <- list(); errors <- list()
for (set_name in names(sets)) {
  x <- sets[[set_name]]$data
  if (length(unique(x$external_id)) < 5L) stop("Too few participants for 5-fold CV")
  for (fold in 1:5) {
    tr <- x[fold_by_id[x$external_id] != fold, , drop = FALSE]
    te <- x[fold_by_id[x$external_id] == fold, , drop = FALSE]
    if (!nrow(tr) || !nrow(te)) stop("Empty training or test fold in ", set_name, "/", fold)
    for (name in sets[[set_name]]$names) {
      result <- tryCatch({
        f <- fit_model(name, tr)
        pred <- predict_model(f, name, te)
        data.frame(dataset = set_name, model = name, fold = fold,
                   record_id = te$record_id, external_id = te$external_id,
                   week = te$week, n = te$n, observed = te$q2, predicted = pred)
      }, error = function(e) {
        errors[[length(errors)+1L]] <<- data.frame(dataset = set_name, model = name,
          fold = fold, message = conditionMessage(e)); NULL
      })
      if (!is.null(result)) all_preds[[length(all_preds)+1L]] <- result
    }
  }
  for (name in sets[[set_name]]$names) {
    result <- tryCatch({
      f <- fit_model(name, x)
      data.frame(dataset = set_name, model = name, N = nrow(x),
                 AIC = AIC(f$fit) + if (name == "opt_k") 4 else 0,
                 a = f$a, b = f$b,
                 successful_starts = if (name == "opt_k") f$successful_starts else NA_integer_)
    }, error = function(e) {
      errors[[length(errors)+1L]] <<- data.frame(dataset = set_name, model = name,
        fold = 0L, message = conditionMessage(e)); NULL
    })
    if (!is.null(result)) fits[[length(fits)+1L]] <- result
  }
}
if (length(errors)) {
  write(do.call(rbind, errors), "fit_errors.csv")
  stop("At least one model failed. Inspect fit_errors.csv; no partial CV summary was issued")
}
p <- do.call(rbind, all_preds)
for (s in names(sets)) for (name in sets[[s]]$names) {
  actual <- p[p$dataset == s & p$model == name, "record_id"]
  expected <- sets[[s]]$data$record_id
  if (length(actual) != length(expected) || !setequal(actual, expected) || anyDuplicated(actual))
    stop("Incomplete out-of-fold predictions for ", s, "/", name)
}
write(p, "oof_predictions.csv")
write(do.call(rbind, fits), "full_sample_fits.csv")
summary_rows <- lapply(split(p, list(p$dataset, p$model), drop = TRUE), function(z)
  data.frame(dataset = z$dataset[1], model = z$model[1], t(metric(z$observed,z$predicted))))
write(do.call(rbind, summary_rows), "cv_metrics.csv")

# Week-level MOS is a separate target, based on held-out participant ratings.
key <- interaction(p$dataset, p$model, p$week, drop = TRUE)
week_rows <- lapply(split(p, key), function(z) data.frame(dataset=z$dataset[1],
  model=z$model[1], week=z$week[1], n_votes=nrow(z),
  actual_mos=mean(z$observed), predicted_mos=mean(z$predicted)))
week_mos <- do.call(rbind, week_rows)
write(week_mos, "cv_week_mos.csv")
write(do.call(rbind, lapply(split(week_mos, list(week_mos$dataset,week_mos$model), drop=TRUE),
  function(z) data.frame(dataset=z$dataset[1], model=z$model[1],
                         t(metric(z$actual_mos,z$predicted_mos))))), "cv_week_metrics.csv")

# Paired uncertainty of RMSE/MAE difference (opt_k minus comparator), resampling
# participants; the CV predictions are held fixed. This is conditional on the
# one chosen CV split and does not account for retraining instability.
if (B > 0L) {
  set.seed(20260926)
  comparisons <- list()
  for (s in names(sets)) {
    z <- p[p$dataset == s, ]
    ref <- z[z$model == "opt_k", c("record_id","external_id","observed","predicted")]
    names(ref)[4] <- "reference"
    for (name in setdiff(sets[[s]]$names, "opt_k")) {
      alt <- z[z$model == name, c("record_id","predicted")]
      names(alt)[2] <- "comparator"
      joined <- merge(ref, alt, by = "record_id", sort = FALSE)
      if (nrow(joined) != nrow(ref)) stop("Unpaired prediction rows")
      delta <- function(v, metric_name) {
        a <- v$observed-v$reference; b <- v$observed-v$comparator
        if (metric_name == "RMSE") sqrt(mean(a*a))-sqrt(mean(b*b))
        else mean(abs(a))-mean(abs(b))
      }
      participants <- unique(joined$external_id)
      boot <- replicate(B, {
        sampled <- sample(participants, length(participants), replace = TRUE)
        v <- do.call(rbind, lapply(sampled, function(id) joined[joined$external_id == id,]))
        c(RMSE=delta(v,"RMSE"), MAE=delta(v,"MAE"))
      })
      for (met in c("RMSE","MAE")) comparisons[[length(comparisons)+1L]] <- data.frame(
        dataset=s, comparator=name, metric=met, estimate=delta(joined,met),
        lower=unname(quantile(boot[met,],.025)), upper=unname(quantile(boot[met,],.975)), B=B)
    }
  }
  write(do.call(rbind, comparisons), "paired_metric_differences.csv")

  # Nonparametric participant-cluster bootstrap; refit opt_k for each resample.
  # Keep duplicate selected clusters with their entire participant histories.
  for (s in names(sets)) {
    x <- sets[[s]]$data
    participants <- unique(x$external_id)
    estimates <- vector("list", B)
    for (i in seq_len(B)) {
      sampled <- sample(participants, length(participants), replace = TRUE)
      bx <- do.call(rbind, lapply(sampled, function(id) x[x$external_id == id, , drop=FALSE]))
      estimates[[i]] <- tryCatch({
        f <- fit_model("opt_k", bx)
        data.frame(rep=i, a=f$a, b=f$b, t(kum_weights(f$a,f$b)))
      }, error=function(e) NULL)
    }
    good <- Filter(Negate(is.null), estimates)
    if (length(good) < .9*B) stop("Too many failed bootstrap fits for ",s)
    z <- do.call(rbind, good)
    names(z)[4:9] <- paste0("day",1:6)
    write(z, paste0("bootstrap_",s,".csv"))
    orig <- fits[[which(vapply(fits, function(f) f$dataset==s && f$model=="opt_k", logical(1)))]]
    original <- c(a=orig$a,b=orig$b,setNames(kum_weights(orig$a,orig$b),paste0("day",1:6)))
    write(data.frame(dataset=s, parameter=names(original), estimate=as.numeric(original),
      lower=vapply(names(original),function(k) unname(quantile(z[[k]],.025)),numeric(1)),
      upper=vapply(names(original),function(k) unname(quantile(z[[k]],.975)),numeric(1)),
      successful=length(good), attempted=B), paste0("parameter_weights_ci_",s,".csv"))
  }
}
cat("Done. Cohort: ", cohort, "; participants: ", length(unique(d$external_id)),
    "; weekly ratings: ", nrow(d), ". Outputs in: ", normalizePath(out), "\n", sep="")
