#!/usr/bin/env Rscript
# Compare grouped-CV predictions trained with all eligible participants versus
# refitted after excluding the six participants listed in manuscript screening.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L)
  stop("Usage: Rscript screening_compare.R FULL_OOF_CSV CLEAN_OOF_CSV OUTPUT_DIR")
out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
full <- read.csv(args[1], stringsAsFactors=FALSE)
clean <- read.csv(args[2], stringsAsFactors=FALSE)
required <- c("dataset","model","fold","record_id","external_id","week","observed","predicted")
if (!all(required %in% names(full)) || !all(required %in% names(clean)))
  stop("Missing out-of-fold prediction columns")
key <- function(x) paste(x$dataset, x$model, x$record_id, sep="|")
if (anyDuplicated(key(full)) || anyDuplicated(key(clean))) stop("Duplicate prediction keys")
idx <- match(key(clean), key(full))
if (anyNA(idx) || any(full$fold[idx] != clean$fold) ||
    any(full$external_id[idx] != clean$external_id) ||
    any(full$week[idx] != clean$week) ||
    any(full$observed[idx] != clean$observed))
  stop("Full and screened-out predictions do not match on observations and folds")
screened_ids <- c("f1308", "bb8d0", "343ab", "9d7d1", "3fe63", "57ea9")
remaining <- intersect(unique(clean$external_id), screened_ids)
if (length(remaining)) stop("Screened participants remain in clean predictions: ",
  paste(remaining,collapse=", "), ".")
clean$predicted_full <- full$predicted[idx]
metric <- function(y, pred) c(RMSE=sqrt(mean((y-pred)^2)),
  MAE=mean(abs(y-pred)), PLCC=cor(y,pred),
  SRCC=suppressWarnings(cor(y,pred,method="spearman")))
splits <- split(clean, list(clean$dataset,clean$model), drop=TRUE)
rows <- lapply(splits, function(z) {
  a <- metric(z$observed,z$predicted_full)
  b <- metric(z$observed,z$predicted)
  data.frame(dataset=z$dataset[1], model=z$model[1],
    N=nrow(z), n_participants=length(unique(z$external_id)),
    train_all_RMSE=a["RMSE"], train_clean_RMSE=b["RMSE"], delta_RMSE=b["RMSE"]-a["RMSE"],
    train_all_MAE=a["MAE"], train_clean_MAE=b["MAE"], delta_MAE=b["MAE"]-a["MAE"],
    train_all_PLCC=a["PLCC"], train_clean_PLCC=b["PLCC"],
    train_all_SRCC=a["SRCC"], train_clean_SRCC=b["SRCC"])
})
write.csv(do.call(rbind,rows), file.path(out,"screening_matched_cv.csv"), row.names=FALSE)
cat("Matched CV comparison saved in: ", normalizePath(out), "\n", sep="")
