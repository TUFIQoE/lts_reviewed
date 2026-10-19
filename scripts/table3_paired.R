#!/usr/bin/env Rscript
# Recalculate manuscript Table 3 with participant-matched weekly ratings.
# Usage: Rscript table3_paired.R INPUT_CSV OUTPUT_DIR
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: Rscript table3_paired.R INPUT_CSV OUTPUT_DIR")
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- read.csv(args[1], stringsAsFactors=FALSE)
required <- c("external_id","week","n","q2")
if (!all(required %in% names(x))) stop("Missing columns: ",paste(setdiff(required,names(x)),collapse=", "))
x <- x[x$n==7L & !is.na(x$q2),required,drop=FALSE]
if (anyDuplicated(x[c("external_id","week")])) stop("Duplicate participant-week rows")
wide <- reshape(x[c("external_id","week","q2")],
                idvar="external_id", timevar="week", direction="wide")
pair_sets <- list(
  list(c(4,3)), list(c(7,6)), list(c(9,10)), list(c(12,13)),
  list(c(25,20)), list(c(16,23)), list(c(17,22)), list(c(21,15)),
  list(c(4,3),c(25,20)), list(c(7,6),c(16,23)),
  list(c(9,10),c(17,22)), list(c(12,13),c(21,15))
)
row_label <- c("4 vs 3","7 vs 6","9 vs 10","12 vs 13",
               "25 vs 20","16 vs 23","17 vs 22","21 vs 15",
               "4 vs 3 AND 25 vs 20","7 vs 6 AND 16 vs 23",
               "9 vs 10 AND 17 vs 22","12 vs 13 AND 21 vs 15")
all_rows <- list(); participant_rows <- list()
for (i in seq_along(pair_sets)) {
  comps <- pair_sets[[i]]
  diffs <- sapply(comps, function(pair) {
    left <- paste0("q2.", pair[1]); right <- paste0("q2.", pair[2])
    if (!all(c(left,right) %in% names(wide))) stop("Missing week columns")
    wide[[left]] - wide[[right]]
  })
  if (is.null(dim(diffs))) diffs <- matrix(diffs,ncol=1)
  # In pooled rows 9-12, require both paired comparisons for each participant,
  # then give nature and slow-motion equal weight within that participant.
  valid <- rowSums(is.na(diffs))==0L
  delta <- rowMeans(diffs[valid,,drop=FALSE])
  if (length(delta)<2L || sd(delta)==0) stop("Insufficient matched data in row ",i)
  test <- t.test(delta,mu=0,alternative="two.sided")
  all_rows[[i]] <- data.frame(row=i, comparison=row_label[i],
    analysis=if(i<=8L) "single paired comparison" else "complete-case mean of two paired differences",
    n_participants=length(delta), mean_difference=mean(delta),
    sd_difference=sd(delta), t=unname(test$statistic), df=unname(test$parameter),
    p=unname(test$p.value), ci_lower=test$conf.int[1], ci_upper=test$conf.int[2],
    cohen_dz=mean(delta)/sd(delta))
  participant_rows[[i]] <- data.frame(row=i, comparison=row_label[i],
    external_id=wide$external_id[valid], difference=delta)
}
tab <- do.call(rbind,all_rows)
tab$p_holm_12 <- p.adjust(tab$p,method="holm")
tab$p_holm_primary_8 <- NA_real_
tab$p_holm_primary_8[1:8] <- p.adjust(tab$p[1:8],method="holm")
write.csv(tab,file.path(out,"table3_paired_results.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,participant_rows),
          file.path(out,"table3_participant_differences.csv"),row.names=FALSE)
cat("Recalculated 12 Table 3 rows. Rows 9-12 use participants with both paired comparisons.\n")
cat("Rows 9-12 are exploratory: consider omitting them from the main table.\n")
