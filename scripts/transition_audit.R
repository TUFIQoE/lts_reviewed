#!/usr/bin/env Rscript
# Descriptive audit of designed quality transitions and odd/even predictions.
# Usage: Rscript transition_audit.R INPUT_CSV OOF_PREDICTIONS_CSV OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript transition_audit.R INPUT_CSV OOF_PREDICTIONS_CSV OUTPUT_DIR")
d <- read.csv(args[1], stringsAsFactors = FALSE)
p <- read.csv(args[2], stringsAsFactors = FALSE)
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
days <- c("mo","tu","we","th","fr","sa")
d <- d[d$n == 7L & !(d$external_id == "343ab" & d$week == 14L), , drop=FALSE]
if (anyNA(d[days])) stop("Incomplete daily sequence in n==7 weeks")
if (length(unique(d$week)) != 25L) stop("Expected 25 weeks")
variation <- aggregate(d[days], by=list(week=d$week),
                       FUN=function(z) length(unique(z)))
if (any(as.matrix(variation[days]) != 1L)) stop("Daily quality differs between participants within a week")
weeks <- aggregate(d[days], by=list(week=d$week), FUN=function(z) z[1])
weeks <- weeks[order(weeks$week), , drop=FALSE]
q <- as.matrix(weeks[days]); storage.mode(q) <- "double"
gap <- q[, 2:6, drop=FALSE] - q[, 1:5, drop=FALSE]
names_gap <- c("change_12","change_23","change_34","change_45","change_56")
colnames(gap) <- names_gap
audit <- cbind(weeks, as.data.frame(gap))
audit$big_between_pairs <- as.integer(abs(gap[,2]) > 0.5) + as.integer(abs(gap[,4]) > 0.5)
audit$big_within_pairs <- rowSums(abs(gap[,c(1,3,5),drop=FALSE]) > 0.5)
audit$total_big_changes <- audit$big_between_pairs + audit$big_within_pairs
audit$max_between_change <- pmax(abs(gap[,2]),abs(gap[,4]))
audit$total_between_change <- abs(gap[,2]) + abs(gap[,4])
write.csv(audit, file.path(out,"week_transition_design.csv"), row.names=FALSE)

# Paired prediction errors, averaged by the 25 actual stimulus sequences.
# Each participant-week contributes once to each model; users are held out in
# the original CV. This is descriptive, since change positions were designed.
odd <- p[p$dataset=="complete" & p$model=="odd",
         c("record_id","week","observed","predicted")]
even <- p[p$dataset=="complete" & p$model=="even",
          c("record_id","week","observed","predicted")]
names(odd)[4] <- "odd_pred"; names(even)[4] <- "even_pred"
z <- merge(odd, even[,c("record_id","even_pred")], by="record_id", sort=FALSE)
if (nrow(z)!=nrow(odd) || anyDuplicated(z$record_id)) stop("Incomplete paired odd/even CV predictions")
z$odd_squared_error <- (z$observed-z$odd_pred)^2
z$even_squared_error <- (z$observed-z$even_pred)^2
z$odd_minus_even_squared_error <- z$odd_squared_error-z$even_squared_error
err <- aggregate(z[c("odd_squared_error","even_squared_error","odd_minus_even_squared_error")],
                 by=list(week=z$week), FUN=mean)
err$n_votes <- as.integer(table(z$week)[as.character(err$week)])
audit <- merge(audit,err,by="week",sort=TRUE)
write.csv(audit,file.path(out,"week_transition_and_cv_errors.csv"),row.names=FALSE)
grouped <- aggregate(audit[c("odd_minus_even_squared_error","total_between_change")],
                     by=list(big_between_pairs=audit$big_between_pairs),FUN=mean)
grouped$n_weeks <- as.integer(table(audit$big_between_pairs)[as.character(grouped$big_between_pairs)])
write.csv(grouped,file.path(out,"transition_group_descriptives.csv"),row.names=FALSE)
cat("Audited 25 weekly sequences. Big changes threshold = 0.5 P.1204.4 units.\n")
cat("Large changes between pairs: ",sum(abs(gap[,c(2,4)])>0.5),
    "; within pairs: ",sum(abs(gap[,c(1,3,5)])>0.5),".\n",sep="")
cat("Descriptive associations cannot identify a transition-memory mechanism.\n")
