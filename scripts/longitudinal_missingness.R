#!/usr/bin/env Rscript
# Reproducible descriptive audit and exploratory within-participant time trends.
# Usage: Rscript longitudinal_missingness.R INPUT_CSV OUTPUT_DIR [BOOT_REPS]
# Base R only; the bootstrap samples participants with their entire histories.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 3L)
  stop("Usage: Rscript longitudinal_missingness.R INPUT_CSV OUTPUT_DIR [BOOT_REPS]")
input <- args[1]; out <- args[2]
B <- if (length(args) == 3L) as.integer(args[3]) else 2000L
if (is.na(B) || B < 100L) stop("BOOT_REPS must be at least 100")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_csv <- function(x, name) write.csv(x, file.path(out,name), row.names=FALSE, na="")
d <- read.csv(input, stringsAsFactors=FALSE, check.names=FALSE)
days <- c("mo","tu","we","th","fr","sa")
if (!all(c("external_id","week","n","q2",days) %in% names(d))) stop("Missing source columns")
if (anyDuplicated(d[c("external_id","week")]) ||
    anyNA(d[c("external_id","week","n","q2")]) ||
    !identical(sort(unique(d$week)),1:25)) stop("Unexpected participant-week structure")
bad <- d$external_id == "343ab" & d$week == 14L
if (sum(bad) != 1L || d$q2[bad] != 3L || d$n[bad] != 6L)
  stop("Expected source-log anomaly 343ab/week 14 has changed")
save_csv(d[bad,,drop=FALSE],"excluded_invalid_outcome.csv")
d <- d[!bad,,drop=FALSE]
if (length(unique(d$external_id)) != 38L || nrow(d) != 828L ||
    any(d$n != 1L + rowSums(!is.na(d[days]))))
  stop("Corrected source counts or daily availability have changed")

# The assigned objective quality for each week/day is common to participants.
# A missing daily value indicates no daily observation in the modeling CSV.
profiles <- do.call(rbind,lapply(1:25,function(w) {
  z <- d[d$week == w,,drop=FALSE]
  values <- vapply(days,function(day) {
    u <- unique(z[[day]][!is.na(z[[day]])])
    if (length(u) != 1L) stop("Assigned quality is not unique: week ",w,"/",day)
    u
  },numeric(1))
  data.frame(week=w,t(values),check.names=FALSE)
}))
save_csv(profiles,"assigned_quality_by_week.csv")
d$mean_seen <- rowMeans(d[days],na.rm=TRUE)
d$mean_assigned <- rowMeans(profiles[match(d$week,profiles$week),days])
d$content_block <- as.integer(d$week >= 14L)
grid <- expand.grid(external_id=sort(unique(d$external_id)),week=1:25,
                    KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
key <- function(z) paste(z$external_id,z$week,sep="|")
grid$observed_week <- as.integer(key(grid) %in% key(d))
grid$content_block <- as.integer(grid$week >= 14L)
save_csv(grid,"participant_week_availability.csv")

long <- do.call(rbind,lapply(days,function(day) data.frame(
  external_id=d$external_id,week=d$week,day=day,
  content_block=d$content_block,
  objective_quality=profiles[[day]][match(d$week,profiles$week)],
  daily_observed=as.integer(!is.na(d[[day]])))))
long$low_objective_quality <- as.integer(long$objective_quality < 3)
save_csv(long,"daily_availability.csv")
count_days <- function(z) data.frame(N=nrow(z),missing=sum(z$daily_observed == 0L),
                                    missing_fraction=mean(z$daily_observed == 0L))
by_day <- do.call(rbind,lapply(days,function(day)
  cbind(day=day,count_days(long[long$day == day,,drop=FALSE]))))
by_quality <- do.call(rbind,lapply(0:1,function(q)
  cbind(low_objective_quality=q,count_days(long[long$low_objective_quality == q,,drop=FALSE]))))
by_block <- do.call(rbind,lapply(0:1,function(b)
  cbind(content_block=b,count_days(long[long$content_block == b,,drop=FALSE]))))
save_csv(by_day,"missing_daily_by_day.csv")
save_csv(by_quality,"missing_daily_by_quality.csv")
save_csv(by_block,"missing_daily_by_content_block.csv")

weeks <- do.call(rbind,lapply(1:25,function(w) {
  x <- d[d$week == w,,drop=FALSE]
  data.frame(week=w,content_block=ifelse(w <= 13,"nature","slow_motion"),
    eligible=38L,valid_weekly_ratings=nrow(x),
    absent_or_invalid_weekly_ratings=38L-nrow(x),
    weekly_response_fraction=nrow(x)/38,
    complete_rated_weeks=sum(x$n == 7L),
    incomplete_rated_weeks=sum(x$n < 7L),
    missing_daily_among_rated_weeks=sum(7L-x$n),
    daily_opportunities_among_rated_weeks=6L*nrow(x),
    weekly_MOS=mean(x$q2),
    assigned_quality_mean=mean(profiles[w,days]),
    seen_quality_mean=mean(x$mean_seen))
}))
save_csv(weeks,"weekly_availability_and_mos.csv")
overview <- data.frame(
  measure=c("eligible_participants","weeks","potential_participant_weeks",
    "valid_weekly_ratings","weekly_without_valid_rating","weekly_rows_absent_from_source",
    "weekly_rows_with_invalid_outcome","complete_rated_weeks","incomplete_rated_weeks",
    "daily_opportunities_with_valid_weekly_rating","daily_observed_with_valid_weekly_rating",
    "daily_missing_with_valid_weekly_rating"),
  value=c(38,25,nrow(grid),nrow(d),sum(grid$observed_week == 0L),
          nrow(grid)-nrow(d)-1L,1L,sum(d$n == 7L),sum(d$n < 7L),
          nrow(long),sum(long$daily_observed),sum(long$daily_observed == 0L)))
save_csv(overview,"missingness_overview.csv")

# Within-person (participant fixed-intercept) linear summaries on the 1-5
# rating scale. This is an exploratory trend, not an ordinal mixed model.
# Bootstrap CIs condition on the observed 25-week stimulus schedule; neither
# dropout nor the ordered content switch is assumed ignorable/causal.
set.seed(20260928)
trend <- function(x, response, predictors, label) {
  groups <- split(x,x$external_id)
  pieces <- lapply(groups,function(z) {
    X <- as.matrix(z[predictors]); storage.mode(X) <- "double"
    y <- as.numeric(z[[response]])
    X <- sweep(X,2L,colMeans(X),"-"); y <- y-mean(y)
    list(xx=crossprod(X),xy=crossprod(X,y))
  })
  estimate <- function(indices) {
    xx <- Reduce(`+`,lapply(indices,function(i) pieces[[i]]$xx))
    xy <- Reduce(`+`,lapply(indices,function(i) pieces[[i]]$xy))
    tryCatch(as.vector(solve(xx,xy)),error=function(e) rep(NA_real_,length(predictors)))
  }
  base <- estimate(seq_along(pieces))
  if (anyNA(base)) stop("Singular trend model: ",label)
  bt <- replicate(B,estimate(sample(seq_along(pieces),length(pieces),replace=TRUE)))
  if (is.vector(bt)) bt <- matrix(bt,nrow=length(predictors))
  valid <- colSums(!is.finite(bt)) == 0L
  if (sum(valid) < .9*B) stop("Too many singular bootstrap samples: ",label)
  bt <- bt[,valid,drop=FALSE]
  ci <- t(apply(bt,1L,quantile,probs=c(.025,.975)))
  data.frame(model=label,predictor=predictors,N=nrow(x),
    n_participants=length(groups),estimate=base,lower=ci[,1],upper=ci[,2],
    successful=ncol(bt),attempted=B)
}
fits <- list(
  trend(d,"q2","week","raw_week_with_participant_intercepts"),
  trend(d,"q2",c("week","mean_seen"),"week_adjusted_for_seen_quality"),
  trend(d,"q2",c("week","mean_assigned"),"week_adjusted_for_assigned_quality"),
  trend(d,"q2",c("week","mean_seen","content_block"),"week_quality_and_content_block"),
  trend(d[d$n == 7L,,drop=FALSE],"q2",c("week","mean_seen"),"complete_weeks_only"),
  trend(d[d$week <= 13L,,drop=FALSE],"q2",c("week","mean_seen"),"nature_block_only"),
  trend(d[d$week >= 14L,,drop=FALSE],"q2",c("week","mean_seen"),"slow_motion_block_only"),
  trend(grid,"observed_week","week","weekly_response_fraction"))
save_csv(do.call(rbind,fits),"exploratory_within_participant_trends.csv")
cat("Done: ",nrow(d)," valid weekly ratings; ",
    sum(grid$observed_week == 0L)," weeks without valid outcome; ",
    sum(long$daily_observed == 0L)," missing daily records among rated weeks.\n",
    "Outputs in ",normalizePath(out),"\n",sep="")
