# Weekly retrospective video QoE: analysis and reproduction

R code, modeling data, and reference results for **How was your week?
Modelling multi-episodic quality of experience judgements in a long-term study**.

The study follows 38 participants over 25 weeks. Participants were scheduled
to view six videos per week and provide one retrospective weekly quality
rating on a 1–5 scale. The analysis compares seven prediction models on
complete viewing histories and three models that accommodate missing days.

## Requirements

- R 4.1 or later, with `Rscript` available in the terminal.
- The CRAN package `ordinal`.

Install the dependency once in the R console:

```r
install.packages("ordinal")
```

All other scripts use base R. Run the commands below from the repository
root. `Rscript ...` commands belong in a system terminal, not at the R `>` prompt.

## Repository contents

| Path | Purpose |
| --- | --- |
| `data/phase4_current_results_p1204.csv` | Input participant–week data. |
| `data/folds.csv` | Fixed participant assignments to the five CV folds. |
| `scripts/` | Analysis, plotting, workflow, and result-comparison scripts. |
| `results/main/` | Reference primary analysis and bootstrap results. |
| `results/screened_out/` | Reference predictions after excluding six flagged participants. |
| `results/screening/` | Matched comparison of the two training cohorts. |
| `results/paired_tests/` | Paired week comparisons and participant differences. |
| `results/transitions/` | Descriptive quality-transition audit. |
| `results/longitudinal/` | Weekly availability/MOS and longitudinal trend results. |
| `figures/` | Reference weekly MOS figure in PDF and PNG. |

Reruns write to `recomputed/`, keeping the supplied reference results intact.
Repeated runs using the same output directory overwrite that rerun's outputs.

## Reproduce the analysis

### 1. Check that the workflow runs

```sh
Rscript scripts/run_all.R quick
```

This fits the models and performs grouped cross-validation, paired tests,
screening sensitivity, transition auditing, and plotting. It skips the
model bootstrap and uses only 100 resamples for the longitudinal trends.
Outputs are written to `recomputed/quick/`. This mode checks execution;
its bootstrap intervals are not intended for reporting.

### 2. Run the full analysis

```sh
Rscript scripts/run_all.R full
```

This uses 1000 participant-bootstrap resamples for the primary model analysis
and 5000 for the longitudinal trends. The screened-out sensitivity analysis
uses zero model-bootstrap resamples. The full model bootstrap can take
substantial time because the Kumaraswamy model is optimized from multiple
starting points in each resample.

Outputs are written to `recomputed/full/`. The driver also records
`sessionInfo.txt` and `run_settings.txt` to document the software environment
and run settings. If any step fails, the workflow stops.

### 3. Compare the rerun with the reference results

```sh
Rscript scripts/check_results.R recomputed/full
```

This compares the principal metric, parameter, interval, paired-test,
screening, transition, and trend tables. Rows are matched by their identifiers;
numeric comparisons use an absolute tolerance of `1e-6`. It reports each
comparison and exits with an error if a table is missing or differs.

Small numerical differences may occur across R or `ordinal` versions.
Inspect the reported columns and `sessionInfo.txt` before accepting differences.
An alternative absolute tolerance can be supplied explicitly:

```sh
Rscript scripts/check_results.R recomputed/full 0.00001
```

## Run individual steps

These commands use `recomputed/manual/` and should be executed in order.

### 1. Fit and validate the primary models

```sh
Rscript scripts/model_revision.R data/phase4_current_results_p1204.csv recomputed/manual/main 1000 all data/folds.csv
```

The third argument is the number of bootstrap resamples; use `0` to skip
the bootstrap. The fixed folds preserve the reference participant split.

### 2. Fit the screened-out sensitivity analysis

```sh
Rscript scripts/model_revision.R data/phase4_current_results_p1204.csv recomputed/manual/screened_out 0 screened_out data/folds.csv
```

### 3. Compare both training cohorts on the same held-out outcomes

```sh
Rscript scripts/screening_compare.R recomputed/manual/main/oof_predictions.csv recomputed/manual/screened_out/oof_predictions.csv recomputed/manual/screening
```

### 4. Calculate paired week comparisons

```sh
Rscript scripts/table3_paired.R data/phase4_current_results_p1204.csv recomputed/manual/paired_tests
```

### 5. Audit quality transitions and odd/even prediction errors

```sh
Rscript scripts/transition_audit.R data/phase4_current_results_p1204.csv recomputed/manual/main/oof_predictions.csv recomputed/manual/transitions
```

### 6. Calculate missingness summaries and longitudinal trends

```sh
Rscript scripts/longitudinal_missingness.R data/phase4_current_results_p1204.csv recomputed/manual/longitudinal 5000
```

### 7. Generate the weekly MOS figure

```sh
Rscript scripts/plot_cv_weekly_mos.R recomputed/manual/main/cv_week_mos.csv recomputed/manual/figures
```

To regenerate only the figure from the supplied reference results:

```sh
Rscript scripts/plot_cv_weekly_mos.R results/main/cv_week_mos.csv recomputed/figures
```

The plot script saves PDF and PNG files without refitting models. Labels,
colors, line styles, axis ranges, and dimensions are configured at its top.
In RStudio, set the working directory to the repository root, edit the
configuration if required, and source this plot script.

## Script guide

| Script | What it does |
| --- | --- |
| `model_revision.R` | Fits cumulative-link models, performs participant-grouped CV, calculates prediction metrics, and bootstraps Kumaraswamy parameters/weights and paired metric differences. Supports primary and screened-out cohorts. |
| `screening_compare.R` | Matches predictions by dataset, model, record, participant, outcome, and fold, then compares training with and without flagged participants. |
| `table3_paired.R` | Calculates 12 participant-paired comparisons on complete weeks, confidence intervals, Cohen's dz, and Holm-adjusted p values. |
| `transition_audit.R` | Reconstructs the 25 objective quality profiles, counts day-to-day transitions, and summarizes matched odd/even held-out errors. |
| `longitudinal_missingness.R` | Summarizes weekly/daily availability and estimates exploratory within-participant linear trends with participant-bootstrap intervals. |
| `plot_cv_weekly_mos.R` | Draws observed weekly MOS and mean out-of-fold predictions for all rated and complete weeks. |
| `run_all.R` | Executes the workflow in dependency order and records the environment/settings. |
| `check_results.R` | Compares a full rerun with the supplied reference tables. |

## Input data and analysis populations

The input table contains 829 rows and the following analysis fields:

| Field | Meaning |
| --- | --- |
| `external_id` | Participant identifier. |
| `week` | Study week, 1–25. |
| `q2` | Retrospective rating on the 1–5 response scale, subject to the source-log correction below. |
| `mo`, `tu`, `we`, `th`, `fr`, `sa` | Objective P.1204.4 quality for each daily video; missing values denote unavailable daily records. |
| `n` | Source count including daily video records and the weekly response; a complete week has `n = 7`, representing six daily videos plus the response. |

Source-log auditing identified one response for participant `343ab`, week 14,
as belonging to a daily stimulus rather than the dedicated weekly rating.
The analysis scripts account for this correction; the primary modeling script
saves the excluded record for inspection. Do not delete it from the input CSV:
the scripts verify the expected anomaly before applying the correction.

The corrected primary sample contains 828 ratings from 38 participants,
including 552 complete weeks. The six participants flagged by exploratory
screening remain in the primary analysis. The sensitivity analysis excludes
`f1308`, `bb8d0`, `343ab`, `9d7d1`, `3fe63`, and `57ea9`, leaving 709 ratings
from 32 participants, including 492 complete weeks.

The scripts contain checks specific to this study's IDs, counts, and schedule.
Adapt these checks and model assumptions before applying the code to another dataset.

## Models and outputs

Seven models are compared on the same 552 complete-week outcomes:
`all` (six daily predictors), `agg` (three day-pair means), `odd` (days 1, 3, 5),
`even` (days 2, 4, 6), `mean`, `min`, and `opt_k` (Kumaraswamy-weighted quality).
Only `mean`, `min`, and `opt_k` are compared on all 828 rated weeks.
Mean/min use observed days; opt_k renormalizes their weights.

| Output | Interpretation |
| --- | --- |
| `cv_metrics.csv` | RMSE, MAE, Pearson correlation (PLCC), and Spearman correlation (SRCC) for individual out-of-fold ratings. |
| `oof_predictions.csv` | Predictions made by models trained without the participant being evaluated. |
| `cv_week_mos.csv`, `cv_week_metrics.csv` | Predictions and observed ratings averaged within study week; a separate target from individual-rating accuracy. |
| `full_sample_fits.csv` | Descriptive in-sample AIC and Kumaraswamy shape estimates. AIC includes the two optimized shape parameters. |
| `bootstrap_*.csv` | Participant-bootstrap refits of shape parameters and six base day weights. |
| `parameter_weights_ci_*.csv` | Full-sample estimates and percentile 95% intervals for parameters/weights, with successful and attempted resample counts. |
| `paired_metric_differences.csv` | opt_k minus comparator; negative RMSE/MAE differences favor opt_k. Intervals resample participants with out-of-fold predictions held fixed. |
| `screening_matched_cv.csv` | Screened-out training minus full-cohort training, evaluated on the same sensitivity-sample outcomes. |
| `table3_paired_results.csv` | Paired tests; `p_holm_12` adjusts all 12 comparisons. The separate eight-comparison adjustment is also supplied. |
| `exploratory_within_participant_trends.csv` | Participant fixed-intercept linear summaries and participant-bootstrap intervals, including adjustment for quality and ordered content blocks. |

Diagnostic CSV files are generated when applicable. `fit_errors.csv` records
model failures; the modeling script then stops rather than issuing a partial
CV summary. Use a new output directory for an independent run so old diagnostic
files cannot be mistaken for new errors.

Compare models within the same analysis population. Cross-validation is grouped
by participant, not by week. It evaluates new participants under the observed
stimulus schedule, not new schedules or future weeks. Content order and
transition positions limit causal interpretations of fitted weights and trends.
The transition threshold of 0.5 P.1204.4 units is descriptive. Missing-day
renormalization does not establish that missingness is random.
