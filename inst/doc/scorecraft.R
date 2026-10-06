## ----include = FALSE----------------------------------------------------------
has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
has_db     <- requireNamespace("DBI", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = has_glmnet)
old_dt <- data.table::setDTthreads(2)

## ----presets------------------------------------------------------------------
library(scorecraft)
library(data.table)
scr_presets()[, c("preset", "target_max", "min_votes", "corr_cutoff", "iv_min")]

## ----config-------------------------------------------------------------------
cfg <- scr_config("moderate", objective = "risk", verbose = FALSE, nthread = 1,
                  use_glmnet = TRUE, use_ranger = FALSE, use_lightgbm = FALSE,
                  xgb_rounds = 60, n_boot = 20)
cfg

## ----select-------------------------------------------------------------------
res <- scr_select(scr_demo, "default", config = cfg, drop = c("id", "churn"),
                  date_col = "ref_date")
res

## ----funnel-------------------------------------------------------------------
table(scr_funnel(res, cols = "all")$exit_stage)
head(scr_funnel(res, only_selected = TRUE)[, .(feature, total_iv, iv_holdout, ks, psi, psi_flag_adjusted)])

## ----scorecard----------------------------------------------------------------
sc <- scr_scorecard(res)
sc
all(sc$sign_check$coef > 0)
scr_score_metrics(sc)[, .(sample, auc, auc_lo, auc_hi, ks, gini)]

## ----cutoff-------------------------------------------------------------------
ct <- scr_cutoff(sc, n_cuts = 8)
ct
ct$table[cut == ct$cuts[4], .(sample, cut, pct_safe, event_rate_safe, events_avoided_pct, ks_at_cut)]

## ----strategy-----------------------------------------------------------------
st <- scr_strategy(sc, revenue_good = 1080, loss_bad = 4500)
st
st$table[, .(band, event_rate = round(event_rate, 4), decision, cum_pct = round(cum_pct, 3), cum_profit)]

## ----reject-------------------------------------------------------------------
set.seed(11)
train <- scr_demo[res$split$train_idx, ]
old_score <- scr_apply(sc, train)$score + stats::rnorm(nrow(train), sd = 20)
declined <- old_score < stats::quantile(old_score, 0.30) | stats::runif(nrow(train)) < 0.05
ttd <- rbind(scr_demo[res$split$holdout_idx, ], train[declined, ])
booked <- rep(c(TRUE, FALSE), c(length(res$split$holdout_idx), sum(declined)))
rj <- scr_reject(sc, population = ttd, accepted = booked)
rj
rj$coverage[, .(band, n_dev, n_unknown, coverage = round(coverage, 2), coverage_flag)]

## ----apply--------------------------------------------------------------------
new <- head(scr_demo, 5)
scr_apply(sc, new, what = "all")[, .(prob = round(prob, 4), score = round(score, 2), score_points,
                                     vl_score_01_woe = round(vl_score_01_woe, 3), vl_score_01_points)]

## ----reasons------------------------------------------------------------------
scr_reasons(sc, new, k = 2)

## ----sql----------------------------------------------------------------------
sql <- scr_sql(sc, table = "prd.customers", dialect = "databricks",
               what = "all", keep_columns = c("id", "ref_date"))
length(sql)
cat(grep("^-- (Scale|score =)", sql, value = TRUE), sep = "\n")
i <- max(which(sql == "SELECT"))
cat(sql[i:(i + 6)], sep = "\n")

## ----duckdb, eval = has_glmnet && has_db && requireNamespace("duckdb", quietly = TRUE), message = FALSE----
con <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
DBI::dbWriteTable(con, "scr_demo", scr_demo)
got <- DBI::dbGetQuery(con, paste(scr_sql(sc, table = "scr_demo", dialect = "duckdb", what = "all",
                                          keep_columns = "id"), collapse = "\n"))
DBI::dbDisconnect(con, shutdown = TRUE)
got <- got[order(got$id), ]
exp <- scr_apply(sc, scr_demo, what = "all")
all.equal(got$score, exp$score)
identical(as.numeric(got$score_points), as.numeric(exp$score_points))
all(vapply(sc$features, function(f) isTRUE(all.equal(got[[paste0(f, "_points")]], exp[[paste0(f, "_points")]])) &&
             isTRUE(all.equal(got[[paste0(f, "_woe")]], exp[[paste0(f, "_woe")]])), logical(1)))

## ----monitor------------------------------------------------------------------
mo <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 20)
mo

## ----monitor-csi--------------------------------------------------------------
mo$csi[variable == "vl_late", .(period, csi = round(csi, 4), flag_fixed, flag_adjusted,
                                points_shift = round(points_shift, 2))]

## ----plan---------------------------------------------------------------------
plan <- sc$monitoring_plan
plan[plan$item != "threshold_source", ]
plan$value[plan$item == "min_events_per_period"] <- "90"
mo2 <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 20, plan = plan)
data.table(period = mo$vintage$period, events = mo$vintage$events,
           status_default = mo$vintage$status, status_edited = mo2$vintage$status)

## ----export, eval = has_glmnet && requireNamespace("openxlsx", quietly = TRUE), message = FALSE----
out <- file.path(tempdir(), "scorecraft-get-started")
files_res <- scr_export(res, out, stamp = FALSE)$files
files_sc  <- scr_export(sc, out, stamp = FALSE, cutoff = ct, strategy = st, reject = rj,
                        monitor = mo)$files
basename(unlist(c(files_res, files_sc)))
mo3 <- scr_monitor(sc, scr_demo, date_col = "ref_date", target = "default", n_boot = 20,
                   plan = files_sc$strategy)
identical(mo3$psi$flag_adjusted, mo$psi$flag_adjusted)

## ----db, eval = has_glmnet && has_db && requireNamespace("RSQLite", quietly = TRUE)----
con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
d <- scr_demo
d$ref_date <- as.character(d$ref_date)
DBI::dbWriteTable(con, "dtm", d)
nrow(scr_fetch(con, "dtm", sample_frac = 0.5, seed = 42))
nrow(scr_fetch(con, "dtm", max_rows = 1000))

## ----run, eval = has_glmnet && has_db && requireNamespace("RSQLite", quietly = TRUE)----
rs <- scr_run(con, "dtm", targets = c("default", "churn"), config = cfg, drop = "id",
              date_col = "ref_date")
rs
c(split = rs$default$split$method, cutoff = rs$default$split$cutoff,
  table = rs$default$config$sql_table)
DBI::dbDisconnect(con)

## ----include = FALSE, eval = TRUE---------------------------------------------
data.table::setDTthreads(old_dt)

