# ============================================================================ #
# sql.R - production SQL engine
# ============================================================================ #
# Two blocks, in this order: a pre-processing CTE (training median,
# special-population flags, categorical COALESCE) and the WOE/BIN
# transformation emitted by the binning engine from the authoritative cut
# points. For the scorecard, a third block composes the score from the WOE
# columns and the whole points from the bin index. Covers exactly the
# shortlist.
# ============================================================================ #

#' Production SQL
#'
#' Code ready to run in the database, covering **exactly** the approved
#' variables (or those of the scorecard), in blocks in this order:
#'
#' \enumerate{
#'   \item CTE `base_scr`: reproduces the Stage 1 pre-processing -
#'     imputation of missing and sentinel by the **training** median,
#'     special-population flags, `COALESCE` of the categorical missing.
#'   \item The WOE/BIN transformation, emitted by
#'     [OptimalBinningWoE::obwoe_sql()] from the authoritative cut points
#'     with full precision.
#'   \item (Scorecard) CTE `woe_scr` with WOE and bin index, followed by the
#'     final `SELECT` with `score` (exact, `a + b * logit`), `<f>_points`
#'     per variable and `score_points` (whole points).
#' }
#'
#' The order matters: without the first block, the WOE would be applied to
#' data different from what was binned. Column names are quoted with the
#' dialect's delimiters only when they are not plain identifiers or are
#' reserved words, the same rule [OptimalBinningWoE::obwoe_sql()] applies,
#' so every block names a column the same way. A row whose value falls in
#' no fitted bin (a category never seen on train) takes a WOE of 0 and the
#' points of a WOE of 0, in the SQL as in [scr_apply()]. The score computed by the SQL
#' matches [scr_apply()] numerically, by an automated test that runs both
#' paths.
#'
#' @section IRB models:
#'
#' `scr_pd` wraps the scorecard SQL in a common table expression and adds a
#' `CASE` on the score cut points that yields `grade` and `pd_final`.
#' `scr_lgd` chains the driver bins of both stages, the logits, the pool
#' `CASE` and the floored result. `scr_ead` computes the utilisation and the
#' undrawn amount, assigns the pool from the frozen cut points and applies
#' the greatest of the model, the drawn amount and the standardised floor.
#' `scr_capital` carries the constants of every pool (PD, LGD, `k`, risk
#' weight) in a `pool_params` table joined on segment and grade, so no
#' normal quantile is evaluated at run time; `level` chooses the exposure
#' or the portfolio output.
#'
#' @param x An object from [scr_select()], [scr_scorecard()], [scr_pd()],
#'   [scr_lgd()], [scr_ead()] or [scr_capital()].
#' @param table Source table name, written verbatim (it may be qualified,
#'   `schema.table`, and is never quoted: pass only a trusted name). `NULL`
#'   uses `config$sql_table`.
#' @param dialect Dialect (`"ansi"`, `"databricks"`, `"spark"`, `"hive"`,
#'   `"mysql"`, `"mariadb"`, `"sqlserver"`, `"bigquery"`, `"postgres"`,
#'   `"oracle"`, `"snowflake"`, `"redshift"`, `"duckdb"`, `"sqlite"`).
#'   `NULL` uses `config$sql_dialect`.
#' @param file Path to write to. `NULL` (default) returns the lines.
#' @param ... Passed on to the methods.
#' @param output For `scr_result`: `"woe"`, `"bin"` or `"both"`. `NULL` uses
#'   `config$sql_output`.
#' @param what For `scr_scorecard`: `"score"` (default: the three blocks, with
#'   the points per variable, the exact score and the whole-points score),
#'   `"woe"` (the WOE/BIN SQL of the scorecard variables only) or `"all"`
#'   (the three blocks plus, for every variable, its bin label, WOE and
#'   points side by side: the deployment layout that reports the band of
#'   each variable next to the score).
#' @param keep_columns For `scr_scorecard`: key columns carried untransformed
#'   into the output (for example the customer identifier and the reference
#'   date); `NULL` uses `config$sql_keep_columns`.
#'
#' @return A character vector with the SQL (invisibly, when `file` is given).
#'
#' @family production
#' @examples
#' cfg <- scr_config(verbose = FALSE, nthread = 1, use_ranger = FALSE,
#'                   xgb_rounds = 60, n_boot = 20)
#' res <- scr_select(scr_demo, "default", config = cfg, drop = "id",
#'                   date_col = "ref_date")
#' cat(head(scr_sql(res, table = "prd.customers", dialect = "databricks"), 20), sep = "\n")
#' sc <- scr_scorecard(res)
#' cat(tail(scr_sql(sc), 12), sep = "\n")
#' @order 1
#' @export
scr_sql <- function(x, table = NULL, dialect = NULL, file = NULL, ...) UseMethod("scr_sql")

#' @rdname scr_sql
#' @export
scr_sql.scr_result <- function(x, table = NULL, dialect = NULL, file = NULL, output = NULL, ...) {
  .scr_local_verbose(x)
  cfg <- x$config
  if (!is.null(table)) cfg$sql_table <- table
  if (!is.null(dialect)) cfg$sql_dialect <- dialect
  if (!is.null(output)) cfg$sql_output <- output
  out <- if (is.null(table) && is.null(dialect) && is.null(output)) x$sql else
    build_sql_woe(x$fit, x$triage$ledger, scr_selected(x), cfg, x$target, provenance = .provenance_line(x$lab))
  .sql_out(out, file)
}

#' @rdname scr_sql
#' @export
scr_sql.scr_scorecard <- function(x, table = NULL, dialect = NULL, file = NULL, what = c("score", "woe", "all"),
                                  keep_columns = NULL, ...) {
  .scr_local_verbose(x)
  what <- match.arg(what)
  cfg <- x$config
  if (!is.null(table)) cfg$sql_table <- table
  if (!is.null(dialect)) cfg$sql_dialect <- dialect
  if (!is.null(keep_columns)) {
    if (!is.character(keep_columns) || anyNA(keep_columns) || any(!nzchar(keep_columns))) {
      stop("scr_sql(): `keep_columns` must be a character vector of column names.", call. = FALSE)
    }
    cfg$sql_keep_columns <- keep_columns
  }
  cached <- is.null(table) && is.null(dialect) && is.null(keep_columns)
  out <- switch(what,
    woe   = build_sql_woe(x$fit, x$ledger, x$features, cfg, x$target, provenance = .provenance_line(x$lab)),
    score = if (cached) x$sql else { y <- x; y$config <- cfg; build_sql_score(y) },
    all   = { y <- x; y$config <- cfg; build_sql_score(y, what = "all") })
  .sql_out(out, file)
}

#' @keywords internal
#' @noRd
.sql_out <- function(out, file) {
  if (!is.null(file)) {
    writeLines(out, con = file)
    msg("SQL written to %s", file)
    return(invisible(out))
  }
  out
}

#' Quote a column identifier exactly as OptimalBinningWoE::obwoe_sql() does
#'
#' The WOE block emitted by the engine (`quote_identifiers = "auto"`) reads
#' the columns of `base_scr` and names its outputs with this rule: a name is
#' quoted, with the dialect's delimiters, only when it is not a plain
#' `[A-Za-z_][A-Za-z0-9_]*` token or is a reserved word. The blocks written
#' here must follow the very same rule, or a variable called `order`, `user`
#' or `my var` is written under one name and read under another (and on
#' case-folding engines quoting changes the name itself).
#' @keywords internal
#' @noRd
.sql_q <- function(x, dialect) {
  reserved <- c(
    "all", "alter", "and", "any", "as", "asc", "between", "by", "case", "cast",
    "check", "column", "create", "cross", "current", "current_date",
    "current_time", "current_timestamp", "default", "delete", "desc",
    "distinct", "drop", "else", "end", "except", "exists", "false", "for",
    "foreign", "from", "full", "grant", "group", "having", "in", "index",
    "inner", "insert", "intersect", "into", "is", "join", "left", "like",
    "limit", "not", "null", "on", "or", "order", "outer", "primary",
    "references", "right", "select", "set", "some", "table", "then", "to",
    "true", "union", "unique", "update", "user", "using", "values", "view",
    "when", "where", "with")
  q <- if (dialect %in% c("mysql", "mariadb", "spark", "hive", "databricks", "bigquery")) c("`", "`")
       else if (identical(dialect, "sqlserver")) c("[", "]") else c("\"", "\"")
  one <- function(part) {
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", part) && !tolower(part) %in% reserved) return(part)
    paste0(q[1], gsub(q[2], paste0(q[2], q[2]), part, fixed = TRUE), q[2])
  }
  # dotted names are quoted part by part, as the engine does
  vapply(as.character(x), function(nm) {
    parts <- strsplit(nm, ".", fixed = TRUE)[[1L]]
    if (!length(parts)) parts <- nm
    paste(vapply(parts, one, character(1)), collapse = ".")
  }, character(1), USE.NAMES = FALSE)
}

#' Block 1: lines of the SELECT of the pre-processing CTE
#' @keywords internal
#' @noRd
.sql_preprocess_lines <- function(ledger, features, cfg) {
  sources <- unique(vapply(features, function(f) {
    r <- ledger[output == f & kind == "num_flag", source]
    if (length(r)) r[1] else f
  }, character(1)))
  dl <- cfg$sql_dialect
  lines <- character()
  for (s in cfg$sql_keep_columns) lines <- c(lines, sprintf("    %s", .sql_q(s, dl)))
  sp <- cfg$special_values
  for (s in sources) {
    imp <- ledger[kind == "num_impute" & source == s]
    flg <- ledger[kind == "num_flag"   & source == s]
    coa <- ledger[kind == "cat_coalesce" & source == s]
    qs <- .sql_q(s, dl)
    if (s %in% features) {
      if (nrow(imp)) {
        cond <- sprintf("%s IS NULL", qs)
        if (length(sp)) cond <- sprintf("%s OR %s IN (%s)", cond, qs, paste(vapply(sp, .sql_num, character(1)), collapse = ", "))
        lines <- c(lines, sprintf("    CASE WHEN %s THEN %s ELSE %s END AS %s", cond, .sql_num(imp$impute_value[1]), qs, qs))
      } else if (nrow(coa)) {
        lines <- c(lines, sprintf("    COALESCE(%s, %s) AS %s", qs, .sql_str("MISSING"), qs))
      } else {
        lines <- c(lines, sprintf("    %s", qs))
      }
    }
    if (nrow(flg) && flg$output[1] %in% features) {
      whens <- sprintf("WHEN %s IS NULL THEN %s", qs, .sql_str("MISSING"))
      for (v in sp) whens <- c(whens, sprintf("WHEN %s = %s THEN %s", qs, .sql_num(v), .sql_str(paste0("S", v))))
      lines <- c(lines, sprintf("    CASE %s ELSE %s END AS %s", paste(whens, collapse = " "), .sql_str("REGULAR"),
                                .sql_q(flg$output[1], dl)))
    }
  }
  lines
}

#' Selection SQL: pre-processing CTE + WOE/BIN from the engine
#' @keywords internal
#' @noRd
build_sql_woe <- function(fit, ledger, features, cfg, target = NULL, provenance = NULL) {
  if (!length(features)) return("-- no approved variable: nothing to generate")
  cte <- "base_scr"
  lines <- .sql_preprocess_lines(ledger, features, cfg)
  woe_sql <- OptimalBinningWoE::obwoe_sql(
    obj = fit, table = cte, features = features, output = cfg$sql_output, style = "select",
    dialect = cfg$sql_dialect, keep_columns = if (length(cfg$sql_keep_columns)) cfg$sql_keep_columns else NULL,
    digits = NULL, comment = TRUE, bin_separator = cfg$bin_separator)
  c("-- =============================================================",
    sprintf("-- scorecraft | target: %s | %d approved variables | dialect: %s", .sql_cmt(target %||% "?"), length(features), cfg$sql_dialect),
    sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "-- Block 1 (CTE base_scr): Stage 1 pre-processing - imputation of missing",
    "--   and sentinel values by the TRAINING median, special-population flags.",
    "-- Block 2: WOE/BIN transformation emitted by OptimalBinningWoE::obwoe_sql().",
    if (!is.null(provenance)) paste0("-- ", .sql_cmt(provenance)),
    "-- =============================================================", "",
    sprintf("WITH %s AS (", cte), "  SELECT", paste(lines, collapse = ",\n"),
    sprintf("  FROM %s", cfg$sql_table), ")", "", as.character(woe_sql)) |> .sql_lines()
}

#' Text safe inside a `--` comment: a line break in a target or variable
#' name would otherwise end the comment and turn the rest into SQL
#' @keywords internal
#' @noRd
.sql_cmt <- function(x) gsub("[\r\n]+", " ", as.character(x))

#' One element per line, so head()/tail()/grep() work line by line
#' @keywords internal
#' @noRd
.sql_lines <- function(x) unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))

#' Extract the column expressions of a SELECT emitted by obwoe_sql
#' @keywords internal
#' @noRd
.sql_select_exprs <- function(sql) {
  l <- unlist(strsplit(as.character(sql), "\n", fixed = TRUE))
  i0 <- which(trimws(l) == "SELECT")[1]
  i1 <- max(which(grepl("^FROM ", trimws(l))))
  body <- l[(i0 + 1L):(i1 - 1L)]
  body <- body[nzchar(trimws(body))]
  body[length(body)] <- sub(",\\s*$", "", body[length(body)])
  body
}

#' Scorecard SQL: pre-processing CTE + WOE/index CTE + score
#' @keywords internal
#' @noRd
build_sql_score <- function(sc, what = c("score", "all")) {
  what <- match.arg(what)
  cfg <- sc$config
  feats <- sc$features
  dl <- cfg$sql_dialect
  keep <- .sql_q(cfg$sql_keep_columns, dl)
  # the names of the WOE, index, bin and points columns, quoted as the engine quotes them
  qn <- function(f, suffix) .sql_q(paste0(f, suffix), dl)
  lines <- .sql_preprocess_lines(sc$ledger, feats, cfg)
  woe_x <- OptimalBinningWoE::obwoe_sql(obj = sc$fit, table = "base_scr", features = feats, output = "woe",
                                        style = "select", dialect = cfg$sql_dialect, digits = NULL,
                                        comment = FALSE, bin_separator = cfg$bin_separator)
  idx_x <- OptimalBinningWoE::obwoe_sql(obj = sc$fit, table = "base_scr", features = feats, output = "index",
                                        style = "select", dialect = cfg$sql_dialect, digits = NULL,
                                        comment = FALSE, bin_separator = cfg$bin_separator)
  exprs <- c(if (length(keep)) paste0(keep, ","), .sql_select_exprs(woe_x))
  exprs[length(exprs)] <- paste0(exprs[length(exprs)], ",")
  exprs <- c(exprs, .sql_select_exprs(idx_x))
  if (identical(what, "all")) {
    # the bin labels, for a deployment that reports the band of every variable
    bin_x <- OptimalBinningWoE::obwoe_sql(obj = sc$fit, table = "base_scr", features = feats, output = "bin",
                                          style = "select", dialect = cfg$sql_dialect, digits = NULL,
                                          comment = FALSE, bin_separator = cfg$bin_separator)
    exprs[length(exprs)] <- paste0(exprs[length(exprs)], ",")
    exprs <- c(exprs, .sql_select_exprs(bin_x))
  }

  al <- sc$alignment
  base_raw <- al$a + al$b * unname(sc$coef["(Intercept)"])
  score_terms <- c(.sql_num(base_raw), vapply(feats, function(f)
    sprintf("%s * %s", .sql_num(al$b * unname(sc$coef[f])), qn(f, "_woe")), character(1)))
  # ELSE: a row in no fitted bin (NULL index) has WOE 0 and takes the points
  # of WOE 0, exactly as scr_apply()
  unb <- .sql_num(.sc_unbinned_points(sc))
  pts_cases <- vapply(feats, function(f) {
    p <- sc$points[variable == f]
    sprintf("    CASE %s %s ELSE %s END AS %s", qn(f, "_idx"),
            paste(sprintf("WHEN %d THEN %s", as.integer(p$bin_id), .sql_num(p$points)), collapse = " "), unb,
            qn(f, "_points"))
  }, character(1))
  # the subselect computes the points columns once; the outer SELECT exposes
  # them by name and sums them into score_points
  per_feature <- if (identical(what, "all")) {
    # bin label, WOE and points of every variable, side by side
    paste0(unlist(lapply(feats, function(f) sprintf("    %s", qn(f, c("_bin", "_woe", "_points"))))), ",")
  } else {
    paste0(vapply(feats, function(f) sprintf("    %s", qn(f, "_points")), character(1)), ",")
  }
  final <- c(
    "SELECT",
    if (length(keep)) sprintf("    %s,", keep),
    sprintf("    %s AS score,", paste(score_terms, collapse = "\n      + ")),
    per_feature,
    sprintf("    %s AS score_points", paste(c(.sql_num(sc$base_points), qn(feats, "_points")), collapse = " + ")),
    "FROM (", "  SELECT", "    *,", paste0("  ", pts_cases, collapse = ",\n"), "  FROM woe_scr", ") pts;")

  c("-- =============================================================",
    sprintf("-- scorecraft | scorecard of target: %s | %d variables | dialect: %s", .sql_cmt(sc$target), length(feats), cfg$sql_dialect),
    sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("-- Scale: %s points at odds %s:1 (%s), PDO %s | %s", format(sc$scale$base_score), format(sc$scale$base_odds),
            sc$odds_orientation, format(sc$scale$pdo), sc$direction),
    sprintf("-- score = %s + %s * logit | base_points = %s", .sql_num(al$a), .sql_num(al$b), format(sc$base_points)),
    "-- Block 1 (CTE base_scr): pre-processing frozen on train.",
    if (identical(what, "all")) "-- Block 2 (CTE woe_scr): WOE, bin index and bin label, emitted by OptimalBinningWoE::obwoe_sql()."
    else "-- Block 2 (CTE woe_scr): WOE and bin index, emitted by OptimalBinningWoE::obwoe_sql().",
    if (identical(what, "all")) "-- Block 3: bin label, WOE and points of every variable, exact score (from the WOE) and whole points (from the bin index)."
    else "-- Block 3: exact score (from the WOE) and whole points (from the bin index).",
    if (!is.null(sc$lab)) paste0("-- ", .sql_cmt(.provenance_line(sc$lab))),
    "-- =============================================================", "",
    "WITH base_scr AS (", "  SELECT", paste(lines, collapse = ",\n"), sprintf("  FROM %s", cfg$sql_table), "),",
    "woe_scr AS (", "  SELECT", exprs, "  FROM base_scr", ")", "", final) |> .sql_lines()
}
