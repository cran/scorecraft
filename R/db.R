# ============================================================================ #
# db.R - database layer (ODBC/DSN, or any DBI driver)
# ============================================================================ #

#' Connect to a database (ODBC DSN or any DBI driver)
#'
#' With `dsn`, a thin wrapper around [DBI::dbConnect()] over [odbc::odbc()]
#' with one deliberate choice: `bigint = "numeric"`. Under the \pkg{odbc}
#' default a BIGINT column arrives as `integer64`, and `is.numeric()` of an
#' `integer64` is `FALSE`: typing would treat the column as a categorical of
#' very high cardinality. With `driver`, any DBI driver object is accepted
#' (e.g. `RSQLite::SQLite()`, `duckdb::duckdb()`), which is how the database
#' path is tested without a DSN.
#'
#' @param dsn Name of the DSN configured on the system. Ignored when
#'   `driver` is given.
#' @param driver Optional DBI driver object, used instead of ODBC.
#' @param timeout Connection timeout in seconds (ODBC only).
#' @param ... Extra arguments passed on to [DBI::dbConnect()] (for example
#'   `dbname = ":memory:"` for SQLite).
#'
#' @return A DBI connection. Close it with [DBI::dbDisconnect()].
#'
#' @family database
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE)
#' con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
#' d <- scr_demo; d$ref_date <- as.character(d$ref_date)   # SQLite has no Date type
#' DBI::dbWriteTable(con, "dtm", d)
#' dt <- scr_fetch(con, "dtm", sample_frac = 0.5, seed = 42)
#' nrow(dt)
#' DBI::dbDisconnect(con)
#' @export
scr_connect <- function(dsn = NULL, driver = NULL, timeout = 20, ...) {
  if (!requireNamespace("DBI", quietly = TRUE)) stop("scr_connect() needs the 'DBI' package.", call. = FALSE)
  if (!is.null(driver)) return(DBI::dbConnect(driver, ...))
  if (is.null(dsn)) stop("scr_connect(): give a `dsn` or a DBI `driver`.", call. = FALSE)
  if (!requireNamespace("odbc", quietly = TRUE)) stop("scr_connect() with a DSN needs the 'odbc' package.", call. = FALSE)
  DBI::dbConnect(odbc::odbc(), dsn = dsn, timeout = timeout, bigint = "numeric", ...)
}

#' Fetch a table with reproducible server-side sampling
#'
#' Sampling happens on the server, not in R: pulling a million rows to
#' discard ninety per cent of them pays the network cost twice. `max_rows` is
#' a memory guard: when it binds, the requested fraction is reduced on the
#' server and the reduction is reported. The random expression follows the
#' connection class (`rand(seed)` on Spark/Databricks/MySQL, `random()` on
#' PostgreSQL/DuckDB, an integer-modulo expression on SQLite, which has no
#' seedable `random()`); pass `sample_expr` to override it.
#'
#' @param con A DBI connection, from [scr_connect()].
#' @param table Qualified table name.
#' @param sample_frac Fraction of rows to fetch, in (0, 1].
#' @param seed Seed of the server-side random function, where supported.
#' @param max_rows Row cap. `NULL` switches it off.
#' @param sample_expr Optional SQL expression yielding a uniform number in
#'   `[0, 1)`, used as `WHERE <sample_expr> <= sample_frac`.
#' @param verbose `TRUE`/`FALSE` to echo (or not) the query for this call;
#'   `NULL` follows [scr_verbose()].
#'
#' @return A `data.table` with the fetched table.
#'
#' @family database
#' @examplesIf requireNamespace("RSQLite", quietly = TRUE)
#' con <- scr_connect(driver = RSQLite::SQLite(), dbname = ":memory:")
#' d <- scr_demo; d$ref_date <- as.character(d$ref_date)
#' DBI::dbWriteTable(con, "dtm", d)
#' nrow(scr_fetch(con, "dtm", sample_frac = 0.5, seed = 42))
#' nrow(scr_fetch(con, "dtm", max_rows = 1000))
#' DBI::dbDisconnect(con)
#' @export
scr_fetch <- function(con, table, sample_frac = 1.0, seed = NULL, max_rows = NULL, sample_expr = NULL,
                      verbose = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) stop("scr_fetch() needs the 'DBI' package.", call. = FALSE)
  if (!is.null(verbose)) {
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("scr_fetch(): `verbose` must be TRUE, FALSE or NULL.", call. = FALSE)
    old <- scr_verbose(verbose); on.exit(scr_verbose(old), add = TRUE)
  }
  .scr_num1(sample_frac, "sample_frac", lower = 0, upper = 1, open_lower = TRUE)
  if (!is.null(max_rows) && is.finite(max_rows)) {
    n_tab <- tryCatch(as.numeric(DBI::dbGetQuery(con, DBI::SQL(sprintf(
      "select count(1) as n from %s", table)))$n[1]), error = function(e) NA_real_)
    if (is.finite(n_tab) && n_tab * sample_frac > max_rows) {
      new <- max_rows / n_tab
      msg("  cap of %s rows: fraction reduced from %.4f to %.4f (table has %s)",
          n_fmt(max_rows), sample_frac, new, n_fmt(n_tab))
      sample_frac <- new
    }
  }
  query <- if (sample_frac > 0 && sample_frac < 1) {
    expr <- sample_expr %||% .sample_expr(con, seed)
    # full precision: under a tight `max_rows` on a large table the fraction
    # can be below 1e-6, which "%.6f" would print as 0 (and fetch nothing)
    sprintf("select * from %s where %s <= %s", table, expr, .sql_num(sample_frac))
  } else sprintf("select * from %s", table)
  msg("SQL: %s", query)
  dt <- DBI::dbGetQuery(con, statement = DBI::SQL(query))
  data.table::setDT(dt)
  dt[]
}

#' Uniform [0, 1) expression for the connection's dialect
#' @keywords internal
#' @noRd
.sample_expr <- function(con, seed = NULL) {
  cls <- class(con)[1]
  if (grepl("SQLite", cls, ignore.case = TRUE)) {
    # SQLite random() is a signed 64-bit integer; no seed is available
    "((abs(random()) % 1000000) / 1000000.0)"
  } else if (grepl("duckdb|Pq|Postgres|Redshift", cls, ignore.case = TRUE)) {
    "random()"
  } else if (is.null(seed)) {
    "rand()"
  } else {
    sprintf("rand(%d)", as.integer(seed))
  }
}
