# ============================================================================ #
# utils.R - internal utilities: verbosity, parallelism, formatting
# ============================================================================ #
# Nothing here knows about the pipeline. Pure functions, reusable and testable
# in isolation.
# ============================================================================ #

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# -- Verbosity -------------------------------------------------------------- #

.scr_env <- new.env(parent = emptyenv())
.scr_env$verbose <- TRUE

#' Switch progress messages on or off
#'
#' Large tables take tens of minutes, and the pipeline reports every stage as
#' it runs: stage name, input and output counts, elapsed time. Messages go
#' through [message()] as single lines, with no progress bar that redraws
#' itself, so that a scheduled job (`Rscript` in batch) produces a readable
#' log. [suppressMessages()] works too; this function exists to switch them
#' off persistently, without wrapping every call. The `verbose` key of
#' [scr_config()] has the same effect per run.
#'
#' @param on `TRUE` to switch on, `FALSE` to switch off, `NULL` (default) to
#'   query the current state only.
#'
#' @return The verbosity state in force *before* the call, invisibly, so that
#'   `old <- scr_verbose(FALSE); ...; scr_verbose(old)` restores it.
#'
#' @family configuration
#' @examples
#' old <- scr_verbose(FALSE)   # silence, keeping the previous state
#' scr_verbose(old)            # restore
#' scr_verbose()               # query
#' @export
scr_verbose <- function(on = NULL) {
  old <- .scr_env$verbose
  if (!is.null(on)) .scr_env$verbose <- isTRUE(on)
  invisible(old)
}

#' Follow the `verbose` key of an object's configuration until the caller exits
#'
#' For the methods that act on a fitted object (export, SQL to a file, the
#' classing lab): the object's `config$verbose` rules, as it does in the
#' function that fitted it. Objects without a configuration leave the
#' session setting of [scr_verbose()] as it is.
#' @keywords internal
#' @noRd
.scr_local_verbose <- function(x, frame = parent.frame()) {
  cfg <- if (is.list(x)) x$config %||% (if (is.list(x$result)) x$result$config) else NULL
  if (is.null(cfg) || is.null(cfg$verbose)) return(invisible(NULL))
  old <- scr_verbose(isTRUE(cfg$verbose))
  do.call(base::on.exit, list(call("scr_verbose", old), add = TRUE), envir = frame)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
msg <- function(fmt, ...) {
  if (isTRUE(.scr_env$verbose)) {
    message(if (length(list(...))) sprintf(fmt, ...) else fmt)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
msg_stage <- function(n, title) {
  msg("\n===== STAGE %s: %s =====", n, title)
}

#' Time an expression, reporting the elapsed time at the end
#' @keywords internal
#' @noRd
time_it <- function(label, expr) {
  t0  <- Sys.time()
  res <- eval.parent(substitute(expr))
  msg("  %s - ok (%.2fs)", label, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  res
}

# -- Parallelism by column -------------------------------------------------- #

#' Parallel lapply with a serial fallback
#'
#' None of the binning algorithms runs in parallel, but columns are
#' embarrassingly parallel. The backend is chosen by
#' `getOption("scorecraft.parallel")`: `"fork"`
#' ([parallel::mclapply()], the default on unix), `"psock"` (a
#' [parallel::makeCluster()] cluster, the default on Windows and the path
#' exercised by the tests on every platform, since PSOCK workers share no
#' memory and therefore prove that every closure serialises), or `"serial"`.
#' Serial when `nthread <= 1` or there are fewer than two items.
#'
#' Every worker returns a sealed envelope (value, error message, warnings),
#' so that in the parent:
#' * an error in any worker is re-thrown with the failing item, instead of
#'   becoming a silent `try-error` inside the list;
#' * a `NULL` can only mean that the worker died (out of memory, a signal)
#'   and is reported as such, instead of surfacing later as a subscript
#'   error; a worker legitimately returning `NULL` still does;
#' * warnings raised in a worker are re-raised here, not lost in a child's
#'   stderr;
#' * a `data.table` returned by a worker is re-allocated
#'   ([data.table::setalloccol()]), so that `:=` on it works without the
#'   "invalid .internal.selfref" warning that follows serialisation.
#'
#' Under fork the number of workers is also capped by the memory available
#' (see `.scr_fork_cap()`): a forked worker starts as a copy-on-write clone,
#' but R's garbage collector writes to every object it marks, so in practice
#' each worker ends up owning a copy of most of the parent heap. Twenty
#' workers over a 9 GB parent exhaust a 60 GB workstation in under two
#' minutes; the cap keeps the run alive at the cost of fewer workers.
#' Forked children use a single data.table thread (data.table does this by
#' itself); PSOCK workers are set to one thread explicitly, since `k`
#' workers times the default thread count oversubscribe the machine.
#' `R CMD check --as-cran` limits every package to two processes through
#' `_R_CHECK_LIMIT_CORES_`, which is honoured.
#' @keywords internal
#' @noRd
.scr_lapply <- function(X, FUN, nthread = 1L, ..., fork_only = FALSE) {
  n <- length(X)
  k <- min(as.integer(nthread %||% 1L), n)
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(tolower(lim), "false")) k <- min(k, 2L)
  backend <- .scr_backend()
  # fork_only: FUN closes over a large table; under fork it is shared
  # copy-on-write, under PSOCK it would be serialised to every worker
  if (n < 2L || is.na(k) || k <= 1L || identical(backend, "serial") ||
      (isTRUE(fork_only) && !identical(backend, "fork"))) return(lapply(X, FUN, ...))

  wrap <- .scr_envelope(FUN)
  if (identical(backend, "fork")) {
    k <- .scr_fork_cap(k)
    if (k <= 1L) return(lapply(X, FUN, ...))
    out <- suppressWarnings(parallel::mclapply(X, wrap, ..., mc.cores = k, mc.preschedule = TRUE))
  } else {
    cl <- parallel::makeCluster(k)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L))
    out <- tryCatch(parallel::parLapply(cl, X, wrap, ...), error = function(e)
      stop("parallel worker failed (PSOCK cluster): ", conditionMessage(e),
           "; lower `nthread` or set options(scorecraft.parallel = \"serial\").", call. = FALSE))
  }
  .scr_unwrap(out, n)
}

#' Wrap `FUN` so that a worker always returns an envelope, never a condition
#'
#' The closure's environment holds only `FUN`, so a PSOCK worker receives
#' `FUN` and nothing else from the caller's frame.
#' @keywords internal
#' @noRd
.scr_envelope <- function(FUN) {
  wrap <- function(x, ...) {
    warns <- character()
    val <- withCallingHandlers(
      tryCatch(list(ok = TRUE, value = FUN(x, ...)),
               error = function(e) list(ok = FALSE, error = conditionMessage(e))),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    val$warnings <- warns
    val
  }
  environment(wrap) <- list2env(list(FUN = FUN), parent = environment(.scr_envelope))
  wrap
}

#' Open the envelopes: dead workers, errors, warnings, then the values
#' @keywords internal
#' @noRd
.scr_unwrap <- function(out, n) {
  dead <- if (length(out) != n) seq_len(n) else which(vapply(out, function(o) is.null(o) || !is.list(o), logical(1)))
  if (length(dead)) {
    stop("parallel worker died without returning a result (item(s) ", lst(dead),
         "): killed by the system, most likely out of memory. Lower `nthread`, or set ",
         "options(scorecraft.parallel = \"serial\").", call. = FALSE)
  }
  fails <- which(!vapply(out, function(o) isTRUE(o$ok), logical(1)))
  if (length(fails)) {
    i <- fails[1]
    stop("parallel worker failed (item ", i, "): ", out[[i]]$error, call. = FALSE)
  }
  warns <- unique(unlist(lapply(out, `[[`, "warnings")))
  for (w in warns) warning(w, call. = FALSE)
  lapply(out, function(o) {
    v <- o$value
    if (data.table::is.data.table(v)) data.table::setalloccol(v) else v
  })
}

#' Resolve the parallel backend from the option and the platform
#' @keywords internal
#' @noRd
.scr_backend <- function() {
  b <- getOption("scorecraft.parallel", NULL)
  if (is.null(b)) return(if (.Platform$OS.type == "unix") "fork" else "psock")
  b <- match.arg(as.character(b), c("fork", "psock", "serial"))
  if (identical(b, "fork") && .Platform$OS.type != "unix") "psock" else b
}

#' Cap the number of fork workers by the memory available
#'
#' Budget: `getOption("scorecraft.fork_mem_fraction", 0.75)` of the memory the
#' kernel reports as available, divided by the resident size of this process
#' (the worst case: every worker duplicating the parent). At least one worker
#' is always returned; `Inf` (or any non-finite fraction) disables the cap.
#' Only Linux exposes both numbers cheaply (`/proc/meminfo` and
#' `/proc/self/statm`); elsewhere the cap is a no-op.
#' @keywords internal
#' @noRd
.scr_fork_cap <- function(k) {
  k <- as.integer(k)
  frac <- suppressWarnings(as.numeric(getOption("scorecraft.fork_mem_fraction", 0.75)))
  if (length(frac) != 1L || !is.finite(frac) || frac <= 0) return(k)
  avail <- .scr_mem_available()
  rss   <- .scr_rss()
  if (!is.finite(avail) || !is.finite(rss) || rss <= 0) return(k)
  cap <- max(1L, as.integer(floor(frac * avail / rss)))
  if (cap < k) {
    msg("  fork workers capped at %d of %d (resident %.1f GB, %.1f GB available)",
        cap, k, rss / 1024^3, avail / 1024^3)
  }
  min(k, cap)
}

#' Memory available to new allocations, in bytes (Linux only, else `NA`)
#' @keywords internal
#' @noRd
.scr_mem_available <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  ln <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
  ln <- ln[startsWith(ln, "MemAvailable:")]
  if (!length(ln)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(regmatches(ln[1], regexpr("[0-9]+", ln[1]))))
  if (!is.finite(kb)) NA_real_ else kb * 1024
}

#' Resident set size of this process, in bytes (Linux only, else `NA`)
#' @keywords internal
#' @noRd
.scr_rss <- function() {
  if (!file.exists("/proc/self/statm")) return(NA_real_)
  st <- tryCatch(scan("/proc/self/statm", what = numeric(), quiet = TRUE), error = function(e) numeric())
  if (length(st) < 2L || !is.finite(st[2])) return(NA_real_)
  st[2] * 4096
}

#' Split a vector into `k` chunks of similar size, preserving order
#' @keywords internal
#' @noRd
.scr_chunks <- function(x, k) {
  k <- max(1L, min(as.integer(k), length(x)))
  if (k == 1L) return(list(x))
  unname(split(x, cut(seq_along(x), k, labels = FALSE)))
}

# -- Sampling --------------------------------------------------------------- #

#' Seed the random number generator for the calling function only
#'
#' `set.seed()` inside a package function would replace the user's random
#' stream for the rest of the session (CRAN policy: a package must not
#' change the global RNG state behind the user's back). This sets `seed`
#' and registers, in the caller's frame, an exit hook that puts
#' `.Random.seed` (and with it the RNG kind) back as it was, or removes it
#' if it did not exist. `NULL` does nothing: the caller then draws from the
#' user's stream, as documented.
#' @keywords internal
#' @noRd
.scr_local_seed <- function(seed, frame = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  genv <- globalenv()
  old  <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  restore <- function() {
    if (is.null(old)) {
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", old, envir = genv)
    }
  }
  # the closure itself goes in the call, and hooks run last-in first-out, so
  # several calls in one frame unwind to the state before the first
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE), envir = frame)
  set.seed(seed)
  invisible(NULL)
}

#' Keep the user's random stream as it is at this point
#'
#' Snapshots `.Random.seed` and registers, in the caller's frame, an exit
#' hook that restores it. Used around bootstrap loops that re-seed once per
#' replicate (so that results do not depend on the backend or the number of
#' workers): the replicate seeds are drawn from the user's stream before the
#' guard, which is therefore advanced by exactly those draws and no more.
#' @keywords internal
#' @noRd
.scr_rng_guard <- function(frame = parent.frame()) {
  genv <- globalenv()
  old  <- if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  restore <- function() {
    if (is.null(old)) {
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", old, envir = genv)
    }
  }
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE), envir = frame)
  invisible(NULL)
}

#' Indices of a subsample stratified by the target
#'
#' Used only for the classifiers: permutation importance and cross-validation
#' grow with the number of rows, and the *ordering* of importance stabilises
#' long before the whole table is used.
#' @keywords internal
#' @noRd
subsample_stratified <- function(y, max_n, seed = NULL) {
  n <- length(y)
  if (!is.finite(max_n) || n <= max_n) return(seq_len(n))
  .scr_local_seed(seed)
  frac <- max_n / n
  idx  <- integer(0)
  for (lv in unique(y)) {
    pos <- which(y == lv)
    k   <- max(1L, min(length(pos), floor(frac * length(pos))))
    idx <- c(idx, if (length(pos) == 1L) pos else sample(pos, k))
  }
  sort(idx)
}

# -- Ranking and formatting ------------------------------------------------- #

#' Rank percentile from 0 to 1, 1 being the best
#' @keywords internal
#' @noRd
rank_pct <- function(x) {
  ok <- is.finite(x)
  out <- rep(0, length(x))
  if (!any(ok)) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- if (sum(ok) == 1L) 1 else (r - 1) / (sum(ok) - 1)
  out
}

#' @keywords internal
#' @noRd
fmt_pct <- function(x, dig = 1) {
  ifelse(is.na(x), "-", sprintf(paste0("%.", dig, "f%%"), 100 * x))
}

#' Integer with thousands separator
#' @keywords internal
#' @noRd
n_fmt <- function(x) {
  formatC(as.integer(x), big.mark = ",", format = "d")
}

#' Truncate a vector of names for a readable log
#' @keywords internal
#' @noRd
lst <- function(x, n = 8) {
  if (!length(x)) return("(none)")
  if (length(x) <= n) paste(x, collapse = ", ")
  else paste0(paste(x[seq_len(n)], collapse = ", "), sprintf(" ... (+%d)", length(x) - n))
}

#' Validate a finite numeric scalar
#' @keywords internal
#' @noRd
.scr_num1 <- function(x, name, lower = -Inf, upper = Inf, open_lower = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be a single finite number.", call. = FALSE)
  }
  if (x < lower || (open_lower && x <= lower) || x > upper) {
    stop("`", name, "` is outside the allowed range (", if (open_lower) "(" else "[",
         lower, ", ", upper, "]).", call. = FALSE)
  }
  x
}

#' Format numbers for SQL without losing precision or falling into scientific notation
#' @keywords internal
#' @noRd
.sql_num <- function(x) {
  vapply(as.double(x), function(v) {
    if (!is.finite(v)) return("NULL")
    s <- format(v, digits = 15L, scientific = FALSE, trim = TRUE)
    if (!isTRUE(all.equal(as.numeric(s), v, tolerance = 0)))
      s <- format(v, digits = 17L, scientific = FALSE, trim = TRUE)
    s
  }, character(1), USE.NAMES = FALSE)
}

#' Quote a string literal for SQL
#'
#' A single quote is always doubled (SQL-92). The backslash is an escape
#' character only in MySQL/MariaDB, Spark/Hive/Databricks and BigQuery; in
#' ANSI SQL, PostgreSQL (standard_conforming_strings), SQL Server, Oracle,
#' SQLite and DuckDB it is an ordinary character, and doubling it there
#' would change the literal. `dialect = NULL` keeps the historical
#' behaviour (backslash doubled); callers that know the dialect should pass
#' it. `NA` becomes `NULL`.
#' @keywords internal
#' @noRd
.sql_str <- function(x, dialect = NULL) {
  x <- as.character(x)
  bs_escape <- is.null(dialect) ||
    dialect %in% c("mysql", "mariadb", "spark", "hive", "databricks", "bigquery")
  s <- if (bs_escape) gsub("\\\\", "\\\\\\\\", x) else x
  out <- paste0("'", gsub("'", "''", s, fixed = TRUE), "'")
  out[is.na(x)] <- "NULL"
  out
}
