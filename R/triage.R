# ============================================================================ #
# triage.R - Stage 1: descriptive triage
# ============================================================================ #
# Binning is not the bottleneck; this stage does not try to guess predictive
# power. It does what the binning engine does not: kills the structurally
# useless, applies a cheap noise floor, and resolves sentinels and missing
# values by decomposition (a categorical flag of its own + the training
# median). After it the data has ZERO missing and ZERO sentinels - a verified
# invariant.
# ============================================================================ #

#' Stage 1: descriptive triage and sentinel resolution
#'
#' Profiles every candidate **on the training rows only**, decides its fate
#' and materialises the clean data for train and hold-out with the same
#' values (training median, `"MISSING"` level). A sentinel or missing mass
#' with weight (`special_min_share`) and signal (`special_min_woe`) becomes a
#' categorical flag `<column><flag_suffix>`, which the engine bins and emits
#' in SQL natively.
#'
#' Failures at this stage: `CONSTANT`, `NEAR_CONSTANT`, `TOO_MANY_MISSING`,
#' `HIGH_CARDINALITY`, `NO_SIGNAL` (coarse IV below `min_iv_quick`) and
#' `DUPLICATE_OF:<column>`. A failed numeric whose special-population flag
#' survives carries the suffix `;RESCUED_AS_FLAG`.
#'
#' @param split An object from [scr_split()].
#' @param config An object from [scr_config()].
#'
#' @return An `scr_triage` object with `profile` (one row per candidate and
#'   derived flag), `ledger` (the source of truth of the pre-processing the
#'   SQL reproduces), `keep`, `derived`, `clean` (target + survivors + flags,
#'   with no `NA`) and the originating `split`.
#'
#' @family stages
#' @examples
#' sp <- scr_split(scr_demo, "default", date_col = "ref_date", drop = "id")
#' tr <- scr_triage(sp, scr_config(verbose = FALSE))
#' tr
#' table(tr$profile$triage_reason)
#' @export
scr_triage <- function(split, config = scr_config()) {
  if (!inherits(split, "scr_split")) stop("`split` must come from scr_split().", call. = FALSE)
  check_config(config, "scr_triage")
  old <- scr_verbose(isTRUE(config$verbose)); on.exit(scr_verbose(old), add = TRUE)

  msg_stage(1, "descriptive triage")
  plan  <- time_it("profile + heuristics", triage_plan(split$data, split$target, split$cols,
                                                        split$train_idx, config))
  clean <- time_it("materialisation (imputation + flags)", apply_triage(split$data, split$target, plan, config))
  dead <- plan$profile[triage_status == "drop", .N, by = triage_reason][order(-N)]
  for (i in seq_len(nrow(dead))) msg("  failed for %-22s %d", dead$triage_reason[i], dead$N[i])
  msg("  survived: %d (%d derived from special populations)", length(plan$keep), length(plan$derived))
  if (!length(plan$keep)) stop("Triage left no candidate - review the configuration limits.", call. = FALSE)

  structure(list(profile = plan$profile, ledger = plan$ledger, keep = plan$keep,
                 keep_num = plan$keep_num, keep_cat = plan$keep_cat, derived = plan$derived,
                 clean = clean, split = split[setdiff(names(split), "data")], config = config),
            class = c("scr_triage", "list"))
}

#' @export
print.scr_triage <- function(x, ...) {
  cat(sprintf("<scr_triage> target \"%s\" | %d candidates -> %d survivors (%d derived)\n",
              x$split$target, nrow(x$profile[is.na(derived_from)]), length(x$keep), length(x$derived)))
  m <- x$profile[triage_status == "drop", .N, by = triage_reason][order(-N)]
  for (i in seq_len(nrow(m))) cat(sprintf("  %-24s %d\n", m$triage_reason[i], m$N[i]))
  invisible(x)
}

#' Profile and decide the fate of each candidate (training statistics ONLY)
#' @keywords internal
#' @noRd
triage_plan <- function(dt, target, cols, train_idx, cfg) {

  y  <- dt[[target]][train_idx]
  n  <- length(train_idx)
  sp <- cfg$special_values

  # One candidate per task, parallel by column. The training slice is
  # taken inside the worker, so under fork nothing is copied up front.
  feats <- cols$features
  out <- .scr_lapply(feats, function(f) {
    .triage_one(f, dt[[f]][train_idx], y, f %in% cols$var_num, sp, cfg, n)
  }, nthread = cfg$nthread, fork_only = TRUE)
  rows   <- lapply(out, `[[`, "row")
  ledger <- unlist(lapply(out, `[[`, "ledger"), recursive = FALSE)
  fps    <- vapply(out, `[[`, character(1), "fp")

  profile <- data.table::rbindlist(rows)
  ledger  <- if (length(ledger)) data.table::rbindlist(ledger) else
    data.table::data.table(kind = character(), source = character(), output = character(),
                           impute_value = numeric(), specials = character(), na_fill = character())
  names(fps) <- cols$features

  # -- profile rows of the derived flags ---------------------------------- #
  flags <- ledger[kind == "num_flag"]
  if (nrow(flags)) {
    extra <- vector("list", nrow(flags)); fps_fl <- character(nrow(flags))
    for (j in seq_len(nrow(flags))) {
      src <- flags$source[j]
      p   <- profile[feature == src]
      g   <- .flag_levels(dt[[src]][train_idx], sp)
      nd  <- data.table::uniqueN(g)
      fps_fl[j] <- paste("F", p$n_missing, p$n_special, paste(sort(unique(g)), collapse = ","), sep = "|")
      extra[[j]] <- data.table::data.table(
        feature = flags$output[j], derived_from = src, type = "categorical",
        n_missing = 0L, pct_missing = 0, n_special = p$n_special, pct_special = p$pct_special,
        n_distinct = nd, pct_mode = max(tabulate(match(g, unique(g)))) / n,
        iv_quick = scr_iv(g, y), woe_special = p$woe_special, decomposition = "derived",
        triage_status = if (nd < 2L) "drop" else "keep",
        triage_reason = if (nd < 2L) "CONSTANT" else "OK")
    }
    names(fps_fl) <- flags$output
    fps <- c(fps, fps_fl)
    profile <- data.table::rbindlist(list(profile, data.table::rbindlist(extra)))
    dead <- profile[derived_from %in% flags$source & triage_status == "drop", feature]
    if (length(dead)) ledger <- ledger[!(kind == "num_flag" & output %in% dead)]
    rescued <- profile[!is.na(derived_from) & triage_status == "keep" &
                       derived_from %in% profile[triage_status == "drop", feature], derived_from]
    if (length(rescued)) profile[feature %in% rescued, triage_reason := paste0(triage_reason, ";RESCUED_AS_FLAG")]
  }

  # -- exact duplicates among survivors (flags included) ------------------- #
  if (isTRUE(cfg$check_duplicates)) {
    alive <- profile[triage_status == "keep", feature]
    orig  <- profile[!is.na(derived_from), .(feature, derived_from)]
    value <- function(f) {
      s <- orig[feature == f, derived_from]
      if (length(s)) .flag_levels(dt[[s[1]]][train_idx], sp) else dt[[f]][train_idx]
    }
    dups <- .find_duplicates(value, alive, fps[alive])
    if (length(dups)) {
      dd <- data.table::data.table(feature = names(dups), of = unlist(dups, use.names = FALSE))
      profile[dd, on = "feature", `:=`(triage_status = "drop", triage_reason = paste0("DUPLICATE_OF:", i.of))]
      ledger <- ledger[!(output %in% names(dups))]
    }
  }

  keep <- profile[triage_status == "keep", feature]
  list(profile = profile, ledger = ledger, keep = keep,
       keep_num = intersect(keep, cols$var_num),
       keep_cat = c(intersect(keep, cols$var_cat), profile[!is.na(derived_from) & triage_status == "keep", feature]),
       derived  = profile[!is.na(derived_from) & triage_status == "keep", feature])
}

#' Profile one candidate on its training slice: verdict, ledger rows, fingerprint
#' @keywords internal
#' @noRd
.triage_one <- function(f, x, y, num, sp, cfg, n) {
  n_miss <- sum(is.na(x))
  n_spec <- if (num && length(sp)) sum(x %in% sp, na.rm = TRUE) else 0L
  share  <- (n_miss + n_spec) / n

  if (num) {
    reg      <- !is.na(x) & !(x %in% sp)
    xr       <- x[reg]
    u        <- unique(xr)
    n_dist   <- length(u)
    pct_mode <- if (length(xr)) max(tabulate(match(xr, u), nbins = n_dist)) / length(xr) else 1
    imput    <- if (length(xr)) stats::median(xr) else NA_real_
    fp       <- paste("N", n_dist, n_miss, n_spec, signif(mean(xr), 12), signif(stats::sd(xr), 12),
                      signif(suppressWarnings(min(xr)), 12), signif(suppressWarnings(max(xr)), 12), sep = "|")
  } else {
    xc       <- if (n_miss) ifelse(is.na(x), "MISSING", x) else x
    u        <- unique(xc)
    n_dist   <- length(u)
    cnt      <- tabulate(match(xc, u), nbins = n_dist)
    ord      <- order(-cnt)[seq_len(min(5L, n_dist))]
    pct_mode <- max(cnt) / n
    imput    <- NA_real_
    fp       <- paste("C", n_dist, n_miss, paste(u[ord], collapse = ","), paste(cnt[ord], collapse = ","), sep = "|")
  }

  iv_q  <- .quick_iv(x, y, num, sp, cfg)
  # sentinels are numeric-only: for text, `x %in% sp` would coerce -999 to
  # "-999" and count a legitimate level as special
  w_sp  <- if (share > 0) woe_subpop(if (num) is.na(x) | x %in% sp else is.na(x), y) else NA_real_

  decomp <- "-"
  if (num && share >= cfg$special_min_share && share <= 1 - cfg$special_min_share &&
      is.finite(w_sp) && abs(w_sp) >= cfg$special_min_woe) {
    decomp <- "flag"
  } else if (num && (n_miss + n_spec) > 0) {
    decomp <- "impute"
  } else if (!num && n_miss > 0) {
    decomp <- "coalesce"
  }

  status <- "keep"; reason <- "OK"
  if (num) {
    if (n_dist < 2L)                      { status <- "drop"; reason <- "CONSTANT" }
    else if (share > cfg$max_missing)      { status <- "drop"; reason <- "TOO_MANY_MISSING" }
    else if (pct_mode > cfg$near_constant) { status <- "drop"; reason <- "NEAR_CONSTANT" }
    else if (iv_q < cfg$min_iv_quick)      { status <- "drop"; reason <- "NO_SIGNAL" }
  } else {
    if (n_dist < 2L)                        { status <- "drop"; reason <- "CONSTANT" }
    else if (n_dist > cfg$max_cat_levels)   { status <- "drop"; reason <- "HIGH_CARDINALITY" }
    else if (pct_mode > cfg$near_constant)  { status <- "drop"; reason <- "NEAR_CONSTANT" }
    else if (iv_q < cfg$min_iv_quick)       { status <- "drop"; reason <- "NO_SIGNAL" }
  }

  row <- data.table::data.table(
    feature = f, derived_from = NA_character_,
    type = if (num) "numeric" else "categorical",
    n_missing = n_miss, pct_missing = n_miss / n, n_special = n_spec, pct_special = n_spec / n,
    n_distinct = n_dist, pct_mode = pct_mode, iv_quick = iv_q, woe_special = w_sp,
    decomposition = decomp, triage_status = status, triage_reason = reason)

  # Derivation ledger: UNCONDITIONAL for every survivor, even when training
  # saw neither a missing value nor a sentinel. Hold-out and production may;
  # this way R and SQL always agree.
  led <- list()
  if (num && status == "keep" && is.finite(imput)) {
    led[[length(led) + 1L]] <- data.table::data.table(
      kind = "num_impute", source = f, output = f, impute_value = imput,
      specials = paste(sp, collapse = ","), na_fill = NA_character_)
  }
  if (num && identical(decomp, "flag")) {
    led[[length(led) + 1L]] <- data.table::data.table(
      kind = "num_flag", source = f, output = paste0(f, cfg$flag_suffix), impute_value = NA_real_,
      specials = paste(sp, collapse = ","), na_fill = NA_character_)
  }
  if (!num && status == "keep") {
    led[[length(led) + 1L]] <- data.table::data.table(
      kind = "cat_coalesce", source = f, output = f, impute_value = NA_real_,
      specials = NA_character_, na_fill = "MISSING")
  }
  list(row = row, ledger = led, fp = fp)
}

#' Materialise the triage plan (same transformation on train and hold-out)
#' @keywords internal
#' @noRd
apply_triage <- function(dt, target, plan, cfg) {
  sp  <- cfg$special_values
  out <- data.table::data.table(.tmp = seq_len(nrow(dt)))
  data.table::set(out, j = target, value = dt[[target]])

  imput <- plan$ledger[kind == "num_impute"]
  flags <- plan$ledger[kind == "num_flag"]
  coal  <- plan$ledger[kind == "cat_coalesce"]
  # set() cannot add a column beyond the over-allocated slots (1024 by
  # default): reserve room for every output column up front
  data.table::setalloccol(out, length(plan$keep) + nrow(flags) + 2L + getOption("datatable.alloccol", 1024L))

  # ledger rows matched once: a lookup (a data.table subset) per column is
  # O(p^2) on a wide table
  feats  <- setdiff(plan$keep, plan$derived)
  imp_v  <- imput$impute_value[match(feats, imput$source)]
  coal_f <- feats %in% coal$source
  flag_o <- flags$output[match(feats, flags$source)]

  for (i in seq_along(feats)) {
    f <- feats[i]
    x <- dt[[f]]
    if (is.numeric(x)) {
      v <- imp_v[i]
      if (is.finite(v)) {
        bad <- is.na(x) | x %in% sp
        if (any(bad)) x[bad] <- v
      }
      data.table::set(out, j = f, value = x)
    } else {
      if (coal_f[i] && anyNA(x)) x[is.na(x)] <- "MISSING"
      data.table::set(out, j = f, value = as.character(x))
    }
    if (!is.na(flag_o[i])) data.table::set(out, j = flag_o[i], value = .flag_levels(dt[[f]], sp))
  }
  orphans <- flags[!source %in% plan$keep]
  for (j in seq_len(nrow(orphans))) {
    data.table::set(out, j = orphans$output[j], value = .flag_levels(dt[[orphans$source[j]]], sp))
  }
  out[, .tmp := NULL]

  missing <- setdiff(plan$keep, names(out))
  if (length(missing)) stop("apply_triage(): the plan promises column(s) the data lacks: ", lst(missing), call. = FALSE)
  leftovers <- names(out)[vapply(out, anyNA, logical(1))]
  if (length(leftovers)) {
    stop("apply_triage(): NA left in ", lst(leftovers),
         " - the Stage 1 invariant (zero missing) was violated.", call. = FALSE)
  }
  out[]
}

#' Levels of the special-population flag
#' @keywords internal
#' @noRd
.flag_levels <- function(x, sp) {
  g <- rep("REGULAR", length(x))
  if (length(sp)) {
    hit <- !is.na(x) & x %in% sp
    if (any(hit)) g[hit] <- paste0("S", x[hit])
  }
  g[is.na(x)] <- "MISSING"
  g
}

#' Coarse IV: quantiles for a numeric, levels for a categorical
#' @keywords internal
#' @noRd
.quick_iv <- function(x, y, num, sp, cfg) {
  if (num) {
    g   <- .flag_levels(x, sp)
    reg <- g == "REGULAR"
    if (sum(reg) > 1L) {
      br <- unique(stats::quantile(x[reg], probs = seq(0, 1, length.out = cfg$quick_iv_groups + 1L),
                                   na.rm = TRUE, type = 1L))
      g[reg] <- if (length(br) >= 3L) paste0("Q", cut(x[reg], br, include.lowest = TRUE, labels = FALSE)) else "R"
    }
  } else {
    g <- ifelse(is.na(x), "MISSING", as.character(x))
    if (data.table::uniqueN(g) > 50L) {
      tb <- table(g); rare <- names(tb)[tb < 0.005 * length(g)]
      if (length(rare)) g[g %in% rare] <- "__RARE__"
    }
  }
  scr_iv(g, y)
}

#' Exact duplicates: cheap fingerprint first, exact comparison only on collisions
#' @keywords internal
#' @noRd
.find_duplicates <- function(value, feats, fp) {
  if (length(feats) < 2L) return(list())
  stopifnot(length(fp) == length(feats))
  dups <- list()
  for (grp in split(feats, fp)) {
    if (length(grp) < 2L) next
    vals <- lapply(grp, value); names(vals) <- grp
    for (a in seq_along(grp)) {
      if (grp[a] %in% names(dups)) next
      for (b in seq_along(grp)) {
        if (b <= a || grp[b] %in% names(dups)) next
        if (identical(vals[[a]], vals[[b]])) dups[[grp[b]]] <- grp[a]
      }
    }
  }
  dups
}
