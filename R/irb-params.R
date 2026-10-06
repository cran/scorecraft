# ============================================================================ #
# irb-params.R - parameter tables of the IRB framework, by preset
# ============================================================================ #
# The package is a technical tool: a regime is a table of numbers (floors,
# correlations, supervisory values, standardised risk weights) selected by a
# preset, never prose. Every table is a data.table the user may edit and pass
# back; the functions that consume it record `params_modified` in the ledger.
# ============================================================================ #

#' IRB parameter tables by framework preset
#'
#' Returns the numeric tables the internal ratings-based (IRB) functions
#' read: probability of default (PD) floors, loss given default (LGD) input
#' floors for own estimates,
#' supervisory LGD values of the foundation approach, standardised credit
#' conversion factors (CCF), asset-correlation parameters of the risk-weight
#' function, maturity rules, the output floor and the standardised risk
#' weights used for the floor comparison. Three presets ship: `"bcb"`
#' (Brazil, BCB Resolutions 303/2023 and 229/2022), `"basel3_final"` (the
#' consolidated Basel Framework in force from 2023) and `"crr3"` (the EU
#' text applicable from 2025). The presets differ in a handful of cells, all
#' visible with `print()`; users who need another jurisdiction edit the
#' tables and pass the object to the functions that take `params`.
#'
#' @section Regulatory texts:
#' The presets encode a reading of the texts below at the time of the
#' release. Regulation changes and is interpreted by each supervisor: before
#' any regulatory use, check every table against the texts in force for the
#' jurisdiction and the portfolio, and edit the tables where they differ.
#' The package implements the calculations; it does not give regulatory
#' advice.
#'
#' @references
#' Basel Committee on Banking Supervision. *The Basel Framework*, chapters
#' CRE31 to CRE36 (IRB approach) and RBC20 (output floor).
#'
#' Regulation (EU) 2024/1623 (CRR3), amending Regulation (EU) No 575/2013.
#'
#' European Banking Authority (2017). *Guidelines on PD estimation, LGD
#' estimation and the treatment of defaulted exposures*, EBA/GL/2017/16.
#'
#' European Banking Authority (2019). *Guidelines for the estimation of LGD
#' appropriate for an economic downturn*, EBA/GL/2019/03.
#'
#' Banco Central do Brasil. BCB Resolution 229/2022 and BCB Resolution
#' 303/2023.
#'
#' International Accounting Standards Board (2014). *IFRS 9 Financial
#' Instruments*.
#'
#' @param framework `"bcb"`, `"basel3_final"` or `"crr3"`.
#'
#' @return An object of class `scr_irb_params`: a list with `framework`,
#'   `source` (one line), `pd_floor`, `lgd_floor`, `lgd_firb`, `ccf_sa`,
#'   `ccf_floor_fraction`, `correlation`, `scaling_factor`, `confidence`,
#'   `m_default`, `m_range`, `output_floor`, `sa_rw` and `modified` (logical,
#'   set by the functions that receive the object when its tables were
#'   edited).
#'
#' @family irb-parameters
#' @examples
#' p <- scr_irb_params("bcb")
#' p
#' p$pd_floor
#' p2 <- p; p2$pd_floor$floor[p2$pd_floor$asset_class == "retail_other"] <- 0.001
#' @export
scr_irb_params <- function(framework = c("bcb", "basel3_final", "crr3")) {
  framework <- match.arg(framework)
  bcb <- identical(framework, "bcb")

  pd_floor <- data.table::data.table(
    asset_class = c("corporate", "bank", "sovereign", "retail_mortgage",
                    "qrre_transactor", "qrre_revolver", "retail_other"),
    floor = c(5e-4, 5e-4, NA_real_, 5e-4, 5e-4, 1e-3, 5e-4))

  lgd_floor <- data.table::data.table(
    asset_class    = c("corporate", "retail_mortgage", "qrre", "retail_other"),
    unsecured      = c(0.25, NA_real_, 0.50, 0.30),
    financial      = c(0.00, NA_real_, NA_real_, 0.00),
    receivables    = c(0.10, NA_real_, NA_real_, 0.10),
    real_estate    = c(0.10, 0.05, NA_real_, 0.10),
    other_physical = c(0.15, NA_real_, NA_real_, 0.15))

  lgd_firb <- if (bcb) {
    data.table::data.table(
      claim = c("senior_unsecured", "priority_claim", "subordinated",
                "secured_financial", "secured_receivables", "secured_real_estate", "secured_other"),
      lgd   = c(0.75, 0.45, 0.75, 0.00, 0.20, 0.20, 0.25))
  } else {
    data.table::data.table(
      claim = c("senior_unsecured_corporate", "senior_unsecured_fi", "subordinated",
                "secured_financial", "secured_receivables", "secured_real_estate", "secured_other"),
      lgd   = c(0.40, 0.45, 0.75, 0.00, 0.20, 0.20, 0.25))
  }

  ccf_sa <- data.table::data.table(
    item = c("uncond_cancellable", "commitment", "nif_ruf", "direct_substitute"),
    ccf  = c(0.10, 0.40, 0.50, 1.00))

  correlation <- list(
    corporate       = c(lo = 0.12, hi = 0.24, k = 50),
    hvcre           = c(lo = 0.12, hi = 0.30, k = 50),
    retail_mortgage = 0.15,
    qrre            = 0.04,
    retail_other    = c(lo = 0.03, hi = 0.16, k = 35),
    fi_multiplier   = 1.25,
    sme             = if (bcb) list(lo = 15, hi = 300, adj = 0.04, unit = "BRL m")
                      else     list(lo = 5,  hi = 50,  adj = 0.04, unit = "EUR m"))

  sa_rw <- data.table::rbindlist(list(
    data.table::data.table(asset_class = "retail_other", sub_class = c("standard", "transactor", "non_granular"),
                           ltv_lo = NA_real_, ltv_hi = NA_real_, rw = c(0.75, 0.45, 1.00)),
    data.table::data.table(asset_class = "retail_mortgage", sub_class = "standard",
                           ltv_lo = c(0, 0.5, 0.6, 0.8, 0.9, 1.0), ltv_hi = c(0.5, 0.6, 0.8, 0.9, 1.0, Inf),
                           rw = c(0.20, 0.25, 0.30, 0.40, 0.50, 0.70)),
    data.table::data.table(asset_class = "retail_mortgage", sub_class = "income_producing",
                           ltv_lo = c(0, 0.5, 0.6, 0.8, 0.9, 1.0), ltv_hi = c(0.5, 0.6, 0.8, 0.9, 1.0, Inf),
                           rw = c(0.30, 0.35, 0.45, 0.60, 0.75, 1.05)),
    data.table::data.table(asset_class = "corporate",
                           sub_class = c("rated_aaa_aa", "rated_a", "rated_bbb", "rated_bb", "rated_below_bb",
                                         "unrated", "investment_grade_unrated", "sme"),
                           ltv_lo = NA_real_, ltv_hi = NA_real_,
                           rw = c(0.20, 0.50, 0.75, 1.00, 1.50, 1.00, 0.65, 0.85)),
    data.table::data.table(asset_class = "defaulted", sub_class = c("provision_lt_20", "provision_ge_20", "mortgage"),
                           ltv_lo = NA_real_, ltv_hi = NA_real_, rw = c(1.50, 1.00, 1.00))
  ), use.names = TRUE)

  source <- switch(framework,
    bcb          = "BCB Resolutions 303/2023 (IRB) and 229/2022 (standardised); values as tables, editable",
    basel3_final = "Basel Framework CRE20/CRE30-36 as in force from 2023; values as tables, editable",
    crr3         = "Regulation (EU) 575/2013 as amended by 2024/1623; values as tables, editable")

  structure(list(
    framework = framework, source = source,
    pd_floor = pd_floor, lgd_floor = lgd_floor, lgd_firb = lgd_firb, ccf_sa = ccf_sa,
    ccf_floor_fraction = 0.5, correlation = correlation,
    scaling_factor = 1, confidence = 0.999, m_default = 2.5, m_range = c(1, 5),
    output_floor = 0.725, sa_rw = sa_rw, modified = FALSE
  ), class = c("scr_irb_params", "list"))
}

#' @export
print.scr_irb_params <- function(x, ...) {
  cat(sprintf("<scr_irb_params> framework: %s%s\n", x$framework, if (isTRUE(x$modified)) " (modified)" else ""))
  cat("  ", x$source, "\n", sep = "")
  pf <- x$pd_floor
  cat("  PD floors:  ", paste(sprintf("%s %s", pf$asset_class, ifelse(is.na(pf$floor), "none", fmt_pct(pf$floor, 2))),
                              collapse = " | "), "\n")
  lf <- x$lgd_floor
  cat("  LGD floors (unsecured): ", paste(sprintf("%s %s", lf$asset_class, ifelse(is.na(lf$unsecured), "n/a", fmt_pct(lf$unsecured, 0))),
                                         collapse = " | "), "\n")
  cat(sprintf("  F-IRB LGD: %s\n", paste(sprintf("%s %s", x$lgd_firb$claim, fmt_pct(x$lgd_firb$lgd, 0)), collapse = " | ")))
  cat(sprintf("  CCF (standardised): %s | own-estimate floor %s of the standardised value\n",
              paste(sprintf("%s %s", x$ccf_sa$item, fmt_pct(x$ccf_sa$ccf, 0)), collapse = " | "),
              fmt_pct(x$ccf_floor_fraction, 0)))
  co <- x$correlation
  cat(sprintf("  correlation: corporate %.2f-%.2f (k=%d) | mortgage %.2f | QRRE %.2f | other retail %.2f-%.2f (k=%d) | FI x%.2f | SME adj %.2f (%s %g-%g)\n",
              co$corporate[["lo"]], co$corporate[["hi"]], as.integer(co$corporate[["k"]]), co$retail_mortgage, co$qrre,
              co$retail_other[["lo"]], co$retail_other[["hi"]], as.integer(co$retail_other[["k"]]), co$fi_multiplier,
              co$sme$adj, co$sme$unit, co$sme$lo, co$sme$hi))
  cat(sprintf("  confidence %.3f | scaling factor %g | M default %.1f in [%g, %g] | output floor %s | SA risk weights: %d rows\n",
              x$confidence, x$scaling_factor, x$m_default, x$m_range[1], x$m_range[2], fmt_pct(x$output_floor, 1), nrow(x$sa_rw)))
  invisible(x)
}

#' Validate a params object and detect user edits
#' @keywords internal
#' @noRd
.check_params <- function(params, fn) {
  if (!inherits(params, "scr_irb_params")) {
    stop(sprintf("%s(): `params` must come from scr_irb_params().", fn), call. = FALSE)
  }
  ref <- scr_irb_params(params$framework)
  same <- isTRUE(all.equal(unclass(params)[setdiff(names(ref), "modified")],
                           unclass(ref)[setdiff(names(ref), "modified")], check.attributes = FALSE))
  params$modified <- !same
  params
}

#' PD floor for an asset class (NA when the framework has none)
#' @keywords internal
#' @noRd
.pd_floor_of <- function(params, asset_class) {
  f <- params$pd_floor$floor[match(asset_class, params$pd_floor$asset_class)]
  if (anyNA(match(asset_class, params$pd_floor$asset_class))) {
    stop("unknown `asset_class`: ", lst(setdiff(unique(asset_class), params$pd_floor$asset_class)),
         ". See scr_irb_params()$pd_floor.", call. = FALSE)
  }
  f
}
