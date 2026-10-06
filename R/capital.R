# ============================================================================ #
# capital.R - expected loss, IRB risk weights, standardised comparison, capital
# ============================================================================ #
# Closed-form arithmetic on per-exposure vectors: the asymptotic single risk
# factor function of the IRB approach, the input floors, the standardised
# lookup for the output floor, the EL-versus-provisions comparison and a
# fixed sensitivity grid. Everything is vectorised; the parameter tables come
# from scr_irb_params() and are never hard-coded here. The Vasicek helpers
# live in this file and are shared with the PD module.
# ============================================================================ #

.irb_classes <- c("corporate", "corporate_sme", "bank", "sovereign", "hvcre",
                  "retail_mortgage", "qrre_revolver", "qrre_transactor", "retail_other")
.irb_wholesale <- c("corporate", "corporate_sme", "bank", "sovereign", "hvcre")

#' Expected loss per exposure
#'
#' The primitive every other function of the module uses: `pd * lgd * ead`
#' for performing exposures and `elbe * ead` for defaulted ones, `elbe`
#' being the best estimate of expected loss. A defaulted exposure without
#' `elbe` uses `lgd` (PD equal to one). Arguments are recycled to a common
#' length.
#'
#' @param pd,lgd,ead Numeric vectors: probability of default, loss given
#'   default (decimals) and exposure at default (currency).
#' @param defaulted Optional 0/1 or logical vector.
#' @param elbe Optional vector with the best estimate of expected loss of
#'   the defaulted rows (decimal of `ead`); ignored on performing rows.
#'
#' @return A numeric vector with the expected loss in currency.
#'
#' @family irb-capital
#' @examples
#' scr_el(c(0.01, 0.02), 0.45, c(1000, 2000))
#' scr_el(0.02, 0.45, 1000, defaulted = TRUE, elbe = 0.6)
#' @export
scr_el <- function(pd, lgd, ead, defaulted = NULL, elbe = NULL) {
  n <- max(length(pd), length(lgd), length(ead), length(defaulted), length(elbe))
  pd <- rep_len(as.double(pd), n); lgd <- rep_len(as.double(lgd), n); ead <- rep_len(as.double(ead), n)
  d <- if (is.null(defaulted)) rep(FALSE, n) else rep_len(isTRUE_vec(defaulted), n)
  e <- if (is.null(elbe)) rep(NA_real_, n) else rep_len(as.double(elbe), n)
  e[is.na(e)] <- lgd[is.na(e)]
  ifelse(d, e * ead, pd * lgd * ead)
}

# -- Vasicek helpers --------------------------------------------------------- #

#' Vasicek conditional PD given a systematic factor
#'
#' `N((G(pd) - sqrt(rho) z) / sqrt(1 - rho))`: the point-in-time PD implied
#' by a through-the-cycle PD when the systematic factor takes the value `z`
#' (positive `z` is a good year). Shared with the PD module.
#' @keywords internal
#' @noRd
.vasicek_pit <- function(pd, z, rho) {
  z <- as.double(z); rho <- as.double(rho)
  if (any(!is.na(rho) & (rho < 0 | rho >= 1))) stop("`rho` must be in [0, 1).", call. = FALSE)
  if (!is.matrix(pd)) pd <- as.double(pd)     # a matrix keeps its shape (n x T, z recycled by row)
  out <- stats::pnorm((stats::qnorm(pd) - sqrt(rho) * z) / sqrt(1 - rho))
  out[pd <= 0] <- 0; out[pd >= 1] <- 1
  out
}

#' Stressed PD of the one-factor model
#'
#' The conditional PD at confidence `q`: `N((G(pd) + sqrt(rho) G(q)) /
#' sqrt(1 - rho))`. With `q = 0.999` and the regulatory correlation it is
#' the PD inside the risk-weight function, so that the capital requirement
#' of a retail exposure is `lgd * (scr_pd_stress(pd, r, 0.999) - pd)`. Used
#' by the sensitivity grid of [scr_capital()] and by the scenario engine of
#' [scr_ecl()]. Arguments are recycled.
#'
#' @param pd Numeric vector of unconditional PDs.
#' @param rho Asset correlation in `[0, 1)`.
#' @param q Confidence level in `(0, 1)`; `0.5` returns the median-year PD.
#'
#' @return A numeric vector of conditional PDs.
#'
#' @references
#' Vasicek, O. (2002). The distribution of loan portfolio value. *Risk*,
#' 15(12), 160-162. Gordy, M. B. (2003). A risk-factor model foundation
#' for ratings-based bank capital rules. *Journal of Financial
#' Intermediation*, 12(3), 199-232.
#'
#' @seealso [scr_pd_pit_ttc()], the same bridge with the systematic factor
#'   given as a value of `z` rather than a quantile `q`.
#' @family irb-capital
#' @examples
#' scr_pd_stress(0.02, rho = 0.15, q = c(0.5, 0.95, 0.99, 0.999))
#' @export
scr_pd_stress <- function(pd, rho, q) {
  q <- as.double(q)
  if (any(is.na(q) | q <= 0 | q >= 1)) stop("`q` must be in (0, 1).", call. = FALSE)
  if (any(!is.na(pd) & (pd < 0 | pd > 1))) stop("`pd` must be in [0, 1].", call. = FALSE)
  .vasicek_pit(pd, -stats::qnorm(q), rho)
}

# -- Correlation, maturity and the K function ------------------------------- #

#' Asset correlation of the risk-weight function
#' @keywords internal
#' @noRd
.irb_correlation <- function(pd, asset_class, sales = NULL, fi = FALSE, params) {
  co <- params$correlation
  w <- function(k) (1 - exp(-k * pd)) / (1 - exp(-k))
  r <- rep(NA_real_, length(pd))
  wc <- asset_class %in% c("corporate", "corporate_sme", "bank", "sovereign")
  if (any(wc)) { cc <- co$corporate; r[wc] <- (cc[["lo"]] * w(cc[["k"]]) + cc[["hi"]] * (1 - w(cc[["k"]])))[wc] }
  hv <- asset_class == "hvcre"
  if (any(hv)) { cc <- co$hvcre; r[hv] <- (cc[["lo"]] * w(cc[["k"]]) + cc[["hi"]] * (1 - w(cc[["k"]])))[hv] }
  r[asset_class == "retail_mortgage"] <- co$retail_mortgage
  r[asset_class %in% c("qrre_revolver", "qrre_transactor")] <- co$qrre
  ro <- asset_class == "retail_other"
  if (any(ro)) { cc <- co$retail_other; r[ro] <- (cc[["lo"]] * w(cc[["k"]]) + cc[["hi"]] * (1 - w(cc[["k"]])))[ro] }
  sme <- asset_class == "corporate_sme"
  if (any(sme)) {
    s <- co$sme
    S <- if (is.null(sales)) rep(NA_real_, length(pd)) else rep_len(as.double(sales), length(pd))
    # the firm-size adjustment needs reported sales (CRE31.9, CRR Art. 153(4)):
    # without them the obligor is treated as a large corporate (no adjustment)
    S[is.na(S)] <- s$hi
    S <- pmin(s$hi, pmax(s$lo, S))
    r[sme] <- r[sme] - s$adj * (1 - (S[sme] - s$lo) / (s$hi - s$lo))
  }
  f <- rep_len(isTRUE_vec(fi), length(pd)) & asset_class %in% .irb_wholesale
  if (any(f)) r[f] <- pmin(0.999, r[f] * co$fi_multiplier)
  r
}

#' Lowest PD at which the maturity slope b(PD) is evaluated
#'
#' `1 - 1.5 b(PD)` reaches zero at PD of about 2.9e-6, so the maturity
#' adjustment explodes and changes sign there, and K stops being increasing
#' in PD below about 1e-5 for any M in `[0, 5]` and R in `[0.03, 0.5]`. Only
#' rows without a PD floor (sovereigns, or `apply_floors` without `"pd"`)
#' reach that region; b is held at its value at this PD below it, which
#' leaves every PD at or above 0.001 % unchanged.
#' @keywords internal
#' @noRd
.irb_b_pd_min <- 1e-5

#' Capital requirement K of the asymptotic single risk factor model
#' @keywords internal
#' @noRd
.irb_k_core <- function(pd, lgd, r, ma, params) {
  g <- stats::qnorm(params$confidence)
  cond <- stats::pnorm(stats::qnorm(pd) / sqrt(1 - r) + sqrt(r / (1 - r)) * g)
  cond[pd >= 1] <- 1
  (lgd * cond - pd * lgd) * ma * params$scaling_factor
}

#' Vectorised risk-weight engine with switchable floors and shocks
#'
#' `floors` is a subset of `c("pd", "lgd", "m")`; `r_mult` scales the
#' correlation (sensitivity grid). Returns the table of [scr_irb_rw()].
#' @keywords internal
#' @noRd
.irb_rw <- function(pd, lgd, ead = 1, m = NULL, asset_class, sales = NULL, fi = FALSE, defaulted = NULL,
                    elbe = NULL, params, approach = "airb", floors = c("pd", "lgd", "m"), r_mult = 1,
                    collateral = NULL, secured_share = NULL, claim = NULL) {
  n <- max(length(pd), length(lgd), length(ead), length(asset_class), length(m), length(defaulted), length(elbe),
           length(sales), length(fi), length(collateral), length(secured_share))
  pd <- rep_len(as.double(pd), n); lgd <- rep_len(as.double(lgd), n); ead <- rep_len(as.double(ead), n)
  ac <- rep_len(as.character(asset_class), n)
  bad <- setdiff(unique(ac), .irb_classes)
  if (length(bad) || anyNA(ac)) stop("unknown `asset_class`: ", lst(c(bad, if (anyNA(ac)) "NA")),
                                     ". Use one of ", lst(.irb_classes, 20), ".", call. = FALSE)
  if (anyNA(pd) || any(pd < 0 | pd > 1)) stop("`pd` must be in [0, 1] without missing values.", call. = FALSE)
  if (anyNA(lgd) || any(lgd < 0)) stop("`lgd` must be non-negative without missing values.", call. = FALSE)
  if (anyNA(ead) || any(ead < 0)) stop("`ead` must be non-negative without missing values.", call. = FALSE)
  d <- if (is.null(defaulted)) rep(FALSE, n) else rep_len(isTRUE_vec(defaulted), n)
  e <- if (is.null(elbe)) rep(NA_real_, n) else rep_len(as.double(elbe), n)
  m <- if (is.null(m)) rep(NA_real_, n) else rep_len(as.double(m), n)
  if (any(!is.na(m) & m < 0)) stop("`m` must be non-negative (years).", call. = FALSE)
  firb <- identical(approach, "firb")
  hit <- c(pd_floor = 0L, lgd_floor = 0L, m_floor = 0L, m_cap = 0L)

  # PD floor (defaulted rows: PD = 1)
  pd_used <- pd
  if ("pd" %in% floors) {
    cls <- ac; cls[cls %in% c("corporate_sme", "hvcre")] <- "corporate"
    fl <- .pd_floor_of(params, cls)
    low <- !d & !is.na(fl) & pd_used < fl
    hit[["pd_floor"]] <- sum(low)
    pd_used[low] <- fl[low]
  }
  pd_used[d] <- 1

  # F-IRB: the supervisory LGD of the claim type replaces the caller's value
  lgd_used <- lgd
  if (firb && !is.null(claim)) {
    cl <- rep_len(as.character(claim), n)
    tb <- params$lgd_firb
    badc <- setdiff(unique(cl[!is.na(cl)]), tb$claim)
    if (length(badc)) stop("unknown `claim`: ", lst(badc), ". Use one of ", lst(tb$claim, 10), " (params$lgd_firb).", call. = FALSE)
    sup <- tb$lgd[match(cl, tb$claim)]
    lgd_used[!is.na(sup)] <- sup[!is.na(sup)]
  }
  # LGD floor: own estimates only (A-IRB); the unsecured column unless a
  # collateral column (and optionally a secured share) is given
  if ("lgd" %in% floors && !firb) {
    lf <- params$lgd_floor
    cls <- ac
    cls[cls %in% c("corporate_sme", "hvcre")] <- "corporate"
    cls[cls %in% c("qrre_revolver", "qrre_transactor")] <- "qrre"
    row <- match(cls, lf$asset_class)
    unsec <- lf$unsecured[row]
    col <- if (is.null(collateral)) rep("unsecured", n) else rep_len(as.character(collateral), n)
    badc <- setdiff(unique(col[!is.na(col)]), setdiff(names(lf), "asset_class"))
    if (length(badc)) stop("unknown `collateral`: ", lst(badc), ". Use one of ", lst(setdiff(names(lf), "asset_class")), ".", call. = FALSE)
    sec <- rep(NA_real_, n)
    for (cc in unique(col[!is.na(col)])) {
      ix <- which(col == cc & !is.na(row))
      sec[ix] <- lf[[cc]][row[ix]]
    }
    # mortgages have a real-estate floor only: fall back to it when unsecured is absent
    fb <- is.na(unsec) & !is.na(row)
    unsec[fb] <- lf$real_estate[row[fb]]
    share <- if (is.null(secured_share)) ifelse(col == "unsecured", 0, 1) else rep_len(as.double(secured_share), n)
    share[is.na(share)] <- 0
    fl <- ifelse(is.na(sec), unsec, (1 - share) * ifelse(is.na(unsec), 0, unsec) + share * sec)
    low <- !d & !is.na(fl) & lgd_used < fl
    hit[["lgd_floor"]] <- sum(low)
    lgd_used[low] <- fl[low]
  }

  # correlation, maturity adjustment, K
  r <- .irb_correlation(pd_used, ac, sales, fi, params)
  if (!identical(r_mult, 1)) r <- pmin(0.999, r * r_mult)
  whole <- ac %in% .irb_wholesale
  b <- rep(NA_real_, n); ma <- rep(1, n)
  mm <- rep(NA_real_, n)                     # maturity is reported for wholesale rows only
  if (any(whole)) {
    mm[whole] <- if (firb) params$m_default else m[whole]
    mm[whole & is.na(mm)] <- params$m_default
    if ("m" %in% floors && !firb) {
      lo <- whole & mm < params$m_range[1]; hi <- whole & mm > params$m_range[2]
      hit[["m_floor"]] <- sum(lo); hit[["m_cap"]] <- sum(hi)
      mm[lo] <- params$m_range[1]; mm[hi] <- params$m_range[2]
    }
    pdw <- pmax(pd_used, .irb_b_pd_min)
    b[whole] <- (0.11852 - 0.05478 * log(pdw[whole]))^2
    ma[whole] <- (1 + (mm[whole] - 2.5) * b[whole]) / (1 - 1.5 * b[whole])
  }
  m <- mm
  k <- .irb_k_core(pd_used, lgd_used, r, ma, params)
  if (any(d)) {
    ed <- e; ed[is.na(ed)] <- lgd[is.na(ed)]
    k[d] <- if (firb) 0 else pmax(0, lgd[d] - ed[d])
    r[d] <- NA_real_; b[d] <- NA_real_; ma[d] <- NA_real_
  }
  out <- data.table::data.table(pd_used = pd_used, lgd_used = lgd_used, m = m, r = r, b = b, ma = ma,
                                k = k, rw = 12.5 * k, rwa = 12.5 * k * ead)
  data.table::setattr(out, "floors_hit", data.table::data.table(floor = names(hit), n = unname(hit)))
  out
}

#' IRB risk weight of one or many exposures
#'
#' The asymptotic single risk factor function, vectorised over exposures:
#' PD floors by asset class, LGD input floors for own estimates
#' (`approach = "airb"`; the unsecured column of `params$lgd_floor` unless
#' `collateral` names another column, blended with `secured_share`), the
#' asset correlation of the class (with the firm-size adjustment of
#' `corporate_sme` from `sales` and the multiplier for large or unregulated
#' financial institutions when `fi` is `TRUE`), the maturity adjustment for
#' wholesale classes only (`m` clipped to `params$m_range`, `params$m_default`
#' when missing or under the foundation approach) and
#'
#' \deqn{K = \left[LGD \cdot N\left(\frac{G(PD) + \sqrt{R}\,G(0.999)}{\sqrt{1-R}}\right) - PD \cdot LGD\right] \cdot MA \cdot s}
#'
#' with `s = params$scaling_factor` and, for wholesale classes,
#' \eqn{MA = (1 + (M - 2.5)\,b) / (1 - 1.5\,b)},
#' \eqn{b = (0.11852 - 0.05478 \ln PD)^2} (`MA = 1` for retail). Below
#' `PD = 1e-5`, reachable only without a PD floor (sovereigns), `b` is
#' held at its value at `1e-5`: the regulatory `b` makes `1 - 1.5 b` vanish
#' near `PD = 2.9e-6`, where the adjustment would explode and change sign.
#' Defaulted rows carry `K = max(0, LGD -
#' ELBE)` under `"airb"` and zero under `"firb"`; a missing `elbe` is taken
#' equal to `lgd`. `RW = 12.5 K` and `RWA = RW * ead`.
#'
#' @inheritParams scr_el
#' @param m Effective maturity in years, non-negative (wholesale classes
#'   only, ignored and reported as `NA` on retail rows; `NULL` or `NA` uses
#'   `params$m_default`).
#' @param asset_class One of `"corporate"`, `"corporate_sme"`, `"bank"`,
#'   `"sovereign"`, `"hvcre"`, `"retail_mortgage"`, `"qrre_revolver"`,
#'   `"qrre_transactor"`, `"retail_other"`; a scalar or a vector.
#' @param sales Annual sales of `corporate_sme` obligors, in the unit of
#'   `params$correlation$sme`, clipped to its bounds. A missing value takes
#'   the upper bound, i.e. no firm-size adjustment: the adjustment requires
#'   reported sales (Basel Framework CRE31.9; CRR Article 153(4)).
#' @param fi Logical: regulated financial institution above the size
#'   threshold, or unregulated one (correlation multiplier).
#' @param params An [scr_irb_params()] object.
#' @param approach `"airb"` (own LGD, floored) or `"firb"` (supervisory LGD
#'   supplied by the caller, no LGD floor, maturity fixed).
#' @param apply_floors `TRUE` (all input floors), `FALSE` (none) or a subset
#'   of `c("pd", "lgd", "m")`.
#' @param collateral Optional column of `params$lgd_floor` naming the
#'   collateral type of each exposure (`"financial"`, `"receivables"`,
#'   `"real_estate"`, `"other_physical"`); `NULL` means unsecured.
#' @param secured_share Optional secured share in `[0, 1]` blending the
#'   unsecured and the collateral floors.
#' @param claim Under `"firb"`, an optional claim type per exposure naming a
#'   row of `params$lgd_firb` (for example `"senior_unsecured"` or
#'   `"subordinated"`); the supervisory LGD of that row replaces `lgd`.
#'   `NULL` keeps the caller's `lgd`. The row names differ by preset
#'   (`"senior_unsecured"` under `"bcb"`, `"senior_unsecured_corporate"` and
#'   `"senior_unsecured_fi"` otherwise); see `params$lgd_firb`.
#'
#' @return A `data.table` with one row per exposure: `pd_used`, `lgd_used`
#'   (after floors; PD one on defaulted rows), `m` (after clipping;
#'   `params$m_default` under `"firb"`; `NA` on retail rows), `r`,
#'   `b`, `ma`, `k`, `rw`, `rwa`; attribute `floors_hit` counts the rows
#'   where each floor was binding.
#'
#' @references
#' Basel Committee on Banking Supervision (2023). *The Basel Framework*,
#' CRE31 (IRB approach: risk-weight functions) and CRE32 (risk components).
#' BCBS (2005). *An explanatory note on the Basel II IRB risk weight
#' functions*.
#'
#' @family irb-capital
#' @examples
#' scr_irb_rw(0.01, 0.45, m = 2.5, asset_class = "corporate")
#' scr_irb_rw(c(0.01, 0.02), c(0.20, 0.80), asset_class = c("retail_mortgage", "qrre_revolver"))
#' r <- scr_irb_rw(1e-4, 0.5, asset_class = "retail_other")
#' attr(r, "floors_hit")
#' # foundation approach: the supervisory LGD of the claim type
#' scr_irb_rw(0.01, lgd = 0, m = 2.5, asset_class = "corporate", approach = "firb",
#'            claim = "senior_unsecured")
#' @export
scr_irb_rw <- function(pd, lgd, ead = 1, m = NULL, asset_class, sales = NULL, fi = FALSE, defaulted = NULL,
                       elbe = NULL, params = scr_irb_params("bcb"), approach = c("airb", "firb"),
                       apply_floors = TRUE, collateral = NULL, secured_share = NULL, claim = NULL) {
  approach <- match.arg(approach)
  params <- .check_params(params, "scr_irb_rw")
  floors <- .floors_of(apply_floors)
  .irb_rw(pd, lgd, ead, m, asset_class, sales, fi, defaulted, elbe, params, approach, floors,
          collateral = collateral, secured_share = secured_share, claim = claim)
}

#' @keywords internal
#' @noRd
.floors_of <- function(apply_floors) {
  all <- c("pd", "lgd", "m")
  if (is.logical(apply_floors)) return(if (isTRUE(apply_floors)) all else character())
  bad <- setdiff(apply_floors, all)
  if (length(bad)) stop("`apply_floors` must be TRUE, FALSE or a subset of ", lst(all), ".", call. = FALSE)
  as.character(apply_floors)
}

# -- Standardised risk weight ------------------------------------------------ #

#' Standardised risk weight of an exposure
#'
#' Lookup in `params$sa_rw`: regulatory retail (`"retail_other"`,
#' `"qrre_*"`: 75 %, or the transactor weight), residential mortgages by
#' loan-to-value band (a missing LTV takes the highest band), corporates by
#' external rating bucket (`"AAA"` to `"AA-"`, `"A"`, `"BBB"`, `"BB"`, below;
#' `NA` is unrated; `"IG"` marks an unrated investment-grade obligor where
#' ratings are not used) or the SME weight, banks and sovereigns through
#' the corporate rating rows, and defaulted exposures by the specific
#' provision ratio (or the mortgage row). Arguments are recycled.
#'
#' @inheritParams scr_irb_rw
#' @param ltv Loan-to-value at origination, decimal (mortgages).
#' @param rating External rating string (corporates), `NA` when unrated.
#' @param transactor Logical: revolving facility repaid in full every month.
#' @param provision_ratio Specific provisions over the outstanding amount
#'   (defaulted rows).
#' @param sme Logical: corporate small or medium enterprise (also implied by
#'   `asset_class = "corporate_sme"`).
#' @param granular Logical (scalar or per exposure): whether the retail
#'   exposure belongs to a granular regulatory retail pool; `FALSE` applies
#'   the non-granular retail weight.
#'
#' @return A numeric vector of standardised risk weights (decimals).
#'
#' @references
#' Basel Committee on Banking Supervision (2023). *The Basel Framework*,
#' CRE20 (standardised approach: individual exposures).
#'
#' @family irb-capital
#' @examples
#' scr_sa_rw(c("retail_other", "retail_mortgage", "corporate"), ltv = c(NA, 0.55, NA),
#'           rating = c(NA, NA, "A+"))
#' scr_sa_rw("retail_other", defaulted = TRUE, provision_ratio = c(0.1, 0.3))
#' @export
scr_sa_rw <- function(asset_class, ltv = NULL, rating = NULL, transactor = NULL, defaulted = NULL,
                      provision_ratio = NULL, sme = NULL, granular = TRUE, params = scr_irb_params("bcb")) {
  params <- .check_params(params, "scr_sa_rw")
  n <- max(length(asset_class), length(ltv), length(rating), length(transactor), length(defaulted),
           length(provision_ratio), length(sme), length(granular))
  gr <- rep_len(isTRUE_vec(granular), n)
  ac <- rep_len(as.character(asset_class), n)
  bad <- setdiff(unique(ac), .irb_classes)
  if (length(bad) || anyNA(ac)) stop("unknown `asset_class`: ", lst(c(bad, if (anyNA(ac)) "NA")), ".", call. = FALSE)
  ltv <- if (is.null(ltv)) rep(NA_real_, n) else rep_len(as.double(ltv), n)
  rt <- if (is.null(rating)) rep(NA_character_, n) else rep_len(as.character(rating), n)
  tr <- if (is.null(transactor)) rep(FALSE, n) else rep_len(isTRUE_vec(transactor), n)
  d <- if (is.null(defaulted)) rep(FALSE, n) else rep_len(isTRUE_vec(defaulted), n)
  pr <- if (is.null(provision_ratio)) rep(NA_real_, n) else rep_len(as.double(provision_ratio), n)
  sm <- (if (is.null(sme)) rep(FALSE, n) else rep_len(isTRUE_vec(sme), n)) | ac == "corporate_sme"
  tab <- params$sa_rw
  look <- function(cls, sub) tab$rw[match(paste(cls, sub), paste(tab$asset_class, tab$sub_class))]
  rw <- rep(NA_real_, n)

  ret <- ac %in% c("retail_other", "qrre_revolver", "qrre_transactor")
  rw[ret] <- look("retail_other", ifelse(!gr[ret], "non_granular",
                                          ifelse(tr[ret] | ac[ret] == "qrre_transactor", "transactor", "standard")))
  mo <- ac == "retail_mortgage"
  if (any(mo)) {
    bands <- tab[tab$asset_class == "retail_mortgage" & tab$sub_class == "standard", ]
    l <- ltv[mo]; l[is.na(l)] <- Inf
    # first band with lo < ltv <= hi (loop over the few bands, not the rows);
    # ltv <= 0 takes the first band, anything unmatched the last
    idx <- rep(NA_integer_, length(l))
    for (i in seq_len(nrow(bands))) idx[is.na(idx) & l > bands$ltv_lo[i] & l <= bands$ltv_hi[i]] <- i
    idx[is.na(idx)] <- ifelse(l[is.na(idx)] <= 0, 1L, nrow(bands))
    rw[mo] <- bands$rw[idx]
  }
  wh <- ac %in% .irb_wholesale
  if (any(wh)) {
    sub <- .sa_rating_bucket(rt[wh])
    sub[sm[wh] & sub == "unrated"] <- "sme"
    rw[wh] <- look("corporate", sub)
  }
  if (any(d)) {
    sub <- ifelse(ac[d] == "retail_mortgage", "mortgage", ifelse(!is.na(pr[d]) & pr[d] >= 0.2, "provision_ge_20", "provision_lt_20"))
    rw[d] <- look("defaulted", sub)
  }
  rw
}

#' @keywords internal
#' @noRd
.sa_rating_bucket <- function(rating) {
  r <- toupper(trimws(as.character(rating)))
  r[is.na(r) | !nzchar(r)] <- "UNRATED"
  core <- gsub("[+-]", "", r)
  ifelse(core %in% c("AAA", "AA"), "rated_aaa_aa",
  ifelse(core == "A", "rated_a",
  ifelse(core == "BBB", "rated_bbb",
  ifelse(core == "BB", "rated_bb",
  ifelse(core %in% c("B", "CCC", "CC", "C", "D", "SD", "RD"), "rated_below_bb",
  ifelse(core %in% c("IG", "INVESTMENT_GRADE"), "investment_grade_unrated", "unrated"))))))
}

# -- Portfolio capital ------------------------------------------------------- #

#' Expected loss, risk-weighted assets and capital of a portfolio
#'
#' Runs [scr_irb_rw()] on every exposure, aggregates by segment, compares
#' the IRB result with the standardised approach for the output floor,
#' reconciles regulatory expected loss with the provision stock
#' (shortfall deducted from capital; excess eligible as tier 2 up to 0.6 %
#' of the IRB risk-weighted assets), measures the impact of each input
#' floor, runs a fixed sensitivity grid and reports the name concentration
#' of the book. The parameter tables come from `params`; the object records
#' whether they were edited.
#'
#' @section Inputs:
#'
#' `x` is either a table of exposures, the remaining arguments naming its
#' columns, or a list `list(pd = , lgd = , ead = , data = )` whose elements
#' are fitted models with an [scr_apply()] method (the PD, LGD and EAD
#' objects of the IRB modules) and `data` the table to apply them to. In
#' the list form each model present fills the corresponding vector from the
#' columns `pd_final`, `lgd_final` and `ead_predicted` of its `scr_apply()`
#' output, and the provenance is written to the ledger; elements that are
#' `NULL` fall back to the named columns of `data`.
#'
#' `asset_class` is a column name of `x` or a single class applied to every
#' row. Segment means are weighted by EAD. The sensitivity grid shocks the
#' PD (x1.10, x1.25, x1.50), the LGD (+5 percentage points), the EAD
#' (+10 %), removes the input floors, scales the correlation (x1.25) and
#' stresses the PD with the one-factor model at `q = 0.95` and `0.99`
#' ([scr_pd_stress()], the stressed PD then re-entering the function so
#' that the correlation follows it).
#'
#' @param x A table of exposures (`data.frame` or `data.table`) or the list
#'   form described above.
#' @param pd,lgd,ead Column names of the probability of default, loss given
#'   default and exposure at default.
#' @param segment Optional column name of the reporting segment.
#' @param asset_class A column name or a single asset class (see
#'   [scr_irb_rw()]).
#' @param m,defaulted,elbe,provisions,ltv,rating,sales,fi,transactor,grade,id
#'   Optional column names: effective maturity, default flag, best estimate
#'   of expected loss, provision stock, loan-to-value, external rating,
#'   annual sales, financial-institution flag, transactor flag, PD grade
#'   (defines the SQL pools together with `segment`) and exposure
#'   identifier.
#' @param claim Optional column name: the claim type of each exposure under
#'   the foundation approach (a row of `params$lgd_firb`); the supervisory
#'   LGD then replaces `lgd`.
#' @param granular `TRUE`, `FALSE` or a column name: whether the retail
#'   exposures belong to a granular regulatory retail pool (the standardised
#'   comparison uses the non-granular weight otherwise).
#' @param params An [scr_irb_params()] object; defaults to the preset of
#'   `config$framework`.
#' @param config An [scr_config()] object (`capital_approach`,
#'   `capital_target_ratio`, `capital_output_floor`, `capital_sensitivity`,
#'   `nthread`, `verbose`).
#' @param keep_rows Keep the per-exposure table in the object.
#'
#' @return An object of class `scr_capital`: a list with `exposures`
#'   (per-exposure table, only with `keep_rows = TRUE`), `segments` (the
#'   reconciliation table: `segment`, `n`, `ead`, `pd_mean`, `lgd_mean`,
#'   `m`, `r_mean`, `k_mean`, `rw`, `rwa_irb`, `rwa_sa`, `irb_sa_ratio`,
#'   `el`, `provisions`, `shortfall_excess`), `pools` (one row per
#'   segment and grade with the constants the SQL emits), `totals` (`n`,
#'   `ead`, `el`, `rwa_irb`, `rwa_sa`, `irb_sa_ratio`, `output_floor`,
#'   `rwa_floor`, `rwa_reported`, `floor_binding`, `headroom`, `density`,
#'   `target_ratio`, `capital`, `provisions`, `shortfall`, `excess`,
#'   `tier2_addback`, `tier2_cap`, `hhi`, `n_eff`, `max_share`,
#'   `granular`), `floors` (`floor`, `n_hit`, `ead_hit`, `delta_rwa`),
#'   `sensitivity` (`shock`, `rwa`, `delta`, `delta_pct`),
#'   `concentration` (share of EAD and RWA by segment), `framework`,
#'   `approach`, `params`, `config`, `ledger`, `model_card` and, after
#'   [scr_export()], `files`.
#'   `segments` and `totals` also carry `n_defaulted`; `totals` also
#'   `el_rate` and `rwa_irb_no_floors`; `concentration` has `segment`, `n`,
#'   `ead`, `rwa`, `ead_share`, `rwa_share` and `hhi_contribution`; `columns`
#'   records the column names the SQL reads.
#'
#' @references
#' Basel Committee on Banking Supervision (2023). *The Basel Framework*,
#' CRE31, CRE35 (treatment of expected losses and provisions), RBC20
#' (output floor).
#'
#' @family irb-capital
#' @examples
#' cfg <- scr_config(verbose = FALSE)
#' cap <- scr_capital(scr_demo_portfolio, segment = "segment", asset_class = "asset_class",
#'                    m = "m", defaulted = "defaulted", elbe = "elbe", provisions = "provision",
#'                    ltv = "ltv", rating = "rating", sales = "sales", transactor = "transactor",
#'                    grade = "grade", id = "id", config = cfg)
#' cap
#' cap$segments[, c("segment", "n", "rw", "irb_sa_ratio")]
#' cap$floors
#' @export
scr_capital <- function(x, pd = "pd", lgd = "lgd", ead = "ead", segment = NULL, asset_class = config$asset_class,
                        m = NULL, defaulted = NULL, elbe = NULL, provisions = NULL, ltv = NULL, rating = NULL,
                        sales = NULL, fi = NULL, transactor = NULL, grade = NULL, id = NULL,
                        claim = NULL, granular = TRUE,
                        params = scr_irb_params(config$framework), config = scr_config(), keep_rows = FALSE) {
  check_config(config, "scr_capital")
  cfg <- config
  old <- scr_verbose(isTRUE(cfg$verbose)); on.exit(scr_verbose(old), add = TRUE)
  params <- .check_params(params, "scr_capital")
  approach <- cfg$capital_approach
  t0 <- Sys.time()
  pd_used <- lgd_used <- rwa <- rwa_sa <- el <- provision <- k <- r <- NULL
  irb_sa_ratio <- shortfall_excess <- rwa_irb <- ead_share <- rwa_share <- hhi_contribution <- rw_sa <- NULL

  # -- resolve the inputs --------------------------------------------------- #
  prov <- list(pd = "column", lgd = "column", ead = "column")
  if (is.data.frame(x)) {
    dt <- data.table::as.data.table(x)
  } else if (is.list(x)) {
    if (is.null(x$data)) stop("scr_capital(): the list form needs `data`, the table to apply the models to.", call. = FALSE)
    dt <- data.table::as.data.table(x$data)
    want <- c(pd = "pd_final", lgd = "lgd_final", ead = "ead_predicted")
    for (nm in names(want)) {
      obj <- x[[nm]]
      if (is.null(obj)) next
      ap <- data.table::as.data.table(scr_apply(obj, dt))
      if (!want[[nm]] %in% names(ap)) {
        stop(sprintf("scr_capital(): scr_apply() of `%s` (class %s) returned no column `%s`.", nm, class(obj)[1], want[[nm]]), call. = FALSE)
      }
      if (nrow(ap) != nrow(dt)) stop(sprintf("scr_capital(): scr_apply() of `%s` returned %d rows for %d exposures.", nm, nrow(ap), nrow(dt)), call. = FALSE)
      data.table::set(dt, j = paste0(".", nm), value = ap[[want[[nm]]]])
      assign(nm, paste0(".", nm))
      prov[[nm]] <- sprintf("scr_apply(%s)$%s", class(obj)[1], want[[nm]])
    }
  } else stop("scr_capital(): `x` must be a table of exposures or a list with `data`.", call. = FALSE)
  if (!nrow(dt)) stop("scr_capital(): no exposure.", call. = FALSE)

  col <- function(nm, what) {
    if (is.null(nm)) return(NULL)
    if (!is.character(nm) || length(nm) != 1L) stop(sprintf("scr_capital(): `%s` must be a column name.", what), call. = FALSE)
    if (!nm %in% names(dt)) stop(sprintf("scr_capital(): column `%s` (%s) not found.", nm, what), call. = FALSE)
    dt[[nm]]
  }
  v_pd <- col(pd, "pd"); v_lgd <- col(lgd, "lgd"); v_ead <- col(ead, "ead")
  v_claim <- col(claim, "claim")
  v_gr <- if (is.character(granular)) col(granular, "granular") else rep_len(isTRUE_vec(granular), nrow(dt))
  ac_src <- if (length(asset_class) == 1L && asset_class %in% names(dt)) "column"
            else if (all(asset_class %in% .irb_classes) && length(asset_class) == 1L) "constant"
            else stop("scr_capital(): `asset_class` must be a column of `x` or one of ", lst(.irb_classes, 20), ".", call. = FALSE)
  v_ac <- if (identical(ac_src, "column")) as.character(dt[[asset_class]]) else rep(asset_class, nrow(dt))
  v_m <- col(m, "m"); v_def <- col(defaulted, "defaulted"); v_elbe <- col(elbe, "elbe")
  v_prov <- col(provisions, "provisions"); v_ltv <- col(ltv, "ltv"); v_rat <- col(rating, "rating")
  v_sales <- col(sales, "sales"); v_fi <- col(fi, "fi"); v_tr <- col(transactor, "transactor")
  v_grade <- col(grade, "grade"); v_id <- col(id, "id")
  v_seg <- if (is.null(segment)) rep("portfolio", nrow(dt)) else as.character(col(segment, "segment"))
  v_seg[is.na(v_seg)] <- "NA"
  n <- nrow(dt)
  d <- if (is.null(v_def)) rep(FALSE, n) else isTRUE_vec(v_def)
  msg("  capital: %s exposures | %s | %s | %s", n_fmt(n), params$framework, approach,
      if (identical(ac_src, "column")) sprintf("%d asset classes", length(unique(v_ac))) else asset_class)

  # -- IRB, floors, SA, EL -------------------------------------------------- #
  run <- function(pd_v = v_pd, lgd_v = v_lgd, ead_v = v_ead, floors = c("pd", "lgd", "m"), r_mult = 1) {
    .irb_rw(pd_v, lgd_v, ead_v, v_m, v_ac, v_sales, v_fi %||% FALSE, v_def, v_elbe, params, approach, floors, r_mult,
            claim = v_claim)
  }
  irb <- run()
  fh <- attr(irb, "floors_hit")
  floors <- data.table::rbindlist(lapply(c("pd", "lgd", "m"), function(f) {
    without <- run(floors = setdiff(c("pd", "lgd", "m"), f))
    binding <- switch(f, pd = irb$pd_used != without$pd_used, lgd = irb$lgd_used != without$lgd_used,
                      m = !is.na(irb$m) & !is.na(without$m) & irb$m != without$m)
    data.table::data.table(floor = paste0(f, "_floor"), n_hit = sum(binding), ead_hit = sum(v_ead[binding]),
                           delta_rwa = sum(irb$rwa) - sum(without$rwa))
  }))
  no_floor <- run(floors = character())
  el_v <- scr_el(irb$pd_used, irb$lgd_used, v_ead, d, v_elbe)
  use_sa <- isTRUE(cfg$capital_output_floor)
  rw_sa <- if (use_sa) {
    pr <- if (is.null(v_prov)) NULL else ifelse(v_ead > 0, v_prov / v_ead, NA_real_)
    scr_sa_rw(v_ac, v_ltv, v_rat, v_tr, d, pr, v_ac == "corporate_sme", granular = v_gr, params = params)
  } else rep(NA_real_, n)
  ex <- data.table::data.table(id = if (is.null(v_id)) seq_len(n) else v_id, segment = v_seg,
                               grade = if (is.null(v_grade)) NA_character_ else as.character(v_grade),
                               asset_class = v_ac, pd = as.double(v_pd), lgd = as.double(v_lgd), ead = as.double(v_ead),
                               defaulted = as.integer(d), elbe = if (is.null(v_elbe)) NA_real_ else as.double(v_elbe),
                               provision = if (is.null(v_prov)) NA_real_ else as.double(v_prov))
  ex <- cbind(ex, irb[, c("pd_used", "lgd_used", "m", "r", "b", "ma", "k", "rw", "rwa")])
  ex[, `:=`(rw_sa = rw_sa, rwa_sa = rw_sa * ead, el = el_v)]

  # -- segments, pools, totals ---------------------------------------------- #
  wmean <- function(v, w) { ok <- !is.na(v); if (!any(ok)) NA_real_ else if (sum(w[ok]) > 0) sum(v[ok] * w[ok]) / sum(w[ok]) else mean(v[ok]) }
  segs <- ex[, list(n = .N, ead = sum(ead), pd_mean = wmean(pd_used, ead), lgd_mean = wmean(lgd_used, ead),
                    m = wmean(m, ead), r_mean = wmean(r, ead), k_mean = wmean(k, ead),
                    rw = if (sum(ead) > 0) sum(rwa) / sum(ead) else NA_real_, rwa_irb = sum(rwa),
                    rwa_sa = if (use_sa) sum(rwa_sa) else NA_real_, el = sum(el),
                    provisions = if (is.null(v_prov)) NA_real_ else sum(provision, na.rm = TRUE),
                    n_defaulted = sum(defaulted)), by = "segment"]
  segs[, irb_sa_ratio := ifelse(!is.na(rwa_sa) & rwa_sa > 0, rwa_irb / rwa_sa, NA_real_)]
  segs[, shortfall_excess := provisions - el]
  data.table::setcolorder(segs, c("segment", "n", "ead", "pd_mean", "lgd_mean", "m", "r_mean", "k_mean", "rw",
                                  "rwa_irb", "rwa_sa", "irb_sa_ratio", "el", "provisions", "shortfall_excess", "n_defaulted"))
  data.table::setorder(segs, -rwa_irb)
  pools <- .capital_pools(ex, has_grade = !is.null(v_grade))

  tot_ead <- sum(ex$ead); rwa_irb <- sum(ex$rwa); rwa_sa_t <- if (use_sa) sum(ex$rwa_sa) else NA_real_
  of <- if (use_sa) params$output_floor else NA_real_
  rwa_floor <- if (use_sa) of * rwa_sa_t else NA_real_
  rwa_rep <- if (use_sa) max(rwa_irb, rwa_floor) else rwa_irb
  prov_t <- if (is.null(v_prov)) NA_real_ else sum(v_prov, na.rm = TRUE)
  el_t <- sum(ex$el)
  shortfall <- if (is.na(prov_t)) NA_real_ else max(0, el_t - prov_t)
  excess <- if (is.na(prov_t)) NA_real_ else max(0, prov_t - el_t)
  t2cap <- 0.006 * rwa_irb
  share <- if (tot_ead > 0) ex$ead / tot_ead else rep(0, n)
  hhi <- sum(share^2)
  totals <- list(
    n = n, n_defaulted = sum(d), ead = tot_ead, el = el_t, el_rate = if (tot_ead > 0) el_t / tot_ead else NA_real_,
    rwa_irb = rwa_irb, rwa_irb_no_floors = sum(no_floor$rwa), rwa_sa = rwa_sa_t,
    irb_sa_ratio = if (use_sa && rwa_sa_t > 0) rwa_irb / rwa_sa_t else NA_real_,
    output_floor = of, rwa_floor = rwa_floor, rwa_reported = rwa_rep,
    floor_binding = if (use_sa) rwa_floor > rwa_irb else FALSE,
    headroom = if (use_sa) rwa_irb - rwa_floor else NA_real_,
    density = if (tot_ead > 0) rwa_rep / tot_ead else NA_real_,
    target_ratio = cfg$capital_target_ratio, capital = cfg$capital_target_ratio * rwa_rep,
    provisions = prov_t, shortfall = shortfall, excess = excess,
    tier2_addback = if (is.na(excess)) NA_real_ else min(excess, t2cap), tier2_cap = t2cap,
    hhi = hhi, n_eff = if (hhi > 0) 1 / hhi else NA_real_, max_share = if (n) max(share) else NA_real_,
    granular = if (n) max(share) <= 0.002 else NA)
  conc <- ex[, list(n = .N, ead = sum(ead), rwa = sum(rwa)), by = "segment"]
  conc[, `:=`(ead_share = if (tot_ead > 0) ead / tot_ead else NA_real_,
              rwa_share = if (rwa_irb > 0) rwa / rwa_irb else NA_real_)]
  conc[, hhi_contribution := ead_share^2]
  data.table::setorder(conc, -ead)

  # -- sensitivity ----------------------------------------------------------- #
  sens <- NULL
  if (isTRUE(cfg$capital_sensitivity)) {
    base_pd <- irb$pd_used; base_r <- irb$r
    shocks <- list(
      base = list(), pd_x1.10 = list(pd = pmin(1, v_pd * 1.10)), pd_x1.25 = list(pd = pmin(1, v_pd * 1.25)),
      pd_x1.50 = list(pd = pmin(1, v_pd * 1.50)), lgd_plus_5pp = list(lgd = v_lgd + 0.05),
      ead_x1.10 = list(ead = v_ead * 1.10), no_floor = list(floors = character()), r_x1.25 = list(r_mult = 1.25),
      vasicek_q0.95 = list(pd = ifelse(d, v_pd, scr_pd_stress(base_pd, ifelse(is.na(base_r), 0, base_r), 0.95))),
      vasicek_q0.99 = list(pd = ifelse(d, v_pd, scr_pd_stress(base_pd, ifelse(is.na(base_r), 0, base_r), 0.99))))
    vals <- .scr_lapply(shocks, function(s) {
      sum(run(pd_v = s$pd %||% v_pd, lgd_v = s$lgd %||% v_lgd, ead_v = s$ead %||% v_ead,
              floors = if (is.null(s$floors)) c("pd", "lgd", "m") else s$floors, r_mult = s$r_mult %||% 1)$rwa)
    }, nthread = cfg$nthread)
    rw_v <- unlist(vals)
    sens <- data.table::data.table(shock = names(shocks), rwa = rw_v, delta = rw_v - rw_v[1],
                                   delta_pct = if (rw_v[1] > 0) (rw_v - rw_v[1]) / rw_v[1] else NA_real_)
  }

  # -- ledger and model card -------------------------------------------------- #
  n_m_default <- if (is.null(v_m)) sum(v_ac %in% .irb_wholesale) else sum(v_ac %in% .irb_wholesale & is.na(v_m))
  n_elbe_missing <- if (is.null(v_elbe)) sum(d) else sum(d & is.na(v_elbe))
  ledger <- data.table::rbindlist(list(
    .cap_ledger("framework", sprintf("%s | approach %s | params %s", params$framework, approach,
                                     if (isTRUE(params$modified)) "modified by the user" else "preset")),
    .cap_ledger("inputs", sprintf("pd: %s | lgd: %s | ead: %s | asset_class: %s", prov$pd, prov$lgd, prov$ead,
                                  if (identical(ac_src, "column")) sprintf("column %s", asset_class) else sprintf("constant %s", asset_class))),
    .cap_ledger("floors", sprintf("pd floor on %d rows, lgd floor on %d rows, maturity clipped on %d rows",
                                  fh$n[fh$floor == "pd_floor"], fh$n[fh$floor == "lgd_floor"],
                                  fh$n[fh$floor == "m_floor"] + fh$n[fh$floor == "m_cap"])),
    if (n_m_default > 0) .cap_ledger("maturity", sprintf("M = %.1f (params$m_default) on %d wholesale rows without maturity", params$m_default, n_m_default)),
    if (n_elbe_missing > 0) .cap_ledger("elbe", sprintf("ELBE taken equal to LGD (K = 0) on %d defaulted rows without elbe", n_elbe_missing)),
    .cap_ledger("output_floor", if (use_sa) sprintf("%s of the standardised RWA: %s", fmt_pct(of, 1),
                                                    if (totals$floor_binding) "BINDING" else "not binding") else "not computed (capital_output_floor = FALSE)"),
    .cap_ledger("provisions", if (is.na(prov_t)) "no provision column: EL comparison skipped"
                              else sprintf("EL %s vs provisions %s: shortfall %s, excess %s, tier 2 add-back %s",
                                           format(round(el_t)), format(round(prov_t)), format(round(shortfall)),
                                           format(round(excess)), format(round(totals$tier2_addback))))
  ), use.names = TRUE, fill = TRUE)
  model_card <- list(
    package = sprintf("scorecraft %s", as.character(utils::packageVersion("scorecraft"))),
    fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    framework = params$framework, approach = approach, params_modified = isTRUE(params$modified),
    asset_classes = paste(sort(unique(v_ac)), collapse = ", "),
    pd_source = prov$pd, lgd_source = prov$lgd, ead_source = prov$ead,
    n_exposures = n, n_defaulted = sum(d), n_segments = nrow(segs), n_pools = nrow(pools),
    ead = tot_ead, el = el_t, rwa_irb = rwa_irb, rwa_sa = rwa_sa_t, rwa_reported = rwa_rep,
    output_floor = of, output_floor_binding = totals$floor_binding, density = totals$density,
    target_ratio = cfg$capital_target_ratio, capital = totals$capital,
    shortfall = shortfall, excess = excess, hhi = hhi, n_eff = totals$n_eff,
    scaling_factor = params$scaling_factor, confidence = params$confidence)
  msg("  capital: EAD %s | RWA %s (density %s) | EL %s | floors hit %d (%.2fs)", format(round(tot_ead)),
      format(round(rwa_rep)), fmt_pct(totals$density, 1), format(round(el_t)), sum(floors$n_hit),
      as.numeric(difftime(Sys.time(), t0, units = "secs")))

  structure(list(
    exposures = if (isTRUE(keep_rows)) ex[] else NULL, segments = segs[], pools = pools, totals = totals,
    floors = floors[], sensitivity = sens, concentration = conc[],
    framework = params$framework, approach = approach, params = params, config = cfg,
    columns = list(pd = pd, lgd = lgd, ead = ead, segment = segment, grade = grade, id = id, defaulted = defaulted,
                   elbe = elbe, asset_class = if (identical(ac_src, "column")) asset_class else NULL),
    ledger = ledger, model_card = model_card, files = NULL
  ), class = c("scr_capital", "list"))
}

#' @keywords internal
#' @noRd
.cap_ledger <- function(action, detail, reason = NA_character_) {
  data.table::data.table(action = action, detail = detail, reason = reason, date = format(Sys.Date()))
}

#' Pool constants for the SQL: one row per segment (and grade)
#'
#' A pool is homogeneous when its performing rows share one PD, one LGD and
#' one K. Otherwise the constants are the EAD-weighted PD and the
#' EL-consistent LGD (`sum(pd lgd ead) / sum(pd ead)`) with the EAD-weighted
#' K, so that EL and RWA are exact at pool level and approximate per row.
#' @keywords internal
#' @noRd
.capital_pools <- function(ex, has_grade) {
  defaulted <- ead <- pd_used <- lgd_used <- k <- r <- b <- ma <- NULL
  by <- if (has_grade) c("segment", "grade") else "segment"
  perf <- ex[defaulted == 0L]
  wm <- function(v, w) if (data.table::uniqueN(v) == 1L) v[1] else if (all(is.na(v))) NA_real_ else if (sum(w) > 0) sum(v * w) / sum(w) else mean(v)
  p <- perf[, list(n = .N, ead = sum(ead),
                   pd = wm(pd_used, ead),
                   lgd = if (data.table::uniqueN(lgd_used) == 1L) lgd_used[1]
                         else if (sum(pd_used * ead) > 0) sum(pd_used * lgd_used * ead) / sum(pd_used * ead) else mean(lgd_used),
                   r = wm(r, ead), b = wm(b, ead), ma = wm(ma, ead), k = wm(k, ead),
                   homogeneous = data.table::uniqueN(pd_used) == 1L && data.table::uniqueN(lgd_used) == 1L &&
                                 data.table::uniqueN(k) == 1L), by = by]
  keys <- unique(ex[, by, with = FALSE])
  p <- merge(keys, p, by = by, all.x = TRUE)
  p[is.na(n), `:=`(n = 0L, ead = 0, homogeneous = TRUE)]
  p[, rw := 12.5 * k]
  data.table::setorderv(p, by)
  p[]
}

#' @export
print.scr_capital <- function(x, ...) {
  t <- x$totals
  cat(sprintf("<scr_capital> %s | %s | %s exposures in %d segments%s\n", x$framework, x$approach, n_fmt(t$n),
              nrow(x$segments), if (isTRUE(x$params$modified)) " | params modified" else ""))
  cat(sprintf("  EAD %s | EL %s (%s) | RWA IRB %s | density %s | capital (%s) %s\n",
              format(round(t$ead), big.mark = ","), format(round(t$el), big.mark = ","), fmt_pct(t$el_rate, 2),
              format(round(t$rwa_irb), big.mark = ","), fmt_pct(t$density, 1), fmt_pct(t$target_ratio, 1),
              format(round(t$capital), big.mark = ",")))
  if (!is.na(t$rwa_sa)) {
    cat(sprintf("  standardised RWA %s | IRB/SA %.3f | output floor %s: %s (headroom %s)\n",
                format(round(t$rwa_sa), big.mark = ","), t$irb_sa_ratio, fmt_pct(t$output_floor, 1),
                if (isTRUE(t$floor_binding)) "BINDING" else "not binding", format(round(t$headroom), big.mark = ",")))
  }
  if (!is.na(t$provisions)) {
    cat(sprintf("  provisions %s vs EL: shortfall %s | excess %s | tier 2 add-back %s (cap %s)\n",
                format(round(t$provisions), big.mark = ","), format(round(t$shortfall), big.mark = ","),
                format(round(t$excess), big.mark = ","), format(round(t$tier2_addback), big.mark = ","),
                format(round(t$tier2_cap), big.mark = ",")))
  }
  f <- x$floors
  cat(sprintf("  floors: %s | HHI %.5f (n_eff %s, max share %s)\n",
              paste(sprintf("%s %d rows, RWA %+s", sub("_floor$", "", f$floor), f$n_hit, format(round(f$delta_rwa), big.mark = ",")), collapse = " | "),
              t$hhi, format(round(t$n_eff)), fmt_pct(t$max_share, 2)))
  s <- utils::head(x$segments, 5)
  cat("  top segments by RWA:\n")
  for (i in seq_len(nrow(s))) {
    cat(sprintf("    %-22s n %-7s EAD %-14s RW %6s  RWA %-14s%s\n", substr(s$segment[i], 1, 22), n_fmt(s$n[i]),
                format(round(s$ead[i]), big.mark = ","), fmt_pct(s$rw[i], 1), format(round(s$rwa_irb[i]), big.mark = ","),
                if (is.na(s$irb_sa_ratio[i])) "" else sprintf("  IRB/SA %.2f", s$irb_sa_ratio[i])))
  }
  if (!is.null(x$sensitivity)) {
    ss <- x$sensitivity[x$sensitivity$shock != "base", ]
    ss <- ss[order(-abs(ss$delta)), ][seq_len(min(3L, nrow(ss))), ]
    cat(sprintf("  sensitivity: %s\n", paste(sprintf("%s %+.1f%%", ss$shock, 100 * ss$delta_pct), collapse = " | ")))
  }
  invisible(x)
}

# -- SQL --------------------------------------------------------------------- #

#' @rdname scr_sql
#' @param level For `scr_capital`: `"exposure"` (default, one row per
#'   exposure with `el`, `k`, `rw`, `rwa`) or `"portfolio"` (the aggregate
#'   by segment).
#' @export
scr_sql.scr_capital <- function(x, table = NULL, dialect = NULL, file = NULL, level = c("exposure", "portfolio"), ...) {
  .scr_local_verbose(x)
  level <- match.arg(level)
  cfg <- x$config
  if (!is.null(table)) cfg$sql_table <- table
  if (!is.null(dialect)) cfg$sql_dialect <- dialect
  dialect <- match.arg(cfg$sql_dialect, c("ansi", "databricks", "spark", "hive", "mysql", "mariadb", "sqlserver",
                                          "bigquery", "postgres", "oracle", "snowflake", "redshift", "duckdb", "sqlite"))
  cols <- x$columns
  p <- x$pools
  has_grade <- !is.null(cols$grade)
  has_seg <- !is.null(cols$segment)
  from_dual <- if (identical(dialect, "oracle")) " FROM DUAL" else ""
  lit <- function(v) if (is.numeric(v)) .sql_num(v) else .sql_str(as.character(v))
  seg_vals <- if (has_seg) lit(p$segment) else NULL
  rows <- vapply(seq_len(nrow(p)), function(i) {
    parts <- c(if (has_seg) sprintf("%s AS segment", seg_vals[i]),
               if (has_grade) sprintf("%s AS grade", lit(p$grade)[i]),
               sprintf("%s AS pd", .sql_num(p$pd[i])), sprintf("%s AS lgd", .sql_num(p$lgd[i])),
               sprintf("%s AS r", .sql_num(p$r[i])), sprintf("%s AS k", .sql_num(p$k[i])),
               sprintf("%s AS rw", .sql_num(p$rw[i])))
    sprintf("  SELECT %s%s", paste(parts, collapse = ", "), from_dual)
  }, character(1))
  join <- if (has_seg || has_grade) {
    on <- c(if (has_seg) sprintf("e.%s = p.segment", cols$segment), if (has_grade) sprintf("e.%s = p.grade", cols$grade))
    sprintf("  JOIN pool_params p ON %s", paste(on, collapse = " AND "))
  } else "  CROSS JOIN pool_params p"
  ead <- sprintf("e.%s", cols$ead)
  lgd_col <- sprintf("e.%s", cols$lgd)
  elbe_expr <- if (is.null(cols$elbe)) lgd_col else sprintf("COALESCE(e.%s, %s)", cols$elbe, lgd_col)
  k_def <- if (identical(x$approach, "firb")) "0" else
    sprintf("(CASE WHEN %s - %s > 0 THEN %s - %s ELSE 0 END)", lgd_col, elbe_expr, lgd_col, elbe_expr)
  if (is.null(cols$defaulted)) {
    el_expr <- sprintf("p.pd * p.lgd * %s", ead)
    k_expr <- "p.k"
  } else {
    dcol <- sprintf("e.%s", cols$defaulted)
    el_expr <- sprintf("CASE WHEN %s = 1 THEN %s * %s ELSE p.pd * p.lgd * %s END", dcol, elbe_expr, ead, ead)
    k_expr <- sprintf("CASE WHEN %s = 1 THEN %s ELSE p.k END", dcol, k_def)
  }
  keep <- c(cols$id, cols$segment, cols$grade)
  keep <- keep[!duplicated(keep)]
  sel <- c(sprintf("    e.%s", keep), sprintf("    %s AS ead", ead), sprintf("    %s AS el", el_expr), sprintf("    %s AS k", k_expr))
  final <- if (identical(level, "exposure")) {
    c("SELECT", paste0("    ", c(keep, "ead", "el", "k"), ","), "    12.5 * k AS rw,", "    12.5 * k * ead AS rwa",
      "FROM exposure_capital;")
  } else {
    grp <- if (has_seg) cols$segment else NULL
    c("SELECT", if (has_seg) sprintf("    %s,", grp), "    COUNT(*) AS n,", "    SUM(ead) AS ead,", "    SUM(el) AS el,",
      "    SUM(12.5 * k * ead) AS rwa,",
      "    CASE WHEN SUM(ead) > 0 THEN SUM(12.5 * k * ead) / SUM(ead) ELSE NULL END AS density",
      "FROM exposure_capital", if (has_seg) sprintf("GROUP BY %s", grp), ";")
  }
  approx <- !all(p$homogeneous)
  c("-- =============================================================",
    sprintf("-- scorecraft | IRB capital | %s | %s | %d pools | dialect: %s", x$framework, x$approach, nrow(p), dialect),
    sprintf("-- Generated on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "-- CTE pool_params: PD, LGD, correlation, K and RW per pool, computed in R",
    "--   (floors applied; no normal quantile needed at run time).",
    "-- CTE exposure_capital: el = pd * lgd * ead, rwa = 12.5 * k * ead per exposure;",
    if (!is.null(cols$defaulted)) sprintf("--   rows with %s = 1 use ELBE * ead and K = max(0, LGD - ELBE).", cols$defaulted),
    if (approx) "-- NOTE: at least one pool is not homogeneous: its constants are EAD-weighted, so",
    if (approx) "--   EL and RWA are exact per pool and approximate per exposure. Pass `grade` for finer pools.",
    "-- =============================================================", "",
    "WITH pool_params AS (", paste(rows, collapse = "\n  UNION ALL\n"), "),",
    "exposure_capital AS (", "  SELECT", paste(sel, collapse = ",\n"), sprintf("  FROM %s e", cfg$sql_table), join, ")", "",
    final) |> .sql_lines() |> .sql_out(file)
}

# -- Export ------------------------------------------------------------------ #

#' @rdname scr_export
#' @export
scr_export.scr_capital <- function(x, dir, stamp = TRUE, ...) {
  .scr_local_verbose(x)
  .need_openxlsx()
  out_dir <- .export_dir(dir, stamp)
  tag <- .file_tag(x$framework)
  t <- x$totals
  p <- x$params
  cfg_rows <- c(list(framework = p$framework, source = p$source, approach = x$approach, params_modified = isTRUE(p$modified),
                     confidence = p$confidence, scaling_factor = p$scaling_factor, m_default = p$m_default,
                     m_range = paste(p$m_range, collapse = "-"), output_floor = p$output_floor,
                     target_ratio = t$target_ratio, fi_multiplier = p$correlation$fi_multiplier,
                     sme_adjustment = sprintf("%g (%s %g-%g)", p$correlation$sme$adj, p$correlation$sme$unit, p$correlation$sme$lo, p$correlation$sme$hi)),
                stats::setNames(as.list(p$pd_floor$floor), paste0("pd_floor_", p$pd_floor$asset_class)),
                stats::setNames(as.list(p$lgd_floor$unsecured), paste0("lgd_floor_unsecured_", p$lgd_floor$asset_class)))
  bridge <- data.frame(step = c("rwa_irb_no_floors", "input_floors", "rwa_irb", "rwa_sa", "output_floor_rwa", "rwa_reported"),
                       value = c(t$rwa_irb_no_floors, t$rwa_irb - t$rwa_irb_no_floors, t$rwa_irb, t$rwa_sa, t$rwa_floor, t$rwa_reported),
                       stringsAsFactors = FALSE)
  elp <- data.frame(item = c("el", "provisions", "shortfall", "excess", "tier2_cap", "tier2_addback"),
                    value = c(t$el, t$provisions, t$shortfall, t$excess, t$tier2_cap, t$tier2_addback), stringsAsFactors = FALSE)
  sheets <- list(
    "Capital_Summary"         = .kv_table(t),
    "Capital_Config"          = .kv_table(cfg_rows),
    "Segments_Reconciliation" = x$segments,
    "Pools"                   = x$pools,
    "Floors_Impact"           = x$floors,
    "Output_Floor_Bridge"     = bridge,
    "Sensitivity"             = x$sensitivity,
    "Concentration"           = x$concentration,
    "EL_vs_Provisions"        = elp,
    "Exposures"               = if (is.null(x$exposures)) NULL else if (nrow(x$exposures) > 1e6) utils::head(x$exposures, 1e6) else x$exposures,
    "Model_Card"              = .kv_table(x$model_card),
    "Decision_Ledger"         = x$ledger)
  files <- list(xlsx = .scr_write_xlsx(sheets, file.path(out_dir, sprintf("capital_%s.xlsx", tag))),
                sql = file.path(out_dir, sprintf("sql_capital_%s.sql", tag)))
  writeLines(scr_sql(x), files$sql)
  for (f in files) msg("  %s", f)
  x$files <- files
  invisible(x)
}

# NSE column names used in data.table expressions of this file
utils::globalVariables(c(
  "rw"
))
