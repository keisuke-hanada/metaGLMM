#'likelihood function
#'
#'@param formula formula of meta-analysis.
#'@param data data of meta-analysis.
#'@param vi variances in each study.
#'@param ni subjects in each study.
#'@param beta a vector of parameter.
#'@param tau2 a value of between-study variance.
#'@param family a family objects for models.
#'@param tau2.min minimum value of tau2 > 0.
#'@param rstdnorm random numbers. The default is 5000 quasi-random values by Sobol' quasi-random sequences
#'
#'@return value of likelihood function
#'
#'@export
likelihood <- function(formula, data, vi, ni, beta, tau2, family=binomial(link="logit"),
                       tau2.min=1e-6, re_group=NULL, trt=NULL,
                       rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000)) {

  if (anyNA(re_group)) stop("re_group must not contain NA.")

  use_random_slope <- !is.null(re_group)
  if (use_random_slope) {
    if (is.null(trt)) stop("re_group is specified, so trt must be provided (0/1).")
    if (!(trt %in% names(data))) stop("trt was not found in data/model.frame.")
    Z <- data[[trt]]
    if (is.logical(Z)) Z <- as.numeric(Z)
    if (is.factor(Z)) {
      if (nlevels(Z) != 2L) stop("trt factor must have 2 levels.")
      Z <- as.numeric(Z == levels(Z)[2L])
    }
    if (anyNA(Z) || !all(Z %in% c(0,1))) stop("trt must be coded as 0/1 (or logical / 2-level factor).")
  } else {
    Z <- rep.int(1, nrow(mf))
  }

  ## model
  yk <- model.response(data)
  X <- model.matrix(formula, data)
  strata <- length(yk)
  tau2.len <- length(tau2)
  tau2 <- max(tau2, tau2.min)
  n.monte <- length(rstdnorm)


  ## model spesific function
  fam_name <- family$family
  if (fam_name == "binomial") {

    clink <- function(x){log(x/(1-x))}
    a <- function(x){numeric(strata)+1}
    b <- function(x){log(1+exp(x))}

  }else if (fam_name == "poisson") {

    clink <- function(x){log(x)}
    a <- function(x){numeric(strata)+1}
    b <- function(x){exp(x)}

  }else if (fam_name == "Gamma") {

    clink <- function(x){-1/x}
    a <- function(x){x}
    b <- function(x){-log(-x)}

  }else if (fam_name == "gaussian") {

    clink <- function(x){x}
    a <- function(x){x}
    b <- function(x){x^2/2}

  }


  etak <- X %*% beta
  a_phik <- a(vi)

  total <- 0.0
  log_n_monte <- log(n.monte)
  rst_tau2 <- rstdnorm * sqrt(tau2)

  if (is.null(re_group)) {
    re_group <- seq_len(strata)
  }
  idx_list <- split(seq_len(strata), re_group)

  for (g in idx_list) {

    ## accumulate arm/row contributions for each MC draw
    z_mc_g <- numeric(n.monte)

    for (j in g) {

      if (use_random_slope) {
        eta_mc <- etak[j] + rst_tau2 * Z[j]
      } else {
        eta_mc <- etak[j] + rst_tau2
      }

      mu_mc    <- family$linkinv(eta_mc)
      theta_mc <- clink(mu_mc)
      b_mc     <- b(theta_mc)

      yj     <- yk[j]
      nj     <- ni[j]
      aj     <- a_phik[j]
      factor <- nj / aj

      z_mc_g <- z_mc_g + factor * ( yj * theta_mc - b_mc )
    }

    logsum <- logSumExp_simple(z_mc_g)
    total <- total + ( - (logsum - log_n_monte) )
  }


  return(total)

}



logSumExp_simple <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}



make_ll_fun <- function(formula, data, vi, ni, tau2, family, tau2_var=FALSE,
                        tau2.min = 1e-6, re_group=NULL, trt=NULL,
                        rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000),
                        ...) {

  if (anyNA(re_group)) stop("re_group must not contain NA.")

  ## ---- cache model objects ONCE ----
  # data is assumed to be a model.frame OR a data.frame that can be used in model.frame
  if (!inherits(data, "data.frame")) {
    stop("data must be a data.frame or model.frame.")
  }

  # If data is not a model.frame, build it once here
  # (If data is already a model.frame, this is cheap and consistent.)
  mf <- model.frame(formula, data)


  use_random_slope <- !is.null(re_group)
  if (use_random_slope) {
    if (is.null(trt)) stop("re_group is specified, so trt must be provided (0/1).")
    if (!(trt %in% names(mf))) stop("trt was not found in data/model.frame.")
    Z <- mf[[trt]]
    if (is.logical(Z)) Z <- as.numeric(Z)
    if (is.factor(Z)) {
      if (nlevels(Z) != 2L) stop("trt factor must have 2 levels.")
      Z <- as.numeric(Z == levels(Z)[2L])
    }
    if (anyNA(Z) || !all(Z %in% c(0,1))) stop("trt must be coded as 0/1 (or logical / 2-level factor).")
  } else {
    Z <- rep.int(1, nrow(mf))
  }

  yk <- model.response(mf)
  X  <- model.matrix(formula, mf)

  strata <- length(yk)
  n.monte <- length(rstdnorm)
  log_n_monte <- log(n.monte)

  fam_name <- family$family

  # Model-specific a(), clink(), b() defined ONCE
  if (fam_name == "binomial") {
    clink <- function(mu) log(mu/(1-mu))              # logit
    a_fun <- function(v) rep.int(1, strata)
    b_fun <- function(theta) log1p(exp(theta))        # stable-ish
  } else if (fam_name == "poisson") {
    clink <- function(mu) log(mu)
    a_fun <- function(v) rep.int(1, strata)
    b_fun <- function(theta) exp(theta)
  } else if (fam_name == "Gamma") {
    clink <- function(mu) -1/mu
    a_fun <- function(v) v
    b_fun <- function(theta) -log(-theta)
  } else if (fam_name == "gaussian") {
    clink <- function(mu) mu
    a_fun <- function(v) v
    b_fun <- function(theta) theta^2/2
  } else {
    stop("Unsupported family: ", fam_name)
  }

  a_phik <- a_fun(vi)              # length strata
  factor <- ni / a_phik            # length strata


  ## --- group indices fixed here (default: each row is its own group) ---
  if (is.null(re_group)) {
    re_group_use <- seq_len(strata)
  } else {
    re_group_use <- re_group
  }
  idx_list <- split(seq_len(strata), re_group_use)


  # local fast log-sum-exp
  logSumExp_simple_local <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }

  ## ---- create ll with explicit formals for mle2 ----
  trm <- terms(formula, data = mf)
  vars <- attr(trm, "term.labels")
  has_int <- attr(trm, "intercept") == 1
  if (has_int) vars <- c("(Intercept)", vars)
  if (tau2_var) vars <- c(vars, "tau2")
  K <- length(vars)

  arg_list <- vector("list", K)
  names(arg_list) <- vars

  ll <- function() NULL
  formals(ll) <- arg_list

  if (tau2_var) {
    # beta = c((Intercept), ...)
    beta_call <- as.call(c(as.name("c"), lapply(vars[-K], as.name)))
    body(ll) <- bquote({
      beta <- .(beta_call)

      tau2_use <- max(tau2, .(tau2.min))
      rst_tau2 <- .(rstdnorm) * sqrt(tau2_use)

      etak <- as.vector(.(X) %*% beta)

      total <- 0.0
      ll_vec <- numeric(length(.(idx_list)))
      g_id <- 0L
      for (g in .(idx_list)) {
        g_id <- g_id + 1L
        z_mc_g <- numeric(.(n.monte))
        for (j in g) {

          if (.(use_random_slope)) {
            eta_mc <- etak[j] + rst_tau2 * .(Z)[j]
          } else {
            eta_mc <- etak[j] + rst_tau2
          }

          mu_mc  <- .(family)$linkinv(eta_mc)
          theta_mc <- .(clink)(mu_mc)

          z_mc_g <- z_mc_g + .(factor)[j] * ( .(yk)[j] * theta_mc - .(b_fun)(theta_mc) )
        }
        ll_g <- ( - (.(logSumExp_simple_local)(z_mc_g) - .(log_n_monte)) )
        ll_vec[g_id] <- ll_g
        total <- total + ll_g
      }
      attr(total, "ll_vec") <- ll_vec
      total
    })
  } else {
    beta_call <- as.call(c(as.name("c"), lapply(vars, as.name)))
    body(ll) <- bquote({
      beta <- .(beta_call)

      tau2_use <- max(.(tau2), .(tau2.min))
      rst_tau2 <- .(rstdnorm) * sqrt(tau2_use)

      etak <- as.vector(.(X) %*% beta)

      total <- 0.0
      ll_vec <- numeric(length(.(idx_list)))
      g_id <- 0L
      for (g in .(idx_list)) {
        g_id <- g_id + 1L
        z_mc_g <- numeric(.(n.monte))
        for (j in g) {

          if (.(use_random_slope)) {
            eta_mc <- etak[j] + rst_tau2 * .(Z)[j]
          } else {
            eta_mc <- etak[j] + rst_tau2
          }

          mu_mc  <- .(family)$linkinv(eta_mc)
          theta_mc <- .(clink)(mu_mc)
          z_mc_g <- z_mc_g + .(factor)[j] * ( .(yk)[j] * theta_mc - .(b_fun)(theta_mc) )
        }
        ll_g <- ( - (.(logSumExp_simple_local)(z_mc_g) - .(log_n_monte)) )
        ll_vec[g_id] <- ll_g
        total <- total + ll_g
      }
      total
    })
  }

  ll
}



#' Fast version of make_ll_fun: MC/QMC integration and log-mean-exp in Rcpp,
#' while link function / clink / b_fun remain in R (works for non-canonical links).
#'
#' Requirements:
#' - Rcpp must be installed.
#' - This function compiles a small C++ routine once per R session (cached in .GlobalEnv).
#'
#' Fast log-likelihood builder: closed-form for Gaussian(identity),
#' otherwise call a compiled Rcpp core (mc_nll_groups) that supports flexible linkinv.
#'
#' Prerequisite:
#' - A compiled Rcpp function `mc_nll_groups()` must be available in the namespace.
#'   (e.g., placed under src/ and exported via Rcpp attributes, or sourceCpp-ed once.)
#'
make_ll_fun.fast <- function(formula, data, vi, ni, tau2, family, tau2_var = FALSE,
                             tau2.min = 1e-6, re_group = NULL, trt = NULL,
                             ghq_Q = 30L,
                             ...) {

  if (anyNA(re_group)) stop("re_group must not contain NA.")
  if (!inherits(data, "data.frame")) stop("data must be a data.frame or model.frame.")

  mf <- model.frame(formula, data)

  ## --- random slope coding (Z) ---
  use_random_slope <- !is.null(re_group)
  if (use_random_slope) {
    if (is.null(trt)) stop("re_group is specified, so trt must be provided (0/1).")
    if (!(trt %in% names(mf))) stop("trt was not found in data/model.frame.")
    Z0 <- mf[[trt]]
    if (is.logical(Z0)) Z0 <- as.numeric(Z0)
    if (is.factor(Z0)) {
      if (nlevels(Z0) != 2L) stop("trt factor must have 2 levels.")
      Z0 <- as.numeric(Z0 == levels(Z0)[2L])
    }
    if (anyNA(Z0) || !all(Z0 %in% c(0, 1))) {
      stop("trt must be coded as 0/1 (or logical / 2-level factor).")
    }
  } else {
    Z0 <- rep.int(1, nrow(mf))
  }

  yk0 <- model.response(mf)
  X0  <- model.matrix(formula, mf)
  strata <- length(yk0)

  ## --- group id ---
  if (is.null(re_group)) re_group_use <- seq_len(strata) else re_group_use <- re_group

  ## --- family / link ---
  fam_name  <- family$family
  link_name <- family$link

  ## closed-form only for Gaussian(identity)
  use_closed_form_gaussian <- identical(fam_name, "gaussian") && identical(link_name, "identity")

  ## ---- a_fun / family_id (canonical part fixed as you said) ----
  if (fam_name == "gaussian") {
    a_fun <- function(v) v
    family_id <- 4L
  } else if (fam_name == "binomial") {
    a_fun <- function(v) rep.int(1, strata)
    family_id <- 1L
  } else if (fam_name == "poisson") {
    a_fun <- function(v) rep.int(1, strata)
    family_id <- 2L
  } else if (fam_name == "Gamma") {
    a_fun <- function(v) v
    family_id <- 3L
  } else {
    stop("Unsupported family: ", fam_name)
  }

  a_phik  <- a_fun(vi)

  if (fam_name == "gaussian") {
    factor0 <- 1 / a_phik
  } else {
    factor0 <- ni / a_phik
  }


  ## ----- link_id mapping (make.link list) -----
  link_map <- c(
    "logit"="1","probit"="2","cauchit"="3","cloglog"="4","identity"="5",
    "log"="6","sqrt"="7","1/mu^2"="8","inverse"="9"
  )
  link_id <- if (!is.null(link_map[[link_name]])) as.integer(link_map[[link_name]]) else 0L

  ## ----- reorder rows by group to build gptr -----
  ord <- order(re_group_use, seq_along(re_group_use))  # stable within-group
  re_sorted <- re_group_use[ord]
  yk     <- as.numeric(yk0[ord])
  Z      <- as.numeric(Z0[ord])
  factor <- as.numeric(factor0[ord])
  X      <- X0[ord, , drop = FALSE]

  rle_g <- rle(re_sorted)
  lens  <- rle_g$lengths
  gptr  <- as.integer(c(0L, cumsum(lens)))  # 0-based boundaries

  ## --- GHQ nodes/weights (only if needed) ---
  if (!use_closed_form_gaussian) {
    if (!requireNamespace("fastGHQuad", quietly = TRUE)) {
      stop("fastGHQuad is required for fast=TRUE GHQ. Please install fastGHQuad.")
    }
    rule  <- fastGHQuad::gaussHermiteData(as.integer(ghq_Q))
    ghq_x <- as.numeric(rule$x)
    ghq_w <- as.numeric(rule$w)
  } else {
    ghq_x <- ghq_w <- NULL
  }

  ## --- linkinv callback (only used if link_id==0 in C++) ---
  linkinv_fun <- family$linkinv

  ## ---- create ll with explicit formals for mle2 ----
  trm <- terms(formula, data = mf)
  vars <- attr(trm, "term.labels")
  has_int <- attr(trm, "intercept") == 1
  if (has_int) vars <- c("(Intercept)", vars)
  if (tau2_var) vars <- c(vars, "tau2")
  K <- length(vars)

  arg_list <- vector("list", K)
  names(arg_list) <- vars
  ll <- function() NULL
  formals(ll) <- arg_list

  ## ---- helpers (R-side closed-form uses pi, log) ----
  # (keep as constants in the closure)
  pi_const <- base::pi

  if (tau2_var) {
    beta_call <- as.call(c(as.name("c"), lapply(vars[-K], as.name)))

    body(ll) <- bquote({
      beta <- .(beta_call)
      tau2_use <- max(tau2, .(tau2.min))
      etak <- as.vector(.(X) %*% beta)

      if (.(use_closed_form_gaussian)) {
        total <- 0.0
        G <- length(.(gptr)) - 1L
        ll_vec <- numeric(G)

        for (g in seq_len(G)) {
          a <- .(gptr)[g] + 1L
          b <- .(gptr)[g + 1L]
          idx <- a:b

          fj <- .(factor)[idx]
          yj <- .(yk)[idx]
          et <- etak[idx]
          zj <- .(Z)[idx]

          # C0 + B0*v + A0*v^2 with prior N(0,tau2)
          C0 <- sum(fj * (yj * et - 0.5 * et * et))
          B0 <- sum(fj * zj * (yj - et))
          A0 <- -0.5 * sum(fj * zj * zj)

          A  <- A0 - 0.5 / tau2_use  # A < 0

          logLg <- C0 - (B0 * B0) / (4 * A) +
            0.5 * log(.(pi_const) / (-A)) -
            0.5 * log(2 * .(pi_const) * tau2_use)

          ll_g <- -logLg
          ll_vec[g] <- ll_g
          total <- total + ll_g
        }

        attr(total, "ll_vec") <- ll_vec
        total

      } else {

        if (!exists("ghq_nll_groups_fast", mode = "function")) {
          stop("ghq_nll_groups_fast() was not found. Compile and load src/ghq_integral.cpp first.")
        }

        out <- ghq_nll_groups_fast(
          etak        = etak,
          y           = .(yk),
          factor      = .(factor),
          Z           = .(Z),
          gptr        = .(gptr),
          tau2        = tau2_use,
          family_id   = .(family_id),
          link_id     = .(link_id),
          ghq_x       = .(ghq_x),
          ghq_w       = .(ghq_w),
          linkinv_fun = .(linkinv_fun)
        )

        total <- out$nll_total
        attr(total, "ll_vec") <- out$nll_group
        total
      }
    })

  } else {

    beta_call <- as.call(c(as.name("c"), lapply(vars, as.name)))

    body(ll) <- bquote({
      beta <- .(beta_call)
      tau2_use <- max(.(tau2), .(tau2.min))
      etak <- as.vector(.(X) %*% beta)

      if (.(use_closed_form_gaussian)) {
        total <- 0.0
        G <- length(.(gptr)) - 1L

        for (g in seq_len(G)) {
          a <- .(gptr)[g] + 1L
          b <- .(gptr)[g + 1L]
          idx <- a:b

          fj <- .(factor)[idx]
          yj <- .(yk)[idx]
          et <- etak[idx]
          zj <- .(Z)[idx]

          C0 <- sum(fj * (yj * et - 0.5 * et * et))
          B0 <- sum(fj * zj * (yj - et))
          A0 <- -0.5 * sum(fj * zj * zj)
          A  <- A0 - 0.5 / tau2_use

          logLg <- C0 - (B0 * B0) / (4 * A) +
            0.5 * log(.(pi_const) / (-A)) -
            0.5 * log(2 * .(pi_const) * tau2_use)

          total <- total + (-logLg)
        }
        total

      } else {

        if (!exists("ghq_nll_groups_fast", mode = "function")) {
          stop("ghq_nll_groups_fast() was not found. Compile and load src/ghq_integral.cpp first.")
        }

        out <- ghq_nll_groups_fast(
          etak        = etak,
          y           = .(yk),
          factor      = .(factor),
          Z           = .(Z),
          gptr        = .(gptr),
          tau2        = tau2_use,
          family_id   = .(family_id),
          link_id     = .(link_id),
          ghq_x       = .(ghq_x),
          ghq_w       = .(ghq_w),
          linkinv_fun = .(linkinv_fun)
        )
        out$nll_total
      }
    })
  }

  ll
}


