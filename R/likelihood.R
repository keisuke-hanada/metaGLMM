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
                       tau2.min=1e-6, re_group=NULL,
                       rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000)) {

  if (anyNA(re_group)) stop("re_group must not contain NA.")

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

      if (j==g[1]) {
        eta_mc <- etak[j] + rst_tau2
      } else {
        eta_mc <- etak[j]
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



make_ll_fun_old <- function(formula, data, vi, ni, tau2, family, tau2_var=FALSE,
                        rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000), ...) {
  trm <- terms(formula, data = data)
  vars <- attr(trm, "term.labels")
  has_int <- attr(trm, "intercept") == 1
  if (has_int) vars <- c("(Intercept)", vars)
  if (tau2_var) vars <- c(vars, "tau2")
  K <- length(vars)

  arg_list <- vector("list", K)
  # names(arg_list) <- paste0("b", seq_len(K))
  names(arg_list) <- vars

  ll <- function() NULL
  formals(ll) <- arg_list

  if (tau2_var) {

    beta_call <- as.call(c(as.name("c"),
                           lapply(names(arg_list[-K]), as.name)))
    tau2_call <- as.name("tau2")

    body(ll) <- bquote({
      beta <- .(beta_call)

      likelihood(formula=.(formula), data=.(data), vi=.(vi), ni=.(ni),
                 beta=.(beta_call), tau2=.(tau2_call), family=.(family),
                 rstdnorm=.(rstdnorm))
    })


  }else {

    beta_call <- as.call(c(as.name("c"),
                           lapply(names(arg_list), as.name)))

    body(ll) <- bquote({
      beta <- .(beta_call)

      likelihood(formula=.(formula), data=.(data), vi=.(vi), ni=.(ni),
                 beta=.(beta_call), tau2=.(tau2), family=.(family),
                 rstdnorm=.(rstdnorm))
    })

  }



  return(ll)
}

make_ll_fun <- function(formula, data, vi, ni, tau2, family, tau2_var=FALSE,
                        tau2.min = 1e-6, re_group=NULL,
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
      for (g in .(idx_list)) {
        z_mc_g <- numeric(.(n.monte))
        for (j in g) {
          eta_mc <- etak[j] + rst_tau2
          mu_mc  <- .(family)$linkinv(eta_mc)
          theta_mc <- .(clink)(mu_mc)
          z_mc_g <- z_mc_g + .(factor)[j] * ( .(yk)[j] * theta_mc - .(b_fun)(theta_mc) )
        }
        total <- total + ( - (.(logSumExp_simple_local)(z_mc_g) - .(log_n_monte)) )
      }
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
      for (g in .(idx_list)) {
        z_mc_g <- numeric(.(n.monte))
        for (j in g) {

          if (j==g[1]) {
            eta_mc <- etak[j] + rst_tau2
          } else {
            eta_mc <- etak[j]
          }

          mu_mc  <- .(family)$linkinv(eta_mc)
          theta_mc <- .(clink)(mu_mc)
          z_mc_g <- z_mc_g + .(factor)[j] * ( .(yk)[j] * theta_mc - .(b_fun)(theta_mc) )
        }
        total <- total + ( - (.(logSumExp_simple_local)(z_mc_g) - .(log_n_monte)) )
      }
      total
    })
  }

  ll
}


