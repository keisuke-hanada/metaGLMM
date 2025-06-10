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
                       tau2.min=1e-6,
                       rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000)) {

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
  for (j in seq_len(strata)) {
    etaj   <- etak[j]
    eta_mc <- etaj + rst_tau2

    mu_mc  <- family$linkinv(eta_mc)
    theta_mc <- clink(mu_mc)

    b_mc   <- b(theta_mc)

    yj     <- yk[j]
    nj     <- ni[j]
    aj     <- a_phik[j]
    factor <- nj / aj

    z_mc   <- factor * ( yj * theta_mc - b_mc )

    logsum <- logSumExp_simple(z_mc)
    total <- total + ( - (logsum - log_n_monte) )
  }

  return(total)

}



logSumExp_simple <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}



make_ll_fun <- function(formula, data, vi, ni, tau2, family, tau2_var=FALSE,
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
