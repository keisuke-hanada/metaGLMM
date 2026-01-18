#' meta-analysis using generalized linear mixed effects model with aggregate data
#'
#'@param formula formula of meta-analysis.
#'@param data data of meta-analysis.
#'@param vi variances in each study.
#'@param ni subjects in each study.
#'@param tau2 a value of between-study variance.
#'@param family a family objects for models.
#'@param tau2_var If TRUE (default), tau2 is estimated. If FALSE, assume tau2 is true value.
#'@param rstdnorm random numbers. The default is 5000 quasi-random values by Sobol' quasi-random sequences
#'@param skip.hessian Bypass Hessian calculation?
#'@param start_beta start values of beta.
#'@param start_tau2 start value of tau2.
#'
#'@return object of meta-analysis under generalized linear mixed effects model
#'
#'@export
metaGLMM <- function(formula, data, vi, ni, tau2, family, tau2_var=TRUE, re_group=NULL, trt=NULL,
                 rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000),
                 skip.hessian = FALSE, fast=FALSE,
                 start_beta=NULL, start_tau2=NULL, ...){

  ## set initial value
  trm <- terms(formula, data = data)
  vars <- attr(trm, "term.labels")
  has_int <- attr(trm, "intercept") == 1
  if (has_int) vars <- c("(Intercept)", vars)

  ## initial beta
  if (is.null(start_beta)) {
    vals <- rnorm(length(vars)) - 2
    names(vals) <- vars
  } else {
    # accept either named numeric vector or list
    if (is.list(start_beta)) start_beta <- unlist(start_beta, use.names = TRUE)

    if (is.null(names(start_beta))) {
      stop("start_beta must be a *named* vector (or a named list).")
    }
    if (!all(vars %in% names(start_beta))) {
      miss <- setdiff(vars, names(start_beta))
      stop(paste0("start_beta is missing: ", paste(miss, collapse = ", ")))
    }
    vals <- as.numeric(start_beta[vars])
    names(vals) <- vars
  }
  beta_init <- as.list(vals)
  names(beta_init) <- vars


  ## define log likelihood
  if (fast) {
    ll <- make_ll_fun.fast(formula=formula, data=data, vi=vi, ni=ni,
                      tau2=tau2, family=family, tau2_var=tau2_var,
                      rstdnorm=rstdnorm, re_group=re_group, trt=trt)

  } else {
    ll <- make_ll_fun(formula=formula, data=data, vi=vi, ni=ni,
                      tau2=tau2, family=family, tau2_var=tau2_var,
                      rstdnorm=rstdnorm, re_group=re_group, trt=trt)
  }


  ## maximum likelihood estimation
  if (tau2_var) {
    lower <- as.list(numeric(length(beta_init))-Inf)
    upper <- as.list(numeric(length(beta_init))+Inf)
    names(lower) <- names(upper) <- vars
    lower$tau2 <- 0
    upper$tau2 <- Inf

    if (is.null(start_tau2)) start_tau2 <- 1
    start_list <- c(beta_init, tau2 = start_tau2)

    fit <- bbmle::mle2(ll, start=start_list,
                       skip.hessian=skip.hessian, method= "L-BFGS-B",
                lower = lower, upper = upper,
                control = list(maxit=2e4), ...)

  }else {
    fit <- bbmle::mle2(ll, start=beta_init,
                       skip.hessian=skip.hessian,
                control = list(reltol=1e-12, maxit=2e4), ...)
  }


  return(fit)
}
#' confidence interval for metaGLMM object.
#'
#'@param object object of metaGLMM function.
#'@param parm a vector of parameter names for estimating confidence interval.
#'@param level a value of significance level. The default is 0.95.
#'@param method method of estimating confidence interval. The default is 'PLSBC'.
#'@param renge.c a value of search renge. The default is 30.
#'
#'@return a matrix of confidence interval using profile likelihood with simple Bartlett correction.
#'
#'@export
ci_metaGLMM <- function(object, parm=names(coef(object))[-length(coef(object))],
                        method="PLSBC", level=0.95, renge.c=30){

  if (method=="PLSBC") {
    res <- confint_SBC(object, parm=parm, level=level, renge.c=renge.c)
  } else if (method=="PL") {
    res <- confint_PL(object, parm=parm, level=level, renge.c=renge.c)
  } else if (method=="AN") {
    res <- confint_AN(object, parm=parm, level=level)
  } else {
    message("method should be 'AN', PL', or 'PLSBC'.")
  }

  return(res)

}
#' confidence interval for metaGLMM object using normal approximation.
#'
#'@param object object of metaGLMM function.
#'@param parm a vector of parameter names for estimating confidence interval.
#'@param level a value of significance level. The default is 0.95.
#'
#'@return a matrix of confidence interval using normal approximation.
#'
#'@export
confint_AN <- function(object, parm=names(coef(object))[-length(coef(object))], level=0.95){

  lower <- coef(object)[parm] - qnorm(1-(1-level)/2) * sqrt(diag(vcov(object))[parm])
  upper <- coef(object)[parm] + qnorm(1-(1-level)/2) * sqrt(diag(vcov(object))[parm])

  val <- cbind(lower, upper)

  return(val)

}

#' confidence interval for metaGLMM object using profile likelihood method.
#'
#'@param object object of metaGLMM function.
#'@param parm a vector of parameter names for estimating confidence interval.
#'@param level a value of significance level. The default is 0.95.
#'@param renge.c a value of search renge. The default is 30.
#'
#'@return a matrix of confidence interval using profile likelihood method.
#'
#'@export
confint_PL <- function(object, parm=names(coef(object))[-length(coef(object))], level=0.95, renge.c=30){

  vars <- names(coef(object))
  lower <- upper <- numeric(length(parm))
  names(lower) <- names(upper) <- parm

  for(i in 1:length(parm)){

    var_name <- parm[i]
    init <- numeric(length(vars)-1)

    mll <- function(var, init){

      names(init) <- vars[vars!=var_name]
      args_list <- as.list(init)
      value <- c(var, args_list)
      names(value)[1] <- var_name
      do.call(object@minuslogl, value)

    }


    pl <- function(var, init){

      if (length(init)==1) {

        renge_tau2 <- renge.c*diag(vcov(object))["tau2"]*qnorm(1-(1-level)/2)
        x <- optimize(mll, var=var, interval=c(max(0,-renge_tau2+coef(object)[2]), renge_tau2+coef(object)[2]))
        tau2h <- x$minimum
        return_val <- 2*(x$objective + logLik(object)) - qchisq(0.95, df=1)

      } else {

        x <- optim(par=init, mll, var=var)
        tau2h <- x$par[length(vars)-1]
        return_val <- 2*(x$value + logLik(object)) - qchisq(0.95, df=1)

      }

      return(return_val)

    }

    renge <- renge.c*diag(vcov(object))[var_name]*qnorm(1-(1-level)/2)
    cil <- try(uniroot(pl, init=init, interval=c(-renge+coef(object)[var_name], coef(object)[var_name]))$root)
    if (class(cil)[1]=="try-error") {
      lower[i] <- NA
    } else {
      lower[i] <- cil
    }

    ciu <- try(uniroot(pl, init=init, interval=c(coef(object)[var_name], renge+coef(object)[var_name]))$root)
    if (class(ciu)[1]=="try-error") {
      upper[i] <- NA
    } else {
      upper[i] <- ciu
    }

  }

  val <- cbind(lower, upper)

  return(val)
}
#' confidence interval for metaGLMM object using profile likelihood with simple Bartlett correction.
#'
#'@param object object of metaGLMM function.
#'@param parm a vector of parameter names for estimating confidence interval.
#'@param level a value of significance level. The default is 0.95.
#'@param renge.c a value of search renge. The default is 30.
#'
#'@return a matrix of confidence interval using profile likelihood with simple Bartlett correction.
#'
#'@export
confint_SBC <- function(object, parm=names(coef(object))[-length(coef(object))], level=0.95, renge.c=30){

  vars <- names(coef(object))
  lower <- upper <- numeric(length(parm))
  names(lower) <- names(upper) <- parm
  vi <- environment(object)$vi

  for(i in 1:length(parm)){

    var_name <- parm[i]
    init <- numeric(length(vars)-1)

    mll <- function(var, init){

      names(init) <- vars[vars!=var_name]
      args_list <- as.list(init)
      value <- c(var, args_list)
      names(value)[1] <- var_name
      do.call(object@minuslogl, value)

    }


    plsbc <- function(var, init){

      if (length(init)==1) {

        renge_tau2 <- renge.c*diag(vcov(object))["tau2"]*qnorm(1-(1-level)/2)
        x <- optimize(mll, var=var, interval=c(max(0,-renge_tau2+coef(object)[2]), renge_tau2+coef(object)[2]))
        tau2h <- x$minimum
        correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
        return_val <- 2*(x$objective + logLik(object))/(1+2*correct) - qchisq(0.95, df=1)

      } else {

        x <- optim(par=init, mll, var=var)
        tau2h <- x$par[length(vars)-1]
        correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
        return_val <- 2*(x$value + logLik(object))/(1+2*correct) - qchisq(0.95, df=1)

      }

      return(return_val)

    }


    renge <- renge.c*diag(vcov(object))[var_name]*qnorm(1-(1-level)/2)
    cil <- try(uniroot(plsbc, init=init, interval=c(-renge+coef(object)[var_name], coef(object)[var_name]))$root)
    if (class(cil)[1]=="try-error") {
      lower[i] <- NA
    } else {
      lower[i] <- cil
    }

    ciu <- try(uniroot(plsbc, init=init, interval=c(coef(object)[var_name], renge+coef(object)[var_name]))$root)
    if (class(ciu)[1]=="try-error") {
      upper[i] <- NA
    } else {
      upper[i] <- ciu
    }

  }

  val <- cbind(lower, upper)

  return(val)
}
#' confidence interval for metaGLMM object using profile likelihood with Godambe calibrated simple Bartlett correction.
#'
#'@param object object of metaGLMM function.
#'@param parm a vector of parameter names for estimating confidence interval.
#'@param level a value of significance level. The default is 0.95.
#'@param renge.c a value of search renge. The default is 30.
#'
#'@return a matrix of confidence interval using profile likelihood with simple Bartlett correction.
#'
#'@export
confint_GSBC <- function(object, parm=names(coef(object))[-length(coef(object))], level=0.95, renge.c=30){

  vars <- names(coef(object))
  lower <- upper <- numeric(length(parm))
  names(lower) <- names(upper) <- parm
  vi <- environment(object)$vi

  for(i in 1:length(parm)){

    var_name <- parm[i]
    # init <- numeric(length(vars)-1)
    init <- coef(object)[names(coef(object))!=var_name]

    mll <- function(var, init){

      names(init) <- vars[vars!=var_name]
      args_list <- as.list(init)
      value <- c(var, args_list)
      names(value)[1] <- var_name
      do.call(object@minuslogl, value)

    }


    plsbc <- function(var, init){

      if (length(init)==1) {

        renge_tau2 <- renge.c*diag(vcov(object))["tau2"]*qnorm(1-(1-level)/2)
        x <- optimize(mll, var=var, interval=c(max(0,-renge_tau2+coef(object)[2]), renge_tau2+coef(object)[2]))
        tau2h <- x$minimum
        correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
        r_tau2 <- calc_r_tau2(var=var, init=x$par, tau2h=tau2h, mll=mll,
                              tau2.min=1e-4, diff_step=1e-2, cap=50)
        if (!is.finite(r_tau2)) r_tau2 <- 1
        correct <- correct * r_tau2
        return_val <- 2*(x$objective + logLik(object))/(1+2*correct) - qchisq(0.95, df=1)

      } else {

        x <- optim(par=init, mll, var=var)
        tau2h <- x$par[length(vars)-1]
        correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
        r_tau2 <- calc_r_tau2(var=var, init=x$par, tau2h=tau2h, mll=mll,
                              tau2.min=1e-4, diff_step=1e-2, cap=10)
        if (!is.finite(r_tau2)) r_tau2 <- 1
        correct <- correct * r_tau2
        return_val <- 2*(x$value + logLik(object))/(1+2*correct) - qchisq(0.95, df=1)

      }

      return(return_val)

    }


    renge <- renge.c*diag(vcov(object))[var_name]*qnorm(1-(1-level)/2)
    cil <- try(uniroot(plsbc, init=init, interval=c(-renge+coef(object)[var_name], coef(object)[var_name]))$root)
    if (class(cil)[1]=="try-error") {
      lower[i] <- NA
    } else {
      lower[i] <- cil
    }

    ciu <- try(uniroot(plsbc, init=init, interval=c(coef(object)[var_name], renge+coef(object)[var_name]))$root)
    if (class(ciu)[1]=="try-error") {
      upper[i] <- NA
    } else {
      upper[i] <- ciu
    }

  }

  val <- cbind(lower, upper)

  return(val)
}
#' profile likelihood function for metaGLMM object.
#'
#'@param value a vector of parameter values.
#'@param object object of metaGLMM function.
#'@param parm parameter name for profile likelihood function.
#'
#'@return a values of profile likelihood function.
#'
#'@export
profile_ll <- function(value, object, parm){

  vars <- names(coef(object))
  var_name <- parm
  init <- numeric(length(vars)-1)

  mll <- function(var, init){

    names(init) <- vars[vars!=var_name]
    args_list <- as.list(init)
    value <- c(var, args_list)
    names(value)[1] <- var_name
    do.call(object@minuslogl, value)

  }


  pl <- function(var, init){

    if (length(init)==1) {

      renge_tau2 <- renge.c*diag(vcov(object))["tau2"]*qnorm(1-(1-level)/2)
      x <- optimize(mll, var=var, interval=c(max(0,-renge_tau2+coef(object)[2]), renge_tau2+coef(object)[2]))
      tau2h <- x$minimum
      return_val <- 2*(x$objective + logLik(object))

    } else {

      x <- optim(par=init, mll, var=var)
      tau2h <- x$par[length(vars)-1]
      return_val <- 2*(x$value + logLik(object))

    }

    return(return_val)

  }


  plvals <- sapply(value, pl, init=init)
  return(plvals)
}
#' profile likelihood with simple Bartlett correction for metaGLMM object.
#'
#'@param value a vector of parameter values.
#'@param object object of metaGLMM function.
#'@param parm parameter name for profile likelihood function.
#'
#'@return a values of profile likelihood function.
#'
#'@export
profile_ll_sbc <- function(value, object, parm){

  vars <- names(coef(object))
  var_name <- parm
  init <- numeric(length(vars)-1)
  vi <- environment(object)$vi

  mll <- function(var, init){

    names(init) <- vars[vars!=var_name]
    args_list <- as.list(init)
    value <- c(var, args_list)
    names(value)[1] <- var_name
    do.call(object@minuslogl, value)

  }


  plsbc <- function(var, init){

    if (length(init)==1) {

      renge_tau2 <- renge.c*diag(vcov(object))["tau2"]*qnorm(1-(1-level)/2)
      x <- optimize(mll, var=var, interval=c(max(0,-renge_tau2+coef(object)[2]), renge_tau2+coef(object)[2]))
      tau2h <- x$minimum
      correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
      return_val <- 2*(x$objective + logLik(object))/(1+2*correct)

    } else {

      x <- optim(par=init, mll, var=var)
      tau2h <- x$par[length(vars)-1]
      correct <- sum(1/(vi+tau2h)^3) / (sum(1/(vi+tau2h)) * sum(1/(vi+tau2h)^2))
      return_val <- 2*(x$value + logLik(object))/(1+2*correct)

    }

    return(return_val)

  }

  plvals <- sapply(value, plsbc, init=init)
  return(plvals)
}


#' Godambe calibration
calc_r_tau2 <- function(var, init, tau2h, mll, tau2.min = 1e-6,
                        diff_step = 1e-3, cap = 50, smooth = 1) {

  # init must be a named numeric vector, with possible element "tau2"
  if (is.null(names(init)) || !("tau2" %in% names(init))) {
    return(1)
  }

  t0 <- max(tau2h, tau2.min)
  h  <- diff_step * (1 + t0)

  tp <- t0 + h
  tm <- t0 - h
  if (tm <= tau2.min) {
    # one-sided difference to avoid boundary stickiness
    tm <- t0
    tp <- t0 + 2*h
  }
  tp <- max(tp, tau2.min)
  tm <- max(tm, tau2.min)

  eval_mll_tau2 <- function(tau2_val) {
    init2 <- init
    init2["tau2"] <- tau2_val
    val <- mll(var = var, init = init2)
    llv <- attr(val, "ll_vec")
    list(mll = val, llvec = llv)
  }

  e0 <- eval_mll_tau2(t0)
  ep <- eval_mll_tau2(tp)
  em <- eval_mll_tau2(tm)

  # basic finiteness checks
  if (!is.finite(e0$mll) || !is.finite(ep$mll) || !is.finite(em$mll)) return(1)
  if (is.null(ep$llvec) || is.null(em$llvec)) return(1)
  if (any(!is.finite(ep$llvec)) || any(!is.finite(em$llvec))) return(1)

  denom <- (tp - tm)
  if (denom <= 0) return(1)

  # H: curvature of minuslogl
  # If one-sided (tm==t0), this becomes a forward-difference curvature proxy.
  if (tm == t0) {
    # forward second difference: f(t0+2h) - 2f(t0+h) + f(t0) over h^2
    # Here ep is at t0+2h and em is at t0 (because tm=t0). We still need f(t0+h):
    eh <- eval_mll_tau2(t0 + h)
    if (!is.finite(eh$mll) || is.null(eh$llvec) || any(!is.finite(eh$llvec))) return(1)
    H <- (ep$mll - 2*eh$mll + e0$mll) / (h^2)
    # J uses scores from llvec: (f(t0+2h)-f(t0))/2h is too rough; use (f(t0+h)-f(t0))/h
    s_g <- (eh$llvec - em$llvec) / h
    J <- sum(s_g^2)
  } else {
    H <- (ep$mll - 2*e0$mll + em$mll) / ((denom/2)^2)
    s_g <- (ep$llvec - em$llvec) / denom
    J <- sum(s_g^2)
  }

  if (!is.finite(H) || H <= 0) return(1)
  if (!is.finite(J) || J <= 0) return(1)

  r_raw <- H / J
  if (!is.finite(r_raw)) return(1)

  # optional smoothing toward 1 (smooth=1 means no smoothing)
  r <- 1 + smooth * (r_raw - 1)

  # enforce conservativeness and cap
  r <- max(1, min(cap, r))
  cat("tau2h=", t0, " H=", H, " J=", J, " r=", r, "\n")

  r
}


confint_PL_count <- function(object, parm=names(coef(object))[-length(coef(object))],
                             level=0.95, renge.c=30){

  vars <- names(coef(object))
  lower <- upper <- numeric(length(parm))
  names(lower) <- names(upper) <- parm

  crit <- qchisq(level, df=1)

  for(i in seq_along(parm)){

    var_name <- parm[i]
    init <- numeric(length(vars)-1)

    ctr_pl  <- 0L
    ctr_mll <- 0L

    mll <- function(var, init){
      ctr_mll <<- ctr_mll + 1L
      names(init) <- vars[vars!=var_name]
      args_list <- as.list(init)
      value <- c(var, args_list)
      names(value)[1] <- var_name
      do.call(object@minuslogl, value)
    }

    pl <- function(var, init){
      ctr_pl <<- ctr_pl + 1L

      if (length(init)==1) {
        se_tau2 <- sqrt(diag(vcov(object))["tau2"])
        tau2hat <- unname(coef(object)["tau2"])
        renge_tau2 <- renge.c*se_tau2*qnorm(1-(1-level)/2)
        x <- optimize(mll, var=var, interval=c(max(0, tau2hat-renge_tau2), tau2hat+renge_tau2), init=init)
        return( 2*(x$objective + logLik(object)) - crit )
      } else {
        x <- optim(par=init, fn=function(p) mll(var, p), method="Nelder-Mead")
        return( 2*(x$value + logLik(object)) - crit )
      }
    }

    renge <- renge.c*sqrt(diag(vcov(object))[var_name])*qnorm(1-(1-level)/2)
    mid <- unname(coef(object)[var_name])

    cil <- try(uniroot(pl, init=init, interval=c(mid-renge, mid))$root, silent=TRUE)
    lower[i] <- if (inherits(cil,"try-error")) NA_real_ else cil

    ciu <- try(uniroot(pl, init=init, interval=c(mid, mid+renge))$root, silent=TRUE)
    upper[i] <- if (inherits(ciu,"try-error")) NA_real_ else ciu

    message(sprintf("[%s] pl calls=%d, mll calls=%d", var_name, ctr_pl, ctr_mll))
  }

  cbind(lower, upper)
}



