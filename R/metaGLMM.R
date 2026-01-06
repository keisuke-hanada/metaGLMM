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
metaGLMM <- function(formula, data, vi, ni, tau2, family, tau2_var=TRUE,
                 rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000),
                 skip.hessian = FALSE,
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
  ll <- make_ll_fun(formula=formula, data=data, vi=vi, ni=ni,
                    tau2=tau2, family=family, tau2_var=tau2_var,
                    rstdnorm=rstdnorm)


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



