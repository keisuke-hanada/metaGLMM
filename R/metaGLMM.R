#' meta-analysis using generalized linear mixed effects model with aggregate data
#'
#'@param formula formula of meta-analysis.
#'@param data data of meta-analysis.
#'@param vi variances in each study.
#'@param ni subjects in each study.
#'@param tau2 a value of between-study variance.
#'@param family a family objects for models.
#'@param tau2_var If FALSE (default), tau2 is estimated. If TRUE, assume tau2 is true value.
#'@param rstdnorm random numbers. The default is 5000 quasi-random values by Sobol' quasi-random sequences
#'
#'@return object of meta-analysis under generalized linear mixed effects model
#'
#'@export
metaGLMM <- function(formula, data, vi, ni, tau2, family, tau2_var=FALSE,
                 rstdnorm=qnorm((qrng::sobol(5000, d=1, scrambling=1)*(5000-1) + 0.5) / 5000)){

  ## set initial value
  trm <- terms(formula, data = data)
  vars <- attr(trm, "term.labels")
  has_int <- attr(trm, "intercept") == 1
  if (has_int) vars <- c("(Intercept)", vars)
  vals <- rnorm(length(vars)) - 2
  beta_init <- as.list(vals)
  names(beta_init) <- vars


  ## define log likelihood
  ll <- make_ll_fun(formula=formula, data=data, vi=vi, ni=ni, tau2=tau2, family=family, tau2_var=tau2_var, rstdnorm=rstdnorm)


  ## maximum likelihood estimation
  if (tau2_var) {
    lower <- as.list(numeric(length(beta_init))-Inf)
    upper <- as.list(numeric(length(beta_init))+Inf)
    names(lower) <- names(upper) <- vars
    lower$tau2 <- 0
    upper$tau2 <- Inf


    fit <- bbmle::mle2(ll, start=c(beta_init, tau2=1), skip.hessian=FALSE, method= "L-BFGS-B",
                lower = lower, upper = upper,
                control = list(maxit=2e4))

  }else {
    fit <- bbmle::mle2(ll, start=beta_init, skip.hessian=FALSE,
                control = list(reltol=1e-12, maxit=2e4))
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
    res <- confint_SBC(object, parm=names(coef(object))[-length(coef(object))], level=0.95, renge.c=30)
  } else if (method=="PL") {
    res <- confint_PL(object, parm=names(coef(object))[-length(coef(object))], level=0.95, renge.c=30)
  }else {
    message("method should be 'PL' or 'PLSBC'.")
  }

  return(res)

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
