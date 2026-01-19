library(metaGLMM)
library(dplyr)

# source(file.path("R", "likelihood.R"))
# source(file.path("R", "metaGLMM.R"))
# Rcpp::sourceCpp("src/mc_integral.cpp")
# Rcpp::sourceCpp("src/ghq_integral.cpp")

set.seed(1234)

n <- 30
nk <- numeric(n) + 30
zk <- rnorm(n)
beta <- c(1,0)
tau2 <- 3
vk <- (rchisq(n, df=1) - 1) * sqrt(tau2/2)
# vk <- rnorm(n, sd=sqrt(tau2))
v2yk <- 50/nk
yk <- beta[1] + beta[2]*zk + rnorm(n, sd=sqrt(v2yk)) + vk

dat3 <- model.frame(formula=yk ~ 1 + zk,
                    data.frame(yk, zk)
                    )



ma.grma <- metaGLMM(formula=yk ~ 1 + zk, data=dat3, vi=v2yk, ni=nk, tau2=NA,
                    family=gaussian(link="identity"), tau2_var=TRUE, fast=FALSE)
summary(ma.grma)

system.time(
  ci.pl <- confint_PL(ma.grma, parm="zk")
)
ci.pl

system.time(
  ci.plsbc <- confint_SBC(ma.grma, parm="zk")
)
ci.plsbc

system.time(
  ci.plgsbc <- confint_GSBC(ma.grma, parm="zk")
)
ci.plgsbc


ma.grma <- metaGLMM(formula=yk ~ 1 + zk, data=dat3, vi=v2yk, ni=nk, tau2=NA,
                    family=gaussian(link="identity"), tau2_var=TRUE, fast=TRUE)
summary(ma.grma)



system.time(
ci.pl <- confint_PL(ma.grma, parm="zk")
)
ci.pl

system.time(
ci.plsbc <- confint_SBC(ma.grma, parm="zk")
)
ci.plsbc

system.time(
ci.plgsbc <- confint_GSBC(ma.grma, parm="zk")
)
ci.plgsbc





set.seed(1234)

n <- 10
nk <- numeric(n) + 30
zk <- rnorm(n)
beta <- c(-2,0.5)
tau2 <- 0
# vk <- (rchisq(n, df=1) - 1) * sqrt(tau2/2)
vk <- rnorm(n, sd=sqrt(tau2))

eta <- beta[1] + beta[2]*zk + vk
pk <- exp(eta) / (1+exp(eta))
yk <- rbinom(n, size=nk, prob=pk) / nk
yk[yk==0] <- 1e-10 / nk[yk==0]
thetahk <- log((nk*yk)/(nk-nk*yk))
v2yk <- 1/(nk*yk) + 1/(nk-nk*yk)
dat3 <- model.frame(formula=yk ~ 1 + zk,
                    data.frame(yk, zk)
)




ma.grma <- metaGLMM(formula=yk ~ 1 + zk, data=dat3, vi=v2yk, ni=nk, tau2=NA,
                    family=binomial(link="logit"), tau2_var=TRUE, fast=FALSE)
summary(ma.grma)

system.time(
  ci.pl <- confint_PL(ma.grma, parm="zk")
)
ci.pl




ma.grma <- metaGLMM(formula=yk ~ 1 + zk, data=dat3, vi=v2yk, ni=nk, tau2=NA, tau2_param="log_tau2",
                    family=binomial(link="logit"), tau2_var=TRUE, fast=TRUE,
                    start_beta=coef(ma.grma)[-3], start_tau2 = 1)
summary(ma.grma)
exp(coef(ma.grma))


system.time(
  ci.pl <- confint_PL(ma.grma, parm="zk")
)
ci.pl

system.time(
  ci.pl <- confint_SBC(ma.grma, parm="zk")
)
ci.pl

system.time(
  ci.pl <- confint_GSBC(ma.grma, parm="zk")
)
ci.pl




