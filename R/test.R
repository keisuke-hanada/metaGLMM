# library(metaGLMM)
# library(dplyr)
#
# source(file.path("R", "likelihood.R"))
# source(file.path("R", "metaGLMM.R"))
#
# set.seed(1234)
#
# n <- 10
# nk <- numeric(n) + 30
# zk <- rnorm(n)
# beta <- c(1,0)
# tau2 <- 3
# vk <- (rchisq(n, df=1) - 1) * sqrt(tau2/2)
# # vk <- rnorm(n, sd=sqrt(tau2))
# v2yk <- 50/nk
# yk <- beta[1] + beta[2]*zk + rnorm(n, sd=sqrt(v2yk)) + vk
#
# dat3 <- model.frame(formula,
#                     data.frame(yk, zk)
#                     )
#
#
#
# ma.grma <- metaGLMM(formula=yk ~ 1 + zk, data=dat3, vi=v2yk, ni=nk, tau2=NA,
#                     family=gaussian(link="identity"), tau2_var=TRUE)
# summary(ma.grma)
#
#
#
# system.time(
# ci.pl <- confint_PL(ma.grma, parm="zk")
# )
# ci.pl
#
# system.time(
# ci.plsbc <- confint_SBC(ma.grma, parm="zk")
# )
# ci.plsbc
#
# system.time(
# ci.plgsbc <- confint_GSBC(ma.grma, parm="zk")
# )
# ci.plgsbc
#
#
#
