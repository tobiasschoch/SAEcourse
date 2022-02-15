# NOTE: Change the following path to YOUR COMPUTER!
setwd("C:/My/code/SAEcourse/application")
library(sae)
#-------------------------------------------------------------------------------
# Preparations
#-------------------------------------------------------------------------------
# Survey data (individual level)
datLCS <- read.table("datLCS.txt", header = TRUE, sep = "\t", dec = ",")
# Auxiliary data (area level)
auxLCS <- read.table("auxLCS.txt", header = TRUE, sep = "\t", dec = ",")

#-------------------------------------------------------------------------------
# TASK 1: Direct estimators of average income
#-------------------------------------------------------------------------------
# Split datLCS into a list by domain (dom)
datLCS_dom <- split(datLCS, datLCS$dom)
# Hajek estimator of the mean
sapply(datLCS_dom, function(u) weighted.mean(u$income, u$w))
# Hajek estimator of mean and (approx.) variance
hajek <- function(x, w)
{
    avg <- weighted.mean(x, w)                          # Hajek estimator
    ni <- length(w)                                     # sample size in area i
    vi <- sum(w * (w - 1) * (x - avg)^2) / sum(w)^2     # variance
    c(avg = avg, vi = vi, ni = ni)
}
res <- sapply(datLCS_dom, function(u) hajek(u$income, u$w))
direct <- as.data.frame(t(res))
direct$dom <- as.numeric(rownames(direct))
# Coefficient of variation (in %)
direct$cv <- 100 * sqrt(direct$vi) / direct$avg

#-------------------------------------------------------------------------------
# TASK 2: Generalized variance function
#-------------------------------------------------------------------------------
est <- lm(log(vi) ~ avg * ni, data = direct)
p <- predict(est)
# Compute error varianc
sigma_e2 <- sum(residuals(est)^2) / df.residual(est)
# Variance of the direct estimator (using generalized variance function)
direct$vi_gvf <- exp(p + sigma_e2 / 2)
# Add auxiliary data to the data.frame 'direct'
direct <- cbind(direct, auxLCS[order(auxLCS$dom), ])

#-------------------------------------------------------------------------------
# TASK 3: Fay-Herriot model
#-------------------------------------------------------------------------------
# Estimate Fay-Herriot model (REML estimate of variance)
m <- eblupFH(avg ~ Mwork + Mnowork + Minact, vardir = vi_gvf, data = direct)
# Estimate Fay-Herriot model (without Mwork, because it is not significantly
# different from zero at 10% level of significance)
m <- eblupFH(avg ~ Mnowork + Minact, vardir = vi_gvf, data = direct)
# Estimate the second-order analytical approximation to the MSE
m <- mseFH(avg ~ Mnowork + Minact, vardir = vi_gvf, data = direct)
# Scatter plot of EBLUP against Hajek estimator
plot(direct$avg, m$est$eblup[, 1])
abline(0, 1)
# Line plot MSE and GVF variance
plot(direct$vi_gvf, type = "b")
lines(m$mse, type = "b", pch = 19)
legend("topleft", legend = c("vi_gvf", "MSE"), pch = c(1, 19), lwd = c(1, 1))
