#===============================================================================
# 0 Preparations
#===============================================================================
# Location
setwd("C:/My/code/SAEcourse")

# Load the data
source("data.R")
head(dat, 3)

#===============================================================================
# 1 Direct estimator y_i
#===============================================================================
# Plot coefficient of variation of the direct estimator
plot(sqrt(dat$vi) / dat$yi, cex.lab = 1.3, cex.axis = 1.2, xlab = "states",
    ylab = "coefficient of variation (in %)", pch = 19)
grid(col = "grey65")

# QQ-plot of direct estimator yi
qqnorm(dat$yi)
qqline(dat$yi)

#===============================================================================
# 2 R package sae
#===============================================================================
# Package
library(sae)

# Utility functions
source("methods.R")

# Formula of the model
f <- yi ~ prIRS + nfIRS + prCensus

# Fay-Herriot model
m <- eblupFH(f, vardir = vi, data = dat)
mysummary(m)
m$fit$goodness
# a different model (without explanatory variable 'prIRS')
eblupFH(yi ~ nfIRS + prCensus, vardir = vi, data = dat)$fit$goodness

# Fit of the same model using the 'myFH' function => returns the standardized
# residuals as an additional slot of 'fit'
tmp <- myFH(f, vardir = vi, data = dat)

# QQ-plot of the standardized resiudals
qqnorm(tmp$fit$residuals)
qqline(tmp$fit$residuals)

# Plot EBLUP vs. direct estimator (with 45 degree-line)
plot(dat$yi, m$eblup[,1], xlab = "direct", ylab = "EBLUP")
abline(0, 1)

# Mean square error (MSE) estimation
# 1) Analytic second-order approximation of MSE (Prasad & Rao, 1990)
mse_analytic <- mseFH(f, vardir = vi, data = dat)
# 2) Bootstrap estimation of MSE with B = 500 replications
mse_boostrap <- mseFH(f, vardir = vi, B = 500, data = dat)

# Plot of analytic MSE and variance of direct estimator by area
orange <- rgb(212, 97, 10, maxColorValue = 255)
blue <- rgb(51, 51, 153, maxColorValue = 255)
plot(dat$vi, type = "b", xlab = "area", ylab = "uncertainty (MSE and variance)",
    col = blue, ylim = c(0, max(dat$vi)), lwd = 2, cex.axis = 1.2,
    cex.lab = 1.3)
grid(col = "grey65")
lines(mse_analytic$mse, col = orange, type = "b", pch = 19, lwd = 2)
legend("bottomright", legend = c("direct estimator", "EBLUP"),
    col = c(blue, orange), pch = c(1, 19), lwd = c(2, 2), cex = 1.2)

# Normal-theory 95% confidence intervals of the EBLUP (using estimated MSE)
# We generate a new data.frame
est <- data.frame(yi = dat$yi, EBLUP = m$eblup[, 1])
# Alpha: 5% level of significance (= 95% confidence level)
alpha <- 0.05
# Add the CI to the new data.frame
est$ci_low <- est$EBLUP - sqrt(mse_analytic$mse) * qnorm(1 - alpha/2)
est$ci_high <- est$EBLUP + sqrt(mse_analytic$mse) * qnorm(1 - alpha/2)
head(est, 3)

#===============================================================================
# 3 R package emdi
#===============================================================================
# Package
library(emdi)

# Fay-Herriot model (fitted by REML) without MSE
m <- fh(yi ~ prIRS + nfIRS + prCensus, vardir = "vi", combined_data = dat,
    method = "reml")

# Retrieve data.frame with direct estimator and EBLUP
head(estimators(m), 3)

# Summary method
summary(m)

# Diagnostic plots (QQ- and density plots of residuals and random effects)
plot(m)

# Fay-Herriot model (fitted by REML) with Prasad-Rao second order analytic
# approximation of MSE
m <- fh(fixed = yi ~ prIRS + nfIRS + prCensus, vardir = "vi",
    combined_data = dat, method = "reml"
    MSE = TRUE, mse_type = "analytical")

# Plot of direct estimator vs. EBLUP and MSE by area
compare_plot(m)

# Miscellaneous
# map_plot(m) # not shown (maps and shape files must be loaded before call)

# Excel output
write.excel(m)
