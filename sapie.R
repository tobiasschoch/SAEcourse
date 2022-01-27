setwd("C:/My/Albanien_INSTAT/SALSTAT/saipe")
load("saipe.RData")
dat <- SAIPE2005[, c(1:4, 7, 8)]
names(dat) <- c("prIRS", "nfIRS", "prCensus", "yi", "vi", "state")
library(sae)
# 1. Utility function: shows estimated coefficients and significance
mysummary <- function(x)
{
    mat <- x$fit$estcoef
    colnames(mat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
    printCoefmat(mat, has.Pvalue = TRUE)
}

# 2. The function 'myFH', below, is identical to 'eblupFH' except that it adds
#    the slot 'residuals' to the return list
myFH <- eblupFH
body(myFH)[[14]] <- substitute({result$fit$estcoef <- coef;
    result$fit$residuals <- as.numeric(resid)} / (vardir + variance))
# From here on, we can use 'myFH' in place of 'eblupFH' (if we like)

#===============================================================================
# direct estimator
#===============================================================================
plot(sqrt(dat$vi) / dat$yi, cex.lab = 1.3, cex.axis = 1.2, xlab = "states",
    ylab = "coefficient of variation (in %)", pch = 19)
grid(col = "grey65")

#===============================================================================
# SAE
#===============================================================================
# data
head(dat, 3)

# QQ-plot of yi
qqnorm(dat$yi)
qqline(dat$yi)

#
m <- eblupFH(yi ~ prIRS + nfIRS + prCensus, vardir = vi, data = dat)
mysummary(m)


myFH(yi ~ prIRS + nfIRS + prCensus, vardir = vi, data = dat)

