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

