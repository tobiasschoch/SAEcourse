# 1. Utility function: shows estimated coefficients and significance
mysummary <- function(x)
{
    mat <- x$fit$estcoef
    colnames(mat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
    printCoefmat(mat, has.Pvalue = TRUE)
}

# 2. The function 'myFH', below, is identical to 'eblupFH' except that it adds
#    the slot 'residuals' to the return list
if ("sae" %in%.packages()) {
    myFH <- eblupFH
    body(myFH)[[14]] <- substitute({result$fit$estcoef <- coef;
        result$fit$residuals <- as.numeric(resid)} / (vardir + variance))
} else {
    cat("The 'sae' package must be loaded!",
        " First, load the package, then run this code snipped again\n")
}
# From here on, we can use 'myFH' in place of 'eblupFH' (if we like)
