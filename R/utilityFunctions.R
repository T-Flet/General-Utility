#### Generic Useful Functions

# Make a Cross-Validation table and produce a summary of it
xTWithSumm <- function(real, predicted) {
    xtab <- table(real, predicted)
    list(xtab = xtab, summary = xTSummary(xtab, levels(real)))
}

# Summary of a cross-validation table, with various True/False Positive/Negatives info
xTSummary <- function(xtab, labels = NULL) {
    d <- diag(xtab)
    t <- sum(xtab)
    rs <- rowSums(xtab)
    cs <- colSums(xtab)

    summary <- list(
        PerClass = tibble(
            Total = rs,
            Proportion = rs / t,
            # True Positive: of class X identified as X
            TP = d,
            TPPerClass = d / rs,
            # False Positive: of class Not-X identified as X 
            FP = cs - d,
            FPPerNotClass = (cs - d) / (t - rs),
            # True Negative: of class Not-X identified as Not-X
            TN = t - rs - cs + d, # Not row or column of X
            TNPerNotClass = (t - rs - cs + d) / (t - rs),
            # False Negative: of class X identified as Not-X
            FN = rs - d,
            FNPerClass = (rs - d) / rs
        ),
        Total = tibble(
            TP = sum(d),
            TPRate = sum(d) / t,
            All3Others = t - sum(d),
            All3OthersRate = (t - sum(d)) / t
        )
    )

    if (!is.null(labels)) {
        summary$PerClass = summary$PerClass %>% add_column(Class = labels, .before = 1)
    }

    summary
}



library(tidyverse)
## Runner function for multiple tuning parameter values
# Inputs:
#   n - number of tests to perform for each tuning parameters' combination
#   target - result to score against when comparing results
#   methodF - function wrapper around the method to test, with only the parameters to tune as arguments
#   params - (possibly named) list of tuning parameter values to test (all combinations will be produced)
# Output: a tibble of tuning parameters and result values, in ascending order of mean deviation from target
# NOTE: the functions in stochGradDesc.R are typical arguments to this function
tuningHelper <- function(n = 100, target, methodF, params) {
    pCombs <- do.call(crossing, params) %>% simplify2array
    map(1:nrow(pCombs), function(i) {
        map(1:n, function(unused) {
            out <- methodF(pCombs[i,])
            out$theta[nrow(out$theta),] - target
        }) %>% simplify2array %>% rowMeans %>% c(pCombs[i,], .)
    }) %>% simplify2array %>% t %>% as.tibble %>%
        mutate(mean = .[, 2:ncol(.)] %>% abs %>% rowMeans) %>%
        arrange(mean)
}



## Functions to go between Factor and Int
numToInt <- function(x) {
    as.integer(unlist(round(x)))
}

numToFac <- function(x, levels) {
    factor(levels[numToInt(x)], levels)
}

# Simple check:
#all(data == toClassFac(toClassInt(data)))



# Partition a data frame into either a given number of random subsets
# or into random subsets of given proportional sizes and names
# Set the argument 'parts' to an integer for the former and to a named list of proportions adding up to 1 for the latter
splitBy <- function(df, parts) {
    labels <- if (length(parts) == 1) {
        cut(seq(nrow(df)), parts, labels = 1:parts)
    } else {
        cut(seq(nrow(df)), nrow(df) * cumsum(c(0, parts)), labels = names(parts))
    }
    split(df, sample(labels))
}



# Plot the useful residual plots for a model on a pre-determined grid
# (Histogram, vs reponse and vs any explanatory of interest, provided in the ...)
library(tidyverse)
library(rlang)
library(gridExtra)
library(moderndive)
residPlotGrid <- function(mod, nrow, ncol, response, ...) {
    rExp <- sym(paste(response, "_hat", sep = ""))
    eExps <- ensyms(...)

    regPs <- get_regression_points(mod)
    baseGgp <- ggplot(regPs)

    hist <- baseGgp + geom_histogram(aes(residual), color = "white") + labs(x = "Residual")

    resp <- baseGgp + geom_point(aes(!!rExp, residual)) +
        labs(x = "Fitted Value", y = "Residual") +
        geom_hline(yintercept = 0, col = "blue", size = 1)

    expls <- map(eExps, function(eExp) {
        baseGgp +
        (if (is.factor(pluck(regPs, as_string(eExp))))
            geom_boxplot(aes(!!eExp, residual))
        else
            geom_point(aes(!!eExp, residual))
        ) +
        labs(x = as_string(eExp), y = "Residual") +
        geom_hline(yintercept = 0, col = "blue", size = 1)
    })

    do.call(grid.arrange, c(list(nrow = nrow, ncol = ncol, hist, resp), expls))
}
# E.g.: residPlotGrid(lm(Y ~ X1 + X2 + X3, data), 3, 2, "Y", "X1", "X2", "X3")