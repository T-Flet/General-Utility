#### Generic Useful Functions

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
