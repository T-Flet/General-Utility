#### Generic Useful Functions


library(purrr)
# NOTE: Need to add the appropriate empty set manually, e.g. c(c(''), powerset(xs))
powerset <- function(xs) { map(1:(length(xs)), ~ combn(xs, ., simplify = F)) %>% unlist(recursive = F) }


# Plot a variables' graph by their correlation (above a given threshold)
library(corrr)
cor_network <- function(df, threshold = 0.5, upper_threshold = 1, repel = T, curved = T) {
  cors <- cor(df, use = 'pairwise.complete.obs')
  cond <- (threshold < abs(cors) & abs(cors) < upper_threshold) & cors != 1
  cor_mat <- cors[rowSums(cond) > 0, colSums(cond) > 0, drop = F]
  print(paste(nrow(cor_mat), 'variables with correlations between', threshold, 'and', upper_threshold))
  network_plot(as_cordf(cor_mat), min_cor = threshold, repel = repel, curved = curved, colours = c('tomato', 'white', 'royalblue'))
}
cor_hist <- function(df) {
  cors <- cor(df, use = 'pairwise.complete.obs')
  diag(cors) <- NA
  tibble(Correlation = as.vector(cors)) %>% ggplot() +
    geom_histogram(aes(Correlation), bins = 50) +
    ylab('Count')
}
# cor_hist(covid$base %>% select(where(is.numeric)))


# Testing Homoscedasticity, Multivariate Normality, and Missing Completely at Random
# This function does tests for specific selections of variables otherwise could do a simple TestMCARNormality(df) for all together
# p_vals_only returns only the two tests' p-values for REJECTING MCAR (Hawkins' p-value assumes normality to reject MCAR); see TestMCARNormality documentation
#   Typical step after computation for p_vals_only: MCAR variables are RES %>% filter(hawkins > 0.05 | non_parametric > 0.05)
library(MissMech)
test_normal_MCAR <- function(df, always_include, include_individually, p_vals_only = F) {
  f <- function(col) { TestMCARNormality(df %>% select(all_of(c(always_include, col)))) }
  if (!p_vals_only) { g <- f }
  else { g  <- function(col) { x <- f(col); list(hawkins = x$pvalcomb, non_parametric = x$pnormality) } }
  
  x <- map(setNames(nm = include_individually), g)
  
  if (p_vals_only) { bind_cols(names(x), bind_rows(x)) %>% rename(Variable_Set = `...1`) } else { x }
}



# Run multiple and select the best of dineof data imputations (in preparation for eof of some variety)
# Dineof imputes data by a form of iterated eof
# Comparison with other missing-data EOFs: https://menugget.blogspot.com/2014/09/pca-eof-for-data-with-missing-values.html
best_dineof <- function(predictor_df, times) {
  best <- NULL
  scores <- list()
  for (i in 1:times) {
    done <- dineof(predictor_df %>% as.matrix())
    scores[[as.character(i)]] <- paste(done$n.eof, tail(done$RMS, 1))
    if (is.null(best) || (tail(done$RMS, 1) < tail(best$RMS, 1))) { best <- done }
  }
  list(best = best, scores = scores)
}

# Multiple points of analysis of princiapl component (PC) or empirical orthogonal function (EOF) analysis:
#   - Number of required components to reach a given explained variance threshold
#   - Bar-plot of component variances
#   - Vertical aligned barplots of loadings of the first few components
#   - Vertical barplots of most contributive variables in the first few components
library(tidytext)
eof_analysis <- function(variances, rotation, var_names, threshold = 0.9, extra_variances_to_plot = 2, comps_to_plot = 4, top_comp_variables = 8) {
  res <- list()
  
  cumulatives <- cumsum(variances / sum(variances))
  res$required <- which(cumulatives > threshold)[1]
  
  cs <- fct_inorder(paste0('C', 1:nrow(rotation)))
  threshold_cs <- c(rep.int(T, res$required), rep.int(F, ncol(rotation) - res$required))
  
  res$variance_plot <- tibble(Component = cs, Variance = variances, Below_Threshold = threshold_cs) %>%
    head(res$required + extra_variances_to_plot) %>% ggplot() +
    geom_col(aes(Component, Variance, fill = Below_Threshold), show.legend = F)
  
  lps <- eof_loadings_plots(rotation, var_names, cs, comps_to_plot = comps_to_plot, top_comp_variables = top_comp_variables)
  res$loadings_plot <- lps$loadings_plot
  res$highest_loadings_plot <- lps$highest_loadings_plot

  res
}

# Further points of analysis of princiapl component (PC) or empirical orthogonal function (EOF) analysis:
#   - Perform a glm and a refined glm (by step_back from iterated_models.R, which NEEDS to be imported)
#   - Produce the same loadings plots as eof_analysis for the refined model's highest p-value eofs
eof_glm_analysis <- function(resp_name, resp_v, rotated_data, required, variances, rotation, var_names, comps_to_plot = 4, top_comp_variables = 8, threshold = 0.1) {
  res <- list()

  df <- as_tibble(rotated_data[,1:required]) %>%
    setNames(paste0('C', 1:required)) %>%
    mutate(!!ensym(resp_name) := resp_v, .before = 1)
  
  res$simple_glm <- glm(as.formula(paste(resp_name, '~ .')), df, family = gaussian())
  res$step_back_glm <- stepLeastSigGLM(res$simple_glm, threshold = threshold, significance_statistic = 't', test = 'F')$mod
  
  rel_vars <- tibble(varName = paste0('C', 1:nrow(rotation)), Rel_Var = variances / sum(variances))
  coeffs <- summary(res$step_back_glm)$coefficients
  res$ord_cs <- as_tibble(coeffs) %>%
    add_column(varName = coeffs %>% rownames, .before = 1) %>%
    filter(varName != '(Intercept)') %>%
    arrange(`Pr(>|t|)`) %>%
    left_join(rel_vars, by = 'varName')

  lps <- eof_loadings_plots(rotation, var_names, cs_subset = res$ord_cs$varName, comps_to_plot = comps_to_plot, top_comp_variables = top_comp_variables)
  res$loadings_plot <- lps$loadings_plot
  res$highest_loadings_plot <- lps$highest_loadings_plot
  
  res
}

# EOF loadings plots as used and described in eof_analysis and eof_glm_analysis
eof_loadings_plots <- function(rotation, var_names, cs = fct_inorder(paste0('C', 1:nrow(rotation))), cs_subset = cs, comps_to_plot = 4, top_comp_variables = 8) {
  res <- list()
  
  df <- rotation %>% t() %>% as_tibble() %>% setNames(var_names) %>%
    mutate(Component = cs, .before = 1) %>%
    pivot_longer(-Component)
  
  res$loadings_plot <- df %>%
    filter(Component %in% cs_subset[1:comps_to_plot]) %>%
    mutate(name = factor(name, levels = sort(unique(name), T))) %>%
    ggplot() + geom_col(aes(value, name, fill = name), show.legend = F) +
      facet_wrap(~ Component, nrow = 1) +
      labs(x = NULL, y = NULL)
  
  res$highest_loadings_plot <- df %>%
    mutate(Sign = factor(sign(value))) %>%
    filter(Component %in% cs_subset[1:comps_to_plot]) %>%
    group_by(Component) %>%
    top_n(top_comp_variables, abs(value)) %>%
    ungroup() %>%
    mutate(name = reorder_within(name, abs(value), Component)) %>%
    ggplot() + geom_col(aes(abs(value), name, fill = Sign), show.legend = F) +
      facet_wrap(~ Component, scales = 'free_y') +
      scale_y_reordered() + scale_fill_manual(values = c('tomato', 'royalblue')) +
      labs(x = NULL, y = NULL)
  
  res
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
numToInt <- function(x) { as.integer(unlist(round(x))) }

numToFac <- function(x, levels) { factor(levels[numToInt(x)], levels) }

# as.numeric is guaranteed to return the level number even if the level is an integer itself
facToInt <- function(x) { as.numeric(x) } # Old: setNames(c(1:length(levels(x))), nm = levels(x))[x] # Older: as.numeric(levels(x))[x]

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



# Compute the overlap of two intervals
  # IMPORTANT NOTE: This is not a vectorised function, therefore DO
  #   EITHER: rowwise %>% mutate(... intervalOverlap(...) ...) %>% ungroup
  #   OR: mutate(... Vectorize(intervalOverlap)(...) ...)
  # If this is not done it will return the result of the first evaluation and that will be it for all rows
intervalOverlap <- function(a, b) {
  max(0, min(a[2], b[2]) - max(a[1], b[1]))
}
