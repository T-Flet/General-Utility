library(stats)
library(tidyverse)



### Joint PCA ###


# Number of PCA components required to achieve varProp cumulative variance proportion of PCA result full_PCA
cumVar_index <- function(full_PCA, varProp) {
    cumVarProp <- summary(full_PCA)$importance[3,]
    which(cumVarProp > varProp)[1]
}


# Perform joint PCA of training and test data, producing multiple pairs of PCAd datasets according to the
# cumulative variance proportions specified in varProps.
# The options are whether to standardise the data before PCA, whether to plot the components and print the summary
# and whether to include the PCA result object in the output
jointPCA <- function(train, test, varProps = c(0.9), center = T, scale = T, verbose = T, return_full_PCA_data = F) {
    varProps <- sort(varProps)
    names(varProps) <- varProps

    both <- bind_rows(train, test)
    #all.equal(both[1:nrow(train),], train)
    #all.equal(both[(nrow(train) + 1):nrow(both),], test)

    both_PCA <- prcomp(both, center = center, scale. = scale)

    if (verbose) { print(summary(both_PCA)); plot(both_PCA) }

    res <- list(varSets = list())
    if (return_full_PCA_data) { res$full_PCA <- both_PCA }

    train_PCAd <- both_PCA$x[1:nrow(train),]
    test_PCAd <- both_PCA$x[(nrow(train) + 1):nrow(both),]
    
    res$varSets <- map(varProps, function(varProp) {
        ind <- cumVar_index(both_PCA, varProp)
        list(train = train_PCAd[, 1:ind], test = test_PCAd[, 1:ind])
    })
    
    res
}



### Generic Empirical Orthogonal Functions (EOF) Analysis ###


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


