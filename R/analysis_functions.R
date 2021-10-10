library(tidyverse)



### Factors ###


# Produce summary statistics of the response variable over the levels of a given factor
df_factor_lv_stats <- function(df, response_name, factor_name) {
  df %>%
    rename(Response = !!sym(response_name), Level = !!sym(factor_name)) %>%
    group_by(Level) %>%
    summarise(across(Response, list(
      Mean = ~ mean(., na.rm = F), SD = ~ sd(., na.rm = F),
      Q1 = ~ quantile(., 0.25, na.rm = F),
      Median = ~ median(., na.rm = F),
      Q3 = ~ quantile(., 0.75, na.rm = F)), .names = '{.fn}'))
}


# For each level of a factor, return its response variable's summary statistics along with the given model's coefficient details
glm_factor_coefs <- function(mod, factor_name) {
  method <- paste(deparse(mod$call), collapse = '') %>% str_extract('^([^\\(]+)')
  data <- mod[[ifelse(method %in% c('polr', 'multinom'), 'model', 'data')]]
  
  lvs <- levels(data[[factor_name]])
  lvs_long <- map(lvs, ~ paste0(factor_name, .)) %>% unlist()

  coefs <- summary(mod)$coefficients
  coefs <- as_tibble(coefs) %>%
    add_column(Var = rownames(coefs), .before = 1)

  intcp <- coefs %>% filter(Var == '(Intercept)')
  lvs_coefs <- coefs %>%
    filter(Var %in% lvs_long) %>%
    mutate(Estimate = Estimate + intcp$Estimate)#, `Std. Error` = `Std. Error` + intcp$`Std. Error`)

  res <- bind_rows(intcp, lvs_coefs) %>%
    mutate(Var = lvs) %>%
    select(-`t value`)
  colnames(res) <- c('Level', 'Coef', 'SE', 'P Val')
  
  factor_stats <- df_factor_lv_stats(data, as.character(mod$formula[[2]]), factor_name)
  inner_join(factor_stats, res, by = 'Level')
}



### Prediction ###


# Produce predictions with s.e. levels (put through the link function) from a glm and new data
predict_response <- function(mod, new_data) {
  method <- paste(deparse(mod$call), collapse = '') %>% str_extract('^([^\\(]+)')
  if (method == 'multinom' | method == 'polr') {
    ppreds <- predict(mod, newdata = new_data, type = 'probs')
    tib <- if (is.null(nrow(ppreds))) { as_tibble_row(ppreds) } else { as_tibble(ppreds) }
    tib <- rowid_to_column(tib)
    tib %>%
      gather(Pred, Prob, -1) %>% # Assume response is the first column
      group_by(rowid) %>% filter(Prob == max(Prob)) %>% ungroup() %>%
      arrange(rowid) %>%
      left_join(tib, by = 'rowid') %>% select(-rowid)
  } else {
    preds <- predict(mod, newdata = new_data, type = 'link', se.fit = T)
    critval <- qnorm(0.975) # ~ 95% CI
    tibble(Pred = preds$fit,
           Upper = preds$fit + (critval * preds$se.fit),
           Lower = preds$fit - (critval * preds$se.fit)) %>%
      mutate(across(everything(), mod$family$linkinv)) %>%
      mutate(TempUpper = Upper) %>% # These final lines are because the inverse link may not be monotonely increasing
      mutate(Upper = if_else(TempUpper >= Lower, TempUpper, Lower), Lower = if_else(TempUpper >= Lower, Lower, TempUpper)) %>%
      select(-TempUpper)
  }
}


# Produce predicted traces for each non-factor covariate in a model
predict_covariate_traces <- function(mod, new_data, n_points = 500) {
  terms <- attr(terms(formula(mod)), which = 'term.labels')
  empty_df <- new_data %>% select(!!terms) %>% slice(1) %>%
    mutate(across(terms, ~ ifelse(is.factor(.x), levels(.x)[1], mean(.x))))
  empty_df <- empty_df[rep(1, n_points),]
  res_dfs <- list()
  for (trm in terms) {
    if (is.factor(new_data[[trm]])) { # The baseline level gets the remainder of extra instances
      lvs <- levels(new_data[[trm]])
      range_col <- c(rep(lvs[1], n_points %% length(lvs)), rep(lvs, each = n_points %/% length(lvs)))
    } else { range_col <- seq(min(new_data[[trm]]), max(new_data[[trm]]), length.out = n_points) }
    res_dfs[[trm]] <- predict_response(mod, empty_df %>% mutate(!!trm := range_col)) %>% mutate(!!trm := range_col)
  }
  res_dfs
}



### Residuals ###


# Plot the useful residual plots for a model on a pre-determined grid
# (Histogram, vs response and vs any explanatory of interest, provided in the ...)
library(gridExtra)
library(moderndive)
residPlotGrid <- function(mod, nrow, ncol, response, ...) {
  rExp <- sym(paste(response, '_hat', sep = ''))
  eExps <- ensyms(...)

  regPs <- get_regression_points(mod)
  baseGgp <- ggplot(regPs)

  hist <- baseGgp + geom_histogram(aes(residual), color = 'white') + labs(x = 'Residual')

  resp <- baseGgp + geom_point(aes(!!rExp, residual)) +
    labs(x = 'Fitted Value', y = 'Residual') +
    geom_hline(yintercept = 0, col = 'blue', size = 1)

  expls <- map(eExps, function(eExp) {
    baseGgp +
    (if (is.factor(pluck(regPs, as_string(eExp)))) geom_boxplot(aes(!!eExp, residual))
     else geom_point(aes(!!eExp, residual)) ) +
    labs(x = as_string(eExp), y = 'Residual') +
    geom_hline(yintercept = 0, col = 'blue', size = 1)
  })

  do.call(grid.arrange, c(list(nrow = nrow, ncol = ncol, hist, resp), expls))
}
# E.g.: residPlotGrid(lm(Y ~ X1 + X2 + X3, data), 3, 2, 'Y', 'X1', 'X2', 'X3')



### Correlation ###


# Plot a variables' graph by their correlation (above a given threshold)
# Variables are positioned by multidimensional scaling on correlation magnitude,
# therefore placement and clustering reflect overall correlation profile similarity
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



### Heteroskedasticity ###


# Testing Homoscedasticity, Multivariate Normality, and Missing Completely at Random
# This function does tests for specific selections of variables otherwise could do a simple TestMCARNormality(df) for all together
# p_vals_only returns only the two tests' p-values for REJECTING MCAR (Hawkins' p-value assumes normality to reject MCAR); see TestMCARNormality documentation
#   Typical step after computation for p_vals_only: MCAR variables are RES %>% filter(hawkins > 0.05 | non_parametric > 0.05)
library(MissMech)
test_normal_MCAR <- function(df, always_include, include_individually, p_vals_only = F) {
  f <- function(col) { TestMCARNormality(df %>% select(all_of(c(always_include, col)))) }
  if (!p_vals_only) g <- f
  else g  <- function(col) { x <- f(col); list(hawkins = x$pvalcomb, non_parametric = x$pnormality) }
  
  x <- map(setNames(nm = include_individually), g)
  
  if (p_vals_only) bind_cols(names(x), bind_rows(x)) %>% rename(Variable_Set = `...1`) else x
}



### Cross-Validation ###


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

  if (!is.null(labels)) summary$PerClass = summary$PerClass %>% add_column(Class = labels, .before = 1)

  summary
}


