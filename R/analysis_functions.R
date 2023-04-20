library(moments)
library(entropy)
library(fitdistrplus)
library(mclust)
library(metRology)
library(GGally)
library(plotly)
library(tidyverse)



### Summary ###


# Produce summary statistics for a list of vectors (a df also works as input; columns are taken to be said vectors)
  # The distributions argument should be F, NULL or a list of the outputs of the distribution_identifier function for each input variable:
  #   for F no additional column is added; for NULL the distribution_identifier function is called for each variable and
  #   a column for the str value of the best distribution is added; same if an output is provided directly
batch_summaries <- function(data, distributions = NULL) {
  data <- map(data, ~ discard(., is.na)) # Better here than for each case
  
  res <- keep(data, is.numeric) %>% map(function(v) {list(
    Mean = mean(v), SD = sd(v), Var = var(v),
    Skewness = skewness(v), Kurtosis = kurtosis(v),
    Q1 = quantile(v, 0.25), Median = median(v), Q3 = quantile(v, 0.75)
  )})
  res <- bind_rows(res) %>% add_column(Name = names(res), .before = 1)
  if (nrow(res) != 0) res <- rowwise(res) %>%
    mutate(Outliers = data[[Name]] %>% keep(function(x) (x < Q1-1.5*(Q3-Q1)) || (x > Q3+1.5*(Q3-Q1))) %>% length())
  
  factor_res <- discard(data, is.numeric) %>% map(function(f) {
    freqs <- table(f) / length(f)
    mode_lvl <- which.max(freqs)
    bit_entropy <- entropy.empirical(freqs, unit = 'log2')
    list(N_Levels = length(freqs), Mode = names(freqs[mode_lvl]), Mode_Freq = unname(freqs[mode_lvl]),
         Bit_Entropy = bit_entropy, Scaled_Entropy = bit_entropy / log2(length(freqs)))#, value_proportions = as.list(freqs))
  })
  factor_res <- bind_rows(factor_res) %>% add_column(Name = names(factor_res), .before = 1)
  res <- bind_rows(res, factor_res)

  if (isFALSE(distributions)) res else {
    if (is.null(distributions)) {
      distributions <- keep(data, is.numeric) %>% map(distribution_identifier)
      distributions <- c(distributions, discard(data, is.numeric) %>% map(distribution_identifier))
    }
    res %>% add_column(Distribution = map(distributions, ~ .[[1]]$str) %>% unlist())
  }
}


# Generate the data and optionally plot a distribution
#   NOTE: the range argument should be made to be integer for discrete distributions,
#           otherwise a straight line at 0 is likely to be produced
get_density <- function(name_with_no_d, params, range = seq(-10, 10, length.out = 100), plot = F) {
  pdf <- function(x) do.call(paste0('d', name_with_no_d), c(list(x = x), params))
  df <- tibble(x = range, pdf = pdf(x))
  if (plot) ggplot(df) + geom_line(aes(x, pdf)) + ggtitle(paste0(name_with_no_d, '(', paste(map(params, ~ round(., 3)), collapse = ', '), ')'))
  else df
}


# Generate tags defining the domain of the given data
domain_tags <- function(xs) {
  xs <- xs %>% discard(is.na)
  uniq <- unique(xs)
  poss_tags <- list()
  if (is.numeric(xs)) {
    poss_tags <- c(poss_tags, list(
      `Non-Negative` = min(xs) >= 0,
      `Non-Positive` = max(xs) <= 0,
      `Non-Zero` = all(xs != 0),
      `[0,1]` = 0 <= min(xs) & max(xs) <= 1
    ))
    poss_tags <- c(poss_tags, if (all(as.integer(xs) == xs)) { list(Integer = T) } else { list(Real = T) })
  } else {
    poss_tags <- c(poss_tags, list(
      Factor = T,
      Binary = length(uniq) == 2
    ))
    if (is.logical(xs)) {
      poss_tags <- c(poss_tags, list(Boolean = T))
    }
  }
  poss_tags %>% keep(identity) %>% names()
}


# Identify and fit a distribution to the given data (domain tags are generated if not known)
#   NOTE: by default if the best identified distribution is a uniform one then a finite Gaussian
#         mixture model is attempted (and discarded if only one is the best fit).
#         This is noteworthy because the object type of the distr field of one element of the
#         might be different from all others, breaking functions equipped to only handle fitdist
#         outputs. This behaviour can be prevented with try_gauss_mix_on_unif = F.
#   Also note that if the given data is a factor no traces are returned even for density_traces = T.
#   Also also note that without unlimited_gauss_mix_n = T mixtures will only be accepted if between 2 and 4 inclusive.
distribution_identifier <- function(xs, tags = domain_tags(xs), forced_dists = NULL, density_traces = F, try_gauss_mix_on_unif = T, unlimited_gauss_mix_n = F) {
  distr_doms <- list(
    beta = c('[0,1]', 'Real'),
    binom = c('Non-Negative', 'Integer'), # No need for Bernoulli and 'Binary' or 'Boolean'
    cauchy = c(),
    chisq = c('Non-Negative'),
    exp = c('Non-Negative'),
    f = c('Non-Negative'),
    gamma = c('Non-Negative', 'Non-Zero'),
    geom = c('Non-Negative', 'Non-Zero', 'Integer'),
    hyper = c('Non-Negative', 'Integer'),
    lnorm = c('Non-Negative'),
    # multinom = c('Non-Negative', 'Integer'), # Here as well
    nbinom = c('Non-Negative', 'Integer'),
    norm = c(),
    pois = c('Non-Negative', 'Integer'),
    t.scaled = c(), # Using t.scaled from metRology because the default t is just bad
    unif = c(),
    weibull = c('Non-Negative')
  )

  xs <- xs %>% discard(is.na)
  if ('Factor' %in% tags) {
    freqs <- table(xs) / length(xs)
    return(list(cat = list(distr = 'No object; use dcat, qcat, rcat if necessary', bic = 0, str = paste0('cat(k = ', length(freqs), ', ps = {', paste(map(freqs, ~ round(.,3)), collapse = ', '), '})'))))
  } else {
    if ('Non-Positive' %in% tags) {
      xs <- -xs
      tags <- tags %>% keep(~ . != 'Non-Positive') %>% append('Non-Negative')
    }
    distributions <- map(distr_doms, ~ setdiff(., tags)) %>%
      keep(~ length(.) == 0) %>% names() %>% setNames(nm = .)
    
    # Starting values or fixed arguments may be required; reasonable results can be obtained even from arbitrary starting values
    params <- list(
      ############## TO DO: for binom, remove fix.arg from here and instead try a few fixed values
      # in the res_distrs map
      #   (by checking distr == 'binom' and then trying, say, from max(xs) to round(1.5 * max(xs))
      binom = list(fix.arg = list(size = max(xs)), start = list(prob = 0.5), method = 'mse'),
      multinom = list(start = list(size = max(xs), prob = 1 / max(xs))),
      chisq = list(start = list(df = 3)),
      f = list(start = list(df1 = 3, df2 = 3)),
      hyper = list(start = list(m = max(xs), n = max(xs) / 2, k = 1.5 * max(xs))),
      t.scaled = list(start = list(df = 3, mean = mean(xs), sd = sd(xs)))
    )
  }
  
  if (!is.null(forced_dists)) distributions <- forced_dists %>% setNames(nm = .)

  # Produce a model name with parameters for individual fitdistr outputs
  get_distr_str <- function(distr) {
    all_params <- if ('fix.arg' %in% names(distr)) c(distr$fix.arg, distr$estimate) else distr$estimate
    paste0(distr$distname, '(',
      paste(map(names(all_params), function(nm)
        paste(nm, '=', round(all_params[[nm]], 3))), collapse = ', ' ), ')')
  }
  
  # Simple distributions
  res_distrs <- map(distributions, function(distr) { tryCatch(
      do.call(fitdist, c(list(data = xs, distr = distr), if (distr %in% names(params)) params[[distr]] else list())),
      error = function(e) as.character(e) )}) %>%
      discard(~ length(.) == 1) %>% discard(~ any(is.na(.$estimate))) %>%
    map(function(res) list(distr = res, bic = res$bic, str = get_distr_str(res))) %>% sort_by_elem(2)

  # Gaussian finite mixture
  if (try_gauss_mix_on_unif && (!'Factor' %in% tags) && (names(res_distrs)[1] == 'unif')) {
    finite_mix <- densityMclust(xs, plot = F)
    if (finite_mix$G > 1 && (unlimited_gauss_mix_n || finite_mix$G <= 4))
      res_distrs <- c(res_distrs, list('gauss_mix' = 
        list(distr = finite_mix, bic = finite_mix$bic, str = get_mclust_str(finite_mix), n_mix = finite_mix$G)
      )) %>% sort_by_elem(2)
  }
   
  if (density_traces) map(res_distrs, function(rd) {
    data_range <- if ((!'n_mix' %in% names(rd)) && 'Integer' %in% distr_doms[[rd$distr$distname]]) seq(min(xs), max(xs)) else seq(min(xs), max(xs), length.out = 100)
    if ('n_mix' %in% names(rd)) c(rd, tibble(x = data_range, pdf = predict.densityMclust(finite_mix, newdata = data_range)))
    else c(rd, get_density(rd$distr$distname, {
        if ('fix.arg' %in% names(rd$distr)) c(rd$distr$fix.arg, rd$distr$estimate) else rd$distr$estimate
      }, data_range)) })
  else res_distrs
}
# distribution_identifier(rbinom(50, 10, 0.75), forced_dists = c('binom'), try_gauss_mix_on_unif = T)
# map(list(test = rbinom(50, 10, 0.75)), ~ distribution_identifier(., forced_dists = c('binom'), density_traces = T, try_gauss_mix_on_unif = F))
# distribution_identifier(rhyper(50, 10, 10, 10), forced_dists = c('hyper'))
# map(list(test = rhyper(50, 10, 10, 10)), ~ distribution_identifier(., forced_dists = c('hyper'), density_traces = T, try_gauss_mix_on_unif = F))
# distribution_identifier(iris$Sepal.Length, forced_dists = c('norm'))


# Produce a model name with parameters for a given (mclust) Gaussian finite mixture fit
get_mclust_str <- function(mod) {
  props <- if ('pro' %in% names(mod$parameters)) mod$parameters$pro else rep(1/mod$G, mod$G)
  paste(map(seq(mod$G), function(i) paste0(
    round(props[i], 3), ' ', 'norm(mean = ', round(mod$parameters$mean[i], 3), ', variance = ', round(mod$parameters$variance$sigmasq[i], 3), ')'
  )), collapse = ' + ') }
# get_mclust_str(densityMclust(c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 2, sd = 0.5))))


# Python-compatible Plotly JSON for a given ggplot, e.g. ggpairs(df)
  # Although https://github.com/plotly/plotly.R/issues/2009, an unrecognised 'key' field can still occur,
  # therefore use this on the Python side: plotly.io.from_json(plotly_json, skip_invalid = True)
ggplotly_json <- function(gg_plot, even_safer = T) {
  pl <- ggplotly(gg_plot)
  (if (even_safer) plotly:::to_JSON(plotly_build(pl)$x[c('data', 'layout', 'config')]) else plotly_json(pl, F, F)) %>%
    str_replace_all(fixed('transparent'), '#ffffff') %>%
    str_replace_all(fixed(',"frame":null'), '') %>%
    str_replace_all(fixed('"family":"",'), '') %>%
    str_replace_all(fixed(',"linetype":[]'), '')
}



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


# Produce predictions with s.e. levels (put through the link function) from a glm and new data.
# The remaining arguments trigger additional outputs beside the raw prediction (and interval) for non-categorical models:
#   bounds for a given probability (to be interpreted as central or not)
#   or probabilities of values being above or below a given bound.
# Specifically:
#   If a custom_prob is given, then the Custom_Lower and Custom_Upper columns are added
#     - if prob_as_central, then custom_prob is centred and the two columns are the bounds of the corresponding CI
#     - if not, then the two columns are each a bound for custom_prob of the values being above or below them
#   If a custom_bound is given then the Prob_Below_Bound and Prob_Above_Bound columns are added
predict_response <- function(mod, new_data, custom_prob = NULL, prob_as_central = F, custom_bound = NULL) {
  method <- paste(deparse(mod$call), collapse = '') %>% str_extract('^([^\\(]+)')
  if (method == 'multinom' | method == 'polr') {
    ppreds <- predict(mod, newdata = new_data, type = 'probs')
    if (length(mod$lev) == 2) { # multinom returns a bare list of probs if only 2 levels; polr works only from 3 up
      tib <- tibble(XX = ppreds) %>% mutate(YY = 1 - XX)
      colnames(tib) <- mod$lev
    } else {
      tib <- if (is.null(nrow(ppreds))) { as_tibble_row(ppreds) } else { as_tibble(ppreds) }
    }
    tib <- rowid_to_column(tib)
    tib %>%
      gather(Pred, Prob, -1) %>% # Assume response is the first column
      group_by(rowid) %>% filter(Prob == max(Prob)) %>% ungroup() %>%
      arrange(rowid) %>%
      left_join(tib, by = 'rowid') %>% select(-rowid)
  } else {
    preds <- predict(mod, newdata = new_data, type = 'link', se.fit = T)
    res <- tibble(Pred = preds$fit, # Central 95% CI below, i.e. same as preds$fit +/- (qnorm(0.975) * preds$se.fit)
      Lower = qnorm(0.975, preds$fit, preds$se.fit, lower.tail = F),
      Upper = qnorm(0.975, preds$fit, preds$se.fit))
    if (!is.null(custom_prob)) {
      custom_prob <- if (prob_as_central) (1 + custom_prob) / 2 else custom_prob # I.e. 1 - (1 - custom_prob) / 2 if it is meant to be central
      res <- res %>% mutate( # Custom bounds
        Custom_Lower = qnorm(custom_prob, preds$fit, preds$se.fit, lower.tail = F),
        Custom_Upper = qnorm(custom_prob, preds$fit, preds$se.fit))
    }
    
    res <- res %>%
      mutate(across(everything(), mod$family$linkinv)) %>%
      mutate(TempLower = Lower) %>% # These final lines are because the inverse link may not be monotonically increasing
      mutate(Lower = if_else(TempLower >= Upper, Upper, TempLower),
             Upper = if_else(TempLower >= Upper, TempLower, Upper)) %>%
      select(-TempLower)
    if (!is.null(custom_prob)) res <- res %>% mutate(TempLower = Custom_Lower) %>% # same as above
      mutate(Custom_Lower = if_else(TempLower >= Custom_Upper, Custom_Upper, TempLower),
             Custom_Upper = if_else(TempLower >= Custom_Upper, TempLower, Custom_Upper)) %>%
      select(-TempLower)

    if (!is.null(custom_bound)) {
      linear_bound <- mod$family$linkfun(custom_bound)
      res <- res %>% mutate(
               Prob_Below_Bound = pnorm(linear_bound, preds$fit, preds$se.fit)) %>%
        mutate(Prob_Above_Bound = 1 - Prob_Below_Bound) %>% # These lines below are because the inverse link may not be monotonically increasing
        mutate(TempLinkIsInverting = xor(linear_bound > preds$fit, custom_bound > Pred)) %>%
        mutate(Prob_Below_Bound = if_else(TempLinkIsInverting, 1 - Prob_Below_Bound, Prob_Below_Bound),
               Prob_Above_Bound = if_else(TempLinkIsInverting, 1 - Prob_Above_Bound, Prob_Above_Bound)) %>%
        select(-TempLinkIsInverting)
    }
    
    res
  }
}
# tibble(iris) %>% head()
# predict_response(glm(Sepal.Length ~ ., data = iris), new_data = tibble(iris) %>% head(),
#                  custom_prob = 0.9, prob_as_central = F, custom_bound = 5)
# predict_response(glm(Sepal.Length ~ ., data = iris, family = Gamma()), new_data = tibble(iris) %>% head(),
#                  custom_prob = 0.9, prob_as_central = F, custom_bound = 5)


# Produce df of new evenly spaced observations in the range of a given one
#   NOTE: n_points is just for reference, and possibly MANY more points may be generated
#           (because the n-th root of the number of numerical columns is ROUNDED UP)
get_range_grid <- function(df, n_points = 1000) {
  axes <- list()

  factors <- keep(colnames(df), ~ is.factor(df[[.x]]))
  if (length(factors) > 0) {
    for (fct in factors) { # Factors get all their levels as grid projection regardless of their number
      axes[[fct]] <- levels(df[[fct]])
      n_points <- n_points / length(axes[[fct]])
    }
    numericals <- setdiff(colnames(df), factors)
  } else numericals <- colnames(df)

  axis_length <- ceiling(n_points ^ (1 / length(numericals)))
  for (num in numericals) axes[[num]] <- seq(min(df[[num]]), max(df[[num]]), length.out = axis_length)

  do.call(expand_grid, axes)
}


# Produce predictions at (new) evenly spaced observations in the range of a given model's data
get_grid_predictions <- function(mod, response, n_points = 1000, keep_covariate_name = F) {
  df <- if (response %in% colnames(mod$data)) mod$data %>% select(-!!response) else mod$data
  df_grid <- get_range_grid(df, n_points = n_points)
  predictions <- predict_response(mod, df_grid)

  res <- map(setNames(nm = colnames(df_grid)), ~ predictions %>% mutate(Covariate = df_grid[[!!.x]]))
  if (keep_covariate_name) res <- map(setNames(nm = names(res)), ~ res[[.x]] %>% rename(!!.x := Covariate))
  res
}
# mod <- glm(Petal.Length ~ Species + Sepal.Length + Sepal.Length:Species + I(Sepal.Length^2) + Sepal.Width + Petal.Width + Sepal.Width:I(Petal.Width^2), iris, family = Gamma())
# get_grid_predictions(mod, 'Petal.Length', n_points = 10000, keep_covariate_name = F)$Sepal.Length
# get_grid_predictions(mod, 'Petal.Length', n_points = 1000, keep_covariate_name = T)$Sepal.Length


# Produce predicted traces for each covariate in a model by setting all other covariates to their mean (median level for factors)
get_covariate_traces <- function(mod, data = NULL, n_points = 200, keep_covariate_name = F) {
  response <- as.character(formula(mod)[[2]])
  if (is.null(data)) data <- mod$data
  df <- if (response %in% colnames(data)) data %>% select(-!!response) %>% as_tibble() else data
  empty_df <- summarise(df, across(.fns = ~ ifelse(is.factor(.x), levels(.x)[ceiling(length(levels(.x)) / 2)], mean(.x, na.rm = T)))) %>%
    slice(rep(1, n_points))

  res <- map(setNames(nm = colnames(df)), function(cn) {
    range_col <- if (is.factor(df[[cn]])) levels(df[[cn]]) else range_col <- seq(min(df[[cn]], na.rm = T), max(df[[cn]], na.rm = T), length.out = n_points)
    predict_response(mod, empty_df %>% head(length(range_col)) %>% mutate(!!cn := range_col)) %>% mutate(Covariate := range_col)
  })
  if (keep_covariate_name) res <- map(setNames(nm = names(res)), ~ res[[.x]] %>% rename(!!.x := Covariate))
  res
}
# mod <- glm(Petal.Length ~ Species + Sepal.Length + Sepal.Length:Species + I(Sepal.Length^2) + Sepal.Width + Petal.Width + Sepal.Width:I(Petal.Width^2), iris, family = Gamma())
# get_covariate_traces(mod, n_points = 100, keep_covariate_name = F)
# get_covariate_traces(mod, n_points = 50, keep_covariate_name = T)



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


