library(glmnet)
library(vip)
library(doParallel)
doParallel::registerDoParallel()

library(tidyverse)
library(tidymodels)


# Use notes:
#   - The mixture and levels arguments should be both integers or 'tune' and (lower than otherwise) integer or 'tune' and integer 2-vector (with 2nd value lower than the 1st)
#   - df should, at least for now, be inputed as.data.frame or as.matrix to avoid workflows bug
#   - df should not contain any variables beside response and predictors (i.e. remove other identifiers)
#   - drop_na() is recommended for the df
# NOTE: This function avoids the currently buggy recipe -> workflow pattern for single model fit
#         Bug with temporary solution: https://community.rstudio.com/t/issue-with-tidymodels-workflows-and-fitting-xgboost-models/75015/11
regreg_grid <- function(fmla, df, family = 'gaussian', mixture = 'tune', levels = c(50, 5),
                        rec = recipe(fmla, data = df) %>% update_role(all_nominal(), has_type('logical'), new_role = 'ID'),
                        boots_n = 10, boots_strata = NULL, seed = 42) {
  res <- list()
  df <- as.data.frame(df) # Because of bug in description
  
  boots <- bootstraps(df, times = boots_n, strata = all_of(boots_strata))
  set.seed(seed)
  res$fit_grid <- if (mixture == 'tune') { tune_grid(
    linear_reg(penalty = tune(), mixture = tune()) %>% set_engine('glmnet', family = family),
    preprocessor = rec, resamples = boots, grid = grid_regular(penalty(), mixture(), levels = levels)
  )} else { tune_grid(
    linear_reg(penalty = tune(), mixture = mixture) %>% set_engine('glmnet', family = family),
    preprocessor = rec, resamples = boots, grid = grid_regular(penalty(), levels = levels)
  )}
  
  res$misc <- list(fmla = fmla, df = df, family = family, mixture = mixture)
  res$misc$resp_and_preds <- rec$term_info %>% filter(role == 'predictor' | role == 'outcome') %>% pluck('variable')

  res
}


# Follow-up to tidy_regreg: after the parameter grid has been explored
# Use notes:
#   - Can use select_best instead of select_by_one_std_err for the numerically-best result
regreg_select <- function(res, selection_crit = select_by_one_std_err, metric = 'rmse', top_n_vars = 20) {
  # Best model
  res$best_params <- res$fit_grid %>% selection_crit(metric)
  set_mix <- if (res$misc$mixture == 'tune') { res$best_params$mixture } else { res$misc$mixture }
  
  # Parameter(s) change plot(s)
  base_metrics <- if (res$misc$mixture == 'tune') { list(
    penalty = collect_metrics(res$fit_grid) %>% filter(mixture == set_mix),
    mixture = collect_metrics(res$fit_grid) %>% filter(penalty == res$best_params$penalty)
  )} else { list(penalty = collect_metrics(res$fit_grid)) }
  
  res$penalty_plot <- ggplot(base_metrics$penalty, aes(penalty, mean, color = .metric)) +
    geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), alpha = 0.5) +
    geom_line(size = 1.5) +
    facet_wrap(~ .metric, scales = 'free', nrow = 2) +
    scale_x_log10() +
    theme(legend.position = 'none') +
    ggtitle(paste('penalty change at', set_mix, 'mixture'))
  
  if (res$misc$mixture == 'tune') {
    res$mixture_plot <- ggplot(base_metrics$mixture, aes(mixture, mean, color = .metric)) +
      geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), alpha = 0.5) +
      geom_line(size = 1.5) +
      facet_wrap(~ .metric, scales = 'free', nrow = 2) +
      scale_x_log10() +
      theme(legend.position = 'none') +
      ggtitle(paste('mixture change at', res$best_params$penalty, 'penalty'))
  }
  
  # Single stand-alone fit
  res$best_fit <- linear_reg(penalty = res$best_params$penalty, mixture = set_mix) %>% set_engine('glmnet', family = res$misc$family) %>%
    fit(res$misc$fmla, data = res$misc$df %>% select(any_of(res$misc$resp_and_preds)))
  
  # Variable importance
  
  # Temporary vi replacement due to this vip bug: https://github.com/koalaverse/vip/issues/103
  vi_glmnet <- function(model_fit, lambda) {
    cs <- coef(model_fit$fit, s = lambda)
    tibble(Variable = rownames(cs), Importance = abs(as.vector(cs)), Sign = if_else(as.vector(cs) >= 0, 'POS', 'NEG')) %>%
      arrange(desc(Importance))
  }
  
  res$best_fit_tibble <- res$best_fit %>% tidy() %>% filter(estimate != 0)
  
  res$importance_plot <- res$best_fit %>%
    vi_glmnet(lambda = res$best_params$penalty) %>% # Would just be vi without the above bug
    mutate(Variable = fct_reorder(Variable, Importance)) %>%
    head(min(nrow(res$best_fit_tibble) - 1, top_n_vars)) %>%
    ggplot() + geom_col(aes(x = Importance, y = Variable, fill = Sign)) +
      scale_x_continuous(expand = c(0, 0)) + scale_fill_manual(values = c('tomato', 'royalblue')) +
      labs(y = NULL)
  
  res
}

