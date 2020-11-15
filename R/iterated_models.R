


# Extract glm model variable names lower than a given p-value
getModTerms <- function(glmMod, threshold) {
  coeffs <- summary(glmMod)$coefficients
  as_tibble(coeffs) %>%
    add_column(varName = coeffs %>% rownames, .before = 1) %>%
    filter(varName != "(Intercept)") %>%
    filter(`Pr(>|z|)` < threshold) %>%
    pull(varName)
}



#### Backward steps ####


# GLM selction by removal of least siginificant term steps
# significance_statistic needs to match the last column of a glm summary; it is usually z, but for, say, the gaussian family it is t
# See anova.glm documentation for most appropriate statistical test (again, usually Chisq, but F for gaussian)
stepLeastSigGLM <- function(glmMod, threshold = 0.2, data = glmMod$data, significance_statistic = 'z', test = 'Chisq', test_threshold = threshold, nested = F) {
  discarded <- c()
  failed_discards <- c()
  mod <- glmMod
  maxIter <- summary(mod)$coefficients %>% length - 1
  
  p_val_col <-  paste0('Pr(>|', significance_statistic, '|)')
  
  for (i in 1:maxIter) {
    coeffs <- summary(mod)$coefficients
    terms <- as_tibble(coeffs) %>%
      add_column(varName = coeffs %>% rownames, .before = 1) %>%
      filter(varName != '(Intercept)') %>%
      arrange(desc(!!ensym(p_val_col)))
    
    for (r in 1:nrow(terms)) {
      if (terms[r,][[p_val_col]] < threshold) { break }
      if (r > 1) { failed_discards <- append(failed_discards, terms[r-1,]$varName) }
      new_mod <- update(mod, formula. = as.formula(paste('~ . -', terms[r,]$varName)), data = data)
      test_res <- anova(new_mod, mod, test = test)
      if (test_res[[paste0('Pr(>', test, ')')]][2] > test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, terms[r,]$varName)
        break
      }
    }
  }
  
  if (!nested) { # Deal with Intercept
    coeffs <- summary(mod)$coefficients
    terms <- as_tibble(coeffs) %>%
      add_column(varName = coeffs %>% rownames, .before = 1) %>%
      filter(varName == '(Intercept)')
    
    if (nrow(terms) != 0 & terms[1,][[p_val_col]] > threshold) {
      new_mod <- update(mod, formula. = ~ . + 0)
      test_res <- anova(new_mod, mod, test = test)
      if (test_res[[paste0('Pr(>', test, ')')]][2] < test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, '(Intercept)')
      }
    }
    redo <- stepLeastSigGLM(mod, threshold = threshold, data = data, significance_statistic = significance_statistic, test = test, test_threshold = test_threshold, nested = T)
    discarded <- c(discarded, redo$discarded)
    mod <- redo$mod
  }
  
  if (!nested & (length(failed_discards) > 0)) { print(paste('Failed discards by', test, 'test:', failed_discards)) }
  
  list(mod = mod, discarded = discarded)
}


# Same as above but for Proportional Odds Logistic Regression
#  NOTE: might want to give the input polr model a start argument; see start_vals function inside below
stepLeastSigPOLR <- function(polrMod, data, threshold = 0.2, test_threshold = threshold, nested = F) {
  discarded <- c()
  failed_discards <- c()
  mod <- polrMod
  maxIter <- tidy(mod) %>% filter(coef.type == 'coefficient') %>% nrow() - 1
  n_levels <- length(mod$lev)
  start_vals <- function(n_predictors) { c(rep.int(0, n_predictors), 0:(n_levels-2) - (n_levels-2)/2) }
  
  for (i in 1:maxIter) {
    terms <- tidy(mod) %>% filter(coef.type == 'coefficient') %>%
      rename(varName = term) %>%
      mutate(p_val = pt(-abs(statistic), nrow(mod$fitted.values) - (n() + (n_levels - 1)) - 1)) %>%
      arrange(desc(p_val)) # Above: degrees of freedom are #observations - #coefficients - 1, and #coefficients is #covariates + #intercepts, where #intercepts is the number of response levels - 1
    
    for (r in 1:nrow(terms)) {
      if (terms[r,]$p_val < threshold) { break }
      if (r > 1) { failed_discards <- append(failed_discards, terms[r-1,]$varName) }
      new_mod <- update(mod, formula. = as.formula(paste('~ . -', terms[r,]$varName)), data = data, start = start_vals(nrow(terms) - 1))
      test_res <- anova(new_mod, mod, test = 'Chisq')
      if (test_res$`Pr(Chi)`[2] > test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, terms[r,]$varName)
        break
      }
    }
  }
  
  if (!nested & (length(failed_discards) > 0)) { print(paste('Failed discards by', test, 'test:', failed_discards)) }
  
  list(mod = mod, discarded = discarded)
}


# GAM selction by removal of least siginificant term steps
stepLeastSigGAM <- function(gamMod, threshold = 0.2, data = gamMod$data, nested = F) {
  discarded <- c()
  mod <- gamMod
  maxIter <- summary(mod)$p.pv %>% length - 1
  for (i in 1:maxIter) {
    terms <- tibble(varName = summary(mod)$p.pv %>% names, pVals = summary(mod)$p.pv) %>%
      filter(varName != "(Intercept)") %>%
      arrange(desc(pVals))
    
    if (terms[1,]$pVals < threshold) { break } else {
      discarded <- append(discarded, terms[1,]$varName)
      mod <- update(mod, formula. = as.formula(paste("~ . -", terms[1,]$varName)))
    }
  }
  
  if (!nested) { # Deal with Intercept
    terms <- tibble(varName = summary(mod)$p.pv %>% names, pVals = summary(mod)$p.pv) %>%
      filter(varName == "(Intercept)")
    
    if (nrow(terms) != 0 & terms[1,]$pVals > threshold) {
      discarded <- append(discarded, "(Intercept)")
      mod <- update(mod, formula. = ~ . + 0)
    }
    redo <- stepLeastSigGAM(mod, threshold = 0.1, data = data, nested = T)
    mod <- redo$mod
    discarded <- c(discarded, redo$discarded)
  }
  
  list(mod = mod, discarded = discarded)
}

# GAM selction by removal of least siginificant term steps for GAMs with 2 formulae (e.g. ziplss)
stepLeastSigGAM2 <- function(gamMod, threshold = 0.2, modUpdater = NA, nested = F) {
  discarded <- c()
  mod <- gamMod
  maxIter <- summary(mod)$p.pv %>% length - 1
  for (i in 1:maxIter) {
    terms <- tibble(varName = summary(mod)$p.pv %>% names, pVals = summary(mod)$p.pv) %>%
      filter(!str_starts(varName, "\\(Intercept\\)")) %>%
      arrange(desc(pVals))
    
    if (terms[1,]$pVals < threshold) { break } else {
      discarded <- append(discarded, terms[1,]$varName)
      if (str_ends(terms[1,]$varName, ".1")) {
        form <- update.formula(summary(mod)$formula[[2]], paste("~ . -", str_remove(terms[1,]$varName, "\\.1")))
        mod <- modUpdater(summary(mod)$formula[[1]], form)
      } else {
        form <- update.formula(summary(mod)$formula[[1]], paste(". ~ . -", terms[1,]$varName))
        mod <- modUpdater(form, summary(mod)$formula[[2]])
      }
    }
  }
  
  if (!nested) { # Deal with Intercept
    for (i in 1:2) {
      terms <- tibble(varName = summary(mod)$p.pv %>% names, pVals = summary(mod)$p.pv) %>%
        filter(str_starts(varName, "\\(Intercept\\)")) %>%
        arrange(desc(pVals))
      
      if (nrow(terms) != 0 & terms[1,]$pVals < threshold) {
        redo <- stepLeastSigGAM2(mod, threshold = 0.1, modUpdater = modUpdater, nested = T)
        mod <- redo$mod
        discarded <- c(discarded, redo$discarded)
      } else {
        discarded <- append(discarded, terms[1,]$varName)
        if (str_ends(terms[1,]$varName, ".1")) {
          form <- update.formula(summary(mod)$formula[[2]], ~ . + 0)
          mod <- modUpdater(summary(mod)$formula[[1]], form)
        } else {
          form <- update.formula(summary(mod)$formula[[1]], . ~ . + 0)
          mod <- modUpdater(form, summary(mod)$formula[[2]])
        }
        redo <- stepLeastSigGAM2(mod, threshold = c(threshold,0.1)[i], modUpdater = modUpdater, nested = T)
        mod <- redo$mod
        discarded <- c(discarded, redo$discarded)
      }
    }
  }
  
  list(mod = mod, discarded = discarded)
}



#### Forward steps ####

# GLM selction by stepwise variable insertion followed by possible removal of least siginificant term if below threshold
# significance_statistic needs to match the last column of a glm summary; it is usually z, but for, say, the gaussian family it is t
# See anova.glm documentation for most appropriate statistical test (again, usually Chisq, but F for gaussian)
#   NOTE: check that the formula in the starting model does not contain any of the variables to add
stepLeastSigGLM_forward <- function(glmMod, variables, threshold = 0.2, data = glmMod$data, test_forward = T, significance_statistic = 'z', test = 'Chisq', test_threshold = threshold, nested = F) {
  discarded <- c()
  failed_discards <- c()
  mod <- glmMod
  
  p_val_col <-  paste0('Pr(>|', significance_statistic, '|)')
  
  for (v in variables) {
    new_mod <- update(mod, formula. = as.formula(paste('~ . +', v)), data = data)
    coeffs <- summary(new_mod)$coefficients
    terms <- as_tibble(coeffs) %>%
      add_column(varName = coeffs %>% rownames, .before = 1) %>%
      filter(varName != '(Intercept)') %>%
      arrange(desc(!!ensym(p_val_col)))
    
    if (terms[1,][[p_val_col]] < threshold) {
      if (test_forward) {
        test_res <- anova(mod, new_mod, test = test)
        if (test_res[[paste0('Pr(>', test, ')')]][2] < test_threshold) { mod <- new_mod }
        else { discarded <- append(discarded, terms[1,]$varName) }
      } else { mod <- new_mod }
    } else if (terms[1,]$varName == v) {
      discarded <- append(discarded, terms[1,]$varName)
    } else { # Perform iterated BACKWARD steps if the first one is not the last added variable
      redo <- stepLeastSigGLM(new_mod, threshold = threshold, data = data, significance_statistic = significance_statistic, test = test, test_threshold = test_threshold, nested = T)
      discarded <- c(discarded, redo$discarded)
      mod <- redo$mod
    }
  }
  
  if (!nested) { # Deal with Intercept
    coeffs <- summary(mod)$coefficients
    terms <- as_tibble(coeffs) %>%
      add_column(varName = coeffs %>% rownames, .before = 1) %>%
      filter(varName == '(Intercept)')
    
    if (nrow(terms) != 0 & terms[1,][[p_val_col]] > threshold) {
      new_mod <- update(mod, formula. = ~ . + 0)
      test_res <- anova(new_mod, mod, test = test)
      if (test_res[[paste0('Pr(>', test, ')')]][2] < test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, '(Intercept)')
      }
    }
    redo <- stepLeastSigGLM(mod, threshold = threshold, data = data, significance_statistic = significance_statistic, test = test, test_threshold = test_threshold, nested = T)
    discarded <- c(discarded, redo$discarded)
    mod <- redo$mod
  }
  
  list(mod = mod, discarded = unique(discarded))
}


# Run multiple and select the lowest AIC run of stepLeastSigGLM_forward
best_steps <- function(glmMod, variables, seeds = c(1,2,7,17,42,73,100,117,1000,1000000), score_f = AIC, threshold = 0.2, data = glmMod$data, test_forward = T, significance_statistic = 'z', test = 'Chisq', test_threshold = threshold, nested = F) {
  best <- NULL
  scores <- list()
  for (seed in seeds) {
    set.seed(seed)
    vs_to_use <- sample(variables)
    done <- stepLeastSigGLM_forward(glmMod, variables, threshold = threshold, data = data, test_forward = test_forward, significance_statistic = significance_statistic, test = test, test_threshold = test_threshold)
    scores[[as.character(seed)]] <- score_f(done$mod)
    if (is.null(best) || (score_f(done$mod) < score_f(best$mod))) { best <- done }
  }
  list(best = best, scores = scores)
}






##########################################################
## TO BE GENERALISED
##########################################################

# Score a model based on prediction correctness
scoreModel <- function(mod, name, responseThreshold = 0.5, data = mod$data) {
  predsDf <- data %>%
    select(t, storm) %>%
    mutate(pred = predict(mod, newdata = data, type = "response")) %>%
    mutate(prediction = if_else(pred > responseThreshold, T, F)) %>%
    filter(!is.na(pred)) %>%
    mutate(correct = storm == prediction)
  
  metricsDf <- tibble(
    model = name,
    Trate = accuracy(predsDf$storm, predsDf$prediction),
    TPrate = recall(predsDf$storm, predsDf$prediction),
    precision = precision(predsDf$storm, predsDf$prediction),
    f1 = 2 * TPrate * precision / (TPrate + precision),
    naCount = nrow(data) - nrow(predsDf),
    naCountRate = (nrow(data) - nrow(predsDf)) / nrow(data),
    missedStorms = sum(fieldForStorms$storm) - sum(predsDf$storm),
    missedStormsRate = (sum(data$storm) - sum(predsDf$storm)) / sum(data$storm))
  
  list(predsDf = predsDf, metricsDf = metricsDf)
}