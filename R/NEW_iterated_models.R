library(MASS)
library(nnet)
library(broom)
library(tidyverse)



#### Backward steps ####


# GLM selction by removal of least siginificant term steps
# significance_statistic needs to match the last column of a glm summary; it is usually z, but for, say, the gaussian family it is t
# See anova.glm documentation for most appropriate statistical test (again, usually Chisq, but F for gaussian)
#   NOTES:
#     - If there are factor covariates no other covariate name may start with their names
#     - The nested argument is relevant for forwards refining (which has nested backwards); all it does is make the intercept non-removable
#     - The term_deps_generated_from_model argument has such verbose name to ensure that a raw dependencies object
#         (as in, generated from gen_dependent_terms and not from gen_dependent_terms_from_model) is NOT passed in
refine_glm_backwards <- function(gmod, threshold = 0.1, test_threshold = threshold,
                                 data = gmod$data, term_deps_generated_from_model = NULL, nested = F) {
  test <- if_else(gmod$family$family %in% c('binomial', 'poisson'), 'Chisq', 'F')
  
  term_details <- get_terms_details(gmod)
  identified <- setNames(map(term_details, ~ .x$identified), map(term_details, ~ .x$original)) # original -> identified
  deps <- if (is.null(term_deps_generated_from_model)) gen_dependent_terms_from_model(mod, term_details) else term_deps_generated_from_model
  
  discarded <- c()
  failed_discards <- c()
  mod <- gmod
  changes_made <- T
  while (changes_made) {
    changes_made <- F
    removable <- removable_terms(deps, discarded)
    removable <- tidy(mod) %>%
      filter(term %in% names(identified)) %>% mutate(id = identified[term]) %>% filter(id %in% removable) %>%
      group_by(id) %>% arrange(p.value) %>% slice_head() %>% ungroup() %>% # Keep most significant factor level
      arrange(desc(p.value))
    
    if (nrow(removable) == 0) break
    for (r in 1:nrow(removable)) {
      if (removable[r,]$p.value < threshold) break # I.e. all remaining are significant
      if (r > 1) failed_discards <- append(failed_discards, removable[r-1,]$id)
      new_mod <- update(mod, formula. = as.formula(paste('~ . -', removable[r,]$id)), data = data)
      test_res <- anova(new_mod, mod, test = test)
      if (test_res[2, ncol(test_res)] > test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, removable[r,]$id)
        changes_made <- T
        break
      }
    }
  }
  
  if (!nested) { # Deal with Intercept
    removable <- tidy(mod) %>% filter(term == '(Intercept)')

    if (nrow(removable) != 0 & removable[1,]$p.value > threshold) {
      new_mod <- update(mod, formula. = ~ . + 0)
      test_res <- anova(new_mod, mod, test = test)
      if (test_res[2, ncol(test_res)] > test_threshold) {
        mod <- new_mod
        discarded <- append(discarded, '(Intercept)')
      }
    }
    redo <- refine_glm_backwards(mod, threshold = threshold, test_threshold = test_threshold, data = data, term_deps_generated_from_model = deps, nested = T)
    discarded <- append(discarded, redo$discarded)
    mod <- redo$mod
  }
  
  if (!nested & (length(failed_discards) > 0)) { print(paste('Failed discards by', test, 'test:', failed_discards)) }
  
  list(mod = mod, discarded = discarded %>% unname() %>% unlist())
}


## Example from already-existing formula
mod <- glm(Petal.Length ~ Species + Sepal.Length + Sepal.Length:Species + I(Sepal.Length^2) + Sepal.Width + Petal.Width + Sepal.Width:I(Petal.Width^2), iris, family = Gamma())
summary(mod)
ref <- refine_glm_backwards(mod)
summary(ref$mod)
ref$discarded

## Example from a full generated formula and raw dependencies caveat
raw_deps <- gen_dependent_terms(c('Sepal.Length', 'Sepal.Width', 'Petal.Width'), c('Species'), max_power = 3, interactions = T)
raw_deps
mod_full <- glm(term_deps_to_formula('Petal.Length', raw_deps, removed = c()), iris, family = Gamma())
summary(mod_full)
ref_WRONG <- refine_glm_backwards(mod_full, term_deps_generated_from_model = raw_deps) #### DO NOT PASS IN RAW DEPENDENCIES
ref_full <- refine_glm_backwards(mod_full)
summary(ref_WRONG$mod) # Much worse than ref_full
summary(ref_full$mod)
ref_WRONG$discarded # Many fewer than ref_full
ref_full$discarded

# mod2 <- polr_tolerant(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, iris)
# tidy_polr_tolerant(mod2)
# mod3 <- multinom(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, iris)
# tidy(mod3)





