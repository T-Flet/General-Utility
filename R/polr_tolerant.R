library(tidyverse)



# Version of polr which manually sets default starting values if no good ones are generated
  # NOTE: if the model produced by this function has a start value (i.e. if it failed without) then
  #         applying broom::tidy with p.vals = T to it will fail since this start will be one element
  #         too long because tidy.polr uses MASS:dropterm, which compares it with all models with 1
  #         fewer variables. The solution is to manually amend the saved call before dropterm:
  #           MODEL$call$start <- MODEL$call$start[2:length(MODEL$call$start)]
  #         This is handled by tidy_polr_tolerant
polr_tolerant <- function(fmla, data, method = 'logistic') {
  data_name <- match.call()$data
  
  tryCatch( # Fit a polr model handling fit failure gracefully
    tryCatch( # Try solving the issue by setting bland starting values
      {
        mod <- polr(fmla, data, Hess = T, method = method)
        mod$call <- call('polr', formula = fmla, data = data_name, Hess = T, method = method)
        mod # ^^ R scoping shenanigans because some call fields are not retrievable outside this function
      },
      error = function(e) {
        response <- all.vars(fmla)[attr(terms(fmla), which = 'response')]
        terms <- attr(terms(fmla), which = 'term.labels')
        n_levels <- length(levels(data[[response]]))
        start_vals <- c(rep.int(0, length(terms)), 1:(n_levels - 1))
        mod <- polr(fmla, data, Hess = T, method = method, start = start_vals)
        mod$call <- call('polr', formula = fmla, data = data_name, Hess = T, method = method, start = start_vals)
        mod # ^^ R scoping shenanigans because some call fields are not retrievable outside this function
      }
    ),
    error = function(e) { warning(paste('A polr model could not be fit:\n', e), call. = F) }
  )
}


library(lmtest)
tidy_polr_tolerant <- function(pt_mod, dropdown_p_values = F) {
  if (dropdown_p_values) { # These are the (expensive and limited) p.values recommended by polr's author
    pt_mod$call$start <- pt_mod$call$start[2:length(pt_mod$call$start)] # Remove one pre-set starting value since dropdown tries to use them
    res <- tidy(pt_mod, p.values = T)
    pt_mod$call$start <- c(0, pt_mod$call$start)
  } else { # These are typical p.values
    cts <- coeftest(pt_mod)
    res <- tidy(pt_mod, p.values = F)
    res <- res %>% mutate(p.value = c(cts[,ncol(cts)], rep(NA, nrow(res) - nrow(cts))))
  }
  res
}


update_polr_tolerant <- function(pt_mod, term, data = pt_mod$model, add = F) {
  # Add or remove one pre-set starting value since dropdown tries to use them
  pt_mod$call$start <- if (add) c(0, pt_mod$call$start) else pt_mod$call$start[2:length(pt_mod$call$start)]
  res <- update(pt_mod, formula. = as.formula(paste('~ .', ifelse(add, '+', '-'), term)), data = pt_mod$model)
  pt_mod$call$start <- if (add) pt_mod$call$start[2:length(pt_mod$call$start)] else c(0, pt_mod$call$start)
  res
}


