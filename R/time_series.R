
library(zoo)
library(forecast)
library(stinepack)

library(smooth)
library(mgcv)

library(broom)
library(tidyverse)

# Requires the match_approx function from generic_functions.R for the fill_if_regular function


t_x_to_ts <- function(t_or_ts, x = NULL) if (is.null(x)) zoo(coredata(t_or_ts), order.by = time(t_or_ts)) else zoo(x, order.by = t_or_ts)
ts_to_df <- function(xs) tibble(t = time(xs), x = coredata(xs))
set_frequency <- function(xs, freq) zoo(coredata(xs), order.by = time(xs), frequency = freq)


# Minor modification of zoo::is.regular.zoo to return the regular interval (NA otherwise) instead of a boolean
get_regular <- function(xs, strict = F, tolerance = 0.001) {
  xs <- as.zoo(xs)
  delta <- suppressWarnings(try(diff(as.numeric(index(xs))), silent = T))
  if (inherits(delta, 'try-error') || anyNA(delta) || length(delta) < 1) NA
  else if (strict) if (isTRUE(all.equal(delta, rep.int(delta[1], length(delta)), tolerance = tolerance))) delta[1] else NA
  else {
    delta <- unique(delta)
    multiples <- isTRUE(all.equal(delta/min(delta), round(delta/min(delta)), tolerance = tolerance))
    if (multiples || isTRUE(all.equal(delta, round(delta), tolerance = tolerance))) min(delta) else NA
  }
}


# Fill a regular series (interval provided or identified) with entries (Stineman interpolations or NAs) for the missing times
fill_if_regular <- function(xs, reg_interval = NULL, interpolate = T) {
  if (is.null(reg_interval)) reg_interval <- get_regular(xs)
  if (is.na(reg_interval)) {
    print('Series was not regular; returning original')
    xs
  } else {
    good_times <- seq(time(xs)[1], time(xs)[length(xs)], by = reg_interval)
    good_is <- match_approx(good_times, time(xs)) # vv the if_else is necessary because good_times are only approximately equal to the real times
    res <- zoo(coredata(xs)[good_is], order.by = if_else(is.na(time(xs)[good_is]), good_times, time(xs)[good_is]))
    if (interpolate) na.stinterp(res) else res
  }
}


# Minor modification of forecast::findfrequency
#   Allows arbitrary detrending (intent: be slightly more than linear but slightly less than required to identify a period
#   Allows toggling whether to only use contiguous intervals from the given series to compute the period
#   Avoids periods greater than half the dataset
dominant_frequency <- function(xs, spline_k = 3, no_na_interval_only = F, avoid_longer_than_half = T, plot_trend = F) {
  mod <- gam(x ~ s(t, k = spline_k, bs = 'tp'), data = ts_to_df(coredata(xs)), method = 'REML')
  xs <- as.ts(residuals(mod))

  n.freq <- 500
  period <- 1L
  spec <- spec.ar(c(if (no_na_interval_only) na.contiguous(xs) else xs), plot = F, n.freq = n.freq)
  if (max(spec$spec) > 10) { # Arbitrary threshold chosen by trial and error.
    candidate <- floor(1 / spec$freq[which.max(spec$spec)] + 0.5)
    if (candidate == Inf || (avoid_longer_than_half && candidate > length(xs) / 2)) { # Find next local maximum
      j <- which(diff(spec$spec) > 0)
      if (length(j) > 0) {
        nextmax <- j[1] + which.max(spec$spec[(j[1] + 1):n.freq])
        if (nextmax < length(spec$freq)) period <- floor(1 / spec$freq[nextmax] + 0.5)
      }
    } else period <- candidate
  }

  if (plot_trend) plot(mod, residuals = T)

  if (avoid_longer_than_half && period > length(xs) / 2) 1L else as.integer(period)
}


# Lightly detrend and compute acf or pacf for each of the given variables
#   NOTE: Factor variables are not detrended
acf_pacf <- function(df_or_list, lag.max = 10, detrend = T, spline_k = 3, partial = F) {
  map(df_or_list, function(xs) {
    if (is.numeric(xs)) {
      mod <- gam(x ~ s(t, k = spline_k, bs = 'tp'), data = tibble(t = time(xs), x = xs), method = 'REML')
      xs <- residuals(mod)
    }
    f <- if (partial) pacf else acf
    as.vector(f(xs, lag.max = lag.max, plot = F)$acf)
  }) %>% bind_rows()
}
# acf_pacf(tibble(iris), partial = T)


# Further automation of smooth::auto.adam, with optional frequency detection
#   The freq argument should be either NULL (for auto-detection) or a positive integer
#   Use ... to pass arguments to auto.adam; already passed-in ones are model (ets_type), distribution, orders, occurrence and initial
auto_adam <- function(xs, freq = NULL, type = c('both', 'ets', 'arima'), max_arima = c(2, 0, 2),
                      ets_type = NULL,
                      distribution = c('dnorm', 'dlaplace', 'ds', 'dgnorm', 'dlnorm', 'dinvgauss', 'dgamma'),
                      initial = c('optimal', 'backcasting'),
                      check_double_freq = T,
                      ...) {
  type <- match.arg(type)
  initial <- match.arg(initial)

  if (is.null(freq)) {
    freq <- dominant_frequency(xs)
    if (freq == 1) print('No frequency detected')
    else print(paste('Detected dominant frequency of', freq, 'observations'))
  } else if (freq != round(freq) || freq <= 0) stop('The freq argument should be either NULL (for auto-detection) or a positive integer')
  else freq <- round(freq)
  xs <- set_frequency(xs, freq)

  orders <- if (type == 'ets') list(ar = 0, i = 0, ma = 0) else {
    if (freq != 1) max_arima <- map(max_arima, ~ c(.x, .x))
    list(ar = max_arima[[1]], i = max_arima[[2]], ma = max_arima[[3]], select = T)
  }

  call_mod <- function(ets_str, used_freq = freq) auto.adam(data = xs, model = ets_str, lags = c(used_freq), distribution = distribution, orders = orders, occurrence = 'auto', initial = initial, ...)
  pick_mod <- function(ets_str) {
    base_freq_mod <- call_mod(ets_str)
    if (check_double_freq && freq != 1 && 2 * freq < length(xs) / 2) {
      tryCatch({ # Try-ing because of https://github.com/config-i1/smooth/issues/199
        double_freq_mod <- call_mod(ets_str, 2 * freq)
        if (BICc(double_freq_mod) < BICc(base_freq_mod)) double_freq_mod else base_freq_mod
      }, error = function(e) base_freq_mod)
    } else base_freq_mod
  }

  if (!is.null(ets_type)) mod <- pick_mod(ets_type)
  else if (type == 'arima') mod <- pick_mod('NNN')
  else if (type == 'both' && frequency(xs) != 1) { # This is because adam has difficulty with seasonal arima and because Z-containing types may not be put in a vector of models to test
    Amod <- pick_mod('ZZA')
    Mmod <- pick_mod('ZZM')
    IC_comparison <- summary(Amod)$ICs < summary(Mmod)$ICs # Might as well go by majority rather than pick one
    mod <- if (sum(IC_comparison) >= 2) Amod else Mmod
  } else mod <- pick_mod('ZZZ')

  mod
}
# AA <- auto_adam(AirPassengers)


# tidy-like output for adam models
tidy_adam <- function(mod) {
  smry <- summary(mod)
  if (is.null(smry$coefficients)) tibble(term = character(), estimate = numeric(), std.error = numeric(), lower = numeric(), upper = numeric(), significance = logical())
  else as_tibble(smry$coefficients) %>%
      rename(estimate = Estimate, std.error = `Std. Error`, lower = `Lower 2.5%`, upper = `Upper 97.5%`) %>%
      mutate(term = rownames(smry$coefficients), .before = 1) %>%
      mutate(significance = smry$significance)
}


# Dataframe of actuals, fitted values, residuals and all states (e.g. level, trend, seasonal, ARIMA components, ...)
#   NOTE: cut_non_actual_times removes the initial rows for which there is no true data (states columns may have values, but the others have NAs)
adam_all_states <- function(mod, cut_non_actual_times = T) {
  data <- cbind(actuals(mod), fitted(mod), mod$states, residuals(mod))
  colnames(data) <- c('actuals', 'fitted', colnames(mod$states), 'residuals')
  if (cut_non_actual_times) data <- data[which(!is.na(data$actuals))[1]:nrow(data)]
  as_tibble(data) %>% mutate(Time = time(data), .before = 1)
}


# Standardised-output (as in with the others in this package) prediction function for adam models
# The last 3 arguments trigger additional outputs beside the raw prediction (and interval):
#   bounds for a given probability (to be interpreted as central or not)
#   or probabilities of values being above or below a given bound.
# Specifically:
#   If a custom_prob is given, then the Custom_Lower and Custom_Upper columns are added
#     - if prob_as_central, then custom_prob is centred and the two columns are the bounds of the corresponding CI
#     - if not, then the two columns are each a bound for custom_prob of the values being above or below them
#   If a custom_bound is given then the Prob_Below_Bound and Prob_Above_Bound columns are added
# Since the returned probabilities for a given custom_bound input are unfortunately approximated by
#   sampling the density at each prediction, said quantile probabilities can be provided manually
#   to increase or decrease the approximation; they should only be values above 0.5 (NOT INCLUDED)
predict_adam <- function(mod, new_data, interval = NULL, only_new_data = F, draw_plot = F,
                         custom_prob = NULL, prob_as_central = F, custom_bound = NULL,
                         # vv Could replace 0.1 with 0.075 for one more pair of values
                         quantiles_probs_for_approximation = c(seq(0.55, 0.9, 0.1), 0.99)) {
  diff_f <- function(a, b) if (is.numeric(a) && is.numeric(b)) a - b else difftime(a, b)
  horizon <- ceiling(as.numeric(diff_f(max(new_data), max(time(mod$fitted)))) / as.numeric(diff_f(time(mod$fitted)[2], time(mod$fitted)[1])))

  if (horizon < 1) {
    print('The requested prediction times are in the series past; returning the modelled values covering the requested interval')
    refit <- reapply(mod)
    confidence <- apply(refit$refitted, 1, function(xs) { c(quantile(xs, 0.025, na.rm = T), quantile(xs, 0.975, na.rm = T)) })
    ok_index <- (index(refit$fitted) >= min(new_data)) & (index(refit$fitted) <= max(new_data))
    res <- refit$fitted[ok_index]
    tibble(Time = time(res), Pred = as.vector(res), Upper = confidence[2, ok_index], Lower = confidence[1, ok_index])
  } else {
    if (is.null(interval)) {
      ets_type <- str_match(mod$model, 'ETS\\((\\w+)\\)') # vv recommended in adam vignette
      interval <- if (!is.na(ets_type[2]) && ets_type[2][1] == 'M') 'simulated' else 'prediction'
    }

    preds <- forecast(mod, h = horizon, interval = interval)
    if (draw_plot) plot(preds)

    res <- tibble(Time = time(preds$mean), Pred = as.vector(preds$mean), Lower = as.vector(preds$lower), Upper = as.vector(preds$upper))

    if (!is.null(custom_prob)) {
      central_prob <- if (prob_as_central) custom_prob else (2 * custom_prob - 1) # I.e. 1 - (1 - custom_prob) * 2 if it is not meant to be central
      custom_preds <- forecast(mod, h = horizon, interval = interval, level = central_prob)
      res <- res %>% mutate(Custom_Lower = as.vector(custom_preds$lower), Custom_Upper = as.vector(custom_preds$upper))
    }
    
    if (!is.null(custom_bound)) {
      ## This is a hard task without replacing the forecast.adam implementation to accommodate this output.
      ##  (Or without compying every if statement in it to then compute ad-hoc solutions)
  
      ## Brute-force optimisation approach: very expensive even for just the last prediction
      # bound_is_upper <- custom_bound >= preds$mean[horizon]
      # bound_exceeds <- if (bound_is_upper) custom_bound > preds$upper[horizon] else custom_bound < preds$lower[horizon]
      # to_minimise <- if (bound_is_upper) {
      #        function(par) { forecast(mod, h = horizon, interval = interval, level = par[1])$upper[horizon] - custom_bound } }
      # else { function(par) { forecast(mod, h = horizon, interval = interval, level = par[1])$lower[horizon] - custom_bound } }
      # prob <- optim(par = c(0.95), fn = to_minimise, method = 'Brent',
      #       lower = if (bound_exceeds) 0.95 else 0.01,
      #       upper = if (bound_exceeds) 1 else 0.95,
      #       control = list(maxit = 10))$par
      # prob <- (1 + prob) / 2 # To non-central probability
      # res <- res %>%
      #   mutate(Prob_Below_Bound = c(rep(NA, horizon - 1), prob), Prob_Above_Bound = c(rep(NA, horizon - 1), 1 - prob)) %>%
      #   mutate(Prob_Below_Bound = if_else(custom_bound >= Pred, Prob_Below_Bound, 1 - Prob_Below_Bound),
      #          Prob_Above_Bound = if_else(custom_bound >= Pred, Prob_Above_Bound, 1 - Prob_Above_Bound))
      
      ## Assume Normality approach: very inaccurate, even splitting upper and lower sides as separate Normals
      ### STRONG APPROXIMATION ASSUMPTION: approximate either side of the error as a normal, ###
      ###   i.e. preds$upper/lower == preds$mean +/- qnorm(0.975) * SD                       ###
      # res <- res %>%
      #   mutate(TempSD = ((Upper - Lower) / 2) / qnorm(0.975)) %>%
      #   mutate(Prob_Below_Bound = pnorm(custom_bound, Pred, TempSD, lower.tail = custom_bound >= Pred)) %>%
      #   mutate(Prob_Above_Bound = 1 - Prob_Below_Bound) %>%
      #   select(-TempSD) %>%
      #   mutate(Prob_Below_Bound = if_else(custom_bound >= Pred, Prob_Below_Bound, 1 - Prob_Below_Bound),
      #          Prob_Above_Bound = if_else(custom_bound >= Pred, Prob_Above_Bound, 1 - Prob_Above_Bound))
      
      ## Density approximation approach (sample a few points and then estimate linearly between closest ones) acceptable
      probs <- quantiles_probs_for_approximation
      if (!is.null(custom_prob)) probs <- probs[probs != custom_prob]
      fss <- probs %>% as.list() %>% map(function(prob) {
        new_preds <- forecast(mod, h = horizon, interval = interval, level = 2 * prob - 1) # Do not want central probs
        named <- list()
        named[[as.character(1 - prob)]] <- as.vector(new_preds$lower) ; named[[as.character(prob)]] <- as.vector(new_preds$upper)
        named
      }) %>% flatten()
      fss[['0.05']] <- res$Lower ; fss[['0.95']] <- res$Upper ; probs <- c(probs, 0.05, 0.95)
      if (!is.null(custom_prob)) { fss[[as.character(1 - custom_prob)]] <- res$Custom_Lower ; fss[[as.character(custom_prob)]] <- res$Custom_Upper }
      
      fss <- do.call('cbind', fss[order(as.numeric(names(fss)))]) # Matrix with ordered columns
      out_probs <- map(1:nrow(fss), function(i) {
        more_than_this_prob <- length(which(fss[i,] < custom_bound))
        if (more_than_this_prob == 0) as.numeric(colnames(fss)[1])
        else if (more_than_this_prob == ncol(fss)) as.numeric(colnames(fss)[ncol(fss)])
        else (as.numeric(colnames(fss)[more_than_this_prob + 1]) + as.numeric(colnames(fss)[more_than_this_prob])) / 2
      }) %>% unlist()

      res <- res %>% mutate(Prob_Below_Bound = out_probs, Prob_Above_Bounds = 1 - out_probs)
    }
    
    if (only_new_data) {
      temp_res <- res %>% filter(Time %in% new_data)
      if (nrow(temp_res) != nrow(new_data)) {
        print('The requested prediction times do not line up with the model time unit; returning lined-up predictions covering the requested interval')
        res %>% filter((Time >= min(new_data)) & (Time <= max(new_data)))
      } else temp_res
    } else res
  }
}
# AA <- auto_adam(AirPassengers)#, distribution = c('dlaplace'))
# BB <- predict_adam(AA, yearmon(c(1965.25, 1965.5)),
#                    custom_prob = 0.999, prob_as_central = F, custom_bound = 550)



## Missing values functions

# AA <- read_csv('https://r-data.pmagunia.com/system/files/datasets/dataset-58057.csv')
# AA <- AA[! index(AA) %in% sample.int(nrow(AA), nrow(AA) %/% 10),]
# AA <- t_x_to_ts(AA$time, AA$AirPassengers)
# is.regular(AA)
# get_regular(AA)
# filled_AA <- fill_if_regular(AA)
# length(filled_AA) > length(AA)
# all(filled_AA == AA) # True because comparing only on same indices
# plot(AA)
# plot(AirPassengers)
# plot(filled_AA)


## Model-related functions

# # AA <- t_x_to_ts(USAccDeaths)
# # AA <- t_x_to_ts(lynx)
# AA <- t_x_to_ts(AirPassengers)
#
# ZZ <- auto_adam(AA, distribution = 'dgamma', type = 'arima')#, initial = 'backcasting')
# summary(ZZ)
# tidy_adam(ZZ)
#
# adam_all_states(ZZ, F) %>% head(15)
# adam_all_states(ZZ, T) %>% head(15)
#
# par(mfcol = c(3, 4))
# plot(ZZ, c(1:11))
# plot(ZZ, which = 12)
# par(mfcol = c(1, 1))
#
#
# # # USAccDeaths
# # predict_adam(ZZ, yearmon(c(1975.25, 1978.25))) # apply yearmon to new_data for USAccDeaths and AirPassnegers
# # predict_adam(ZZ, yearmon(1980), draw_plot = T)
# # # lynx
# # predict_adam(ZZ, c(1950.25, 1955.25)) # apply yearmon to new_data for USAccDeaths and AirPassnegers
# # predict_adam(ZZ, 1970, draw_plot = T)
# # AirPassengers
# predict_adam(ZZ, yearmon(c(1950.25, 1955.25))) # apply yearmon to new_data for USAccDeaths and AirPassnegers
# predict_adam(ZZ, yearmon(1970), draw_plot = T)


