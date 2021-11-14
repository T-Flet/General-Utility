
library(zoo)
library(forecast)

library(smooth)
library(mgcv)

library(broom)
library(tidyverse)


t_x_to_ts <- function(t_or_ts, x = NULL) if (is.null(x)) zoo(coredata(t_or_ts), order.by = time(t_or_ts)) else zoo(x, order.by = t_or_ts)
ts_to_df <- function(xs) tibble(t = time(xs), x = coredata(xs))
set_frequency <- function(xs, freq) zoo(coredata(xs), order.by = time(xs), frequency = freq)


# Minor modification of forecast::findfrequency
#   Allows arbitrary detrending (intent: be slightly more than linear but slightly less than required to identify a period
#   Allows toggling whether to only use contiguous intervals from the given series to compute the period
dominant_frequency <- function(xs, spline_k = 3, no_na_interval_only = F, plot_trend = F) {
  mod <- gam(x ~ s(t, k = spline_k, bs = 'tp'), data = ts_to_df(coredata(xs)), method = 'REML')
  xs <- as.ts(residuals(mod))
  
  n.freq <- 500
  period <- 1L
  spec <- spec.ar(c(if (no_na_interval_only) na.contiguous(xs) else xs), plot = F, n.freq = n.freq)
  if (max(spec$spec) > 10) { # Arbitrary threshold chosen by trial and error.
    candidate <- floor(1 / spec$freq[which.max(spec$spec)] + 0.5)
    if (candidate == Inf) { # Find next local maximum
      j <- which(diff(spec$spec) > 0)
      if (length(j) > 0) {
        nextmax <- j[1] + which.max(spec$spec[(j[1] + 1):n.freq])
        if (nextmax < length(spec$freq)) period <- floor(1 / spec$freq[nextmax] + 0.5)
      }
    } else period <- candidate
  }
  
  if (plot_trend) plot(mod, residuals = T)
  
  as.integer(period)
}


# Further automation of smooth::auto.adam, with optional frequency detection
#   The freq argument should be either NULL (for auto-detection) or a positive integer
#   Use ... to pass arguments to auto.adam; already passed-in ones are model (ets_type), distribution, orders, occurrence and initial
auto_adam <- function(xs, freq = NULL, type = c('both', 'ets', 'arima'), max_arima = c(2, 0, 2),
                      ets_type = NULL,
                      distribution = c('dnorm', 'dlaplace', 'ds', 'dgnorm', 'dlnorm', 'dinvgauss', 'dgamma'),
                      initial = c('optimal', 'backcasting'), ...) {
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
  
  # The following two lines are here instead of a simple auto.adam-callling function until this issue is resolved: https://github.com/config-i1/smooth/issues/188
  call_mod <- function(ets_str) do.call(if (length(distribution) == 1) adam else auto.adam, c(list(data = xs, model = ets_str, distribution = distribution, orders = orders, occurrence = 'auto', initial = initial), list(...)))
  if (type == 'arima' && (length(distribution) != 1 || distribution != 'dnorm')) stop('At this time ARIMA-only models can only be fit with the "dnorm" distribution argument; see https://github.com/config-i1/smooth/issues/188')
  
  if (!is.null(ets_type)) mod <- call_mod(ets_type)
  else if (type == 'arima') mod <- call_mod('NNN')
  else if (type == 'both' && frequency(xs) != 1) { # This is because adam has difficulty with seasonal arima and because Z-containing types may not be put in a vector of models to test 
    Amod <- call_mod('ZZA')
    Mmod <- call_mod('ZZM')
    IC_comparison <- summary(Amod)$ICs < summary(Mmod)$ICs # Might as well go by majority rather than pick one
    mod <- if (sum(IC_comparison) >= 2) Amod else Mmod
  } else mod <- call_mod('ZZZ')
  
  mod
}


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
  # as_tibble(mod$states) %>% mutate(Time = time(mod$states), .before = 1)
  data <- cbind(actuals(mod), fitted(mod), mod$states, residuals(mod))
  colnames(data) <- c('actuals', 'fitted', colnames(mod$states), 'residuals')
  if (cut_non_actual_times) data <- data[which(!is.na(data$actuals))[1]:nrow(data)]
  as_tibble(data) %>% mutate(Time = time(data), .before = 1)
}


# Standardised-output (as in with the others in this package) prediction function for adam models
predict_adam <- function(mod, new_data, interval = NULL, only_new_data = F, draw_plot = F) {
  diff_f <- function(a, b) if (is.numeric(a)) a - b else difftime(a, b)
  horizon <- ceiling(as.numeric(diff_f(max(new_data), max(time(mod$fitted)))) / as.numeric(diff_f(time(mod$fitted)[2], time(mod$fitted)[1])))
  
  if (horizon < 1) {
    print('The requested prediction times are in the series past; returning the modelled values covering the requested interval')
    res <- mod$fitted[(index(mod$fitted) >= min(new_data)) & (index(mod$fitted) <= max(new_data))]
    tibble(Time = time(res), Pred = as.vector(res), Upper = NA, Lower = NA)
  } else {
    if (is.null(interval)) {
      ets_type <- str_match(mod$model, 'ETS\\((\\w+)\\)') # vv recommended in adam vignette
      interval <- if (!is.na(ets_type[2]) && ets_type[2][1] == 'M') 'simulated' else 'prediction'
    }
    
    preds <- forecast(mod, h = horizon, interval = interval)
    if (draw_plot) plot(preds)
    
    res <- tibble(Time = time(preds$mean), Pred = as.vector(preds$mean), Upper = as.vector(preds$upper), Lower = as.vector(preds$lower))
    
    if (only_new_data) {
      temp_res <- res %>% filter(Time %in% new_data)
      if (nrow(temp_res) != nrow(new_data)) {
        print('The requested prediction times do not line up with the model time unit; returning lined-up predictions covering the requested interval')
        res %>% filter((Time >= min(new_data)) & (Time <= max(new_data)))
      } else temp_res
    } else res
  }
}



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


