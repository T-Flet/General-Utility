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


