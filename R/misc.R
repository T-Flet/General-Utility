library(tidyverse)



library(metap)
# Adjust and combine a set of p-values without any assumption of independence/correlation pattern or evidence distribution
#   Methods achieving this: Bonferroni-Holm for adjustment and logistic quantiles for combining (see Loughin 2004)
adjust_combine_pvs <- function(ps) logitp(p.adjust(ps, 'holm'))


# Compute the overlap of two intervals
  # IMPORTANT NOTE: This is not a vectorised function, therefore DO
  #   EITHER: rowwise %>% mutate(... intervalOverlap(...) ...) %>% ungroup
  #   OR: mutate(... Vectorize(intervalOverlap)(...) ...)
  # If this is not done it will return the result of the first evaluation and that will be it for all rows
intervalOverlap <- function(a, b) max(0, min(a[2], b[2]) - max(a[1], b[1]))


# Area under the tail on the opposite side of 0 from the mean of the difference between two normal distributions
#   I.e. the probability that it is not the sign it is
p_normal_diff_other_sign <- function(mu1, sd1, mu2, sd2) {
  mu <- mu1 - mu2
  sd <- sqrt(sd1^2 + sd2^2)
  pnorm(0, mean = mu, sd = sd, lower.tail = mu >= 0)
}


