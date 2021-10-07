library(rlang)
library(coda)


# Metropolis-Hastings Sampler
#   Generates a MH sample with the given parameters
#   Inputs: the number of required samples, the proposal density sampler and density functions,
#           the 0th sample and the log-likelihood function with only one argument (the proposed value)
#   Output: a LIST with the samples and the acceptance rate of the MH sampler
# Note:  
#   The required input is the log-likelihood instead of the likelihood because using the latter in
#   a MH chain can in many cases lead to numerical instabilities as the numbers may be quite small.
#   Doing calculations on the log scale and transforming back at the end can avoid underflow problems
mhSampler <- function(iterations = 10000, r_prop = rnorm(1, 0, 1), d_prop = function(s_to, s_from) { dnorm(s_to, 0, 1) }, s0 = 0, log_like) {
    sExpr <- enquo(r_prop)
    ss <- matrix(0, nrow = iterations, ncol = length(s0))
    colnames(ss) <- names(s0)
    ss[1,] <- as.vector(s0)
    rejRate <- 0

    for (i in 1:(iterations - 1)) {
        ss[i + 1,] <- as.vector(eval_tidy(sExpr))

        accRatio <- exp(log_like(ss[i + 1,]) + log(d_prop(ss[i,], ss[i + 1,]))
                        - log_like(ss[i,]) - log(d_prop(ss[i + 1,], ss[i,])))

        if (accRatio < 1) {
            if (rbinom(1, 1, accRatio) == 0) {
                ss[i + 1,] <- ss[i,]
                rejRate <- rejRate + 1
            }
        }
    }

    list(samples = mcmc(ss), accRate = 1 - rejRate / iterations)
}


# Perform burn-in and/or thinning on a matrix (or vector) of samples
#   Note: burn = 0 and thin = 1 for identity
burnNthin <- function(ss, burn = floor(nrow(as.matrix(ss)) / 10), thin = 2) {
    mcmc(
        as.matrix(ss)[(burn + 1):nrow(ss),][seq(1, nrow(ss) - burn, thin),],
        start = burn + 1,
        thin = thin
    )
}


# Metropolis-Hastings sampler with burn-in and thinning yielding the required final number of samples
mhSampler_burnNthin <- function(final_iterations, burn = floor(final_iterations / 10), thin = 2, ...) {
    mh_samples <- mhSampler(iterations = burn + final_iterations * thin, ...)
    mh_samples$samples <- burnNthin(mh_samples$samples, burn, thin)
    mh_samples$accRate <- 1 - as.numeric(rejectionRate(mh_samples$samples)[1])
    mh_samples
}


