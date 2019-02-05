#### Stochastic Gradient Descent Methods

## Cost and Gradient functions for linear regression (simplest arguments for below methods)
lmf <- function(y, X, theta, ndata) { sum((X %*% theta - y) ^ 2) / (2 * ndata) }
lmgrad <- function(y, X, theta, ndata) { t(X) %*% (X %*% theta - y) / ndata }



### Gradient descent (GD) given a cost function f and its gradient grad
gd <- function(f, grad, y, X, theta0, npars, ndata, a, niters) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  
  theta[1, ] <- theta0
  cost[1] <- f(y, X, theta0, ndata)
  
  for (i in 2:niters) {
    theta[i, ] <- theta[i-1, ]-a*grad(y, X, theta[i-1, ], ndata)
    cost[i] <- f(y, X, theta[i, ], ndata)
  }
  
  return(list(theta=theta, cost=cost))
}



### Stochastic gradient descent (SGD) given a cost function f and its gradient grad
sgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples) {
    theta <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        theta[i,] <- theta[i - 1,] - a * grad(y[sIDs], X[sIDs,], theta[i - 1,], nsubsamples)
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}



### Stochastic gradient descent with momentum (MSGD) given a cost function f and its gradient grad
msgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0) {
    theta <- m <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    m[1,] <- m0
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        m[i,] <- b * m[i - 1,] + (1 - b) * grad(y[sIDs], X[sIDs,], theta[i - 1,], nsubsamples)
        theta[i,] <- theta[i - 1,] - a * m[i,]
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}



### Stochastic gradient descent with Nesterov accelerated gradient (NAGSGD) given a cost function f and its gradient grad
nagsgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0) {
    theta <- m <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    m[1,] <- m0
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        m[i,] <- b * m[i - 1,] + (1 - b) * grad(y[sIDs], X[sIDs,], theta[i - 1,] - a * b * m[i - 1,], nsubsamples)
        theta[i,] <- theta[i - 1,] - a * m[i,]
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}



### AdaGrad given a cost function f and its gradient grad
adagrad <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, epsilon, G0) {
    theta <- g <- diagG <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    g[1,] <- diagG[1,] <- G0 # Forcing this equality so that it will be true for all diagG that they are the sum of squared gs
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        theta[i,] <- theta[i - 1,] - a * g[i - 1,] / sqrt(diagG[i - 1,] + epsilon)
        g[i,] <- grad(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
        diagG[i,] <- diagG[i - 1,] + g[i,] ^ 2 # Works along with earlier commented line, removing the need to recompute all squares and sums
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}



### RMSProp given a cost function f and its gradient grad
rmsprop <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, c, epsilon, v0) {
    theta <- g <- v <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    v[1,] <- v0
    g[1,] <- grad(y[sIDs], X[sIDs,], theta0, nsubsamples)
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        v[i,] <- c * v[i - 1,] + (1 - c) * g[i - 1,] ^ 2
        theta[i,] <- theta[i - 1,] - a * g[i - 1,] / sqrt(v[i,] + epsilon)
        g[i,] <- grad(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}



### Adam given a cost function f and its gradient grad
adam <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, c, epsilon, m0, v0) {
    theta <- g <- m <- v <- matrix(data = NA, nrow = niters, ncol = npars)
    cost <- vector(mode = "numeric", length = niters)

    sIDs <- sample(1:ndata, nsubsamples)

    theta[1,] <- theta0
    m[1,] <- m0
    v[1,] <- v0
    g[1,] <- grad(y[sIDs], X[sIDs,], theta0, nsubsamples)
    cost[1] <- f(y[sIDs], X[sIDs,], theta0, nsubsamples)

    for (i in 2:niters) {
        sIDs <- sample(1:ndata, nsubsamples)
        m[i,] <- b * m[i - 1,] + (1 - b) * g[i - 1,]
        v[i,] <- c * v[i - 1,] + (1 - c) * g[i - 1,] ^ 2
        m_hat <- m[i,] / (1 - b ^ i)
        v_hat <- v[i,] / (1 - c ^ i)
        theta[i,] <- theta[i - 1,] - a * m_hat / (sqrt(v_hat) + epsilon)
        g[i,] <- grad(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
        cost[i] <- f(y[sIDs], X[sIDs,], theta[i,], nsubsamples)
    }

    return(list(theta = theta, cost = cost))
}
