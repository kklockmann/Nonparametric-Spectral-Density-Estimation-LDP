library(extraDistr)

#truncates data at +-tau and adds laplace noise Lap(fac * tau / alpha)
truncLap <- function(data, alpha, tau, fac) {
  p <- length(data)
  Lap.noise <- rlaplace(p, 0, fac * tau / alpha)
  Z <- data
  Z[Z > tau] <- tau
  Z[Z < -tau] <- -tau
  Z <- Z + Lap.noise
  return(Z)
}

##########################################################
#Functions for covariance coefficient estimation
#with sequentially interactive (SI) mechanism
##########################################################

#j=index of covariance coefficient (0: variance 1: simga_1, ...)
#tau, tau.tilde = truncation parameters
priv.data.SI.acf <- function(data, j, alpha, tau, tau.tilde) {
  p <- length(data)
  if (j != 0) {
    Z <- truncLap(data, alpha, tau, fac = 4)
    Z.bar <- data[(j + 1):p] * Z[1:(p - j)]
    Z.bar <- truncLap(Z.bar, alpha, tau.tilde, fac = 4)
    return(Z.bar)
  }
  if (j == 0) {
    Z <- data^2
    Z.bar <- truncLap(Z, alpha, tau, fac = 2)
    return(Z.bar)
  }
}

#estimates covariance coefficient from SI-privatized data
sigma.hat.SI <- function(data) {
  return(mean(data))
}

############################################################
#Functions for pointwise spectral density estimation
#with sequentially interactive (SI) mechanism
############################################################

#omega = frequency to be estimated, value in [0,pi]
#K = number of Fourier coefficients to be included in the estimate
#    e.g. K=0: only sigma_0 included
#    e.g. K=1: only sigma_(-1),sigma_0,sigma_1 included
#tau, tau.tilde = truncation parameters
priv.data.SI.sdf0 <- function(data, omega, alpha, tau, tau.tilde, K) {
  Z <- truncLap(data, alpha, tau, fac = 4)
  n <- length(data)
  V <- rep(0 + 0i, n)
  
  aw <- ak_weights(K)
  k_vals <- aw$k
  a_vals <- aw$a
  
  lap_scale <- 4 * tau.tilde / alpha
  
  for (i in (K + 1):n) {
    s <- 0 + 0i
    
    for (idx in seq_along(k_vals)) {
      k <- k_vals[idx]
      
      if (k == 0)
        next
      
      j <- i - abs(k)
      if (j >= 1) {
        s <- s +
          a_vals[idx] *
          data[i] * Z[j] *
          exp(-1i * omega * k)
      }
    }
    
    V[i] <- data[i]^2 + s
  }
  if (max(abs(Im(V))) < 0.01) {
    V = Re(V)
  } else {
    cat("imgainary number in V")
  }
  Z.tilde <- truncLap(V, alpha, tau.tilde, fac = 4)
  return(list(Z = Z, Z.tilde = Z.tilde))
}

ak_weights <- function(K) {
  k <- -K:K
  a <- numeric(length(k))
  
  for (j in seq_along(k)) {
    kk <- abs(k[j])
    if (kk <= K / 2) {
      a[j] <- 1
    } else {
      a[j] <- 2 * (1 - kk / K)
    }
  }
  
  return(list(k = k, a = a))
}

#estimates pointwise spectral density from SI-privatized data
sdf0.hat.SI <- function(data.priv) {
  res = mean(data.priv$Z.tilde) / (2 * pi)
  return(res)
}



############################################################
#Functions for global spectral density estimation
#no privacy
############################################################

#omega = frequency to be estimated, value in [0,pi]
#K = number of Fourier coefficients to be included in the estimate
sdf.hat.global <- function(data, K, omegas) {
  p = length(data)
  sigma.hat = c()
  for (j in 0:K) {
    X1 = data[1:(p - j)]
    X2 = data[(j + 1):p]
    sigma.hat <- c(sigma.hat, mean(X1 * X2))
  }
  ks = (-K):K
  sdf.res = c()
  for (o in omegas) {
    sdf.res = c(sdf.res, sum(sigma.hat[abs(ks) + 1] *
                               exp(-1i * ks * o)))
  }
  
  sdf.res = sdf.res / (2 * pi)
  if (max(abs(Im(sdf.res))) < 0.01) {
    sdf.res = Re(sdf.res)
  } else {
    cat("imgainary number in sdf.res")
  }
  return(sdf.res)
}




##########################################################
#Functions for covariance coefficient estimation
#with non-interactive (NI) mechanism
#based on (Kroll, 2024)
##########################################################

priv.data.NI <- function(data, alpha, tau) {
  p <- length(data)
  
  Lap.noise <- rlaplace(p, 0, 2 * tau / alpha)
  
  Z <- data
  Z[Z > tau] <- tau
  Z[Z < -tau] <- -tau
  Z <- Z + Lap.noise
  return(Z)
}


sigma.hat.NI <- function(data, j, alpha, tau, debias = FALSE) {
  p = length(data)
  Z1 = data[1:(p - j)]
  Z2 = data[(j + 1):p]
  if (debias) {
    sigma.hat <- sum(Z1 * Z2) / (p - j)
  } else{
    sigma.hat <- sum(Z1 * Z2) / p
  }
  
  if (j == 0) {
    sigma.hat = sigma.hat - 8 * tau^2 / alpha^2
  }
  return(sigma.hat)
}


##########################################################
#Functions for global spectral density estimation
#with non-interactive (NI) mechanism
#based on (Kroll, 2024)
##########################################################
#omegas = vector of frequencies to be estimated, value in [0,pi]
#K = number of Fourier coefficients to be included in the estimate

sdf.hat.global.NI <- function(data, K, omegas, alpha, tau, debias = FALSE) {
  jj = 0:K
  sigmas = sapply(jj, function(j) {
    sigma.hat.NI(data, j, alpha, tau, debias)
  })
  
  sdf.res = c()
  for (o in omegas) {
    val <- 0 + 0i
    for (k in-K:K) {
      val <- val +
        sigmas[abs(k) + 1] *
        exp(-1i * k * o)
    }
    sdf.res = c(sdf.res, val)
  }
  sdf0 = sdf.res / (2 * pi)
  if (max(abs(Im(sdf0))) < 0.01) {
    sdf0 = Re(sdf0)
  } else {
    cat("imgainary number in sdf0")
  }
  return(sdf0)
}