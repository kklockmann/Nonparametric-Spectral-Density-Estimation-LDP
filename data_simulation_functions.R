##############################################################
#EXAMPLES OF STATIONARY GAUSSIAN PROCESSES
##############################################################
library(MASS) #mvrnorm
library(dtt)

#GLOBAL VARIABLES

#INPUT
## p: length of vector
## n: number of iid. samples
## sd: standard deviation

#OUTPUT
## Y: pxn dimensional data matrix
## sdf: spectral density function of Y
## acf: auto-covariance function of Y

phi <- 0.8
n <- 10  # e.g., for a 10×10 correlation matrix

#EX1: ARMA(0.8) process
example1 <- function(p, n, sd) {
  x = seq(0, 1, length = p)
  pp = 1
  coef = c(0.8)
  R.mat <- outer(1:p, 1:p, function(i, j) coef^abs(i - j)) #alternativ command
  Sigma = sd^2 * R.mat
  acf = Sigma[1, ]
  sdf = acf2sdf(acf)
  Y = matrix(mvrnorm(n, mu = numeric(p), Sigma = Sigma), n, p)
  return(list(Y = Y, sdf = sdf, acf = acf))
}

#EX2: Hölder continuos sdf
example2 <- function(p, n, sd, gamma) {
  x = seq(0, 1, length = p)
  sdf = (abs(cos(pi * x))^gamma + 0.45)
  #normalizing factor such that integral of sdf = 1
  C = mean(abs(cos(pi * seq(0, 1, length = 100000)))^gamma + 0.45)
  sdf = sd^2 * sdf / C
  acf = sdf2acf(sdf)
  Sigma = toeplitz(acf)
  Y = matrix(mvrnorm(n, mu = numeric(p), Sigma = Sigma), n, p)
  return(list(Y = Y, sdf = sdf, acf = acf))
}

#EX3: acf with polynomial decay
#gamma = decay of acf
example3 <- function(p, n, sd, gamma) {
  x = seq(0, 1, length = p)
  acf = c(sd^2, sd^2 / (1 + seq(1, p - 1))^gamma)
  Sigma = toeplitz(acf)
  sdf = acf2sdf(acf)
  Y = matrix(mvrnorm(n, mu = numeric(p), Sigma = Sigma), n, p)
  return(list(Y = Y, sdf = sdf, acf = acf))
}

acf2sdf <- function(x) {
  N = length(x)
  y = 2 * dtt(x, type = "dct", variant = 1) #*sqrt(2/(N-1))
  #y=Re(fft(c(x,x[(N-1):2]))) #the same
  return (y)
}

sdf2acf <- function(x) {
  N = length(x)
  y = dtt(x, type = "dct", variant = 1) / (N - 1)
  return (y)
}
