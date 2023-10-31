# This file contains functions to estimate the simulation study performance measures

### Functions to estimate performance measures of simulation study

# bias estimate
bias <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(theta.hat)/nsim - theta
  return(est)
}

# Monte Carlo SE of bias estimate
bias.mcse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  mcse <- sqrt(1/(nsim*(nsim-1))*tmp)
  return(mcse)
}

# coverage estimate
coverage <- function(theta.hat.low, theta.hat.upp, theta) {
  nsim <- length(theta.hat.low)
  est <- sum(ifelse(theta>=theta.hat.low & theta<=theta.hat.upp,1,0))/nsim
  return(est)
}

# Monte Carlo SE of coverage estimate
coverage.mcse <- function(coverage, nsim) {
  mcse <- sqrt((coverage*(1-coverage))/nsim)
  return(mcse)
}

# MSE estimate
mse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum((theta.hat-theta)^2)/nsim
  return(est)
}

# Monte Carlo SE of MSE estimate
mse.mcse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  tmp <- (theta.hat-theta)^2
  mse.est <- sum(tmp)/nsim
  mcse <- sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# Monte Carlo SE of any continuous performance metric
mcse.estimate <- function(perf.measure) {
  nsim <- length(perf.measure)
  perf.measure.mean <- sum(perf.measure)/nsim
  mcse <- sqrt(sum((perf.measure-perf.measure.mean)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# Empirical standard error 
empse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  est <- sqrt(tmp/(nsim-1))
  return(est)
}

# EmpSE MCSE
empse.mcse <- function(empse, nsim) {
  mcse <- empse/(sqrt(2*(nsim-1)))
  return(mcse)
}