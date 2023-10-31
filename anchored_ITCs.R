# This file performs the anchored indirect comparison methodologies to analyze the simulated data

rm(list=ls()) # clear

# Set working directory
# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/conditional_marginal_effect_modifiers")
load(file="sim_settings.RData") # load simulation study settings

# load packages
library("copula") # to simulate individual-level covariates in parametric G-computation 
library("boot") # for non-parametric bootstrap
library("parallel") # to detect cores
library("doSNOW") # for parallel cluster

set.seed(444) # set seed for reproducibility

# settings for the covariate-adjusted indirect comparison methods
resamples <- 500 # number of resamples for non-parametric bootstrap
n.star <- 1000 # number of simulated covariate profiles for parametric G-computation

# simulated patient-level (S=1) and aggregate-level (S=2) datasets for all scenarios
IPD.S1.all <- vector(mode="list", n.scenarios)
ALD.S2.all <- vector(mode="list", n.scenarios)

# load simulated datasets
for (i in 1:n.scenarios) {
  file.id <- paste0("setting", scenarios$setting[i], "OM", scenarios$OM[i])  
  load(paste0("Data/IPD_S1_", file.id, ".RData")) # load index (S=1) patient-level data
  load(paste0("Data/ALD_S2_", file.id, ".RData")) # load target (S=2) aggregate-level data
  IPD.S1.all[[i]] <- IPD.S1
  ALD.S2.all[[i]] <- ALD.S2
}
  
# Anchored indirect treatment comparison methodologies

## MAIC (matching-adjusted indirect comparison) 
maic.wrapper <- function(data.S1, data.S2, resamples, outcome.model="linear") {
  # Inputs
  # data.S1: index trial (S=1) individual-level data 
  # data.s2: target trial (S=2) aggregate-level data
  # resamples: number of resamples for non-parametric bootstrap
  # outcome.model: "linear", "log-linear" or "logistic"
  #
  # non-parametric bootstrap with replacement over S=1 data
  maic.boot <- function(data.S1, indices) {
    dat.S1 <- data.S1[indices,]  
    # weight estimation uses method of moments as per Signorovitch et al. approach
    # center S=1 covariates on S=2 means
    dat.S1$x1 <- dat.S1$x1 - data.S2$mean.x1
    dat.S1$x2 <- dat.S1$x2 - data.S2$mean.x2
    dat.S1$x3 <- dat.S1$x3 - data.S2$mean.x3
    alpha <- rep(1, 3) # arbitrary starting point for the optimizer
    # objective function to be minimized for standard method of moments
    Q <- function(alpha, X) {
      return(sum(exp(X %*% alpha)))
    }
    z <- as.matrix(dat.S1[,c("x1","x2","x3")]) # matrix of centered covariates
    # objective function minimized using BFGS
    Q.min <- optim(fn=Q, X=z, par=alpha, method="BFGS")
    # finite solution is the logistic regression parameters
    hat.alpha <- Q.min$par
    log.hat.w <- z %*% hat.alpha
    dat.S1$hat.w <- exp(log.hat.w) # estimated weights
    if (outcome.model=="linear") {
      # fit weighted simple linear regression of outcome on treatment to S=1  
      outcome.fit <- lm(y~trt, weights=hat.w, data=dat.S1)
    } else if (outcome.model=="log-linear") {
      # fit weighted simple Poisson regression of outcome on treatment to S=1
      outcome.fit <- glm(y~trt, weights=hat.w, data=dat.S1, family=poisson(link="log"))
    } else if (outcome.model=="logistic") {
      # fit weighted simple logistic regression of outcome on treatment to S=1
      outcome.fit <- glm(y~trt, weights=hat.w, data=dat.S1, family=quasibinomial(link="logit"))
    }  
    # treatment coefficient of fitted model is marginal treatment effect estimate for A vs. C
    delta.AC <- summary(outcome.fit)$coef[2]
    return(delta.AC)    
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.S1, statistic=maic.boot, R=resamples)
  # bootstrap mean of marginal treatment effect estimate for A vs. C
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of marginal treatment effect estimate for A vs. C
  hat.var.Delta.AC <- var(boot.object$t)
  # published marginal B vs. C treatment effect estimate
  hat.Delta.BC <- with(data.S2, hat.Delta)
  # published standard error estimate for marginal B vs. C treatment effect
  hat.SE.Delta.BC <- with(data.S2, SE.hat.Delta)
  # combine relative effects and standard errors for A vs. B indirect comparison
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC 
  hat.var.Delta.AB <- hat.var.Delta.AC + (hat.SE.Delta.BC)^2    
  list(hat.Delta.AB, hat.var.Delta.AB)
}

## Parametric G-computation with maximum likelihood estimation
param.gcomp.wrapper <- function(data.S1, data.S2, n.star, resamples, setting=1, outcome.model="linear") { 
  # Inputs
  # data.S1: index trial (S=1) individual-level data 
  # data.s2: target trial (S=2) aggregate-level data
  # n.star: number of simulated individual-level covariate profiles for S=2 
  # resamples: number of resamples for non-parametric bootstrap
  # outcome.model: "linear", "log-linear" or "logistic"
  # setting=1: homogeneity working model
  # setting=2 or setting=3: heterogeneity working model 
  #
  # matrix of pairwise correlations between S=1 subject-level covariates
  rho <- cor(data.S1[,c("x1","x2","x3")])
  # individual-level covariate simulation for S=2 using multivariate Gaussian copula
  cop <- normalCopula(param=c(rho[1,2],rho[1,3],rho[2,3]), dim=3, dispstr="un") # S=1 pairwise correlations
  # simulate S=2 covariates from copula with normally-distributed marginals
  mvd <- mvdc(copula=cop, margins=c("norm", "norm", "norm"),
              # marginal distribution means and standard deviations from S=2 publication
              paramMargins=list(list(mean=data.S2$mean.x1, sd=data.S2$sd.x1),
                                list(mean=data.S2$mean.x2, sd=data.S2$sd.x2),       
                                list(mean=data.S2$mean.x3, sd=data.S2$sd.x3)))
  # data frame of S=2 simulated covariate profiles
  x.star <- as.data.frame(rMvdc(n.star, mvd))
  colnames(x.star) <- c("x1", "x2", "x3")  
  param.gcomp.boot <- function(data.S1, indices) {
    dat.S1 <- data.S1[indices,]
    if (outcome.model=="linear") {
      # covariate-adjusted linear regression fitted to S=1 subject-level data using maximum-likelihood      
      if (setting==1) { # homogeneity working model (assumed correctly specified)
        mod <- lm(y~trt+x1+x2+x3, data=dat.S1)  
      } else { # heterogeneity working model (assumed correctly specified)
        mod <- lm(y~trt+x1*trt+x2+x3, data=dat.S1)
      }
    } else if (outcome.model=="log-linear") {
      # covariate-adjusted Poisson regression fitted to S=1 subject-level data using maximum-likelihood
      if (setting==1) { # homogeneity working model
        mod <- glm(y~trt+x1+x2+x3, data=dat.S1, family=poisson(link="log")) 
      } else { # heterogeneity working model
        mod <- glm(y~trt+x1*trt+x2+x3, data=dat.S1, family=poisson(link="log"))
      }
    } else if (outcome.model=="logistic") {
      # covariate-adjusted logistic regression fitted to S=1 subject-level data using maximum-likelihood
      if (setting==1) { # homogeneity working model
        mod <- glm(y~trt+x1+x2+x3, data=dat.S1, family=binomial(link="logit"))  
      } else { # heterogeneity working model
        mod <- glm(y~trt+x1*trt+x2+x3, data=dat.S1, family=binomial(link="logit"))
      }
    }
    # counterfactual datasets with fixed covariates
    data.1 <- data.0 <- x.star
    data.1$trt <- 1 # intervene so that everyone receives A
    data.0$trt <- 0 # intervene so that everyone receives C
    # predict conditional outcome mean for simulated S=2 profiles under treatment A
    mu1 <- predict(mod, newdata=data.1, type="response")
    # marginal mean outcome on natural scale for treatment A
    hat.mu1 <- mean(mu1)
    # predict conditional outcome mean for simulated S=2 profiles under treatment C
    mu0 <- predict(mod, newdata=data.0, type="response")
    # marginal mean outcome on natural scale for treatment C
    hat.mu0 <- mean(mu0)
    if (outcome.model=="linear") {
      # covariate-adjusted marginal mean difference point estimate for A vs. C      
      delta.AC <- hat.mu1-hat.mu0
    } else if (outcome.model=="log-linear") {
      # covariate-adjusted marginal log risk ratio point estimate for A vs. C
      delta.AC <- log(hat.mu1)-log(hat.mu0)
    } else if (outcome.model=="logistic") {
      # covariate-adjusted marginal log odds ratio point estimate for A vs. C
      delta.AC <- log((hat.mu1/(1-hat.mu1))/(hat.mu0/(1-hat.mu0)))  
    }
    return(delta.AC)  
  }    
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.S1, statistic=param.gcomp.boot, R=resamples)
  # bootstrap mean of marginal treatment effect estimate for A vs. C
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of marginal treatment effect estimate for A vs. C
  hat.var.Delta.AC <- var(boot.object$t)
  # published marginal B vs. C treatment effect estimate
  hat.Delta.BC <- with(data.S2, hat.Delta)
  # published standard error estimate for marginal B vs. C treatment effect
  hat.SE.Delta.BC <- with(data.S2, SE.hat.Delta)
  # combine relative effects and standard errors for A vs. B indirect comparison
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC 
  hat.var.Delta.AB <- hat.var.Delta.AC + (hat.SE.Delta.BC)^2    
  list(hat.Delta.AB, hat.var.Delta.AB)  
}  

# Bucher method (unadjusted ITC: does not adjust for covariates)
bucher.wrapper <- function(data.S1, data.S2, outcome.model="linear") {
  if (outcome.model=="linear") {
    # fit simple linear regression of outcome on treatment to S=1  
    outcome.fit <- lm(y~trt, data=data.S1)
  } else if (outcome.model=="log-linear") {
    # fit simple Poisson regression of outcome on treatment to S=1
    outcome.fit <- glm(y~trt, data=data.S1, family=poisson(link="log"))
  } else if (outcome.model=="logistic") {
    # fit simple logistic regression of outcome on treatment to S=1
    outcome.fit <- glm(y~trt, data=data.S1, family=binomial(link="logit"))
  }  
  # treament coefficient is point estimate of marginal treatment effect for A vs. C
  hat.Delta.AC <- summary(outcome.fit)$coef[2]
  # variance estimate for marginal treatment effect of A vs. C
  hat.var.Delta.AC <- vcov(outcome.fit)[2,2]
  # published marginal B vs. C treatment effect estimate
  hat.Delta.BC <- with(data.S2, hat.Delta)
  # published standard error estimate for marginal B vs. C treatment effect
  hat.SE.Delta.BC <- with(data.S2, SE.hat.Delta)
  # combine relative effects and standard errors for A vs. B indirect comparison
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC 
  hat.var.Delta.AB <- hat.var.Delta.AC + (hat.SE.Delta.BC)^2 
  list(hat.Delta.AB, hat.var.Delta.AB)
}
  
# set up cluster for parallel computing
num.cores <- detectCores()-1 # leave one core available
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
pb <- txtProgressBar(max=n.sim, style=3) # progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# function to combine lists in parallelization
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# For a given scenario/method, run anchored indirect comparison methods in parallel for all replicates
for(i in 1:n.scenarios) {
  IPD.S1 <- IPD.S1.all[[i]]
  ALD.S2 <- ALD.S2.all[[i]]
  file.id <- paste0("setting", scenarios$setting[i], "OM", scenarios$OM[i])  
  # matching-adjusted indirect comparison (MAIC)
  maic.results <- foreach(j=1:n.sim, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list()), .options.snow=opts,
                          .packages=c("boot")) %dopar% {
                            results <- maic.wrapper(data.S1=IPD.S1[[j]], data.S2=ALD.S2[[j]],
                                                    resamples=resamples,
                                                    outcome.model=scenarios$OM[i])
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  save(means, file=paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file.id, ".RData"))
  # parametric G-computation with maximum-likelihood estimation
  gcomp.results <- foreach(j=1:n.sim, .combine='comb', .multicombine=TRUE,
                           .init=list(list(),list()), .options.snow=opts,
                           .packages=c("boot", "copula")) %dopar% {
                                results <- param.gcomp.wrapper(data.S1=IPD.S1[[j]], data.S2=ALD.S2[[j]],
                                                               n.star=n.star, resamples=resamples,
                                                               setting=scenarios$setting[i],
                                                               outcome.model=scenarios$OM[i])
                                return(results)
                           }
  close(pb)
  means <- unlist(gcomp.results[[1]])
  variances <- unlist(gcomp.results[[2]])
  save(means, file=paste0("Results/GComp/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/GComp/variances_", file.id, ".RData"))
  # Bucher method (unadjusted anchored indirect treatment comparison)
  bucher.results <- foreach(j=1:n.sim, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts) %dopar% {
                              results <- bucher.wrapper(data.S1=IPD.S1[[j]], data.S2=ALD.S2[[j]],
                                                        outcome.model=scenarios$OM[i])
                              return(results)
                            }
  close(pb)
  means <- unlist(bucher.results[[1]])
  variances <- unlist(bucher.results[[2]])
  save(means, file=paste0("Results/Bucher/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/Bucher/variances_", file.id, ".RData"))
}

stopCluster(cluster)