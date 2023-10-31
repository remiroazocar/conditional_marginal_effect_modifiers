# This file specifies the simulation setup, generates the simulation study data
# and computes the true marginal estimands for active treatment vs. control

rm(list=ls()) # clear

# Set working directory
# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/conditional_marginal_effect_modifiers")

# Load packages
library("dplyr") # for data manipulation
library("MASS") # to simulate covariates from a multivariate normal distribution

set.seed(444) # set random seed for reproducibility

# define some simulation study settings
n.sim <- 500 # number of Monte Carlo replicates (simulated datasets) per scenario/method
n <- 5000 # number of subjects per trial
settings <- 1:3
outcome.models <- c("linear", "log-linear", "logistic")
# setting/model combinations for each simulation scenario
scenarios <- expand.grid(setting=settings, OM=outcome.models)
n.scenarios <- nrow(scenarios) # number of scenarios
save(n.sim, scenarios, n.scenarios, file="sim_settings.RData") # save simulation study settings

# function to generate simulated datasets for S=1 and S=2
sim.data <- function(n, setting=1, study=1, outcome.model="linear") {
  # Inputs
  # n: the number of subjects in each study
  # setting=1: treatment effect homogeneity, imbalanced covariate means and uncorrelated covariates
  # setting=2: treatment effect heterogeneity, balanced covariate means and different correlation structures
  # setting=3: treatment effect heterogeneity, imbalanced covariate means and different correlation structures
  # study=1 or study=2
  # outcome.model is "linear", "log-linear" or "logistic"
  #
  # treatment assignment indicator (1:1 allocation ratio in randomization)
  t <- c(rep(1, n/2), rep(0, n/2)) # 1: active; 0: control
  b0 <- -1 # intercept of outcome-generating model
  bX <- 1 # main (prognostic) coefficient for each covariate
  bXT <- 0.5 # interaction coefficient for conditional effect measure modifier
  bT <- 1.05 # main treatment coefficient (conditional trt effect at X=0) as per Austin and Stafford (2008)
  if (setting==1) {
    mean_X <- ifelse(study==1, 0, -1.4) # marginal mean for covariates depending on study
    rho <- cbind(c(1,0,0),c(0,1,0),c(0,0,1)) # pairwise correlation matrix for both studies in setting=1
    # three uncorrelated normally-distributed continuous baseline covariates
    x1 <- rnorm(n=n, mean=mean_X, sd=1)
    x2 <- rnorm(n=n, mean=mean_X, sd=1)
    x3 <- rnorm(n=n, mean=mean_X, sd=1)
    # treatment effect homogeneity at the individual level
    LP <- b0+bX*x1+bX*x2+bX*x3+bT*t # linear predictor (homogeneity)
  } else if (setting==2 | setting==3) {
    mean_X <- ifelse(study==1 & setting==3, 0, -1.4) # marginal mean for covariates depending on study/setting
    sd.vec <- rep(1, 3) # vector of standard deviations for each covariate
    # set correlation matrix depending on study
    if (study==1) {  
      rho <- cbind(c(1,0,0),c(0,1,0),c(0,0,1))
    } else if (study==2) {
      rho <- cbind(c(1,0,.4),c(0,1,.4),c(.4,.4,1))
    }
    cor2cov <- function(R, S) {
      # compute covariance matrix from correlation matrix R and vector
      # of standard deviations S
      sweep(sweep(R, 1, S, "*"),2,S,"*")
    }
    # covariance matrix required as input for mvrnorm
    cov.mat <- cor2cov(rho, sd.vec) 
    # three possibly correlated covariates generated from multivariate normal
    x <- as.data.frame(MASS::mvrnorm(n=n, mu=rep(mean_X,3), Sigma=cov.mat))
    x1 <- x[,1]
    x2 <- x[,2]
    x3 <- x[,3]
    # treatment effect heterogeneity at the individual level
    LP <- b0+bX*x1+bXT*x1*t+bX*x2+bX*x3+bT*t # linear predictor (heterogeneity)
  }
  if (outcome.model=="linear") {
    # continuous outcomes generated using linear model with error terms from a standard normal distribution
    eps <- rnorm(n, 0, 1) # 
    y <- LP + eps
  } else if (outcome.model=="log-linear") {
    # count outcomes generated from Poisson distribution with mean from log-linear model
    y <- rpois(n,exp(LP))
  } else if (outcome.model=="logistic") {
    # binary outcomes generated from Bernoulli distribution with mean from logistic model 
    y <- rbinom(n,1,exp(LP)/(1+exp(LP)))
  }
  # patient-level dataset
  data <- data.frame(x1=x1,x2=x2,x3=x3,trt=t,y=y)
  return(data)
}

# generate simulated datasets for each scenario
for (i in 1:n.scenarios) {
  print(i)
  # generate subject-level covariates, treatment and outcome for index trial S=1
  IPD.S1 <- replicate(n=n.sim, expr=sim.data(n=n, setting=scenarios$setting[i], study=1, 
                                             outcome.model=scenarios$OM[i]),
                      simplify=FALSE)
  # generate subject-level covariates, treatment and outcome for target trial S=2
  IPD.S2 <- replicate(n=n.sim, expr=sim.data(n=n, setting=scenarios$setting[i], study=2, 
                                             outcome.model=scenarios$OM[i]),
                      simplify=FALSE)
  # summarize the subject-level data for S=2 as aggregate-level data  
  if (scenarios$OM[i]=="linear") { # continuous outcomes, linear model, mean difference
    ALD.S2 <- lapply(1:n.sim, function(j) {
      as.data.frame(cbind(
        # summarize covariates as means with standard deviations
        summarise(IPD.S2[[j]], mean.x1=mean(x1), mean.x2=mean(x2), mean.x3=mean(x3),
                  sd.x1=sd(x1), sd.x2=sd(x2), sd.x3=sd(x3),
                  # point estimate for marginal mean difference for B vs. C and standard error
                  hat.Delta=summary(lm(y~trt, data=IPD.S2))$coef[2],
                  SE.hat.Delta=sqrt(vcov(lm(y~trt, data=IPD.S2))[2,2]))))
    } )
  } else if (scenarios$OM[i]=="log-linear") { # count outcomes, log-linear model, log risk ratio
    ALD.S2 <- lapply(1:n.sim, function(j) {
      as.data.frame(cbind(
        # summarize covariates as means with standard deviations
        summarise(IPD.S2[[j]], mean.x1=mean(x1), mean.x2=mean(x2), mean.x3=mean(x3),
                  sd.x1=sd(x1), sd.x2=sd(x2), sd.x3=sd(x3),
                  # point estimate for marginal log risk ratio for B vs. C and standard error
                  hat.Delta=summary(glm(y~trt, data=IPD.S2, family=poisson(link="log")))$coef[2],
                  SE.hat.Delta=sqrt(vcov(glm(y~trt, data=IPD.S2, family=poisson(link="log")))[2,2]))))
    } )    
  } else if (scenarios$OM[i]=="logistic") { # binary outcomes, logistic model, log odds ratio 
    ALD.S2 <- lapply(1:n.sim, function(j) {
      as.data.frame(cbind(
        # summarize covariates as means with standard deviations
        summarise(IPD.S2[[j]], mean.x1=mean(x1), mean.x2=mean(x2), mean.x3=mean(x3),
                  sd.x1=sd(x1), sd.x2=sd(x2), sd.x3=sd(x3),
                  # point estimate for marginal log risk ratio for B vs. C and standard error
                  hat.Delta=summary(glm(y~trt, data=IPD.S2, family=binomial(link="logit")))$coef[2],
                  SE.hat.Delta=sqrt(vcov(glm(y~trt, data=IPD.S2, family=binomial(link="logit")))[2,2]))))
    } )       
  }
  file.id <- paste0("setting", scenarios$setting[i], "OM", scenarios$OM[i]) 
  save(IPD.S1, file=paste0("Data/IPD_S1_", file.id, ".RData"))
  save(IPD.S2, file=paste0("Data/IPD_S2_", file.id, ".RData"))
  save(ALD.S2, file=paste0("Data/ALD_S2_", file.id, ".RData"))  
} 

# compute true marginal estimand for active treatment vs. control in each study/scenario
compute.truth <- function(n, setting, study, outcome.model) {
  dataset <- sim.data(n=n, setting=setting, study=study, outcome.model=outcome.model)   
  mu1 <- mean(dataset[which(dataset$trt==1),]$y)
  mu0 <- mean(dataset[which(dataset$trt==0),]$y)
  if (outcome.model=="linear") {
    truth.marginal <- mu1-mu0 # true marginal mean difference for active vs. control  
  } else if (outcome.model=="log-linear") {
    truth.marginal <- log(mu1)-log(mu0) # true marginal log risk ratio for active vs. control
  } else if (outcome.model=="logistic") {
    truth.marginal <- log((mu1/(1-mu1))/(mu0/(1-mu0))) # true marginal log odds ratio for active vs. control
  }
  return(truth.marginal)
}

# true marginal mean differences for active treatment versus control in each study/setting
true.MD.S1.setting1 <- compute.truth(n=50000000, setting=1, study=1, outcome.model="linear")
true.MD.S2.setting1 <- compute.truth(n=50000000, setting=1, study=2, outcome.model="linear")
true.MD.S1.setting2 <- compute.truth(n=50000000, setting=2, study=1, outcome.model="linear")
true.MD.S2.setting2 <- compute.truth(n=50000000, setting=2, study=2, outcome.model="linear")
true.MD.S1.setting3 <- compute.truth(n=50000000, setting=3, study=1, outcome.model="linear")
true.MD.S2.setting3 <- compute.truth(n=50000000, setting=3, study=2, outcome.model="linear")

# true marginal log risk ratio for active treatment versus control in each study/setting
true.log.RR.S1.setting1 <- compute.truth(n=50000000, setting=1, study=1, outcome.model="log-linear")
true.log.RR.S2.setting1 <- compute.truth(n=50000000, setting=1, study=2, outcome.model="log-linear")
true.log.RR.S1.setting2 <- compute.truth(n=50000000, setting=2, study=1, outcome.model="log-linear")
true.log.RR.S2.setting2 <- compute.truth(n=50000000, setting=2, study=2, outcome.model="log-linear")
true.log.RR.S1.setting3 <- compute.truth(n=50000000, setting=3, study=1, outcome.model="log-linear")
true.log.RR.S2.setting3 <- compute.truth(n=50000000, setting=3, study=2, outcome.model="log-linear")

# true marginal log odds ratio for active treatment versus control in each study/setting
true.log.OR.S1.setting1 <- compute.truth(n=50000000, setting=1, study=1, outcome.model="logistic")
true.log.OR.S2.setting1 <- compute.truth(n=50000000, setting=1, study=2, outcome.model="logistic")
true.log.OR.S1.setting2 <- compute.truth(n=50000000, setting=2, study=1, outcome.model="logistic")
true.log.OR.S2.setting2 <- compute.truth(n=50000000, setting=2, study=2, outcome.model="logistic")
true.log.OR.S1.setting3 <- compute.truth(n=50000000, setting=3, study=1, outcome.model="logistic")
true.log.OR.S2.setting3 <- compute.truth(n=50000000, setting=3, study=2, outcome.model="logistic")
