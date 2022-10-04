# Set working directory
# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/conditional_marginal_effect_modifiers")

rm(list=ls()) # clear 

# Load packages
library("doSNOW") # for parallel cluster
library("parallel") # to detect cores
library("boot") # for non-parametric bootstrap with replacement
library("ggplot2") # plot simulation study results
library("ggridges") # plot simulation study results
library("gridExtra") # plot simulation study results

set.seed(555) # set random seed for reproducibility

#### SIMULATE DATA #####

n.sim <- 2000 # number of Monte Carlo replicates (simulated datasets)

# function to generate data for simulation study
sim.data <- function(n, mean_X, outcome.model="logistic") {
  # Inputs: number of subjects in each trial and mean of each baseline covariate 
  #
  # three normally-distributed continuous baseline covariates
  x1 <- rnorm(n=n, mean=mean_X, sd=1)
  x2 <- rnorm(n=n, mean=mean_X, sd=1)
  x3 <- rnorm(n=n, mean=mean_X, sd=1)
  # treatment assignment indicator (1:1 allocation ratio in randomization)
  t <- c(rep(1, n/2), rep(0, n/2)) # 1: active; 0: control
  b0 <- -1 # intercept of outcome-generating model
  bX <- 1 # main (prognostic) coefficient for each covariate
  bT <- 1.048569 # main treatment coefficient as per Austin and Stafford (2008)
  if (outcome.model=="logistic") {
    # binary outcomes generated using logistic model
    y <- rbinom(n,1,exp(b0+bX*x1+bX*x2+bX*x3+bT*t)/(1+exp(b0+bX*x1+bX*x2+bX*x3+bT*t)))
  } else if (outcome.model=="linear") {
    # continuous outcomes generated using linear model
    LP <- b0+bX*x1+bX*x2+bX*x3+bT*t # linear predictor
    eps <- rnorm(n, 0, 1) # normally-distributed error terms
    y <- LP + eps
  }
  # patient-level dataset
  data <- data.frame(x1=x1,x2=x2,x3=x3,trt=t,y=y)
  return(data)
}

#### COVARIATE ADJUSTMENT METHODOLOGIES ####

### MAIC (matching-adjusted indirect comparison) 
maic.wrapper <- function(data.S1, data.S2, resamples, outcome.model="logistic") {
  # Inputs: patient-level datasets and number of resamples in non-parametric bootstrap
  #
  # non-parametric bootstrap with replacement over S=1 data
  maic.boot <- function(data.S1, indices) {
    dat.S1 <- data.S1[indices,]  
    # weight estimation uses method of moments as per Signorovitch et al. approach
    # center S=1 covariates on S=2 means
    dat.S1$x1 <- dat.S1$x1 - mean(data.S2$x1) 
    dat.S1$x2 <- dat.S1$x2 - mean(data.S2$x2)
    dat.S1$x3 <- dat.S1$x3 - mean(data.S2$x3)
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
    dat.S1$hat.w <- exp(log.hat.w)
    if (outcome.model=="logistic") {
      # fit weighted logistic regression to S=1 IPD  
      outcome.fit <- glm(y~trt, weights=hat.w, data=dat.S1, family=quasibinomial)
    } else if (outcome.model=="linear") {
      # fit weighted linear regression to S=1 IPD
      outcome.fit <- lm(y~trt, weights=hat.w, data=dat.S1)
    }
    # treatment coefficient of fitted model is marginal treatment effect estimate for A vs. C
    # marginal log odds ratio for logistic, marginal mean difference for linear
    delta.AC <- summary(outcome.fit)$coef[2]
    return(delta.AC)    
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.S1, statistic=maic.boot, R=resamples)
  # bootstrap mean of marginal treatment effect estimate for A vs. C
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of marginal treatment effect estimate for A vs. C
  hat.var.Delta.AC <- var(boot.object$t)
  # analysis of S=2 IPD
  if (outcome.model=="logistic") {
    outcome.fit <- glm(y~trt, data=data.S2, family=binomial)
  } else if (outcome.model=="linear") {
    outcome.fit <- lm(y~trt, data=data.S2)
  }
  # marginal treatment effect estimate for B vs. C
  hat.Delta.BC <- summary(outcome.fit)$coef[2]
  # variance estimate for marginal treatment effect of B vs. C
  hat.var.Delta.BC <- vcov(outcome.fit)[2,2]
  # indirect treatment comparison (A vs. B)
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
  list(hat.Delta.AB, hat.var.Delta.AB)
}

### Parametric G-computation with maximum likelihood estimation
param.gcomp.wrapper <- function(data.S1, data.S2, resamples, outcome.model="logistic") { 
  # Inputs: patient-level datasets and number of resamples in non-parametric bootstrap
  # 
  # non-parametric bootstrap with replacement over S=1 data
  param.gcomp.boot <- function(data.S1, indices) {
    dat.S1 <- data.S1[indices,]
    if (outcome.model=="logistic") {
      # logistic regression outcome model fitted to S=1 IPD with maximum-likelihood
      mod <- glm(y~trt+x1+x2+x3, data=dat.S1, family=binomial)  
    } else if (outcome.model=="linear") {  
      # linear regression outcome model fitted to S=1 IPD with maximum-likelihood
      mod <- lm(y~trt+x1+x2+x3, data=dat.S1)  
    }
    # potential (counterfactual) study S=2 under treatment A
    data.S2$trt <- 1 # 
    # predict outcomes on natural scale for each individual under A
    mu1 <- predict(mod, newdata=data.S2, type="response")
    # marginal mean outcome on natural scale under treatment A
    hat.mu1 <- mean(mu1)
    # potential (counterfactual) study S=2 under treatment C
    data.S2$trt <- 0 
    # predict outcomes on natural scale for each individual under C
    mu0 <- predict(mod, newdata=data.S2, type="response")
    # marginal mean outcome on natural scale under C
    hat.mu0 <- mean(mu0)
    if (outcome.model=="logistic") {
      # marginal log odds ratio estimate
      delta.AC <- log((hat.mu1/(1-hat.mu1))/(hat.mu0/(1-hat.mu0)))      
    } else if (outcome.model=="linear") {
      # marginal mean difference estimate
      delta.AC <- hat.mu1-hat.mu0
    }
    return(delta.AC)  
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.S1, statistic=param.gcomp.boot, R=resamples)
  # bootstrap mean of marginal treatment effect estimate for A vs. C
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of marginal treatment effect estimate for A vs. C
  hat.var.Delta.AC <- var(boot.object$t)
  # analysis of S=2 IPD
  if (outcome.model=="logistic") {
    outcome.fit <- glm(y~trt, data=data.S2, family=binomial)
  } else if (outcome.model=="linear") {
    outcome.fit <- lm(y~trt, data=data.S2)
  }
  # marginal treatment effect estimate for B vs. C
  hat.Delta.BC <- summary(outcome.fit)$coef[2]
  # variance estimate for marginal treatment effect of B vs. C
  hat.var.Delta.BC <- vcov(outcome.fit)[2,2]
  # indirect treatment comparison (A vs. B)
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
  list(hat.Delta.AB, hat.var.Delta.AB)
} 

### Bucher method (unadjusted: does not adjust for covariates)
bucher.wrapper <- function(data.S1, data.S2, outcome.model="logistic") {
  if (outcome.model=="logistic") {
    # simple logistic regressions of outcome on treatment fitted to IPD
    outcome.fit.S1 <- glm(y~trt, data=data.S1, family=binomial)
    outcome.fit.S2 <- glm(y~trt, data=data.S2, family=binomial)
  } else if (outcome.model=="linear") {
    # simple linear regressions of outcome on treatment fitted to IPD
    outcome.fit.S1 <- lm(y~trt, data=data.S1)
    outcome.fit.S2 <- lm(y~trt, data=data.S2)
  }
  # treatment coefficient is marginal effect estimate for A vs. C
  hat.Delta.AC <- summary(outcome.fit.S1)$coef[2]
  # variance estimate for marginal effect of A vs. C
  hat.var.Delta.AC <- vcov(outcome.fit.S1)[2,2]
  # treatment coefficient is marginal effect estimate for B vs. C
  hat.Delta.BC <- summary(outcome.fit.S2)$coef[2]
  # variance estimate for marginal effect of B vs. C
  hat.var.Delta.BC <- vcov(outcome.fit.S2)[2,2]
  # indirect treatment comparison (A vs. B)
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
  list(hat.Delta.AB, hat.var.Delta.AB)  
} 

#### SIMULATED DATASETS USING LOGISTIC MODEL (MAIN TEXT SIMULATION STUDY) ####

# simulate IPD for index trial (S=1)
IPD.S1 <- replicate(n=n.sim, expr=sim.data(n=10000, mean_X=0, outcome.model="logistic"), 
                    simplify=FALSE)

# simulate IPD for competitor trial (S=2)
IPD.S2 <- replicate(n=n.sim, expr=sim.data(n=10000, mean_X=-1.4, outcome.model="logistic"),
                    simplify=FALSE)

#### ANALYSIS OF SIMULATED DATASETS ####
## logistic outcome model (simulation study in main text)

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

## run indirect treatment comparison methods for all simulation replicates

# Matching-adjusted indirect comparison
maic.results <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                        .init=list(list(), list()), .options.snow=opts,
                        .packages=c("boot")) %dopar% {
                          results <- maic.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                  resamples=1000, outcome.model="logistic")
                          return(results)
                        }
close(pb)
maic.mean.log.OR <- unlist(maic.results[[1]])
maic.var.log.OR <- unlist(maic.results[[2]])
save(maic.mean.log.OR, file=paste0("Results/MAIC/means.RData"))
save(maic.var.log.OR, file=paste0("Results/MAIC/variances.RData"))
  
# Parametric G-computation using maximum likelihood estimation
gcomp.results <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                         .init=list(list(), list()), .options.snow=opts,
                         .packages=c("boot")) %dopar% {
                           results <- param.gcomp.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                          resamples=1000, outcome.model="logistic")
                           return(results)
                        }
close(pb)
gcomp.mean.log.OR <- unlist(gcomp.results[[1]])
gcomp.var.log.OR <- unlist(gcomp.results[[2]])
save(gcomp.mean.log.OR, file=paste0("Results/GComp/means.RData"))
save(gcomp.var.log.OR, file=paste0("Results/GComp/variances.RData"))

# Bucher method (unadjusted anchored indirect treatment comparison)
bucher.results <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list()), .options.snow=opts) %dopar% {
                            results <- bucher.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                      outcome.model=="logistic")
                            return(results)
                          }
close(pb)
bucher.mean.log.OR <- unlist(bucher.results[[1]])
bucher.var.log.OR <- unlist(bucher.results[[2]])
save(bucher.mean.log.OR, file=paste0("Results/Bucher/means.RData"))
save(bucher.var.log.OR, file=paste0("Results/Bucher/variances.RData"))

### COMPUTE TRUE MARGINAL ESTIMANDS ### 
## logistic outcome model (simulation study in main text)

# true marginal estimand for active treatment vs. control in S=1
truth.S1 <- sim.data(n=10000000, mean_X=0)
p1.S1 <- mean(truth.S1[which(truth.S1$trt==1),]$y)
p0.S1 <- mean(truth.S1[which(truth.S1$trt==0),]$y)
# true marginal log OR in S=1 for active vs. control
true.log.OR.S1 <- log((p1.S1/(1-p1.S1))/(p0.S1/(1-p0.S1))) 

# true marginal estimand for active treatment vs. control in S=2
truth.S2 <- sim.data(n=10000000, mean_X=-1.4)
p1.S2 <- mean(truth.S2[which(truth.S2$trt==1),]$y)
p0.S2 <- mean(truth.S2[which(truth.S2$trt==0),]$y)
# true marginal log OR in S=2 for active vs. control
true.log.OR.S2 <- log((p1.S2/(1-p1.S2))/(p0.S2/(1-p0.S2))) 

# true marginal log odds ratio for A vs. B in S=2 is zero
Delta.AB <- 0 

### PROCESS RESULTS, COMPUTE AND PLOT SIMULATION STUDY METRICS ###
## logistic outcome model (simulation study in main text)

source("functions.R") # load functions to compute performance measures

# to store performance measures
simulation.metrics <- as.data.frame(matrix(nrow=3, ncol=11)) # 3 rows: 3 methods
colnames(simulation.metrics) <- c("Method", "Bias", "Bias.MCSE", "LCI", "LCI.MCSE", 
                                  "UCI", "UCI.MCSE", "Cov", "Cov.MCSE", "MSE", "MSE.MCSE")

# to store all A vs. B marginal log odds ratio point estimates
ate.table <- as.data.frame(matrix(nrow=3*n.sim, ncol=2))
colnames(ate.table) <- c("Method", "ATE")

# function that computes performance metrics for a given method
process.metrics <- function(means, variances, truth) {
  bias.metric <- bias(means, truth)
  bias.metric.mcse <- bias.mcse(means)
  # construct Wald-type interval estimates using normal distribution
  lci <- means + qnorm(0.025)*sqrt(variances)
  uci <- means + qnorm(0.975)*sqrt(variances)
  lci.mean <- mean(lci)
  lci.mcse <- mcse.estimate(lci)
  uci.mean <- mean(uci)
  uci.mcse <- mcse.estimate(uci)
  cov <- coverage(lci, uci, truth)
  cov.mcse <- coverage.mcse(cov, length(means))
  mse.metric <- mse(means, truth) 
  mse.metric.mcse <- mse.mcse(means, truth) 
  list(bias.metric, bias.metric.mcse, lci.mean, lci.mcse, uci.mean, 
       uci.mcse, cov, cov.mcse, mse.metric, mse.metric.mcse)
} 

### Matching-adjusted indirect comparison (MAIC)
load(paste0("Results/MAIC/means.RData"))
load(paste0("Results/MAIC/variances.RData"))
simulation.metrics[1,1] <- "MAIC"
maic.metrics <- process.metrics(maic.mean.log.OR, maic.var.log.OR, Delta.AB)
simulation.metrics[1,2:11] <- unlist(maic.metrics)
ate.table[1:n.sim,1] <- "MAIC"
ate.table[1:n.sim,2] <- maic.mean.log.OR

### Parametric G-computation
load(paste0("Results/GComp/means.RData"))
load(paste0("Results/GComp/variances.RData"))
simulation.metrics[2,1] <- "G-computation"
gcomp.metrics <- process.metrics(gcomp.mean.log.OR, gcomp.var.log.OR, Delta.AB)
simulation.metrics[2,2:11] <- unlist(gcomp.metrics)
ate.table[(n.sim+1):(n.sim*2),1] <- "G-computation"
ate.table[(n.sim+1):(n.sim*2),2] <- gcomp.mean.log.OR

### Bucher method (indirect comparison without covariate adjustment)
load(paste0("Results/Bucher/means.RData"))
load(paste0("Results/Bucher/variances.RData"))
simulation.metrics[3,1] <- "Bucher"
bucher.metrics <- process.metrics(bucher.mean.log.OR, bucher.var.log.OR, Delta.AB)
simulation.metrics[3,2:11] <- unlist(bucher.metrics)
ate.table[((n.sim*2)+1):(n.sim*3),1] <- "Bucher"
ate.table[((n.sim*2)+1):(n.sim*3),2] <- bucher.mean.log.OR

## Function generates ridgeline plot for point estimates
plot.results <- function() {
  ridge.plot <- ggplot(ate.table, aes(x=ATE, y=Method, fill=Method)) +
    geom_density_ridges(alpha=0.65) +
    geom_vline(xintercept=0, linetype="dashed", color ="red") +
    scale_x_continuous(limits=c(-1, 1)) +
    scale_y_discrete(limits=c("MAIC", "G-computation", "Bucher")) +
    theme_classic() + 
    theme(legend.position = "none", 
          axis.text.y = element_text(color="grey20", size=9, face ="plain"),
          axis.text.x = element_text(color="grey20", size=7, face="plain"),
          plot.title = element_text(color="grey20", size=10, face ="plain"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(color="grey20", size=9, face="plain")) +
    scale_fill_brewer(palette="Dark2") + xlab("Point estimates")
  # a ridgeline plot with the spread of point estimates is returned
  return(ridge.plot)
}  

## Function generates results table with performance metrics 
table.results <- function() {
  display.table <- cbind(Method=simulation.metrics$Method,
                         Bias=paste0(format(round(simulation.metrics$Bias,digits=3),nsmall=3)," (",
                                     format(round(simulation.metrics$Bias.MCSE,digits=3),nsmall=3),")"),
                         Coverage=paste0(format(round(simulation.metrics$Cov,digits=3),nsmall=3)," (",
                                    format(round(simulation.metrics$Cov.MCSE,digits=3),nsmall=3),")"),
                         MSE=paste0(format(round(simulation.metrics$MSE,digits=3),nsmall=3)," (",
                                    format(round(simulation.metrics$MSE.MCSE,digits=3),nsmall=3),")"))
  table.grob <- tableGrob(display.table, theme=ttheme_minimal(base_size=7))
  # a table with the performance measures and corresponding MCSEs is returned
  return(table.grob)
}

ridge.plot <- plot.results() # ridgeline plot
table.grob <- table.results() # table of results

ridge.grid <- arrangeGrob(ridge.plot, table.grob, ncol=2, widths=c(0.8,1.2))

ggsave(file="Figure1.pdf", plot=ridge.grid, width=170, height=45, units="mm", dpi = 300)

#### SIMULATED DATASETS USING LINEAR MODEL (APPENDIX A SIMULATION STUDY) ####

# simulate IPD for index trial (S=1)
IPD.S1 <- replicate(n=n.sim, expr=sim.data(n=10000, mean_X=0, outcome.model="linear"), 
                    simplify=FALSE)

# simulate IPD for competitor trial (S=2)
IPD.S2 <- replicate(n=n.sim, expr=sim.data(n=10000, mean_X=-1.4, outcome.model="linear"),
                    simplify=FALSE)

#### ANALYSIS OF SIMULATED DATASETS ####
## linear outcome model (simulation study in Appendix A)

## run indirect treatment comparison methods for all simulation replicates

# Matching-adjusted indirect comparison
maic.results.lr <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                           .init=list(list(), list()), .options.snow=opts,
                           .packages=c("boot")) %dopar% {
                             results <- maic.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                     resamples=1000, outcome.model="linear")
                             return(results)
                           }
close(pb)
maic.mean.MD <- unlist(maic.results.lr[[1]])
maic.var.MD <- unlist(maic.results.lr[[2]])
save(maic.mean.MD, file=paste0("Results/MAIC/means_lr.RData"))
save(maic.var.MD, file=paste0("Results/MAIC/variances_lr.RData"))

# Parametric G-computation using maximum likelihood estimation
gcomp.results.lr <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts,
                            .packages=c("boot")) %dopar% {
                              results <- param.gcomp.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                             resamples=1000, outcome.model="linear")
                              return(results)
                            }
close(pb)
gcomp.mean.MD <- unlist(gcomp.results.lr[[1]])
gcomp.var.MD <- unlist(gcomp.results.lr[[2]])
save(gcomp.mean.MD, file=paste0("Results/GComp/means_lr.RData"))
save(gcomp.var.MD, file=paste0("Results/GComp/variances_lr.RData"))

# Bucher method (unadjusted anchored indirect treatment comparison)
bucher.results.lr <- foreach(i=1:n.sim, .combine='comb', .multicombine=TRUE,
                             .init=list(list(), list()), .options.snow=opts) %dopar% {
                               results <- bucher.wrapper(data.S1=IPD.S1[[i]], data.S2=IPD.S2[[i]],
                                                         outcome.model="linear")
                               return(results)
                             }
close(pb)
bucher.mean.MD <- unlist(bucher.results.lr[[1]])
bucher.var.MD <- unlist(bucher.results.lr[[2]])
save(bucher.mean.MD, file=paste0("Results/Bucher/means_lr.RData"))
save(bucher.var.MD, file=paste0("Results/Bucher/variances_lr.RData"))

stopCluster(cluster)

### PROCESS RESULTS, COMPUTE AND PLOT SIMULATION STUDY METRICS ###
## linear outcome model (simulation study in Appendix A)

# true marginal log relative risk for A vs. B in S=2 is zero
Delta.AB <- 0 

# to store performance measures
simulation.metrics <- as.data.frame(matrix(nrow=3, ncol=11)) # 3 rows: 3 methods
colnames(simulation.metrics) <- c("Method", "Bias", "Bias.MCSE", "LCI", "LCI.MCSE", 
                                     "UCI", "UCI.MCSE", "Cov", "Cov.MCSE", "MSE", "MSE.MCSE")

# to store all A vs. B marginal log relative risk point estimates
ate.table <- as.data.frame(matrix(nrow=3*n.sim, ncol=2))
colnames(ate.table) <- c("Method", "ATE")

### Matching-adjusted indirect comparison (MAIC)
load(paste0("Results/MAIC/means_lr.RData"))
load(paste0("Results/MAIC/variances_lr.RData"))
simulation.metrics[1,1] <- "MAIC"
maic.metrics.lr <- process.metrics(maic.mean.MD, maic.var.MD, Delta.AB)
simulation.metrics[1,2:11] <- unlist(maic.metrics.lr)
ate.table[1:n.sim,1] <- "MAIC"
ate.table[1:n.sim,2] <- maic.mean.MD

### Parametric G-computation
load(paste0("Results/GComp/means_lr.RData"))
load(paste0("Results/GComp/variances_lr.RData"))
simulation.metrics[2,1] <- "G-computation"
gcomp.metrics.lr <- process.metrics(gcomp.mean.MD, gcomp.var.MD, Delta.AB)
simulation.metrics[2,2:11] <- unlist(gcomp.metrics.lr)
ate.table[(n.sim+1):(n.sim*2),1] <- "G-computation"
ate.table[(n.sim+1):(n.sim*2),2] <- gcomp.mean.MD

### Bucher method (indirect comparison without covariate adjustment)
load(paste0("Results/Bucher/means_lr.RData"))
load(paste0("Results/Bucher/variances_lr.RData"))
simulation.metrics[3,1] <- "Bucher"
bucher.metrics.lr <- process.metrics(bucher.mean.MD, bucher.var.MD, Delta.AB)
simulation.metrics[3,2:11] <- unlist(bucher.metrics.lr)
ate.table[((n.sim*2)+1):(n.sim*3),1] <- "Bucher"
ate.table[((n.sim*2)+1):(n.sim*3),2] <- bucher.mean.MD

ridge.plot.lr <- plot.results() # ridgeline plot
table.grob.lr <- table.results() # table of results

ridge.grid.lr <- arrangeGrob(ridge.plot.lr, table.grob.lr, ncol=2, widths=c(0.8,1.2))

ggsave(file="Figure2.pdf", plot=ridge.grid.lr, width=170, height=45, units="mm", dpi = 300)
