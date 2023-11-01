## This file processes the results of the simulation study and computes and graphs the relevant performance metrics

rm(list=ls()) # clear

# Set working directory
# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/conditional_marginal_effect_modifiers")
source("functions.R") # load functions to compute performance measures
load(file="sim_settings.RData") # load simulation study settings

# load packages required to plot the simulation study results
library("ggplot2")
library("ggridges")
library("gridExtra")

# simulation study scenario settings
n.sim <- 500 # number of Monte Carlo replicates (simulated datasets) per scenario/method
settings <- 1:3
outcome.models <- c("linear", "log-linear", "logistic")
# setting/model combinations for each simulation scenario
scenarios <- expand.grid(setting=settings, OM=outcome.models)
n.scenarios <- nrow(scenarios) # number of scenarios

Delta.AB <- 0 # true value of marginal A vs. B treatment effect in S=2 is zero

# data frame to store performance measures (3 methods, each row repeated 3 times)
simulation.metrics <- scenarios[rep(seq_len(n.scenarios), each=3), ]
metrics.names <- c("Method", "Bias", "Bias.MCSE", "LCI", "LCI.MCSE", 
                  "UCI", "UCI.MCSE", "Cov", "Cov.MCSE", "MSE", "MSE.MCSE")
simulation.metrics[metrics.names] <- NA

# data frame to store all A vs. B marginal treatment effect point estimates
ate.table <- scenarios[rep(seq_len(n.scenarios), each=3*n.sim),]
ate.table["Method"] <- NA
ate.table["ATE"] <- NA
  
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

j <- 1 # row counter for simulation metrics
k <- 1 # row counter for ATEs

for (i in 1:n.scenarios) {
  file.id <- paste0("setting", scenarios$setting[i], "OM", scenarios$OM[i])   
  # Bucher method (indirect treatment comparison without covariate adjustment)
  load(paste0("Results/Bucher/means_", file.id, ".RData"))
  load(paste0("Results/Bucher/variances_", file.id, ".RData"))  
  simulation.metrics[j,3] <- "Bucher"
  bucher.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:13] <- unlist(bucher.metrics)
  ate.table[k:(k+n.sim-1),3] <- "Bucher"
  ate.table[k:(k+n.sim-1),4] <- means
  j <- j+1
  k <- k+n.sim   
  # parametric G-computation
  load(paste0("Results/GComp/means_", file.id, ".RData"))
  load(paste0("Results/GComp/variances_", file.id, ".RData"))  
  simulation.metrics[j,3] <- "G-computation"
  gcomp.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:13] <- unlist(gcomp.metrics)
  ate.table[k:(k+n.sim-1),3] <- "G-computation"
  ate.table[k:(k+n.sim-1),4] <- means
  j <- j+1
  k <- k+n.sim  
  # matching-adjusted indirect comparison (MAIC)
  load(paste0("Results/MAIC/means_", file.id, ".RData"))
  load(paste0("Results/MAIC/variances_", file.id, ".RData"))
  simulation.metrics[j,3] <- "MAIC"
  maic.metrics <- process.metrics(means, variances, Delta.AB) 
  simulation.metrics[j,4:13] <- unlist(maic.metrics)
  ate.table[k:(k+n.sim-1),3] <- "MAIC"
  ate.table[k:(k+n.sim-1),4] <- means
  j <- j+1
  k <- k+n.sim
}

# function generates ridgeline plot of point estimates for a specific scenario
plot.results <- function(scenario) {
  i <- scenario
  OM.ates <- ate.table[ate.table$OM==scenarios$OM[i], ]
  scenario.ates <- subset(OM.ates, setting==scenarios$setting[i])
  ridge.plot <- ggplot(scenario.ates, aes(x=ATE, y=Method, fill=Method)) +
    geom_density_ridges(alpha=0.65) +
    geom_vline(xintercept=0, linetype="dashed", color ="red") +
    scale_x_continuous(limits=c(-1.2, 1.2)) +
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

# function generates results table with performance metrics for a specific scenario
table.results <- function(scenario) {
  i <- scenario
  OM.metrics <- simulation.metrics[simulation.metrics$OM==scenarios$OM[i], ]
  scenario.metrics <- subset(OM.metrics, setting==scenarios$setting[i])
  display.table <- cbind(Method=scenario.metrics$Method,
                         Bias=paste0(format(round(scenario.metrics$Bias,digits=3),nsmall=3)," (",
                                     format(round(scenario.metrics$Bias.MCSE,digits=3),nsmall=3),")"),
                         Coverage=paste0(format(round(scenario.metrics$Cov,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$Cov.MCSE,digits=3),nsmall=3),")"),
                         MSE=paste0(format(round(scenario.metrics$MSE,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$MSE.MCSE,digits=3),nsmall=3),")"))
  table.grob <- tableGrob(display.table, theme=ttheme_minimal(base_size=7))
  # a table with the performance measures and corresponding MCSEs is returned
  return(table.grob)
}

# ridgeline plot for each simulation scenario

## continuous outcome, linear outcome model and mean difference as summary effect measure

# setting 1: treatment effect homogeneity: imbalanced means and uncorrelated covariates
ridge.plot.S1.MD <- plot.results(scenario=1) + ggtitle(expression("Treatment effect homogeneity: imbalanced means, uncorrelated covariates"))

# setting 2: treatment effect heterogeneity: balanced means and differences in correlation structures
ridge.plot.S2.MD <- plot.results(scenario=2) + ggtitle(expression("Treatment effect heterogeneity: balanced means, different correlation structures"))

# setting 3: treatment effect heterogeneity: imbalanced means and differences in correlation structures
ridge.plot.S3.MD <- plot.results(scenario=3) + ggtitle(expression("Treatment effect heterogeneity: imbalanced means, different correlation structures"))

## count outcome, log-linear outcome model and log risk ratio as summary effect measure

# setting 1: treatment effect homogeneity: imbalanced means and uncorrelated covariates
ridge.plot.S1.log.RR <- plot.results(scenario=4) + ggtitle(expression("Treatment effect homogeneity: imbalanced means, uncorrelated covariates"))
                                                                      
# setting 2: treatment effect heterogeneity: balanced means and differences in correlation structures
ridge.plot.S2.log.RR <- plot.results(scenario=5) + ggtitle(expression("Treatment effect heterogeneity: balanced means, different correlation structures"))
  
# setting 3: treatment effect heterogeneity: imbalanced means and differences in correlation structures
ridge.plot.S3.log.RR <- plot.results(scenario=6) + ggtitle(expression("Treatment effect heterogeneity: imbalanced means, different correlation structures"))

## binary outcome, logistic outcome model and log odds ratio as summary effect measure

ridge.plot.S1.log.OR <- plot.results(scenario=7) + ggtitle(expression("Treatment effect homogeneity: imbalanced means, uncorrelated covariates"))
  
ridge.plot.S2.log.OR <- plot.results(scenario=8) + ggtitle(expression("Treatment effect heterogeneity: balanced means, different correlation structures"))

ridge.plot.S3.log.OR <- plot.results(scenario=9) + ggtitle(expression("Treatment effect heterogeneity: imbalanced means, different correlation structures"))

# table of results for each scenario
table.grob.S1.MD <- table.results(scenario=1)  
table.grob.S2.MD <- table.results(scenario=2)  
table.grob.S3.MD <- table.results(scenario=3)  
table.grob.S1.log.RR <- table.results(scenario=4)  
table.grob.S2.log.RR <- table.results(scenario=5)  
table.grob.S3.log.RR <- table.results(scenario=6)  
table.grob.S1.log.OR <- table.results(scenario=7)  
table.grob.S2.log.OR <- table.results(scenario=8)  
table.grob.S3.log.OR <- table.results(scenario=9)  

ridge.grid.MD <- arrangeGrob(ridge.plot.S1.MD, table.grob.S1.MD, ridge.plot.S2.MD, table.grob.S2.MD,
                             ridge.plot.S3.MD, table.grob.S3.MD, ncol=2, widths=c(0.8,1.2)) 
ggsave(file="Results/Figure1.pdf", plot=ridge.grid.MD, width=170, height=135, units="mm", dpi = 300)

ridge.grid.log.RR <- arrangeGrob(ridge.plot.S1.log.RR, table.grob.S1.log.RR, ridge.plot.S2.log.RR, table.grob.S2.log.RR,
                                 ridge.plot.S3.log.RR, table.grob.S3.log.RR, ncol=2, widths=c(0.8,1.2)) 
ggsave(file="Results/Figure2.pdf", plot=ridge.grid.log.RR, width=170, height=135, units="mm", dpi = 300)

ridge.grid.log.OR <- arrangeGrob(ridge.plot.S1.log.OR, table.grob.S1.log.OR, ridge.plot.S2.log.OR, table.grob.S2.log.OR,
                                 ridge.plot.S3.log.OR, table.grob.S3.log.OR, ncol=2, widths=c(0.8,1.2)) 
ggsave(file="Results/Figure3.pdf", plot=ridge.grid.log.OR, width=170, height=135, units="mm", dpi = 300)


# Save simulation study performance metrics
write.csv(simulation.metrics, "Results/performance_metrics.csv", row.names = FALSE)
