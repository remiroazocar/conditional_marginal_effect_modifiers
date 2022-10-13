# Purely prognostic variables may modify marginal treatment effects for non-collapsible effect measures: Code

### Antonio Remiro-Azócar

### *remiroantonio@gmail.com*

### *2022*

This repository contains the R code used for my paper [Purely prognostic variables may modify marginal treatment effects for non-collapsible effect measures][1]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location where the download was made, and open and run the `main.R` script. This script specifies the settings of the simulation study, generates the data, performs the indirect treatment comparison methods (saving the point estimates and variances to the `"./Results/"` subdirectory), computes the relevant performance metrics, and graphs the results of the simulation study. The `functions.R` script contains user-defined functions to evaluate the performance measures of interest. 

In the simulation study, the `doSNOW` package is used to parallelize the performance of the methods, distributing the tasks to different cores of the computer. 

The code was prepared in `RStudio` using `R` version `4.1.1` in a Windows architecture, with a 64-bit operating system. The following packages and versions were used:

* `boot 1.3.28` required for use of the non-parametric bootstrap in covariate-adjusted indirect comparisons 
* `doSNOW 1.0.19` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `ggplot2 3.3.5` to plot the simulation study results (Figure 1 in the article)
* `ggridges 0.5.3` to plot the simulation study results (Figure 1 in the article)
* `gridExtra 2.3` to plot the simulation study results (Figure 1 in the article)
* `parallel 4.1.1` to detect the number of CPU cores

[1]: https://doi.org/10.48550/arXiv.2210.01757
