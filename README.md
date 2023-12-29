# Transportability of model-based estimands in evidence synthesis: Code

### Antonio Remiro-Azócar

### *remiroantonio@gmail.com*

### *2023*

This repository contains the R code used for my paper [Transportability of model-based estimands in evidence synthesis][1]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location where the download was made.  To run the pipeline, the user should then open and run the scripts in the following order. 

|       Script        | Explanation                                                  |
| :-----------------: | ------------------------------------------------------------ |
|    `gen_data.R`     | Specifies the settings of the simulation study and saves them in `"./sim_settings.RData"`. Generates the data for the simulation study according to the settings, saving the data to the `"./Data/"` subdirectory. |
|  `anchored_ITCs.R`  | Performs the anchored indirect comparison methods on the simulated data (saving the point estimates and variances to the `"./Results/"` subdirectory) |
| `process_results.R` | Processes the results of the simulation study and computes and graphs the relevant performance metrics (saved to the `"./Results/"` subdirectory) |

The `functions.R` script contains user-defined functions to evaluate the performance measures of interest. The file `./Results/performance_metrics.csv` has recorded the key performance measures associated with each simulation scenario/setting, as presented in Figure 1, Figure 2 and Figure 3 of the paper. 

In the simulation study, the `doSNOW` package is used to parallelize the performance of the anchored indirect comparison methods, distributing the tasks to different cores of the computer. 

The code was prepared in `RStudio` using `R` version `4.1.1` in a Windows architecture, with a 64-bit operating system. The following packages and versions were used:

* `boot 1.3.28` required for use of the non-parametric bootstrap in covariate-adjusted indirect comparisons 
* `copula 1.0.1` to generate individual-level covariates in the covariate simulation step of parametric G-computation, drawing from a multivariate Gaussian copula   
* `doSNOW 1.0.19` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `dplyr 1.0.7` for data manipulation when aggregating the subject-level data for the target study
* `ggplot2 3.3.5` to plot the simulation study results (Figure 1, Figure 2 and Figure 3 in the article)
* `ggridges 0.5.3` to plot the simulation study results (Figure 1, Figure 2 and Figure 3 in the article)
* `gridExtra 2.3` to plot the simulation study results (Figure 1, Figure 2 and Figure 3 in the article)
* `MASS 7.3-55` to simulate covariates in some data-generating processes of the simulation study, drawing from a multivariate normal distribution
* `parallel 4.1.1` to detect the number of CPU cores

[1]: https://doi.org/10.48550/arXiv.2210.01757
