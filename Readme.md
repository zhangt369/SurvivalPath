## Brief Introduction
The survival path mapping approach has been proposed for dynamic prognostication of cancer 
patients using time-series survival data. The SurvivalPath R package was developed to facilitate 
building personalized survival path models. The package contains
functions to convert time-series data by into time-slices data by fixed interval based on
time information of input medical records. After the pre-processing of data, under a userdefined 
parameters on covariates, significance level, minimum bifurcation sample size and
number of time slices for analysis, survival paths can be computed using the main function,
which can be visualized as a tree diagram, with important parameters annotated. The
package also includes function for analyzing the connections between exposure/treatment
and node transitions, and function for screening patient subgroup with specific features,
which can be used for further exploration analysis. The SurvivalPath R package is freely
available from github, https://github.com/zhangt369/SurvivalPath.

**Installation**

- 1.Download and install the package “SurvivalPath” in the CRAN. 
- 2.Download the latest release (tar.gz or zip file) at “https://github.com/zhangt369/SurvivalPath” and install the package through local path through the command:
    install.packages("~/Yourpath/filename", repos = NULL, type = "source")).
