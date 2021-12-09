# Possibly non-invertible SVARMA Models: The Normalised Canonical WHF Parametrisation 

This repository contains the R-package and scripts for the empirical application in the associated article.
The directory [svarmawhf](./svarmawhf) contains the R-package. 
The file [Package Structure](./scriptsP_whf_revision_sgt_bq/z_package_structure.html) (or the [pdf version](./scriptsP_whf_revision_sgt_bq/z_package_structure.pdf))is a good starting point to get familiar with the important functions and the call structure.
All functions are documented. 
The documentation can be accessed as usual with `help(<function_name>)` or `?<function_name>`.
The package builds on the R-packages **RLDM** and **rationalmatrices**, authored jointly with Wolfgang Scherrer.
Since these packages might change, the parts which are necessary for the analysis in the associated article are extracted to R files whose names start with **~/svarmawhf/R/zz_ratmat_** and  **~/svarmawhf/R/zz_rldm_**. 

## Installation

By opening the file **~/svarmawhf.Rproj**, the working directory of R is set to the project root. 
This is necessary that all scripts run as intended and that the package compiles.
In R-Studio, use the shortcuts *cmd+shift+d* (or **devtools::document()**) to generate the documentation, and *cmd+shift+b* (or **devtools::build()**) to build the package, see [http://r-pkgs.had.co.nz/].

# Steps in the Analysis

All files which are relevant for the analysis of the empirical application are contained in the directory [scriptsP_whf_revision_sgt_bq](./scriptsP_whf_revision_sgt_bq).
In the following, we will describe the content of the files in this directory.
Since all results are also uploaded into the [folder](./local_data), it is not necessary to run any script.

The analysis is separated into different steps.
Each file starts with an overview of the main conclusion, i.e. it is not necessary to go through the files in detail if one is only interested in the main take-aways.

## Data Preparation

The Rmarkdown file [Data Preparation](./scriptsP_whf_revision_sgt_bq/z_data_preparation.html) (or the [pdf version](./scriptsP_whf_revision_sgt_bq/z_data_preparation.pdf)) loads and transforms the data of [Blanchard and Quah](./local_data/g_gmr_bq/bqdata.csv) in the same way as this is done in GMR.
Moreover, visualisations of intermediary data transformation steps are shown and described.

Eventually, the data are saved and serve as input to the main script which we will describe next.

## Evaluating all Combinations of Integer-Valued Parameters

The [main script](./scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R) performs the main work of estimating the normalised canonical WHF model for all combinations of $\(p,q,\kappa, k\)$.
This script is run on the HPC of the University of Helsinki. 
The [equivalent script](./scriptsP_whf_revision_sgt_bq/ab_Rscript_local_bq.R) produces the same results but is intended for local use. 
It can be adjusted to use MS Azure as parallel backend with the R-package [doAzureparallel](https://github.com/Azure/doAzureParallel).

Eventually, the results are saved as rds-files and later extracted. 
For the particular steps taken, see the documentation package function [create_results()](./svarmawhf/R/ab_create_results.R).
In essence, we create a data frame of all combinations of integer-valued parameters, separate the data frame for parallelisation and perform the optimisation steps for each set of integer-valued parameters.

### SLURM script

The [SLURM script](./scriptsP_whf_revision_sgt_bq/b_slurm) is uploaded to the HPC cluster of the University of Helsinki which in turn calls the [main script](./scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R).

Since evaluation of our model for different integer-valued parameters is embarrassingly parallel, we use array-jobs to perform the calculations.
In our case, the evaluation of all models takes about 15 minutes.

## Analysis of Errors in Optimization

The Rmarkdown file [Error Analysis](./scriptsP_whf_revision_sgt_bq/c_error_analysis.html) analyses the convergence properties of the Nelder-Mead and BFGS optimizations for the Gaussian, Laplace, SGT densities respectively.

The main take aways are that the optimisation works well except for rather high MA orders (larger than 5).

## Model Selection

The Rmarkdown file [Model Selection](./scriptsP_whf_revision_sgt_bq/d_modelselection_aic_bic.html) generates AIC and BIC values for all integer-valued parameters and plots them.

This is a preliminary step for final model selection which is based on both AIC/BIC values and the independence properties of the residuals.

## Residual Checks

Since the residuals should be independent and non-Gaussian, we perform the Shapiro-Wilk and Jarque-Bera tests as wells as the Ljung-Box test in the Rmarkdown file [Residual Checks](./scriptsP_whf_revision_sgt_bq/e_residualcheck.html).

Together with the AIC/BIC values, we decide on the model with integer-valued parameters $(p,q,\kappa, k) = (1,2,1,0)$.

## Impulse Response Functions

In the Rmarkdown file [IRFs](./scriptsP_whf_revision_sgt_bq/f_irfs.html), we obtain the IRFs for some of the best models (with respect to AIC/BIC/Shapiro-Wilk/Jarque-Bera/Ljung-Box).
We use either estimated values for the static shock transmission matrix $B$ or impose the long-run restriction suggested by Blanchard and Quah (1989).

In addition, we generate the IRFs of Blanchard and Quah (1989) and GMR.

## Robust Standard Errors

In the Rmarkdown file [Robust Standard Errors](./scriptsP_whf_revision_sgt_bq/h_bestmodel_robust.html), we obtain different estimates for the standard deviations of the system and noise parameters.
In particular, we compare estimates under the assumption of correct model specification with robust ones.

Under correct specification, the information matrix matrix may be calculated as the 

* outer product gradient (OPG) of the scores
* the Hessian obtained as output for `stats::optim()`
* the analytic version of the Hessian obtained from `numDeriv::hessian()`

Allowing for misspecification, the robust standard errors are obtained from the covariance matrix $\Omega^{-1}\left(\hat{\theta_T}\right) I\left(\hat{\theta_T}\right) \Omega^{-1}\left(\hat{\theta_T}\right)$.
Again, the Hessian may be calculated as output of `stats::optim()` or in an analytic way from `numDeriv::hessian()`.

First, we investigate the OPG and the different versions of the Hessians quantitatively and by plotting the values.
Subsequently, we compare the diagonal elements of the OPG and the different versions of the Hessian.
Next, we evaluate the implied standard deviations for the system and noise parameters quantitatively and with plots.
Finally, we provide the results, i.e. each part of the parameter estimates together with their robust and non-robust standard errors.

# Specials

In order to investigate the [SGT distributions](./flexdashboards/sgt_dashboard.Rmd) and [Gaussian mixtures](./flexdashboards/mixtures_dashboard.Rmd), we have created two dashboards in the subdirectory **flexdashboards**.
They may also be directly accessed via the links [https://funber.shinyapps.io/sgt_dashboard/](https://funber.shinyapps.io/sgt_dashboard/) and [https://funber.shinyapps.io/mixtures_dashboard/](https://funber.shinyapps.io/mixtures_dashboard/)

