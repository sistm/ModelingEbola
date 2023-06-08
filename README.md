# ODEmodeling_AbPrediction_Ebola
R and mlxtran codes for generating results about the evaluation and prediction of the long-term antibody response induced by Ad26/MVA vaccine against Ebola.


In this project, we were interested in evaluating the persistence of the long-term immune response induced by the two-dose heterologous Ad26.ZEBOV, MVA-BN-Filo vaccine regimen against Ebola.
We used a mechanistic model based on ordinary differential equations describing the dynamics of binding antibodies and short- and long-lived antibody-secreting cells (ASCs) to model the humoral response from 7 days after the 2nd vaccination. Initially estimated on data from Phase I clinical trials with a follow-up period of 1 year (Pasin et al., JVI 2019), first we assessed the robustness of the model using Phase II data with extended follow-up of 2 years. Then, we updated model parameter estimated using all Phase I and Phase II data.  

In this folder, we gathered R and Mlxtran codes used to perform this work.


## Content of folders
**R codes:**
* *ProfileLikelihood.R*: Code used to generate profile likelihood for the fixed parameter defining the decay rate of long-lived ASCs ($\delta_L$)
* *Model_Averaging.R*: Code used to perform a model averaging approach on our model to integrate model uncertainty in the value of $\delta_L$ in the calculation of the confidence intervals of model parameter estimates.
* *MCCV.R*: Code used to perform Monte Carlo cross-validation (MCCV) to evaluate the ability of the model to predict unseen data.

**Mlxtran codes [code used by Monolix software]:**
* *MlxtranCode_ODEmodel_validationPrediction.txt*: Code describing the mechanistic model used in the validation step of the work (model developed by Pasin et al.)
* *MlxtranProject_validationPrediction.txt*: Extraction of the .mlxtran file obtained for the Monolix project in the validation step. This file describe parameter transformation and covariate model.
* *MlxtranCode_ODEmodel_ReEstimation.txt*: Code describing the second version the mechanisitic model integrating an adjustment for laboratory effects in the observation model. Model used in the re-estimation step of the work.
*  *MlxtranProject_ReEstimation_FinalModel.txt*: Extraction of the .mlxtran file obtained for the Monolix project in the re-estimation step. This file describe parameter transformation and covariate model.
*  *MlxranModel_Individual_Simulation.txt*: Code used by Simulx in the MCCV to simulate individual model parameters and calculate the resulting individual antibody dynamics.
