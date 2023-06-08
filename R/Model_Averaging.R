# ------------------------- #
# OBJECTIVE: To Perform model averaging approach on model fixed parameters to account for their uncertainty on the estimation of other model parameters
# Author: Marie Alexandre
# Date: 2023/06/08

# Related paper: Alexandre et al. (2023) - "Evaluation and prediction of the long-term humoral immune response induced by the two-dose heterologous Ad26.ZEBOV,MVA-BN-Filo vaccine regimen against Ebola"

# R version: 4.2.1
# Monolix version: >= 2019R1


# Description: In this file, we want to perform a model averaging approach on model parameter estimations in order to take into account uncertainty related to fixed parameters (delta L here)
#              The code is splitted into two parts: 
#                     1) Estimation of candidate models involved in model averaging
#                     2) Averaging of model parameters estimated for the different candidate models
# ------------------------- #




rm(list=ls())




# --- LIBRARIES --- ####
# - Packages for Monolix
library(lixoftConnectors)
# - Package for parallel calculation
library(parallel) ; library(snow) ; library(doSNOW)
# - Package for plots
library(ggplot2)
# ---------------- #


# --- FUNCTIONS --- ####
# - PART 1
Profile_likelihood_function <- function(Init_project,result_folder,parameter,tested_values,final_project_names, Nb_Cores,monolix_path){
  require(parallel) ; require(snow) ; require(doSNOW) ; require(foreach)  # packages for parallel calculation
  
  N <- length(tested_values)
  
  # -- Arguments parallel calculation
  NCores <- min(detectCores()-1,Nb_Cores)
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  pb <- txtProgressBar(max=N,style=3)         # We add a progress bar
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress=progress)
  print("cluster")
  clusterExport(Cluster,c("Init_project","result_folder","parameter","tested_values","final_project_names","monolix_path"))
  
  print(paste("Parallel calculation - Profile likelihood ",parameter,sep=""))
  
  # n <- 1
  Results <- foreach(n=seq(1,N),.combine = "rbind",.errorhandling = "remove",.packages = c("doParallel","foreach","doSNOW","lixoftConnectors"),.options.snow = opts)%dopar%{
    initializeLixoftConnectors(software = "monolix", path=monolix_path)
    
    # load of initial project 
    loadProject(Init_project)
    
    # Modification of the value of the tested parameter
    Initial_Pop_params <- getPopulationParameterInformation() 
    Initial_Pop_params$initialValue[which(Initial_Pop_params$name == paste(parameter,"pop",sep="_"))] <- tested_values[n]
    Initial_Pop_params$method[which(Initial_Pop_params$name == paste(parameter,"pop",sep="_"))] <- "FIXED"
    setPopulationParameterInformation(Initial_Pop_params)
    
    # Remove "plots" from the scenario
    scenario <- getScenario()
    scenario$tasks["plots"] <- FALSE
    setScenario(scenario)
    
    # Save project
    New_project_name <- paste(result_folder,final_project_names[n],sep="/")
    saveProject(projectFile = New_project_name)
    
    # Run the project
    runScenario()
    saveProject(projectFile = New_project_name)
    
    # Extraction of LL value 
    LL <- as.numeric(-0.5*getEstimatedLogLikelihood()$importanceSampling['-2LL'])
    res <- data.frame(param=tested_values[n],LL=LL)
    return(res)
  }
  close(pb)
  stopCluster(Cluster)
  
  return(Results)
}

# - PART 2
Model_averaging_weights_FUNCTION <- function(values){
  # Function of weight calculation that is usually applied for AIC or BIC values
  min_value <- min(values,na.rm=TRUE)
  delta_value <- values - min_value
  weights <- sapply(seq(1,length(values)),function(i) exp(-delta_value[i]/2)/sum(exp(-delta_value/2)))
  return(weights)
}
Parameter_distribution_TABLE <- function(parameter_distributions){
  
  table_results <- NULL
  
  # > a. Distribution of the decay rate of antibodies ####
  # delta Ab in Women (ref)
  deltaAb_Mean_Women <- parameter_distributions$mean[which(parameter_distributions$parameter == "delta_Ab_pop")]
  
  deltaAb_Sd_Women <- parameter_distributions$sd[which(parameter_distributions$parameter == "delta_Ab_pop")]
  deltaAb_CImin_Women <- exp(log(deltaAb_Mean_Women) - 1.96*deltaAb_Sd_Women/deltaAb_Mean_Women)
  deltaAb_CImax_Women <- exp(log(deltaAb_Mean_Women) + 1.96*deltaAb_Sd_Women/deltaAb_Mean_Women)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="deltaAb_Women",Mean=format(deltaAb_Mean_Women,digits = 3),
                                    CI=paste("[",format(deltaAb_CImin_Women,digits=3)," ; ",format(deltaAb_CImax_Women,digits=3),"]",sep="")))
  
  # delta Ab in Men
  beta_deltaAb_Mean_Men <- parameter_distributions$mean[which(parameter_distributions$parameter == "beta_delta_Ab_Sex_M")]
  deltaAb_Mean_Men <- deltaAb_Mean_Women*exp(beta_deltaAb_Mean_Men)
  
  beta_deltaAb_Sd_Men <- parameter_distributions$sd[which(parameter_distributions$parameter == "beta_delta_Ab_Sex_M")]
  deltaAb_CImin_Men <- exp(log(deltaAb_Mean_Women) + beta_deltaAb_Mean_Men - 1.96*(sqrt((beta_deltaAb_Sd_Men)^2 + (deltaAb_Sd_Women/deltaAb_Mean_Women)^2)))
  deltaAb_CImax_Men <- exp(log(deltaAb_Mean_Women) + beta_deltaAb_Mean_Men + 1.96*(sqrt((beta_deltaAb_Sd_Men)^2 + (deltaAb_Sd_Women/deltaAb_Mean_Women)^2)))
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="deltaAb_Men",Mean=format(deltaAb_Mean_Men,digits = 3),
                                    CI=paste("[",format(deltaAb_CImin_Men,digits=3)," ; ",format(deltaAb_CImax_Men,digits=3),"]",sep="")))
  
  
  # > b. Distribution of the half life of antibodies ####
  # Half life in Women
  HLAb_Mean_Women <- log(2)/deltaAb_Mean_Women
  
  HLAb_CImin_Women <- log(2)/deltaAb_CImax_Women
  HLAb_CImax_Women <- log(2)/deltaAb_CImin_Women
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="HL_Ab_Women",Mean=format(HLAb_Mean_Women,digits = 3),
                                    CI=paste("[",format(HLAb_CImin_Women,digits=3)," ; ",format(HLAb_CImax_Women,digits=3),"]",sep="")))
  
  # Half life in Men
  HLAb_Mean_Men <- log(2)/deltaAb_Mean_Men
  
  HLAb_CImin_Men <- log(2)/deltaAb_CImax_Men
  HLAb_CImax_Men <- log(2)/deltaAb_CImin_Men
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="HL_Ab_Men",Mean=format(HLAb_Mean_Men,digits = 3),
                                    CI=paste("[",format(HLAb_CImin_Men,digits=3)," ; ",format(HLAb_CImax_Men,digits=3),"]",sep="")))
  
  
  # > c. Distribution of the decay rate of SL ASCs ####
  deltaS_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "delta_S_pop")]
  
  deltaS_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "delta_S_pop")]
  deltaS_CImin <- exp(log(deltaS_Mean) - 1.96*deltaS_Sd/deltaS_Mean)
  deltaS_CImax <- exp(log(deltaS_Mean) + 1.96*deltaS_Sd/deltaS_Mean)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="deltaS",Mean=format(deltaS_Mean,digits = 3),
                                    CI=paste("[",format(deltaS_CImin,digits=3)," ; ",format(deltaS_CImax,digits=3),"]",sep="")))
  
  # > d. Distribution of the half life of SL ASCs ####
  HLSL_Mean <- log(2)/deltaS_Mean
  
  HLSL_CImin <- log(2)/deltaS_CImax
  HLSL_CImax <- log(2)/deltaS_CImin
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="HL_SL",Mean=format(HLSL_Mean,digits = 3),
                                    CI=paste("[",format(HLSL_CImin,digits=3)," ; ",format(HLSL_CImax,digits=3),"]",sep="")))
  
  
  # > e. Distribution of the decay rate of the LL ASCs ####
  deltaL_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "delta_L_pop")]
  
  deltaL_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "delta_L_pop")]
  deltaL_CImin <- parameter_distributions$CImin_total[which(parameter_distributions$parameter == "delta_L_pop")]
  deltaL_CImax <- parameter_distributions$CImax_total[which(parameter_distributions$parameter == "delta_L_pop")]
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="deltaL",Mean=format(deltaL_Mean,digits = 3),
                                    CI=paste("[",format(deltaL_CImin,digits=3)," ; ",format(deltaL_CImax,digits=3),"]",sep="")))
  
  # > f. Distribution of the half life of SL ASCs (in years) ####
  HLLL_Mean <- log(2)/(deltaL_Mean*365)
  
  HLLL_CImin <- log(2)/deltaL_CImax
  HLLL_CImax <- log(2)/deltaL_CImin
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="HL_LL",Mean=format(HLLL_Mean,digits = 3),
                                    CI=paste("[",format(HLLL_CImin,digits=3)," ; ",format(HLLL_CImax,digits=3),"]",sep="")))
  
  # > g. Distribution of the influx of SL ASCs ####
  # For subjects with the mean age (mean ~ 31.3 years)
  PhiS_Mean_meanAge <- parameter_distributions$mean[which(parameter_distributions$parameter == "phi_S_pop")]
  
  PhiS_Sd_meanAge <- parameter_distributions$sd[which(parameter_distributions$parameter == "phi_S_pop")]
  PhiS_CImin_meanAge <- exp(log(PhiS_Mean_meanAge) - 1.96*PhiS_Sd_meanAge/PhiS_Mean_meanAge)
  PhiS_CImax_meanAge <- exp(log(PhiS_Mean_meanAge) + 1.96*PhiS_Sd_meanAge/PhiS_Mean_meanAge)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="PhiS_meanAge",Mean=format(PhiS_Mean_meanAge,digits = 3),
                                    CI=paste("[",format(PhiS_CImin_meanAge,digits=3)," ; ",format(PhiS_CImax_meanAge,digits=3),"]",sep="")))
  
  # Fold change (+ 1 Year)
  beta_PhiS_Mean_meanAge <-  parameter_distributions$mean[which(parameter_distributions$parameter == "beta_phi_S_CenteredAge")]
  FC_PhiS_mean_1year <- exp(beta_PhiS_Mean_meanAge)
  
  beta_PhiS_Sd_meanAge <-  parameter_distributions$sd[which(parameter_distributions$parameter == "beta_phi_S_CenteredAge")]
  FC_PhiS_CImin_1year <- exp(beta_PhiS_Mean_meanAge - 1.96*beta_PhiS_Sd_meanAge)
  FC_PhiS_CImax_1year <- exp(beta_PhiS_Mean_meanAge + 1.96*beta_PhiS_Sd_meanAge)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="FC_PhiS_1year",Mean=format(FC_PhiS_mean_1year,digits = 3),
                                    CI=paste("[",format(FC_PhiS_CImin_1year,digits=3)," ; ",format(FC_PhiS_CImax_1year,digits=3),"]",sep="")))
  
  # > h. Distribution of the influx of the LL ASCs ####
  # PhiL in Africa (ref)
  PhiL_Mean_Africa <-  parameter_distributions$mean[which(parameter_distributions$parameter == "phi_L_pop")]
  
  PhiL_Sd_Africa <-  parameter_distributions$sd[which(parameter_distributions$parameter == "phi_L_pop")] 
  PhiL_CImin_Africa <- exp(log(PhiL_Mean_Africa) - 1.96*PhiL_Sd_Africa/PhiL_Mean_Africa)
  PhiL_CImax_Africa <- exp(log(PhiL_Mean_Africa) + 1.96*PhiL_Sd_Africa/PhiL_Mean_Africa)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="PhiL_Africa",Mean=format(PhiL_Mean_Africa,digits = 3),
                                    CI=paste("[",format(PhiL_CImin_Africa,digits=3)," ; ",format(PhiL_CImax_Africa),"]",sep="")))
  
  # PhiL in Europe
  beta_PhiL_Mean_Europe <- parameter_distributions$mean[which(parameter_distributions$parameter == "beta_phi_L_Continent_Europe")]
  PhiL_Mean_Europe <- PhiL_Mean_Africa*exp(beta_PhiL_Mean_Europe)
  
  beta_PhiL_Sd_Europe <- parameter_distributions$sd[which(parameter_distributions$parameter == "beta_phi_L_Continent_Europe")]
  PhiL_CImin_Europe <- exp(log(PhiL_Mean_Africa) + beta_PhiL_Mean_Europe - 1.96*sqrt(beta_PhiL_Sd_Europe^2 + (PhiL_Sd_Africa/PhiL_Mean_Africa)^2))
  PhiL_CImax_Europe <- exp(log(PhiL_Mean_Africa) + beta_PhiL_Mean_Europe + 1.96*sqrt(beta_PhiL_Sd_Europe^2 + (PhiL_Sd_Africa/PhiL_Mean_Africa)^2))
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="PhiL_Europe",Mean=format(PhiL_Mean_Europe,digits = 3),
                                    CI=paste("[",format(PhiL_CImin_Europe,digits=3)," ; ",format(PhiL_CImax_Europe),"]",sep="")))
  
  # > i. Scaling factors of lab effects ####
  # Laboratory FOCUS
  alphaFocus_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "alpha_focus_pop")]
  
  alphaFocus_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "alpha_focus_pop")]
  alphaFocus_CImin <- exp(log(alphaFocus_Mean) - 1.96*alphaFocus_Sd/alphaFocus_Mean)
  alphaFocus_CImax <- exp(log(alphaFocus_Mean) + 1.96*alphaFocus_Sd/alphaFocus_Mean)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="alpha_Focus",Mean=format(alphaFocus_Mean,digits = 3),
                                    CI=paste("[",format(alphaFocus_CImin,digits=3)," ; ",format(alphaFocus_CImax),"]",sep="")))
  
  # Laboratory Q2solution
  alphaQ2_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "alpha_Q2sol_pop")]
  
  alphaQ2_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "alpha_Q2sol_pop")]
  alphaQ2_CImin <- exp(log(alphaQ2_Mean) - 1.96*alphaQ2_Sd/alphaQ2_Mean)
  alphaQ2_CImax <- exp(log(alphaQ2_Mean) + 1.96*alphaQ2_Sd/alphaQ2_Mean)
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="alpha_Q2sol",Mean=format(alphaQ2_Mean,digits = 3),
                                    CI=paste("[",format(alphaQ2_CImin,digits=3)," ; ",format(alphaQ2_CImax),"]",sep="")))
  
  # > j. Distribution of the random effects ####
  # Random effects on delta Ab
  omega_deltaAb_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "omega_delta_Ab")]
  
  omega_deltaAb_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "omega_delta_Ab")]
  omega_deltaAb_CImin <- omega_deltaAb_Mean -1.96*omega_deltaAb_Sd
  omega_deltaAb_CImax <- omega_deltaAb_Mean +1.96*omega_deltaAb_Sd
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="omega_deltaAb",Mean=format(omega_deltaAb_Mean,digits = 3),
                                    CI=paste("[",format(omega_deltaAb_CImin,digits=3)," ; ",format(omega_deltaAb_CImax),"]",sep="")))
  
  # Random effects on PhiS
  omega_PhiS_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "omega_phi_S")]
  
  omega_PhiS_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "omega_phi_S")]
  omega_PhiS_CImin <- omega_PhiS_Mean -1.96*omega_PhiS_Sd
  omega_PhiS_CImax <- omega_PhiS_Mean +1.96*omega_PhiS_Sd
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="omega_PhiS",Mean=format(omega_PhiS_Mean,digits = 3),
                                    CI=paste("[",format(omega_PhiS_CImin,digits=3)," ; ",format(omega_PhiS_CImax),"]",sep="")))
  
  # Random effects on PhiL
  omega_PhiL_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "omega_phi_L")]
  
  omega_PhiL_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "omega_phi_L")]
  omega_PhiL_CImin <- omega_PhiL_Mean -1.96*omega_PhiL_Sd
  omega_PhiL_CImax <- omega_PhiL_Mean +1.96*omega_PhiL_Sd
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="omega_PhiL",Mean=format(omega_PhiL_Mean,digits = 3),
                                    CI=paste("[",format(omega_PhiL_CImin,digits=3)," ; ",format(omega_PhiL_CImax),"]",sep="")))
  
  # > k. Distribution of the error model ####
  sigma_Ab_Mean <- parameter_distributions$mean[which(parameter_distributions$parameter == "a")]
  
  sigma_Ab_Sd <- parameter_distributions$sd[which(parameter_distributions$parameter == "a")]
  sigma_Ab_CImin <- sigma_Ab_Mean - 1.96*sigma_Ab_Sd
  sigma_Ab_CImax <- sigma_Ab_Mean + 1.96*sigma_Ab_Sd
  
  table_results <- rbind(table_results,
                         data.frame(Parameter="sigma_Ab",Mean=format(sigma_Ab_Mean,digits = 3),
                                    CI=paste("[",format(sigma_Ab_CImin,digits=3)," ; ",format(sigma_Ab_CImax),"]",sep="")))
  
  
  return(table_results)
}
# ---------------- #



# We assume that an initial model has been build
Project_Folder <-  "Project_Folder_To_fill"
Project_Name <- "Project_Name_To_fill"      # Name of monolix project on which we want applied the model averaging approach

Results_Folder <-  "PROFILE LIKELIHOOD"
dir.create(paste(Project_Folder,Results_Folder,sep="/"),recursive = TRUE)


Tested_LL_HalfLife <- seq(1,40) # in years
Tested_delta_L <- log(2)/(Tested_LL_HalfLife*365.25) # in 1/day

# Names of monolix projects generated in "Results_Folder"
Final_project_names <- paste("PL_LLHalfLife_",Tested_LL_HalfLife,"years",sep="")


# ----- #
# -- PART 1: Estimation of candidate models -- ####
# ----- #
# Here, we include as potential candidate models all models involved in the profile likelihood of the fixed parameter delta L
# We applied the same code than in the R/ProfileLikelihood.R" file


# -- Run models' estimation ####
# \!/\!/ Can take long time according to the number of tested values and number of cores \!/\!/ 
# Initialization of lixoftconnectors 
monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify according to the Monolix version used
initializeLixoftConnectors(software = "monolix", path=monolix_path)

Init_project <- paste(Project_Folder,paste(Project_Name,"mlxtran",sep="."),sep="/")

Profile_likelihood_function(Init_project = Init_project, result_folder = paste(Project_Folder,Results_Folder,sep="/"),
                            parameter = "delta_L", tested_values = Tested_delta_L,final_project_names = Final_project_names,
                            Nb_Cores = 19,monolix_path = monolix_path)





# ----- #
# -- PART 2: Model averaging approach -- ####
# ----- #


# -- 1. Download of model estimation results ####
# We download LL-based criteria & population parameter estimates for each estimated model

Model_estimations <- NULL
LL_criteria <- NULL
for(i in 1:length(Tested_delta_L)){
  # i <- 1
  
  project_folder <- strsplit(Final_project_names[i],split = ".mlxtran",fixed=TRUE)[[1]]
  LLcriteria_result <- read.table(file=paste(Project_Folder,Results_Folder,project_folder,"LogLikelihood","logLikelihood.txt",sep="/"),sep=",",dec=".",stringsAsFactors = FALSE,header = TRUE)
  
  # \!/\!/ Code written for 2019R2 version. May require adaptation for newer monolix version \!/\!/ 
    LL_criteria <- rbind(LL_criteria,data.frame(deltaL = Tested_delta_L[i],LL_HalfLife=Tested_LL_HalfLife[i],
                                              LL=-0.5*as.numeric(LLcriteria_result$importanceSampling[which(LLcriteria_result$criteria == "-2LL")]),
                                              AIC=as.numeric(LLcriteria_result$importanceSampling[which(LLcriteria_result$criteria == "AIC")]),
                                              BIC=as.numeric(LLcriteria_result$importanceSampling[which(LLcriteria_result$criteria == "BIC")]),
                                              BICc=as.numeric(LLcriteria_result$importanceSampling[which(LLcriteria_result$criteria == "BICc")])))
  
  PopParams_results <- read.table(file=paste(Project_Folder,Results_Folder,project_folder,"populationParameters.txt",sep="/"),sep=",",dec=".",stringsAsFactors = FALSE,header = TRUE)
  PopParams_results <- cbind(deltaL = Tested_delta_L[i],LL_HalfLife=Tested_LL_HalfLife[i],PopParams_results)
  
  Model_estimations <- rbind(Model_estimations,PopParams_results)
}


# -- 2. Extraction of models to average and MA weight calculation ####
# We include in the set of candidate models for model averaging only those with a difference of AIC with the minimum AIC lower than 7  (Burhnam, 2004)
min_AIC <- min(LL_criteria$AIC,na.rm = T)
LL_criteria$Delta_AIC <- LL_criteria$AIC - min_AIC
Selected_models <- subset(LL_criteria,Delta_AIC<=7)

# Calculation of model averaging weights on these selected models
Selected_models$weights <- Model_averaging_weights_FUNCTION(values=Selected_models$AIC)

# Plot of the distribution of weights 
Weight_distribution_plot <- ggplot(data=Selected_models) + 
  geom_col(aes(x=as.factor(LL_HalfLife),y=weights),fill="gray30",alpha=1.0,col="black") + 
  coord_cartesian(ylim=c(0,0.25),expand = 0) + 
  xlab("LL ASCs Half-life (years)") + ylab("Model averaging weights") + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.5),
        axis.text.y =  element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        panel.grid.major.y = element_line(color="gray80",size=0.25),
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(color="white",size=0.8),
        legend.title = element_text(color="black",size=7,face="bold"),
        legend.text = element_text(color="gray20",size=7),
        legend.key.width = unit(1.0,"cm"), legend.key.height = unit(0.4,"cm"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        strip.background = element_rect(fill="navyblue",color="white"),
        strip.text = element_text(color="white",size=10,face="bold",margin=margin(t=0,r=0,b=0,l=0)),
        legend.position = c(0.6,0.6), 
        legend.direction = "vertical",
        legend.justification = "center",
        legend.box = "horizontal")



# -- 3. Calculation of averaged population parameters ####
Selected_Model_estimation <- subset(Model_estimations,LL_HalfLife %in% unique(Selected_models$LL_HalfLife))
Population_parameters <- unique(Model_estimations$parameter)
Parameters_distribution <- data.frame(parameter=Population_parameters,transformation=NA)
log_params <- c("delta_S_pop","phi_S_pop","phi_L_pop","delta_Ab_pop","alpha_focus_pop","alpha_Q2sol_pop")

Parameters_distribution$transformation <- ifelse(Parameters_distribution$parameter %in% log_params,"LN",
                                                 ifelse(Parameters_distribution$parameter =="delta_L_pop","Fixed","N"))


# Calculation of distribution of model parameters (we assume here a diagonal FIM - Otherwise, a bootstrap approach is required to sample population parameter from the FIM)
Averaged_PopParams <- NULL
for(p in 1:length(Population_parameters)){
  
  parameter <- Population_parameters[p]
  parameter_distribution <- Parameters_distribution$transformation[which(Parameters_distribution$parameter == parameter)]
  
  selected_model_estimation_param <- Selected_Model_estimation[Selected_Model_estimation$parameter==parameter,]
  parameter_weights <- setNames(Selected_models$weights,paste("HL =",selected_model_estimation_param$LL_HalfLife))
  parameter_values <- setNames(selected_model_estimation_param$value,paste("HL =",selected_model_estimation_param$LL_HalfLife))
  parameter_sd <- setNames(selected_model_estimation_param$se_sa,paste("HL =",selected_model_estimation_param$LL_HalfLife))
  
  # Calculation of mean and variances (inter and intra model)
  Averaged_mean <- as.numeric(parameter_values %*% parameter_weights)
  
  Averaged_intramodel_variance <- as.numeric(parameter_sd^2 %*% parameter_weights)
  Averaged_intermodel_variance <- as.numeric((parameter_values - Averaged_mean)^2 %*% parameter_weights)
  Averaged_total_variance <- Averaged_intramodel_variance + Averaged_intermodel_variance
  
  # calculation of confidence intervals according to the distribution of parameters
  if(parameter_distribution == "LN"){
    Averaged_CImin <- exp(log(Averaged_mean) - 1.96*sqrt(Averaged_total_variance)/Averaged_mean)
    Averaged_CImax <- exp(log(Averaged_mean) + 1.96*sqrt(Averaged_total_variance)/Averaged_mean)
  }else{
    Averaged_CImin <- Averaged_mean - 1.96*sqrt(Averaged_total_variance)
    Averaged_CImax <- Averaged_mean + 1.96*sqrt(Averaged_total_variance)
  }
  
  Averaged_PopParams <- rbind(Averaged_PopParams,
                              data.frame(parameter=parameter,mean=Averaged_mean,sd=sqrt(Averaged_total_variance),
                                         CImin=Averaged_CImin,CImax=Averaged_CImax))
}




# Results for averaged estimation (model averaging)
Distribution_Table_MA <- Parameter_distribution_TABLE(parameter_distributions = Averaged_PopParams)


# Results for all values of LL ASCs HL - distinct models
All_parameter_distribution <- Selected_Model_estimation
colnames(All_parameter_distribution) <- c("deltaL","LL_HalfLife","parameter","mean","sd","rse")

Distribution_Table_Allvalues <- setNames(lapply(unique(All_parameter_distribution$LL_HalfLife), function(halflife) Parameter_distribution_TABLE(parameter_distributions = subset(All_parameter_distribution,LL_HalfLife==halflife))),paste("HL =",unique(All_parameter_distribution$LL_HalfLife)))


