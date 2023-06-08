# ------------------------- #
# OBJECTIVE: Perform Monte-carlo cross-validation (MCCV)
# Author: Marie Alexandre
# Date: 2023/06/08

# Related paper: Alexandre et al. (2023) - "Evaluation and prediction of the long-term humoral immune response induced by the two-dose heterologous Ad26.ZEBOV,MVA-BN-Filo vaccine regimen against Ebola"

# R version: 4.2.1
# Monolix version: 2019R1 or 2019R2  (use of package mlxR which is not suitable for monolix version >= 2020)
# ------------------------- #




rm(list=ls())
`%notin%` <- Negate(`%in%`)



# --- LIBRARIES --- ####
# - Packages for Monolix
library(lixoftConnectors)
library(mlxR)  # used in PART 2 to simulate individual parameters (percent of coverage)

# - Package for parallel calculation
library(parallel) ; library(snow) ; library(doSNOW)

library(plyr) ; library(dplyr)

# - Package for plots (PART 3)
library(ggplot2) ;  library(lemon)
library(gridExtra) ; library(cowplot) ; library(grid)
# ---------------- #



# --- FUNCTIONS --- ####
# - PART 2
Root_Mean_Squared_Error <- function(observed_data,predicted_data,censure=NULL){
  if(!is.null(censure) & sum(censure) == 0){
    censure <- NULL
  }
  
  if(is.null(censure)){ # No censored data
    res <- sqrt(sum( (predicted_data-observed_data)^2 )/length(observed_data))
  }else{
    # we split data according to the censure
    obs_data_nocens <- observed_data[which(censure == 0)]
    pred_data_nocens <- predicted_data[which(censure == 0)]
    square_diff_nocens <- (pred_data_nocens-obs_data_nocens)^2
    
    obs_data_cens <- observed_data[which(censure == 1)]
    pred_data_cens <- predicted_data[which(censure == 1)]
    square_diff_cens <- NULL
    for(i in 1:length(obs_data_cens)){
      obs_i <- obs_data_cens[i]
      pred_i <- pred_data_cens[i]
      if(pred_i <= obs_i){
        square_diff_cens <- c(square_diff_cens,0)
      }else{
        square_diff_cens <- c(square_diff_cens,(pred_i-obs_i)^2)
      }
    }
    
    square_diff <- c(square_diff_nocens,square_diff_cens)
    res <- sqrt(sum(square_diff)/length(observed_data))
  }
  return(res)
}
# ---------------- #




# Initialization of lixoftconnectors 
monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify according to the Monolix version used
initializeLixoftConnectors(software = "monolix", path=monolix_path)


# We assume that the model on which we apply MCCV approach has already been estimated (only loaded)
Project_Folder <-  "Project_Folder_To_fill"
Project_Name <- "Project_Name_To_fill"      # Name of monolix project on which we want applied the MCCV approach  (without mlxtran extension)

Results_Folder <-  "CROSS VALIDATION"  # General folder in which we gather results obtained for different train-test split percentages (1 folder for each percentage)
dir.create(paste(Project_Folder,Results_Folder,sep="/"),recursive = TRUE)






# Calculation taking time, we didn't apply for loop but tested each percentage separately on calculation server
MC_CV_folder <- "MCCV_60"   # Folder gathering results for a given train-test split percentage (e.g 60% here)
dir.create(paste(Results_Folder,MC_CV_folder,sep="/"),recursive = TRUE)

Project_Name_init <- paste(Project_Folder,"/",Project_Name,".mlxtran",sep="")
New_Project_Name <- paste(Results_Folder,"/",MC_CV_folder,"/","Model_Init.mlxtran",sep="")

## Load the model 
loadProject(Project_Name_init)
saveProject(New_Project_Name)
runScenario()
saveProject()

Data <- read.csv2(file=getData()$dataFile,header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
Trials <- unique(Data$trial)
Subject_distribution <- unique(Data[,c("ID","Country","trial","Continent","Region")])

Nb_loop <- 100  # Number of iterations
Training_prop <- 0.60 # proportion of the dataset (keeping ratio between groups) selected for training 





# ----- #
# -- PART 1: MCCV approach - parallel calculation -- ####
# ----- #
Nb_Cores <- 19
NCores <- min(detectCores()-1,Nb_Cores)
Cluster <- parallel::makeCluster(NCores,type="SOCK")
registerDoSNOW(Cluster)
pb <- txtProgressBar(max=Nb_loop,style=3)         # We add a progress bar
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)
clusterExport(Cluster,c("Nb_loop","Results_Folder","MC_CV_folder","Training_prop","Data","Trials","Subject_distribution",
                        "New_Project_Name","monolix_path"))



Results <- foreach(n=seq(1,Nb_loop),.errorhandling = "remove",.packages = c("doParallel","foreach","doSNOW","lixoftConnectors"),.options.snow = opts)%dopar%{
  initializeLixoftConnectors(software = "monolix", path = monolix_path)
  
  dir.create(paste(Results_Folder,"/",MC_CV_folder,"/","Simu_",n,sep=""))
  tmp_folder <- paste(Results_Folder,"/",MC_CV_folder,"/","Simu_",n,sep="")
  
  # > 1.a Training dataset creation ####
  # For each group, we have to keep the same proportion of subject
  Training_subjects <- setNames(lapply(Trials,function(trial,train_prop,subjects){
    subjects_trial <- subjects$ID[which(subjects$trial == trial)]
    Nb_subj <- round(train_prop*length(subjects_trial)) # number of subject to select
    selected_subj <- sample(x=subjects_trial,size=Nb_subj,replace=FALSE)   # random sampling
    return(selected_subj)
  },train_prop=Training_prop,subjects=Subject_distribution),Trials)
  Training_dataset <- do.call("rbind",lapply(Trials,function(trial,data,training_subj){
    return(data[which(data$ID %in% training_subj[[trial]]),])
  },data=Data,training_subj=Training_subjects))
  
  # > 1.b Testing dataset creation ####
  Testing_subjects <- setNames(lapply(Trials,function(trial,training_subj,subjects){
    subjects_trial <- subjects$ID[which(subjects$trial == trial)]
    selected_subj <- setdiff(subjects_trial,training_subj[[trial]])
    return(selected_subj)
  },training_subj=Training_subjects,subjects=Subject_distribution),Trials)
  Testing_dataset <- do.call("rbind",lapply(Trials,function(trial,data,testing_subj){
    return(data[which(data$ID %in% testing_subj[[trial]]),])
  },data=Data,testing_subj=Testing_subjects))
  
  # Record of data
  write.csv2(Training_dataset,file=paste(tmp_folder,"Training_dataset.csv",sep="/"),row.names = FALSE)
  write.csv2(Testing_dataset,file=paste(tmp_folder,"Testing_dataset.csv",sep="/"),row.names = FALSE)
  
  
  # load of the initial model
  loadProject(New_Project_Name)
  MC_CV_Training_modelName <- paste(tmp_folder,"Model_Estimation.mlxtran",sep="/")
  saveProject(MC_CV_Training_modelName)
  
  # > 1.c Estimation of the model on training dataset ####
  header <- getData()$header ; headerType <- getData()$headerTypes
  observationName <- getData()$observationNames ; observationType <- getData()$observationTypes
  setData(dataFile=paste(tmp_folder,"Training_dataset.csv",sep="/"),
          headerTypes=headerType,observationTypes=observationType)
  saveProject()
  runScenario()
  saveProject()
  
  # > 1.d Estimation of EBEs on testing dataset ####
  MC_CV_Testing_modelName <- paste(tmp_folder,"Model_Test.mlxtran",sep="/")
  saveProject(MC_CV_Testing_modelName)
  Population_Estimates <- getEstimatedPopulationParameters()
  
  # Population parameters are fixed at values estimated on training data
  Initial_Pop_Estimates <- getPopulationParameterInformation()
  Initial_Pop_Estimates$initialValue <- sapply(seq(1,nrow(Initial_Pop_Estimates)),function(i){
    # i <- 1
    if(Initial_Pop_Estimates$method[i] == "FIXED"){
      return(as.numeric(Initial_Pop_Estimates$initialValue[i]))
    }else{
      return(as.numeric(Population_Estimates[which(names(Population_Estimates) == Initial_Pop_Estimates$name[i])]))
    }
  })
  Initial_Pop_Estimates$method <- "FIXED"
  setPopulationParameterInformation(Initial_Pop_Estimates)
  
  # Modification of the dataset
  header <- getData()$header ; headerType <- getData()$headerTypes
  observationName <- getData()$observationNames ; observationType <- getData()$observationTypes
  setData(dataFile=paste(tmp_folder,"Testing_dataset.csv",sep="/"),
          headerTypes=headerType,observationTypes=observationType)
  saveProject()
  runScenario()
  saveProject()
}
close(pb)
stopCluster(Cluster)





# ----- #
# -- PART 2: RMSE and percent of coverage estimation - parallel calculation -- ####
# ----- #
initializeLixoftConnectors(software = "simulx", path=monolix_path,force = TRUE) # we need to switch from monolix to simulx
initMlxR(path=monolix_path)

Nb_loop <- 100 # Number of iterations
Nb_Cores <- 19 # Number of cores to use for parallel calculation

file_individual_simulation <- "Mlxtran/MlxranModel_Individual_Simulation.txt"  # File used by simulx function to simulate individual parameters and the resulting trajectories to estimated percent of coverage
 

# We merge results for different values of train-test split percentages 
MC_CV_folders <- c(paste("MCCV_",seq(20,80,by=10),sep=""))
results_files <- c(paste("MCCV_",seq(20,80,by=10),"training.csv",sep=""))


# - Evaluation of the MCCV performances
for(m in 1:length(MC_CV_folders)){
  # m <- 1
  MC_CV_folder <- MC_CV_folders[m] 
  print(MC_CV_folder)
  
  # Parameters for parallel calculation
  NCores <- min(detectCores()-1,Nb_Cores)
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  pb <- txtProgressBar(max=Nb_loop,style=3)         # We add a progress bar
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress=progress)
  clusterExport(Cluster,c("Nb_loop","Results_Folder","Root_Mean_Squared_Error","MC_CV_folder","file_individual_simulation","monolix_path"))
  
  
  MC_CV_criteria <- NULL
  MC_CV_criteria <- foreach(n=seq(1,Nb_loop),.errorhandling = "remove",.packages = c("doParallel","foreach","doSNOW","mlxR","plyr","dplyr"),.options.snow = opts)%dopar%{

    initMlxR(path=monolix_path)

    tmp_folder <- paste(Results_Folder,"/",MC_CV_folder,"/","Simu_",n,sep="")
    
    Training_dataset <- read.csv2(file=paste(tmp_folder,"Training_dataset.csv",sep="/"),header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
    Testing_dataset <- read.csv2(file=paste(tmp_folder,"Testing_dataset.csv",sep="/"),header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
    
    # > 1.e Estimation of the Root Mean Square Error (RMSE) ####
    # We estimate the RMSE on the testing dataset
    Predicted_data <- read.csv(file=paste(tmp_folder,"Model_Test","predictions.txt",sep="/"),header=TRUE,stringsAsFactors = FALSE)
    # Addition of information about censure
    Predicted_data$Cens <- Testing_dataset$Cens
    
    if(sum(Predicted_data$Cens == 0)){
      RMSE_SAEM <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_SAEM)
      RMSE_mean <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_mean)
      RMSE_mode <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_mode)
    }else{
      RMSE_SAEM <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_SAEM,censure=Predicted_data$Cens)
      RMSE_mean <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_mean,censure=Predicted_data$Cens)
      RMSE_mode <- Root_Mean_Squared_Error(observed_data=Predicted_data$logAb,predicted_data=Predicted_data$indivPred_mode,censure=Predicted_data$Cens)
    }
    
    # > 1.f Estimation of the mean Individual Weighted Residual (mean IWRES) ####
    mean_IWRES_SAEM <- mean(Predicted_data$indWRes_SAEM,na.rm=TRUE)
    mean_IWRES_mean <- mean(Predicted_data$indWRes_mean,na.rm=TRUE)
    mean_IWRES_mode <- mean(Predicted_data$indWRes_mode,na.rm=TRUE)
    
    
    # > 1.g Estimation of the percentage of coverage ####
    parameter_names <- c("delta_L","delta_S","phi_L","phi_S","delta_Ab","alpha_focus","alpha_Q2sol")
    # Download of individual parameter distribution
    Individual_parameters <- read.csv(file=paste(tmp_folder,"Model_Test","IndividualParameters","estimatedIndividualParameters.txt",sep="/"),header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
    # Download of population parameter
    Population_parameters <- read.csv(file=paste(tmp_folder,"Model_Test","populationParameters.txt",sep="/"),header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
    Testing_subjects <- unique(Individual_parameters$id)
    Time_simulation <- seq(0,720,by=10)
    
    Distribution_simu <- NULL
    for(id in Testing_subjects){
      parameters_id <- Individual_parameters[which(Individual_parameters$id == id),]
      
      # Selection of parameters (mode)
      params_value <- parameters_id[,which(names(parameters_id) %in% paste(parameter_names,"_mode",sep=""))]
      names(params_value) <- paste(unlist(strsplit(names(params_value),split="_mode",fixed=TRUE)),"pop",sep="_")
      params_sd <- parameters_id[,which(names(parameters_id) %in% paste(parameter_names,"_sd",sep=""))]
      names(params_sd) <- paste("omega",unlist(strsplit(names(params_sd),split="_sd",fixed=TRUE)),sep="_")
      params <- as.numeric(c(params_value,params_sd,Population_parameters$value[which(Population_parameters$parameter == "a")]))
      names(params) <- c(names(params_value),names(params_sd),"a")
      
      params["omega_phi_L"] <- 0 #Population_parameters$value[which(Population_parameters$parameter == "omega_phi_L")]
      params["omega_phi_S"] <- 0 #Population_parameters$value[which(Population_parameters$parameter == "omega_phi_S")]
      params["omega_delta_Ab"] <- 0 #Population_parameters$value[which(Population_parameters$parameter == "omega_delta_Ab")]
      
      names(params)[which(names(params) == "alpha_focus_pop")] <- "alpha_focus"
      names(params)[which(names(params) == "alpha_Q2sol_pop")] <- "alpha_Q2sol"
      
      
      # definition of the output of the model
      individual_time <- sort(unique(c(Predicted_data$time[which(Predicted_data$id == id)],Time_simulation)))
      f <- list(name="log_AbValue",time=individual_time)
      y <- list(name="y",time=individual_time)
      g <- list(size=100,level='individual')
      # definition of the regressors
      reg_AbTi <- list(name="Ab_Ti",time=individual_time,value=rep(unique(Predicted_data$Ab_Ti[which(Predicted_data$id == id)]),length(individual_time)))
      reg_Ti <- list(name="Ti",time=individual_time,value=rep(unique(Predicted_data$Ti[which(Predicted_data$id == id)]),length(individual_time)))
      reg_Test <- list(name="Test",time=individual_time,value=rep(unique(Predicted_data$Test[which(Predicted_data$id == id)]),length(individual_time)))
      
      simulations <- simulx(model = file_individual_simulation,
                            parameter=params,regressor=list(reg_AbTi,reg_Ti,reg_Test),group=g,
                            output=list(f,y,list(name="delta_L"),list(name="delta_S"),list(name="phi_L"),list(name="phi_S"),list(name="delta_Ab")))
      
      simulated_parameters <- simulations$parameter
      simulated_dynamics <- simulations$log_AbValue
      simulated_noised_dynamics <- simulations$y
      
      distribution_noised_dynamics <- ddply(.data=simulated_noised_dynamics,.variables = .(time),summarize,
                                            Mean=mean(y,na.rm=TRUE),Median=median(y,na.rm=TRUE),
                                            ICMIN=quantile(y,na.rm=TRUE,probs = c(0.025)),ICMAX=quantile(y,na.rm=TRUE,probs = c(0.975)))
      
      # Estimation of the percentage of data within the 95% CI
      distribution_observedtime <- distribution_noised_dynamics[which(distribution_noised_dynamics$time %in% Predicted_data$time[which(Predicted_data$id == id)]),]
      distribution_observedtime$Obs <- Predicted_data$logAb[which(Predicted_data$id == id)]
      distribution_observedtime$in.CI <- sapply(seq(1,nrow(distribution_observedtime)),function(i){
        1*between(distribution_observedtime$Obs[i],left=distribution_observedtime$ICMIN[i],right=distribution_observedtime$ICMAX[i])
      })
      
      Distribution_simu <- rbind(Distribution_simu,cbind(id=id,distribution_observedtime))
    }
    Coverage_perc <- sum(Distribution_simu$in.CI)*100/nrow(Distribution_simu)
    
    tmp_mccv_criteria <- data.frame(iter=n,RMSE_mean=RMSE_mean,RMSE_mode=RMSE_mode,RMSE_SAEM=RMSE_SAEM,
                                    mean_IWRES_SAEM=mean_IWRES_SAEM,mean_IWRES_mean=mean_IWRES_mean,
                                    mean_IWRES_mode=mean_IWRES_mode,Coverage=Coverage_perc)
    write.csv2(tmp_mccv_criteria,file=paste(tmp_folder,"mccv_criteria.csv",sep="/"),row.names = FALSE)
    return(tmp_mccv_criteria)
  }
  close(pb)
  stopCluster(Cluster)
  
  MC_CV_criteria <- do.call("rbind",MC_CV_criteria)
  write.csv2(MC_CV_criteria,file=paste(Results_Folder,MC_CV_folder,results_files[m],sep="/"),row.names = FALSE)
}





# ----- #
# -- PART 3: Plot of results -- ####
# ----- #
MC_CV_folders <- c(paste("MCCV_",seq(20,80,by=10),sep=""))
results_files <- c(paste("MCCV_",seq(20,80,by=10),"training.csv",sep=""))


# - Download of results
MCCV_results <- NULL
for(m in 1:length(MC_CV_folders)){
  # m <- 1
  MC_CV_folder <- MC_CV_folders[m]
  Prop <- as.numeric(substr(MC_CV_folder,start=nchar(MC_CV_folder)-1,stop=nchar(MC_CV_folder)))
  tmp_results <- read.csv2(file=paste(Results_Folder,MC_CV_folder,results_files[m],sep="/"),header=TRUE,stringsAsFactors = FALSE)  
  MCCV_results <- rbind(MCCV_results,cbind(Prop=Prop,tmp_results))
}



# - Estimation of the distribution of the coverage and of the RMSE
MCCV_RMSE_distribution <- ddply(.data=MCCV_results,.variables = .(Prop),summarise,Mean=mean(RMSE_mode,na.rm=TRUE),
                                ICMIN=quantile(RMSE_mode,na.rm=TRUE,probs=c(0.025)),ICMAX=quantile(RMSE_mode,na.rm=TRUE,probs=c(0.975)))
MCCV_Coverage_distribution <- ddply(.data=MCCV_results,.variables = .(Prop),summarise,Mean=mean(Coverage,na.rm=TRUE),
                                    ICMIN=quantile(Coverage,na.rm=TRUE,probs=c(0.025)),ICMAX=quantile(Coverage,na.rm=TRUE,probs=c(0.975)))
MCCV_IWRES_distribution <- ddply(.data=MCCV_results,.variables = .(Prop),summarise,Mean=mean(mean_IWRES_mode,na.rm=TRUE),
                                 ICMIN=quantile(mean_IWRES_mode,na.rm=TRUE,probs=c(0.025)),ICMAX=quantile(mean_IWRES_mode,na.rm=TRUE,probs=c(0.975)))

MCCV_criteria_distribution <- rbind(cbind(Criteria="RMSE",MCCV_RMSE_distribution),
                                    cbind(Criteria="Coverage",MCCV_Coverage_distribution),
                                    cbind(Criteria="IWRES",MCCV_IWRES_distribution))


# - Plot
Criteria_labels <- setNames(c("Root Mean Square Error (RMSE)","Coverage (%)"),c("RMSE","Coverage"))

# RMSE
MCCV_RMSE_Plot <- ggplot(data=subset(MCCV_criteria_distribution,Criteria=="RMSE")) + 
  geom_line(aes(x=Prop,y=Mean),size=1.2,color="black") + 
  geom_ribbon(aes(x=Prop,ymin=ICMIN,ymax=ICMAX),color="black",linetype="dashed",size=0.75,fill=NA) + 
  facet_rep_wrap(.~Criteria,repeat.tick.labels = TRUE,scales = "free",labeller=labeller(Criteria=Criteria_labels)) + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=1.5),
        axis.text.y =  element_text(color="black",size=10),
        axis.text.x =  element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10,face="bold"),
        legend.background = element_rect(color="white",size=0.8),
        legend.title = element_text(color="black",size=11,face="bold"),
        legend.text = element_text(color="black",size=10),
        legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.2,"cm"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        strip.background = element_rect(fill="white",color="white"),
        strip.text = element_text(color="black",size=10,face="bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box = "horizontal") +
  xlab("Percentage of subjects used for training datasets") + 
  ylab("Value calculated on testing datasets") + 
  scale_x_continuous(expand=c(0,0),breaks = seq(0,100,by=10),limits=c(20,82)) + 
  scale_y_continuous(limits=c(0.05,0.125),expand = c(0,0.001),breaks = seq(0,0.2,0.015))


# Percent coverage
y_breaks <- data.frame(y=c(seq(90,100,by=2),95),color=c(rep("black",length(seq(90,100,by=2))),"red"))

MCCV_Coverage_plot <- ggplot(data=subset(MCCV_criteria_distribution,Criteria=="Coverage")) + 
  geom_hline(yintercept = 95,linetype="dotted",color="red",size=1) + 
  geom_line(aes(x=Prop,y=Mean),color="black",size=1.2) + 
  geom_ribbon(aes(x=Prop,ymin=ICMIN,ymax=ICMAX),color="black",linetype="dashed",size=0.75,fill=NA) + 
  facet_rep_wrap(.~Criteria,repeat.tick.labels = TRUE,scales = "free",labeller=labeller(Criteria=Criteria_labels)) + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=1.5),
        axis.text.y =  element_text(color=y_breaks$color,size=10),
        axis.text.x =  element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10,face="bold"),
        legend.background = element_rect(color="white",size=0.8),
        legend.title = element_text(color="black",size=11,face="bold"),
        legend.text = element_text(color="black",size=10),
        legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.2,"cm"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        strip.background = element_rect(fill="white",color="white"),
        strip.text = element_text(color="black",size=10,face="bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box = "horizontal") +
  xlab("Percentage of subjects used for training datasets") + 
  ylab("Value calculated on testing datasets") + 
  scale_x_continuous(expand=c(0,0),breaks = seq(0,100,by=10),limits=c(20,82)) + 
  scale_y_continuous(limits=c(90,100),expand = c(0,0),breaks = y_breaks$y) 



# Merge of the 2 plots
grid <- grid.arrange(MCCV_RMSE_Plot + theme(axis.title=element_blank()),
                     MCCV_Coverage_Plot + theme(axis.title = element_blank()),
                     ncol=2,left=textGrob("Value on testing datasets",rot = 90, gp = gpar(fontsize = 11, fontface = 'bold')),
                     bottom=textGrob("Percent of subjects used within training datasets",vjust = 0, gp = gpar(fontsize = 11, fontface = 'bold')))