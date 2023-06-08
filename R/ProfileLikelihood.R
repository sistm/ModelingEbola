# ------------------------- #
# OBJECTIVE: Perform profile likelihood on model parameters
# Author: Marie Alexandre
# Date: 2023/06/08

# Related paper: Alexandre et al. (2023) - "Evaluation and prediction of the long-term humoral immune response induced by the two-dose heterologous Ad26.ZEBOV,MVA-BN-Filo vaccine regimen against Ebola"

# R version: 4.2.1
# Monolix version: >= 2019R1
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
# ---------------- #




# Initialization of lixoftconnectors 
monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify according to the Monolix version used
initializeLixoftConnectors(software = "monolix", path=monolix_path)




# We assume that an initial model has been build
Project_Folder <-  "Project_Folder_To_fill"
Project_Name <- "Project_Name_To_fill"


Results_Folder <-  "PROFILE LIKELIHOOD"
dir.create(paste(Project_Folder,Results_Folder,sep="/"),recursive = TRUE)




# ----- #
# -- Profile likelihood on delta L -- ####
# ----- #
Tested_LL_HalfLife <- seq(1,40) # in years
Tested_delta_L <- log(2)/(Tested_LL_HalfLife*365.25) # in 1/day


# Names of monolix projects generated in "Results_Folder"
Final_project_names <- paste("PL_LLHalfLife_",Tested_LL_HalfLife,"years",sep="")




# -- 1. Run profile Likelihood estimation ####
# \!/\!/ Can take long time according to the number of tested values and number of cores \!/\!/ 
Init_project <- paste(Project_Folder,paste(Project_Name,"mlxtran",sep="."),sep="/")

Profile_likelihood_function(Init_project = Init_project, result_folder = paste(Project_Folder,Results_Folder,sep="/"),
                            parameter = "delta_L", tested_values = Tested_delta_L,final_project_names = Final_project_names,
                            Nb_Cores = 19,monolix_path = monolix_path)






# -- 2. Plot of profile Likelihood #### 
# Extraction of the results 
Profile_likelihood_results <- NULL
for(i in 1:length(Tested_delta_L)){
  
  project_folder <- strsplit(Final_project_names[i],split = ".mlxtran",fixed=TRUE)[[1]]
  tmp_result <- read.table(file=paste(Project_Folder,Results_Folder,project_folder,"LogLikelihood","logLikelihood.txt",sep="/"),sep=",",dec=".",stringsAsFactors = FALSE)
  
  # \!/\!/ Code written for 2019R2 version. May require adaptation for newer monolix version \!/\!/ 
  Profile_likelihood_results <- rbind(Profile_likelihood_results,data.frame(deltaL = Tested_delta_L[i],LL_HalfLife=Tested_LL_HalfLife[i],
                                                                            LL=-0.5*as.numeric(tmp_result$V2[which(tmp_result$V1 == "-2LL")]),
                                                                            AIC=as.numeric(tmp_result$V2[which(tmp_result$V1 == "AIC")]),
                                                                            BIC=as.numeric(tmp_result$V2[which(tmp_result$V1 == "BIC")]),
                                                                            BICc=as.numeric(tmp_result$V2[which(tmp_result$V1 == "BICc")])))
}



X_labels <- data.frame(x=c(seq(1,14,by=2),seq(15,40,by=5)),color="black",stringsAsFactors = FALSE)
X_labels$color[which(X_labels$x == 15)] <- "red"
X_labels$color[which(X_labels$x == 5)] <- "blue"


Profile_Likelihood_deltaL_PAPER <- ggplot(data=Profile_likelihood_results,aes(x=LL_HalfLife,y=LL)) +
  
  geom_vline(aes(xintercept=15,color="new",linetype="new"),size=0.75,key_glyph = "path") + 
  geom_vline(aes(xintercept=5,color="old",linetype="old"),size=0.75,key_glyph = "path") + 
  geom_line(size=0.5,color="black") + 
  
  scale_color_manual(name="Lower bound of LL ASCs Half-life",breaks=c("old","new"),values=c("blue","red"),
                     labels=c("Previous estimate (Phase I data)","New estimate (Phase I-II data)")) + 
  scale_linetype_manual(name="Lower bound of LL ASCs Half-life",breaks=c("old","new"),values=c("dotted","dashed"),
                        labels=c("Previous estimate (Phase I data)","New estimate (Phase I-II data)")) + 
  
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.5),
        axis.text.y =  element_text(color="black",size=8),
        axis.text.x = element_text(color=X_labels$color,size=8),
        axis.title = element_text(color="black",size=9,face="bold"),
        panel.grid.minor = element_blank(),
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
        legend.box = "horizontal") +
  
  scale_x_continuous(name= "LL ASCs Half-life (years)",breaks=X_labels$x,limits=c(0.5,43),expand = c(0,0)) +
  scale_y_continuous(name="Log-Likelihood",breaks=seq(-0,500,by=25),expand = c(0,0)) +
  coord_cartesian(ylim=c(0,NA))
