
library(SurvivalPath)
library(readxl)
library(dplyr)
library(survival)
library(rms)
library(survivalROC)
library(survival)
library(survminer)
library(treeio)
library(ggtree)

data("dataset")

dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)

time <- list()

status <- list()

tsdata <- list()

tsid <- list()

treatment <- list()

for (i in 1:10){

  data <- dataset[dataset['time_slice']==i,]

  time <- c(time,list(data['OStime_new']))

  status <- c(status,list(data['Status_of_death']))

  tsid <- c(tsid,list(data['ID']))

  c_data <- subset(data, select = c('Age','Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)','New Lesion','Vascular Invasion','Local Lymph Node Metastasis',
                                    'Distal Metastasis','Ascites','Massive Ascites','Moderate or Mild Ascites'
                                    ,'Infected with Hepatitis','ALB',"TBLT",'PT',"AFP"))

  tsdata <- c(tsdata,list(c_data))

  c_treatment <- subset(data, select = c("Resection"))

  treatment <- c(treatment,list(c_treatment))
}


# cutoff
tsdata <- classifydata(time,status,tsdata,tsid,predict.time=365*5)

#conduct survival path mapping on whole dataset
result <- survivalpath(time, status, tsdata[[1]], tsid, time_slices = 8, treatments = treatment, p.value=0.05, degreeofcorrelation=0.7, rates=365)

#draw survival paths
mytree <- result$tree
ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = T, )+
  theme_tree2() +
  geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),hjust=0.6, vjust=-1.85 ,size=3.0)+
  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=F)+
  labs(size= "Nitrogen",
       x = "TimePoints",
       y = "Survival",
       subtitle = "node_name/sample number/Median survival time/Survival rate",
       title = "Survival Tree") +
  theme(legend.title=element_blank(),legend.position = c(0.1,0.9))

#draw KM plots
treepoints = c(16,23)
plotKM(result$data, treepoints,mytree,risk.table=T)

treepoints = c(40,42)
compareTreatmentPlans(result$data, treepoints,mytree,dataset,"TargetedTherapy")

#Crosstabs for treatment and node change
treepoint=40
A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Resection")
mytable <- xtabs(~ `Resection`+treepoint, data=A)
prop.table(mytable,1)

#Generate subgroup data based on personalized design
varname=list('Amount of Hepatic Lesions')
varvalue=list(1)
df <- matchsubgroup(time,status,tsdata[[1]],tsid,
                    varname=varname ,varvalue=varvalue)

#conduct survival path mapping based on personalized requirement
result <- survivalpath(df$time,df$status,df$timeslicedata,df$tspatientid,time_slices=9, degreeofcorrelation=0.7, Rates=365)
