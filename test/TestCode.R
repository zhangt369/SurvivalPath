
library(SurvivalPath)

library(readxl)
library(dplyr)
library(survival)
library(rms)
library(survivalROC)
library(survival)
library(survminer)#包中需要
library(treeio) #包中需要
library(ggtree) #包中不需要
#X2021data <- read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/20210704data.xlsx")
X2021data <-  read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/20210704data.xlsx",
                           col_types = c("text", "numeric", "date",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric",
                                         "numeric", "date", "numeric"))
data("dataset")

dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=1.5)

time <- list()

status <- list()

tsdata <- list()

tsid <- list()

treatment <- list()

for (i in 1:10){

  data <- dataset[dataset['timenode']==i,]

  time <- c(time,list(data['OStime_new']))

  status <- c(status,list(data['Status_of_death'])) #Status_of_death

  tsid <- c(tsid,list(data['ID']))

  c_data <- subset(data, select = c('Age','Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)','New Lesion','Vascular Invasion','Local Lymph Node Metastasis',
                                    'Distal Metastasis','Ascites','Massive Ascites','Moderate or Mild Ascites'
                                    ,'Infected with Hepatitis','ALB',"TBLT",'PT',"AFP"))

  tsdata <- c(tsdata,list(c_data))

  c_treatment <- subset(data, select = c("Resection"))

  treatment <- c(treatment,list(c_treatment))
}


# cutoff
tsdata <- classifydata(time,status,tsdata,tsid,cutoff=365*1)

#varname=list('Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)') ,varvalue=list(1,1)
#df <- matchsubgroup(time,status,tsdata[[1]],tsid,varname=
#                       list('Amount of Hepatic Lesions') ,varvalue=list(1))

#ggtree
#result <- survivalpath(df$time,df$status,df$tsdata,df$tsid,time_slices=7)
#df <- survivalpath(time,status,tsdata[[1]],tsid,time_slices = 9)


result <- survivalpath(time,status,tsdata[[1]],tsid,time_slices=10)

mytree <- result$tree
ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = T, )+
  theme_tree2() +
  geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),hjust=0.6, vjust=-1.85 ,size=3.0)+
  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=F)+
  #guides(color=guide_legend(title="node name/sample number/Median survival time/Survival rate")) +
  labs(size= "Nitrogen",
       x = "TimePoints",
       y = "Survival",
       subtitle = "node_name/sample number/Median survival time/Survival rate",
       title = "Survival Tree") +
  theme(legend.title=element_blank(),legend.position = c(0.1,0.9))

treepoints = c(16,23)
plotKM(result$data, treepoints,mytree,risk.table=T)

treepoints = c(17,22)
compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Treatment2")

treepoint=16
A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Treatment2")
mytable <- xtabs(~ `Treatment2`+treepoint, data=A)
prop.table(mytable,1)

