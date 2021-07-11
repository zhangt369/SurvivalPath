
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

NPC <- read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/NPC.xlsx", sheet = "Sheet2",
                       col_types = c("numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric"))
# 容错时间
#dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=1.5)
#write.csv(dataset,"D:/Code/R/SurvivalPath0225/Survivaldata/generator/mid4_4.csv", row.names = FALSE)

time <- list()

status <- list()

tsdata <- list()

tsid <- list()

treatment <- list()

for (i in 1:8){

  data <- NPC[NPC['timenode']==i,]

  time <- c(time,list(data['time']))

  status <- c(status,list(data['status']))

  tsid <- c(tsid,list(data['ID']))

  c_data <- subset(data, select = c("gender",	"age",	"weight",	"BMI",	"Tstage",	 "Nstage",	 "Mstage",	"TNM",
                                    "hypertension",	"diabetes-mellitu",	"WBC",	"RBC",	"HB",	"TPROT",	"ALB",	"CRP"))

  tsdata <- c(tsdata,list(c_data))

  c_treatment <- subset(data, select = c("Induced-chemotherapy-regimen"))

  treatment <- c(treatment,list(c_treatment))
}


# 输出cutoff值
tsdata <- classifydata(time,status,tsdata,tsid,cutoff=365*1)

result <- survivalpath(time,status,tsdata[[1]],tsid,time_slices = 7,treatments = treatment,p.value=0.05)

mytree <- result$tree
ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = T, )+
  theme_tree2() +
  xlim(c("放疗前","放疗后一个月","放疗后两个月"))+
  geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),hjust=0.6, vjust=-1.85 ,size=3.0)+
  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=F)+
  #guides(color=guide_legend(title="节点名称/样本量/中位生存时间/一年生存率")) +
  labs(size= "Nitrogen",
       x = "治疗节点",
       y = "鼻咽癌人群治疗后进展状态",
       subtitle = "节点名称/样本量/中位生存时间/一年生存率",
       title = "放疗生存进展树") +
  theme(legend.title=element_blank(),legend.position = c(0.1,0.9))

treepoints = c(5,3)
plotKM(result$data, treepoints,mytree,risk.table=T)

treepoints = c(17,22)
compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Treatment2")

treepoint=16
A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Treatment2")
mytable <- xtabs(~ `Treatment2`+treepoint, data=A)
prop.table(mytable,1)

