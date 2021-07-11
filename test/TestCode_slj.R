library(readxl)
library(dplyr)
library(survival)
library(rms)
library(survivalROC)
library(survival)
library(survminer)#??????Ҫ
library(treeio) #??????Ҫ
library(ggtree) #???в???Ҫ

data("dataset")
# ?ݴ?ʱ??
dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=1.5)  #这个函数需要再检查一下，明确数据整理的规则。
#write.csv(dataset,"D:/Code/R/SurvivalPath0225/Survivaldata/generator/mid4_4.csv", row.names = FALSE)
table(dataset$timenode)  #查看时间片的数据

time <- list()

status <- list()

tsdata <- list()

tsid <- list()

treatment <- list()

for (i in 1:10){

  data <- dataset[dataset['time_slice']==i,]  #改成time slice

  time <- c(time,list(data['OStime_new']))

  status <- c(status,list(data['Status_new']))

  tsid <- c(tsid,list(data['ID']))

  c_data <- subset(data, select = c('Age','Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)','New Lesion','Vascular Invasion','Local Lymph Node Metastasis',
                                    'Distal Metastasis','Ascites','Massive Ascites','Moderate or Mild Ascites'
                                    ,'Infected with Hepatitis','ALB',"TBLT",'PT',"AFP"))

  tsdata <- c(tsdata,list(c_data))

  c_treatment <- subset(data, select = c("Treatment"))  #最好不要带2，直接Treatment

  treatment <- c(treatment,list(c_treatment))
}


# 生存ROC对连续性及分类变量找cutoffֵ
tsdata <- classifydata(time,status,tsdata,tsid,cutoff=365*1)

#varname=list('Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)') ,varvalue=list(1,1)
#df <- selectgroup(time,status,tsdata[[1]],tsid,varname=
#                       list('Amount of Hepatic Lesions') ,varvalue=list(1))   #selectgroup名称改matchsubgroup()

#ggtree
#result <- survivalpath(df$time,df$status,df$tsdata,df$tsid,time_slices=7)  #这行代码用于52行


result <- survivalpath(time,status,tsdata[[1]],tsid,time_slices = 10,treatments = treatment,p.value=0.05,degreeofcorrelation=0.7)  #time_slices为允许跑的时间节点的最大值
#(Q)1.主函数输出结果，varname要改成varvalue, Class是否需要舍去，输出结果再加一个输出tree节点（数值）信息的表格，
#(Q)2.主函数中，参数Survival rate要选择率的时间单位
#(S)3.主函数中，共线性的相关变量的排除的相关系数节值如何设置？最好有个Arguement.

mytree <- result$tree    #(Q)1.图像没有画到最后一个时间片，需要再看一下作图的数据
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

#生存曲线
treepoints = c(16,23)
plotKM(result$data, treepoints,mytree,risk.table=T)

#治疗方式生存曲线
treepoints = c(17,22)
compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Treatment")  #治疗方式最好只选一个

treepoint=16
A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Treatment")
mytable <- xtabs(~ `Treatment2`+treepoint, data=A)
prop.table(mytable,1)

