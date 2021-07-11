library(readxl)
library(dplyr)
library(survival)
library(rms)
X2021data <- read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/2021data.xlsx",
                               col_types = c("text", "numeric", "date",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric",
                                                       "date", "numeric"))
# 容错时间
dataset = timedivision(X2021data,"ID","Date",period = 30,low = 15,hight=45)
write.csv(dataset,"D:/Code/R/SurvivalPath0225/Survivaldata/generator/mid4_4.csv", row.names = FALSE)

time <- list()

status <- list()

tsdata <- list()

tsid <- list()

for (i in 1:9){

  data <- dataset[dataset['timenode']==i,]

  time <- c(time,list(data['OS时间']))

  status <- c(status,list(data['Status']))

  tsid <- c(tsid,list(data['ID']))

  c_data <- subset(data, select = c('Age','Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)','New Lesion','Vascular Invasion','Local Lymph Node Metastasis',
                                    'Distal Metastasis','Ascites','Massive Ascites','Moderate or Mild Ascites'
                                    ,'Infected with Hepatitis','ALB',"TBLT",'PT',"AFP"))

  tsdata <- c(tsdata,list(c_data))
}

# 输出cutoff值
tsdata <- classifydata(time,status,tsdata,tsid)
result <- survivalpath(time,status,tsdata[[1]],tsid,timepoints=9)

df = structureResult(result,time,status,tsid)

mytree = df2newick(df,innerlabel = T,survivaltime = 3,period=12)

ggtree(mytree, color="black",linetype=1,size=1.,ladderize = T, )+
  theme_tree2() +
  geom_text2(aes(label=label),hjust=0.6, vjust=-0.5 ,size=2.4)+
  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),hjust=0.6, vjust=-1.65 ,size=2.4)+
  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=F)+
  #guides(color=guide_legend(title="node/sample number/Median survival time/Survival rate")) +
  labs(size= "Nitrogen",
       x = "TimePoints",
       y = "Survival",
       subtitle = "node/sample number/Median survival time/Survival rate",
       title = "Survival Tree") +
  theme(legend.title=element_blank(),legend.position = c(0.1,0.9))
