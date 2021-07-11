# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#'@title  Draw the KM curve of different nodes after the specified treatment plan
#'@description Based on the survival tree, specify nodes and specify treatment methods,
#'draw survival curves to compare prognostic effects.
#'@usage compareTreatmentPlans(
#'df,
#'treepoints,
#'mytree,
#'source,
#'treatment
#')
#'@param df "data" in the return result of the survivalpath function
#'@param treepoints list object;Specify the node for drawing the KM curve, which is in the survival tree
#'@param  mytree  "mytree" in the return result of the survivalpath function
#'@param source Raw data
#'@param treatment Choose treatment
#'@details By specifying the node, using mytree to calculate the previous branch variables,
#'and calculating with the original data, the group of subjects with different treatment plans
#'at different nodes can be screened out. Draw the survival curve of patients based on the selected groups
#'@seealso survminer
#'
#'@export
#'@examples
#'data("dataset")
#'dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=1.5)
#'
#'time <- list()
#'status <- list()
#'tsdata <- list()
#'tsid <- list()
#'
#'treatment <- list()
#'for (i in 1:10){
#'
#'  data <- dataset[dataset['timenode']==i,]
#'
#'  time <- c(time,list(data['OStime_new']))
#'
#'  status <- c(status,list(data['Status_new']))
#'
#'  tsid <- c(tsid,list(data['ID']))
#'
#'  c_data <- subset(data, select = c('Age','Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)','New Lesion','Vascular Invasion','Local Lymph Node Metastasis',
#'                                    'Distal Metastasis','Ascites','Massive Ascites','Moderate or Mild Ascites'
#'                                    ,'Infected with Hepatitis','ALB',"TBLT",'PT',"AFP"))
#'
#'  tsdata <- c(tsdata,list(c_data))
#'
#'  c_treatment <- subset(data, select = c("Treatment2"))
#'
#'  treatment <- c(treatment,list(c_treatment))
#'}
#'
#'tsdata <- classifydata(time,status,tsdata,tsid,cutoff=365*1)
#'
#'result <- survivalpath(time,status,tsdata[[1]],tsid,time_slices = 10,treatments = treatment,p.value=0.05,degreeofcorrelation=0.7)
#'
#'mytree <- result$tree
#'ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = T, )+
#'  theme_tree2() +
#'  geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
#'  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),hjust=0.6, vjust=-1.85 ,size=3.0)+
#'  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
#'  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=F)+
#'  #guides(color=guide_legend(title="node name/sample number/Median survival time/Survival rate")) +
#'  labs(size= "Nitrogen",
#'       x = "TimePoints",
#'       y = "Survival",
#'       subtitle = "node_name/sample number/Median survival time/Survival rate",
#'       title = "Survival Tree") +
#'  theme(legend.title=element_blank(),legend.position = c(0.1,0.9))
#'
#'#plot KM curve
#'treepoints = c(16,23)
#'plotKM(result$data, treepoints,mytree,risk.table=T)
#'#Comparing the efficacy of treatment methods by drawing survival curves
#'treepoints = c(17,22)
#'compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Treatment")
#'

compareTreatmentPlans <- function(df,treepoints,mytree,source,treatment){
  data.df <- data.frame()
  for (d in treepoints) {
    newdf <- getPatients(df,d,mytree,source,treatment)
    newdf$treepoint <- d
    data.df <- rbind(data.df,newdf)
  }
  #print(table(data.df[,3:4]))
  name <- names(data.df)
  fit <- survfit(Surv( time,status)~Treatment2+treepoint, data = data.df)
  ggsurvplot(fit,data = data.df,
             pval = TRUE,
             #conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             #surv.median.line = "hv",
             ggtheme = theme_bw(),
             #palette = c("#E7B800", "#2E9FDF"))
  )

}

getPatients <- function(df,treepoint,mytree,source,treatment){

  node = mytree@phylo[["edge"]][,"node"]
  parent = mytree@phylo[["edge"]][,"parent"]

  tip.point = sort(setdiff(node,parent))
  tip.label = mytree@phylo[["tip.label"]]

  node.point = sort(unique(parent))
  node.label= mytree@phylo[["node.label"]]

  #get root node
  non.tip.point = sort(setdiff(node,tip.point))
  rootnode = setdiff(parent,non.tip.point)

  if(treepoint> max(node) | treepoint< min(node)){
    stop("Out of node range")
  }

  #get treepoint survival path
  path <- c(treepoint)
  current_node <- treepoint
  while (current_node!=rootnode) {
    index <- which(node==current_node)
    current_node <- parent[index]
    path <- c(path,current_node)
  }
  #get  treepoint survival periods
  time_slice <- length(path)
  print(path)

  for (i in time_slice:1){

    if (i==time_slice){
      newdf <- df
    }else{

      if (path[i] %in% node.point){
        variable <- node.label[which(node.point==path[i])]

        if("all" %in% variable){
          next
        }
        varname <- substr(variable, 1,nchar(variable)-2)

        value <- substr(variable, nchar(variable),nchar(variable))
        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_variable",sep = "")]==varname,),]
        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_varname",sep = "")]==value,),]

      }
      if (path[i] %in% tip.point){
        variable <- tip.label[which(tip.point==path[i])]

        if("all" %in% variable){
          next
        }
        varname <- substr(variable, 1,nchar(variable)-2)

        value <- substr(variable, nchar(variable),nchar(variable))

        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_variable",sep = "")]==varname,),]
        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_varname",sep = "")]==value,),]
      }

    }
  }

  result <- newdf[,1:3]
  print(dim(result))

  sourcedata <- source[which(source$timenode==time_slice-1),]
  print(dim(sourcedata))

  sourcedata <- sourcedata[sourcedata$ID %in% result$ID,]

  result <-  subset(sourcedata,select=c("OStime_new","Status_new"))
  names(result) <- c("time","status")
  result <- cbind(result,sourcedata[,treatment])

  return(result)
}






