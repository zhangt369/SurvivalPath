# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#'@title  Compare and Draw the KM curve of any given node
#'@description According to the survival tree, draw the KM curve of the nodes on the survival tree
#'@usage plotKM(
#'df,
#'treepoints,
#'mytree,
#'risk.table=T
#')
#'@param df "data" in the return result of the survivalpath function
#'@param treepoints list object;Specify the node for drawing the KM curve, which is in the survival tree
#'@param  mytree "mytree" in the return result of the survivalpath function
#'@param risk.table Allowed values include:
#'TRUE or FALSE specifying whether to show or not the risk table. Default is FALSE."absolute" or "percentage".
#'Shows the absolute number and the percentage of subjects at risk by time, respectively. "abs_pct" to show both absolute number and percentage.
#'"nrisk_cumcensor" and "nrisk_cumevents". Show the number at risk and, the cumulative number of censoring and events, respectively.
#'@details draw the KM curve
#'@seealso survminer
#'@export
#'@examples data("DTSDHCC")
#'dataset = timedivision(DTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'
#'time <- list()
#'status <- list()
#'tsdata <- list()
#'tsid <- list()
#'
#'treatment <- list()
#'for (i in 1:10){
#'
#'  data <- dataset[dataset['time_slice']==i,]
#'
#'  time <- c(time,list(data['OStime_day']))
#'
#'  status <- c(status,list(data['Status_of_death']))
#'
#'  tsid <- c(tsid,list(data['ID']))
#'
#'  c_data <- subset(data, select = c( "Age", "Amount of Hepatic Lesions", "Largest Diameter of Hepatic Lesions (mm)", "New Lesion",
#'    "Vascular Invasion" ,"Local Lymph Node Metastasis", "Distant Metastasis" , "Child_pugh_score" ,"AFP"))
#'
#'
#'  tsdata <- c(tsdata,list(c_data))
#'
#'  c_treatment <- subset(data, select = c("Resection"))
#'
#'  treatment <- c(treatment,list(c_treatment))
#'}
#'
#'tsdata <- classifydata(time,status,tsdata,tsid,predict.time=365*1)
#'
#'result <- survivalpath(time,status,tsdata[[1]],tsid,time_slices = 10,treatments = treatment,p.value=0.05,
#'degreeofcorrelation=0.7)
#'
#'mytree <- result$tree
#'
#'library(ggtree)
#'library(ggplot2)
#'ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = TRUE )+
#'  theme_tree2() +
#'  geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
#'  geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),
#'  hjust=0.6, vjust=-1.85 ,size=3.0)+
#'  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
#'  geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size%/%200+1,show.legend=FALSE)+
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
#'


plotKM <- function(df,treepoints,mytree,risk.table=T){
  data.df <- data.frame()
  for (d in treepoints) {
    newdf <- findPatients(df,d,mytree)
    newdf$treepoint <- d
    data.df <- rbind(data.df,newdf)
  }
  #data.df$treepoint <- factor(data.df$treepoint)

  fit <- survfit(Surv( time,status) ~ treepoint, data = data.df)
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

findPatients <- function(df,treepoint,mytree){

  node = mytree@phylo[["edge"]][,"node"]
  parent = mytree@phylo[["edge"]][,"parent"]

  tip.point = sort(setdiff(node,parent))
  tip.label = mytree@phylo[["tip.label"]]

  node.point = sort(unique(parent))
  node.label= mytree@phylo[["node.label"]]

  if(treepoint> max(node) | treepoint< min(node)){
    stop("Out of node range")
  }

  index = which(node==treepoint)
  #print(index)
  branch = vector()

  branch = c(branch,node[index],parent[index])
  while ( parent[index] %in% node ) {
    index = which(node==parent[index])
    branch = c(branch,parent[index])
  }

  for (i in length(branch):1){
    if (i==length(branch)){
      newdf <- df
    }else{
      if (branch[i] %in% node.point){
        variable <- node.label[which(node.point==branch[i])]
        varname <- substr(variable, 1,nchar(variable)-2)
        value <- substr(variable, nchar(variable),nchar(variable))

        newdf <- newdf[which(newdf[,paste("time_slices_",length(branch)-i,"_variable",sep = "")]==varname,),]
        newdf <- newdf[which(newdf[,paste("time_slices_",length(branch)-i,"_varname",sep = "")]==value,),]
      }

    }
  }
  result <- newdf[,1:3]
  #result$time_slice <- rep(treepoint,dim(result)[1])
  return(result)
}






