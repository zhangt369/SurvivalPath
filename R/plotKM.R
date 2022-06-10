
#'@title  Compare and Draw the KM curves of any given nodes
#'@description According to the survival path tree, draw the KM curves of the using any nodes on the survival tree
#'@usage plotKM(
#'df,
#'treepoints,
#'mytree,
#'risk.table=TRUE
#')
#'@param df "data" in the returned result of the \code{survivalpath()} function
#'@param treepoints list object;Specify the node for drawing the KM curve, which is in the survival path tree
#'@param  mytree "tree" in the returned result of the \code{survivalpath()} function
#'@param risk.table Logical value. Allowed values include:TRUE or FALSE specifying whether to show the risk table. Default is FALSE.
#'@details Plot survival curves for patients contained in nodes in the survival path tree.
#'@return No return value.
#'@seealso survminer
#'@importFrom ggplot2 theme_bw
#'@export
#'@examples
#'library(dplyr)
#'data("DTSDHCC")
#'id = DTSDHCC$ID[!duplicated(DTSDHCC$ID)]
#'set.seed(123)
#'id = sample(id,500)
#'miniDTSDHCC <- DTSDHCC[DTSDHCC$ID %in% id,]
#'dataset = timedivision(miniDTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'resu <- generatorDTSD(dataset,periodindex="time_slice",IDindex="ID" ,timeindex="OStime_day",
#'  statusindex="Status_of_death",variable =c( "Age", "Amount.of.Hepatic.Lesions",
#'  "Largest.Diameter.of.Hepatic.Lesions",
#'  "New.Lesion","Vascular.Invasion" ,"Local.Lymph.Node.Metastasis",
#'  "Distant.Metastasis" , "Child_pugh_score" ,"AFP"),predict.time=365*1)
#'result <- survivalpath(resu,time_slices =9)
#'
#'mytree <- result$tree
#'
#'library(ggplot2)
#'library(ggtree)
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
#'treepoints = c(14,20)
#'plotKM(result$data, treepoints,mytree,risk.table=T)
#'


plotKM <- function(df,treepoints,mytree,risk.table=TRUE){
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
        newdf <- newdf[which(newdf[,paste("time_slices_",length(branch)-i,"_varvalue",sep = "")]==value,),]
      }

    }
  }
  result <- newdf[,1:3]
  #result$time_slice <- rep(treepoint,dim(result)[1])
  return(result)
}






