# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#'@title  Compare and Draw the KM curve of specified treatment plan or exposure in selected nodes
#'@description Based on the survival tree, specify the node of interest and the treatment methods, draw survival curves to evaluate the impact of treatments or exposure.
#'@usage compareTreatmentPlans(
#'df,
#'treepoints,
#'mytree,
#'source,
#'treatment
#')
#'@param df "data" in the return result of the \code{survivalpath()} function
#'@param treepoints list object;Specify the node for drawing the KM curve, which is displayed in the survival path graphs
#'@param  mytree  "tree" in the return result of the \code{survivalpath()} function
#'@param source Data.frame of time slice data, which could be returned by \code{timedivision()}
#'@param treatment Factor variable in the source data.frame. This argument is to specify the intervention or exposure that of interest at a specific node.
#'@details The function creates survival curves of specified treatment plan or exposure in selected nodes. The results should be interpreted with caution
#'as the effect of covariates have not been adjusted.
#'@return No return value.
#'@seealso survminer
#'@import survival
#'@importFrom stats as.formula
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
#'#Comparing the efficacy of treatment methods by drawing survival curves
#'treepoints = c(14,20)
#'compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Resection")
#'

compareTreatmentPlans <- function(df,treepoints,mytree,source,treatment){
  data.df <- data.frame()
  for (d in treepoints) {
    newdf <- getPatients(df,d,mytree,source,treatment)
    newdf$treepoint <- d
    data.df <- rbind(data.df,newdf)
  }

  formula <- paste("Surv(time,status)~",paste(names(data.df)[3:4],collapse ="+") )

  formula <- as.formula(formula)

  fit <- surv_fit(formula,data=data.df)

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
  #print(path)

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
        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_varvalue",sep = "")]==value,),]

      }
      if (path[i] %in% tip.point){
        variable <- tip.label[which(tip.point==path[i])]

        if("all" %in% variable){
          next
        }
        varname <- substr(variable, 1,nchar(variable)-2)

        value <- substr(variable, nchar(variable),nchar(variable))

        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_variable",sep = "")]==varname,),]
        newdf <- newdf[which(newdf[,paste("time_slices_",time_slice-i,"_varvalue",sep = "")]==value,),]
      }

    }
  }

  result <- newdf[,1:3]
  #print(dim(result))

  sourcedata <- source[which(source$time_slice==time_slice-1),]
  #print(dim(sourcedata))

  sourcedata <- sourcedata[sourcedata$ID %in% result$ID,]

  result <-  subset(sourcedata,select=c("OStime_day","Status_of_death"))
  names(result) <- c("time","status")
  result <- cbind(result,sourcedata[,treatment])
  names(result) <- c("time","status",treatment)
  return(result)
}






