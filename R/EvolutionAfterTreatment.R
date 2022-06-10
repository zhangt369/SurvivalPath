
#'@title Display node transition with specified treatment plan or exposure
#'@description Calculate the number of subjects (proportion) assigned to different sub-nodes after specified treatment plan or exposure in certain node.
#'@usage EvolutionAfterTreatment(
#'df,
#'treepoint,
#'mytree,
#'source,
#'treatment
#')
#'@param df "data" in the return result of the \code{survivalpath()} function
#'@param treepoint list object;Specify the node for drawing the KM curve, the node number is displayed in the survival path tree graph.
#'@param  mytree  "tree" in the return result of the \code{survivalpath()} function
#'@param source Data.frame of time slice data, which could be returned by \code{timedivision()}
#'@param treatment Factor variable in the source data.frame. This argument is to specify the intervention or exposure that of interest at a specific node.
#'@return A data.frame object, whose rows and columns represents the number of subjects in the sub-nodes (in the next time slice) and treatment plan, respectively.
#'@export
#'@import dplyr
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
#'mytree <- result$tree
#'
#'#Draw the survival Path model
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
#'treepoint=15
#'A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Resection")
#'mytable <- xtabs(~ `Resection`+treepoint, data=A)
#'prop.table(mytable,1)
#'

EvolutionAfterTreatment <- function(df,treepoint,mytree,source,treatment){


  node = mytree@phylo[["edge"]][,"node"]
  parent = mytree@phylo[["edge"]][,"parent"]

  index <- which(parent==treepoint)

  parentDF <- getnodePatients(df,treepoint,mytree,source,treatment)
  parentDF$parentnode <- treepoint

  data.df <- data.frame()
  for (d in node[index]) {

    newdf <- getnodePatients(df,d,mytree,source,treatment)
    newdf$treepoint <- d
    #print(dim(newdf))
    data.df <- rbind(data.df,newdf)
  }
  data.df <- subset(data.df,select = -2)

  A<- merge(parentDF,data.df,all = TRUE,by.y = "ID")
  A$treepoint[which(is.na(A$treepoint))] <- "Missing follow-up"

  return(A)
}

getnodePatients <- function(df,treepoint,mytree,source,treatment){

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

  sourcedata <- source[which(source$time_slice==time_slice-1),]
  #print(dim(sourcedata))
  sourcedata <- sourcedata[sourcedata$ID %in% result$ID,]

  result <-  subset(sourcedata,select=c("ID"))
  #names(result) <- c("ID")
  result <- cbind(result,sourcedata[,treatment])
  names(result) <- c("ID",treatment)

  return(result)
}



