
#'@title  Calculate the number of people (proportion) assigned to different branches after a certain
#'node is treated by a specified plan
#'@description Calculate the number of people (proportion) assigned to different branches after a
#'certain node is treated by a specified plan
#'@usage EvolutionAfterTreatment(
#'df,
#'treepoint,
#'mytree,
#'source,
#'treatment
#')
#'@param df "data" in the return result of the survivalpath function
#'@param treepoint  Numerical object;Specify the node for drawing the KM curve, which is in the survival tree
#'@param  mytree  "mytree" in the return result of the survivalpath function
#'@param source Raw data
#'@param treatment Choose treatment
#'@details By specifying the node, using mytree to calculate the relevant branch variables and subsequent nodes,
#'and using the original data to calculate, filter out the group of subjects who branch to the next node after
#'the node adopts different treatment plans. Calculate the number of people assigned to different nodes according
#'to different treatment methods (proportion).
#'@return A dataframe whose rows and columns are the next node and treatment plan
#'@export
#'@examples data("dataset")
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
#'treepoint=16
#'A = EvolutionAfterTreatment(result$data, treepoint,mytree,dataset,"Treatment")
#'mytable <- xtabs(~ `Treatment2`+treepoint, data=A)
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
    print(dim(newdf))
    data.df <- rbind(data.df,newdf)
  }
  data.df <- subset(data.df,select = -2)

  A<- merge(parentDF,data.df,all = T,by.y = "ID")
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

  sourcedata <- source[which(source$timenode==time_slice-1),]
  #print(dim(sourcedata))
  sourcedata <- sourcedata[sourcedata$ID %in% result$ID,]

  result <-  subset(sourcedata,select=c("ID"))
  names(result) <- c("ID")
  result <- cbind(result,sourcedata[,treatment])

  return(result)
}



