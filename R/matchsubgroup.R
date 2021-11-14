
#'@title Screen and collect data of subjects that meet the given conditions
#'@description The data of the subjects that meet the conditions for the first time and the data of
#'subsequent time slices will also be selected. Their data from initial time slice that meet the given conditions to last time slice were compiled into a new time-slice dataset
#'@usage matchsubgroup(
#'time,
#'status,
#'timeslicedata,
#'tspatientid,
#'varname,
#'varvalue
#')
#'@param time list object;Elements sorted by time node, each element is a
#'Dataframe object, representing event time or censoring time for subjects
#'@param status list object;Elements sorted by time node, each element is a
#'Dataframe object, representing status, 1 if death or event, 0 otherwise.
#'@param  timeslicedata list object; Elements sorted by time node, each element is a
#'Dataframe object, representing risk factors for the subject
#'@param tspatientid list object; Elements sorted by time node, each element represents the subjects identification number
#'@param  varname list object;The variable used to screen subjects, and the variable
#'needs to be a variable in timeslicedata
#'@param varvalue  list object;Subjects whose varname variable is equal to varvalue will be selected
#'@details According to the input time, status, variables, subject ID, etc., the data of eligible
#'subjects is screened through specified conditions, and the subject and the variable data of the
#'first and subsequent time slices are re-screened. The final returned result contains four list
#'objects: time, state, timeslicedata, subject ID(tspatientid).
#'@return Returns a list of the following items:
#'time,state,timeslicedata,tspatientid
#'@author Shen Lujun and ZhangTao
#'@export
#'@examples
#'data("DTSDHCC")
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
#'varname=list('Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)')
#'varvalue=list(1,1)
#'df <- matchsubgroup(time,status,tsdata[[1]],tsid,varname=
#'                      list('Amount of Hepatic Lesions') ,varvalue=list(1))

#'#ggtree
#'result <- survivalpath(df$time,df$status,df$timeslicedata,df$tspatientid,time_slices=9)
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

matchsubgroup <- function(time,status,timeslicedata,tspatientid,varname,varvalue){

  d_time <- data.frame()
  d_status <- data.frame()
  d_tsdata <- data.frame()
  d_tsid <- data.frame()
  d_timenode <- data.frame()
  for (i in 1:length(time)) {
    time_slice <- rep(i,dim(time[[i]])[1])
    time_slice <- data.frame(time_slice)
    d_timenode <- rbind(d_timenode,time_slice)
    d_time <- rbind(d_time,time[[i]])
    d_status <- rbind(d_status,status[[i]])
    d_tsdata <- rbind(d_tsdata,timeslicedata[[i]])
    d_tsid <- rbind(d_tsid,tspatientid[[i]])
  }

  dataset <- cbind(d_time,d_status,d_tsdata,d_tsid,d_timenode)
  names(d_tsid) <- "id_id"
  dataset <- cbind(dataset,d_tsid)

  names(d_time) <- "time_time"
  dataset <- cbind(dataset, d_time)

  dataset <- dataset %>% group_by(id_id) %>% arrange(id_id,time_slice)

  data <- dataset

  if (!is.list(varname)) {
    stop('varname should be a list!')
  }
  if (!is.list(varvalue)) {
    stop('varvalue should be a list!')
  }

  if (length(varname)!=length(varvalue)) {
    stop('The length ofvarvalue should equal with the length of varvalue!')
  }

  len <- length(varname)

  if (len==0){
    stop("Variables not specified for selecting data")
  }


  for (i in 1:len) {

    dataset <- dataset[which(dataset[,varname[[i]]]==varvalue[[i]]),]
  }

  ind <- duplicated(dataset$id_id)
  id_index <- subset(dataset,select = c("id_id","time_slice"))
  id_index = data.frame(id_index)
  id_index = id_index[!ind,]

  n_data = data.frame()
  for (i in 1:dim(id_index)[1]) {
    m_data = data[which(data$id_id==id_index$id_id[i]),]
    m_data = m_data[which(m_data$time_slice<=id_index$time_slice[i]),]
    n_data = rbind(n_data,m_data)
  }

  timename = names(time[[1]])
  statusname = names(status[[1]])
  varname = names(timeslicedata[[1]])

  num <- length(time)

  time <- list()

  status <- list()

  timeslicedata <- list()

  tspatientid <- list()

  for (i in 1:num){

    data <- n_data[n_data['time_slice']==i,]
    if(dim(data)[1]==0){
      break()
    }

    time <- c(time,list(data[timename]))

    status <- c(status,list(data[statusname]))

    ID <- data['id_id']
    names(ID) = "ID"

    tspatientid <- c(tspatientid,list(ID))

    c_data <- subset(data, select =varname)

    timeslicedata <- c(timeslicedata,list(c_data))
  }

  result <- list()
  result$time = time
  result$status = status
  result$tspatientid = tspatientid
  result$timeslicedata = timeslicedata
  return(result)
}
