
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#'@title  Build survival path, generate survival tree
#'@description build survival path, plot survival tree
#'@usage survivalpath(
#'time,
#'status,
#'timeslicedata,
#'tspatientid,
#'time_slices,
#'treatments=NULL,
#'num_categories=2,
#'p.value=0.05,
#'minsample = 15,
#'degreeofcorrelation=0.7,
#'rates
#')
#'@param time list object;Elements sorted by time node, each element is a
#'Dataframe object, representing event time or censoring time for subjects
#'@param status  list object;Elements sorted by time node, each element is a
#'Dataframe object, representing status, 1 if death or event, 0 otherwise.
#'@param  timeslicedata list object; Elements sorted by time node, each element is a
#'Dataframe object, representing risk factors for the subject
#'@param tspatientid list object; Elements sorted by time node, each element is a
#'Dataframe object, representing the subjects’ identification number
#'@param time_slices null
#'@param treatments default NULL.It is possible to specify the intervention measures taken
#'by the subjects at different time points. list object
#'@param num_categories default 2. The maximum number of branches that each node can divide
#'@param p.value p.value for univariate selection; variables less than p.value are given meaning
#'@param minsample Minimum sample size for branching
#'@param degreeofcorrelation default 0.7;When the correlation between variables is greater than this value,
#'the variables are considered to be correlated.
#'@param rates This parameter specifies the time point at which each node calculates the survival rate.
#'@details for survivalpath
#'@return Returns a list of the following items:
#' data dataframe:contains the main risk factors and corresponding values of each
#'subject divided at different time points
#' tree survival tree of newick structure
#'@author Shen Lujun and ZhangTao
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


# data.tree

survivalpath <- function(time,status,timeslicedata,tspatientid,time_slices,treatments=NULL,num_categories=2,
                         p.value=0.05,minsample = 15,degreeofcorrelation=0.7,rates=365) {
  l_time <- length(time)
  l_status <- length(status)
  l_timesl <- length(timeslicedata)
  l_id <- length(tspatientid)


  ######################################################################################################
  #format

  if (!is.list(time)) {
    stop('time should be a list of data.frame!')
  }
  if (!is.list(status)) {
    stop('status should be a list of data.frame!')
  }
  if (!is.list(timeslicedata)) {
    stop('timeslicedata should be a list of data.frame!')
  }
  if (!is.list(tspatientid)) {
    stop('tspatientid should be a list of data.frame!')
  }

  #Data preprocessing, class 2 variables are transformed into class 2 variables

  #timeslicedata <- classifydata(time,status,timeslicedata,tspatientid)

  #Judging the number of input data time nodes is consistent
  if (l_time!=l_status | l_time!= l_timesl | l_time != l_id){
    stop(sprintf('The length of time,status,timeslicedata and tspatientid are:%d, %d, %d, %d,
                                but they should be the same length.',l_time,l_status,l_timesl,l_id))
  }

  if(!is.null(treatments)){
    l_tr <- length(treatments)
    if (!is.list(treatments)) {
      stop('treatments should be a list of data.frame!')
    }

    if (l_time!=l_status | l_time!= l_timesl | l_time != l_id |l_time != l_tr){
      stop(sprintf('The length of time,status,timeslicedata , tspatientid,treatments are:%d, %d, %d, %d,%d,
                                but they should be the same length.',l_time,l_status,l_timesl,l_id,l_tr))
    }
  }

  #A patient's unique coding id required
  if (l_id==0){
    stop('The input can\'t be empty')
  }


  #A patient's unique coding id required,judgement
  patientID <- tspatientid[[1]]


  patientID_a <-unique(patientID[duplicated(patientID),])

  lenpid <- length(as.numeric(unlist(patientID_a)))
  if (lenpid!=0){
    stop(paste('The patient ID repeat. Every patient ID should be unique:',patientID_a))
  }
  ###############################################################################################################
  #id in patientID


  #Number of time nodes, equal to the length of input data
  #time_slices <- 9

  #all patients ID
  branchID <- list(tspatientid[[1]])

  #init
  ts_var <- vector("list", time_slices)
  ts_num <- vector("list", time_slices)
  ts_id <- vector("list", time_slices)
  ts_value <- vector("list", time_slices)

  alldata <- list()
  for(i in 1:l_time){
    ts_data <- data.frame( tspatientid[[i]], time[[i]],status[[i]],class=1)
    names(ts_data) <- c("ID","time","status","class")
    if(!is.null(treatments[[i]])){
      ts_data <- data.frame(ts_data,treatments[[i]],timeslicedata[[i]])
    }
    else{
      ts_data <- data.frame(ts_data,timeslicedata[[i]])
    }

    alldata$data <- rbind(alldata$data,ts_data)
  }
  alldata$treatmentsnum <- dim(treatments[[i]])[2]


  for (time_slice in 1:time_slices){
    #print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print(paste('time_slice:',time_slice,sep = ''))

    previousbranchID <- branchID  #list


    #Full sample time at current time node
    current_time <- time[[time_slice]]

    #Full sample status of current time node
    current_status <- status[[time_slice]]

    #Full sample id of current time node
    current_id <- tspatientid[[time_slice]]

    #All sample variables for the current time node: data frames
    current_variables <- timeslicedata[[time_slice]]

    #Print variable name
    #print(paste("datavar:",names(current_variables)))


    #Determine whether the dimensions of the four data are consistent
    if (dim(current_time)[1]!=dim(current_status)[1] | dim(current_time)[1]!= dim(current_variables)[1] |
        dim(current_time)[1] != dim(current_id)[1]){stop(sprintf('The number of the %d-th timeslice of time,status,
                        timeslicedata and tspatientid are:%d, %d, %d, %d, but they should be the same ',time_slice,
                                                                 l_time,l_status,l_timesl,l_id))}


    # id,time,status,varibales
    current_ts_data <- data.frame(current_id,current_time,current_status,current_variables)

    #Dimensional Judgment
    if (dim(current_time)[2]!=1){
      stop(sprintf('the %d-th timeslice of survival time should be unique, but the dim is:%d',2,dim(current_time)[2]))
    }

    if (dim(current_status)[2]!=1){
      stop(sprintf('the %d-th timeslice of survival status should be unique, but the dim is:%d',2,dim(current_status)[2]))
    }

    if (dim(current_id)[2]!=1){
      stop(sprintf('the %d-th timeslice of patient id should be unique, but the dim is:%d',2,dim(current_id)[2]))
    }


    ############################################# Former branch  ##########################################################

    # Number of branches retained at the previous moment previousbranchid-> the id,list() of each branch,data.frame,brnum:number of branches
    brnum <- length(previousbranchID)

    if(brnum<1){
      ts_var = ts_var[-c(time_slice:time_slices)]

      ts_num = ts_num[-c(time_slice:time_slices)]

      ts_id = ts_id[-c(time_slice:time_slices)]

      ts_value = ts_value[-c(time_slice:time_slices)]
      next
    }


    #Current time, variable name of current branch
    ts_br_var <- vector('list',brnum)
    #ts_br_v_names <- vector('list',brnum)


    # The number of branches of the current branch
    ts_br_num <- vector('list',brnum)


    # ??????
    branchID <- list()
    branchID_var_value <- list()

    for (i_id in 1:brnum) {

      #print('########################################')

      #print(paste('current branch:',i_id,sep = ''))


      #ts_br_v_names[[i_id]] <- paste('ts',time_slice,'_br',i_id)

      #id of the previous time node to the current branch
      #print(i_id)


      br_id <- data.frame(previousbranchID[[i_id]])

      nam <- names(current_ts_data)

      names(br_id) <- nam[1]
      #print(names(br_id)[1])

      #Dimensions of current survival data

      #print("pre branch:")
      #print(dim(current_ts_data))


      #To extract current_ts_data data rows in the br_id by name nam[1]"
      current_ts_b_data <- merge(current_ts_data,br_id,by= nam[1],all=FALSE)

      #print("after merge, current branch:")
      #print(dim(current_ts_b_data))


      #id of current time node branches
      current_ts_b_id <- subset(current_ts_b_data, select = nam[1])

      current_names <- isbranch(current_ts_b_id,minsample = minsample) #Empty or stop

      #Judge whether to proceed
      if (!is.null(current_names)) {
        ts_br_var[[i_id]] <- current_names
        ts_br_num[[i_id]] <- 0
        next
      }

      #ts_br_var variable name list, the value is stop, the branch stops the branch

      #Time of branch at current moment
      current_ts_b_time <- current_ts_b_data[,2]

      #Status of branches at the current moment
      current_ts_b_status <- current_ts_b_data[,3]

      #All variables of the current time branch
      current_ts_b_variables <- subset(current_ts_b_data, select = -c(1,2,3))


      #Choose a meaningful variable from all variables, and audition returns cox$variables;cox$pvalue
      cox <- kmselectfeature(current_ts_b_time,current_ts_b_status,current_ts_b_variables,num_categories,p.value)

      #print("cox: ")
      #print(dim(cox$variables))
      #print(length(current_ts_b_time))


      #Determine the number of variables, if the number is 0, the current branch stops branching
      if (dim(cox$variables)[2]==0) {

        ts_br_var[[i_id]] <- "all"

        branch <- list(current_ts_b_id)

        b_varname = "all"

        ts_br_num[[i_id]] <- length(branch)

        branchID <- c(branchID,branch)
        branchID_var_value <- c(branchID_var_value,b_varname)

        next
      }


      #Determine the correlation of variables, group and select the main variables
      ###############################################################################################################
      #print('Determine the correlation of variables, group and select the main variables')

      #Returns the main variable of the current branch

      #degreeofcorrelation=0.7
      #

      mainvar <- selectvariables(current_ts_b_time,current_ts_b_status,cox,degreeofcorrelation)

      #cox.variables <- selectnoncorrvariables(current_ts_b_time,current_ts_b_status,mainvar)

      #print(paste('variables:',names(cox.variables)))

      #Selection of the most important variables

      #print(names(mainvar$variables))
      #print(dim(mainvar$variables))
      #print(names(mainvar$variables))
      #print(dim(mainvar$variables)[2])

      if(dim(mainvar$variables)[2]>1){
        mainvalue <- stepbyselectfeature(current_ts_b_time,current_ts_b_status, mainvar)
        curr_names <- names(mainvalue$variables)
      }else{
        curr_names <- names(mainvar$variables)
      }

      #print("curr_names")
      #print(curr_names)

      ts_br_var[[i_id]]  <- curr_names[1]

      #branch
      ts_b_v_status <- subset(cox$variables, select = curr_names[1])


      v_status <- unique(unlist(ts_b_v_status))
      v_status <- sort(v_status)


      #model <- cph(Surv(current_ts_b_time,current_ts_b_status)~ts_b_v_status)

      num = length(v_status)
      #print(paste('current branch num:',num,sep = ''))



      #go to branch
      branch <- lapply(1:num, branchfun,ts_b_v_status,v_status,current_ts_b_id)

      #####################################
      b_varname = list(v_status)
      #####################################

      ts_br_num[[i_id]] <- length(branch)

      branchID <- c(branchID,branch)
      branchID_var_value <- c(branchID_var_value,b_varname)

      #print(paste('current branch variable:',unlist(ts_br_var[[i_id]]),sep = ''))

    }

    #ts_br_var
    ts_var[[time_slice]] <- ts_br_var

    ts_num[[time_slice]] <- ts_br_num

    #print(paste('current branch num:',unlist(ts_br_num),sep = ''))

    ts_id[[time_slice]] <- branchID

    ts_value[[time_slice]]  <- branchID_var_value

  }

  survivalpathresults <- list(ts_var,ts_num,ts_id,ts_value)

  result  <- list()

  ## (result,time,status,tsid){
  result$data <- structureResult(survivalpathresults,time,status,tspatientid)

  result$tree <-  df2newick(result$data,innerlabel = T, period=rates)

  return(result)
}

#--------------------------------------------------------------------------------------------------------------
pvalue <- function(x,time,status,variables,v_names,num_categories){

  options(warn = 0)

  time = as.numeric(unlist(time))
  status = unlist(status)

  variable <- variables[x]
  #print(names(variable))
  variable <- unlist(variable)
  num = length(unique(variable))
  if (num > num_categories){

    stop(paste('The number of classes in ',x,'-th variable:',v_names[x],
               ' is ',num,',is higher than num_categories:',num_categories,sep = ''))
    return(1000)
  }

  if(num<=1){
    warning(paste('The number of classes in ',x,'-th variable:',v_names[x],
                  ' is ',num,',is lower than 2:',num_categories,sep = ''))
    return(1000)
  }

  else{
    S <-survdiff(Surv(time,status)~variable)
  }

  #print(S)
  p.value <- 1 - pchisq(S$chisq, length(S$n) - 1)
  #print(p.value)

  return(p.value)
}



#------------------------------------------------------------------------------------------------------------

kmselectfeature <- function(time,status,variables,num_categories,p.value){

  #print(paste("*************************kmselectfeature***************************"))

  #Gets the variable name
  v_names <- names(variables)


  #lapply:list  pvalue calculated p values
  var_s <- lapply(1:dim(variables)[2], pvalue, time,status,variables,v_names,num_categories)

  cox_v_pvalue <- unlist(var_s[var_s<p.value])

  #If the variable is num_categories category, the variable is cleared
  cox_v_names <- unlist(v_names[var_s<p.value])

  cox <- list()

  #print(cox$names)

  #Select meaningful variables
  cox$variables <- subset(variables, select = unlist(cox_v_names))

  #p values of meaningful variables
  cox$pvalue <- cox_v_pvalue

  return(cox)

}


#------------------------------------------------------------------------------------------------------------
#Batch cycles, each eliminating a major sample
selectvariables <- function(time,status,coxvar,degreeofcorrelation){

  #print(paste("*************************selectnoncorrvariables***************************"))
  #Calculation of correlation matrix
  corr <- rcorr(as.matrix(coxvar$variables))

  #r correlation coefficient
  r <- corr$r
  rshape <- dim(r)
  #Remove the polygonal line, its own correlation coefficient
  for(ri in 1:rshape[1]){
    r[(ri-1)*rshape[1]+ri] <- 0
  }

  #Remove null variables
  r[is.na(r)] <- 0
  FLAG <- r[which.max(r)]>= degreeofcorrelation

  # Remove the correlated variables
  if(FLAG) {
    coxvar <- pairwise(time, status,coxvar,degreeofcorrelation)
  }

  return(coxvar)
}

#--------------------------------------------------------------

#The correlation matrix is grouped, the most relevant pair of variables is
# selected, and a main variable is filtered by backward method
pairwise <- function(time, status,coxvar,degreeofcorrelation){

  #Calculation of correlation matrix
  corr <- rcorr(as.matrix(coxvar$variables))

  #r correlation coefficient
  r <- corr$r

  rshape <- dim(r)

  #Remove the polygonal line, its own correlation coefficient
  for(ri in 1:rshape[1]){
    r[(ri-1)*rshape[1]+ri] <- 0
  }

  #Remove null variables
  r[is.na(r)] <- 0

  wh <- dim(r)


  #Find the largest pair or pairs and return the corresponding id
  FLAG <- r[which.max(r)]>= degreeofcorrelation

  if(!FLAG){
    return(coxvar)
  }

  max_ind <- which(r==r[which.max(r)],arr.ind = T)
  #print(max_ind)

  #The largest, take all the largest front
  max_ind <- max_ind[1,] #integer

  max <- list()

  #Temporary two pairs of related variables
  max$variables <- subset(coxvar$variables, select = max_ind)

  max$pvalue <- coxvar$pvalue[max_ind]

  #Selection of main variables
  maxvar <- stepbyselectfeature(time, status, max)

  #print(maxvar$index)


  lenvar <- 1:length(max_ind)
  lenindex <- lenvar[-maxvar$index]

  index <- max_ind[lenindex]

  coxvar$variables <- subset(coxvar$variables,select = -index)
  coxvar$pvalue <- coxvar$pvalue[-index]

  return(coxvar)
}




#------------------------------------------------------------------------------------------------------------
#Backward selection of variables
stepbyselectfeature <- function(time,status,coxvar){
  #print(paste("*************************stepbyselectfeature***************************"))

  S <- try({
    #var_ <-data.frame(coxvar$variables)
    S<-coxph(Surv(time,status)~.,coxvar$variables,singular.ok = T)
  }
  ,silent=FALSE)


  # Determines whether the expression in the try statement of the current loop is running correctly
  if('try-error' %in% class(S))
  {
    max <- list()
    max$index <- which.min(coxvar$pvalue)
    max$variables <- subset(coxvar$variables,select = max$index)
    return(max)
  }

  S_step <- step(S, direction = 'backward',trace = 0)

  step_aic <- drop1(eval(S_step))


  max <- list()

  row_names <- row.names(step_aic)
  #print(row_names)

  varname <- names(coxvar$variables)
  #print(step_aic[,"AIC"])
  max$index <- which.max(step_aic[,"AIC"])
  #max$variables <- subset(coxvar$variables,select = max$index)
  max$names <- row_names[max$index]

  max$index <- which(varname==max$names)
  #print(max$index)
  max$variables <- subset(coxvar$variables,select = max$index)
  #print("*************************stepbyselectfeature end***************************")
  #max$

  return(max)

}


#------------------------------------------------------------------------------------------------------------
isbranch <- function(c_ts_b_id,minsample=15){
  if (dim(c_ts_b_id)[1]<minsample){
    current_names <- 'STOP'
  }else{
    return(NULL)
  }
}


#------------------------------------------------------------------------------------------------------------
#branch
branchfun <- function(x,status,value_status,id){
  status_id <- data.frame(status,id)

  index <- value_status[x]
  #print(paste("condition:",x,sep = ""))
  #print(index)

  b <- status_id[status_id[1]==index,2]
  b <- data.frame(b)

  #print(dim(b))
  return(b)
}


library(haven)
library(ggtree)
# Converts list-structured data to a data format that can be converted to-structured data
structureResult <- function(result,time,status,tsid){

  ts_var <- result[[1]]
  ts_num <- result[[2]]
  ts_id <- result[[3]]
  ts_varname <- result[[4]]

  ID = rbind(ts_id[[1]][[1]],ts_id[[1]][[2]])

  ##get ostime status by id

  new_data <- cbind(tsid[[1]],time[[1]],status[[1]])
  names(new_data) <- c("ID","time","status")

  names(ID) <- "ID"

  loc = match(ID$ID,new_data$ID)
  new_data <- new_data[loc,]

  # All branches and sort

  #Full time node of flashback cycle
  #Number of required categories less than 10

  Class =  rep("0", times=dim(new_data)[1])

  for (ts in length(ts_var):1){

    tp_var = vector("list", 0)
    tp_varname = vector("list", 0)


    #Class[idi] = paste("0",Class[idi],sep = "")

    #Grading each id of each time node
    for(idi in 1:dim(new_data)[1]){
      id = new_data$ID[idi]
      #print(id)
      #print(idi)


      tp_var = c(tp_var,list("None"))
      tp_varname = c(tp_varname,list(list("None")))


      isexist = FALSE
      #br: Each branch variable
      #print(paste("length:",length(ts_var[[ts]])))

      new_ts_var_ = ts_var[[ts]]
      new_ts_var_ = new_ts_var_[new_ts_var_!="STOP"]

      brni = 0
      if (length(new_ts_var_)==0){
        next
      }
      for (br in 1:length(new_ts_var_)){
        #print(paste("br:",br))

        if (isexist){break}

        #brn: Number of branches per branch variable
        #print(paste(ts,br,ts_varname[[ts]][[br]]))
        for (brn in 1:length(ts_varname[[ts]][[br]])){

          brni = brni+1
          #Determine whether id branch
          #print(paste("ts_brn:",ts,brni))
          #print(ts_id[[ts]][[brni]])

          isexist = id %in% ts_id[[ts]][[brni]][,1]
          if (isexist){


            tp_var[idi] <- new_ts_var_[[br]]
            tp_varname[idi] <- ts_varname[[ts]][[br]][[brn]]

            Class[idi] = paste(as.character(brn),substr(Class[idi],2,nchar(Class[idi])),sep = "")

            break

          }

        }

      }
      #print("----------------------------------------")
      #print(tp_var[1:10])
      #print(tp_varname[1:10])
      #print(Class[1:10])

    }
    #print("iiiiiiiiiiiiiiiiii")
    #print(length(tp_var))
    #print(length(tp_varname))
    assign(paste("time_slices_",ts,"_variable"),tp_var)
    assign(paste("time_slices_",ts,"_varname"),tp_varname)

    df1 = as.data.frame(matrix(unlist(tp_var),  byrow=1))
    names(df1) = paste("time_slices_",ts,"_variable",sep = "")

    df2 = as.data.frame(matrix(unlist(tp_varname),  byrow=1))
    names(df2) = paste("time_slices_",ts,"_varname",sep = "")


    new_data = cbind(new_data,df1,df2)

  }

  new_data = cbind(new_data,Class)
}


#------------------------------------------------------------------------------------------------------------

## recursion function
traverse <- function(data,a,i,innerl,survivaltime,period){
  res = list()

  desc <- NULL
  size <- NULL
  survival <- NULL
  survivalrate <- NULL
  if(i >4){
    df2 = data[which(as.character(paste(data[,i+1],data[,i+2],sep = "_"))==a),]
    alevelinner <- as.character(unique(paste(df2[,i-1],df2[,i],sep = "_")))
    #print(alevelinner)

    #drop NA data,no show
    if ("None_None" %in% alevelinner) {
      alevelinner <- alevelinner[-which(alevelinner=="None_None")]
      #print(alevelinner)
    }

    desc <- NULL
    size <- NULL
    survival <- NULL
    survivalrate <- NULL

    #
    new_df2 = df2
    #new_df2$time = new_df2$time-survivaltime
    if(length(alevelinner) == 1) {
      il <- NULL; if(innerl==TRUE) il <- a
      if ("all_all" %in% alevelinner){
        (newickout <- paste("(all_all:1)",il,":",1,sep=""))

        (size <-  paste("(",dim(df2)[1],":1)",dim(df2)[1],":",1,sep="") )
        fit <- surv_fit(Surv(time,status)~1,data=df2)
        median <- round(surv_median(fit)$median,1)
        #print(surv_median(fit))
        (survival <-  paste("(",median,":1)",median,":",1,sep="") )
        rate <- round(fit$surv[fit$time>period][1],3)
        (survivalrate <-  paste("(",rate,":1)",rate,":",1,sep="") )
      }
      #else if ( "None_None"  %in% alevelinner){(newickout <- paste(a,":",0.1,sep = "")) }
      else{
        r <- traverse(new_df2,alevelinner,i-2,innerl,survivaltime,period)
        (newickout <- r$newickout) #traverse(df2,alevelinner,i-2,innerl)
        (size <- r$size)
        (survival <- r$survival)
        (survivalrate <- r$survivalrate)

      }
    }
    else if(length(alevelinner) > 1) {


      # for  branch
      for(b in alevelinner){
        r <- traverse(new_df2,b,i-2,innerl,survivaltime,period)
        desc <- c(desc,r$newickout) #c(desc,traverse(df2,b,i-2,innerl))
        size <- c(size,r$size)
        survival <- c(survival,r$survival)
        survivalrate <- c(survivalrate,r$survivalrate)
      }
      il <- NULL; if(innerl==TRUE) il <- a
      #(newickout <- paste("(",paste(desc,collapse=","),")",il,":",i*0.05,sep=""))
      (newickout <- paste("(",paste(desc,collapse=","),")",il,":",1,sep=""))

      (size <- paste("(",paste(size,collapse=","),")",dim(df2)[1],":",1,sep=""))

      fit <- surv_fit(Surv(time,status)~1,data=df2)
      median <- round(surv_median(fit)$median,1)
      #print(surv_median(fit))
      (survival <- paste("(",paste(survival,collapse=","),")",median,":",1,sep=""))
      rate <- round(fit$surv[fit$time>period][1],3)
      (survivalrate <- paste("(",paste(survivalrate,collapse=","),")",rate,":",1,sep=""))
      #(size <- c(size,dim(df2)[1]))
    }
    else{
      (newickout <- paste(a,":",1,sep = ""))
      (size <- paste(dim(df2)[1],":",1,sep = ""))
      fit <- surv_fit(Surv(time,status)~1,data=df2)
      median <- round(surv_median(fit)$median,1)
      #print(surv_median(fit))
      (survival <- paste(median,":",1,sep = ""))
      rate <- round(fit$surv[fit$time>period][1],3)
      (survivalrate <- paste(rate,":",1,sep = ""))
      #(size <- c(size,dim(df2)[1]))
    }

  }
  #else { (newickout <- paste(a,":",0.1,sep = "")) }
  else {
    (newickout <- paste(a,":",1,sep = ""))
    (size <- paste(dim(data)[1],":",1,sep = ""))

    fit <- surv_fit(Surv(time,status)~1,data=data)
    median <- round(surv_median(fit)$median,1)
    #print(surv_median(fit))
    (survival <- paste(median,":",1,sep = ""))
    rate <- round(fit$surv[fit$time>period][1],3)
    (survivalrate <- paste(rate,":",1,sep = ""))
    #(size <- c(size,dim(data)[1]))
  }

  (res$newickout <- newickout)
  (res$size <- size)
  (res$survival <- survival)
  (res$survivalrate <- survivalrate)
  res
}


#result$data,innerlabel = T, period=365
## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE,survivaltime=3,period=365){
  df <- data.frame(df)
  len <- dim(df)[2]+1
  res <- list()
  alevel <- as.character(unique(paste(df[,len-3],df[,len-2],sep = "_")))
  newick <- NULL
  size <- NULL
  survival <- NULL
  survivalrate <- NULL

  # survival time period reduce
  new_df = df
  #new_df$time = new_df$time-survivaltime
  for(x in alevel){
    r <- traverse(new_df,x,len-4,innerlabel,survivaltime,period)
    newick <- c(newick,r$newickout)
    size <- c(size,r$size)
    survival <- c(survival,r$survival)
    survivalrate <- c(survivalrate,r$survivalrate)
  }

  (newick <- paste("(",paste(newick,collapse=","),")ALL;",sep=""))
  (size <- paste("(",paste(size,collapse=","),")",dim(df)[1],";",sep=""))
  fit <- surv_fit(Surv(time,status)~1,data=df)
  median <- round(surv_median(fit)$median,1)
  #print(surv_median(fit))
  rate <- round(fit$surv[fit$time>period][1],3)
  (survival <- paste("(",paste(survival,collapse=","),")",median,";",sep=""))
  (survivalrate <- paste("(",paste(survivalrate,collapse=","),")",rate,";",sep=""))
  res$newick <- newick
  res$size <- size
  res$survival <- survival
  res$survivalrate <- survivalrate


  mytree <- read.tree(text=res$newick)
  sizetree <- read.tree(text=res$size)
  survivaltree <- read.tree(text=res$survival)
  survivalratetree <- read.tree(text=res$survivalrate)
  #plot(mytree,family='STXihei')

  x = as_tibble(mytree)
  sizex = as_tibble(sizetree)
  survivalx = as_tibble(survivaltree)
  survivalratex = as_tibble(survivalratetree)
  #x[which(x$label!="" | x$branch.length==NaN),3]=0.1

  x$size = as.numeric(sizex$label)
  x$survival = as.numeric(survivalx$label)
  x$survivalrate = as.numeric(survivalratex$label)

  mytree1 = as.treedata(x)
}

