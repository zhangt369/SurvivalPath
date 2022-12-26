
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#'@title  survivalpath
#'@description build survival path
#'@usage survivalpath()
#'@param survivaltime survival time
#' survivalresult survival result
#' timepotin the number of timepoint
#'@details for survivalpath
#'@return mm all result
#'@author zs
#'@seealso surv
#'@example example test x=y+7
#'@importFrom rms cph
#'@export
survivalpath <- function(time,status,timeslicedata,tspatientid,num_categories=2) {
  l_time <- length(time)
  l_status <- length(status)
  l_timesl <- length(timeslicedata)
  l_id <- length(tspatientid)


  ######################################################################################################
  #输入数据格式判断

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

  #数据预处理，类别不为2的变量变换为2类变量
  #timeslicedata <- classifydata(time,status,timeslicedata,tspatientid)

  #输入数据时间节点个数是否一致判断
  if (l_time!=l_status | l_time!= l_timesl | l_time != l_id){
         stop(sprintf('The length of time,status,timeslicedata and tspatientid are:%d, %d, %d, %d,
                                but they should be the same length.',l_time,l_status,l_timesl,l_id))
     }

  #要求输入患者唯一编码ID:
  if (l_id==0){
    stop('The input can\'t be empty')
  }


  #要求ID唯一，判断
  patientID <- tspatientid[[1]]


  patientID_a <-unique(patientID[duplicated(patientID),])

  lenpid <- length(as.numeric(unlist(patientID_a)))
  if (lenpid!=0){
    stop(paste('The patient ID repeat. Every patient ID should be unique:',patientID_a))
  }
  ###############################################################################################################
  #id in patientID


  #时间节点个数，等于输入数据的长度
  timepoints <- 9

  #全部患者ID
  branchID <- list(tspatientid[[1]])

  #初始化
  ts_var <- vector("list", timepoints)
  ts_brnum <- vector("list", timepoints)
  ts_brid <- vector("list", timepoints)
  ts_br_varname <- vector("list", timepoints)


  for (timepoint in 1:timepoints){
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print(paste('timepoint:',timepoint,sep = ''))

    previousbranchID <- branchID  #list


    #当前时间节点的全部样本时间
    current_time <- time[[timepoint]]

    #当前时间节点的全部样本状态
    current_status <- status[[timepoint]]

    #当前时间节点的全部样本ID
    current_id <- tspatientid[[timepoint]]

    #当前时间节点的全部样本变量：数据帧
    current_variables <- timeslicedata[[timepoint]]

    #打印变量名称
    print(paste("datavar:",names(current_variables)))


    #判断四个数据的维度是否一致
    if (dim(current_time)[1]!=dim(current_status)[1] | dim(current_time)[1]!= dim(current_variables)[1] |
        dim(current_time)[1] != dim(current_id)[1]){stop(sprintf('The number of the %d-th timeslice of time,status,
                        timeslicedata and tspatientid are:%d, %d, %d, %d, but they should be the same ',timepoint,
                                                                 l_time,l_status,l_timesl,l_id))
    }

    #合并数据：id,time,status,varibales
    current_ts_data <- data.frame(current_id,current_time,current_status,current_variables)

    #维度判断
    if (dim(current_time)[2]!=1){
      stop(sprintf('the %d-th timeslice of survival time should be unique, but the dim is:%d',2,dim(current_time)[2]))
      }

    if (dim(current_status)[2]!=1){
      stop(sprintf('the %d-th timeslice of survival status should be unique, but the dim is:%d',2,dim(current_status)[2]))
    }

    if (dim(current_id)[2]!=1){
      stop(sprintf('the %d-th timeslice of patient id should be unique, but the dim is:%d',2,dim(current_id)[2]))
    }


    ############################################# 前一分支分支 ##########################################################

    #前一时刻保留下来的分支数量：previousbranchID-》每一分支的ID,list()，data.frame,brnum:分支数
    brnum <- length(previousbranchID)


    #当前时刻，当前分支的变量名
    ts_br_var <- vector('list',brnum)
    #ts_br_v_names <- vector('list',brnum)


    # 当前时刻，当前分支的分出的分支数
    ts_br_num <- vector('list',brnum)


    # ？？？？？？？？？？？？？？？？
    branchID <- list()
    branchID_varname <- list()

    for (i_id in 1:brnum) {
      print('########################################')

      print(paste('current branch:',i_id,sep = ''))

      #ts_br_v_names[[i_id]] <- paste('ts',timepoint,'_br',i_id)

      #前一时间节点留存到当前分支的ID
      br_id <- data.frame(previousbranchID[[i_id]])

      nam <- names(current_ts_data)
      print(nam[1])

      names(br_id) <- nam[1]
      print(names(br_id)[1])

      #当前生存数据的维度

      print("pre branch:")
      print(dim(current_ts_data))


      #提取current_ts_data在br_id中有的数据行，按名称“nam[1]”
      current_ts_b_data <- merge(current_ts_data,br_id,by= nam[1],all=FALSE)

      print("after merge, current branch:")
      print(dim(current_ts_b_data))


      #当前时间节点分支的ID
      current_ts_b_id <- subset(current_ts_b_data, select = nam[1])

      current_names <- isbranch(current_ts_b_id) #空值 或者 STOP

      #判断是否继续划分
      if (!is.null(current_names)) {
        ts_br_var[[i_id]] <- current_names
        ts_br_num[[i_id]] <- 0
        next
      }

      #ts_br_var 变量名list中，值为STOP，该分支停止分支

      #当前时刻分支的时间
      current_ts_b_time <- current_ts_b_data[,2]

      #当前时刻分支的状态
      current_ts_b_status <- current_ts_b_data[,3]

      #当前时刻分支的全部变量
      current_ts_b_variables <- subset(current_ts_b_data, select = -c(1,2,3))


      #从全部变量中选择有意义变量,海选 返回cox$variables; cox$pvalue
      cox <- kmselectfeature(current_ts_b_time,current_ts_b_status,current_ts_b_variables,num_categories)

      #print("cox: ")
      #print(dim(cox$variables))
      #print(length(current_ts_b_time))


      #判断变量个数，如果个数为0，当前分支停止分支
      if (dim(cox$variables)[2]==0) {

        ts_br_var[[i_id]] <- "all"

        branch <- list(current_ts_b_id)

        b_varname = "all"

        ts_br_num[[i_id]] <- length(branch)

        branchID <- c(branchID,branch)
        branchID_varname <- c(branchID_varname,b_varname)

        next
      }


      #判断变量的相关性，并分组，选出主要变量
      ###############################################################################################################
      print('判断变量的相关性，并分组，选出主要变量')

      #返回当前分支的主要变量
      mainvar <- selectvariables(current_ts_b_time,current_ts_b_status,cox)

      #cox.variables <- selectnoncorrvariables(current_ts_b_time,current_ts_b_status,mainvar)

      #print(paste('variables:',names(cox.variables)))

      #后退法选择最主要变量

      print(names(mainvar$variables))
      print(dim(mainvar$variables))
      print(mainvar$pvalue)

      mainvalue <- stepbyselectfeature(current_ts_b_time,current_ts_b_status, mainvar)
      #print(mainvalue)

      curr_names <- names(mainvalue$variables)

      #print("curr_names")
      #print(curr_names)

      ts_br_var[[i_id]]  <- curr_names[1]
      #分支
      ts_b_v_status <- subset(cox$variables, select = curr_names)


      v_status <- unique(unlist(ts_b_v_status))


      num = length(v_status)
      #print(paste('current branch num:',num,sep = ''))


      #进行分支
      branch <- lapply(1:num, branchfun,ts_b_v_status,v_status,current_ts_b_id)

      #####################################
      b_varname = list(v_status)
      #####################################

      ts_br_num[[i_id]] <- length(branch)

      branchID <- c(branchID,branch)
      branchID_varname <- c(branchID_varname,b_varname)

      print(paste('current branch variable:',unlist(ts_br_var[[i_id]]),sep = ''))

    }

    #ts_br_var
    ts_var[[timepoint]] <- ts_br_var

    ts_brnum[[timepoint]] <- ts_br_num

    print(paste('current branch num:',unlist(ts_br_num),sep = ''))

    ts_brid[[timepoint]] <- branchID

    ts_br_varname[[timepoint]]  <- branchID_varname

  }

  survivalpathresults <- list(ts_var,ts_brnum,ts_brid,ts_br_varname)

}


################################ Test ######################################
library(haven)
library(rms)
library(Hmisc)
library(survival)

dataset <- read_sav("C:/Users/Administrator/Desktop/SurvivalPath/肝癌R包制作标准数据_new.sav")

# one id to one time and status ----------------------------------------------------------------
time <- list()
status <- list()
tsdata <- list()
tsid <- list()
for (i in 1:9){
  data <- dataset[dataset['TimeSlice']==i,]
  time <- c(time,list(data['OS时间']))
  status <- c(status,list(data['是否死亡']))
  tsid <- c(tsid,list(data['病历号']))
  c_data <- subset(data, select = c('AFP25','AFP200','AFP400','ChildBCvA','是否单个病灶','肝脏病灶最大直径50',
                                    '肝脏病灶最大直径70','肝脏病灶最大直径100','是否无活性病灶转','是否血管侵犯'
                                    ,'是否腹水','是否远处转移',"是否23病灶最大小于3cm",'是否23个最大大于3cm','是否4个结节或者23个结节最大大于3cm转'))

  tsdata <- c(tsdata,list(c_data))
}


result <- survivalpath(time,status,tsdata,tsid)

