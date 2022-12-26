
library(pROC)

classifydata <- function(time,status,timeslicedata,tspatientid){

  varnames <- names(timeslicedata[[1]])

  variables <- timeslicedata[[1]]

  for (v in 1:length(varnames)){

    varn <- varnames[v]

    variable <- variables[varn]

    #print(names(variable))

    variable <- unlist(variable)

    num = length(unique(variable))

    if (num > 1){

      roc <- roc(status[[1]][[1]], variable,main="Confidence intervals", percent=TRUE,ci=TRUE )

      ci <- ci(roc, of="thresholds", thresholds="best")

      thresholds <- attr(ci,"thresholds")[1]

      for (i in 1:length(timeslicedata)){

        # 重新转为data.frame

        timeslicedata[[i]][[v]] <- ifelse(timeslicedata[[i]][[v]] < thresholds,0,100)
      }

    }

  }

  return(timeslicedata)

}


autodividedtimenode <- function()


library(haven)

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

tsdata <- classifydata(time,status,tsdata,tsid)

