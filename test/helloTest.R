################################ Test ######################################
library(haven)
library(rms)
library(Hmisc)
library(survival)

dataset <- read_sav("D:/Code/R/SurvivalPath0225/Survivaldata/R_new.sav")

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

############################## 20190627ceshiresult.R####################################
library(haven)

dataset <- read_sav("D:/Code/R/SurvivalPath0225/Survivaldata/R_new.sav")

ts_var <- result[[1]]
ts_brnum <- result[[2]]
ts_brid <- result[[3]]
ts_br_varname <- result[[4]]

ID = rbind(ts_brid[[1]][[1]],ts_brid[[1]][[2]])

keys = c()
for (i in 1:length(ts_var)){
  for (keyi in c("Ts_","variable_","variable_value_")){
    keys = c(keys,paste(keyi,i,sep = ""))
    print(keys)
  }
}

for (ii in 1:dim(ID)[1]){

  id =  ID[ii,1]

  for (jj in 1:length(ts_var)){


  }

}


###############################StructuralData###########################


library(haven)

dataset <- read_sav("D:/Code/R/SurvivalPath0225/Survivaldata/R_new.sav")

ts_var <- result[[1]]
ts_brnum <- result[[2]]
ts_brid <- result[[3]]
ts_br_varname <- result[[4]]

ID = rbind(ts_brid[[1]][[1]],ts_brid[[1]][[2]])

## 通过id获取 ostime status
dup_dataset = dataset[which(dataset$TimeSlice==1),]
new_data <- subset(dup_dataset,select = c("病历号" ,"OS时间","是否死亡"))
names(new_data) <- c("ID","time","status")

names(ID) <- "ID"

loc = match(ID$ID,new_data$ID)
new_data <- new_data[loc,]

# 全部分支并排序

#倒叙循环全部时间节点
#要求类别个数小于10

Class =  rep("0", times=dim(new_data)[1])

for (ts in length(ts_var):1){

  tp_var = vector("list", 0)
  tp_varname = vector("list", 0)


  #Class[idi] = paste("0",Class[idi],sep = "")

  #对每一个时间节点的每一个ID分级
  for(idi in 1:dim(new_data)[1]){
    id = new_data$ID[idi]
    print(id)
    print(idi)


    tp_var = c(tp_var,list("None"))
    tp_varname = c(tp_varname,list(list("None")))


    isexist = FALSE
    #br: 每一个分支变量
    print(paste("length:",length(ts_var[[ts]])))

    new_ts_var_ = ts_var[[ts]]
    new_ts_var_ = new_ts_var_[new_ts_var_!="STOP"]

    brni = 0
    for (br in 1:length(new_ts_var_)){
      print(paste("br:",br))

      if (isexist){break}

      #brn: 每一个分支变量分支出来的个数
      print(paste(ts,br,ts_br_varname[[ts]][[br]]))
      for (brn in 1:length(ts_br_varname[[ts]][[br]])){

        brni = brni+1
        #判断 id 是否在分支
        print(paste("ts_brn:",ts,brni))
        print(ts_brid[[ts]][[brni]])

        isexist = id %in% ts_brid[[ts]][[brni]][,1]
        if (isexist){


          tp_var[idi] <- new_ts_var_[[br]]
          tp_varname[idi] <- ts_br_varname[[ts]][[br]][[brn]]

          Class[idi] = paste(as.character(brn),substr(Class[idi],2,nchar(Class[idi])),sep = "")

          break

        }

      }

    }
    print("----------------------------------------")
    print(tp_var[1:10])
    print(tp_varname[1:10])
    print(Class[1:10])

  }
  print("iiiiiiiiiiiiiiiiii")
  print(length(tp_var))
  print(length(tp_varname))
  assign(paste("timepoints_",ts,"_variable"),tp_var)
  assign(paste("timepoints_",ts,"_varname"),tp_varname)

  df1 = as.data.frame(matrix(unlist(tp_var),  byrow=1))
  names(df1) = paste("timepoints_",ts,"_variable",sep = "")

  df2 = as.data.frame(matrix(unlist(tp_varname),  byrow=1))
  names(df2) = paste("timepoints_",ts,"_varname",sep = "")


  new_data = cbind(new_data,df1,df2)

}

new_data = cbind(new_data,Class)

#write.csv(new_data,"D:/Code/R/SurvivalPath0225/Survivaldata/generator/mid3.csv", row.names = FALSE)

