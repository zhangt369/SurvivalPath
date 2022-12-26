
#'
#' @param data: data.frame includes all survival data
#' @param variable:  Optional variables used for Cox
#'
#' @return Meaningful variables that can be used for Cox regression

library(survival)

#--------------------------------------------------------------------------------------------------------------
#' @description 计算单变量的p值，并判断变量类别是否为 num_categories，不是就warning
#' @param time 时间
#' @param status 状态
#' @param variables 单变量
#' @param v_names 单变量的名字
#' @param num_categories 离散变量的分类数
#' @return P.value，均为正 如果返回值为1000，意味该变量不是num_categories类 （是否需要区分1类和多类）
pvalue <- function(x,time,status,variables,v_names,num_categories){

  options(warn = 0)

  time = as.numeric(unlist(time))
  status = unlist(status)

  variable <- variables[x]
  #print(names(variable))
  variable <- unlist(variable)
  num = length(unique(variable))
  if (num != num_categories){

    warning(paste('The number of classes in ',x,'-th variable:',v_names[x],
                  ' is ',num,', not equal to num_categories:',num_categories,sep = ''))
    return(1000)
  }

  else
    S <-survdiff(Surv(time,status)~variable)
  #print(S)
  p.value <- 1 - pchisq(S$chisq, length(S$n) - 1)
  #print(p.value)

  return(p.value)
}



#------------------------------------------------------------------------------------------------------------

#' @description 对全部变量，选择出有意义的变量，logrank 检验 p值小于0.05
#' @param time 时间
#' @param status 事件标志
#' @param variables 待验证的全部变量
#' @param num_categories 变量是几分类数据
#' @return cox, cox$names 变量名字；
#' cox$variables 变量，数据类型为data.frame; cox$pvalue 每个变量的P值
kmselectfeature <- function(time,status,variables,num_categories){

  #print(paste("*************************kmselectfeature***************************"))

  v_names <- names(variables) #获取变量名

  #lapply:list， pvalue计算P值
  var_s <- lapply(1:dim(variables)[2], pvalue, time,status,variables,v_names,num_categories)

  cox_v_pvalue <- unlist(var_s[var_s<0.05])
  cox_v_names <- unlist(v_names[var_s<0.05]) #如果变量是num_categories个类别，该变量会被清除

  print(cox_v_pvalue)
  print(cox_v_names)

  cox <- list()

  #print(cox$names)

  #选择有意义变量
  cox$variables <- subset(variables, select = unlist(cox_v_names))

  #有意义变量的P值
  cox$pvalue <- cox_v_pvalue

  return(cox)

}


#------------------------------------------------------------------------------------------------------------
#一批批循环，每个循环剔除一个主要样本

selectvariables <- function(time,status,coxvar){

  #print(paste("*************************selectnoncorrvariables***************************"))


  while (dim(coxvar$variables)[2]>2) {

    coxvar <- pairwise(time, status,coxvar)
  }


  return(coxvar)
}

#--------------------------------------------------------------

#相关矩阵分组,选择最相关的一对变量，使用后退法筛选一个主要变量
pairwise <- function(time, status,coxvar){

  #print("pairwise-------------->")

  #计算相关矩阵
  corr <- rcorr(as.matrix(coxvar$variables))

  #r相关系数
  r <- corr$r


  rshape <- dim(r)


  #去掉多角线，自己相关的系数
  for(ri in 1:rshape[1]){
    r[(ri-1)*rshape[1]+ri] <- 0
  }

  #去掉空变量
  r[is.na(r)] <- 0

  #print(r)

  wh <- dim(r)


  #找出最大的一对或多对，返回对应ID
  max_ind <- which(r==r[which.max(r)])
  #print(max_ind)


  #最大,取全部最大的最前面一个
  max_ind <- max_ind[1] #integer

  max_ind <- c((max_ind-1) %% wh[1] + 1, max_ind %/% wh[1]+1)

  #print(max_ind)


  max <- list()


  #临时两个对相关变量
  max$variables <- subset(coxvar$variables, select = max_ind)

  max$pvalue <- coxvar$pvalue[max_ind]

  #后退法选择主要变量
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
#后退法选择变量
stepbyselectfeature <- function(time,status,coxvar){
  #print(paste("*************************stepbyselectfeature***************************"))


  S <- try({

    S<-coxph(Surv(time,status)~.,coxvar$variables,singular.ok = FALSE)
    }
    ,silent=FALSE)

  #print(S)
  if('try-error' %in% class(S))            # 判断当前循环的try语句中的表达式是否运行正确
  {
    max <- list()
    max$index <- which.min(coxvar$pvalue)
    max$variables <- subset(coxvar$variables,select = max$index)
    return(max)
  }


  #print(S)

  S_step <- step(S, direction = 'backward',trace = 0)

  step_aic <- drop1(eval(S_step))
  #print(step_aic)

  max <- list()

  row_names <- row.names(step_aic)
  #print(row_names)

  varname <- names(coxvar$variables)
  #print(varname)

  max$index <- which.max(step_aic[,"AIC"])
  #print(max$index)

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

#' @description 本函数判断该分支数中样本个数是否达到（单变量分支最低要求）个数12
#' @return 如果达不到分支数，返回“STOP”
isbranch <- function(c_ts_b_id){
  if (dim(c_ts_b_id)[1]<12){
    current_names <- 'STOP'
  }
}


#------------------------------------------------------------------------------------------------------------
#branch
branchfun <- function(x,status,value_status,id){
  status_id <- data.frame(status,id)

  index <- value_status[x]
  print(paste("condition",x,":",sep = ""))
  print(index)

  b <- status_id[status_id[1]==index,2]
  b <- data.frame(b)

  #print(dim(b))
}
