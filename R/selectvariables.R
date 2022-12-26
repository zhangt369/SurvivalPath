#' kmselectfeature
#'
#' @param data: data.frame includes all survival data
#' @param variable:  Optional variables used for Cox
#'
#' @return Meaningful variables that can be used for Cox regression

library(survival)

#--------------------------------------------------------------------------------------------------------------

pandanvariable <- function(x,time,status,variables,v_names,num_categories){

  options(warn = 0)

  time = as.numeric(unlist(time))
  status = unlist(status)

  variable <- variables[x]
  print(names(variable))
  variable <- unlist(variable)
  num = length(unique(variable))
  if (num != num_categories){

    warning(paste('The number of classes in ',x,'-th variable:',v_names[x],
                  ' is ',num,', not equal to num_categories:',num_categories,sep = ''))
    return(-1)
  }

  else
    S <-survdiff(Surv(time,status)~variable)
  print(S)
  p.value <- 1 - pchisq(S$chisq, length(S$n) - 1)
  print(p.value)
  if (p.value >= 0.05)
    return(0)
  else
    return(x)
}



#------------------------------------------------------------------------------------------------------------


#' @return 选择出来的变量
kmselectfeature <- function(time,status,variables,num_categories){

  print(paste("*************************kmselectfeature***************************"))

  v_names <- names(variables)

  var_s <- lapply(1:dim(variables)[2], pandanvariable, time,status,variables,v_names,num_categories)

  cox_v_names <- unlist(var_s[var_s>0])

  cox <- list()

  cox$names <- v_names[cox_v_names]

  #print(cox$names)

  cox$variables <- subset(variables, select = unlist(cox_v_names))

  return(cox)

}


#------------------------------------------------------------------------------------------------------------
#判断是否有相关性比较高的变量
selectnoncorrvariables <- function(time,status,variables,corrvalue = 0.7){

  #print(paste("*************************selectnoncorrvariables***************************"))

  corr <- rcorr(as.matrix(variables))
  r <- corr$r
  print(r)


  r[is.na(r)] <- 0
  r[r==1.0] <- 0

  g <- groups(r,corrvalue)

  gvalue <- as.numeric(unlist(g))

  if (length(g)<1){
    return(variables)
  }
  print(g)
  gmain <- lapply(g,function(x,time,status){
    var <- subset(variables, select=x)
    main <- stepbyselectfeature(time,status,var)
    return(main$index)
  },time,status)


  gmain <- as.numeric(unlist(gmain))


  exist <- list()
  for (i in 1:length(g)){
    exist[[i]] <- g[[i]][gmain[[i]]]
  }

  #print(exist)

  exist <- as.numeric(unlist(exist))

  exist <- gvalue %in% exist

  notmaincorrvariables <- gvalue[!exist]

  #print(paste("notmaincorrvariables:",notmaincorrvariables))

  allvalueindex <- 1:dim(variables)[2]

  correxist <- allvalueindex %in% notmaincorrvariables

  #print(paste("correxist",correxist))

  noncorrvariablesindex <- allvalueindex[!correxist]

  #print(names(variables))
  #print(paste("noncorrvariablesindex:",noncorrvariablesindex))

  noncorrvariables <- subset(variables, select=noncorrvariablesindex)
  #print(names(noncorrvariables))

  return(noncorrvariables)

}

#--------------------------------------------------------------
#相关矩阵分组
groups <- function(r,corrvalue){
  wh <- dim(r)

  group <- list()
  allvalues <- c()
  ind <- 1

  for (i in 2:wh[2]-1){

    if (i %in% allvalues)
      next

    group[ind] <- c(i)

    for (j in (i+1):wh[2]){

      if (r[i,j]>corrvalue){
        group[[ind]] <- c(group[[ind]],j)
      }
    }

    allvalues <- c(allvalues,group[[ind]])

    if(length(group[[ind]])!=1){

      ind <- 1+ind
    }

    group[ind] <- NULL



  }
  #print(allvalues)
  #print(paste('length:',length(group)))
  #print(group)

  group

}




#------------------------------------------------------------------------------------------------------------
#后退法选择变量
stepbyselectfeature <- function(time,status,variables){
  #print(paste("*************************stepbyselectfeature***************************"))

  S <- cph(Surv(time,status)~.,variables,singular.ok = FALSE)
  #print(S)

  S_step <- step(S, direction = 'backward',trace = 0)

  step_aic <- drop1(eval(S_step))
  #print(step_aic)

  max <- list()

  row_names <- row.names(step_aic)

  varname <- names(variables)

  max$index <- which.max(step_aic[,"AIC"])

  max$names <- row_names[max$index]

  max$index <- which(varname==max$names)
  #print("*************************stepbyselectfeature end***************************")


  return(max)

}


#------------------------------------------------------------------------------------------------------------

#' @describeIn
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

  b <- status_id[status_id[1]==index,2]
  b <- data.frame(b)

  #print(dim(b))
  }

