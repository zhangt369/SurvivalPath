
#'
#' @param data: data.frame includes all survival data
#' @param variable:  Optional variables used for Cox
#'
#' @return Meaningful variables that can be used for Cox regression

library(survival)

#--------------------------------------------------------------------------------------------------------------
#' @description To calculate the p value of a single variable and to determine whether the variable category is num_categories, warning or not
#' @param time
#' @param status
#' @param variables
#' @param v_names variable name
#' @param num_categories Classification of discrete variables
#' @return P.value£¬All positive if the return value is 1000, meaning that the variable is not a num_categories class (need to distinguish between class 1 and class multiple)
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

#' @description To all variables, select meaningful variables, logrank test p value is less than 0.05
#' @param time time
#' @param status event flag
#' @param variables All variables to be validated
#' @param num_categories A variable is a sort of data
#' @return cox, cox$names  variable name
#' cox$variables Variable, data type is the p value data.frame;cox$pvalue each variable
kmselectfeature <- function(time,status,variables,num_categories){

  #print(paste("*************************kmselectfeature***************************"))

  #Gets the variable name
  v_names <- names(variables)

  #lapply:list  pvalue calculated p values
  var_s <- lapply(1:dim(variables)[2], pvalue, time,status,variables,v_names,num_categories)

  cox_v_pvalue <- unlist(var_s[var_s<0.05])

  #If the variable is num_categories category, the variable is cleared
  cox_v_names <- unlist(v_names[var_s<0.05])

  print(cox_v_pvalue)
  print(cox_v_names)

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

selectvariables <- function(time,status,coxvar){

  #print(paste("*************************selectnoncorrvariables***************************"))


  while (dim(coxvar$variables)[2]>2) {

    coxvar <- pairwise(time, status,coxvar)
  }


  return(coxvar)
}

#--------------------------------------------------------------

#The correlation matrix is grouped, the most relevant pair of variables is
# selected, and a main variable is filtered by backward method
pairwise <- function(time, status,coxvar){

  #print("pairwise-------------->")

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

  #print(r)

  wh <- dim(r)


  #Find the largest pair or pairs and return the corresponding id
  max_ind <- which(r==r[which.max(r)])
  #print(max_ind)


  #The largest, take all the largest front
  max_ind <- max_ind[1] #integer

  max_ind <- c((max_ind-1) %% wh[1] + 1, max_ind %/% wh[1]+1)

  #print(max_ind)


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

    S<-coxph(Surv(time,status)~.,coxvar$variables,singular.ok = FALSE)
    }
    ,silent=FALSE)

  #print(S)
  # Determines whether the expression in the try statement of the current loop is running correctly
  if('try-error' %in% class(S))
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

#' @description This function determines whether the number of samples in the branch number meets (the minimum
#' requirement for univariate branches)12
#' @return current_names If the number of branches is not reached, return stop"
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
