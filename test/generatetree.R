library(ggplot2)
library(readxl)
library(data.tree)
library(ggtree)
library(dplyr)
library(treeio)
library(survival)
library(survminer)


## recursion function
traverse <- function(data,a,i,innerl,survivaltime,period){
  res = list()
  print(i)
  print(a)
  desc <- NULL
  size <- NULL
  survival <- NULL
  survivalrate <- NULL
  if(i >4){
    df2 = data[which(as.character(paste(data[,i+1],data[,i+2],sep = "_"))==a),]
    alevelinner <- as.character(unique(paste(df2[,i-1],df2[,i],sep = "_")))
    print(alevelinner)

    #drop NA data,no show
    if ("None_None" %in% alevelinner) {
      alevelinner <- alevelinner[-which(alevelinner=="None_None")]
      print(alevelinner)
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
        print(surv_median(fit))
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
      print(surv_median(fit))
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
      print(surv_median(fit))
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
    print(surv_median(fit))
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



## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE,survivaltime=3,period=12){
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
  print(surv_median(fit))
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
