library(ggplot2)
library(readxl)
library(data.tree)
library(ggtree)
library(dplyr)
library(treeio)
library(survival)
library(survminer)

data <- read.csv("D:/Code/R/SurvivalPath0225/Survivaldata/generator/mid.csv")

## recursion function
traverse <- function(data,a,i,innerl){
  res = list()
  print(i)
  print(a)
  desc <- NULL
  size <- NULL
  if(i >2){
    df2 = data[which(as.character(paste(data[,i+1],data[,i+2],sep = "_"))==a),]
    alevelinner <- as.character(unique(paste(df2[,i-1],df2[,i],sep = "_")))
    print(alevelinner)


    if ("None_None" %in% alevelinner) {
      alevelinner <- alevelinner[-which(alevelinner=="None_None")]
      print(alevelinner)
    }

    desc <- NULL
    size <- NULL
    if(length(alevelinner) == 1) {
      il <- NULL; if(innerl==TRUE) il <- a
      if ("all_all" %in% alevelinner){
        (newickout <- paste("(all_all:0.5)",il,":",0.5,sep=""))

        (size <-  paste("(",dim(df2)[1],":0.5)",dim(df2)[1],":",0.5,sep="") )
      }
      #else if ( "None_None"  %in% alevelinner){(newickout <- paste(a,":",0.1,sep = "")) }
      else{
        r <- traverse(df2,alevelinner,i-2,innerl)
        (newickout <- r$newickout) #traverse(df2,alevelinner,i-2,innerl)
        (size <- r$size)

      }
    }
    else if(length(alevelinner) > 1) {


      # loop brach
      for(b in alevelinner){
        r <- traverse(df2,b,i-2,innerl)
        desc <- c(desc,r$newickout) #c(desc,traverse(df2,b,i-2,innerl))
        size <- c(size,r$size)
      }
      il <- NULL; if(innerl==TRUE) il <- a
      #(newickout <- paste("(",paste(desc,collapse=","),")",il,":",i*0.05,sep=""))
      (newickout <- paste("(",paste(desc,collapse=","),")",il,":",0.5,sep=""))

      (size <- paste("(",paste(size,collapse=","),")",dim(df2)[1],":",0.5,sep=""))
      #(size <- c(size,dim(df2)[1]))
    }
    else{
      (newickout <- paste(a,":",0.3,sep = ""))
      (size <- paste(dim(df2)[1],":",0.3,sep = ""))
      #(size <- c(size,dim(df2)[1]))
    }

  }
  #else { (newickout <- paste(a,":",0.1,sep = "")) }
  else {
    (newickout <- paste(a,":",0.5,sep = ""))
    (size <- paste(dim(data)[1],":",0.5,sep = ""))
    #(size <- c(size,dim(data)[1]))
  }

  (res$newickout <- newickout)
  (res$size <- size)
  res
}



## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  res <- list()
  alevel <- as.character(unique(paste(df[,19],df[,20],sep = "_")))
  newick <- NULL
  size <- NULL
  for(x in alevel){
    r <- traverse(df,x,18,innerlabel)
    newick <- c(newick,r$newickout)
    size <- c(size,r$size)
  }
  (newick <- paste("(",paste(newick,collapse=","),")ALL;",sep=""))
  (size <- paste("(",paste(size,collapse=","),")",dim(df)[1],";",sep=""))
  res$newick <- newick
  res$size <- size
  res
}

df = data

mk = df2newick(df,innerlabel = T)

mytree <- read.tree(text=mk$newick)
sizetree <- read.tree(text=mk$size)
#plot(mytree,family='STXihei')

x = as_tibble(mytree)
sizex = as_tibble(sizetree)
#x[which(x$label!="" | x$branch.length==NaN),3]=0.1


x$size = as.numeric(sizex$label)

mytree1 = as.treedata(x)

p1 <- ggtree(mytree1, color="black",linetype=1,size=1.,ladderize = T, )+
  theme_tree2() +
  geom_text2(aes(label=label),hjust=0.6, vjust=-0.5 ,size=2.0)+
  geom_text2(aes(label=paste(size,sep = "/")),hjust=0.6, vjust=-1.65 ,size=2.4)+
  #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
  geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$size%/%35+1)+
  labs(size= "Nitrogen",
       x = "TimePoints",
       y = "Survival",
       title = "Survival Tree")

p1


