% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareTreatmentPlans.R
\name{compareTreatmentPlans}
\alias{compareTreatmentPlans}
\title{Compare and Draw the KM curve of specified treatment plan or exposure in selected nodes}
\usage{
compareTreatmentPlans(
df,
treepoints,
mytree,
source,
treatment
)
}
\arguments{
\item{df}{"data" in the return result of the \code{survivalpath()} function}

\item{treepoints}{list object;Specify the node for drawing the KM curve, which is displayed in the survival path graphs}

\item{mytree}{"tree" in the return result of the \code{survivalpath()} function}

\item{source}{Data.frame of time slice data, which could be returned by \code{timedivision()}}

\item{treatment}{Factor variable in the source data.frame. This argument is to specify the intervention or exposure that of interest at a specific node.}
}
\value{
No return value.
}
\description{
Based on the survival tree, specify the node of interest and the treatment methods, draw survival curves to evaluate the impact of treatments or exposure.
}
\details{
The function creates survival curves of specified treatment plan or exposure in selected nodes. The results should be interpreted with caution
as the effect of covariates have not been adjusted.
}
\examples{
library(dplyr)
data("DTSDHCC")
id = DTSDHCC$ID[!duplicated(DTSDHCC$ID)]
set.seed(123)
id = sample(id,500)
miniDTSDHCC <- DTSDHCC[DTSDHCC$ID \%in\% id,]
dataset = timedivision(miniDTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
resu <- generatorDTSD(dataset,periodindex="time_slice",IDindex="ID" ,timeindex="OStime_day",
 statusindex="Status_of_death",variable =c( "Age", "Amount.of.Hepatic.Lesions",
 "Largest.Diameter.of.Hepatic.Lesions",
 "New.Lesion","Vascular.Invasion" ,"Local.Lymph.Node.Metastasis",
 "Distant.Metastasis" , "Child_pugh_score" ,"AFP"),predict.time=365*1)
result <- survivalpath(resu,time_slices =9)

mytree <- result$tree

library(ggplot2)
library(ggtree)
ggtree(mytree, color="black",linetype=1,size=1.2,ladderize = TRUE )+
 theme_tree2() +
 geom_text2(aes(label=label),hjust=0.6, vjust=-0.6 ,size=3.0)+
 geom_text2(aes(label=paste(node,size,mytree@data$survival,mytree@data$survivalrate,sep = "/")),
 hjust=0.6, vjust=-1.85 ,size=3.0)+
 #geom_point2(aes(shape=isTip, color=isTip), size=mytree1@data$os/40)+
 geom_point2(aes(shape=isTip, color=isTip), size=mytree@data$size\%/\%200+1,show.legend=FALSE)+
 #guides(color=guide_legend(title="node name/sample number/Median survival time/Survival rate")) +
 labs(size= "Nitrogen",
      x = "TimePoints",
      y = "Survival",
      subtitle = "node_name/sample number/Median survival time/Survival rate",
      title = "Survival Tree") +
 theme(legend.title=element_blank(),legend.position = c(0.1,0.9))

#Comparing the efficacy of treatment methods by drawing survival curves
treepoints = c(14,20)
compareTreatmentPlans(result$data, treepoints,mytree,dataset,"Resection")

}
\seealso{
survminer
}
