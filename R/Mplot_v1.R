library(readxl)
library(data.tree)

data <- read_excel("C:/Users/Administrator/Desktop/SurvivalPath/0829/result.xlsx",sheet = "Sheet4")

data$pathString = paste("ALL",
                        paste(data$timepoints_1_variable,data$timepoints_1_varname,sep="_"),
                        paste(data$timepoints_2_variable,data$timepoints_2_varname,sep="_"),
                        paste(data$timepoints_3_variable,data$timepoints_3_varname,sep="_"),
                        paste(data$timepoints_4_variable,data$timepoints_4_varname,sep="_"),
                        paste(data$timepoints_5_variable,data$timepoints_5_varname,sep="_"),
                        paste(data$timepoints_6_variable,data$timepoints_6_varname,sep="_"),
                        paste(data$timepoints_7_variable,data$timepoints_7_varname,sep="_"),
                        paste(data$timepoints_8_variable,data$timepoints_8_varname,sep="_"),
                        paste(data$timepoints_9_variable,data$timepoints_9_varname,sep="_"),
                                                    sep = "/")

datatree = as.Node(data)

Prune(datatree, function(x) x$name != "None_None")

datatree1 = as.phylo.Node(datatree)
library(ggtree)
ggtree(datatree1)+   geom_tiplab() +
  geom_cladelabel(node=5, label="Some random clade",
                  color="red2", offset=-.3, align=TRUE) +
  geom_cladelabel(node=30, label="A different clade",
                  color="blue", offset=-.3, align=TRUE) +
  theme_tree2() +
  #  xlim(0, 1) +
  theme_tree()

