
library(haven)

dataset <- read_sav("C:/Users/Administrator/Desktop/SurvivalPath/肝癌R包制作标准数据_new.sav")

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
