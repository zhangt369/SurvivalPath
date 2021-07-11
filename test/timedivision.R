

library(readxl)
library(dplyr)



#dataset <- read_excel("C:/Users/Administrator/Desktop/SurvivalPath/时间序列逻辑清洗.xls",sheet="总数据库")
dataset <- read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/2020-10-03 新数据分类后.xlsx",sheet="Sheet1")


dataset[["住院日期"]] <- as.Date(dataset[["影像检查/住院日期"]])


test <- dataset %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期)

test_1 <- dataset %>% group_by(病历号) %>% filter(n() == 1) %>% arrange(病历号,住院日期)

test1 <- dataset %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)


while(dim(test)[1] != dim(test1)[1]){
  print("---------------------------------------------")
  print(dim(test)[1])

  test1 <- test %>% mutate(time3 = round(as.numeric((住院日期 - lag(住院日期) )/30))) %>% filter( is.na(time3) == TRUE | time3 ==1 ) %>% arrange(病历号,desc(住院日期)) %>% ungroup


  test <- test1 %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期)

}

print("test")
test2 <- test %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)

test1 <- test %>% group_by(病历号) %>% arrange(病历号,住院日期) %>% mutate(time4 = 1) %>% mutate(timenode = cumsum(time4)) %>% ungroup


test1 <- bind_rows(test_1,test1)


write.csv(test1, "D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/test2021.csv")




