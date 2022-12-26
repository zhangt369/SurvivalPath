

library(readxl)
library(dplyr)



dataset <- read_excel("C:/Users/Administrator/Desktop/SurvivalPath/时间序列逻辑清洗.xls",sheet="总数据库")

k=10

dataset[["住院日期"]] <- as.Date(dataset[["住院日期"]])

#test <- dataset %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期) 
#test0 <- test %>% slice(1) %>% ungroup # 按病历号和 住院日期排序
#test <- test %>% ungroup
#test0 <- test0 %>% ungroup
#a <-bind_rows(test,test0)

test <- dataset %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期) 

test_1 <- dataset %>% group_by(病历号) %>% filter(n() == 1) %>% arrange(病历号,住院日期) 

test1 <- dataset %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)


while(dim(test)[1] != dim(test1)[1]){
  print("---------------------------------------------")
  print(dim(test)[1])
  
  test1 <- test %>% mutate(time3 = round(as.numeric((住院日期 - lag(住院日期) )/30))) %>% filter( is.na(time3) == TRUE | time3 ==1 ) %>% arrange(病历号,desc(住院日期)) %>% ungroup 
  #计算二次入院和前一次入院的时间差
  
  test <- test1 %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期) 
  
}

print("test")
test2 <- test %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)

test1 <- test %>% group_by(病历号) %>% arrange(病历号,住院日期) %>% mutate(time4 = 1) %>% mutate(timenode = cumsum(time4)) %>% ungroup 


test1 <- bind_rows(test_1,test1)


write.xlsx(test1, "C:/Users/Administrator/Desktop/SurvivalPath/test1.xlsx",sheetName = "data",append = TRUE)




