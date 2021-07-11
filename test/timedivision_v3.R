

library(readxl)
library(dplyr)



#dataset <- read_excel("C:/Users/Administrator/Desktop/SurvivalPath/时间序列逻辑清洗.xls",sheet="总数据库")
dataset <- read_excel("D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/2020-10-03 新数据分类后.xlsx",sheet="Sheet1")


dataset[["住院日期"]] <- as.Date(dataset[["影像检查/住院日期"]])

#test <- dataset %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期)
#test0 <- test %>% slice(1) %>% ungroup # 按病历号和 住院日期排序
#test <- test %>% ungroup
#test0 <- test0 %>% ungroup
#a <-bind_rows(test,test0)

test <- dataset %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期) #6890

test_1 <- dataset %>% group_by(病历号) %>% filter(n() == 1) %>% arrange(病历号,住院日期) #85

test1 <- dataset %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期) #1146

#去重

print("----------------------dup-----------------------")
print(dim(test)[1])

test1 <- test %>% mutate(time2 = round(as.numeric(住院日期 - min(住院日期, na.rm = TRUE))/30))%>%
  mutate(time3 = time2 - lag(time2) ) %>%

  filter( is.na(time3) == TRUE | time3 >=1 ) %>%
  arrange(病历号,desc(住院日期)) %>%
  ungroup

print(dim(test1)[1])
#计算二次入院和前一次入院的时间差
test <- test1 %>% group_by(病历号) %>%   arrange(病历号,住院日期)



print("---------------------------------------------")
print(dim(test)[1])

test1 <- test %>% mutate(time4 = round(as.numeric(住院日期 - min(住院日期, na.rm = TRUE))/30))%>%
  mutate(time5 = time4 - lag(time4) ) %>%
  filter( is.na(time5) == TRUE | time5 ==1 ) %>%
  arrange(病历号,desc(住院日期)) %>%
  ungroup
#计算二次入院和前一次入院的时间差
test <- test1 %>% group_by(病历号)  %>% arrange(病历号,住院日期)




print("test")
test2 <- test %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)

test1 <- test %>% group_by(病历号) %>% arrange(病历号,住院日期) %>% mutate(time6 = 1) %>% mutate(timenode = cumsum(time6)) %>% ungroup

test1 <- bind_rows(test_1,test1)

write.csv(test1, "D:/Code/R/SurvivalPath0225/Survivaldata/SurvivalPath/test2.csv")

A <- test1 %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期) #1146
B <- test1 %>% group_by(病历号) %>% filter(n() > 1) %>% arrange(病历号,住院日期)
C <- B %>% group_by(病历号) %>% slice(1) %>% arrange(病历号,住院日期)
