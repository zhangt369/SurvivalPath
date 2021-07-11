loc<-"http://archive.ics.uci.edu/ml/machine-learning-databases/"
ds<-"breast-cancer-wisconsin/breast-cancer-wisconsin.data"
url<-paste(loc,ds,sep="")
data<-read.table(url,sep=",",header=F,na.strings="?")
names(data)<-c("编号","肿块厚度","肿块大小","肿块形状","边缘黏附","单个表皮细胞大小","细胞核大小","染色质","细胞核常规","有丝分裂","类别")
#print(data)
data$类别[data$类别==2]<-"良性"
data$类别[data$类别==4]<-"恶性"
#print(data)
data<-data[-1] #删除第一列元素#
#print(data)
set.seed(1234) #随机抽样设置种子
train<-sample(nrow(data),0.7*nrow(data)) #抽样函数，第一个参数为向量，nrow()返回行数 后面的是抽样参数前
tdata<-data[train,] #根据抽样参数列选择样本，都好逗号是选择行
vdata<-data[-train,] #删除抽样行


library(rpart)
dtree<-rpart(类别~.,data=tdata,method="class", parms=list(split="information"))
printcp(dtree)

library(rattle)
asRules(dtree)


opar<-par(no.readonly = T)
par(mfrow=c(1,2))
library(rpart.plot)
#png(file = "./R/tree1.png")
rpart.plot(dtree,branch=1,type=2, fallen.leaves=T,cex=0.8, sub="剪枝前")

par(opar)
dev.off()
