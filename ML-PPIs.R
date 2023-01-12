###机器学习导引课程论文代码###
###李炳华 2022317110014###


### 清楚环境设置路径
rm(list = ls())
setwd("F:/李炳华专属文档/甲流实验/R-test")

#安装protr包并加载
install.packages("protr")
library(protr)

# 导入FASTA文件
#导入蛋白质序列文件
#导入正样本,使用protcheck()函数进行氨基酸类型健全性检查，并移除非标准的序列
pro <- readFASTA("0.3pro.fasta")
length(pro)
pro <- pro[(sapply(pro, protcheck))]
length(pro)
#根据正样本的比例导入负样本(1:1、1:3、1:5)，这里只做1:1
allpro <- readFASTA("all_new.fasta")
length(allpro)
allpro <- allpro[(sapply(allpro, protcheck))]
length(allpro)
restpro <- setdiff(allpro, pro)
nonpro <- sample(restpro, length(pro)*1, replace=FALSE, prob=NULL)
#nonpro <- sample(restpro, length(pro)*3, replace=FALSE, prob=NULL)
#nonpro <- sample(restpro, length(pro)*5, replace=FALSE, prob=NULL)

# 计算氨基酸组成、伪氨基酸组成和联合三元组法的描述符
x1 <- t(sapply(pro, extractCTriad))
# x1 <- t(sapply(pro, extractAAC))
# x1 <- t(sapply(pro, extractPAAC))
# x1 <- t(sapply(pro, extractCTriad))
x2 <- t(sapply(nonpro, extractCTriad))
# x2 <- t(sapply(nonpro, extractAAC))
# x2 <- t(sapply(nonpro, extractPAAC))
# x2 <- t(sapply(nonpro, extractCTriad))
x <- rbind(x1, x2)



# 设置分类标签
labels <- as.factor(c(rep(0, length(nonpro)), rep(1, length(pro))))

set.seed(1001)

# 划分75%的训练集和25%的测试集
tr.idx <- c(
  sample(1:nrow(x1), round(nrow(x1) * 0.75)),
  sample(nrow(x1) + 1:nrow(x2), round(nrow(x2) * 0.75))
)
te.idx <- setdiff(1:nrow(x), tr.idx)

x.tr <- x[tr.idx, ]
x.te <- x[te.idx, ]
y.tr <- labels[tr.idx]
y.te <- labels[te.idx]

### 随机森林分类模型
install.packages("randomForest")
library("randomForest")
rf.fit <- randomForest(x.tr, y.tr, cv.fold = 5)
rf.fit
# 在测试集中进行预测
rf.pred <- predict(rf.fit, x.te)
table(y.te,rf.pred)
# 输出评价参数
# install.packages("caret")
# library("caret")
# confusionMatrix(rf.pred,y.te,positive = '1')

### SVM模型
install.packages("e1071")
library("e1071")

# SVM_radial
svmr.fit <- svm(x.tr, y.tr, kernel='radial')
svmr.pred <- predict(svmr.fit, x.te)
table(y.te,svmr.pred)

# SVM_linear
svml.fit <- svm(x.tr, y.tr, kernel='linear')
svml.pred <- predict(svml.fit, x.te)
table(y.te,svml.pred)

# SVM_polynomia
svmp.fit <- svm(x.tr, y.tr, kernel='polynomia')
svmp.pred <- predict(svmp.fit, x.te)
table(y.te,svmp.pred)


### KNN模型
install.packages('class')
library("class")
knn.pred <- knn(x.tr, x.te, y.tr)
table(y.te,knn.pred)


### 朴素贝叶斯模型
install.packages("e1071")
library("e1071")
nb.fit <- naiveBayes(x.tr, y.tr)
nb.pred <- predict(nb.fit, x.te)
table(y.te,nb.pred)

