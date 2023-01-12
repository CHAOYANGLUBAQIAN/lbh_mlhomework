###����ѧϰ�����γ����Ĵ���###
###����� 2022317110014###


### �����������·��
rm(list = ls())
setwd("F:/�����ר���ĵ�/����ʵ��/R-test")

#��װprotr��������
install.packages("protr")
library(protr)

# ����FASTA�ļ�
#���뵰���������ļ�
#����������,ʹ��protcheck()�������а��������ͽ�ȫ�Լ�飬���Ƴ��Ǳ�׼������
pro <- readFASTA("0.3pro.fasta")
length(pro)
pro <- pro[(sapply(pro, protcheck))]
length(pro)
#�����������ı������븺����(1:1��1:3��1:5)������ֻ��1:1
allpro <- readFASTA("all_new.fasta")
length(allpro)
allpro <- allpro[(sapply(allpro, protcheck))]
length(allpro)
restpro <- setdiff(allpro, pro)
nonpro <- sample(restpro, length(pro)*1, replace=FALSE, prob=NULL)
#nonpro <- sample(restpro, length(pro)*3, replace=FALSE, prob=NULL)
#nonpro <- sample(restpro, length(pro)*5, replace=FALSE, prob=NULL)

# ���㰱������ɡ�α��������ɺ�������Ԫ�鷨��������
x1 <- t(sapply(pro, extractCTriad))
# x1 <- t(sapply(pro, extractAAC))
# x1 <- t(sapply(pro, extractPAAC))
# x1 <- t(sapply(pro, extractCTriad))
x2 <- t(sapply(nonpro, extractCTriad))
# x2 <- t(sapply(nonpro, extractAAC))
# x2 <- t(sapply(nonpro, extractPAAC))
# x2 <- t(sapply(nonpro, extractCTriad))
x <- rbind(x1, x2)



# ���÷����ǩ
labels <- as.factor(c(rep(0, length(nonpro)), rep(1, length(pro))))

set.seed(1001)

# ����75%��ѵ������25%�Ĳ��Լ�
tr.idx <- c(
  sample(1:nrow(x1), round(nrow(x1) * 0.75)),
  sample(nrow(x1) + 1:nrow(x2), round(nrow(x2) * 0.75))
)
te.idx <- setdiff(1:nrow(x), tr.idx)

x.tr <- x[tr.idx, ]
x.te <- x[te.idx, ]
y.tr <- labels[tr.idx]
y.te <- labels[te.idx]

### ���ɭ�ַ���ģ��
install.packages("randomForest")
library("randomForest")
rf.fit <- randomForest(x.tr, y.tr, cv.fold = 5)
rf.fit
# �ڲ��Լ��н���Ԥ��
rf.pred <- predict(rf.fit, x.te)
table(y.te,rf.pred)
# ������۲���
# install.packages("caret")
# library("caret")
# confusionMatrix(rf.pred,y.te,positive = '1')

### SVMģ��
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


### KNNģ��
install.packages('class')
library("class")
knn.pred <- knn(x.tr, x.te, y.tr)
table(y.te,knn.pred)


### ���ر�Ҷ˹ģ��
install.packages("e1071")
library("e1071")
nb.fit <- naiveBayes(x.tr, y.tr)
nb.pred <- predict(nb.fit, x.te)
table(y.te,nb.pred)
