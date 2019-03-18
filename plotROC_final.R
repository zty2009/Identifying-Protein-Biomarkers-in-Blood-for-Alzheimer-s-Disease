data=read.table('C:/Users/zty20/Documents/R/Protein_AD/cross_ad_protein_new.txt',stringsAsFactors=FALSE)

feature=read.table('C:/Users/zty20/Documents/R/Protein_AD/feature42.txt',stringsAsFactors=FALSE)
AD=feature[39,1:1481]
AD=t(AD)
A=rep(AD,184)
data1=as.matrix(data)
data1=as.vector(t(data1))
library(pROC)
#############cross validation roc############
	r=plot.roc(A,data1)
	roc1<-roc(A,data1)
	auc1=roc1$auc

#############cross validation roc for InBioMap############
data=read.table('C:/Users/zty20/Documents/R/Protein_AD/inbiomap_crossvalidation.txt',stringsAsFactors=FALSE)
network_data=data[,3:1483]
data1=as.matrix(network_data)
data1=as.vector(t(data1))
	r1=plot.roc(A,data1)
	roc2<-roc(A,data1)
	auc2=roc2$auc
#############roc for initialvalue 0.6############
data=read.table('C:/Users/zty20/Documents/R/Protein_AD/result_0.6.txt',stringsAsFactors=FALSE)
r=plot.roc(as.factor(AD),as.numeric(data[1,]))
	roc3<-roc(as.factor(AD),as.numeric(data[1,]))
	auc3=roc1$auc
	
	#############roc for InBioMap############
data=read.table('C:/Users/zty20/Documents/R/Protein_AD/inbiomap_randomwalk_result.txt',stringsAsFactors=FALSE)
r3=plot.roc(as.factor(AD),as.numeric(data[1,]))
	roc4<-roc(as.factor(AD),as.numeric(data[1,]))
	auc4=roc1$auc