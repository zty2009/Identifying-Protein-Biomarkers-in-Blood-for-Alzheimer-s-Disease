rm(list = ls())
library(lars)
library(glmnet)
library(ridge)
data=read.table('C:/Users/zty20/Documents/R/Protein_AD/feature42.txt',stringsAsFactors=FALSE)
x1=as.matrix(data[,1:1481])
y=as.numeric(unlist(data[,1482]))
adpro=which(x1[nrow(x1),]==1)

for (lpl in 1: length(adpro)) {
x=x1
x[nrow(x1),adpro[lpl]]=0
a=data.frame(x,y)
ad=which(x[nrow(x),]!=0)
ppp=0.6
w=matrix((1-ppp)/(1481-length(ad)),1481,1)
w[ad]=ppp/length(ad)
new_w=w
T = 1 
Tmin = 1e-8 
k = 100  
delta = 0.9  
result=matrix(0,k,1)

calculate <- function(initial_w,x)
{
sum_fw=0
for (i in 1:nrow(x)){
se=(t(initial_w)%*%(1-x[i,])*(1-t(initial_w)%*%(1-x[i,]))/sum(x[i,]))^0.5
fw=abs((t(initial_w)%*%x[i,]-y[i])/se)
sum_fw=sum_fw+fw
} 
return(sum_fw)
}


for (j in 1:k){
if (j%%4==1){
result[j]=calculate(w,x)	
X=matrix(0,nrow(x),ncol(x))
	for (p in 1:39){
		X[p,]=w*x[p,]
	}
	
laa = lars(X, y, type = "lar",use.Gram=FALSE)    
r=summary(laa)
best=which(r$Cp==min(abs(r$Cp)))-1
index=which(laa$beta[best,]!=0)


lp=1:1481
ziji=setdiff(lp,index)
loss=min(w[ziji])/2
new_w[index]=w[index]+loss
sumloss=loss*length(index)
new_w[ziji]=new_w[ziji]-sumloss/length(ziji)

result[j+1]=calculate(new_w,x)
if (result[j]<result[j+1]){
	prop=1/(1+exp(-sum(abs((w-new_w)))/T))
	random=runif(1,min=0,max=1)
	if (random>prop){
		w=new_w
	}
}
if(result[j]>result[j+1]){
w=new_w
}

}


else{
result[j]=calculate(w,x)	
X=matrix(0,nrow(x),ncol(x))
	for (p in 1:39){
		X[p,]=w*x[p,]
	}

r=glmnet(X,y,alpha=0.5)
best=which(r$dev.ratio==max(abs(r$dev.ratio)))
index=which(r$beta[,best]!=0)

lp=1:1481
ziji=setdiff(lp,index)
loss=min(w[ziji])/2
new_w[index]=w[index]+loss
sumloss=loss*length(index)
new_w[ziji]=new_w[ziji]-sumloss/length(ziji)

result[j+1]=calculate(new_w,x)
if (result[j]<result[j+1]){
	prop=1/(1+exp(-sum(abs((w-new_w)))/T))
	random=runif(1,min=0,max=1)
	if (random>prop){
		w=new_w
	}
}
if(result[j]>result[j+1]){
w=new_w
}

}
}
write.table(file='cross_ad_protein_new.txt',t(w),append=T,col.names=F,row.names=F)


}



mymod=linearRidge(y~.,a)
r=summary(mymod)
pvalue=r[[6]]$summary16$coefficients[,5]
index=which(pvalue<0.05)-1







