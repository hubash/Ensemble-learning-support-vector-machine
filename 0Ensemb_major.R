
#rm(list=ls (all= TRUE) ) 

#bootstrap
bss<-function(data,fn,p){
  bsslist <- list()
  ind=sample(1:round(10*p),nrow(data),replace = TRUE)
  temp=ind==1
  Trainsetx=data[!temp,1:(fn)]
  Trainsety=as.factor(data[!temp,(fn+1)])
  Testsetx=data[temp,1:(fn)]
  Testsety=as.factor(data[temp,(fn+1)])
  bsslist [[length(bsslist)+1]] <- Trainsetx 
  bsslist [[length(bsslist)+1]] <- Trainsety 
  bsslist [[length(bsslist)+1]] <- Testsetx 
  bsslist [[length(bsslist)+1]] <- Testsety 
  names(bsslist)<-c("Trx" ,"Try" ,"Tex" ,"Tey")
  return(bsslist)
}




library("kernlab")

# library("R.matlab")
# data <- readMat("SVMDATA1412.mat")
# nfac=30
# w<-data$WisDdata30
# #w=w[,-32]
# w[w[,31]==0,31]=-1

Repperf<-data.frame(matrix(0, nrow=nrep, ncol=6), row.names = NULL)
colnames(Repperf) = c("Acc", "Spec","Sens","Fscore","Auc","TrainTime")
for (rep in 1:nrep){
  #Cross validation
  w[,ncol(w)]=sample(1:Ncv,nrow(w),replace = TRUE)
  
  crossperfmajor=data.frame(matrix(0, nrow=Ncv, ncol=6), row.names = NULL)
  colnames(crossperfmajor) = c("Acc", "Spec","Sens","Fscore","Auc","TrainTime")
  
  ##----Start cross v--
  for (crv in 1:Ncv)
  {

      
    
temp=w[,ncol(w)]==crv
w1=w[!temp,]
if (balance_flag=='in'){
  print(balance_flag)
data_major <- w1[w1[,nfac+1] == 1, ]
data_min <- w1[w1[,nfac+1] == -1, ]
indices <- sample(nrow(data_major), nrow(data_min))
temp_data <- rbind(data_major[indices,], data_min)
indices <- sample(nrow(temp_data), nrow(temp_data))
w1<-temp_data[indices,]
}

w2=w[temp,]
ptm <- proc.time()  
#Bootstrap
bagginglist <- list()
bagginglistacc <- list()

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "rbfdot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "polydot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "vanilladot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "tanhdot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "laplacedot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "anovadot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 
#svp <- ksvm(Trainsetx, Trainsety, type = "C-svc", kernel = "splinedot",cross = 0,prob.model = TRUE)
#result=predict(svp,Testsetx)
#sum(result==Testsety)/nrow1(Testsety)
bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "rbfdot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "polydot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "vanilladot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "tanhdot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "laplacedot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

bres<-bss(w1,nfac,0.3)
svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "anovadot",cross = 0,prob.model = TRUE)
bagginglist [[length(bagginglist)+1]] <- svp 

ptm <- proc.time()-ptm   

#Predict
baggingPredi=matrix(0, nrow=nrow(w2), ncol=length(bagginglist)+1)
for (i in 1:length(bagginglist)){
  baggingPredi[,i]<- predict(bagginglist[[i]],w2[,1:nfac])
}
baggingPredi[baggingPredi==1]=-1
baggingPredi[baggingPredi==2]=1

#Majority
for (i in 1:nrow(baggingPredi)){
  res=as.numeric(names(which.max(table(baggingPredi[i,1:length(bagginglist)]))))
  baggingPredi[i,length(bagginglist)+1]=res
}


result=baggingPredi[,length(bagginglist)+1]

Acc=sum(result==w2[,nfac+1])/nrow(w2)
Spec=sum((result==-1)& (w2[,nfac+1]==-1))/sum((w2[,nfac+1]==-1))
Sens=sum((result==1)& (w2[,nfac+1]==1))/sum((w2[,nfac+1]==1))
Fscore=2*(Spec*Sens)/(Spec+Sens)

# ROC area under the curve
library(ROCR)
pred<-prediction(as.numeric(result),w2[,nfac+1])
auc.tmp <- performance(pred,"auc")
auc <- as.numeric(auc.tmp@y.values)


crossperfmajor$Acc[crv]=Acc
crossperfmajor$Spec[crv]=Spec
crossperfmajor$Sens[crv]=Sens
crossperfmajor$Fscore[crv]=Fscore
crossperfmajor$Auc[crv]=auc
crossperfmajor$TrainTime[crv]=ptm[3]
}
Repperf$Acc[rep]=mean(crossperfmajor$Acc)
Repperf$Spec[rep]=mean(crossperfmajor$Spec)
Repperf$Sens[rep]=mean(crossperfmajor$Sens)
Repperf$Fscore[rep]=mean(crossperfmajor$Fscore)
Repperf$Auc[rep]=mean(crossperfmajor$Auc)
Repperf$TrainTime[rep]=mean(crossperfmajor$TrainTime)
}
Repperfmajor<-Repperf

#ROC curve
pred<-prediction(as.numeric(result),w2[,nfac+1])
ROCmajor <- performance(pred,measure="tpr", x.measure="fpr")
#plot(perf)
