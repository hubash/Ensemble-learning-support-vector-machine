
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
  names(bsslist)<-c("Trx","Try","Tex","Tey")
  return(bsslist)
}

getperf<-function(svp,data){

  tempper <- data.frame(matrix(0, nrow=1, ncol=5), row.names = NULL)
  colnames(tempper) = c("Acc", "Spec","Sens","Fscore","Auc")
  result=predict(svp,data$Tex)
  
  tempper$Acc=sum(result==data$Tey)/nrow(data$Tex)
  tempper$Spec=sum((result==-1)& (data$Tey==-1))/sum((data$Tey==-1))
  tempper$Sens=sum((result==1)& (data$Tey==1))/sum((data$Tey==1))
  tempper$Fscore=2*(tempper$Spec*tempper$Sens)/(tempper$Spec+tempper$Sens)
  
  library(ROCR)
  pred<-prediction(as.numeric(result),data$Tey)
  auc.tmp <- performance(pred,"auc")
  tempper$Auc <- as.numeric(auc.tmp@y.values)
  return(tempper)
}




#Breast Cancer Data Set 286 
#  Breast Cancer Wisconsin (Diagnostic) 569 
#   Breast Cancer Wisconsin (Original) 699 
#  Breast Cancer Wisconsin (Prognostic) 198 

library("kernlab")

# library("R.matlab")
# setwd("C:/360/SameStep/ResearchD/Paper/ExpertSystempaper/ExpertSystemCode")
# data <- readMat("SVMDATA1412.mat")
# nfac=30
# w<-data$WisDdata30
# #w=w[,-32]
# w[w[,31]==0,31]=-1

#Ncv=5
#nrep=5
Repperf<-data.frame(matrix(0, nrow=nrep, ncol=6), row.names = NULL)
colnames(Repperf) = c("Acc", "Spec","Sens","Fscore","Auc","TrainTime")
for (rep in 1:nrep){
  #Cross validation
  w[,ncol(w)]=sample(1:Ncv,nrow(w),replace = TRUE)
  
  crossperfwtauc=data.frame(matrix(0, nrow=Ncv, ncol=6), row.names = NULL)
  colnames(crossperfwtauc) = c("Acc", "Spec","Sens","Fscore","Auc","TrainTime")
  
  ##----Start cross v--
  for (crv in 1:Ncv)
  {
    temp=w[,ncol(w)]==crv
    w1=w[!temp,]
    
    if (balance_flag=='in'){
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
    bagginglistacc <- data.frame(matrix(0, nrow=12, ncol=5), row.names = NULL)
    colnames(bagginglistacc) = c("Acc", "Spec","Sens","Fscore","Auc")
    
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "rbfdot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp
    bagginglistacc[1,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "polydot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
    bagginglistacc[2,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "vanilladot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
    bagginglistacc[3,]=getperf(svp,bres)
    
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "tanhdot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
    bagginglistacc[4,]=getperf(svp,bres)
    
    
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "laplacedot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[5,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "C-svc", kernel = "anovadot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[6,]=getperf(svp,bres)
    #svp <- ksvm(Trainsetx, Trainsety, type = "C-svc", kernel = "splinedot",cross = 0,prob.model = TRUE)
    #result=predict(svp,Testsetx)
    #sum(result==Testsety)/nrow1(Testsety)
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "rbfdot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[7,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "polydot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[8,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "vanilladot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[9,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "tanhdot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[10,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "laplacedot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[11,]=getperf(svp,bres)
    
    bres<-bss(w1,nfac,0.3)
    svp <- ksvm(bres$Trx, bres$Try, type = "nu-svc", kernel = "anovadot",cross = 0,prob.model = TRUE)
    bagginglist [[length(bagginglist)+1]] <- svp 
  bagginglistacc[12,]=getperf(svp,bres)
  ptm <- proc.time()-ptm 
    
    #Predict
    baggingPredi=matrix(0, nrow=nrow(w2), ncol=length(bagginglist)+1)
    for (i in 1:length(bagginglist)){
      baggingPredi[,i]<- predict(bagginglist[[i]],w2[,1:nfac])
    }
    baggingPredi[baggingPredi==1]=-1
    baggingPredi[baggingPredi==2]=1
    
    ##Weighted average (Acc*Auc)
    for (i in 1:nrow(baggingPredi)){
      #temp=bagginglistacc$Acc*bagginglistacc$Auc
      temp=bagginglistacc$Auc
      res=sum(baggingPredi[i,1:length(bagginglist)]*temp)/sum(temp)
      if(res<0){ 
        baggingPredi[i,length(bagginglist)+1]=-1
      }
      else
      {
        baggingPredi[i,length(bagginglist)+1]=1
      }
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
    
    crossperfwtauc$Acc[crv]=Acc
    crossperfwtauc$Spec[crv]=Spec
    crossperfwtauc$Sens[crv]=Sens
    crossperfwtauc$Fscore[crv]=Fscore
    crossperfwtauc$Auc[crv]=auc
    crossperfwtauc$TrainTime[crv]=ptm[3]
  }
  Repperf$Acc[rep]=mean(crossperfwtauc$Acc)
  Repperf$Spec[rep]=mean(crossperfwtauc$Spec)
  Repperf$Sens[rep]=mean(crossperfwtauc$Sens)
  Repperf$Fscore[rep]=mean(crossperfwtauc$Fscore)
  Repperf$Auc[rep]=mean(crossperfwtauc$Auc)
  Repperf$TrainTime[rep]=mean(crossperfwtauc$TrainTime)
}
Repperfwtave2<-Repperf

#ROC curve
pred<-prediction(as.numeric(result),w2[,nfac+1])
ROCwtave2<- performance(pred,measure="tpr", x.measure="fpr")
#Sensitivity/specificity plots: measure="sens", x.measure="spec".
#perf <- performance(pred,measure="sens", x.measure="spec")
#plot(perf)

#ROC curve
#pred<-prediction(as.numeric(result),w2[,nfac+1])
#perf <- performance(pred,measure="tpr", x.measure="fpr")
#plot(perf)
