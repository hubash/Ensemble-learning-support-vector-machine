rm(list=ls (all= TRUE) ) 

setwd("*****")

Ncv=10
nrep=5
#dataflag='SEER'
#if(dataflag=='SEER'){
#  data <- read.csv(file="SEER19732014.csv", header=TRUE, sep=",")}else{library("R.matlab")  
#    data <- readMat("SVMDATA1412.mat")}

#balance_flag='Out'
#  Breast Cancer Wisconsin (Diagnostic) 569 (wdbc) 1995

w<-data$WisDdata30
nfac=30
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
for(i in 1:nfac){w[,i]=range01(w[,i])}
###

#nx=22
#p1<-hist(w[w[,nfac+1] == 1, nx],prob = TRUE, col="green",border="green")#X15Num_nodes   X13Tumor_size  X12Age
#p1$counts=p1$counts/sum(p1$counts)
#lines(density(data[data[,nfac+1] == 1, 'X13Tumor_size']))
#p2<-hist(w[w[,nfac+1] == -1, nx],prob = TRUE, col="blue",border="blue")
#p2$counts=p2$counts/sum(p2$counts)

#plot( p1, col=rgb(0,0,1,1/4))
#plot( p2, col=rgb(1,0,0,1/4), add=T)

#nx=22
#plot(density(w[w[,nfac+1] == 1, nx]), col="blue", lwd=2)
#datpoint=density(w[w[,nfac+1] == 1, nx])
#x=datpoint$x
#y=datpoint$y
#xy=data.frame(x,y)
#write.csv(xy,'xy1.csv' )

#write.csv(w[,c(nx,nfac+1)],'22colandlabel.csv' )

#lines(density(w[w[,nfac+1] == -1, nx]), col="green", lwd=2)

# Breast Cancer Wisconsin (Original) 699 wobc 1992 WBC

# w<-na.omit(data$Wis9)
# nfac=9
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# for(i in 1:nfac){w[,i]=range01(w[,i])}


# Breast Cancer Wisconsin (Prognostic) 198 wpbc 1995 
#1-- nonrecur  -1--recur

#SEER data
nfac=14
data_major <- data[data[,nfac+1] == 1, ]
data_min <- data[data[,nfac+1] == -1, ]
indices <- sample(nrow(data_major), nrow(data_min))
temp_data <- rbind(data_major[indices,], data_min)
indices <- sample(nrow(temp_data), nrow(temp_data))
data1<-temp_data[indices,]
rm(data_major,data_min,data)
data1[,ncol(data1)+1]=0
dataNum <- matrix(data = NA, nrow = dim(data1)[1], ncol = dim(data1)[2])
for (i in 1:dim(data1)[2]) {
  dataNum[,i] <- c(as.numeric(data1[[i]]))
}
w<-dataNum

#################################
# p1<-hist(data[data[,nfac+1] == 1, 'X15Num_nodes'],prob = TRUE, col="green",border="green")#X15Num_nodes   X13Tumor_size  X12Age
# p1$counts=p1$counts/sum(p1$counts)
# #lines(density(data[data[,nfac+1] == 1, 'X13Tumor_size']))
# p2<-hist(data[data[,nfac+1] == -1, 'X15Num_nodes'],prob = TRUE, col="blue",border="blue")
# p2$counts=p2$counts/sum(p2$counts)
# 
# plot( p1, col=rgb(0,0,1,1/4))
# plot( p2, col=rgb(1,0,0,1/4), add=T)

source("0Ensemb_major.R")
source("0Ensemb_max.R")
source("0Ensemb_min.R")
source("0Ensemb_weightav.R")
source("0Ensemb_weightav_2.R")
source("0Ensemb_sbest.R")


plot(ROCmajor)
plot(ROCmax,add = TRUE)
plot(ROCmin,add = TRUE)
plot(ROCwtave,add = TRUE)
plot(ROCwtave2,add = TRUE)
plot(ROCsbest,add = TRUE)

colMeans(Repperfmajor)
colMeans(Repperfmax)
colMeans(Repperfmin)
colMeans(Repperfsbest)
colMeans(Repperfwtave)
colMeans(Repperfwtave2)



boxplot(Repperfmajor$Acc,Repperfmax$Acc,Repperfmin$Acc,Repperfwtave$Acc,Repperfwtave2$Acc,Repperfsbest$Acc)

#save.image(paste("wpbc","RData", sep="."))
