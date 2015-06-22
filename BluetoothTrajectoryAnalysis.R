rm(list=ls())
library("ggplot2")
library("gridExtra")
library("VIM")
library("sm")
library("mice")
library("mitools")

setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
source("MixtureModels_Univariate.R")
source("MixtureModels_Multivariate.R")
source("HiddenmarkovModel_Scaled.R")

####################################### 
## West Speedway section - Westbound ##
#######################################
data = read.csv('Bluetooth data\\Trajectory_WestSpeedway_WB_150512-150609.csv', header=T)

# subset none NA dat
trajectory = subset(data, !is.na(data$Camp2Mnt) & !is.na(data$Mnt2Park) & !is.na(data$Park2Euclid))
hist(trajectory$Camp2Mnt,breaks=100)
hist(trajectory$Mnt2Park,breaks=100)
hist(trajectory$Park2Euclid,breaks=100)

# Fit Univeraite GMM
GaussianMixture(trajectory$Camp2Mnt,3)
GaussianMixture(trajectory$Mnt2Park,3)
GaussianMixture(trajectory$Park2Euclid,3)

# Fit Multivariate GMM
GassuianMixture_Multi(trajectory[,c(3:5)],3)

# Fit HMM
pi.init = c(0.59006003,0.07263568,0.33730429);  # first link
mu.init = list()
mu.init[[1]] = c(95.81823,277.82411,56.34372);
mu.init[[2]] = c(63.17064,27.72879,180.58148);
mu.init[[3]] = c(222.2542,301.1632,160.0822);

sigma.init = list()
sigma.init[[1]] = c(798.1444,16666.0754,144.9083);
sigma.init[[2]] = c(540.51626,61.95351,7966.37868);
sigma.init[[3]] = c(303.0292,27072.0995,130.9301);

D = length(mu.init);
A.init = list();
for(i in 1:(D-1))
{
  K1 = length(mu.init[[i]]);
  K2 = length(mu.init[[i+1]]);
  A.init[[i]] = matrix(0,nrow=K1,ncol=K2,byrow=TRUE);
  
  #generating K1xK2 random numbers with summation of 1  
  for(j in 1:K1)
  {
    rnd = runif(K2);
    rnd = rnd / sum(rnd);
    A.init[[i]][j,] = rnd;
  }  
}

HiddenMarkovModel(trajectory[,c(3:5)], mu.init, sigma.init, pi.init, A.init);

######################################### 
## Middle Speedway section - Eastbound ##
#########################################
data = read.csv('Bluetooth data\\Trajectory_MidSpeedway_EB_140610-140725.csv', header=T)
colnames(data) = c('Date','TripID','Link1','Link2','Link3','Link4','Link5','Link6');
data$sumNA = apply(data,MARGIN=1,FUN=function(x) length(x[is.na(x)]));
aggr(data[,c(3:8)]);
data.lessNA = subset(data,sumNA<=3);
aggr(data.lessNA[,c(3:8)]);

nrows = nrow(data.lessNA);
ncols = ncol(data.lessNA);

# imputing missing data first
ini<-mice(data.lessNA,max=0,print=FALSE)
meth<-ini$meth
post<-ini$post
pred<-ini$pred

data.imp<-mice(data.lessNA,method="pmm",seed=30031,maxit=20)
stripplot(data.imp,pch=20,cex=1.2)
densityplot(data.imp,scales=list(x=list(relation="free")),layout=c(3,2))

data.all = complete(data.imp,"long")
trajectory = aggregate(data.all[,5:10],list(data.all$.id),mean);

hist(trajectory[,2],breaks=100,main="Alvernon to Columbus")
hist(trajectory[,3],breaks=100,main="Columbus to Swan")
hist(trajectory[,4],breaks=100,main="Swan to Rosemont")
hist(trajectory[,5],breaks=100,main="Rosemont to Craycroft")
hist(trajectory[,6],breaks=100,main="Craycroft to Wilmot")
hist(trajectory[,7],breaks=100,main="Wilmot to Kolb")

trajectory = subset(trajectory,trajectory[,2]<182 & trajectory[,3]<180 & trajectory[,4]<180 & trajectory[,5]<178 & trajectory[,6]<371 & trajectory[,7]<356);
summary(trajectory)

# fit univariate GMM
GaussianMixture(trajectory[,2],3) # Alvernon to Columbus
GaussianMixture(trajectory[,3],3) # Columbus to Swan
GaussianMixture(trajectory[,4],3) # Swan to Rosemont
GaussianMixture(trajectory[,5],3) # Rosemont to Craycroft
GaussianMixture(trajectory[,6],3) # Craycroft to Wilmot
GaussianMixture(trajectory[,7],3) # Wilmot to Kolb

# fit multivariate GMM
GassuianMixture_Multi(trajectory[,c(2:7)],3)

# Fit HMM
pi.init = c(0.6547958,0.2325503,0.1126539);
mu.init = list()
mu.init[[1]] = c(53.40069,66.85142,115.34070);
mu.init[[2]] = c(57.47536,74.91055,125.52219);
mu.init[[3]] = c(97.66305,51.14884,62.75426);
mu.init[[4]] = c(72.91124,160.31854,91.51193);
mu.init[[5]] = c(132.9261,188.7008,107.4272);
mu.init[[6]] = c(186.92614,153.90338,97.21446);


sigma.init = list()
sigma.init[[1]] = c(26.33393,135.17557,981.08500);
sigma.init[[2]] = c(67.8512,162.9658,631.5054);
sigma.init[[3]] = c(901.17347,31.00597,99.20501);
sigma.init[[4]] = c(177.6018,101.4874,420.1962);
sigma.init[[5]] = c(162.66843,3082.46684,65.88103);
sigma.init[[6]] = c(2530.5552,200.2076,242.4268);

D = length(mu.init);
A.init = list();
for(i in 1:(D-1))
{
  K1 = length(mu.init[[i]]);
  K2 = length(mu.init[[i+1]]);
  A.init[[i]] = matrix(0,nrow=K1,ncol=K2,byrow=TRUE);
  
  #generating K1xK2 random numbers with summation of 1  
  for(j in 1:K1)
  {
    rnd = runif(K2);
    rnd = rnd / sum(rnd);
    A.init[[i]][j,] = rnd;
  }  
}

HiddenMarkovModel(trajectory[,c(2:7)], mu.init, sigma.init, pi.init, A.init);

