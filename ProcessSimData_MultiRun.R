library("mixdist")
library("mixtools")

rm(list=ls())
SimData<-read.csv(file.choose(), header=T)

tt.WB = data.frame();
tt.EB = data.frame();

for(run in 1:10)
{
  WB.total<-SimData[which((SimData$No.==1) & (SimData$Run == run)),]
  EB.total<-SimData[which((SimData$No.==7) & (SimData$Run == run)),]
  
  #Eastbound
  for(i in 1:nrow(WB.total))
  {
    veh.no = WB.total[i,4] # get vehicle no
    trajectory = SimData[which((SimData$veh==veh.no) & (SimData$Run==run)),];
    trajectory = trajectory[order(trajectory$No.),];
    newrow = c(veh.no,trajectory[1,6],trajectory[2,6],trajectory[3,6],trajectory[4,6],trajectory[5,6],trajectory[6,6]);
    tt.WB = rbind(tt.WB,newrow);
  }

  #Westboudn
  for(i in 1:nrow(EB.total))
  {
    veh.no = EB.total[i,4] # get vehicle no
    trajectory = SimData[which((SimData$veh==veh.no) & (SimData$Run==run)),];
    trajectory = trajectory[order(trajectory$No.),];
    newrow = c(veh.no,trajectory[1,6],trajectory[2,6],trajectory[3,6],trajectory[4,6],trajectory[5,6],trajectory[6,6]);
    tt.EB = rbind(tt.EB,newrow);
  }
}
names(tt.WB) = c('veh','tt1','tt2','tt3','tt4','tt5','tt6');
names(tt.EB) = c('veh','tt7','tt8','tt9','tt10','tt11','tt12');


load("LinkTravelTime.RData")
#WB Total
hist(tt.WB$tt1,breaks=100)
mix.fit1<-normalmixEM(tt.WB$tt1,lambda=0.333,mu=c(150,250,350),sigma=c(50,50,50),maxit=3000)
plot(mix.fit1,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit1)

#WB ->Campbell
hist(tt.WB$tt2,breaks=100)
mix.fit2<-normalmixEM(tt.WB$tt2,lambda=0.333,mu=c(20,50,80),sigma=c(10,10,10),maxit=3000)
plot(mix.fit2,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit2)

#WB Campbell->Cherry
hist(tt.WB$tt3,breaks=100)
mix.fit3<-normalmixEM(tt.WB$tt3,lambda=0.333,mu=c(20,40,80),sigma=c(10,10,10),maxit=3000)
plot(mix.fit3,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit3)

#WB Cherry->Mountain
hist(tt.WB$tt4,breaks=100)
mix.fit4<-normalmixEM(tt.WB$tt4,lambda=0.333,mu=c(20,40,80),sigma=c(10,10,10),maxit=3000)
plot(mix.fit4,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit4)

#WB Mountain->Park
hist(tt.WB$tt5,breaks=100)
mix.fit5<-normalmixEM(tt.WB$tt5,lambda=0.333,mu=c(30,36,40),sigma=c(5,5,5),maxit=3000)
plot(mix.fit5,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit5)

#WB Park->Euclid
hist(tt.WB$tt6,breaks=100)
mix.fit6<-normalmixEM(tt.WB$tt6,lambda=0.333,mu=c(30,36,40),sigma=c(10,10,10),maxit=3000)
plot(mix.fit6,density=TRUE,breaks=100,cex.axis=1.2,cex.lab=1.2,cex.main=1.4,main2="Mixure Gaussian model",xlab2="Travel Time")
summary(mix.fit6)

lambda=c(1/3,1/3,1/3)
mix.out = mvnormalmixEM(tt.WB[,c(2:7)], arbmean = FALSE, lambda=lambda, epsilon = 1e-02)




#EB Total
hist(tt.EB$tt7,breaks=100)
hist(tt.EB$tt8,breaks=100)
hist(tt.EB$tt9,breaks=100)
hist(tt.EB$tt10,breaks=100)
hist(tt.EB$tt11,breaks=100)
