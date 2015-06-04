library("mixdist")
library("mixtools")

rm(list=ls())
SimData<-read.csv(file.choose(), header=T)
WB.total<-SimData[which(SimData$No.==1),]
EB.total<-SimData[which(SimData$No.==7),]

#Eastbound
tt.WB = data.frame();
for(i in 1:nrow(WB.total))
{
  veh.no = WB.total[i,3] # get vehicle no
  trajectory = SimData[which(SimData$veh==veh.no),];
  if(nrow(trajectory) == 6)
  {
    trajectory = trajectory[order(trajectory$No.),];
    newrow = c(veh.no,trajectory[1,5],trajectory[2,5],trajectory[3,5],trajectory[4,5],trajectory[5,5],trajectory[6,5]);  
    tt.WB = rbind(tt.WB,newrow);
  }  
}
names(tt.WB) = c('veh','tt1','tt2','tt3','tt4','tt5','tt6');

#Westboudn
tt.EB = data.frame();
for(i in 1:nrow(EB.total))
{
  veh.no = EB.total[i,3] # get vehicle no
  trajectory = SimData[which(SimData$veh==veh.no),];
  if(nrow(trajectory) == 6)
  {
    trajectory = trajectory[order(trajectory$No.),];    
    newrow = c(veh.no,trajectory[1,5],trajectory[2,5],trajectory[3,5],trajectory[4,5],trajectory[5,5],trajectory[6,5]);  
    tt.EB = rbind(tt.EB,newrow);
  }    
}
names(tt.EB) = c('veh','tt7','tt8','tt9','tt10','tt11','tt12');

