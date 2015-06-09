setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
load('LinkTT_10runs.RData');

source("GA_TravelTimeAllocation.R");
source('HMM_FindMostprobableSequence.R');
source('Parameters.R');

row.num = sample(1:nrow(tt.WB),1000);
dataset = tt.WB[row.num,];
proportion = c(0.192,0.229,0.199,0.189,0.191); #used for Benchmark method

x.benchmark = list();
MAPE.benchmark = data.frame();
MAE.benchmark = data.frame();

x.UniGMM = list();
MAE.UniGMM = data.frame();  
MAPE.UniGMM = data.frame();

x.MultiGMM = list();
MAE.MultiGMM = data.frame();  
MAPE.MultiGMM = data.frame();

x.HMM = list();
MAE.HMM = data.frame();  
MAPE.HMM = data.frame();

for(n in 1:1000)
{  
  cat('n = ',n,'\n');
  
  y.obs = dataset[n,2];
  x.obs = dataset[n,c(3:7)];
  
  cat('Obs.:            ',x.obs$tt2,x.obs$tt3,x.obs$tt4,x.obs$tt5,x.obs$tt6,'\n');
  
  # Benchmark method
  x = proportion * y.obs;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('Benchmark:       ',x,'\n');
  
  x.benchmark[[n]] = x;
  MAE.benchmark = rbind(MAE.benchmark,MAE);
  MAPE.benchmark = rbind(MAPE.benchmark,MAPE);
  
  # Univariate Gaussian Mixture Model
  out.Uni = GA.MainFunction(mu,sigma,Pi,NULL,y.obs,'Univariate');
  x = out.Uni$solution;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('Univariate GMM:  ',x,'\n');
  
  x.UniGMM[[n]] = x;
  MAE.UniGMM = rbind(MAE.UniGMM,MAE);
  MAPE.UniGMM = rbind(MAPE.UniGMM,MAPE);  

  # Multivariate Gaussian Mixture Model
  out.Multi = GA.MainFunction(mu2,sigma2,Pi2,NULL,y.obs,'Multivariate');
  x = out.Multi$solution;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('Multivariate GMM:',x,'\n');
  
  x.MultiGMM[[n]] = x;
  MAE.MultiGMM = rbind(MAE.MultiGMM,MAE);
  MAPE.MultiGMM = rbind(MAPE.MultiGMM,MAPE);
  
  # Hidden Markov Model
  out.HMM = GA.MainFunction(mu3,sigma3,Pi3,A,y.obs,'HMM');
  x = out.HMM$solution;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('HMM             :',x,'\n');
  
  x.HMM[[n]] = x;
  MAE.HMM = rbind(MAE.HMM,MAE);
  MAPE.HMM = rbind(MAPE.HMM,MAPE);  
}

# load('ModelComparison4.RData');
# MAPE.byObs.benchmark = rowSums(MAPE.benchmark);
# MAPE.byObs.UniGMM = rowSums(MAPE.UniGMM);
# MAPE.byObs.MultiGMM = rowSums(MAPE.MultiGMM);
# 
# # Comparison by Observation
# # Histogram
# diff1 = MAPE.byObs.benchmark - MAPE.byObs.UniGMM;
# diff2 = MAPE.byObs.benchmark - MAPE.byObs.MultiGMM;
# diff3 = MAPE.byObs.UniGMM - MAPE.byObs.MultiGMM;
# 
# par(mfrow=c(3,1));
# hist(diff1,breaks=50,main="Benchmark vs. Univariate GMM",xlab="MAPE(Benchmark) - MAPE(Univariate GMM)");
# abline(v=0,col="red",lty=3,lwd=3);
# 
# hist(diff2,breaks=50,main="Benchmark vs. Multivariate GMM",xlab="MAPE(Benchmark) - MAPE(Multivariate GMM)");
# abline(v=0,col="red",lty=3,lwd=3);
# 
# hist(diff3,breaks=50,main="Univariate GMM vs. Multivariate GMM",xlab="MAPE(Univariate GMM) - MAPE(Multivariate GMM)");
# abline(v=0,col="red",lty=3,lwd=3);
# 
# #Error Plot
# par(mfrow=c(1,1));
# plot(x = MAPE.byObs.benchmark, y = MAPE.byObs.UniGMM,xlim=c(0,5),ylim=c(0,5));
# lines(x=seq(0,5,1),y=seq(0,5,1));
# 
# plot(x = MAPE.byObs.benchmark, y = MAPE.byObs.MultiGMM,xlim=c(0,5),ylim=c(0,5));
# lines(x=seq(0,5,1),y=seq(0,5,1));
# 
# # Comparison by variable
# #Boxplot
# par(mfrow=c(3,1));
# boxplot(MAPE.benchmark$tt2,MAPE.UniGMM$tt2,MAPE.MultiGMM$tt2,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
# boxplot(MAPE.benchmark$tt3,MAPE.UniGMM$tt3,MAPE.MultiGMM$tt3,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
# boxplot(MAPE.benchmark$tt4,MAPE.UniGMM$tt4,MAPE.MultiGMM$tt4,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
# boxplot(MAPE.benchmark$tt5,MAPE.UniGMM$tt5,MAPE.MultiGMM$tt5,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
# boxplot(MAPE.benchmark$tt6,MAPE.UniGMM$tt6,MAPE.MultiGMM$tt6,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
# 
# MAPE.benchmark$method = 'Benchmark';
# MAPE.UniGMM$method = 'Uni. GMM';
# MAPE.MultiGMM$method = 'Multi. GMM';
# MAPE.all = rbind(MAPE.benchmark,MAPE.UniGMM,MAPE.MultiGMM);
# 
# library(ggplot2)
# library(gridExtra)
# p1 = ggplot(MAPE.all,aes(x=method,y=tt3)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE for tt3") + 
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));                                                                                                      
# 
# p2 = ggplot(MAPE.all,aes(x=method,y=tt4)) + theme_bw() +
#   geom_boxplot() +  ylab("MAPE for tt4") +
#   theme(title=element_text(face="bold",size=16),        
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));
# 
# p3 = ggplot(MAPE.all,aes(x=method,y=tt5)) + theme_bw() +
#   geom_boxplot() + ylab("MAPE for tt5") +
#   theme(title=element_text(face="bold",size=16),        
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));
# 
# p4 = ggplot(MAPE.all,aes(x=method,y=tt6)) + theme_bw() +
#   geom_boxplot() + ylab("MAPE for tt6") +
#   theme(title=element_text(face="bold",size=16),        
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));
#                                                                                                 
# grid.arrange(p1, p2, p3, p4, ncol=2)




