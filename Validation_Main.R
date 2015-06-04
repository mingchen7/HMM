setwd("C:\\Users\\mingchen7\\Google Drive\\4 Courses\\CE560\\Term project\\GA_Travel Time Decomposition\\");
load('LinkTT_10runs.RData');

source("GA_TravelTimeAllocation.R");
source("TravelTimeDecomposition.R");

# Model parameters
# Multivariate - K = 3
# mu2 = matrix(c(64.39,84.24,36.44,53.22,132.87,71.31,81.48,36.88,40.17,96.70,60.67,64.34,49.76,102.36,103.12),ncol=5,byrow=TRUE);
# sigma2 = array(0,dim=c(5,5,3));
# sigma2[,,1] = matrix(c(674.54,-327.20,2.168,-16.58,177.27,-327.20,484.59,4.933,47.057,-63.86,2.168,4.933,9.412,-4.391,16.904,-16.583,47.057,-4.391,80.320,-168.508,177.270,-63.866,16.904,-168.508,1569.619),ncol=5,byrow=TRUE);
# sigma2[,,2] = matrix(c(607.875,-171.285,-1.695,2.683,14.407,-171.285,416.381,1.495,20.225,38.613,-1.695,1.495,14.424,-4.294,-6.834,2.683,20.225,-4.294,18.468,4.507,14.407,38.613,-6.834,4.507,33.326),ncol=5,byrow=TRUE);
# sigma2[,,3] = matrix(c(556.202,-140.974,-4.953,6.083,-4.428,-140.974,371.512,44.138,22.477,-6.369,-4.953,44.138,179.011,-146.354,-15.017,6.0873,22.477,-146.354,200.895,40.886,-4.428,-6.369,-15.017,40.886,170.937),ncol=5,byrow=TRUE)
# Pi2 = c(0.077,0.415,0.506);

# Multivariate - K = 5
mu2 = matrix(c(65.97,75.896,39.951,67.415,169.657,72.181,79.912,35.830,39.328,97.040,36.707,68.785,50.249,100.477,100.357,60.403,86.392,38.954,46.279,94.444,72.074,62.511,49.582,102.560,101.421),ncol=5,byrow=TRUE);
sigma2 = array(0,dim=c(5,5,5));
sigma2[,,1] = matrix(c(614.559,-210.098,-29.355,-181.795,3.713,-210.098,502.008,-37.279,-288.513,23.904,-29.244,-37.279,93.794,100.081,10.500,-181.795,-288.513,100.081,1043.551,-21.571,3.713,23.904,10.500,-21.571,100.225),ncol=5,byrow=TRUE)
sigma2[,,2] = matrix(c(590.936,-138.961,-2.753,9.003,15.033,-138.961,409.852,-2.174,15.753,37.589,-2.753,-2.174,3.704,-1.335,-1.640,9.003,15.753,-1.335,10.261,3.842,15.033,37.589,-1.640,3.842,28.220),ncol=5,byrow=TRUE);
sigma2[,,3] = matrix(c(12.151,6.799,0.763,0.876,1.202,6.799,356.246,31.066,15.291,8.628,0.763,31.066,191.711,-147.743,-10.002,0.876,15.291,-147.743,177.549,2.668,1.202,8.628,-10.002,2.668,25.406),ncol=5,byrow=TRUE);
sigma2[,,4] = matrix(c(673.442,-269.773,17.786,-30.793,12.824,-269.773,427.453,-6.152,30.485,45.030,17.786,-6.152,29.649,-28.981,-10.831,-30.793,30.485,-28.981,77.559,-11.433,12.824,45.030,-10.831,-11.433,56.570),ncol=5,byrow=TRUE);
sigma2[,,5] = matrix(c(416.181,-140.444,-1.533,-9.996,-2.728,-140.444,371.397,49.734,36.651,3.475,-1.533,49.734,170.750,-142.195,-15.053,-9.996,36.651,-142.195,191.507,18.742,-2.728,3.475,-15.053,18.742,22.505),ncol=5,byrow=TRUE);
Pi2 = c(0.056,0.310,0.157,0.141,0.336);

# univariate
Pi = list()
Pi[[1]] = c(0.345,0.238,0.417);
Pi[[2]] = c(0.649,0.149,0.201);
Pi[[3]] = c(0.250,0.166,0.583);
Pi[[4]] = c(0.239,0.280,0.480);
Pi[[5]] = c(0.94395103,0.05604897);

mu = list()
mu[[1]] = c(89.19,39.76,58.95);
mu[[2]] = c(87.01,58.93,43.84);
mu[[3]] = c(58.78,39.60,35.60);
mu[[4]] = c(106.65,81.66,39.64);
mu[[5]] = c(98.76526,170.08265);

sigma = list()
sigma[[1]] = c(102.749,3.165,123.500);
sigma[[2]] = c(252.64,24.81,5.39);
sigma[[3]] = c(59.14,8.61,2.78);
sigma[[4]] = c(34.10,498.01,14.03);
sigma[[5]] = c(36.66009,72.18567);


row.num = sample(1:nrow(tt.WB),1000);
dataset = tt.WB[row.num,];
proportion = c(0.192,0.229,0.199,0.189,0.191);

x.benchmark = list();
MAPE.benchmark = data.frame();
MAE.benchmark = data.frame();

x.UniGMM = list();
MAE.UniGMM = data.frame();  
MAPE.UniGMM = data.frame();

x.MultiGMM = list();
MAE.MultiGMM = data.frame();  
MAPE.MultiGMM = data.frame();

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
  out.Uni = GA.MainFunction(mu,sigma,Pi,y.obs,'Univariate');
  x = out.Uni$solution;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('Univariate GMM:  ',x,'\n');
  
  x.UniGMM[[n]] = x;
  MAE.UniGMM = rbind(MAE.UniGMM,MAE);
  MAPE.UniGMM = rbind(MAPE.UniGMM,MAPE);  

  # Multivariate Gaussian Mixture Model
  out.Multi = GA.MainFunction(mu2,sigma2,Pi2,y.obs,'Multivariate');
  x = out.Multi$solution;
  MAE = abs(x - x.obs);
  MAPE = abs(x - x.obs) / x.obs;
  
  cat('Multivariate GMM:',x,'\n');
  
  x.MultiGMM[[n]] = x;
  MAE.MultiGMM = rbind(MAE.MultiGMM,MAE);
  MAPE.MultiGMM = rbind(MAPE.MultiGMM,MAPE); 
  
}

load('ModelComparison4.RData');
MAPE.byObs.benchmark = rowSums(MAPE.benchmark);
MAPE.byObs.UniGMM = rowSums(MAPE.UniGMM);
MAPE.byObs.MultiGMM = rowSums(MAPE.MultiGMM);

# Comparison by Observation
# Histogram
diff1 = MAPE.byObs.benchmark - MAPE.byObs.UniGMM;
diff2 = MAPE.byObs.benchmark - MAPE.byObs.MultiGMM;
diff3 = MAPE.byObs.UniGMM - MAPE.byObs.MultiGMM;

par(mfrow=c(3,1));
hist(diff1,breaks=50,main="Benchmark vs. Univariate GMM",xlab="MAPE(Benchmark) - MAPE(Univariate GMM)");
abline(v=0,col="red",lty=3,lwd=3);

hist(diff2,breaks=50,main="Benchmark vs. Multivariate GMM",xlab="MAPE(Benchmark) - MAPE(Multivariate GMM)");
abline(v=0,col="red",lty=3,lwd=3);

hist(diff3,breaks=50,main="Univariate GMM vs. Multivariate GMM",xlab="MAPE(Univariate GMM) - MAPE(Multivariate GMM)");
abline(v=0,col="red",lty=3,lwd=3);

#Error Plot
par(mfrow=c(1,1));
plot(x = MAPE.byObs.benchmark, y = MAPE.byObs.UniGMM,xlim=c(0,5),ylim=c(0,5));
lines(x=seq(0,5,1),y=seq(0,5,1));

plot(x = MAPE.byObs.benchmark, y = MAPE.byObs.MultiGMM,xlim=c(0,5),ylim=c(0,5));
lines(x=seq(0,5,1),y=seq(0,5,1));

# Comparison by variable
#Boxplot
par(mfrow=c(3,1));
boxplot(MAPE.benchmark$tt2,MAPE.UniGMM$tt2,MAPE.MultiGMM$tt2,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
boxplot(MAPE.benchmark$tt3,MAPE.UniGMM$tt3,MAPE.MultiGMM$tt3,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
boxplot(MAPE.benchmark$tt4,MAPE.UniGMM$tt4,MAPE.MultiGMM$tt4,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
boxplot(MAPE.benchmark$tt5,MAPE.UniGMM$tt5,MAPE.MultiGMM$tt5,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))
boxplot(MAPE.benchmark$tt6,MAPE.UniGMM$tt6,MAPE.MultiGMM$tt6,ylab="MAPE",xlab="Methods",names=c('Benchmark','Uni GMM','Multi GMM'))

MAPE.benchmark$method = 'Benchmark';
MAPE.UniGMM$method = 'Uni. GMM';
MAPE.MultiGMM$method = 'Multi. GMM';
MAPE.all = rbind(MAPE.benchmark,MAPE.UniGMM,MAPE.MultiGMM);

library(ggplot2)
library(gridExtra)
p1 = ggplot(MAPE.all,aes(x=method,y=tt3)) + theme_bw() +
  geom_boxplot()+ ylab("MAPE for tt3") + 
  theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));                                                                                                      

p2 = ggplot(MAPE.all,aes(x=method,y=tt4)) + theme_bw() +
  geom_boxplot() +  ylab("MAPE for tt4") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));

p3 = ggplot(MAPE.all,aes(x=method,y=tt5)) + theme_bw() +
  geom_boxplot() + ylab("MAPE for tt5") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));

p4 = ggplot(MAPE.all,aes(x=method,y=tt6)) + theme_bw() +
  geom_boxplot() + ylab("MAPE for tt6") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));
                                                                                                
grid.arrange(p1, p2, p3, p4, ncol=2)




