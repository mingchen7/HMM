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





