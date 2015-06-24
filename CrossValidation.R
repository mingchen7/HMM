library("caret")
setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");

source("GA_TravelTimeAllocation.R");
source('HMM_FindMostProbableSequence.R');
source("MixtureModels_Univariate.R");
source("MixtureModels_Multivariate.R");
source("HiddenmarkovModel_Scaled.R");

# sink('Log.txt');

CrossValidation = function(data,nfolds,proportion)
{
  K.uniGMM = 3;
  K.multiGMM = 3;
  # Split the data into ten bins
  nrows = nrow(data);
  ncols = ncol(data);
  folds = createFolds(1:nrows,k=nfolds,list=TRUE,returnTrain=FALSE);
  order = sample(1:10,10);
  
  output = list();
  
  # loop for cross-validation
  for(run in 1:10)
  {
    cat('Current Run: ',run,'\n');
    idx.validation = folds[[order[run]]];
    # specify training dataset
    data.training = data[-idx.validation,];
    data.validation = data[idx.validation,];
    
    ## TRAINING PART ##
    # train Univariate GMM
    model.uniGMM = list();
    for(i in 1:ncols)
    {
      cat('Fitting Univaraite GMM #',i,'\n');
      model.uniGMM[[i]] = GaussianMixture(data.training[,i],K.uniGMM);
      par(mfrow=c(1,1)); #reset plot
    }
    
    # training Multivariate GMM
    cat('Fitting Multivariate GMM\n');
    model.multiGMM = GassuianMixture_Multi(data.training,K.multiGMM);
    
    # training HMM
    pi.init = model.uniGMM[[1]]$Pi;
    pi.uni=list();mu.uni = list();sigma.uni=list();
    for(i in 1:ncols)
    {
      pi.uni[[i]] = model.uniGMM[[i]]$Pi;
      mu.uni[[i]] = model.uniGMM[[i]]$Mu;
      sigma.uni[[i]] = model.uniGMM[[i]]$Sigma;
    }
            
    D = length(mu.uni);
    A.init = list();
    for(i in 1:(D-1))
    {
      K1 = length(mu.uni[[i]]);
      K2 = length(mu.uni[[i+1]]);
      A.init[[i]] = matrix(0,nrow=K1,ncol=K2,byrow=TRUE);          
      for(j in 1:K1)
      {
        rnd = runif(K2);
        rnd = rnd / sum(rnd);
        A.init[[i]][j,] = rnd;
      }  
    }
    cat('Fitting Hidden Markov model\n');
    model.HMM = HiddenMarkovModel(data.training,mu.uni,sigma.uni,pi.init,A.init);
    
    ## VALIDATION PART ##
    x.benchmark = data.frame(); #store estimated x
    MAPE.benchmark = data.frame(); #store MAPE
    MAE.benchmark = data.frame(); #store MAE
    
    x.UniGMM = data.frame();
    MAE.UniGMM = data.frame();  
    MAPE.UniGMM = data.frame();
    
    x.MultiGMM = data.frame();
    MAE.MultiGMM = data.frame();  
    MAPE.MultiGMM = data.frame();
    
    x.HMM = data.frame();
    MAE.HMM = data.frame();  
    MAPE.HMM = data.frame();
    
    nrows.validation = nrow(data.validation);
    for(n in 1:nrows.validation)
    {      
      x.obs = data.validation[n,];
      y.obs = sum(x.obs);      
      cat('n = ',n,'\n');
      print(x.obs);      
      
      # Benchmark method
      x = proportion * y.obs;
      MAE = abs(x - x.obs);
      MAPE = abs(x - x.obs) / x.obs;      
      cat('Benchmark:       ',x,'\n');      
      x.benchmark = rbind(x.benchmark,x);
      MAE.benchmark = rbind(MAE.benchmark,MAE);
      MAPE.benchmark = rbind(MAPE.benchmark,MAPE);
      
      # Univariate Gaussian Mixture Model      
      out.Uni = GA.MainFunction(mu.uni,sigma.uni,pi.uni,NULL,y.obs,'Univariate');
      x = out.Uni$solution;
      MAE = abs(x - x.obs);
      MAPE = abs(x - x.obs) / x.obs;      
      cat('Univariate GMM:  ',x,'\n');
      x.UniGMM = rbind(x.UniGMM,x);      
      MAE.UniGMM = rbind(MAE.UniGMM,MAE);
      MAPE.UniGMM = rbind(MAPE.UniGMM,MAPE);  
      
      # Multivariate Gaussian Mixture Model
      out.Multi = GA.MainFunction(model.multiGMM$Mu,model.multiGMM$Sigma,model.multiGMM$Pi,NULL,y.obs,'Multivariate');
      x = out.Multi$solution;
      MAE = abs(x - x.obs);
      MAPE = abs(x - x.obs) / x.obs;      
      cat('Multivariate GMM:',x,'\n');
      x.MultiGMM = rbind(x.MultiGMM,x);      
      MAE.MultiGMM = rbind(MAE.MultiGMM,MAE);
      MAPE.MultiGMM = rbind(MAPE.MultiGMM,MAPE);
      
      # Hidden Markov model
      out.HMM = GA.MainFunction(model.HMM$Mu,model.HMM$Sigma,model.HMM$Pi,model.HMM$A.transition,y.obs,'HMM');
      x = out.HMM$solution;
      MAE = abs(x - x.obs);
      MAPE = abs(x - x.obs) / x.obs;      
      cat('HMM:             ',x,'\n');
      x.HMM = rbind(x.HMM,x);      
      MAE.HMM = rbind(MAE.HMM,MAE);
      MAPE.HMM = rbind(MAPE.HMM,MAPE);       
    }
    
    #store results
    output[[run]] = list()
    output[[run]]$x.benchmark = x.benchmark;
    output[[run]]$MAE.benchmark = MAE.benchmark;
    output[[run]]$MAPE.benchmark = MAPE.benchmark;
    output[[run]]$x.UniGMM = x.UniGMM;
    output[[run]]$MAE.UniGMM = MAE.UniGMM;
    output[[run]]$MAPE.UniGMM = MAPE.UniGMM;
    output[[run]]$x.MultiGMM = x.MultiGMM;
    output[[run]]$MAE.MultiGMM = MAE.MultiGMM;
    output[[run]]$MAPE.MultiGMM = MAPE.MultiGMM;
    output[[run]]$x.HMM = x.HMM;
    output[[run]]$MAE.HMM = MAE.HMM;
    output[[run]]$MAPE.HMM = MAPE.HMM;
  }  
  return(output);  
}

# load data
load('TrajectoryData_MidSpeedway_EB.RData');
#load('TrajectoryData_WestSpeedway_WB.RData');

#proportion.SimuData = c(0.192,0.229,0.199,0.189,0.191);
#proportion.WestSpeedway_WB = c(0.5405105,0.2807991,0.1786903);
proportion.MidSpeedway_EB = c(0.125725,0.124182,0.124481,0.123237,0.256225,0.24615);
results = CrossValidation(trajectory[,2:7],10,proportion.MidSpeedway_EB);
save(results,file='Result_CrossValidation_MidSpeedway_EB.RData');
#save(results,file='Result_CrossValidation_WestSpeedway_WB.RData');

