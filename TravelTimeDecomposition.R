# Travel Time Decomposition: main function
# Use EM for updating mixing responsibility
# use GA to maximize log-likelihood

TravelTimeDecomposition = function(mu,sigma,Pi,y.obs)
{
  #initiailization
  D = length(mu);

  #Gamma = Pi;
  Gamma = list();
  for(i in 1:D)
  {
    K = length(mu[[i]]);
    Gamma[[i]] = c(0);
    for(k in 1:K)
    {
      Gamma[[i]][k] = 1/K; 
    }    
  }
  
  # solution for observation
  # use Pi to calculate the initial values for x  
  GA.result = GA.MainFunction(mu,sigma,Pi,y.obs,'Univariate');
  x = c();
  x = GA.result$solution;
  LogLikelihood = c(0);
  solutions = list();
  
  LogLikelihood[1] = GA.result$loglikelihood;
  solutions[[1]] = x;
    
  for(iter in 1:20)
  {
    # E-step
    for(i in 1:D)
    {
      K = length(mu[[i]]);
      sum = 0;
      for(k in 1:K)
      {
        sum = sum + Pi[[i]][k] * dnorm(x[i], mu[[i]][k], sd = sqrt(sigma[[i]][k]));
      }
      
      for(k in 1:K)
      {
        Gamma[[i]][k] = (Pi[[i]][k] * dnorm(x[i], mu[[i]][k], sd = sqrt(sigma[[i]][k]))) / sum; 
      }
      
    }
    
    # control of values
    for(i in 1:D)
    {
      Gamma[[i]][Gamma[[i]] < 1e-60] = 1e-15;
      Gamma[[i]][Gamma[[i]] > (1-1e-60)] = 1e-15;
    }
        
    # M-step
    GA.result = GA.MainFunction(mu,sigma,Gamma,y.obs,'Univariate');
    x = GA.result$solution;
    solutions[[iter+1]] = x;
    
    # Evaluate log-likelihood
    LogLikelihood[iter+1] = GA.result$loglikelihood;
            
    print(iter);
#     print(Gamma);
  }

  # plot log-likelihood
  plot(x = 1:length(LogLikelihood), y = LogLikelihood, xlab='iteration',type = 'l');
  
  output = list();
  output$solutions = solutions;
  output$Loglikelihood = LogLikelihood;
  output$Gamma = Gamma;
  return(output);
}

# setwd("C:\\Users\\mingchen7\\Google Drive\\4 Courses\\CE560\\Term project\\GA_Travel Time Decomposition\\");
#load('LinkTT_10runs.RData');
# source("GA_TravelTimeAllocation.R");

## TEST DATA 
# multivariate
# mu = matrix(c(64.39,84.24,36.44,53.22,132.87,71.31,81.48,36.88,40.17,96.70,60.67,64.34,49.76,102.36,103.12),ncol=5,byrow=TRUE);
# sigma = array(0,dim=c(5,5,3))
# sigma[,,1] = matrix(c(674.54,-327.20,2.168,-16.58,177.27,-327.20,484.59,4.933,47.057,-63.86,2.168,4.933,9.412,-4.391,16.904,-16.583,47.057,-4.391,80.320,-168.508,177.270,-63.866,16.904,-168.508,1569.619),ncol=5,byrow=TRUE);
# sigma[,,2] = matrix(c(607.875,-171.285,-1.695,2.683,14.407,-171.285,416.381,1.495,20.225,38.613,-1.695,1.495,14.424,-4.294,-6.834,2.683,20.225,-4.294,18.468,4.507,14.407,38.613,-6.834,4.507,33.326),ncol=5,byrow=TRUE);
# sigma[,,3] = matrix(c(556.202,-140.974,-4.953,6.083,-4.428,-140.974,371.512,44.138,22.477,-6.369,-4.953,44.138,179.011,-146.354,-15.017,6.0873,22.477,-146.354,200.895,40.886,-4.428,-6.369,-15.017,40.886,170.937),ncol=5,byrow=TRUE)
# Gamma = c(0.077,0.415,0.506);

# univariate
# Pi = list()
# Pi[[1]] = c(0.345,0.238,0.417);
# Pi[[2]] = c(0.649,0.149,0.201);
# Pi[[3]] = c(0.250,0.166,0.583);
# Pi[[4]] = c(0.239,0.280,0.480);
# Pi[[5]] = c(0.94395103,0.05604897);
# 
# mu = list()
# mu[[1]] = c(89.19,39.76,58.95);
# mu[[2]] = c(87.01,58.93,43.84);
# mu[[3]] = c(58.78,39.60,35.60);
# mu[[4]] = c(106.65,81.66,39.64);
# mu[[5]] = c(98.76526,170.08265);
# 
# sigma = list()
# sigma[[1]] = c(102.749,3.165,123.500);
# sigma[[2]] = c(252.64,24.81,5.39);
# sigma[[3]] = c(59.14,8.61,2.78);
# sigma[[4]] = c(34.10,498.01,14.03);
# sigma[[5]] = c(36.66009,72.18567);
# 
# out = TravelTimeDecomposition(mu,sigma,Pi,401.3);