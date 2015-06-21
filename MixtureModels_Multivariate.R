library(ggplot2)
library(lubridate)
library(mixtools)
library(MASS)

# mu: 1xM matrix
# covs: MxM matrix
# x: 1xM vector

ProbGaussian = function(mu, covs, X)
{
  sum = 0;
  D = length(mu);
  X = as.matrix(X);
  const = -(D/2)*log(2*pi) - 0.5*log(det(covs));
  Xm = X - mu;
  Xm = matrix(Xm,nrow=1,ncol=D)
  lg.px = const - 0.5 * diag(Xm %*% solve(covs) %*% t(Xm));    
  
  return(exp(lg.px));
}


GassuianMixture_Multi = function(data, numCluster)
{
    #get dimension of variables
    data = as.matrix(data);
    D = ncol(data);
    N = nrow(data);
    # using Kmeans to cluster the data first              
    KmeansCluster = kmeans(data, numCluster, iter.max = 20, algorithm = 'Hartigan-Wong');  
    
    # model parameters: needs to be estimated eventually
    Pi = c();  # 1xK vector
    Mu = matrix(nrow=numCluster,ncol=D);  #KxD matrix
    Sigma = array(0,dim=c(D,D,numCluster)); #DxDxK array
    
    # Initialize pi_k, mu_k, sigma_k
    for(k in 1: numCluster)
    {
        Pi[k] = KmeansCluster$size[k] / N;
        for(var in 1:D)
        {
          Mu[k,var] = mean(data[KmeansCluster$cluster == k,var]);  
        }        
        Sigma[,,k] = cov(data[KmeansCluster$cluster == k,]);
    }
    
    Q = matrix(data = NA, nrow = N, ncol = numCluster); #NxK matrix, responsibility for each observation
    LogLikelihood = c(0);
    
    for(iter in 1: 50)
    {
        # E-Step
        for(n in 1: N)
        {
            SumDenominator = 0;
            x = data[n,];
                        
            for(k in 1: numCluster)
            {
                covs = Sigma[,,k];
                mu = Mu[k,];
                Temp = Pi[k] * ProbGaussian(mu,covs,x);
                SumDenominator = SumDenominator + Temp;
            }
            
            # update Q
            for(k in 1: numCluster)
            {
                covs = Sigma[,,k];
                mu = Mu[k,];
                Q[n, k] = Pi[k] * ProbGaussian(mu,covs,x) / SumDenominator;
            }
        }
        
        # Minor control for numeric issues
        Q[Q < 1e-60] = 1e-15;
        Q[Q > 1-1e-60] = 1e-15;
            
        # M-Step
        # update Pi,Mu, and Sigma
        for(k in 1: numCluster)
        {
            Pi[k] = mean(Q[, k]);            
            Mu[k,] = colSums(matrix(rep(Q[, k],D), nrow=N) * data) / sum(Q[, k]);
            
            Xm = data - matrix(rep(Mu[k,],N), ncol=D, byrow=TRUE);
            Sigma[,,k] = t(Xm * matrix(rep(Q[, k],D), nrow=N)) %*% Xm;
            Sigma[,,k] = Sigma[,,k] / sum(Q[,k]);                        
        }
        
        # Computing log-likelihood
        LogLikelihood[iter] = 0;
        for(n in 1: N)
        {
            likelihood = 0;
            x = data[n,]
            for(k in 1: numCluster)
            {              
                mu = Mu[k,];
                covs = Sigma[,,k];
                likelihood = likelihood + Pi[k] * ProbGaussian(mu,covs,x);
            } 
            
            LL = log(likelihood); 
            LogLikelihood[iter] = LogLikelihood[iter] + LL;
            
        }            
        print(sprintf('iter = %d  loglikelihood = %f', iter, LogLikelihood[iter]));
#         print(Pi);
#         print(Mu);
#         print(Sigma);
    }
            
    plot(x = 1:length(LogLikelihood), y = LogLikelihood, xlab='iteration',type = 'l');
    
    Estimation  = list();
    Estimation$Pi = Pi;
    Estimation$Mu = Mu;
    Estimation$Sigma = Sigma;
    Estimation$loglikelihood = LogLikelihood[length(LogLikelihood)];
    
    return(Estimation);
}

# example: 
# setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
# load('LinkTT_10hours.RData'); 
# GassuianMixture_Multi(tt.WB[,c(3:7)],3);
