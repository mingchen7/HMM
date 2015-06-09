library(ggplot2)
library(lubridate)
library(mixtools)
library(MASS)

GassuianMixture = function(data, numCluster)
{
    # using Kmeans to cluster the data first  
    KmeansCluster = kmeans(data, numCluster, iter.max = 20, algorithm = 'Hartigan-Wong');
    
    Pi = c(); Mu = c(); Sigma = c();
    # Initialize pi_k, mu_k, sigma_k
    for(k in 1: numCluster)
    {
        Pi[k] = KmeansCluster$size[k] / length(data);
        Mu[k] = mean(data[KmeansCluster$cluster == k]);
        Sigma[k] = var(data[KmeansCluster$cluster == k]);
    }
    
    Q = matrix(data = NA, nrow = length(data), ncol = numCluster);
    LogLikelihood = c(0);
    
    for(iter in 1: 80)
    {
        # Maxization Step
        for(n in 1: length(data))
        {
            SumDenominator = 0;
            for(j in 1: numCluster)
            {
                Temp = Pi[j] * dnorm(data[n], Mu[j], sd = sqrt(Sigma[j]));
                SumDenominator = SumDenominator + Temp;
            }
            
            for(k in 1: numCluster)
            {
                Q[n, k] = Pi[k] * dnorm(data[n], Mu[k], sd = sqrt(Sigma[k])) / SumDenominator;
            }
        }
        
        Q[Q < 1e-60] = 1e-15;
        Q[Q > 1-1e-60] = 1e-15;
    
        
        # Expectation Step
        for(k in 1: numCluster)
        {
            Pi[k] = mean(Q[, k]);
            
            Mu[k] = sum(Q[, k] * t(data)) / sum(Q[, k]);
            Sigma[k] = sum(Q[, k] * t((data - Mu[k])^2)) / sum(Q[, k]);
            
        }
        
        # Likelihood
        LogLikelihood[iter] = 0;
        for(n in 1: length(data))
        {
            Temp = 0;
            for(k in 1: numCluster)
            {
                Temp = Pi[k] * dnorm(data[n], Mu[k], sd = sqrt(Sigma[k])) + Temp;
            }
            
            LogTemp = log(Temp); 
            LogLikelihood[iter] = LogLikelihood[iter] + LogTemp;
            
        }
        
        #print(sprintf('Pi = %0.2f, %0.2f   Mu = %0.2f, %0.2f    Sigma = %0.2f, %0.2f    LogLikelihood = %s', 
        #              Pi[1], Pi[2], Mu[1], Mu[2], Sigma[1], Sigma[2],  LogLikelihood[iter]));
        
        print(sprintf('iter = %d  loglikelihood = %f', iter, LogLikelihood[iter]));
#         print(Pi);
#         print(Mu);
#         print(Sigma);
    }
    
    
    #Plotting
    layout(matrix(c(1:4), 2, 2, byrow = TRUE));
    
    MinData = min(data); MaxData = max(data);
    XValue = seq(from = MinData, to = MaxData, by = 0.1);
    YValues = list(); MaxYValue = -1;
    for(k in 1: numCluster)
    {
        YValues[[k]] = dnorm(XValue, Mu[k], sd = sqrt(Sigma[k]));
        MaxYValue = max(MaxYValue, YValues[[k]]);
    }
    
    
    hist(data, nclass = 30);
    
    plot(x = XValue, y = YValues[[1]], type = 'l', ylim = c(0, MaxYValue));
    if(k >= 2)
    {
        for(k in 2: numCluster)
        {
            par(new = TRUE);
            plot(x = XValue, y = YValues[[k]], type = 'l', ylim = c(0, MaxYValue), axes = FALSE);
        }
        par(new = FALSE);
    }
    
    plot(x = 1:length(LogLikelihood), y = LogLikelihood, type = 'l');
    
    Estimation  = list();
    Estimation$Pi = Pi;
    Estimation$Mu = Mu;
    Estimation$Sigma = Sigma;
    Estimation$loglikelihood = LogLikelihood[length(LogLikelihood)];
    
    return(Estimation);
}

# example: 
setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
load('LinkTT_10hours.RData')
GassuianMixture(tt.WB$tt4,2)