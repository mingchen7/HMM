library(ggplot2)
library(lubridate)
library(MASS)

# Evaluate the emission probability
evalEmissionProb = function(obs,mu,sigma,i,k)
{
  # get the parameters of the ith variable
  # mu.i and sigma.i are K length vectors
  mu.i = mu[[i]];
  sigma.i = sigma[[i]];
  
  prob = dnorm(obs[i],mu.i[k],sd=sqrt(sigma.i[k]));
  
  return(prob);
}

# Evaluating the alpha(Zi)  by recursion
# Taking index i representing the index of variable
# obs takes one row of data matrix
evalAlpha = function(obs,pi,mu,sigma,A,i)
{
  alpha.crt = c(0);
  output = list();
  
  # evaluating the 1st variable
  if(i == 1)
  {
    # scaling factor
    c = 0; 
    
    K = length(pi);
    for(k in 1:K)
    {
      alpha.crt[k] = pi[k] * evalEmissionProb(obs,mu,sigma,1,k);
      c = c + alpha.crt[k];
    }
    
    alpha.crt = alpha.crt / c;
    output$alpha = alpha.crt;
    output$c = c;
    return(output);
  }
  else
  {
    c = 0; #scaling factor
    K = length(mu[[i]]);
    A.crt = A[[i-1]]; #get transition matrix
    
    # call recursivelly the get the previous value of alpha    
    alpha.pre = evalAlpha(obs,pi,mu,sigma,A,i-1)$alpha;
    
    for(k in 1:K)
    {
      # get emmission probablity p(xi|zi)
      prob = evalEmissionProb(obs,mu,sigma,i,k);
      K.pre = length(alpha.pre);
      sum = 0;
      tmp = 0;
            
      for(j in 1:K.pre)
      {
        tmp = alpha.pre[j]*A.crt[j,k];
        sum = sum + tmp;        
      }
      
      alpha.crt[k] = prob * sum;      
      c = c + alpha.crt[k];
      
    }
    
    alpha.crt = alpha.crt / c;
    output$alpha = alpha.crt;
    output$c = c;
    return(output);    
  }
}

# Evaluating the Beta(Zi) by recursion
evalBeta = function(obs,mu,sigma,A,i,c)
{
  beta.crt = c(0);
  D = length(mu);
  if(i == D)
  {
    return(rep(1,length(mu[[i]])));
  }
  else
  {
    K = length(mu[[i]]);
    A.crt = A[[i]];
    
    # get the beta values of i+1
    beta.nxt = evalBeta(obs,mu,sigma,A,i+1,c);
    
    for(j in 1:K)
    {
      sum = 0;
      tmp = 0;
      K.nxt = length(beta.nxt);
      
      for(k in 1:K.nxt)
      {
        prob = evalEmissionProb(obs,mu,sigma,i+1,k);
        tmp = beta.nxt[k] * prob * A.crt[j,k] / c[i+1];
        sum = sum + tmp;        
      }
      
      beta.crt[j] = sum;
    }
    
    return(beta.crt);
  }
    
}


# Evaluating the Ksi(Zi-1,Zi) based on alpha(Zi-1) and beta(Zi)
evalKsi = function(obs,mu,sigma,A,i,alpha.pre,beta.crt,c)
{
  D = length(mu);
  
  if((i < 2) | (i > D))
  {
    print("Index i out of range!");
    return(NULL);
  }
  else
  {
    # the result is stored in a K.pre X K matrix
    K.pre = length(mu[[i-1]]);
    K = length(mu[[i]]);    
    ksi = matrix(0,nrow=K.pre,ncol=K);
        
    A.crt = A[[i-1]];    
          
    for(j in 1:K.pre)
    {
      for(k in 1:K)
      {
        prob = evalEmissionProb(obs,mu,sigma,i,k);
        
        # p(X) is not calcualted here as it will be canceled out in the maximization step
        ksi[j,k] = alpha.pre[j] * prob * A.crt[j,k] * beta.crt[k] / c[i];                
      }
    }
    
    return(ksi);
  }  
}

evalCompleteLikelihood = function(data,mu,sigma,pi,A,gamma,ksi)
{
  N = nrow(data);
  D = length(mu);
  
  term1 = 0;
  for(n in 1:N)
  {
    K = length(mu[[1]]);
    for(k in 1:K)
    {
      term1 = term1 + gamma[[1]][n,k] * log(pi[k]);
    }
  }
  
  term2 = 0;
  for(n in 1:N)
  {
    for(i in 2:D)
    {
      K2 = length(mu[[i]]);
      K1 = length(mu[[i-1]]);
      for(k in 1:K2)
      {
        for(j in 1:K1)
        {
          term2 = term2 + ksi[[i-1]][j,k,n] * log(A[[i-1]][j,k]);
        }
      }
    }
  }
  
  term3 = 0;
  for(n in 1:N)
  {
    for(i in 1:D)
    {
      K = length(mu[[i]]);
      for(k in 1:K)
      {
        prob = evalEmissionProb(data[n,],mu,sigma,i,k);
        term3 = term3 + gamma[[i]][n,k] * log(prob);
      }
    }
  }
  
  sum = term1 + term2 + term3;
  return(sum);
}

# Hidden Markov Model for sequential links using Bayesian Networks
# Taking N x D (# of observations x dimensions of variables), Each observatino is treated as independent
# Taking following initial parameters 
# mu: list of D vectors, initial mean for travel times
# sigma: list of D vectors, initial variance for travel times
# pi: Kx1 vector, initial prior probability for traffic states on the first link
# A: list of D-1 KxK matrices, transition probabilities

HiddenMarkovModel = function(data,mu.init,sigma.init,pi.init,A.init)
{
  data = as.matrix(data);
  N = nrow(data);
  D = ncol(data);
  K = c(0);
  
  # initialization
  mu = mu.init;
  sigma = sigma.init;
  pi = pi.init;
  A = A.init;
  
  for(i in 1:D)
  {
    # Number of components for each link
    K[i] = length(mu[[i]]); 
  }  
  
  gamma = list();
  alpha = list();
  beta = list();
  
  for(i in 1:D)
  {
    # dimension: DxNxK
    gamma[[i]] = matrix(0,nrow=N,ncol=K[i]);
    alpha[[i]] = matrix(0,nrow=N,ncol=K[i]);
    beta[[i]] = matrix(0,nrow=N,ncol=K[i]);
  }
  
  # scaling factor ci
  c = matrix(0,nrow=N,ncol=D);
  
  ksi = list();  
  for(i in 1:(D-1))
  {    
    #each observation is represented as a KxK matrix that sum to unity
    ksi[[i]] = array(0,dim=c(K[i],K[i+1],N)); 
  }
  
  LogLikelihood = c(0);
  CompLogLikelihood = c(0);
  
  # main loop for evaluating parameters
  for(iter in 1:100)
  {
    # E-step
    # Evaluate gamma and ksi
    for(n in 1:N)
    {
      obs = data[n,];      
      
      for(i in 1:D)
      {    
        #Evaluate alpha, beta        
        alpha.rst = evalAlpha(obs,pi,mu,sigma,A,i);
        alpha[[i]][n, ] = alpha.rst$alpha;
        c[n,i] = alpha.rst$c;
      }
      
      for(i in 1:D)
      {
        beta[[i]][n, ]= evalBeta(obs,mu,sigma,A,i,c[n,]);                                     
        
        #use alpha and beta to calcualte gamma and ksi
        gamma[[i]][n, ] = alpha[[i]][n, ] * beta[[i]][n, ];        
        
        if(i >= 2)
        {
          ksi[[i-1]][,,n] = evalKsi(obs,mu,sigma,A,i,alpha[[i-1]][n,],beta[[i]][n,],c[n,]);  
        }          
      }
    }
    
    # tricky control
#     for(i in 1:D)
#     {
#       gamma[[i]][gamma[[i]] < 1e-60] = 1e-15;
#       gamma[[i]][gamma[[i]] > 1-1e-60] = 1e-15;      
#       
#       if(i < D)
#       {
#         ksi[[i]][ksi[[i]] < 1e-60] = 1e-15;
#         ksi[[i]][ksi[[i]] > 1-1e-60] = 1e-15;  
#       }      
#     }
              
    
    # M-step
    # update mu, sigma, pi, and A
    for(i in 1:D)
    {
      for(k in 1:K[i])
      {
        mu[[i]][k] = sum(gamma[[i]][,k] * data[,i]) / sum(gamma[[i]][,k]);
        sigma[[i]][k] = sum(gamma[[i]][,k] * t((data[,i]-mu[[i]][k])^2)) / sum(gamma[[i]][,k]);
                                    
        if(i == 1)
        {          
          pi[k] = sum(gamma[[1]][,k]) / sum(gamma[[1]][,]);
        }
        else #i>=2
        {
          for(j in 1:K[i-1])
          {
            # summation over n
            A[[i-1]][j,k] = apply(ksi[[i-1]],1:2,sum)[j,k] / sum(apply(ksi[[i-1]],1:2,sum)[j,]);
          }          
        }
                            
      }
    }
        
    # Evaluate log-likelihood
    LogLikelihood[iter] = 0;
    # CompLogLikelihood[iter] = 0;
    tmp = 0;
    for(n in 1:N)
    { 
      for(i in 1:D)
      {
        tmp = tmp + log(c[n,i]);
      }            
    }
    
    LogLikelihood[iter] = tmp;            
    # CompLogLikelihood = evalCompleteLikelihood(data,mu,sigma,pi,A,gamma,ksi);
        
    print(sprintf('iter = %d  loglikelihood = %f', iter, LogLikelihood[iter]));
    print("Pi:")
    print(pi);
    print("Mu:");
    print(mu);
    print("Sigma:");
    print(sigma);
    print("Transition matrix:");
    print(A);
  }
  
  # plot log-likelihood
  plot(x = 1:length(LogLikelihood), y = LogLikelihood, xlab='iteration',type = 'l');
  
  results  = list();
  results$pi = pi;
  results$mu = mu;
  results$sigma = sigma;
  results$A.transition = A;
  results$loglikelihood = LogLikelihood[length(LogLikelihood)];
  
  return(results);    
}                      

# example: 
setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");

# loading data
load('LinkTT_10hours.RData')
data = tt.WB[,c(3:7)];

# initializing parameters
mu.init = list();
mu.init[[1]] = c(89.19,39.76,58.95);
mu.init[[2]] = c(87.01,58.93,43.84);
mu.init[[3]] = c(58.78,39.60,35.60);
mu.init[[4]] = c(106.65,81.66,39.64);
mu.init[[5]] = c(98.76526,170.08265);

sigma.init = list()
sigma.init[[1]] = c(102.749,3.165,123.500);
sigma.init[[2]] = c(252.64,24.81,5.39);
sigma.init[[3]] = c(59.14,8.61,2.78);
sigma.init[[4]] = c(34.10,498.01,14.03);
sigma.init[[5]] = c(36.66009,72.18567);

pi.init = c(0.345,0.238,0.417);

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

HiddenMarkovModel(data, mu.init, sigma.init, pi.init, A.init);

