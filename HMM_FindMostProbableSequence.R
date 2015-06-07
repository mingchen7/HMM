library(lubridate)
library(MASS)

evalEmissionProb = function(x,mu,sigma,i,k)
{
  # get the parameters of the ith variable
  # mu.i and sigma.i are K length vectors
  mu.i = mu[[i]];
  sigma.i = sigma[[i]];
  
  prob = dnorm(x[i],mu.i[k],sd=sqrt(sigma.i[k]));
  
  return(prob);
}

# recursion function for evalating Omega
evalOmega = function(x,mu,sigma,pi,A,i)
{   
  K = length(mu[[i]]);
  omega = rep(0,K); 
  psi.i = rep(0,K);  
  out = list();
  
  if(i == 1)
  {
    for(k in 1:K)
    {
      prob = evalEmissionProb(x,mu,sigma,i,k);
      omega[k] = log(pi[k]) + log(prob);
    }
    
    # create a lsit to store values of psi
    D = length(mu);
    psi = list();
    for(i in 1:D)
    {    
      psi[i] = rep(0,length(mu[[i]]));
    }
                
    out$psi = psi;
    out$omega = omega;
    return(out);
  }
  else
  {
    for(k in 1:K)
    {
      K.pre = length(mu[[i-1]]);      
      max = 0;
      tmp = 0;
      mark = 0;
      
      rslt = evalOmega(x,mu,sigma,pi,A,i-1);
      psi = rslt$psi;
      
      for(j in 1:K.pre)
      {
        tmp = log(A[[i-1]][j,k]) + rslt$omega[j];
        if(tmp >= max)
        {
          max = tmp;
          mark = j
        }
      }
      
      psi.i[k] = mark;      
      prob = evalEmissionProb(x,mu,sigma,i,k);
      omega[k] = log(prob) + max;      
    }

    psi[[i]] = psi.i;
    out$psi = psi;
    out$omega = omega;
    return(out);
  }
}

# Given a set of observations X, find the msot probable sequence of hidden state Z and the max joint probability
# Taking inputs of observations X and HMM parameters
Viterbi = function(X,mu,sigma,pi,A)
{
  Omega = c(0);
  D = length(mu);
  K = length(mu[[D]]);
    
  loglikelihood.max = 0;
  tmp = 0
  D.mark = 0;
  
  rslt = evalOmega(x,mu,sigma,pi,A,D,psi);
  psi = rslt$psi;
    
  #evalating the joint probablity X by evalating Omega(zD)
  for(k in 1:K)
  {
    tmp = rslt$omega[k];
    if(tmp >= loglikelihood.max)
    {
      loglikelihood.max = tmp;
      D.mark = k;
    }
  }
  
  # Back-tracking procedure
  z = D.mark;
  z.pre= 0;
  pritn(z);
  for(i in (D-1):1)
  {    
    z.pre = psi[[i]][z];
    z = z.pre;    
    print(z);
  }
  
  return(loglikelihood.max);
}