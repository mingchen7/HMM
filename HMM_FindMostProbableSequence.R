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
  omega = rep(0,K);  # to save the results for omega at local step i 
  psi.i = rep(0,K);  # to save the results of psi at local step i
  out = list();
  
  if(i == 1)
  {
    for(k in 1:K)
    {
      prob = evalEmissionProb(x,mu,sigma,i,k);
      omega[k] = log(pi[k]) + log(prob);
    }
    
    # create a list to store values of psi
    D = length(mu);    
    psi = list();
    for(i in 1:(D-1))
    {
      NumComps = length(mu[[i+1]]);
      psi[[i]] = rep(0,NumComps);
    }
                
    out$psi = psi;
    out$omega = omega;
    return(out);
  }
  else # i >= 2
  {
    K.pre = length(mu[[i-1]]);          
    rslt = evalOmega(x,mu,sigma,pi,A,i-1);
    psi = rslt$psi;
    for(k in 1:K)
    {      
      max = -Inf;
      tmp = 0;
      mark = 0;
      
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

    psi[[i-1]] = psi.i;
    out$psi = psi;
    out$omega = omega;
    return(out);
  }
}

# Given a set of observations X, find the msot probable sequence of hidden state Z and the max joint probability
# Taking inputs of observations X and HMM parameters
Viterbi = function(x,mu,sigma,pi,A)
{
  Omega = c(0);
  D = length(mu);
  K = length(mu[[D]]);
    
  loglikelihood.max = -Inf;
  tmp = 0
  D.mark = 0;
  
  rslt = evalOmega(x,mu,sigma,pi,A,D);
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
  # cat('traffic state of link 5 is',z,'\n');
  for(i in (D-1):1)
  {    
    z.pre = psi[[i]][z];
    z = z.pre;    
    # cat('traffic state of link',i,'is',z,'\n');    
  }
  
  return(loglikelihood.max);
}

# this function use enumeration to list all the possible solutions 
# .. and verify the correctness of Viterbi algorithm
Enumeration = function(x,mu,sigma,pi,A)
{
  D = length(mu);
  K = c(0);
  for(i in 1:D)
  {
    K[i] = length(mu[[i]]);
  }
  
  tmp = 0;
  max = -Inf;
  term1 = 0;
  term2 = 0;
  term3 = 0;
  for(z1 in 1:K[1])
  {
    for(z2 in 1:K[2])
    {
      for(z3 in 1:K[3])
      {
        for(z4 in 1:K[4])
        {
          for(z5 in 1:K[5])
          {
            term1 = log(pi[z1]);
            term2 = log(A[[1]][z1,z2]) + log(A[[2]][z2,z3]) + log(A[[3]][z3,z4]) + log(A[[4]][z4,z5]);
            term3 = log(evalEmissionProb(x,mu,sigma,1,z1)) + log(evalEmissionProb(x,mu,sigma,2,z2)) + 
              log(evalEmissionProb(x,mu,sigma,3,z3)) + log(evalEmissionProb(x,mu,sigma,4,z4)) + 
              log(evalEmissionProb(x,mu,sigma,5,z5));              
            
            tmp = term1 + term2 + term3;
            cat('Sequence:',z1,z2,z3,z4,z5,'LL:',tmp,'\n');
            if(tmp >= max)
            {
              max = tmp;              
            }
          }
        }
      }
    }
  }
  
  return(max);
}

# FOR TEST
# setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
# load('LinkTT_10hours.RData');
# source('Parameters.R');
# 
# # x=c(115,77.2,56.2,34.4,78.4);
# x=c(36.37660,98.23747,36.23666,89.58783,100.86144);
# 
# Viterbi(x,mu3,sigma3,Pi3,A);
# Enumeration(x,mu3,sigma3,Pi3,A);
