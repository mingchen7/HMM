# **IMPORTANT**
# Functions to evaluate fitness
# Depending on different models

EvaluationFunction.UniGMM = function(mu,sigma,Gamma,x)
{
  Sum = 0;
  D = length(mu);
  for(i in 1:D)
  {
    K = length(mu[[i]]);
    for(k in 1:K)
    {
      Sum = Sum + Gamma[[i]][k] * dnorm(x[i], mu[[i]][k], sd = sqrt(sigma[[i]][k]));
    }
  }
  
  LogLikelihood = log(Sum);
  return(-LogLikelihood);
}

EvaluationFunction.MultiGMM = function(mu, sigma, Gamma, x)
{
  sum = 0;
  K = nrow(mu);
  for(k in 1:K)
  {
    D = ncol(mu);    
    const = -(D/2)*log(2*pi) - 0.5*log(det(sigma[,,k]));
    Xm = x - mu[k,];
    Xm = matrix(Xm,nrow=1,ncol=D)
    lg.px = const - 0.5 * diag(Xm %*% solve(sigma[,,k]) %*% t(Xm));
    sum = sum + Gamma[k]*exp(lg.px);    
  }
  
  return((-1)*log(sum));
}


EvaluationFunction.HMM = function(mu,sigma,Pi,A,x)
{
  LogLikelihood = -Inf;
  LogLikelihood = Viterbi(x,mu,sigma,Pi,A);
  return(-LogLikelihood);
}

# THIS FUNCTION STARTS THE GA ALGORITHM
# mu: mean
# sigma: variance/covariance
# Gamma: pi
# A:transition matrix for HMM
GA.MainFunction = function(mu,sigma,Gamma,A,y.obs,model)
{
    set.seed(y.obs);    
    iter = 1;    
    tol = 1e-5;
    miniter = 100;
    maxiter = 100;
    objValue = 0;
    
    if((model =='Univariate') | (model == 'HMM'))
    {
      D = length(mu);  
    }
    else
    {
      D = ncol(mu); # get number of links  
    }
    
    InitialGen = matrix(nrow=50, ncol=D);
    
    # initialize a population with size of 50
    for(i in 1:50)
    {
      rnd = runif(D);
      rnd = (rnd / sum(rnd)) * y.obs;
      InitialGen[i,] = rnd;
    }
    
    if(model =='Univariate')
    {
      objValue.new = EvaluationFunction.UniGMM(mu,sigma,Gamma,InitialGen[1,])
    }
    else if(model == 'Multivariate')
    {
      objValue.new = EvaluationFunction.MultiGMM(mu,sigma,Gamma,InitialGen[1,]);
    }      
    else
    {
      objValue.new = EvaluationFunction.HMM(mu,sigma,Gamma,A,InitialGen[1,]);
    }
        
    
    while(((objValue.new - objValue) > tol) | (iter < miniter))
    {
      # Call for generating a new population
      InitialGen = GA.NewGeneration(mu,sigma,Gamma,A,InitialGen,y.obs,D,model);
      objValue = objValue.new;
      
      if(model =='Univariate')
      {
        objValue.new = EvaluationFunction.UniGMM(mu,sigma,Gamma,InitialGen[1,]);
      }
      else if(model == 'Multivariate')
      {
        objValue.new = EvaluationFunction.MultiGMM(mu,sigma,Gamma,InitialGen[1,]);
      }      
      else
      {
        objValue.new = EvaluationFunction.HMM(mu,sigma,Gamma,A,InitialGen[1,]);
      }
      
      iter = iter + 1;
      # print(c(iter,objValue.new));
      
      if(iter > maxiter)
      {
        break;
      }
    }
            
    output = list();
    output$solution = InitialGen[1,];
    output$loglikelihood = -objValue.new;
    return(output);
}


GA.NewGeneration = function(mu,sigma,Gamma,A,originalGen,y.obs,D,model)
{
    NewGeneration = c();
    
    SortedFitness = GA.Evaluate(mu,sigma,Gamma,A,originalGen,model);
    SortedOriginalGen = GA.SortGen(originalGen,SortedFitness$ix);
    
    NewGeneration = GA.Elitism(SortedOriginalGen);    
    
    for(i in 1:((nrow(originalGen))/ 2 - 1))
    {
        SelectedChromos = GA.SelectedChromosomes(SortedOriginalGen);
        
        CrossedChromos = GA.Crossover(SelectedChromos[,c(1:D)],y.obs);
        MutatedChromos = GA.Mutation(CrossedChromos[,c(1:D)],y.obs);
        
        NewGeneration = rbind(NewGeneration, MutatedChromos);
    }
    
    return(NewGeneration);
}

GA.Evaluate = function(mu,sigma,Gamma,A,originalGen,model)
{
    Fitness = c();
    
    for(i in 1: nrow(originalGen))
    {
      if(model =='Univariate')
      {
        Fitness[i] = EvaluationFunction.UniGMM(mu,sigma,Gamma,originalGen[i,]);
      }
      else if(model == 'Multivariate')
      {
        Fitness[i] = EvaluationFunction.MultiGMM(mu,sigma,Gamma,originalGen[i,]);
      }      
      else
      {
        Fitness[i] = EvaluationFunction.HMM(mu,sigma,Gamma,A,originalGen[i,]);
      }
    }
    
    SortedFitness = sort(Fitness, index.return = TRUE);
    #print(SortedFitness$x);
    
    return(SortedFitness);
}

# General GA functions: independent of any model
GA.Elitism = function(sortedGeneration) # from min to max
{
    return(sortedGeneration[c(1,2), ]);
}

GA.SortGen = function(originalGen, sortIndex)
{
    return(originalGen[sortIndex, ]);
}

GA.SelectedChromosomes = function(sortedGeneration)
{
    FirstChromoIndex = GA.RankSelection(sortedGeneration);
    SecondChromoIndex = GA.RankSelection(sortedGeneration);
    
    #print(sprintf('%02d, %02d', FirstChromoIndex, SecondChromoIndex));
    
    return(sortedGeneration[c(FirstChromoIndex, SecondChromoIndex), ]);
}

GA.Mutation = function(crossedChromosomes, y.obs)
{
    L = ncol(crossedChromosomes);  
    subChromosomes = crossedChromosomes[,c(1:L-1)];
    L.sub = ncol(subChromosomes);
    MutationProbability = 0.03;
    RandomProbability = runif(1);
        
    if(RandomProbability <= MutationProbability)
    {
      for(i in 1:2)
      {
        RanMutation = runif(1); # randomly determine which mutation method to use
        if(RanMutation <= 0.7) # uniform mutation
        {
          
          j = sample(1:L.sub,1);
          # calculated upperbound and lowerboudn for gene j
          lowerbnd = 0.1;
          upperbnd = 0.999999*(y.obs - sum(subChromosomes[i,-j]));                    
          
          gene.new = runif(1,lowerbnd,upperbnd);
          subChromosomes[i,j] = gene.new;
        }
        else # boundary mutation
        {
          j = sample(1:L.sub,1);
          lowerbnd = 0.1;
          upperbnd = 0.999999*(y.obs - sum(subChromosomes[i,-j]));
          if(runif(1)>0.5)
          {
            gene.new = upperbnd;
          }
          else
          {
            gene.new = lowerbnd;
          }
          subChromosomes[i,j] = gene.new;
        }                        
      }
      
      # Complete the chromosomes with last element
      crossedChromosomes[,c(1:L.sub)] = subChromosomes[,]
      crossedChromosomes[1,L] = y.obs - sum(subChromosomes[1,]);
      crossedChromosomes[2,L] = y.obs - sum(subChromosomes[2,]);
    }
    
    return(crossedChromosomes);     
}

Generate.Sign = function()
{
    RandomVar = runif(1, min = -1, max = 1);
    
    return(RandomVar / abs(RandomVar));
}

GA.Crossover = function(selectedChromosomes, y.obs)
{        
    L = ncol(selectedChromosomes);
    subChromosomes = selectedChromosomes[,c(1:L-1)]; #get the first 1 to L-1 genes
    L.sub = ncol(subChromosomes);
    cross.pnt = ceiling(2*L.sub/3);    
    CrossoverProbability = 0.5;
    
    RandomProbability = runif(1);
    
    if(RandomProbability <= CrossoverProbability)
    {
#         print("======= Selected chromosomes ===========");
#         print(selectedChromosomes);
        
        RandomCross= runif(1); #determine which crossover method                
        if(RandomCross <= 0.33) #Simple Crossover
        {
#            print("coming in simple cross");
                      
           a= seq(1,0,-0.05); 
           for(i in 1:length(a))
           {
             # make copy to avoid overwrite
             v.copy = subChromosomes[1,];
             w.copy = subChromosomes[2,];
             
             Sv.cross = a[i]*subChromosomes[2,c(cross.pnt:L.sub)] + (1-a[i]) * subChromosomes[1,c(cross.pnt:L.sub)];
             Sw.cross = a[i]*subChromosomes[1,c(cross.pnt:L.sub)] + (1-a[i]) * subChromosomes[2,c(cross.pnt:L.sub)];
             v.copy[c(cross.pnt:L.sub)] = Sv.cross;
             w.copy[c(cross.pnt:L.sub)] = Sw.cross;
             
#              cat('a=',a[i],'\n');
             
             if((sum(v.copy) < y.obs) & (sum(w.copy) < y.obs))
             {
               a.largest = a[i]
               break; # found the largest a
             }                          
           }
           
           #update new chromosomes
           subChromosomes[1,c(cross.pnt:L.sub)] = Sv.cross;
           subChromosomes[2,c(cross.pnt:L.sub)] = Sw.cross;
        }
        else if(RandomCross <=0.66) # single arithmetical crossover
        {
#           print("Coming in single cross");          
          
          j = sample(1:L.sub,1); 
          v.lowerbnd = 0.1;
          v.upperbnd = 0.999999*(y.obs - sum(subChromosomes[1,-j]));
          w.lowerbnd = 0.1;
          w.upperbnd = 0.999999*(y.obs - sum(subChromosomes[2,-j]));
                              
          range = getFeasibleRange(w.lowerbnd,v.lowerbnd,w.upperbnd,v.upperbnd,subChromosomes[1,j],subChromosomes[2,j]);          
                    
#           cat('Bounds:',c(v.upperbnd,w.upperbnd),'\n');
#           cat('Range:',range,'\n');
#           cat('j=',j,'\n');
          
          a = runif(1,range[1],range[2]);          
          vj.new = a*subChromosomes[2,j] + (1-a)*subChromosomes[1,j];
          wj.new = a*subChromosomes[1,j] + (1-a)*subChromosomes[2,j];            
          
#           cat('vj=',vj.new,'wj=',wj.new,'\n');
          
          subChromosomes[1,j] = vj.new;
          subChromosomes[2,j] = wj.new;                    
          
        }
        else # whole arithmetical crossover
        {
#           print("Coming in whole crossover");
          
          a = runif(1);
          tmp = subChromosomes[1,];
          subChromosomes[1,] = a*subChromosomes[2,] + (1-a)*subChromosomes[1,];
          subChromosomes[2,] = a*tmp + (1-a)*subChromosomes[2,];          
        }
        
        # last element of chromosome
#         print(subChromosomes);
#         print((sum(subChromosomes[1,]) < y.obs));
#         print((sum(subChromosomes[2,]) < y.obs));

        selectedChromosomes[,c(1:L.sub)] = subChromosomes[,]
        selectedChromosomes[1,L] = y.obs - sum(subChromosomes[1,]);
        selectedChromosomes[2,L] = y.obs - sum(subChromosomes[2,]);
    }

    return(selectedChromosomes);    
}

getFeasibleRange = function(w.LB,v.LB,w.UB,v.UB,vj,wj)
{
  if( abs(vj - wj) < 1e-6 )
  {
    range.min = 0;
    range.max = 0;
  }  
  else if(vj > wj)
  {
    range.min = max((w.LB-wj)/(vj-wj),(v.UB-vj)/(wj-vj));
    range.max = min((v.LB-vj)/(wj-vj),(w.UB-wj)/(vj-wj));
  }  
  else if(vj < wj)
  {
    range.min = max((v.LB-vj)/(wj-vj),(w.UB-wj)/(vj-wj));
    range.max = min((w.LB-wj)/(vj-wj),(v.UB-vj)/(wj-vj));
  }
  else
  {
    range.min = 0;
    range.max = 0;
  }
  
  # control for the range
  if(range.min > range.max)
  {
    range.min = 0;
    range.max = 0;
  }
#   if(range.min > range.max)
#   {
#     cat('vj=',vj,'\n','wj=',wj,'\n','v.low=',v.LB,'\n','v.up',v.UB,'\n','w.low=',w.LB,'\n','w.up=',w.UB,'\n');
#   }
 
  return(c(range.min,range.max));
}

GA.RankSelection = function(generation)
{
    LabelVector = c();
    for(i in 1: nrow(generation))
    {
        for(j in (sum(1: (i-1)) + 1): sum(1:i))
        {
            LabelVector[j] = i;
        }
    }
    
    RandomVar = runif(1);
    
    By = 1 / sum(1: nrow(generation));
    Index = ceiling(RandomVar / By);
    
    return(nrow(generation) - LabelVector[Index] + 1);
}


# setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM");
# source('Parameters.R');
# source('HMM_FindMostprobableSequence.R');

# examples
# GA.MainFunction(mu,sigma,Pi,NULL,361.3,'Univariate')
# GA.MainFunction(mu2,sigma2,Pi2,NULL,361.3,'Multivariate')
# GA.MainFunction(mu3,sigma3,Pi3,A,361.3,'HMM')