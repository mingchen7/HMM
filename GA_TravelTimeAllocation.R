# Evaluate log-likelihood function
# mu: kxM matrix
# sigma: MxMxk array
# Gamma: 1xM vector
# x: 1xM vector

EvaluationFunction.MultiGMM = function(mu, sigma, Gamma, X)
{
  sum = 0;
  K = nrow(mu);
  for(k in 1:K)
  {
    D = ncol(mu);    
    const = -(D/2)*log(2*pi) - 0.5*log(det(sigma[,,k]));
    Xm = X - mu[k,];
    Xm = matrix(Xm,nrow=1,ncol=D)
    lg.px = const - 0.5 * diag(Xm %*% solve(sigma[,,k]) %*% t(Xm));
    sum = sum + Gamma[k]*exp(lg.px);    
  }
  
  return((-1)*log(sum));
}

EvaluationFunction.UniGMM = function(mu,sigma,Gamma,X)
{
  Sum = 0;
  D = length(mu);
  for(i in 1:D)
  {
    K = length(mu[[i]]);
    for(k in 1:K)
    {
      Sum = Sum + Gamma[[i]][k] * dnorm(X[i], mu[[i]][k], sd = sqrt(sigma[[i]][k]));
    }
  }
  
  LogLikelihood = log(Sum);
  return(-LogLikelihood);
}


GA.MainFunction = function(mu,sigma,Gamma,y.obs,model)
{
    set.seed(y.obs);    
    iter = 1;    
    tol = 1e-5;
    miniter = 100;
    maxiter = 100;
    objValue = 0;
    
    if(model =='Univariate')
    {
      D = length(mu);  
    }
    else
    {
      D = ncol(mu); # get number of links  
    }
    
    InitialGen = matrix(nrow=50, ncol=D);
            
    for(i in 1:2)
    {
      rnd = runif(D);
      rnd = (rnd / sum(rnd)) * y.obs;
      InitialGen[i,] = rnd;
    }
    
    if(model =='Univariate')
    {
      objValue.new = EvaluationFunction.UniGMM(mu, sigma, Gamma, InitialGen[1,])
    }
    else
    {
      objValue.new = EvaluationFunction.MultiGMM(mu, sigma, Gamma, InitialGen[1,]);
    }                
        
    
    while(((objValue.new - objValue) > tol) | (iter < miniter))
    {
      InitialGen = GA.NewGeneration(mu,sigma,Gamma,InitialGen,y.obs,D,model);
      objValue = objValue.new;
      
      if(model =='Univariate')
      {
        objValue.new = EvaluationFunction.UniGMM(mu, sigma, Gamma, InitialGen[1,]);
      }
      else
      {
        objValue.new = EvaluationFunction.MultiGMM(mu, sigma, Gamma, InitialGen[1,]);
      }        
      
      iter = iter + 1;
      #print(c(iter,objValue.new));
      
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


GA.NewGeneration = function(mu,sigma,Gamma,originalGen,y.obs,D,model)
{
    NewGeneration = c(); #matrix(data = -999, nrow = nrow(originalGens), ncol = ncol(originalGens));
    
    SortedFitness = GA.Evaluate(mu,sigma,Gamma,originalGen,model);
    SortedOriginalGen = GA.SortGen(originalGen,SortedFitness$ix);
    
    NewGeneration = GA.Elitism(SortedOriginalGen);
    #D = ncol(mu);
    
    for(i in 1:((nrow(originalGen))/ 2 - 1))
    {
        SelectedChromos = GA.SelectedChromosomes(SortedOriginalGen);
        
        CrossedChromos = GA.Crossover(SelectedChromos[,c(1:D)],y.obs);
        MutatedChromos = GA.Mutation(CrossedChromos[,c(1:D)],y.obs);
        
        NewGeneration = rbind(NewGeneration, MutatedChromos);
    }
    
    return(NewGeneration);
}

GA.Evaluate = function(mu,sigma,Gamma,originalGen,model)
{
    Fitness = c();
    
    for(i in 1: nrow(originalGen))
    {
      if(model =='Univariate')
      {
        Fitness[i] = EvaluationFunction.UniGMM(mu,sigma,Gamma,originalGen[i,]);
      }
      else
      {
        Fitness[i] = EvaluationFunction.MultiGMM(mu,sigma,Gamma,originalGen[i,]);
      }                     
    }
    
    SortedFitness = sort(Fitness, index.return = TRUE);
    #print(SortedFitness$x);
    
    return(SortedFitness);
}


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

# example: GA.RankSelection(InitialGen)
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

## TEST DATA 
# multivariate
# mu = matrix(c(64.39,84.24,36.44,53.22,132.87,71.31,81.48,36.88,40.17,96.70,60.67,64.34,49.76,102.36,103.12),ncol=5,byrow=TRUE);
# sigma = array(0,dim=c(5,5,3))
# sigma[,,1] = matrix(c(674.54,-327.20,2.168,-16.58,177.27,-327.20,484.59,4.933,47.057,-63.86,2.168,4.933,9.412,-4.391,16.904,-16.583,47.057,-4.391,80.320,-168.508,177.270,-63.866,16.904,-168.508,1569.619),ncol=5,byrow=TRUE);
# sigma[,,2] = matrix(c(607.875,-171.285,-1.695,2.683,14.407,-171.285,416.381,1.495,20.225,38.613,-1.695,1.495,14.424,-4.294,-6.834,2.683,20.225,-4.294,18.468,4.507,14.407,38.613,-6.834,4.507,33.326),ncol=5,byrow=TRUE);
# sigma[,,3] = matrix(c(556.202,-140.974,-4.953,6.083,-4.428,-140.974,371.512,44.138,22.477,-6.369,-4.953,44.138,179.011,-146.354,-15.017,6.0873,22.477,-146.354,200.895,40.886,-4.428,-6.369,-15.017,40.886,170.937),ncol=5,byrow=TRUE)
# Gamma = c(0.077,0.415,0.506);
# 
# mu2 = matrix(c(65.97,75.896,39.951,67.415,169.657,72.181,79.912,35.830,39.328,97.040,36.707,68.785,50.249,100.477,100.357,60.403,86.392,38.954,46.279,94.444,72.074,62.511,49.582,102.560,101.421),ncol=5,byrow=TRUE);
# sigma2 = array(0,dim=c(5,5,5));
# sigma2[,,1] = matrix(c(614.559,-210.098,-29.355,-181.795,3.713,-210.098,502.008,-37.279,-288.513,23.904,-29.244,-37.279,93.794,100.081,10.500,-181.795,-288.513,100.081,1043.551,-21.571,3.713,23.904,10.500,-21.571,100.225),ncol=5,byrow=TRUE)
# sigma2[,,2] = matrix(c(590.936,-138.961,-2.753,9.003,15.033,-138.961,409.852,-2.174,15.753,37.589,-2.753,-2.174,3.704,-1.335,-1.640,9.003,15.753,-1.335,10.261,3.842,15.033,37.589,-1.640,3.842,28.220),ncol=5,byrow=TRUE);
# sigma2[,,3] = matrix(c(12.151,6.799,0.763,0.876,1.202,6.799,356.246,31.066,15.291,8.628,0.763,31.066,191.711,-147.743,-10.002,0.876,15.291,-147.743,177.549,2.668,1.202,8.628,-10.002,2.668,25.406),ncol=5,byrow=TRUE);
# sigma2[,,4] = matrix(c(673.442,-269.773,17.786,-30.793,12.824,-269.773,427.453,-6.152,30.485,45.030,17.786,-6.152,29.649,-28.981,-10.831,-30.793,30.485,-28.981,77.559,-11.433,12.824,45.030,-10.831,-11.433,56.570),ncol=5,byrow=TRUE);
# sigma2[,,5] = matrix(c(416.181,-140.444,-1.533,-9.996,-2.728,-140.444,371.397,49.734,36.651,3.475,-1.533,49.734,170.750,-142.195,-15.053,-9.996,36.651,-142.195,191.507,18.742,-2.728,3.475,-15.053,18.742,22.505),ncol=5,byrow=TRUE);
# Pi2 = c(0.056,0.310,0.157,0.141,0.336);

#univariate
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

# GA.MainFunction(mu,sigma,Pi,244.7,'Univariate')
# GA.MainFunction(mu2,sigma2,Pi2,244.7,'Multivariate')