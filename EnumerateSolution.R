EvaluationFunction = function(mu, covs, Gamma, X)
{
  sum = 0;
  K = nrow(mu);
  for(k in 1:K)
  {
    D = ncol(mu);    
    const = -(D/2)*log(2*pi) - 0.5*log(det(covs[,,k]));
    Xm = X - mu[k,];
    Xm = matrix(Xm,nrow=1,ncol=D)
    lg.px = const - 0.5 * diag(Xm %*% solve(covs[,,k]) %*% t(Xm));
    sum = sum + Gamma[k]*exp(lg.px);    
  }
  
  return((-1)*log(sum));
}

mu = matrix(c(75.49,107.10,106.09,112.36,183.08,161.51,402.74,518.00,319.98),ncol=3,byrow=TRUE);
covs = array(0,dim=c(3,3,3))
for(i in 1:3)
{
  covs[,,i] = matrix(c(4131.15,78.38,-135.91,75.38,4879.42,132.92,-135.91,132.92,4806.1859),ncol=3,byrow=TRUE);
}
Gamma = matrix(c(0.716,0.864,0.292,0.251,0.111,0.669,0.032,0.025,0.087),ncol=3,byrow=TRUE)

min = 999;
for(x in 1:600)
{
  for(y in 1:600)
  {
    z = 600 - x -y;
    value = EvaluationFunction(mu,covs,Gamma,c(x,y,z));
    
    if(value < min)
    {
      min = value;
      print(min);
    } 
  }
}