setwd("C:\\Users\\mingchen7\\Documents\\GitHub\\HMM")

library(ggplot2)
library(reshape2)
library(gridExtra)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

AggregateData = function(results){
  MAE.Benchmark = list()
  MAE.UniGMM = list()
  MAE.MultiGMM = list()
  MAE.HMM = list()
  
  MAPE.Benchmark = list()
  MAPE.UniGMM = list()
  MAPE.MultiGMM = list()
  MAPE.HMM = list()
  
  RMSE.Benchmark = list()
  RMSE.UniGMM = list()
  RMSE.MultiGMM = list()
  RMSE.HMM = list()
  
  N.Alg = 3
  
  MAE = matrix(NA,nrow=10,ncol=N.Alg)
  MAPE = matrix(NA,nrow=10,ncol=N.Alg)
  RMSE = matrix(NA,nrow=10,ncol=N.Alg)
  
  for(run in 1:10)
  {
    MAE.Benchmark[[run]] = abs(results[[run]]$x.benchmark - results[[run]]$x.obs)
    MAE.UniGMM[[run]] = abs(results[[run]]$x.UniGMM - results[[run]]$x.obs)
    #MAE.MultiGMM[[run]] = abs(results[[run]]$x.MultiGMM - results[[run]]$x.obs)
    MAE.HMM[[run]] = abs(results[[run]]$x.HMM - results[[run]]$x.obs)
    
    MAPE.Benchmark[[run]] = abs(results[[run]]$x.benchmark - results[[run]]$x.obs) / results[[run]]$x.obs
    MAPE.UniGMM[[run]] = abs(results[[run]]$x.UniGMM - results[[run]]$x.obs) / results[[run]]$x.obs
    #MAPE.MultiGMM[[run]] = abs(results[[run]]$x.MultiGMM - results[[run]]$x.obs) / results[[run]]$x.obs
    MAPE.HMM[[run]] = abs(results[[run]]$x.HMM - results[[run]]$x.obs) / results[[run]]$x.obs
    
    RMSE.Benchmark[[run]] = (results[[run]]$x.benchmark - results[[run]]$x.obs)^2
    RMSE.UniGMM[[run]] = (results[[run]]$x.UniGMM - results[[run]]$x.obs)^2
    #RMSE.MultiGMM[[run]] = (results[[run]]$x.MultiGMM - results[[run]]$x.obs)^2
    RMSE.HMM[[run]] = (results[[run]]$x.HMM - results[[run]]$x.obs)^2  
    
    # calculate mean of MAEs, MAPEs and RMSEs
    MAE[run,1] = mean(as.matrix(MAE.Benchmark[[run]]))
    MAE[run,2] = mean(as.matrix(MAE.UniGMM[[run]]))
    #MAE[run,3] = mean(as.matrix(MAE.MultiGMM[[run]]))
    MAE[run,3] = mean(as.matrix(MAE.HMM[[run]]))
    
    MAPE[run,1] = mean(as.matrix(MAPE.Benchmark[[run]]))
    MAPE[run,2] = mean(as.matrix(MAPE.UniGMM[[run]]))
    #MAPE[run,3] = mean(as.matrix(MAPE.MultiGMM[[run]]))
    MAPE[run,3] = mean(as.matrix(MAPE.HMM[[run]]))
    
    RMSE[run,1] = sqrt(mean(as.matrix(RMSE.Benchmark[[run]])))
    RMSE[run,2] = sqrt(mean(as.matrix(RMSE.UniGMM[[run]])))
    #RMSE[run,3] = sqrt(mean(as.matrix(RMSE.MultiGMM[[run]])))
    RMSE[run,3] = sqrt(mean(as.matrix(RMSE.HMM[[run]])))
  }
  
  MAE = as.data.frame(MAE)
  MAPE = as.data.frame(MAPE)
  RMSE = as.data.frame(RMSE)
  MAE$run = c(1:10)
  MAPE$run = c(1:10)
  RMSE$run = c(1:10)
  
  colnames(MAE) = c("Benchmark","GMM","HMM","run")
  colnames(MAPE) = c("Benchmark","GMM","HMM","run")
  colnames(RMSE) = c("Benchmark","GMM","HMM","run")
  
  out = list()
  out$MAE = MAE
  out$MAPE = MAPE
  out$RMSE = RMSE
  out$MAPE.Benchmark = MAPE.Benchmark
  out$MAPE.UniGMM = MAPE.UniGMM
  out$MAPE.HMM =MAPE.HMM
  return(out)
}

plotME = function(MAE,MAPE,RMSE){
  df.MAE = melt(MAE,id=c("run"),variable.name="Algorithm",value.name="MAE")
  df.MAPE = melt(MAPE,id=c("run"),variable.name="Algorithm",value.name="MAPE")
  df.RMSE = melt(RMSE,id=c("run"),variable.name="Algorithm",value.name="RMSE")
  
#   p1 = ggplot(df.MAE,aes(x=Algorithm,y=MAE)) + theme_bw() +
#     geom_boxplot()+ ylab("MAE") + 
#     theme(title=element_text(face="bold",size=16),      
#           axis.title.y = element_text(face="bold", size=16),
#           axis.text = element_text(colour="#000000", size=16),
#           legend.title = element_text(size=16),
#           legend.text = element_text(size = 16));    
#   
#   p2 = ggplot(df.MAPE,aes(x=Algorithm,y=MAPE)) + theme_bw() +
#     geom_boxplot()+ ylab("MAPE") + 
#     theme(title=element_text(face="bold",size=16),      
#           axis.title.y = element_text(face="bold", size=16),
#           axis.text = element_text(colour="#000000", size=16),
#           legend.title = element_text(size=16),
#           legend.text = element_text(size = 16));    
#   
#   p3 = ggplot(df.RMSE,aes(x=Algorithm,y=RMSE)) + theme_bw() +
#     geom_boxplot()+ ylab("RMSE") + 
#     theme(title=element_text(face="bold",size=16),      
#           axis.title.y = element_text(face="bold", size=16),
#           axis.text = element_text(colour="#000000", size=16),
#           legend.title = element_text(size=16),
#           legend.text = element_text(size = 16));   
#   
#   grid.arrange(p1,p2,p3,ncol=3);
  
  df.MAE.err = summarySE(df.MAE,measurevar="MAE",groupvars=c("Algorithm"))
  df.MAPE.err = summarySE(df.MAPE,measurevar="MAPE",groupvars=c("Algorithm"))
  df.RMSE.err = summarySE(df.RMSE,measurevar="RMSE",groupvars=c("Algorithm"))

  p1 = ggplot(df.MAE.err, aes(x=Algorithm, y=MAE, colour="red",group=1)) + 
    geom_errorbar(aes(ymin=MAE-sd, ymax=MAE+sd), colour="black", width=.1) +
    geom_line(color="red",size=1.2) +
    geom_point(size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Algorithm") + ylab("MAE") + 
    theme_bw() +
    theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.position="none"); 

  p2 = ggplot(df.MAPE.err, aes(x=Algorithm, y=MAPE, colour="blue",group=1)) + 
    geom_errorbar(aes(ymin=MAPE-sd, ymax=MAPE+sd), colour="black", width=.1) +
    geom_line(color="blue",size=1.2) +
    geom_point(size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Algorithm") + ylab("MAE") + 
    theme_bw() +
    theme(title=element_text(face="bold",size=16),      
          axis.title.y = element_text(face="bold", size=16),
          axis.text = element_text(colour="#000000", size=16),
          legend.position="none"); 

  p3 = ggplot(df.RMSE.err, aes(x=Algorithm, y=RMSE, colour="purple",group=1)) + 
    geom_errorbar(aes(ymin=RMSE-sd, ymax=RMSE+sd), colour="black", width=.1) +
    geom_line(color="Purple",size=1.2) +
    geom_point(size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Algorithm") + ylab("MAE") + 
    theme_bw() +
    theme(title=element_text(face="bold",size=16),      
          axis.title.y = element_text(face="bold", size=16),
          axis.text = element_text(colour="#000000", size=16),
          legend.position="none"); 
  
  grid.arrange(p1,p2,p3,ncol=3);
}
  
#load('Results\\NewResult_CrossValidation_WestSpeedway_WB2.RData')
load('Results\\NewResult_CrossValidation_MidSpeedway_EB.RData')
#load('Results\\NewResult_CrossValidation_SimuData.RData')

tmp = AggregateData(results)
MAE = tmp$MAE
MAPE = tmp$MAPE
RMSE = tmp$RMSE

#plotME(MAE,MAPE,RMSE)

MAPE.Benchmark = tmp$MAPE.Benchmark
MAPE.UniGMM = tmp$MAPE.UniGMM
MAPE.HMM = tmp$MAPE.HMM

MAPE.Benchmark.all = data.frame()
MAPE.UniGMM.all = data.frame()
MAPE.HMM.all = data.frame()

for(run in 1:10)
{
#   East Speedway
  colnames(MAPE.Benchmark[[run]]) = c("L4","L5","L6","L7","L8","L9")
  colnames(MAPE.UniGMM[[run]]) = c("L4","L5","L6","L7","L8","L9")
  colnames(MAPE.HMM[[run]]) = c("L4","L5","L6","L7","L8","L9")

#   West Speedway  
#   colnames(MAPE.Benchmark[[run]]) = c("L3","L2","L1")
#   colnames(MAPE.UniGMM[[run]]) = c("L3","L2","L1")
#   colnames(MAPE.HMM[[run]]) = c("L3","L2","L1")
  
  MAPE.Benchmark[[run]]$run = run
  MAPE.UniGMM[[run]]$run = run
  MAPE.HMM[[run]]$run = run
    
  MAPE.Benchmark.all = rbind(MAPE.Benchmark.all,MAPE.Benchmark[[run]])
  MAPE.UniGMM.all = rbind(MAPE.UniGMM.all,MAPE.UniGMM[[run]])
  MAPE.HMM.all = rbind(MAPE.HMM.all,MAPE.HMM[[run]])  
}

MAPE.Benchmark.all$method = "Benchmark"
MAPE.UniGMM.all$method = "GMM"
MAPE.HMM.all$method = "HMM"

MAPE.cmb = rbind(MAPE.Benchmark.all,MAPE.UniGMM.all,MAPE.HMM.all)


# p1 = ggplot(MAPE.cmb,aes(x=method,y=L1)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L1") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));   
# 
# p2 = ggplot(MAPE.cmb,aes(x=method,y=L2)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of L2") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));   
# 
# p3 = ggplot(MAPE.cmb,aes(x=method,y=L3)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L3") + ylim(0,4) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16)); 
# grid.arrange(p1,p2,p3,ncol=3);

# East Speedway
# p1 = ggplot(MAPE.cmb,aes(x=method,y=L4)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L4") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));   
# 
# p2 = ggplot(MAPE.cmb,aes(x=method,y=L5)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of L5") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));   
# 
# p3 = ggplot(MAPE.cmb,aes(x=method,y=L6)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L6") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));  
# 
# p4 = ggplot(MAPE.cmb,aes(x=method,y=L7)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L7") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));  
# 
# p5 = ggplot(MAPE.cmb,aes(x=method,y=L8)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L8") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));  
# 
# p6 = ggplot(MAPE.cmb,aes(x=method,y=L9)) + theme_bw() +
#   geom_boxplot()+ ylab("MAPE of  L9") + ylim(0,1) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16));  

# grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3);

# Density plot
# West Speedway
# MAPE.cmb$sum = MAPE.cmb$L1+MAPE.cmb$L2+MAPE.cmb$L3
# MAPE.cmb$method = factor(MAPE.cmb$method,levels=c('Benchmark','GMM','HMM'));
# ggplot(MAPE.cmb,aes(sum,fill=method)) + geom_density(alpha=0.5) + theme_bw() + 
#   scale_fill_manual(values = c("gray", "red","green")) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16),
#         legend.justification=c(1,1),legend.position=c(1,1)) +   
#   labs(x="Total MAPE Each Observation",y="Density");

# East Speedway
# MAPE.cmb$sum = MAPE.cmb$L4+MAPE.cmb$L5+MAPE.cmb$L6 + MAPE.cmb$L7 + MAPE.cmb$L8 + MAPE.cmb$L9
# MAPE.cmb$method = factor(MAPE.cmb$method,levels=c('Benchmark','GMM','HMM'));
# ggplot(MAPE.cmb,aes(sum,fill=method)) + geom_density(alpha=0.5) + theme_bw() + 
#   scale_fill_manual(values = c("gray", "red","green")) +
#   theme(title=element_text(face="bold",size=16),      
#         axis.title.y = element_text(face="bold", size=16),
#         axis.text = element_text(colour="#000000", size=16),
#         legend.title = element_text(size=16),
#         legend.text = element_text(size = 16),
#         legend.justification=c(1,1),legend.position=c(1,1)) +   
#   labs(x="Total MAPE Each Observation",y="Density");
