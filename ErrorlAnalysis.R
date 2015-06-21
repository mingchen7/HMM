library(ggplot2)
library(gridExtra)
load('ValidationResults.RData');

cat('MAE of Benchmark',sum(MAE.benchmark));
cat('MAE of UniGMM',sum(MAE.UniGMM));
cat('MAE of MultiGMM',sum(MAE.MultiGMM));
cat('MAE of HMM',sum(MAE.HMM));

cat('MAPE of Benchmark',sum(MAPE.benchmark));
cat('MAPE of UniGMM',sum(MAPE.UniGMM));
cat('MAPE of MultiGMM',sum(MAPE.MultiGMM));
cat('MAPE of HMM',sum(MAPE.HMM));

colSums(MAPE.benchmark);
colSums(MAPE.UniGMM);
colSums(MAPE.MultiGMM);
colSums(MAPE.HMM);

# MAE
# MAE.benchmark$sum = MAE.benchmark$tt2 + MAE.benchmark$tt3 + MAE.benchmark$tt4 + MAE.benchmark$tt5 +MAE.benchmark$tt6;
# MAE.UniGMM$sum = MAE.UniGMM$tt2 + MAE.UniGMM$tt3 + MAE.UniGMM$tt4 + MAE.UniGMM$tt5 +MAE.UniGMM$tt6;
# MAE.MultiGMM$sum = MAE.MultiGMM$tt2 + MAE.MultiGMM$tt3 + MAE.MultiGMM$tt4 + MAE.MultiGMM$tt5 +MAE.MultiGMM$tt6;
# MAE.HMM$sum = MAE.HMM$tt2 + MAE.HMM$tt3 + MAE.HMM$tt4 + MAE.HMM$tt5 +MAE.HMM$tt6;

MAE.benchmark$sum = MAE.benchmark[,1] + MAE.benchmark[,2] + MAE.benchmark[,3];
MAE.UniGMM$sum = MAE.UniGMM[,1] + MAE.UniGMM[,2] + MAE.UniGMM[,3];
MAE.MultiGMM$sum = MAE.MultiGMM[,1] + MAE.MultiGMM[,2] + MAE.MultiGMM[,3];
MAE.HMM$sum = MAE.HMM[,1] + MAE.HMM[,2] + MAE.HMM[,3];

MAE.benchmark$method = 'Benchmark';
MAE.UniGMM$method = 'Uni. GMM';
MAE.MultiGMM$method = 'Multi. GMM';
MAE.HMM$method = 'HMM';
MAE.all = rbind(MAE.benchmark,MAE.UniGMM,MAE.MultiGMM,MAE.HMM);

# MAPE
# MAPE.benchmark$sum = MAPE.benchmark$tt2 + MAPE.benchmark$tt3 + MAPE.benchmark$tt4 + MAPE.benchmark$tt5 +MAPE.benchmark$tt6;
# MAPE.UniGMM$sum = MAPE.UniGMM$tt2 + MAPE.UniGMM$tt3 + MAPE.UniGMM$tt4 + MAPE.UniGMM$tt5 +MAPE.UniGMM$tt6;
# MAPE.MultiGMM$sum = MAPE.MultiGMM$tt2 + MAPE.MultiGMM$tt3 + MAPE.MultiGMM$tt4 + MAPE.MultiGMM$tt5 +MAPE.MultiGMM$tt6;
# MAPE.HMM$sum = MAPE.HMM$tt2 + MAPE.HMM$tt3 + MAPE.HMM$tt4 + MAPE.HMM$tt5 +MAPE.HMM$tt6;

MAPE.benchmark$sum = MAPE.benchmark[,1] + MAPE.benchmark[,2] + MAPE.benchmark[,3];
MAPE.UniGMM$sum = MAPE.UniGMM[,1] + MAPE.UniGMM[,2] + MAPE.UniGMM[,3];
MAPE.MultiGMM$sum = MAPE.MultiGMM[,1] + MAPE.MultiGMM[,2] + MAPE.MultiGMM[,3];
MAPE.HMM$sum = MAPE.HMM[,1] + MAPE.HMM[,2] + MAPE.HMM[,3];

MAPE.benchmark$method = 'Benchmark';
MAPE.UniGMM$method = 'Uni. GMM';
MAPE.MultiGMM$method = 'Multi. GMM' 
MAPE.HMM$method = 'HMM';
MAPE.all = rbind(MAPE.benchmark,MAPE.UniGMM,MAPE.MultiGMM,MAPE.HMM);


# Histogram of MAE and MAPE
MAE.all$method = factor(MAE.all$method,levels=c('Benchmark','Uni. GMM','Multi. GMM','HMM'));
h1 = ggplot(MAE.all,aes(sum,fill=method)) + geom_density(alpha=0.5) + theme_bw() + 
     scale_fill_manual(values = c("gray", "blue","red","green")) +
     theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16),
        legend.justification=c(1,1),legend.position=c(1,1)) + 
     xlim(0,800) + 
     labs(x="Sum of MAEs",y="Density");

MAPE.all$method = factor(MAPE.all$method,levels=c('Benchmark','Uni. GMM','Multi. GMM','HMM'));
h2 =ggplot(MAPE.all,aes(sum,fill=method)) + geom_density(alpha=0.5) + theme_bw() + 
    scale_fill_manual(values = c("gray", "blue","red","green")) +
    theme(title=element_text(face="bold",size=16),      
          axis.title.y = element_text(face="bold", size=16),
          axis.text = element_text(colour="#000000", size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size = 16),
          legend.justification=c(1,1),legend.position=c(1,1)) + 
    xlim(0,10) +
    labs(x="Sum of MAPEs",y="Density");

grid.arrange(h1,h2,ncol=2);

# Comparison by variable
p1 = ggplot(MAPE.all,aes(x=method,y=tt2)) + theme_bw() +
  geom_boxplot()+ ylab("MAPE for tt2") + 
  theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));    

p2 = ggplot(MAPE.all,aes(x=method,y=tt3)) + theme_bw() +
  geom_boxplot()+ ylab("MAPE for tt3") + 
  theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));                                                                                                      

p3 = ggplot(MAPE.all,aes(x=method,y=tt4)) + theme_bw() +
  geom_boxplot() +  ylab("MAPE for tt4") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));

p4 = ggplot(MAPE.all,aes(x=method,y=tt5)) + theme_bw() +
  geom_boxplot() + ylab("MAPE for tt5") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));

p5 = ggplot(MAPE.all,aes(x=method,y=tt6)) + theme_bw() +
  geom_boxplot() + ylab("MAPE for tt6") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));
                                                                                                
grid.arrange(p1, p2, p3, p4, p5, ncol=3)


# For Bluetooth data
p1 = ggplot(MAPE.all,aes(x=method,y=Camp2Mnt)) + theme_bw() +
  geom_boxplot()+ ylab("MAPE for Camp2Mnt") + 
  theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));    

p2 = ggplot(MAPE.all,aes(x=method,y=Mnt2Park)) + theme_bw() +
  geom_boxplot()+ ylab("MAPE for Mnt2Park") + 
  theme(title=element_text(face="bold",size=16),      
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));                                                                                                      

p3 = ggplot(MAPE.all,aes(x=method,y=Park2Euclid)) + theme_bw() +
  geom_boxplot() +  ylab("MAPE for Park2Euclid") +
  theme(title=element_text(face="bold",size=16),        
        axis.title.y = element_text(face="bold", size=16),
        axis.text = element_text(colour="#000000", size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16));

grid.arrange(p1, p2, p3,ncol=3)

