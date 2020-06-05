#Loading required packages
library(edgeR)

#Set working directory
setwd('/home/alexadol/TKA3/')

#Open tables with raw counts
res_with <- read.csv('counts_res_with.csv',sep=' ') #with antibiotic
res_wo <- read.csv('counts_res_without.csv',sep=' ') #without antibiotic

#binding two tables in one
res_full <- cbind(res_wo,res_with)
#making matrix from data frame
my_y <- as.matrix((res_full))

#Creating DGEList object (EdgeR), setting group names for samples
#wo - without antibiotic with - with antibiotic 0,10,20,60 - time after treatment
my_y <- DGEList(counts = my_y, group=c(rep('wo_0min',3),rep('wo_10min',3),rep('wo_20min',3),rep('wo_60min',3),rep('with_100_60min',3),rep('with_500_10min',3),rep('with_500_20min',3),rep('with_500_60min',3)))
#Normalization
my_y <- calcNormFactors(my_y)
#Getting counts per million instead of raw counts
my_z <- cpm(my_y, normalized.lib.size=TRUE)

#transpose table
scaledata <- t(my_y)
#select rows without NA
scaledata <- scaledata[complete.cases(scaledata),]
scaledata <- as.matrix(scaledata)
#Get mean of 3 samples in all condition (we have 3 repeats for each condition)
scaledata_mean <- data.frame(ID=my_y[,0], wo_0min_mean=rowMeans(my_y[,1:3]),wo_10min_mean=rowMeans(my_y[,4:6]), wo_20min_mean=rowMeans(my_y[,7:9]),wo_60min_mean=rowMeans(my_y[,10:12]),with_100_60min_mean=rowMeans(my_y[,13:15]),with_500_10_mean=rowMeans(my_y[,16:18]),with_500_20min_mean=rowMeans(my_y[,19:21]),with_500_60min_mean=rowMeans(my_y[,22:24]))


#Identify optimal number of clusters (when K-means is used for clusterization) using Gap-statistic
library(cluster)
set.seed(13)
gap <- clusGap(scaledata_mean, kmeans, 20, B = 100, verbose = interactive())
#Visualize results
plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)

##Identify optimal number of clusters using Within groups sum of squares
wss <- (nrow(scaledata_mean)-1)*sum(apply(scaledata_mean,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata_mean,
                                     centers=i)$withinss)
#Visualize results
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#Implement clusterization using kmeans function (In our case optimal numbers of clusters was determined as 4)
set.seed(20)
kClust <- kmeans(scaledata_mean, centers=4, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
#Prepare results for visualization
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata_mean, kClusters)
kClustcentroids_with <- kClustcentroids[c(1,6,7,8,5),]

#Change names of points in which genes within clusters have similar expression dinamics
rownames(kClustcentroids_with) <- c ('0 minutes','10 minutes novobiocin 500 mkg','20 minutes novobiocin 500 mkg','60 minutes novobiocin 500 mkg','60 minutes novobiocin 100 mkg')

library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids_with)
colnames(Kmolten) <- c('sample','cluster','value')

#Creating plot to assess dynamics of genes within clusters
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time with novobiocin treatment",color = "Cluster")+
  scale_color_manual(values=c('#993404','#E69F00', '#56B4E9','#31a354'))+
  theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1

#Checking if clusters have low correlation level (Less number of clusters should be used if some of them have high correlation values)
cor(kClustcentroids)
#Extract genes from clusters
K1 <- (scaledata_mean[kClusters==4,])
#Get score of genes in cluster (how far it from core, values near 1 identify genes which have dinamics near to core of cluster)
score <- apply(K1, 1, corscore)
score_names <- names(score)
score_names <- as.data.frame(score_names)
#save genes from cluster to table
write.table(score_names,'cluster_up_in_60_wo_ab_names.txt')
