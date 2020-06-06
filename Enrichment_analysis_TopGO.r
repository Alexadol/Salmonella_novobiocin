#set working directory
setwd('/home/alexadol/TKA3')

#packages loading
library(DESeq2)
library(topGO)
library(ggplot2)
library(ALL)
data("geneList") #For topDiffGene function


#open and bind tables with raw counts 
res_with <- read.csv('counts_res_with.csv',sep=' ') #with antibiotic
res_wo <- read.csv('counts_res_without.csv',sep=' ') #without
res_full <- cbind(res_wo,res_with)

#adding FAKE 3 columns with zero (as zero treatment) for making groups equal
fake_zero_treat <- res_full[,c(1,2,3)]
colnames(fake_zero_treat)<- c('fakewith_0_1','fakewith_0_2','fakewith_0_3')
res_full_new <- cbind(fake_zero_treat,res_full)

#Creating of colData object, 3 factors (time 0/10/20/60, treat treated/untreated and concentration 0/100/500)
time <- c(rep(0,6),rep(10,3),rep(20,3),rep(60,6),rep(10,3),rep(20,3),rep(60,3))
treat<-c(rep('treated',3),rep('untreated',12),rep('treated',12))
concentration <- c(rep(0,15),rep(100,3),rep(500,9))
full_coldata <- as.data.frame(cbind(time,treat,concentration))
rownames(full_coldata)<-colnames(res_full_new)

#creating DESEqDataSet object 
dds <- DESeqDataSetFromMatrix(res_full_new,colData = full_coldata,design=~treat*time)
#set reference for treat (untreated)
dds$treat = relevel(dds$treat,"untreated")

#DE analysis
dds <- DESeq(dds)
#Get the names of all possible comparisons (since we have complex design)
resultsNames(dds)

#results for difference in expression between 60 and zero minutes when antibiotic added
res_60_ab <- results(dds, contrast=list(c('time_60_vs_0','treattreated.time60')))
#results for difference in expression between 20 and zero minutes when antibiotic added
res_20_ab <- results(dds, contrast=list(c('time_20_vs_0','treattreated.time20')))
#results for difference in expression between 10 and zero minutes when antibiotic added
res_10_ab <- results(dds, contrast=list(c('time_10_vs_0','treattreated.time10')))

#exclude NA results
res_wo_NA_60 <- res_60_ab[!is.na(res_60_ab$log2FoldChange) &!is.na(res_60_ab$padj) ,]
res_wo_NA_20 <- res_20_ab[!is.na(res_20_ab$log2FoldChange) &!is.na(res_20_ab$padj),]
res_wo_NA_10 <- res_10_ab[!is.na(res_10_ab$log2FoldChange) &!is.na(res_10_ab$padj),]

#select UP and DOWn regulated genes
Up_regulated_60 <- res_wo_NA_60[res_wo_NA_60$log2FoldChange > 1 & res_wo_NA_60$padj < 0.05,]
Down_regulated_60 <- res_wo_NA_60[res_wo_NA_60$log2FoldChange < -1 & res_wo_NA_60$padj < 0.05,]

Up_regulated_20 <- res_wo_NA_20[res_wo_NA_20$log2FoldChange > 1 & res_wo_NA_20$padj < 0.05,]
Down_regulated_20 <- res_wo_NA_20[res_wo_NA_20$log2FoldChange < -1 & res_wo_NA_20$padj < 0.05,]

Up_regulated_10 <- res_wo_NA_10[res_wo_NA_10$log2FoldChange > 1 & res_wo_NA_10$padj < 0.05,]
Down_regulated_10 <- res_wo_NA_10[res_wo_NA_10$log2FoldChange < -1 & res_wo_NA_10$padj < 0.05,]

#open table with gene2GO mapping
map_paste_x <- readMappings('/home/alexadol/TKA3/Salmonella_enterica_14028S_gene_to_GO')

#create named vectors with p-values
#60 minutes
test_num_60_up <- Up_regulated_60$pvalue
names(test_num_60_up)<-rownames(Up_regulated_60)
test_num_60_down <- Down_regulated_60$pvalue
names(test_num_60_down)<-rownames(Down_regulated_60) 
#20 minutes
test_num_20_up <- Up_regulated_20$pvalue
names(test_num_20_up)<-rownames(Up_regulated_20)
test_num_20_down <- Down_regulated_20$pvalue
names(test_num_20_down)<-rownames(Down_regulated_20) 
#10 minutes
test_num_10_up <- Up_regulated_10$pvalue
names(test_num_10_up)<-rownames(Up_regulated_10)
test_num_10_down <- Down_regulated_10$pvalue
names(test_num_10_down)<-rownames(Down_regulated_10) 

#Create topGO object for analysis using named vector with p-values, topDiffGene function (select the differentially expressed genes, at 0.01 significance level) and object with gene2GO mapping
sample_wo_GOdata_60_up <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = test_num_60_up, geneSel = topDiffGenes,nodeSize = 5, annot = annFUN.gene2GO, gene2GO = map_paste_x)
sample_wo_GOdata_60_down <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = test_num_60_down, geneSel = topDiffGenes,nodeSize = 5, annot = annFUN.gene2GO, gene2GO = map_paste_x)

resultKS_60_up <- runTest(sample_wo_GOdata_60_up, algorithm = "classic", statistic = "ks")
resultKS_60_down <- runTest(sample_wo_GOdata_60_down, algorithm = "classic", statistic = "ks")

allGO_60_up=usedGO(sample_wo_GOdata_60_up)
allGO_60_down=usedGO(sample_wo_GOdata_60_down)


sample_wo_GOdata_20_up <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = test_num_20_up, geneSel = topDiffGenes,nodeSize = 1, annot = annFUN.gene2GO, gene2GO = map_paste_x)
sample_wo_GOdata_20_down <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = test_num_20_down, geneSel = topDiffGenes,nodeSize = 1, annot = annFUN.gene2GO, gene2GO = map_paste_x)

resultKS_20_up <- runTest(sample_wo_GOdata_20_up, algorithm = "classic", statistic = "ks")
resultKS_20_down <- runTest(sample_wo_GOdata_20_down, algorithm = "classic", statistic = "ks")

allGO_20_up=usedGO(sample_wo_GOdata_20_up)
allGO_20_down=usedGO(sample_wo_GOdata_20_down)

#HERE WE SHOWED EXAMPLE OF ENRICHMENT ANALYSIS FOR 20 AND 60 MINUTES WHEN ANTIBIOTIC IS ADDED BUT THIS CODE CAN BE EXTENDED FOR OTHER CASES

#create tables with te results of analysis
res_60_up = GenTable(sample_wo_GOdata_60_up, weightFisher=resultKS_60_up, topNodes=length(allGO_60_up)) 
#Take first 20 GO with the highest enrichment score
res_sign_60_up <- res_60_up[1:20,]

res_60_down = GenTable(sample_wo_GOdata_60_down, weightFisher=resultKS_60_down, topNodes=length(allGO_60_down))
res_sign_60_down <- res_60_down[1:20,]

res_20_up = GenTable(sample_wo_GOdata_20_up, weightFisher=resultKS_20_up, topNodes=length(allGO_20_up)) 
res_sign_20_up <- res_20_up[1:20,]

res_20_down = GenTable(sample_wo_GOdata_20_down, weightFisher=resultKS_20_down, topNodes=length(allGO_20_down))
res_sign_20_down <- res_20_down[1:20,]


#create bar plot for enrichment results(one example)
res_sign_20_up$Term <- factor(res_sign_20_up$Term, levels=rev(res_sign_20_up$Term))
ggplot(res_sign_20_up, aes(x=Term, y=-log10(as.numeric(weightFisher)))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value), KS") +
  ggtitle("Upregulated in treatment 20 min vs 0 min")+ coord_flip()+scale_y_continuous(breaks = round(seq(0, max(-log10(as.numeric(res_sign$weightFisher))), by = 2), 1))+
  geom_col(aes(fill = as.numeric(res_sign$Significant)))+
  scale_fill_gradient(aes(colour = "Mapped genes"),high = "#c51b8a",low='cornsilk')+ #COLOR CAN BE CHANGED HERE
  theme(
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=20, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=15,  hjust=1.10),
    axis.text.y=element_text(angle=0, size=14, vjust=0.7),
    axis.title=element_text(size=20, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(0.8, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=14),
    panel.background = element_blank(),#Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5)))+
  geom_hline(yintercept=2, linetype="dashed", 
             color = "gray28", size=1)
                     
#The following code was used to performe a likelihood ratio test, where we remove the treatment-specific differences over time. Genes with small p values from this test are those which, at one or more time points after time 0 showed a treatment-specific effect 
library(DESeq2)
ddsTC <- DESeqDataSetFromMatrix(res_full_new,colData = full_coldata,design=~treat+time+treat:time)
ddsTC$treat = relevel(ddsTC$treat,"untreated")
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ treat + time)
resTC <- results(ddsTC)
resultsNames(ddsTC)
resTC <- resTC[!is.na(resTC$padj),]

test_num_TC <- resTC$pvalue
names(test_num_TC)<-rownames(resTC)
 

sample_wo_GOdata_TC <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = test_num_TC, geneSel = topDiffGenes,nodeSize = 1, annot = annFUN.gene2GO, gene2GO = map_paste_x)

resultKS_TC <- runTest(sample_wo_GOdata_TC, algorithm = "classic", statistic = "ks")

allGO_TC=usedGO(sample_wo_GOdata_TC)

res_TC = GenTable(sample_wo_GOdata_TC, weightFisher=resultKS_TC, topNodes=length(allGO_TC))
#Take first 20 GO with the highest enrichment score
res_sign_TC <- res_TC[1:20,]

res_sign_TC$Term <- factor(res_sign_TC$Term, levels=rev(res_sign_TC$Term))
ggplot(res_sign_TC, aes(x=Term, y=-log10(as.numeric(weightFisher)))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value), KS") +
  ggtitle("Upregulated in treatment 20 min vs 0 min")+ coord_flip()+scale_y_continuous(breaks = round(seq(0, max(-log10(as.numeric(res_sign_TC$weightFisher))), by = 2), 1))+
  geom_col(aes(fill = as.numeric(res_sign_TC$Significant)))+
  scale_fill_gradient(aes(colour = "Mapped genes"),high = "blue",low='cornsilk')+
  theme(
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=20, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=15,  hjust=1.10),
    axis.text.y=element_text(angle=0, size=14, vjust=0.7),
    axis.title=element_text(size=20, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(0.8, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=14),
    panel.background = element_blank(),#Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5)))+
  geom_hline(yintercept=2, linetype="dashed", 
             color = "gray28", size=1)

#Extract gene names from enriched GO
AnnotatedGenes = lapply(res_sign_TC$GO.ID, function(x) as.character(unlist(genesInTerm(object = sample_wo_GOdata_TC, whichGO = x))))
names(AnnotatedGenes) <- res_sign_TC$Term
AnnotatedGenes$`branched-chain amino acid metabolic proc...`


#SOMETIMES IT IS ALSO VERY USEFUL TO  ASSESS AN EXPRESSION DINAMIC BY RAW COUNTS 
#HERE is a short code for taking raw counts data for one gene and for visualization of these data
ssrA<- plotCounts(ddsTC,'STM14_1083', 
                   intgroup = c("time","treat"), returnData = TRUE)
ssrA$time <- as.numeric(as.character(ssrA$time))

p <- ggplot(ssrA,
       aes(x = time, y = count, color = treat, group = treat)) + scale_color_manual(values=c('#2ca25f','#ef6548'))+
  geom_point() + stat_summary(fun.y=mean, geom="line")+ggtitle('SsrA')+theme_classic(base_size = 15)

p



