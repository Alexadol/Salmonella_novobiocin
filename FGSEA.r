#packages loading
library(DESeq2)
library(fgsea)
library(ggplot2)
library(apeglm) #It's needed for lfcshrink

#set working directory
setwd('/home/alexadol/TKA3')

#open and bind tables with raw counts 
res_with <- read.csv('fc_res_with.csv',sep=' ') #with antibiotic
res_wo <- read.csv('fc_res_wo.csv',sep=' ') #without
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
resultsNames(dds)

#Getting shrunk lfc (as example when we want to see difference in 60 minutes when antibiotic( treatment) is added )
resLFC <- lfcShrink(dds, coef = "treattreated.time60", type="apeglm")

#Creating vector with ranked genes (Shrunked LFC is used as rank)
resLFC$name <- rownames(resLFC)
ranks_lfc <- as.data.frame(cbind(resLFC$name,resLFC$log2FoldChange))
ranks_lfc <- as.data.frame(ranks_lfc[complete.cases(ranks_lfc),])
colnames(ranks_lfc) <- c('name','lfc')
ranks_lfc$lfc <- as.numeric(as.character(ranks_lfc$lfc))
ranks_lfc_vec <- setNames(ranks_lfc$lfc, ranks_lfc$name) 


#Open file with custom gene sets in .gmt format
pathways <- gmtPathways('geneset_conv_div_without_frames_with_stream.gmt') 

#fgsea Analysis implementation
fgseaRes <- fgsea(pathways=pathways, stats=ranks_lfc_vec, nperm=1000)

#Barplot for results visualization
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() + scale_fill_manual(values = c('#02818a','#ae017e'))+
  labs(x="Gene set", y="Normalized Enrichment Score",
       title="Gene sets enriched in 20 minutes while treated (GSEA analysis)") + 
  theme_classic(base_size = 15)

#Default Enrichment plot
plotEnrichment(pathways[["diverging_genes_forward"]],ranks_lfc_vec)  # я в основном вот на этот график смотрю, ну и на padj в табличке fgseaRes

#Function for making custom version of Enrichment plot
plotEnr = function (pathway, stats, gseaParam = 1, ticksSize = 0.2)
{   
  LINECOL = "#fe9929"  #Here we can change line color
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "grey",
               linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms),
               colour = "grey", linetype = "dashed") + geom_hline(yintercept = 0,
                                                                  colour = "grey") + geom_line(col=LINECOL) + theme_classic(base_size = 15) +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x,
                                                               y = -diff/2, xend = x, yend = diff/2), size = ticksSize, colour = '#636363') +#here we can change color and tickness of bars which visualize genes from gene set 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "Rank", y = "Enrichment score") + ggtitle(label = 'Converging genes reverse (<)')#Here we can change color
  g
}

#Example of usage 
plotEnr(pathways[["diverging_genes_forward"]],ranks_lfc_vec) 
