library( "DESeq2" )
library(ggplot2)

countData <- read.csv('counts.csv', header = TRUE, sep = ",")
head(countData)
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
metaData
dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~dex, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
#ENTER GENE OF INTEREST
par(mfrow=c(2,3))
####change it .. just for demo
plotCounts(dds, gene="XX", intgroup="dex")
plotCounts(dds, gene="XX", intgroup="dex")
plotCounts(dds, gene="XX", intgroup="dex")
plotCounts(dds, gene="XX", intgroup="dex")
plotCounts(dds, gene="XX", intgroup="dex")
plotCounts(dds, gene="XX", intgroup="dex")
#Next steps in exploring these data...BLAST to database to find associated gene function

#####Log fold change shrinkage
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")


#Volcano Plot

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA

#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup="dex") #using the DESEQ2 plotPCA fxn we can

#MA-plot
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


write.csv(as.data.frame(res),file="condition_treated_results.csv")
