library(tximportData)
library(tximport)
library("DESeq2")
library("vsn")
library("RColorBrewer")


# setting up environment
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis")

setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_6")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_17")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_24")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_39")

setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_23")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_29")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_chloroplast")
setwd("/Users/guritsk1/Desktop/RNA_SEQ_2017/DESeq2_analysis/MAG_ANALYSIS/BIN_algae")



samples <- read.table("samples.txt", header = TRUE)
names(files) <- paste0("sample", 1:6)
files <- file.path("salmon_out", samples$Dir, "quant.sf")
all(file.exists(files))


# loading Salmon data and pre-processing
txi <- tximport(files, type = "salmon", txOut = TRUE)
names(txi)
head(txi$counts)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Time)
keep <- rowSums(counts(ddsTxi)) >= 1
dds <- ddsTxi[keep,]


# DESeq analysis
dds <- DESeq(dds)
res <- results(dds, alpha=0.05) #(FDR=5%)
write.csv(as.data.frame(resOrdered), file="significance_results.csv")
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.01, na.rm=TRUE)

# transform data
vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)


# sample comparison plots
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Time, vsd$Day, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

plotMA(res, ylim=c(-2,2))
heatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
plotCounts(dds, gene=which.min(res$padj), intgroup="Time") # plot expresison of most differential gene
plotPCA(vsd, intgroup="Time")
sum(res$padj < 0.01, na.rm=TRUE)







