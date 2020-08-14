#Differential Gene Expression using DESeq2
#Author: Anurag Kanase
# Version: 1
# Github: dezember
# Date: August 12, 2020
# Project: Precision Drug Prediction

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install("DESeq2")
library("DESeq2")

dat<-read.csv("Features_1.csv",header = T,quote = "",row.names = 1,sep = ",")

# Convert to matrix
dat <- as.matrix(dat)
head(dat)

# Assign condition (first three are WT, next three are mutants)

condition <- factor(c(rep("Healthy",4),rep("Infected",4)))
condition=relevel(condition,ref = "Healthy")


# Create a coldata frame: its rows correspond to columns of dat (i.e., matrix representing the countData)
coldata <- data.frame(row.names=colnames(dat), condition)

head(coldata)

#     condition

#H1   Healthy
#H2   Healthy
#H3   Healthy
#H4   Healthy
#I1  Infected
#I2  Infected

##### DESEq pipeline, first the design and the next step, normalizing to model fitting
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)

dds <- DESeq(dds)

# Plot Dispersions:
png("qc-dis_1.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")

dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Principal Components Analysis
plotPCA(rld)
dev.off()

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diff_exp_results.csv",quote = FALSE,row.names = F)

