
library(umap)
library(Seurat)
library(gplots)
library(limma)
library(ReactomePA)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reactome.db)
library(SingleR)
library(ggplot2)
library(ggfortify)
library(survival)
library(tidyverse)
library(tableone)
library(beeswarm)
library(pROC)


##############################################
a_name = "13_rib_min_scaled_withGN_unscaled_PCAs_QN"

setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#exclude possible spike-ins and blacklist genes
raw.genecounts=raw.genecounts[grep("ERCC",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("UGDH.AS1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM2P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("LOC100131257",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("KCNQ1OT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MALAT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MAB21L3",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("EEF1A1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]


#Scale before quantile normalization
raw.genecounts=t(t(raw.genecounts)/colSums(raw.genecounts))*100000

#quantile normalize
raw.genecounts=round(limma::normalizeQuantiles(raw.genecounts))

#create seurat object;  my first dataset I ever analyzed was from colon... I keep habbit of naming the single cell object this way ;)
colon <- CreateSeuratObject(raw.genecounts, min.cells = 1,min.features = 1,
                            project = "AE")


#select fixed number of variable genes
colon<- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(colon), 35)

#plot variable genes
plot1 <- VariableFeaturePlot(colon)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scale data based on all genes
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = c(all.genes))


#runPCAS

colon <- RunPCA(colon, features = VariableFeatures(object = colon))
VizDimLoadings(colon, dims = 1:2, reduction = "pca")
DimPlot(colon, reduction = "pca")
DimHeatmap(colon, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(colon, dims = 1:9, cells = 500, balanced = TRUE)
colon <- JackStraw(colon, num.replicate = 100)
colon <- ScoreJackStraw(colon, dims = 1:20)
JackStrawPlot(colon, dims = 1:20)
ElbowPlot(colon)


colon <- NormalizeData(object = colon, normalization.method = "LogNormalize", 
                         scale.factor = 10000)

#find clusters and make tSNE

colon <- FindNeighbors(colon, dims = 1:12)
colon <- FindClusters(colon, resolution = 0.5)
head(Idents(colon), 10)
colon <- RunTSNE(object = colon, dims.use = 1:12, do.fast = TRUE)
colon <- RunUMAP(object = colon, dims = 1:12)

#DimPlot(colon)
DimPlot(colon, reduction = "umap")
DimPlot(colon, reduction = "pca")
DimPlot(colon, reduction = "tsne")
TSNEPlot(object = colon)



#find marker genes
colon.markers <- FindAllMarkers(object = colon, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)



write.table(colon.markers,file = paste(a_name,paste("all","_DEGs.txt")),sep = "\t",quote = F)

