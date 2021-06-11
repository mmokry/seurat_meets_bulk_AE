
library(umap)
library(Seurat)
library(limma)



##############################################
setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

#read in read count matrix files
raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log; there are 4096 possible UMIs
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#check
raw.genecounts[1:10,1:10]

#exclude spike-ins
raw.genecounts=raw.genecounts[grep("ERCC",row.names(raw.genecounts),invert=TRUE),]
#exclude unvanted genes - historically those genes have mapping issues
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

#create seurat object
colon <- CreateSeuratObject(raw.genecounts, min.cells = 1,min.features = 1,
                            project = "AE")

#add correlation coefficients into the object 
cor.matrix = cor(log2(raw.genecounts+10))
colon[["mean.cor"]] <- colMeans(cor.matrix)
aa <- cor.matrix[order(row(cor.matrix), -cor.matrix)]
second.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 2])
tenth.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 10])
names(second.best)=row.names(cor.matrix)
colon[["second.best.cor"]] <- second.best
names(tenth.best)=row.names(cor.matrix)
colon[["tenth.best.cor"]] <- tenth.best



#get percentatge of mito genes
colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
percent.mt=PercentageFeatureSet(colon, pattern = "^MT-")
VlnPlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(colon, features = c("second.best.cor", "tenth.best.cor", "mean.cor"), ncol = 3)

colon<- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 5000)
top <- head(VariableFeatures(colon), 35)

#plot variable genes
plot1 <- VariableFeaturePlot(colon)
plot1
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
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


colon.markers <- FindAllMarkers(object = colon, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)


