#BiocManager::install(c("ReactomePA","biomaRt","clusterProfiler"))


library(Seurat)
library(gplots)
library(limma)
library(ReactomePA)
library(dplyr)
library(biomaRt)
library(clusterProfiler)

#a_name = "02_POLII_genes_scaled_QN"
a_name = "02_all_genes_scaled_QN"
#setwd("D:/R/atheroexpress")
setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress")
raw.genecounts = read.table (file = "raw_counts.txtt",header = T,sep = "\t",row.names = 1)

#correct for UMI sampling
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#check header
raw.genecounts[1:10,1:10]

#exclude spike-ins
raw.genecounts=raw.genecounts[grep("ERCC",row.names(raw.genecounts),invert=TRUE),]
#exclude unvanted genes
#raw.genecounts=raw.genecounts[grep("^MT-",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("^RPL",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("^RPS",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("^RP",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("UGDH.AS1",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("PGM2P2",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("LOC100131257",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("KCNQ1OT1",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("MALAT1",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("MAB21L3",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("EEF1A1",row.names(raw.genecounts),invert=TRUE),]
#raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]

#Scale before quantile normalization
raw.genecounts=t(t(raw.genecounts)/colSums(raw.genecounts))*100000

#quantile normalize
raw.genecounts=round(limma::normalizeQuantiles(raw.genecounts))

#create seurat object
colon <- CreateSeuratObject(raw.genecounts, min.cells = 20,min.features = 10000,
                            project = "AE")

#add correlation coefficients onto the object for future excludion of possible duplicates and samples that are large outliers
#cor.matrix = cor(raw.genecounts)
cor.matrix = cor(log2(raw.genecounts+10))
colon[["mean.cor"]] <- colMeans(cor.matrix)
aa <- cor.matrix[order(row(cor.matrix), -cor.matrix)]
second.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 2])
tenth.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 10])
names(second.best)=row.names(cor.matrix)
colon[["second.best.cor"]] <- second.best
names(tenth.best)=row.names(cor.matrix)
colon[["tenth.best.cor"]] <- tenth.best

#pdf(file = "figures_seurat_bulk_add_more_pheno_v3.4_all_genes_scale_before_QN.pdf")
#pdf(file = "figures_seurat_bulk_add_more_pheno_v3.44_all_genes_scale_before_QN.pdf")

#get percentatge of mito genes
colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
percent.mt=PercentageFeatureSet(colon, pattern = "^MT-")
pdf(file = paste(a_name,"_mt_percentage_ngenes_coverage.pdf"),height = 4)
VlnPlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#normalize the counts
colon
colon <- NormalizeData(object = colon, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

#select fixed number of variable genes
colon<- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(colon), 15)

#plot variable genes
pdf(file = paste(a_name,"_variable genes.pdf"),height = 4,width = 7)
plot1 <- VariableFeaturePlot(colon)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

#scale data based on all genes
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = c(all.genes))

#regress out batch - this does not work properly
#f=read.table(file = "metadata_good.txt",header = T,sep = "\t",row.names = 1)
#colon[["seq.batch"]] <- f
#colon <- ScaleData(colon, vars.to.regress = c("nFeature_RNA","seq.batch"))


#runPCAS
pdf(file = paste(a_name,"_PCA_QC.pdf"),height = 4)
colon <- RunPCA(colon, features = VariableFeatures(object = colon))
VizDimLoadings(colon, dims = 1:2, reduction = "pca")
DimPlot(colon, reduction = "pca")
DimHeatmap(colon, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(colon, dims = 1:9, cells = 500, balanced = TRUE)
colon <- JackStraw(colon, num.replicate = 100)
colon <- ScoreJackStraw(colon, dims = 1:20)
JackStrawPlot(colon, dims = 1:20)
ElbowPlot(colon)
dev.off()


#find clusters and make tSNE
pdf(file = paste(a_name,"_tSNE_QC.pdf"),height = 4,width = 4)
colon <- FindNeighbors(colon, dims = 1:16)
colon <- FindClusters(colon, resolution = 0.9)
head(Idents(colon), 10)
colon <- RunTSNE(object = colon, dims.use = 1:16, do.fast = TRUE)
#colon <- RunUMAP(colon, dims = 1:16)
#DimPlot(colon)
DimPlot(colon, reduction = "pca")
TSNEPlot(object = colon)
dev.off()

#find marker genes
colon.markers <- FindAllMarkers(object = colon, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)


####pathway analysis
DEGs_entrez_all=as.character(0)
DEGs_entrez_all=as.list(DEGs_entrez_all)
for (i in c(0:( max(as.numeric(colon$seurat_clusters))-1))){
  print(i)
#extract genes
  DEGs=sapply( strsplit( (colon.markers)[colon.markers$cluster==i,7], "-ENSG"), "[", 2)
  DEGs=paste("ENSG",DEGs,sep = "")
  if (length(DEGs)>500){
    DEGs = DEGs[1:500]
  }

  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  DEGs_entrez <- getBM(
     filters="ensembl_gene_id",
     attributes=c( "entrezgene_id"),
     values=DEGs,
     mart=mart
  )
  DEGs_entrez_all[[as.character(i)]]=as.vector(DEGs_entrez[,1])
  PA <- enrichPathway(gene=as.vector(DEGs_entrez[,1]),pvalueCutoff=0.05, readable=T)
  if (length(as.data.frame(PA)[,1])>1){
    pdf(file = paste(a_name,paste(i,"_Reactome.pdf")),height = 7, width = 13)
    print(barplot(PA, showCategory=10))
    print(barplot(PA, showCategory=15))
    print(barplot(PA, showCategory=20))
    print(barplot(PA, showCategory=50))
    print(dotplot(PA, showCategory=10))
    print(dotplot(PA, showCategory=15))
    print(dotplot(PA, showCategory=20))
    print(dotplot(PA, showCategory=50))
    print(emapplot(PA, showCategory=10))
    print(emapplot(PA, showCategory=15))
    print(emapplot(PA, showCategory=20))
    print(emapplot(PA, showCategory=50))
    print(emapplot(PA, showCategory=200))
    print(cnetplot(PA, categorySize="pvalue",showCategory=8))
    print(cnetplot(PA, categorySize="pvalue",showCategory=20))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_barplot10.pdf")),height = 3, width = 13)
    print(barplot(PA, showCategory=10))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_barplot15.pdf")),height = 4, width = 13)
    print(barplot(PA, showCategory=15))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_barplot20.pdf")),height = 5, width = 13)
    print(barplot(PA, showCategory=20))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_barplot50.pdf")),height = 9, width = 13)
    print(barplot(PA, showCategory=50))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_dotplot10.pdf")),height = 3.5, width = 13)
    print(dotplot(PA, showCategory=10))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_dotplot15.pdf")),height = 4, width = 13)
    print(dotplot(PA, showCategory=15))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_dotplot20.pdf")),height = 5, width = 13)
    print(dotplot(PA, showCategory=20))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_dotplot50.pdf")),height = 9, width = 13)
    print(dotplot(PA, showCategory=50))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot10.pdf")),height = 5, width = 5)
    print(emapplot(PA, showCategory=10))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot15.pdf")),height = 6, width = 6)
    print(emapplot(PA, showCategory=15))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot20.pdf")),height = 7, width = 7)
    print(emapplot(PA, showCategory=20))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot50.pdf")),height = 10, width = 10)
    print(emapplot(PA, showCategory=50))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot100.pdf")),height = 15, width = 15)
    print(emapplot(PA, showCategory=100))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_emapplot200.pdf")),height = 20, width = 20)
    print(emapplot(PA, showCategory=200))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_cneplot10.pdf")),height = 10, width = 10)
    print(cnetplot(PA, categorySize="pvalue",showCategory=10))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_cneplot20.pdf")),height = 14, width = 14)
    print(cnetplot(PA, categorySize="pvalue",showCategory=20))
    dev.off()
    
    pdf(file = paste(a_name,paste(i,"_Reactome_cneplot40.pdf")),height = 18, width = 18)
    print(cnetplot(PA, categorySize="pvalue",showCategory=40))
    dev.off()
    
  }
  write.table(as.data.frame(PA),file = paste(a_name,paste(i,"_Reactome.txt")),sep = "\t",quote = F)
}

#compare_pathways_clusters
cRes <- compareCluster(DEGs_entrez_all[-1], fun="enrichPathway")
pdf(file = paste(a_name,"_COmpare_clusters.pdf"),height = 12,width = 12)
dotplot(cRes,showCategory=30)
dev.off()



pdf(file = paste(a_name,"_VLN_plots_QC.pdf"),height = 3.5,width = 7)
VlnPlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(colon, features = c("second.best.cor", "tenth.best.cor", "mean.cor"), ncol = 3, same.y.lims = TRUE)


dev.off()
pdf(file = paste(a_name,"_VLN_plots_genes.pdf"),height = 3.5,width = 3)
VlnPlot(colon,x.lab.rot=T, features = c("ACTA2-ENSG00000107796"))
VlnPlot(colon,x.lab.rot=T, features = c("MYH11-ENSG00000276480"))
VlnPlot(colon,x.lab.rot=T, features = c("FGF13-ENSG00000129682"))
VlnPlot(colon,x.lab.rot=T, features = c("ACKR1-ENSG00000213088"))
VlnPlot(colon,x.lab.rot=T, features = c("KLF4-ENSG00000136826"))
VlnPlot(colon,x.lab.rot=T, features = c("APOE-ENSG00000130203"))
VlnPlot(colon,x.lab.rot=T, features = c("HBB-ENSG00000244734"))
VlnPlot(colon,x.lab.rot=T, features = c("XIST-ENSG00000229807"))
VlnPlot(colon,x.lab.rot=T, features = c("CXCL12-ENSG00000107562"))
VlnPlot(colon,x.lab.rot=T, features = c("FGF13-ENSG00000129682"))
VlnPlot(colon,x.lab.rot=T, features = c("ACTB-ENSG00000075624"))
VlnPlot(colon,x.lab.rot=T, features = c("CXCL12-ENSG00000107562"))
VlnPlot(colon,x.lab.rot=T, features = c("STARD9-ENSG00000159433"))
VlnPlot(colon,x.lab.rot=T, features = c("MTATP8P2???ENSG00000229604"))
VlnPlot(colon,x.lab.rot=T, features = c("CD74???ENSG00000019582"))
VlnPlot(colon,x.lab.rot=T, features = c("CYSLTR1???ENSG00000173198"))
VlnPlot(colon,x.lab.rot=T, features = c("IFI30???ENSG00000216490"))

dev.off()




pdf(file = paste(a_name,"_Feature_plots_genes.pdf"),height = 3.5,width = 3.5)
FeaturePlot(object = colon, features = "FGF13-ENSG00000129682",cols = rainbow(100))
FeaturePlot(object = colon, features = "ACTA2-ENSG00000107796",cols = rainbow(100))
FeaturePlot(object = colon, features = "HBB-ENSG00000244734",cols = rainbow(100))
FeaturePlot(object = colon, features = "MYH11-ENSG00000276480",cols = rainbow(100))
FeaturePlot(object = colon, features = "CXCL12-ENSG00000107562",cols = rainbow(100))
FeaturePlot(object = colon, features = "ACKR1-ENSG00000213088",cols = rainbow(100))
FeaturePlot(object = colon, features = "KLF4-ENSG00000136826",cols = rainbow(100))
FeaturePlot(object = colon, features = "APOE-ENSG00000130203",cols = rainbow(100))
FeaturePlot(object = colon, features = "XIST-ENSG00000229807",cols = rainbow(100))
FeaturePlot(object = colon, features = "ACTB-ENSG00000075624",cols = rainbow(100))
FeaturePlot(object = colon, features = "MYH11-ENSG00000276480",cols = rainbow(100))
FeaturePlot(object = colon, features = "STARD9-ENSG00000159433",cols = rainbow(100))
FeaturePlot(object = colon, features = "MTATP8P2???ENSG00000229604",cols = rainbow(100))
FeaturePlot(object = colon, features = "CD74???ENSG00000019582",cols = rainbow(100))
FeaturePlot(object = colon, features = "CYSLTR1???ENSG00000173198",cols = rainbow(100))
FeaturePlot(object = colon, features = "IFI30???ENSG00000216490",cols = rainbow(100))

dev.off()


library("dplyr")
top10 <- colon.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf(file = paste(a_name,"_heatmap_top10markers.pdf"),height = 10,width = 10)
DoHeatmap(colon, features = top10$gene) + NoLegend()
dev.off()


#correlate results with clinical data
pdf(file = paste(a_name,"_clinical_data.pdf"),height = 5,width = 5)
f=read.table(file = "clinical_data_good_reduced.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:25)){
  tbl = table(f[names(Idents(colon)),i],as.numeric(Idents(colon)))
#  if (chisq.test(tbl)$p.value<0.1){
    par(mar=c(5,8,4,3))
    barplot(tbl,
        col= c(1:length(tbl[,1])),
        main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value))
        )
  
    par(mar=c(5,8,4,3))
    barplot(scale(tbl,center = FALSE,
        scale = colSums(tbl)),
        col = c(1:length(tbl[,1])),
        main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value))
        )

    par(mar=c(5,8,4,3))
    barplot(scale(tbl,center = FALSE,
        scale = colSums(tbl)),
        col = c(1:length(tbl[,1])),
        main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
        #legend = row.names(tbl)
        )
   legend(1,1,legend=row.names(tbl), col = c(1:length(tbl[,1])),pch = 15)
#  }
}


for (i in c(14:22,24)){
#  if (kruskal.test(f[names(Idents(colon)),i]~as.numeric(Idents(colon)))$p.value<0.1){
  par(mar=c(5,8,4,3))
  boxplot(f[names(Idents(colon)),i]~as.numeric(Idents(colon)),
        col= c(1:6),
        main =  paste(names(f)[i],paste(" p = ",kruskal.test(f[names(Idents(colon)),i]~as.numeric(Idents(colon)))$p.value)),
    )
 # }
}
dev.off()



pdf(file = paste(a_name,"_outcomes_covs.pdf"),height = 5,width = 5)
f=read.table(file = "metadata_check.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:7,10:13)){
#  tbl = table(f[names(Idents(colon)),i],sample(as.numeric(Idents(colon))))
  
  tbl = table(f[names(Idents(colon)),i],as.numeric(Idents(colon)))
  
  par(mar=c(5,8,4,3))
  barplot(tbl,
          col= c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
          
  )
  
  par(mar=c(5,8,4,3))
  barplot(scale(tbl,center = FALSE,
                scale = colSums(tbl)),
          col = c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
          #  legend = row.names(tbl)
  )
  
  
  par(mar=c(5,8,4,3))
  barplot(scale(tbl,center = FALSE,
                scale = colSums(tbl)),
          col = c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
        #  legend = row.names(tbl)
  )
  legend(1,1,legend=row.names(tbl), col = c(1:length(tbl[,1])),pch = 15)
}


f=read.table(file = "metadata_good.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:5)){
  #  tbl = table(f[names(Idents(colon)),i],sample(as.numeric(Idents(colon))))
  
  tbl = table(f[names(Idents(colon)),i],as.numeric(Idents(colon)))
  
  par(mar=c(5,8,4,3))
  barplot(tbl,
          col= c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
          
  )
  
  par(mar=c(5,8,4,3))
  barplot(scale(tbl,center = FALSE,
                scale = colSums(tbl)),
          col = c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
          #  legend = row.names(tbl)
  )
  
  
  par(mar=c(5,8,4,3))
  barplot(scale(tbl,center = FALSE,
                scale = colSums(tbl)),
          col = c(1:length(tbl[,1])),
          main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
          #  legend = row.names(tbl)
  )
  legend(1,1,legend=row.names(tbl), col = c(1:length(tbl[,1])),pch = 15)
}


dev.off()





