#BiocManager::install(c("ReactomePA","biomaRt","clusterProfiler"))

library(umap)
library(Seurat)
library(gplots)
library(limma)
library(ReactomePA)
library(dplyr)
library(biomaRt)#
library(clusterProfiler)
library(org.Hs.eg.db)
library(reactome.db)
library(SingleR)


#reticulate::py_install(packages ='umap-learn')
##############################################

copycat = 1 # copycat object on and off
logNorm_later = 1 # performs log norm after PCA calculatiobs

#stamp for filenames from this analysis
a_name = "11_ribo_min_genes_scaled_withGN_unscaled_PCAs_QN"
setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis2")

raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

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
colon <- CreateSeuratObject(raw.genecounts, min.cells = 20,min.features = 9000,
                            project = "AE")
if (copycat == 1){
  raw.genecountsENS=raw.genecounts
  namesENS=sapply( strsplit( row.names(raw.genecountsENS), "_ENSG"), "[", 2)
  namesENS=paste("ENSG",namesENS,sep = "")
  row.names(raw.genecountsENS)=namesENS
  
  colonENS <- CreateSeuratObject(raw.genecountsENS, min.cells = 20,min.features = 9000,
                              project = "AE")
  }

#add correlation coefficients onto the object for future excludion of possible duplicates and samples that are large outliers
cor.matrix = cor(log2(raw.genecounts+10))
colon[["mean.cor"]] <- colMeans(cor.matrix)
aa <- cor.matrix[order(row(cor.matrix), -cor.matrix)]
second.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 2])
tenth.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 10])
names(second.best)=row.names(cor.matrix)
colon[["second.best.cor"]] <- second.best
names(tenth.best)=row.names(cor.matrix)
colon[["tenth.best.cor"]] <- tenth.best

#add metadata
f=read.table(file = "clinical_data_good_selection.txt",header = T,sep = "\t",row.names = 1)
for (i in names(f)){
  feature = f[,i]
  names(feature)=row.names(f)
  colon[[i]] <- feature
}


#get percentatge of mito genes
colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
percent.mt=PercentageFeatureSet(colon, pattern = "^MT-")
pdf(file = paste(a_name,"_mt_percentage_ngenes_coverage.pdf"),height = 4)
VlnPlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(colon, features = c("second.best.cor", "tenth.best.cor", "mean.cor"), ncol = 3)



dev.off()

#normalize the counts
if (logNorm_later == 0) {
  colon
  colon <- NormalizeData(object = colon, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
  if (copycat == 1){
    colonENS <- NormalizeData(object = colonENS, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  }
}
#select fixed number of variable genes
colon<- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 5000)
if (copycat == 1){
  colonENS<- FindVariableFeatures(colonENS, selection.method = "vst", nfeatures = 5000)
  }

top10 <- head(VariableFeatures(colon), 35)

#plot variable genes
pdf(file = paste(a_name,"_variable genes.pdf"),height = 7,width = 7)
plot1 <- VariableFeaturePlot(colon)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

#scale data based on all genes
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = c(all.genes))

if (copycat ==1){
  all.genesENS <- rownames(colonENS)
  colonENS <- ScaleData(colonENS, features = c(all.genesENS))
  
}

#regress out batch - this does not work properly
#f=read.table(file = "metadata_good.txt",header = T,sep = "\t",row.names = 1)
#colon[["seq.batch"]] <- f
#colon <- ScaleData(colon, vars.to.regress = c("nFeature_RNA","seq.batch"))


#runPCAS
pdf(file = paste(a_name,"_PCA_QC.pdf"),height = 5)
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

if (copycat ==1){
  colonENS <- RunPCA(colonENS, features = VariableFeatures(object = colonENS))

}
if (logNorm_later == 1) {
  colon
  colon <- NormalizeData(object = colon, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  if (copycat == 1){
    colonENS <- NormalizeData(object = colonENS, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
  }
}


#find clusters and make tSNE
pdf(file = paste(a_name,"_tSNE_QC.pdf"),height = 4,width = 4)
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
dev.off()

if (copycat==1){
  colonENS <- FindNeighbors(colonENS, dims = 1:12)
  colonENS <- FindClusters(colonENS, resolution = 0.5)
  colonENS <- RunTSNE(object = colonENS, dims.use = 1:12, do.fast = TRUE)
  colonENS <- RunUMAP(object = colonENS, dims = 1:12)
  
  
}



#find marker genes
colon.markers <- FindAllMarkers(object = colon, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)

c02.markers <- FindMarkers(colon, ident.1 = "0", ident.2 = "2")

write.table(colon.markers,file = paste(a_name,paste("all","_DEGs.txt")),sep = "\t",quote = F)
write.table(c02.markers,file = paste(a_name,paste("c0-c2","_DEGs.txt")),sep = "\t",quote = F)
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
#  dev.off()
#}
  ###################### compare pathway scores between modules
  if (length(as.data.frame(PA)[,1])>1){
    pdf(file = paste(a_name,paste(i,"_Volcano.pdf")),height = 4, width = 4)
    for (j in 1: length(as.data.frame(PA)[,1])){
      pathway=c(as.data.frame(PA)[j,1])
      path_genes <- getBM(
      filters="reactome",
      attributes=c( "ensembl_gene_id"),
      values=pathway,
      mart=mart
      )
      if (length(path_genes[,1])){
      colonENS= AddModuleScore(colonENS, 
                           features = path_genes, 
                           pool = NULL, 
                           nbin = 24, 
                           ctrl = 100,
                           k = FALSE, 
                           assay = NULL, name = as.data.frame(PA)[j,2], seed = 1)
      print(VlnPlot(colonENS, features = c(paste(as.data.frame(PA)[j,2],1,sep = "")), ncol = 3))
      print (paste(j,length(as.data.frame(PA)[,1]),sep = " pathway projection out of "))
      }
 
    }
 


  dev.off()
  }
  
  ############custom pathwaysprojection
#  xx <- as.list(reactomePATHID2EXTID)
#  PAS = c("R-HSA-71403")
  
#  for (j in 1: length(PAS)){
#  print(PAS[j])
    
#    pathway=c(PAS[j])
#    path_genes <- getBM(
#      filters="reactome",
#      attributes=c( "ensembl_gene_id"),
#      values=pathway,
#      mart=mart
#    )
#    if (length(path_genes[,1])){
#      colonENS= AddModuleScore(colonENS, 
##                               features = path_genes, 
#                               pool = NULL, 
#                               nbin = 24, 
#                               ctrl = 100,
#                               k = FALSE, 
#                               assay = NULL, name = as.data.frame(PA)[j,2], seed = 1)
#      print(VlnPlot(colonENS, features = "The citric acid (TCA) cycle and respiratory electron transport1" , ncol = 3))
#      print (paste(j,length(as.data.frame(PA)[,1]),sep = " pathway projection out of "))
#    }
    
#  }
  
  ########################
  
  
  
  
  ####################
  
  
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



genesv=c(
  "ACAT1-ENSG00000075239",
 "ACADSB-ENSG00000196177",
  "KYNU-ENSG00000115919",
  "LPIN1-ENSG00000134324",
  "FGF13-ENSG00000129682",
 
  "ACTA2-ENSG00000107796",
  "MYH11-ENSG00000276480",
  "MYH10-ENSG00000133026",
 
 "CXCL12-ENSG00000107562",
 "NOS1-ENSG00000089250",
 "nFeature_RNA", "nCount_RNA", "percent.mt",
 "second.best.cor", "tenth.best.cor", "mean.cor",
 "VDR-ENSG00000111424",
 "SOD2-ENSG00000112096",
 

  "CD14-ENSG00000170458",
  "C1QA-ENSG00000173372",
  "CD63-ENSG00000135404"  ,   
   "CD74-ENSG00000019582" ,
 "APOE-ENSG00000130203",
 "CD4-ENSG000000106101"
 
# "KLF4-ENSG00000136826"
 
  
  )



for (i in 1:length(genesv)){
print(genesv[i])
pdf(file = paste(paste(a_name,genesv[i]),"_VLN_plots.pdf"),height = 2.5,width = 3.5)
print(VlnPlot(colon,x.lab.rot=T, features = c(genesv[i]),pt.size = 0.7))
dev.off()
}


for (i in 1:length(genesv)){
  print(genesv[i])
  pdf(file = paste(paste(a_name,genesv[i]),"_fature_plots.pdf"),height = 3,width = 3)
  print(FeaturePlot(object = colon, features = c(genesv[i]),cols = topo.colors(100)))
  dev.off()
}




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
VlnPlot(colon,x.lab.rot=T, features = c("MKI67-ENSG00000148773"))
VlnPlot(colon,x.lab.rot=T, features = c("ACTB-ENSG00000075624"))
VlnPlot(colon,x.lab.rot=T, features = c("CXCL12-ENSG00000107562"))
VlnPlot(colon,x.lab.rot=T, features = c("STARD9-ENSG00000159433"))
VlnPlot(colon,x.lab.rot=T, features = c("MTATP8P2-ENSG00000229604"))
VlnPlot(colon,x.lab.rot=T, features = c("CD74-ENSG00000019582"))
VlnPlot(colon,x.lab.rot=T, features = c("CYSLTR1-ENSG00000173198"))
VlnPlot(colon,x.lab.rot=T, features = c("IFI30-ENSG00000216490"))
VlnPlot(colon,x.lab.rot=T, features = c("CSF3R-ENSG00000119535"))
VlnPlot(colon,x.lab.rot=T, features = c("CYSLTR1-ENSG00000173198"))
VlnPlot(colon,x.lab.rot=T, features = c("KIT-ENSG00000157404"))
VlnPlot(colon,x.lab.rot=T, features = c("CCN5-ENSG00000064205"))   #CCN5
VlnPlot(colonENS,x.lab.rot=T, features = c("ENSG00000064205"))   #PDGFRB
VlnPlot(colon,x.lab.rot=T, features = c("CD4-ENSG00000010610")) 
VlnPlot(colon,x.lab.rot=T, features = c("VDR-ENSG00000111424"))
VlnPlot(colon,x.lab.rot=T, features = c("CD3E-ENSG00000198851"))

CD4-ENSG000000106101



dev.off()


pdf(file = paste(a_name,"_Feature_plots_genes.pdf"),height = 3.5,width = 3.5)
FeaturePlot(object = colon, features = "FGF13-ENSG00000129682",cols = topo.colors(100))
FeaturePlot(object = colon, features = "ACTA2-ENSG00000107796",cols = topo.colors(100))
FeaturePlot(object = colon, features = "HBB-ENSG00000244734",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MYH11-ENSG00000276480",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CXCL12-ENSG00000107562",cols = topo.colors(100))
FeaturePlot(object = colon, features = "ACKR1-ENSG00000213088",cols = topo.colors(100))
FeaturePlot(object = colon, features = "KLF4-ENSG00000136826",cols = topo.colors(100))
FeaturePlot(object = colon, features = "APOE-ENSG00000130203",cols = topo.colors(100))
FeaturePlot(object = colon, features = "XIST-ENSG00000229807",cols = topo.colors(100))
FeaturePlot(object = colon, features = "ACTB-ENSG00000075624",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MYH11-ENSG00000276480",cols = topo.colors(100))
FeaturePlot(object = colon, features = "STARD9-ENSG00000159433",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MTATP8P2-ENSG00000229604",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CD74-ENSG00000019582",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CYSLTR1-ENSG00000173198",cols = topo.colors(100))
FeaturePlot(object = colon, features = "IFI30-ENSG00000216490",cols = topo.colors(100))

FeaturePlot(object = colon, features = "CSF3R-ENSG00000119535",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CYSLTR1-ENSG00000173198",cols = topo.colors(100))
FeaturePlot(object = colon, features = "KIT-ENSG00000157404",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CCN5-ENSG00000064205",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MKI67-ENSG00000148773",cols = topo.colors(100))
FeaturePlot(object = colonENS, features = "ENSG00000113721",cols = topo.colors(100))
FeaturePlot(object = colon, features = "VDR-ENSG00000111424",cols = topo.colors(100))



dev.off()


library("dplyr")
top10 <- colon.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

pdf(file = paste(a_name,"_heatmap_top10markers.pdf"),height = 3,width = 8)
DoHeatmap(colon, features = genesv) + NoLegend()
dev.off()
pdf(file = paste(a_name,"_heatmap_fewselectedmarkers.pdf"),height = 10,width = 10)
DoHeatmap(colon, features = top10$gene) + NoLegend()
dev.off()
pdf(file = paste(a_name,"_heatmap_topselectedmarkers.pdf"),height = 10,width = 10)
DoHeatmap(colon, features = c("CCDC144A-ENSG00000170160","CYSLTR1-ENSG00000173198",
                              "KYNU-ENSG00000115919", "ADAM32-ENSG00000275594",  
                              "CYP46A1-ENSG00000036530",   "TMEM212-ENSG00000186329",
                              "SOD2-ENSG00000112096",   
                              
                              "PCSK5-ENSG00000099139", "CALD1-ENSG00000122786",   
                              "MYH11-ENSG00000276480"   ,  "DSTN-ENSG00000125868" ,
                              "ACTA2-ENSG00000107796"  ,    "LMOD1-ENSG00000163431" ,  
                              "MYH10-ENSG00000133026"   , 
                                
                              "IGFBP7-ENSG00000163453"   , "USP8-ENSG00000138592"  ,    
                              "ZNF286A-ENSG00000187607"  ,
                               "ATXN3-ENSG00000066427"   ,  "RPS6KA5-ENSG00000100784" ,  "SLC35E3-ENSG00000175782",  
                              "ORAI2-ENSG00000160991"  ,       
                              
                               "CD81-ENSG00000110651"   ,   "C1QA-ENSG00000173372"    ,  "CD63-ENSG00000135404"  ,   
                               "C1QB-ENSG00000173369"   ,   "PLTP-ENSG00000100979"    ,  "GRN-ENSG00000030582"   ,   
                               "SRGN-ENSG00000122862"   ,   "TYROBP-ENSG00000011600"  ,  "CD74-ENSG00000019582"  ,   
                               "SAT1-ENSG00000130066"      , 
                              
                              "CTSB-ENSG00000164733"   ,  
                              "IFI30-ENSG00000216490"  ,   "FTH1-ENSG00000167996"    ,  "PSAP-ENSG00000197746"    , 
                              "FTL-ENSG00000087086"    ,   "GPNMB-ENSG00000136235"  ,   "TYROBP-ENSG00000011600"  , 
                         
                              "APOE-ENSG00000130203"   ,   "APOC1-ENSG00000130208"  ,  
                              "SRGN-ENSG00000122862"   ,   "PKM-ENSG00000067225" ,     
                              "IGLL5-ENSG00000254709"    )) + NoLegend()
dev.off()

#correlate_data_with_single_cell_deconvolution
#pdf(file = paste(a_name,"_single_cell.pdf"),height = 4,width = 4)
f=read.table(file = "MuSiC.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:7)){

  #  if (kruskal.test(f[names(Idents(colon)),i]~as.numeric(Idents(colon)))$p.value<0.1){
  pdf(file = paste(paste(a_name,names(f)[i]),"_single_cell.pdf"),height = 3.5,width = 3.5)
  par(mar=c(5,8,4,3))
  boxplot(f[names(Idents(colon)),i]~as.numeric(Idents(colon)),
          col= c(1:6),
          main =  paste(names(f)[i],paste(" p = ",kruskal.test(f[names(Idents(colon)),i]~as.numeric(Idents(colon)))$p.value)),
  )
  dev.off()
  # }
}





#correlate results with clinical data
pdf(file = paste(a_name,"_clinical_data.pdf"),height = 5,width = 5)
f=read.table(file = "clinical_data_good_selection.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:18)){
  tbl = table(f[names(Idents(colon)),i],as.numeric(Idents(colon)))

  
  tb=t(tbl)
  for (j in colnames(tb)){
    feature = tb[,j]
    nonfeature = rowSums(tb)-tb[,j]
    tbs = cbind(feature,nonfeature)
    chisq.test(cbind(feature,nonfeature))
    barplot((tbs[,"feature"]/(tbs[,"nonfeature"]+tbs[,"feature"])-mean((tbs[,"feature"]/(tbs[,"nonfeature"]+tbs[,"feature"]))))*100,
          col = c(1:10),ylab = "percentage change",
          main =  paste(names(f)[i],paste(j,paste(" p = ",chisq.test(tbs)$p.value)))
          )
    barplot((tbs[,"feature"]/(tbs[,"nonfeature"]+tbs[,"feature"])*100),
           col = c(1:10),ylab = "percentage",
          main =  paste(names(f)[i],paste(j,paste(" p = ",chisq.test(tbs)$p.value)))
          )
  }
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
   
   par(mar=c(5,8,4,3))
   tbl=t(tbl)
   barplot(scale(tbl,center = FALSE,
                 scale = colSums(tbl)),
           col = c(1:length(tbl[,1])),
           main =  paste(names(f)[i],paste(" p = ",chisq.test(tbl)$p.value)),
           #legend = row.names(tbl)
   )
   legend(1,1,legend=row.names(tbl), col = c(1:length(tbl[,1])),pch = 15)
#  }
}


for (i in c(3,6,19:52)){
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


f=read.table(file = "metadata_batch.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1:3)){
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




f=read.table(file = "metadata_check.txt",header = T,sep = "\t",row.names = 1)

The citric acid (TCA) cycle and respiratory electron transport1

wilcox.test( list(stroke=colonENS@meta.data$"The citric acid (TCA) cycle and respiratory electron transport1"[f$symptoms_inclusion=="stroke"],
              TIA=colonENS@meta.data$"The citric acid (TCA) cycle and respiratory electron transport1"[f$symptoms_inclusion=="TIA"],
              asymptomatic=colonENS@meta.data$"The citric acid (TCA) cycle and respiratory electron transport1"[f$symptoms_inclusion=="ocular"],
              ocular=colonENS@meta.data$"The citric acid (TCA) cycle and respiratory electron transport1"[f$symptoms_inclusion=="asymptomatic"]
))

col = c(1:6))
boxplot( list(stroke=colonENS@meta.data$"Neutrophil degranulation1"[f$symptoms_inclusion=="stroke"],
              TIA=colonENS@meta.data$"Neutrophil degranulation1"[f$symptoms_inclusion=="TIA"],
              asymptomatic=colonENS@meta.data$"Neutrophil degranulation1"[f$symptoms_inclusion=="ocular"],
              ocular=colonENS@meta.data$"Neutrophil degranulation1"[f$symptoms_inclusion=="asymptomatic"]
              ),
         
         col = c(1:6))





## single R
singler <- CreateSinglerObject(counts = seuset@raw.data[,rownames(seuset@meta.data)],       
                               annot = seuset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "Single Cell Plaques", 
                               min.genes = 0,                          # if starting from seurat, set min.genes to 0 
                               technology = "CEL-seq2", 
                               species = "Human", 
                               citation = "",
                               ref.list = list(), 
                               normalize.gene.length = FALSE, 
                               variable.genes = "de",
                               fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
                               do.signatures = TRUE,
                               clusters = seuset@ident,
                               do.main.types = TRUE, 
                               reduce.file.size = TRUE, 
                               numCores = 4
)
pdf(paste(result.folder, "/", Sys.Date(), " SingleR heatmap 16 communities2.pdf", sep = ""))
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)
dev.off()


#####survival
library(ggfortify)
library(survival)


f=read.table(file = "ep_major.txt",header = T,sep = "\t",row.names = 1)
surv=data.frame(Idents(colon))
survv=f[names(Idents(colon)),c(1:2)]
combined=cbind(surv,survv)
combined$ep.major=as.numeric(combined$ep_major)

combined$ep.major[combined$ep.major == 2] = 5
combined$ep.major[combined$ep.major == 1] = 2
combined$ep.major[combined$ep.major == 5] = 1
fit = survfit(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined)
autoplot(fit, conf.int = F)
combined$Idents.colon.[combined$Idents.colon. == 4] = 3
fit = survfit(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined)
autoplot(fit, conf.int = F)

combined$Idents.colon.[combined$Idents.colon. == 4] = 3
combined$Idents.colon.[combined$Idents.colon. == 2] = 1
combined$Idents.colon.[combined$Idents.colon. == 0] = 1


fit = survfit(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined)
autoplot(fit, conf.int = T)






f=read.table(file = "ep_composite.txt",header = T,sep = "\t",row.names = 1)
surv=data.frame(Idents(colon))
survv=f[names(Idents(colon)),c(1:2)]
combined=cbind(surv,survv)
combined$ep.composite=as.numeric(combined$ep_composite)

combined$ep.composite[combined$ep.composite == 2] = 5
combined$ep.composite[combined$ep.composite == 1] = 2
combined$ep.composite[combined$ep.composite == 5] = 1
fit = survfit(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined)
autoplot(fit, conf.int = F)
combined$Idents.colon.[combined$Idents.colon. == 4] = 3
fit = survfit(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined)
autoplot(fit, conf.int = T)


combined$Idents.colon.[combined$Idents.colon. == 4] = 3
combined$Idents.colon.[combined$Idents.colon. == 2] = 1
combined$Idents.colon.[combined$Idents.colon. == 0] = 1


fit = survfit(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_composite_time,ep.composite)~Idents.colon.,data = combined)
autoplot(fit, conf.int = T)



