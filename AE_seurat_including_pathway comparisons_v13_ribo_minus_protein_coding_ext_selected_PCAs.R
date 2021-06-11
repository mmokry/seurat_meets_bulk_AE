#BiocManager::install(c("ReactomePA","biomaRt","clusterProfiler"))

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

#reticulate::py_install(packages ='umap-learn')
##############################################

copycat = 1 # copycat object on and off - it creates seurat objects with same content but different gene names (ensembl, HGCN), as it is difficult to change it in existing seurat object
logNorm_later = 1 # performs log norm after PCA calculations

#stamp for filenames from this analysis
a_name = "14_as13_PCAs_model_QN"

setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log
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
colon <- CreateSeuratObject(raw.genecounts, min.cells = 1,min.features = 1,
                            project = "AE")
if (copycat == 1){
  raw.genecountsENS=raw.genecounts
  namesENS=sapply( strsplit( row.names(raw.genecountsENS), "_ENSG"), "[", 2)
  namesENS=paste("ENSG",namesENS,sep = "")
  row.names(raw.genecountsENS)=namesENS
  
  colonENS <- CreateSeuratObject(raw.genecountsENS, min.cells = 1,min.features = 1,
                              project = "AE")
  raw.genecountsHG=raw.genecounts
  namesHG=sapply( strsplit( row.names(raw.genecountsHG), "_ENSG"), "[", 1)
  row.names(raw.genecountsHG)=namesHG
  
  colonHG <- CreateSeuratObject(raw.genecountsHG, min.cells = 1,min.features = 1,
                                 project = "AE")
  
  }

#add correlation coefficients onto the object for future exclusion of possible duplicates and samples that are large outliers
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

f=read.table(file = "MuSiC.txt",header = T,sep = "\t",row.names = 1)
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
    if (copycat == 1){
      colonHG <- NormalizeData(object = colonHG, normalization.method = "LogNormalize", 
                                scale.factor = 10000)
      
  }
}
#select fixed number of variable genes
colon<- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 5000)
if (copycat == 1){
  colonENS<- FindVariableFeatures(colonENS, selection.method = "vst", nfeatures = 5000)
  colonHG<- FindVariableFeatures(colonHG, selection.method = "vst", nfeatures = 5000)
  }

top10 <- head(VariableFeatures(colon), 35)

#plot variable genes
pdf(file = paste(a_name,"_variable genes.pdf"),height = 5,width = 5)
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
  all.genesHG <- rownames(colonHG)
  colonHG <- ScaleData(colonHG, features = c(all.genesHG))
  
}

#regress out batch - this does not work properly and its probably not needed
#f=read.table(file = "metadata_good.txt",header = T,sep = "\t",row.names = 1)
#colon[["seq.batch"]] <- f
#colon <- ScaleData(colon, vars.to.regress = c("seq.batch"))


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
  colonHG <- RunPCA(colonHG, features = VariableFeatures(object = colonHG))
}
if (logNorm_later == 1) {
  colon
  colon <- NormalizeData(object = colon, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  if (copycat == 1){
    colonENS <- NormalizeData(object = colonENS, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
    colonHG <- NormalizeData(object = colonHG, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
    
  }
}

##############   asymptomatic / TIA+ocular ####
#RNA_PCA = "symptoms_2g ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_6 + PC_7 + PC_18 + PC_20 + PC_22 + PC_24 + PC_26 + PC_28 + PC_31 + PC_34 + PC_45"

#all_parameters = "symptoms_2g ~ macmean0 + calcification + smc + macrophages_location +     smc_macrophages_ratio + PC_1 + PC_2 + PC_3 + PC_4 + PC_6 +     PC_7 + PC_8 + PC_18 + PC_20 + PC_22 + PC_23 + PC_24 + PC_26 +     PC_27 + PC_28 + PC_31 + PC_37 + PC_40 + PC_42 + PC_45"

##############   logistic model no TIA ####
###################################
"symptoms_2g ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_7 + PC_11"

#find clusters and make tSNE
pdf(file = paste(a_name,"_tSNE_QC.pdf"),height = 4,width = 4)
colon <- FindNeighbors(colon, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))
colon <- FindClusters(colon, resolution = 0.5)
head(Idents(colon), 10)
colon <- RunTSNE(object = colon, dims.use = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45), do.fast = TRUE)
colon <- RunUMAP(object = colon, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))

#DimPlot(colon)
DimPlot(colon, reduction = "umap")
DimPlot(colon, reduction = "pca")
DimPlot(colon, reduction = "tsne")
TSNEPlot(object = colon)
dev.off()

if (copycat==1){
  colonENS <- FindNeighbors(colonENS, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))
  colonENS <- FindClusters(colonENS, resolution = 0.5)
  colonENS <- RunTSNE(object = colonENS, dims.use = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45), do.fast = TRUE)
  colonENS <- RunUMAP(object = colonENS, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))
  colonHG <- FindNeighbors(colonHG, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))
  colonHG <- FindClusters(colonHG, resolution = 0.5)
  colonHG <- RunTSNE(object = colonHG, dims.use = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45), do.fast = TRUE)
  colonHG <- RunUMAP(object = colonHG, dims = c(1,2,3,4,6,7,11,18,20,22,23,24,26,27,28,31,34,37,40,42,45))
  
  
}



#find marker genes
colon.markers <- FindAllMarkers(object = colon, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)
if (copycat ==1){
  colonENS.markers <- FindAllMarkers(object = colonENS, only.pos = TRUE, min.pct = 0.2, 
                                  thresh.use = 0.2)
}

#c02.markers <- FindMarkers(colon, ident.1 = "0", ident.2 = "2")

write.table(colon.markers,file = paste(a_name,paste("all","_DEGs.txt")),sep = "\t",quote = F)
#write.table(c02.markers,file = paste(a_name,paste("c0-c2","_DEGs.txt")),sep = "\t",quote = F)


####################
####pathway analysis  sometimes the internet connection hampers... and it fails. There will be a lot of warnings as we use only protein coding genes
####it's based on differential genes identified by Seurat what does not need to be the best choice, not all clusters get ####significant results - but when some pathways are checked manually they come back as cluster specific. I have to think ####about it a bit more...
DEGs_entrez_all=as.character(0)
DEGs_entrez_all=as.list(DEGs_entrez_all)


for (i in c(0:( max(as.numeric(colon$seurat_clusters))-1))){
  print(i)
#extract genes
  DEGs=sapply( strsplit( (colon.markers)[colon.markers$cluster==i,7], "-ENSG"), "[", 2)
  DEGs=paste("ENSG",DEGs,sep = "")
  #use only 500 top
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

  ###################### compare pathway scores between modules
  if (length(as.data.frame(PA)[,1])>1){
    pdf(file = paste(a_name,paste(i,"_Volcano.pdf")),height = 4, width = 4)
    npathways = length(as.data.frame(PA)[,1])
    if(50<npathways){npathways = 50}   #only tot 50 pathways are processed so some interesting ones might be skipped... 
    for (j in 1: npathways){
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
  
  
}


############custom pathwaysprojection
#xx <- as.list(reactomePATHID2EXTID)
PAS = c("R-HSA-170834","R-HSA-69478","R-HSA-69002","R-HSA-68952","R-HSA-69306","R-HSA-68962","R-HSA-383280","R-HSA-350054","R-HSA-212436","R-HSA-977606","R-HSA-71406","R-HSA-71403","R-HSA-5675482","R-HSA-203615","R-HSA-75105","R-HSA-75109","R-HSA-6798695","R-HSA-71240","R-HSA-71403","R-HSA-109581","R-HSA-2559583","R-HSA-77288","R-HSA-211935","R-HSA-77289","R-HSA-70171","R-HSA-71336")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

for (j in 1: length(PAS)){
  print(PAS[j])

  pathway=c(PAS[j])
  path_genes <- getBM(
    filters="reactome",
    attributes=c( "ensembl_gene_id"),
    values=pathway,
    mart=mart
  )
  if (length(path_genes[,1])>1){
    pdf(file = paste(a_name,paste(PAS[j],"_VLN.pdf")),height = 4, width = 4)
    colonENS= AddModuleScore(colonENS, 
                             features = path_genes, 
                             pool = NULL, 
                             nbin = 24, 
                             ctrl = 100,
                             k = FALSE, 
                             assay = NULL, name = as.data.frame(PAS)[j,], seed = 1)
    print(VlnPlot(colonENS, features = paste(as.data.frame(PAS)[j,],"1",sep = "") , ncol = 3))
    print (paste(j,length(as.data.frame(PAS)[j,]),sep = " pathway projection"))
    dev.off()
    pdf(file = paste(a_name,paste(PAS[j],"_fature_plots.pdf")),height = 3,width = 3)
    print(FeaturePlot(object = colonENS, features = c(paste(as.data.frame(PAS)[j,],"1",sep = "")),cols = colorpanel(100,"grey100","grey90","darkgreen")))
    dev.off()
    
  }
  
}

########################

#compare_pathways_clusters
cRes <- compareCluster(DEGs_entrez_all[-1], fun="enrichPathway")
pdf(file = paste(a_name,"_COmpare_clusters.pdf"),height = 12,width = 12)
dotplot(cRes,showCategory=30)
dev.off()

############
############

pdf(file = paste(a_name,"_VLN_plots_QC.pdf"),height = 3.5,width = 7)
VlnPlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(colon, features = c("second.best.cor", "tenth.best.cor", "mean.cor"), ncol = 3, same.y.lims = TRUE)


dev.off()

#################
##############

genesv=c(
  "ACAT1-ENSG00000075239",
 "ACADSB-ENSG00000196177",
  "KYNU-ENSG00000115919",
  "LPIN1-ENSG00000134324",
  "FGF13-ENSG00000129682",
 "NOS1-ENSG00000089250",
 "SOD2-ENSG00000112096",
 
  "ACTA2-ENSG00000107796",
  "MYH11-ENSG00000276480",
 "MYH10-ENSG00000133026",
 "IGFBP7-ENSG00000163453"   ,

 
 "CXCL12-ENSG00000107562",
 
 "nFeature_RNA", "nCount_RNA", "percent.mt",
 "second.best.cor", "tenth.best.cor", "mean.cor",

 
 

  "CD14-ENSG00000170458",
  "C1QA-ENSG00000173372",
  "CD63-ENSG00000135404"  ,   
   "CD74-ENSG00000019582" ,
 "APOE-ENSG00000130203",
  "USP8-ENSG00000138592"  ,    
 "ZNF286A-ENSG00000187607"  ,
 "ATXN3-ENSG00000066427"   ,  "RPS6KA5-ENSG00000100784" ,  "SLC35E3-ENSG00000175782",  
 "SLC35E3-ENSG00000175782",
 "VDR-ENSG00000111424",
  "ORAI2-ENSG00000160991",
 "MYH11..Smooth.Muscle.Cells",
 "CD3.CD8.",
 "CD3.CD4.",
 "CD34.",
  "CD14.CD68."
 
  
  )

pathways = c("Platelet degranulation 1","MHC class II antigen presentation1",
             "Extracellular matrix organization1","Collagen degradation1",
             "Iron uptake and transport1","Neutrophil degranulation1"
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
  print(FeaturePlot(object = colon, features = c(genesv[i]),cols = colorpanel(100,"grey100","grey90","darkgreen")))
  dev.off()
}


for (i in 1:length(pathways)){
  print(pathways[i])
  pdf(file = paste(paste(a_name,pathways[i]),"_fature_plots.pdf"),height = 3,width = 3)
  print(FeaturePlot(object = colonENS, features = c(pathways[i]),cols = colorpanel(100,"grey100","grey90","darkgreen")))
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
VlnPlot(colon,x.lab.rot=T, features = c("SPIC-ENSG00000166211"))



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
FeaturePlot(object = colonENS, features = "Neutrophil degranulation1",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MYH11..Smooth.Muscle.Cells",cols = topo.colors(100))
FeaturePlot(object = colon, features = "VDR-ENSG00000111424",cols = topo.colors(100))

FeaturePlot(object = colon, features = "SPIC-ENSG00000166211",cols = topo.colors(100))





dev.off()


library("dplyr")
top10 <- colon.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf(file = paste(a_name,"_heatmap_top10markers.pdf"),height = 8,width = 8)
DoHeatmap(colon, features = top10$gene) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))
dev.off()
pdf(file = paste(a_name,"_heatmap_fewselectedmarkers.pdf"),height = 3.5,width = 8)
DoHeatmap(colon, features = genesv) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))
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
                              "ORAI2-ENSG00000160991"  ,   "VDR-ENSG00000111424",    
                              
                               "CD81-ENSG00000110651"   ,   "C1QA-ENSG00000173372"    ,  "CD63-ENSG00000135404"  ,   
                               "C1QB-ENSG00000173369"   ,   "PLTP-ENSG00000100979"    ,  "GRN-ENSG00000030582"   ,   
                               "SRGN-ENSG00000122862"   ,   "TYROBP-ENSG00000011600"  ,  "CD74-ENSG00000019582"  ,   
                               "SAT1-ENSG00000130066"      , 
                              
                              "CTSB-ENSG00000164733"   ,  
                              "IFI30-ENSG00000216490"  ,   "FTH1-ENSG00000167996"    ,  "PSAP-ENSG00000197746"    , 
                              "FTL-ENSG00000087086"    ,   "GPNMB-ENSG00000136235"  ,   "TYROBP-ENSG00000011600"  , 
                         
                              "APOE-ENSG00000130203"   ,   "APOC1-ENSG00000130208"  ,  
                              "SRGN-ENSG00000122862"   ,   "PKM-ENSG00000067225" ,     
                              "IGLL5-ENSG00000254709"    )) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))

colon.cluster.averages <- AverageExpression(colon, return.seurat = TRUE)
DoHeatmap(colon.cluster.averages, features = c("CCDC144A-ENSG00000170160","CYSLTR1-ENSG00000173198",
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
                              "ORAI2-ENSG00000160991"  ,   "VDR-ENSG00000111424",    
                              
                              "CD81-ENSG00000110651"   ,   "C1QA-ENSG00000173372"    ,  "CD63-ENSG00000135404"  ,   
                              "C1QB-ENSG00000173369"   ,   "PLTP-ENSG00000100979"    ,  "GRN-ENSG00000030582"   ,   
                              "SRGN-ENSG00000122862"   ,   "TYROBP-ENSG00000011600"  ,  "CD74-ENSG00000019582"  ,   
                              "SAT1-ENSG00000130066"      , 
                              
                              "CTSB-ENSG00000164733"   ,  
                              "IFI30-ENSG00000216490"  ,   "FTH1-ENSG00000167996"    ,  "PSAP-ENSG00000197746"    , 
                              "FTL-ENSG00000087086"    ,   "GPNMB-ENSG00000136235"  ,   "TYROBP-ENSG00000011600"  , 
                              
                              "APOE-ENSG00000130203"   ,   "APOC1-ENSG00000130208"  ,  
                              "SRGN-ENSG00000122862"   ,   "PKM-ENSG00000067225" ,     
                              "IGLL5-ENSG00000254709"    ),draw.lines = F) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))

dev.off()

#correlate_data_with_single_cell_deconvolution, this will be fixed...
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




###################
#correlate results with clinical data   this one will be changed for more proper method
pdf(file = paste(a_name,"_clinical_data.pdf"),height = 5,width = 5)
f=read.table(file = "clinical_data_good_selection.txt",header = T,sep = "\t",row.names = 1)

for (i in c(1,2,4,5,7:18)){
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




## single R
library(SingleR)
singler <- CreateSinglerObject(counts = colonHG@assays$RNA@counts[,rownames(colonHG@meta.data)],       
                               annot = Idents(colon),                                        # if is NULL it takes 10 clusters
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
                               clusters = Idents(colon),
                               do.main.types = TRUE, 
                               reduce.file.size = TRUE, 
                               numCores = 4
)
pdf(file = paste(a_name,"_SingleR_data.pdf"),height = 5,width = 4)
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 25)
SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 25)

dev.off()

pdf(file = paste(a_name,"_SingleR_individual cells_data.pdf"),height = 5,width = 7)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n = 25, clusters = Idents(colon),order.by.clusters = T)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main,top.n = 25, clusters = Idents(colon),order.by.clusters = T)

dev.off()



#####survival
library(ggfortify)
library(survival)


pdf(file = paste(a_name,"_survival_ep_major.pdf"),height = 5,width = 5)
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
autoplot(fit, conf.int = F,ylim = c(0.75,1),xlim = c(0,5))
combined$Idents.colon.[combined$Idents.colon. == 3] = 1
fit = survfit(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined)
autoplot(fit, conf.int = F,ylim = c(0.75,1),xlim = c(0,5))

combined$Idents.colon.[combined$Idents.colon. == 4] = 0
combined$Idents.colon.[combined$Idents.colon. == 2] = 0



fit = survfit(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined,conf.int=TRUE)
survdiff(Surv(ep_major_time,ep.major)~Idents.colon.,data = combined)
autoplot(fit, conf.int = F,ylim = c(0.75,1),xlim = c(0,5))
autoplot(fit, conf.int = T,ylim = c(0.75,1),xlim = c(0,5))




dev.off()

############### GWAS

#CAD=read.table(file = "UKBB_CAD.genes.txt",header = T,sep = "\t",row.names = 1)[,1]


for (thr in c(0.05,0.01,0.001,0.0001,0.00001, 0.000001)){
CAD=read.table(file = "UKBB_CAD.MAGMAGenes.txt",header = T,sep = "\t")
CAD = CAD[CAD[,9]<thr,1]


print(length(CAD))

pdf(file = paste(paste(a_name,thr),"GWAS.pdf"),height = 3,width = 12)
par(mfrow=c(1,5))
for (i in 0:4){
  cluster.genes = colonENS.markers[colonENS.markers[,6] == i,7]
  print(length(cluster.genes))
  if (length(cluster.genes)>250){
    cluster.genes = cluster.genes[1:250]
  }
  print(length(cluster.genes))
  overlap = intersect(CAD,cluster.genes)

  m=0
  c=0
  set.seed(123)
  for (j in 1:10000){
    random.genes = sample(row.names(raw.genecountsENS),length(cluster.genes),replace = T)
    random.overlap = intersect(CAD,random.genes)
    m[j]=length(random.overlap)
    if (length(overlap)>length(random.overlap)){
      c=c+1
    }
  }

  
  (hist.default(m, xlim = c(0,(2*max(m))), breaks = (max(m)+1),col = scales::hue_pal()(5)[i+1],main = paste("p = ",1-c/10000)))
  (points((length(overlap)),0, pch = 20, cex = 4, col = "red"))
  abline(v=(length(overlap)), col = "red",lwd = 4)
}
dev.off()
}
############# baseline table
#install.packages("tidyverse")
library(tidyverse)
library(tableone)
#write.table(data.frame(Idents(colon)),file = "Seurat_clusters_Michal_vXX.txt" ,sep = "\t",quote = F)
seurat_clusters <- read_tsv("Seurat_clusters_Michal_v13.txt")
seurat_PCA <- read_tsv("PCAs_AE_bulk_v13.txt")
scaden_predictions <- read_tsv("scaden_predictions_model.txt")

seurat_clusters  <- seurat_clusters %>% 
  left_join(seurat_PCA, by = "study_number")
seurat_clusters  <- seurat_clusters %>% 
  left_join(scaden_predictions, by = "study_number")

clinical_data <- read_tsv("bulk_RNAseq_clinical_data_MM.txt")

clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number, cluster, sex, smokercurrent, dm_composite, risk614, hypertension1, cad_history, stroke_history, paod, symptoms_4g, stenosis_ipsilateral, stenosis_con_bin, med_statin_derived, med_statin_lld, med_all_antiplatelet, age, bmi, gfr_mdrd, totalchol, triglyceriden, ldl, hdl, plaquephenotype, epmajor_3years, ep_major, time_event_or,fabp4_serum_luminex,	cst3_pg_ug,	cst3_serum_luminex,	hdac9,	vegfa_plasma,	vwf_plasma,	pla2_plasma,	pcsk9_plasma,	gdf15_plasma,	ctni_plasma,	rantes_plasma,	hscrp_plasma,	tat_plasma,	mpo_plasma,	ntprobnp_plasma,	pdgf_bb_plasma,	opg_plasma, il6,symptoms_2g,macmean0,	smcmean0, fat, thrombus,	thrombus_organization, thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	neutrophils,	mast_cells_plaque,	vessel_density, iph_bin,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,CD3CD8,CD14CD68, MYH11SmoothMuscleCells,CD79ABCells,CD3CD4,CD34,KITMastCells) %>% 
  mutate_at(vars(cluster, sex, smokercurrent, dm_composite, risk614, hypertension1, cad_history, stroke_history, paod, symptoms_4g, stenosis_ipsilateral, stenosis_con_bin, med_statin_derived, med_statin_lld, med_all_antiplatelet, plaquephenotype, epmajor_3years, ep_major,symptoms_2g,fat, thrombus,	thrombus_organization,	thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,iph_bin), list(as.factor)) %>%
  print()




listVars <- c("sex", "smokercurrent", "dm_composite","hypertension1", "cad_history", "stroke_history", "paod", "symptoms_4g", "stenosis_ipsilateral", "stenosis_con_bin", "med_statin_derived", "med_statin_lld", "med_all_antiplatelet", "age", "bmi", "gfr_mdrd", "totalchol", "triglyceriden", "ldl", "hdl", "plaquephenotype", "epmajor_3years", "ep_major","time_event_or","fabp4_serum_luminex", "cst3_pg_ug", "cst3_serum_luminex", "hdac9", "vegfa_plasma", "vwf_plasma", 	"pla2_plasma", "pcsk9_plasma", "gdf15_plasma", "ctni_plasma","rantes_plasma", "hscrp_plasma", "tat_plasma", "mpo_plasma", "ntprobnp_plasma", "pdgf_bb_plasma", "opg_plasma", "il6")
catVars <- c("cluster", "sex", "smokercurrent", "dm_composite", "risk614", "hypertension1", "cad_history", "stroke_history", "paod", "symptoms_4g", "stenosis_ipsilateral", "stenosis_con_bin", "med_statin_derived", "med_statin_lld", "med_all_antiplatelet", "plaquephenotype", "epmajor_3years", "ep_major")


table1 <- CreateTableOne(vars = listVars, data = clinical_sel, factorVars = catVars)
table2 <- CreateTableOne(listVars, clinical_sel, catVars, strata = c("cluster"))
beeswarm(clinical_sel$time_event_or ~ clinical_sel$cluster)
beeswarm(clinical_sel$hscrp_plasma ~ clinical_sel$cluster,ylim = c(0,1000))
table_female <- CreateTableOne(listVars, clinical_sel[clinical_sel$sex=="female",], catVars, strata = c("cluster"))
table3 <- CreateTableOne(listVars, clinical_sel, catVars, strata = c("plaquephenotype"))
table4 <- CreateTableOne(listVars, clinical_sel, catVars, strata = c("cluster","plaquephenotype"))
table5 <- CreateTableOne(listVars, clinical_sel, catVars, strata = c("cluster","sex"))





################################33
##############   logistic model####
###################################
#check male female separatelly, 
#take only one category e.g fibrous   and see if RNA adds something.
#repeat even the clustering only on one plaque type
#PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+PC_11+PC_12+PC_13+PC_14+PC15

library(tidyverse)
library(tableone)
#write.table(data.frame(Idents(colon)),file = "Seurat_clusters_Michal_vXX.txt" ,sep = "\t",quote = F)
seurat_clusters <- read_tsv("Seurat_clusters_Michal_v13.txt")
seurat_PCA <- read_tsv("PCAs_AE_bulk_v13.txt")
scaden_predictions <- read_tsv("scaden_predictions_model.txt")

seurat_clusters  <- seurat_clusters %>% 
  left_join(seurat_PCA, by = "study_number")
seurat_clusters  <- seurat_clusters %>% 
  left_join(scaden_predictions, by = "study_number")

clinical_data <- read_tsv("bulk_RNAseq_clinical_data_MM.txt")

clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g, plaquephenotype ,cluster, sex, macmean0,	smcmean0, fat, thrombus,	thrombus_organization, thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50  )  %>%
  mutate_at(vars(symptoms_2g, plaquephenotype ,cluster, sex, fat, thrombus,	thrombus_organization, thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)

pdf(file = paste(a_name,"GLM_SYmptoms_2g_selected_variables.pdf"),height = 7,width = 7)
par(mfrow=c(2,2))


###################HIST
clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g, plaquephenotype, macmean0,	smcmean0, fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin )  %>%
  mutate_at(vars(symptoms_2g, plaquephenotype , fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)
clinical_sel <- clinical_sel %>% 
  drop_na()
summary(clinical_sel)

histology = "symptoms_2g ~ plaquephenotype+ macmean0+	smcmean0+ fat+	calcification+	collagen	+smc	+smc_location	+macrophages+	macrophages_location+	smc_macrophages_ratio+	media+ iph_bin "

model2 <- glm(histology, data = clinical_sel, family = binomial)
summary(model2)
steps=step(model2)
histology = "symptoms_2g ~ macmean0 + fat + calcification + smc"
model2 <- glm(histology, data = clinical_sel, family = binomial)
summary(model2)
predictions <- as.data.frame(1-predict(model2, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.hist = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "histology",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("orange","red"), pch = 20,cex = 0.6)

########RNA

clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g, PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50  )  %>%
  mutate_at(vars(symptoms_2g ), list(as.factor))  %>%
  print()
summary(clinical_sel)


RNA_PCA = "symptoms_2g ~ PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+PC_11+PC_12+PC_13+PC_14+PC_15+PC_16+PC_17+PC_18+PC_19+PC_20+PC_21+PC_22+PC_23+PC_24+PC_25+PC_26+PC_27+PC_28+PC_29+PC_30+PC_31+PC_32+PC_33+PC_34+PC_35+PC_36+PC_37+PC_38+PC_39+PC_40+PC_41+PC_42+PC_43+PC_44+PC_45+PC_46+PC_47+PC_48+PC_49+PC_50"
model3 <- glm(RNA_PCA, data = clinical_sel, family = binomial)
summary(model3)
steps=step(model3)
RNA_PCA = "symptoms_2g ~ PC_1 + PC_3 + PC_7 + PC_14 + PC_16 + PC_20 + PC_21 + 
  PC_23 + PC_24 + PC_26 + PC_28 + PC_30 + PC_34 + PC_41 + PC_44 + 
  PC_45 + PC_49"
model3 <- glm(RNA_PCA, data = clinical_sel, family = binomial)
summary(model3)
predictions <- as.data.frame(1-predict(model3, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.rna = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "RNA",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("lightblue","darkblue"), pch = 20,cex = 0.6)

####ALL


clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g, plaquephenotype, macmean0,	smcmean0, fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin , PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50 )  %>%
  mutate_at(vars(symptoms_2g, plaquephenotype , fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)
clinical_sel <- clinical_sel %>% 
  drop_na()
summary(clinical_sel)

all_parameters = "symptoms_2g ~ plaquephenotype+ macmean0+	smcmean0+ fat+	calcification+	collagen	+smc	+smc_location	+macrophages+	macrophages_location+	smc_macrophages_ratio+	media+ iph_bin+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+PC_11+PC_12+PC_13+PC_14+PC_15+PC_16+PC_17+PC_18+PC_19+PC_20+PC_21+PC_22+PC_23+PC_24+PC_25+PC_26+PC_27+PC_28+PC_29+PC_30+PC_31+PC_32+PC_33+PC_34+PC_35+PC_36+PC_37+PC_38+PC_39+PC_40+PC_41+PC_42+PC_43+PC_44+PC_45+PC_46+PC_47+PC_48+PC_49+PC_50"


model1 <- glm(all_parameters, data = clinical_sel, family = binomial)
summary(model1)
steps=step(model1)
all_parameters = "symptoms_2g ~ fat + calcification + PC_1 + PC_3 + PC_11 + PC_14 + 
    PC_15 + PC_16 + PC_20 + PC_21 + PC_23 + PC_28 + PC_30 + PC_34 + 
    PC_41 + PC_44 + PC_45 + PC_49"
model1 <- glm(all_parameters, data = clinical_sel, family = binomial)
summary(model1)
predictions <- as.data.frame(1-predict(model1, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.all = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "histology + RNA",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("lightgreen","darkgreen"), pch = 20,cex = 0.6)



#model4 <- glm(in_silico_histology, data = clinical_sel, family = binomial)
#summary(model4)
#steps=step(model4)
#predictions <- as.data.frame(1-predict(model1, clinical_sel, type = "response"))
#names(predictions)="mild"
#predictions$severe=1-predictions$mild
#predictions$observed <- clinical_sel$symptoms_2g
#predictions = predictions[complete.cases(predictions),]
#predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
#head(predictions)
#predictions = predictions[complete.cases(predictions),]
#roc.ishistology = roc(predictions$observed, predictions$severe)
#boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "IS histology",ylim = c(0,1))
#beeswarm(severe~observed,data = predictions,add = T, col = c("lightpurple","purple"), pch = 20,cex = 0.6)


plot(roc.all, col = "green", xlim = c(1,0))
lines(roc.rna, col = "blue")
lines(roc.hist, col = "red")
legend(legend = c(paste ("H AUC = ",round(roc.hist$auc,digits = 3)),
                  paste ("R AUC = ",round(roc.rna$auc,digits = 3)),
                  paste ("R+H AUC = ",round(roc.all$auc,digits = 3))
                  ),x=0.6,y=0.3,col = c("red","blue","green"),pch = 19)
dev.off()



################################33
##############   logistic model no TIA ####
###################################
#clinical_sel = clinical_sel[clinical_sel$symptoms_4g!="TIA",]
pdf(file = paste(a_name,"GLM_SYmptoms_2g_noTIA_selected_vars.pdf"),height = 7,width = 7)
par(mfrow=c(2,2))
seurat_clusters <- read_tsv("Seurat_clusters_Michal_v13.txt")
seurat_PCA <- read_tsv("PCAs_AE_bulk_v13.txt")
scaden_predictions <- read_tsv("scaden_predictions_model.txt")

seurat_clusters  <- seurat_clusters %>% 
  left_join(seurat_PCA, by = "study_number")
seurat_clusters  <- seurat_clusters %>% 
  left_join(scaden_predictions, by = "study_number")

clinical_data <- read_tsv("bulk_RNAseq_clinical_data_MM.txt")

clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_4g,symptoms_2g, plaquephenotype ,cluster, sex, macmean0,	smcmean0, fat, thrombus,	thrombus_organization, thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50  )  %>%
  mutate_at(vars(symptoms_2g, symptoms_4g,plaquephenotype ,cluster, sex, fat, thrombus,	thrombus_organization, thrombus_organization_v2,	thrombus_location,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)
clinical_sel = clinical_sel[clinical_sel$symptoms_4g!="TIA",]


###################HIST
clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g,symptoms_4g, plaquephenotype, macmean0,	smcmean0, fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin )  %>%
  mutate_at(vars(symptoms_2g,symptoms_4g, plaquephenotype , fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)
clinical_sel <- clinical_sel %>% 
  drop_na()
clinical_sel = clinical_sel[clinical_sel$symptoms_4g!="TIA",]
summary(clinical_sel)

histology = "symptoms_2g ~ plaquephenotype+ macmean0+	smcmean0+ fat+	calcification+	collagen	+smc	+smc_location	+macrophages+	macrophages_location+	smc_macrophages_ratio+	media+ iph_bin "

model2 <- glm(histology, data = clinical_sel, family = binomial)
summary(model2)
steps=step(model2)
histology = "symptoms_2g ~ plaquephenotype + macmean0 + calcification + smc"
model2 <- glm(histology, data = clinical_sel, family = binomial)
summary(model2)
predictions <- as.data.frame(1-predict(model2, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.hist = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "histology",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("orange","red"), pch = 20,cex = 0.6)

########RNA

clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g,symptoms_4g, PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50  )  %>%
  mutate_at(vars(symptoms_2g ,symptoms_4g), list(as.factor))  %>%
  print()
clinical_sel = clinical_sel[clinical_sel$symptoms_4g!="TIA",]
summary(clinical_sel)


RNA_PCA = "symptoms_2g ~ PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+PC_11+PC_12+PC_13+PC_14+PC_15+PC_16+PC_17+PC_18+PC_19+PC_20+PC_21+PC_22+PC_23+PC_24+PC_25+PC_26+PC_27+PC_28+PC_29+PC_30+PC_31+PC_32+PC_33+PC_34+PC_35+PC_36+PC_37+PC_38+PC_39+PC_40+PC_41+PC_42+PC_43+PC_44+PC_45+PC_46+PC_47+PC_48+PC_49+PC_50"
model3 <- glm(RNA_PCA, data = clinical_sel, family = binomial)
summary(model3)
steps=step(model3)
RNA_PCA = "symptoms_2g ~ PC_1 + PC_2 + PC_3 + PC_6 + PC_7 + PC_10 + PC_12 + 
    PC_13 + PC_14 + PC_17 + PC_18 + PC_19 + PC_20 + PC_21 + PC_22 + 
PC_26 + PC_29 + PC_33 + PC_34 + PC_41 + PC_43 + PC_45 + PC_47"
model3 <- glm(RNA_PCA, data = clinical_sel, family = binomial)
summary(model3)
predictions <- as.data.frame(1-predict(model3, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.rna = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "RNA",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("lightblue","darkblue"), pch = 20,cex = 0.6)

####ALL


clinical_sel <- clinical_data %>% 
  left_join(seurat_clusters, by = "study_number") %>% 
  dplyr::select(study_number,symptoms_2g,symptoms_4g, plaquephenotype, macmean0,	smcmean0, fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media, iph_bin , PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20,PC_21,PC_22,PC_23,PC_24,PC_25,PC_26,PC_27,PC_28,PC_29,PC_30,PC_31,PC_32,PC_33,PC_34,PC_35,PC_36,PC_37,PC_38,PC_39,PC_40,PC_41,PC_42,PC_43,PC_44,PC_45,PC_46,PC_47,PC_48,PC_49,PC_50 )  %>%
  mutate_at(vars(symptoms_2g,symptoms_4g, plaquephenotype , fat,	calcification,	collagen	,smc	,smc_location	,macrophages,	macrophages_location,	smc_macrophages_ratio,	media,	iph_bin  ), list(as.factor))  %>%
  print()
summary(clinical_sel)
clinical_sel <- clinical_sel %>% 
  drop_na()
clinical_sel = clinical_sel[clinical_sel$symptoms_4g!="TIA",]
summary(clinical_sel)

all_parameters = "symptoms_2g ~ plaquephenotype+ macmean0+	smcmean0+ fat+	calcification+	collagen	+smc	+smc_location	+macrophages+	macrophages_location+	smc_macrophages_ratio+	media+ iph_bin+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10+PC_11+PC_12+PC_13+PC_14+PC_15+PC_16+PC_17+PC_18+PC_19+PC_20+PC_21+PC_22+PC_23+PC_24+PC_25+PC_26+PC_27+PC_28+PC_29+PC_30+PC_31+PC_32+PC_33+PC_34+PC_35+PC_36+PC_37+PC_38+PC_39+PC_40+PC_41+PC_42+PC_43+PC_44+PC_45+PC_46+PC_47+PC_48+PC_49+PC_50"


model1 <- glm(all_parameters, data = clinical_sel, family = binomial)
summary(model1)
steps=step(model1)
all_parameters = "symptoms_2g ~ plaquephenotype + fat + calcification + PC_1 + 
    PC_3 + PC_6 + PC_10 + PC_12 + PC_14 + PC_16 + PC_17 + PC_18 + 
    PC_19 + PC_20 + PC_21 + PC_22 + PC_26 + PC_28 + PC_29 + PC_34 + 
    PC_37 + PC_39 + PC_41 + PC_43 + PC_45 + PC_47"
model1 <- glm(all_parameters, data = clinical_sel, family = binomial)
summary(model1)
predictions <- as.data.frame(1-predict(model1, clinical_sel, type = "response"))
names(predictions)="mild"
predictions$severe=1-predictions$mild
predictions$observed <- clinical_sel$symptoms_2g
predictions = predictions[complete.cases(predictions),]
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
predictions = predictions[complete.cases(predictions),]
roc.all = roc(predictions$observed, predictions$severe)
boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "histology + RNA",ylim = c(0,1))
beeswarm(severe~observed,data = predictions,add = T, col = c("lightgreen","darkgreen"), pch = 20,cex = 0.6)



#model4 <- glm(in_silico_histology, data = clinical_sel, family = binomial)
#summary(model4)
#steps=step(model4)
#predictions <- as.data.frame(1-predict(model1, clinical_sel, type = "response"))
#names(predictions)="mild"
#predictions$severe=1-predictions$mild
#predictions$observed <- clinical_sel$symptoms_2g
#predictions = predictions[complete.cases(predictions),]
#predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
#head(predictions)
#predictions = predictions[complete.cases(predictions),]
#roc.ishistology = roc(predictions$observed, predictions$severe)
#boxplot(severe~observed,data = predictions,outline = F, ylab = ("prediction score"),main = "IS histology",ylim = c(0,1))
#beeswarm(severe~observed,data = predictions,add = T, col = c("lightpurple","purple"), pch = 20,cex = 0.6)


plot(roc.all, col = "green", xlim = c(1,0))
lines(roc.rna, col = "blue")
lines(roc.hist, col = "red")
legend(legend = c(paste ("H AUC = ",round(roc.hist$auc,digits = 3)),
                  paste ("R AUC = ",round(roc.rna$auc,digits = 3)),
                  paste ("R+H AUC = ",round(roc.all$auc,digits = 3))
),x=0.6,y=0.3,col = c("red","blue","green"),pch = 19)
dev.off()









##########################PCA loadings
PCA_load = Loadings(colonENS, reduction = "pca")
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "useast")


for (i in 1:100){
  print(i)
  if((i %% 2) == 0){j=i/2
    #extract genes
    DEGs=names(PCA_load[order(PCA_lead[,j],decreasing = F),1][1:100])
  } else{
      j=i/2+0.5
    #extract genes
    DEGs=names(PCA_load[order(PCA_lead[,j],decreasing = T),1][1:100])
    }
  print(j)  
  
  
  DEGs_entrez <- getBM(
  filters="ensembl_gene_id",
  attributes=c( "entrezgene_id"),
  values=DEGs,
  mart=mart
  )
  DEGs_entrez_all[[as.character(i)]]=as.vector(DEGs_entrez[,1])
  PA <- enrichPathway(gene=as.vector(DEGs_entrez[,1]),pvalueCutoff=0.05,minGSSize = 5,maxGSSize = 500, readable=T,pAdjustMethod="none")
  if (length(as.data.frame(PA)[,1])>1){
    pdf(file = paste(a_name,paste(paste(i,j,sep="_"),"PCA_Reactome_100genes.pdf")),height = 7, width = 13)
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
  }
}

#only those used in TIA model
cRes <- compareCluster(DEGs_entrez_all[-1][c(1,2,5,6,11,12,19,20,27,28,31:44,51,52,55:58,67,68,73,74,77,78,81,82,85,86,89,90,93,94)], fun="enrichPathway",minGSSize = 5,maxGSSize = 500,pvalueCutoff=0.05,pAdjustMethod="none")
pdf(file = paste(a_name,"_COmpare_PCAs_pathway_100_genes.pdf"),height = 26,width = 18)
dotplot(cRes,showCategory=30)
dev.off()


###################


library(randomForest)
library(pROC)
set.seed(100)
clinical_sel_model = clinical_sel[100:300,]
output.forest <- randomForest(symptoms_2g ~  cluster+plaquephenotype+macmean0 + smcmean0 , ntree = 500,
                              data = clinical_sel_model, na.action = na.omit,  mtry = 4, importance = TRUE)
print(output.forest)
print(importance(output.forest,type = 2))
varImpPlot(output.forest)
plot(output.forest)
predictions <- as.data.frame(predict(output.forest, clinical_sel[c(1:99,501:600),], type = "prob"))
predictions$observed <- clinical_sel[c(1:99,401:500),]$symptoms_2g
predictions = predictions[complete.cases(predictions),]

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
head(predictions)
roc.all = roc(predictions$observed, predictions$severe,predictions$predicted)
plot(roc.all, col = "gray60")


################
#################
pdf(file = paste(a_name,"dists_symptoms_4g_cluster.pdf"),height = 6,width = 8)
par(mfrow=c(2,3),mar = c(0,0,4,0))
for (i in c(3,4,2,5,1)){
  print (pie(table(clinical_sel$symptoms_4g, clinical_sel$cluster)[,i],main = paste ("cluster ",colnames(table(clinical_sel$symptoms_4g, clinical_sel$cluster))[i],sep = ""),col = c(1,2,3,4),cex = 1.3))
}
dev.off()


pdf(file = paste(a_name,"dists_symptoms_4g_plaque_phenotype.pdf"),height = 3,width = 8)
par(mfrow=c(1,3),mar = c(0,0,4,0))
for (i in 1:3){
  print (pie(table(clinical_sel$symptoms_4g, clinical_sel$plaquephenotype)[,i],main = paste ("cluster ",colnames(table(clinical_sel$symptoms_4g, clinical_sel$plaquephenotype))[i],sep = ""),col = c(1,2,3,4),cex = 1.3))
}
dev.off()

pdf(file = paste(a_name,"dists_symptoms_4g_plaque_phenotyp_clustere.pdf"),height = 7,width = 10)
par(mfrow=c(3,5),mar = c(0,0,4,0))
for (j in 1:3){
 for (i in c(3,4,2,5,1)){


    print (pie(table(clinical_sel$cluster,clinical_sel$symptoms_4g, clinical_sel$plaquephenotype)[i,,j],main =paste("n=",sum(table(clinical_sel$cluster,clinical_sel$symptoms_4g, clinical_sel$plaquephenotype)[i,,j])),col = c(1,2,3,4),cex = 1.3))
    
    
  }
}
dev.off()
pdf(file = paste(a_name,"dists_ep_major3y_plaque_phenotyp_clustere.pdf"),height = 7,width = 10)
par(mfrow=c(3,5),mar = c(0,0,4,0))
for (j in 1:3){
  for (i in c(3,4,2,5,1)){
    
    
    print (pie(table(clinical_sel$cluster,clinical_sel$epmajor_3years, clinical_sel$plaquephenotype)[i,,j],main =paste("n=",sum(table(clinical_sel$cluster,clinical_sel$epmajor_3years, clinical_sel$plaquephenotype)[i,,j])),col = c(1,2,3,4),cex = 1.3))
    
    
  }
}
dev.off()
pdf(file = paste(a_name,"dists_sex_plaque_phenotyp_clustere.pdf"),height = 7,width = 10)
par(mfrow=c(3,5),mar = c(0,0,4,0))
for (j in 1:3){
  for (i in c(3,4,2,5,1)){
    
    
    print (pie(table(clinical_sel$cluster,clinical_sel$sex, clinical_sel$plaquephenotype)[i,,j],main =paste("n=",sum(table(clinical_sel$cluster,clinical_sel$sex, clinical_sel$plaquephenotype)[i,,j])),col = c(1,2,3,4),cex = 1.3))
    
    
  }
}
dev.off()


pdf(file = paste(a_name,"dists_paod_plaque_phenotyp_clustere.pdf"),height = 7,width = 10)
par(mfrow=c(3,5),mar = c(0,0,4,0))
for (j in 1:3){
  for (i in c(3,4,2,5,1)){
    
    
    print (pie(table(clinical_sel$cluster,clinical_sel$paod, clinical_sel$plaquephenotype)[i,,j],main =paste("n=",sum(table(clinical_sel$cluster,clinical_sel$paod, clinical_sel$plaquephenotype)[i,,j])),col = c(1,2,3,4),cex = 1.3))
    
    
  }
}
dev.off()

pdf(file = paste(a_name,"dists_med_statin_derived_plaque_phenotyp_clustere.pdf"),height = 7,width = 10)
par(mfrow=c(3,5),mar = c(0,0,4,0))
for (j in 1:3){
  for (i in c(3,4,2,5,1)){
    
    
    print (pie(table(clinical_sel$cluster,clinical_sel$med_statin_derived, clinical_sel$plaquephenotype)[i,,j],main =paste("n=",sum(table(clinical_sel$cluster,clinical_sel$med_statin_derived, clinical_sel$plaquephenotype)[i,,j])),col = c(1,2,3,4),cex = 1.3))
    
    
  }
}
dev.off()
med_statin_derived


table(clinical_sel$plaquephenotype, clinical_sel$cluster)[,c(3,4,2,5,1)]
############GWAS2

## loop over different clusters and save differential exression data
# create dataframe
differential.expression.clusters <- data.frame()

cluster.number = 5
identity.classes = c(0,1,2,3,4)
# prepare all possible combinations
choices <- choose(cluster.number, 2)        # identify the number of possible pairs
combinations <- combn(cluster.number, 2)    # create those pairs

# pairwise comparison for differential expression
for (choice in 1:choices) {
  # select the cluster numbers to run
  cluster1 <- combinations[,choice][1] # the first index of the combination
  cluster2 <- combinations[,choice][2] # the second index of the combination
  
  # input function must be character
  cluster1 <- identity.classes[cluster1]
  cluster2 <- identity.classes[cluster2]
  
  # run formula
  pairwise.comparison <- FindMarkers(object = colonHG,
                                     ident.1 = cluster1 ,
                                     ident.2 = cluster2,
                                     logfc.threshold = 0.6, #default 0.25
                                     min.pct = 0.1,
                                     test.use = "wilcox") # default wilcox
  
  # add gene names for later
  pairwise.comparison$gene <- rownames(pairwise.comparison)
  # add clusters
  pairwise.comparison$cluster1 <- cluster1
  pairwise.comparison$cluster2 <- cluster2
  # select only values with p adjusted of certain value
  pairwise.comparison <- pairwise.comparison[(pairwise.comparison$p_val_adj <= .05), ]
  # bind to dataframe
  differential.expression.clusters <- rbind(differential.expression.clusters, pairwise.comparison)
  
}

dim(differential.expression.clusters)

## differential expression findAllMarkers
differential.expression.all <- FindAllMarkers(object = colonHG,
                                              logfc.threshold = 0.6,
                                              min.pct = 0.1,
                                              only.pos = FALSE,
                                              test.use = "wilcox")

dim(differential.expression.all)
differential.expression.all <- differential.expression.all[(differential.expression.all$p_val_adj <= 0.05),]
dim(differential.expression.all)

## combine differential expression based on gene names, not rownames
differential.expression.data <- union(differential.expression.all$gene, differential.expression.clusters$gene)
length(differential.expression.data)
length(unique(differential.expression.data))

## save the table with results
# write.table(differential.expression.clusters, "DEG_one_vs_one.txt")
# write.table(differential.expression.all, "DEG_one_vs_all.txt")
# write.table(differential.expression.data, "All_DEG_(gene_names).txt")

#-----------------------------------------------------------------------------------------
### K means
## K means with average expression

# calculate average expression with the genes from differential.expression.data that are present in the object
average.expression <-  AverageExpression(colonHG,
                                         features = differential.expression.data,
                                         return.seurat = TRUE)

# save object
#saveRDS(average.expression, "seuset_average_expression_DEG.RDS")

## K means
geneclusters <- 8

## calculate Kmeans
# option 1 (seurat heatmap)

set.seed(576)
Kmeans.object <- kmeans(average.expression@assays$RNA@scale.data, centers = geneclusters)
Kmeans.object$cluster["MYH11"]
image(Kmeans.object$centers)
heatmap(Kmeans.object$centers, col = colorpanel(100,"grey100","grey90","darkgreen"))
Kmeans.object$cluster[Kmeans.object$cluster==7]



for (thr in c(0.05,0.01,0.001,0.0001,0.00001, 0.000001)){
  CAD=read.table(file = "UKBB_CAD.MAGMAGenes.txt",header = T,sep = "\t")
  CAD = CAD[CAD[,9]<thr,10]
  
  
  print(length(CAD))
  
  pdf(file = paste(paste(a_name,thr),"GWAS_kmean.pdf"),height = 3,width = 12)
  par(mfrow=c(2,4))
  for (i in 1:geneclusters){
    cluster.genes = Kmeans.object$cluster[Kmeans.object$cluster==i]
    print(length(cluster.genes))
    if (length(cluster.genes)>250){
      cluster.genes = cluster.genes[1:250]
    }
    print(length(cluster.genes))
    overlap = intersect(CAD,names(cluster.genes))
    
    m=0
    c=0
    set.seed(123)
    for (j in 1:10000){
      random.genes = sample(row.names(raw.genecountsHG),length(cluster.genes),replace = T)
      random.overlap = intersect(CAD,random.genes)
      m[j]=length(random.overlap)
      if (length(overlap)>length(random.overlap)){
        c=c+1
      }
    }
    
    
    (hist.default(m, xlim = c(0,(2*max(m))), breaks = (max(m)+1),col = scales::hue_pal()(5)[i+1],main = paste("p = ",1-c/10000)))
    (points((length(overlap)),0, pch = 20, cex = 4, col = "red"))
    abline(v=(length(overlap)), col = "red",lwd = 4)
  }
  dev.off()
}

# save object
#saveRDS(average.expression, "seuset_Kmeans_object_average_expression_DEG.RDS")

###################
seuset <- readRDS("temp.seuset.RDS")
pdf(file = paste(a_name,"_Feature_plots_genes.pdf"),height = 3.5,width = 3.5)
FeaturePlot(object = seuset, features = "FGF13",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "ACTA2",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "HBB",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "MYH11",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "CXCL12",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "ACKR1",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "KLF4",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "APOE-ENSG00000130203",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "XIST-ENSG00000229807",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "ACTB-ENSG00000075624",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "MYH11-ENSG00000276480",cols = topo.colors(100))
FeaturePlot(object = colon, features = "STARD9-ENSG00000159433",cols = topo.colors(100))
FeaturePlot(object = colon, features = "MTATP8P2-ENSG00000229604",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CD74-ENSG00000019582",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CYSLTR1-ENSG00000173198",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "IFI30",cols = colorpanel(100,"grey90","darkgreen","black"))

FeaturePlot(object = colon, features = "CSF3R-ENSG00000119535",cols = topo.colors(100))
FeaturePlot(object = colon, features = "CYSLTR1-ENSG00000173198",cols = topo.colors(100))
FeaturePlot(object = colon, features = "KIT-ENSG00000157404",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "USP8",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = colon, features = "MKI67-ENSG00000148773",cols = topo.colors(100))
FeaturePlot(object = colonENS, features = "Neutrophil degranulation1",cols = topo.colors(100))
FeaturePlot(object = seuset, features = "C1QA",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "CD63",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "SOD2", cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "TYROBP",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "CTSB",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "ADAM32", cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "SLC35E3",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "KYNU",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seuset, features = "CTSB", cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = colon, features = "SPIC-ENSG00000166211",cols = topo.colors(100))

SLC35E3



dev.off()




DoHeatmap(seusetild_v3, features = c("CCDC144A","CYSLTR1","KYNU", "ADAM32", "CYP46A1",   "TMEM212","SOD2",   
                                     "PCSK5", "CALD1",   "MYH11"   ,  "DSTN" , "ACTA2"  ,    "LMOD1" , "MYH10", 
                                     "IGFBP7","USP8","ZNF286A","ATXN3", "RPS6KA5", "SLC35E3", "ORAI2" ,"VDR",    
                                     "CD81" ,   "C1QA"    ,  "CD63"  ,  "C1QB"   , "PLTP", "GRN",   
                                     "SRGN"   ,   "TYROBP"  ,  "CD74" ,"SAT1" ,
                                     "CTSB"   ,   "IFI30"  ,   "FTH1"    ,  "PSAP" , "FTL"    ,   "GPNMB"  , "TYROBP",
                                     "APOE"   ,   "APOC1"  ,  
                                     "SRGN"   ,   "PKM" ,     
                                     "IGLL5" ),) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))

seusetild_v3.cluster.averages <- AverageExpression(seusetild_v3, return.seurat = TRUE)
DoHeatmap(seusetild_v3.cluster.averages, features = c("CCDC144A","CYSLTR1","KYNU", "ADAM32", "CYP46A1",   "TMEM212","SOD2",
                                     "PCSK5", "CALD1",   "MYH11"   ,  "DSTN" , "ACTA2"  ,    "LMOD1" , "MYH10", 
                                     "IGFBP7","USP8","ZNF286A","ATXN3", "RPS6KA5", "SLC35E3", "ORAI2" ,"VDR",    
                                     "CD81" ,   "C1QA"    ,  "CD63"  ,  "C1QB"   , "PLTP", "GRN",   
                                     "SRGN"   ,   "TYROBP"  ,  "CD74" ,"SAT1" ,
                                     "CTSB"   ,   "IFI30"  ,   "FTH1"    ,  "PSAP" , "FTL"    ,   "GPNMB"  , "TYROBP",
                                     "APOE"   ,   "APOC1"  ,  
                                     "SRGN"   ,   "PKM" ,     
                                     "IGLL5" ),draw.lines = F) + scale_fill_gradientn(colors = c("grey100","grey90","darkgreen"))


##############
load("all.seur.combined.RData")
seurat_40 = readRDS("seurat_40_pts_20200416.RDS")
seusetild_v3 <- UpdateSeuratObject(object = all.seur.combined)
FeaturePlot(object = seusetild_v3, features = "FGF13",cols = colorpanel(100,"grey90","darkgreen","black"))
FeaturePlot(object = seusetild_v3, features = "KYNU",cols = colorpanel(100,"grey90","darkgreen","black"))


FeaturePlot(object = seusetild_v3, features = "FURIN",cols = colorpanel(100,"grey90","darkgreen","black"))


VlnPlot(object = seuset, features = "ACTA2")
VlnPlot(object = seuset, features = "MYH11")
VlnPlot(object = seuset, features = "MYH10")
VlnPlot(object = seuset, features = "SOD2")
VlnPlot(object = seuset, features = "IL6")
VlnPlot(object = seuset, features = "XBP1")
FeaturePlot(object = colonHG, features = "NOS1")
FeaturePlot(object = colonHG, features = "XBP1")
FeaturePlot(object = colonHG, features = "ACE2")
FeaturePlot(object = colonHG, features = "FURIN")
###############Y chrom
VlnPlot(object = colonHG, features = "ZFY")
VlnPlot(object = colonHG, features = "SRY")
VlnPlot(object = colonHG, features = "TSPY2")
VlnPlot(object = colonHG, features = "AMELY")

VlnPlot(object = colonHG, features = "SRY")
VlnPlot(object = colonHG, features = "UTY")
VlnPlot(object = colonHG, features = "PRY")

VlnPlot(object = colonHG, features = "ACE2")
VlnPlot(object = colonHG, features = "FURIN")


#save.image(file='whole_analysis.RData')
load('whole_analysis.RData')



a <- DotPlot(object = seusetild_v3, features = c("SVEP1"))
a$data


