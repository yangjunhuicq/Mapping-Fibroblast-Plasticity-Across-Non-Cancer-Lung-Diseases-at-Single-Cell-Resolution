args <- commandArgs(TRUE)
if(length(args)!=4)
{
        print("<word_dir> <input rds> <output rds> <groupby>")
        print("Author: Yang Junhui")
        q()
}
library(Seurat)
library(cowplot)
library(harmony)

setwd(args[1])
lung_fib.my.clean=readRDS(file = args[2])
dim(lung_fib.my.clean@meta.data)

lung_fib.my.clean


lung_fib.my.clean[["percent.mt"]] <- PercentageFeatureSet(lung_fib.my.clean, pattern = "^mt-|^MT-")

pdf("lung_fib.my.clean.QC.new_filter.clean.ori.pdf",width=12)
VlnPlot(lung_fib.my.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0,group.by=args[4])
plot1 <- FeatureScatter(lung_fib.my.clean, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by=args[4])
plot2 <- FeatureScatter(lung_fib.my.clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by=args[4])
print(plot1 + plot2)
dev.off()

lung_fib.my.clean <- subset(lung_fib.my.clean, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 15)
lung_fib.my.clean
pdf("lung_fib.my.clean.QC.new_filter.clean.pdf",width=12)
VlnPlot(lung_fib.my.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0,group.by=args[4])
plot1 <- FeatureScatter(lung_fib.my.clean, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by=args[4])
plot2 <- FeatureScatter(lung_fib.my.clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by=args[4])
print(plot1 + plot2)
dev.off()

lung_fib.my.clean <- NormalizeData(lung_fib.my.clean, normalization.method = "LogNormalize", scale.factor = 10000)

lung_fib.my.clean <- FindVariableFeatures(lung_fib.my.clean, selection.method = "vst", nfeatures = 2000)

lung_fib.my.clean


dim(lung_fib.my.clean[['RNA']]@scale.data)

all.genes <- rownames(lung_fib.my.clean)
lung_fib.my.clean=ScaleData(lung_fib.my.clean, features = all.genes)#"percent.mt" , vars.to.regress = "nFeature_RNA"
dim(lung_fib.my.clean[['RNA']]@scale.data)

lung_fib.my.clean=RunPCA(lung_fib.my.clean, features = VariableFeatures(lung_fib.my.clean))#npcs: Total Number of PCs to compute and store (50 by default)


pdf("lung_fib.my.clean.ElbowPlot.new_filter.clean.pdf",width=20)
ElbowPlot(lung_fib.my.clean)#by default, ndims=20, reduction="pca"
ElbowPlot(lung_fib.my.clean,ndims=50)
DimPlot(lung_fib.my.clean, reduction = "pca",group.by=args[4])
DimHeatmap(lung_fib.my.clean, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

lung_fib.my.clean <- RunHarmony(lung_fib.my.clean, "orig.ident") #by default, reduction = "pca"
#unique(lung_fib.my.clean@active.ident)

markers=read.table("markers2",header=F,row.names=1,sep="\t")
human_genes=intersect(rownames(markers),rownames(lung_fib.my.clean))
for(i in c(15,20,25,30))
{
  lung_fib.my.clean <- RunUMAP(lung_fib.my.clean, reduction = "harmony", dims = 1:i, reduction.name=paste("umapbyharmony",i,sep=""), reduction.key=paste("umapbyharmony",i,"_",sep="")) #by default, reduction = "pca"
  pdf(paste("lung_fib.my.clean.new_filter.clean.markers.FeaturePlot.","umapbyharmony",i,".new.pdf",sep=""),width=15)
  for(g in human_genes)
  {print(FeaturePlot(lung_fib.my.clean, features = g, raster=TRUE, reduction=paste("umapbyharmony",i,sep=""), pt.size=0.1))}
  dev.off()
}

lung_fib.my.clean <- FindNeighbors(lung_fib.my.clean, dims = 1:15, reduction = "harmony") #by default, reduction = "pca"
lung_fib.my.clean <- FindClusters(lung_fib.my.clean, resolution = 0.3) #会覆盖harmony之前的聚类结果。Note that 'seurat_clusters' will be overwritten everytime FindClusters is run

unique(lung_fib.my.clean@active.ident)
pdf("lung_fib.my.clean.ElbowPlot.Harmony.new_filter.clean.width.res0.3.pdf",width=15)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=T)
dev.off()

lung_fib.my.clean <- FindClusters(lung_fib.my.clean, resolution = 0.2)
unique(lung_fib.my.clean@active.ident)
pdf("lung_fib.my.clean.ElbowPlot.Harmony.new_filter.clean.width.res0.2.pdf",width=15)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=T)
dev.off()

lung_fib.my.clean <- FindClusters(lung_fib.my.clean, resolution = 0.1)
unique(lung_fib.my.clean@active.ident)
pdf("lung_fib.my.clean.ElbowPlot.Harmony.new_filter.clean.width.res0.1.pdf",width=15)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=T)
dev.off()

lung_fib.my.clean <- FindClusters(lung_fib.my.clean, resolution = 0.05)
unique(lung_fib.my.clean@active.ident)
pdf("lung_fib.my.clean.ElbowPlot.Harmony.new_filter.clean.width.res0.05.pdf",width=15)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=T)
dev.off()


saveRDS(lung_fib.my.clean, file = args[3])
lung_fib.my.clean@commands


lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="RNA_snn_res.0.3")
unique(lung_fib.my.clean@active.ident)
lung_fib.my.clean.markers <- FindAllMarkers(lung_fib.my.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #logfc.threshold 0.25: fc>1.189207 || fc<0.8408964. logfc.threshold 0.5: fc>1.414214 || fc<0.7071068.
write.table(lung_fib.my.clean.markers,file="lung_fib.my.clean.clustering_after_harmony.FindAllMarkers.logfc0.25.har15.res0.3.xls",col.names=T,row.names=T,quote=F,sep="\t")

lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="RNA_snn_res.0.2")
unique(lung_fib.my.clean@active.ident)
lung_fib.my.clean.markers <- FindAllMarkers(lung_fib.my.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #logfc.threshold 0.25: fc>1.189207 || fc<0.8408964. logfc.threshold 0.5: fc>1.414214 || fc<0.7071068.
write.table(lung_fib.my.clean.markers,file="lung_fib.my.clean.clustering_after_harmony.FindAllMarkers.logfc0.25.har15.res0.2.xls",col.names=T,row.names=T,quote=F,sep="\t")

lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="RNA_snn_res.0.1")
unique(lung_fib.my.clean@active.ident)
lung_fib.my.clean.markers <- FindAllMarkers(lung_fib.my.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #logfc.threshold 0.25: fc>1.189207 || fc<0.8408964. logfc.threshold 0.5: fc>1.414214 || fc<0.7071068.
write.table(lung_fib.my.clean.markers,file="lung_fib.my.clean.clustering_after_harmony.FindAllMarkers.logfc0.25.har15.res0.1.xls",col.names=T,row.names=T,quote=F,sep="\t")

lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="RNA_snn_res.0.05")
unique(lung_fib.my.clean@active.ident)
lung_fib.my.clean.markers <- FindAllMarkers(lung_fib.my.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #logfc.threshold 0.25: fc>1.189207 || fc<0.8408964. logfc.threshold 0.5: fc>1.414214 || fc<0.7071068.
write.table(lung_fib.my.clean.markers,file="lung_fib.my.clean.clustering_after_harmony.FindAllMarkers.logfc0.25.har15.res0.05.xls",col.names=T,row.names=T,quote=F,sep="\t")


