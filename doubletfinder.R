
print("##Start:sc-rna-8-doubletfinder")
date()

library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

## load the data
data <- readRDS("GSM5133603.ana.rds")

## pK Identification (no ground-truth)
pdf("GSM5133603.df-sc-rna-8-doubletfinder.pdf")
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line() + theme(axis.text.x = element_text( angle = 90, hjust=1,vjust=1 ))
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK) 
pK <- as.numeric(as.character(pK[[1]]))
cat("pK selected is :",pK,"\n")

## Homotypic Doublet Proportion Estimate
annotations <- data@meta.data[,"RNA_snn_res.0.5"]
homotypic.prop <- modelHomotypic(annotations) 

## Exp doublet number calculation
cell_num <- nrow(data@meta.data)
print("Cell number is :")
cat(cell_num,"\n")

multiplet_rate <- cell_num * 0.0008 + 0.0527
print("Multiplet rate is :")
cat(multiplet_rate,"%\n")

nExp_poi <- round(cell_num * multiplet_rate * 0.01)
print("Expected doublet number is :")
cat(nExp_poi,"\n")

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
print("Expected doublet number adjusted by homotypicc doublet is :")
cat(nExp_poi.adj,"\n")

## Run DoubletFinder with varying classification stringencies 
data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

DimPlot(data, reduction = "umap", group.by= paste0("DF.classifications_0.25_",pK,"_",nExp_poi.adj))
dev.off()

metadata <- data@meta.data
write.table(metadata,"GSM5133603.df-metadata_doubletfinder.txt",quote=F,sep="\t",col.names=NA)

print("##Finished:sc-rna-8-doubletfinder")
date()

