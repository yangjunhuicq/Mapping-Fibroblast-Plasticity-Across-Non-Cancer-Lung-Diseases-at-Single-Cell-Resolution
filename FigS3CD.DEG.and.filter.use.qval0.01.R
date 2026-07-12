args <- commandArgs(TRUE)
if(length(args)<5)
{
        print("<word_dir> <input rds> <group to call DEG> <output prefix> <group1> <group2>")
        print("Author: Yang Junhui")
        q()
}

library(Seurat)
library(cowplot)
library(harmony)

setwd(args[1])
lung_fib.my.clean1=readRDS(file = args[2])
lung_fib.my.clean1

print("the last resolution:")
unique(lung_fib.my.clean1@active.ident)

print("");print("resolution selected:")
lung_fib.my.clean1 <- SetIdent(lung_fib.my.clean1, value=args[3])
table(lung_fib.my.clean1@active.ident)

group1=strsplit(args[5], ",", fixed = TRUE)[[1]]
group2=setdiff(unique(lung_fib.my.clean1@active.ident),args[5])
if(length(args)>5)
{
 group2=strsplit(args[6], ",", fixed = TRUE)[[1]]
}
print(group1)
print(group2)
COPD.markers <- FindMarkers(lung_fib.my.clean1, ident.1 = group1, ident.2 = group2, min.pct = 0.25)
write.table(COPD.markers,file=paste0(args[4],".DEG.xls"),col.names=NA,row.names=T,quote=F,sep="\t")
cat("DEG num: ")
dim(COPD.markers)

COPD.markers=read.table(paste0(args[4],".DEG.xls"),header=T,row.names=1,sep="\t")
perc=rownames(COPD.markers)[which(COPD.markers$p_val_adj<=0.01 & abs(COPD.markers$pct.1-COPD.markers$pct.2)>=0.2)]
log2FC=rownames(COPD.markers)[which(COPD.markers$p_val_adj<=0.01 & abs(COPD.markers$avg_log2FC)>=0.4)]
COPD.markers.fil=COPD.markers[union(perc,log2FC),]
print(dim(COPD.markers.fil))
write.table(COPD.markers.fil,file=paste0(args[4],".DEG.filt.qval0.01.xls"),col.names=NA,row.names=T,quote=F,sep="\t")
if(length(grep('^RPS|^RPL|^MT-|^MTRNR',rownames(COPD.markers.fil),ignore.case=T))>0)
{
COPD.markersuse=rownames(COPD.markers.fil)[-grep('^RPS|^RPL|^MT-|^MTRNR',rownames(COPD.markers.fil),ignore.case=T)]
COPD.markers.fil2=COPD.markers.fil[COPD.markersuse,]
print(dim(COPD.markers.fil2))
write.table(COPD.markers.fil2,file=paste0(args[4],".DEG.filt.qval0.01.noRP_noMT.xls"),col.names=NA,row.names=T,quote=F,sep="\t")
}

