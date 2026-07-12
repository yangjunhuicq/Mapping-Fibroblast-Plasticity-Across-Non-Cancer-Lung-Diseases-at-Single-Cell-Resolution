library(SCP)
library(ggplot2)
library(Seurat)
inp=readRDS('inp.rds')
inp
colnames(inp@meta.data)

table(inp$Disease,useNA='always')
table(inp$celltypeuse,useNA='always')
inp$Disease=factor(inp$Disease,levels=c('Donor','COPD','PF','COVID19'))
inp$celltypeuse=factor(inp$celltypeuse,levels=c('AT1', 'AT2', 'Basal', 'Secretory', 'Ciliated', 'CD8', 'NK', 'CD4', 'B_Plasma', 'Macrophage', 'Dividing', 'Monocyte', 'DC', 'Mast', 'Endothelial', 'SMC_Pericytes', 'Fibroblast', 'Mesothelial'))
table(inp$Disease,useNA='always')
table(inp$celltypeuse,useNA='always')


inpx=read.table('genelist.xls',header=T,row.names=NULL,sep="\t",fill=T)
head(inpx)

for(i in seq(1,ncol(inpx)))
{
 genes=inpx[,i]
 genes=intersect(rownames(inp),genes)
 print(length(genes))
 inp=AddModuleScore(inp,features=list(unique(intersect(rownames(inp),genes))),name=c(paste(colnames(inpx)[i],sep="")))
 print(colnames(inp@meta.data)[ncol(inp@meta.data)])
}
print(colnames(inp@meta.data))



pdf(paste0("output.score.FeatureStatPlot.pdf"),width=6.5,height=3)
for(pathw in colnames(inpx))
{
print(FeatureStatPlot(inp, stat.by = paste0(pathw,"1"), group.by = 'celltypeuse', split.by="Disease", bg.by = 'celltypeuse', add_box = T, box_ptsize = 1, y.min=-0.6, y.max=2.5)+ theme(plot.title = element_text(size = 10), axis.title=element_text(colour = "black", size=7, face='bold'), legend.title = element_text(colour = "black", size=6), legend.text=element_text(colour = "black", size=6), axis.text = element_text(size=6)))
}
dev.off()

write.table(inp@meta.data,file=paste0('inp.rds.add_score2.metadata.xls'),row.names=T,col.names=NA,quote=F,sep='\t')

