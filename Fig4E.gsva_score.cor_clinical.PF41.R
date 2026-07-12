library(GSVA)

out=read.table('GSE47460_GPL14550_429samples.exp.txt.addGenesymbol',header=T,row.names=NULL,sep='\t')
rownames(out)=make.unique(out[,ncol(out)])
out=out[,-c(1,ncol(out))]
head(out[,c(1:3)])
dim(out)

FPKM2=as.matrix(out)
head(FPKM2[,c(1:3)])
dim(FPKM2)

genesets=read.table("PF41.xls",header=T,sep="\t",row.names=NULL)
head(genesets)
genesets_list <- list(PF53 = as.character(genesets[, 1]))
genesets_list

result=gsva(FPKM2,genesets_list,method="gsva",verbose=TRUE)
head(result[,c(1:3)])
dim(result)
write.table(result,file="GSE47460_GPL14550_429samples.gsva_PF41.xls",col.names =NA,row.names = T,quote=F,sep="\t")



library(ggplot2)
library(ggpubr)
out=read.table('GSE47460_GPL14550_429samples.gsva_PF41.xls',header=T,row.names=1,sep='\t')
dim(out)


clinical=read.table('clinicalinfo.xls',header=T,row.names=2,sep='\t')
dim(clinical)
head(clinical[,c(1:3)])
length(intersect(colnames(out),row.names(clinical)))
samples=intersect(colnames(out),row.names(clinical))
clinical=clinical[samples,]
dim(clinical)
Control=rownames(clinical)[which(clinical$disease=='Control')]
COPD=rownames(clinical)[which(clinical$disease=='Chronic Obstructive Lung Disease')]
IPF=rownames(clinical)[which(clinical$disease=='Interstitial lung disease')]
length(Control)
length(COPD)
length(IPF)

clinicaluse=clinical[which(clinical$disease=='Control' | (clinical$disease=='Interstitial lung disease' & clinical$ild_subtype=='2-UIP/IPF')),]
dim(clinicaluse)

clinicaluse=clinical[which(clinical$disease=='Control' | (clinical$disease=='Interstitial lung disease' & clinical$ild_subtype=='2-UIP/IPF')),]
dim(clinicaluse)
for(i in c("age","emphysema","predicted_fev1_prebd","predicted_fvc_prebd","predicted_fev1_postbd","predicted_fvc_postbd","predicted_dlco"))
{
 print(i)
 print(table(clinicaluse[,i],useNA='always'))
 clinicaluse2=clinicaluse[which(clinicaluse[,i]!=''),]
 print(table(clinicaluse2[,i],useNA='always'))
 print(dim(clinicaluse2))
 test_out=c()
 pdf(paste0('gsva_clinical.scatter',i,'.PF41.pdf'))
 for(g in rownames(out))
 {
  data = data.frame(Expression = t(out[g,rownames(clinicaluse2)]), group = clinicaluse2[,i], disease=clinicaluse2[,'disease'])
  print(ggscatter(data, x = "group", y = g, add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = i, ylab = g))
  test_out=c(test_out,g,cor.test(data[,g],data[,'group'])$p.value,cor.test(data[,g],data[,'group'])$estimate[[1]],cor.test(data[,g],data[,'group'],method="spearman")$p.value,cor.test(data[,g],data[,'group'],method="spearman")$estimate[[1]])
 }
 print(dev.off())
 test_out_mat=matrix(test_out,ncol=5,byrow=T)
 colnames(test_out_mat)=c('Gene','pearson.pval','pearson.cor','spearman.pval','spearman.cor')
 write.table(test_out_mat,file=paste0('gsva_clinical.',i,'.result.PF41.xls'),row.names=F,col.names=T,quote=F,sep='\t')
}

