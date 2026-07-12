args <- commandArgs(TRUE)
if(length(args)!=3)
{
        print("<word_dir> <input file: genelist> <output prefix>")
        print("Author: Yang Junhui")
        q()
}

library(msigdbr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

setwd(args[1])
inp=read.table(args[2],header=T,row.names=NULL,sep="\t")
genelist <- unique(inp[,1])

#KEGG
genesets = msigdbr(species = "Homo sapiens", category = "C2")
# View(msigdbr_collections()) #查看msigdbr包中所有的基因集
print(unique(genesets$gs_subcat))  # 有多个数据库来源的基因集可选，这里选用KEGG
genesetsuse <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse$gs_name))) #查看有多少条通路（186个）
print(length(intersect(genelist,unique(genesetsuse$gene_symbol)))) #查看有多少个靶基因有通路信息

res <- enricher(genelist, TERM2GENE = genesetsuse, pvalueCutoff=1, qvalueCutoff = 1)
saveRDS(res,file=paste0(args[3],".enricher_KEGG.rds"))
res <- data.frame(res)
res$set_size=rep(-1,nrow(res))
res$set_in=rep(-1,nrow(res))
res$set_fraction=rep(-1,nrow(res))
for(pat in rownames(res))
{
 setall=unique(genesetsuse$gene_symbol[which(genesetsuse$gs_name==pat)])
 setin=intersect(setall,genelist)
 setfrac=length(setin)/length(setall)
 res[pat,"set_size"]=length(setall)
 res[pat,"set_in"]=length(setin)
 res[pat,"set_fraction"]=setfrac
}
write.table(res, file=paste0(args[3],".enricher_KEGG.xls"), row.names = F, sep="\t", quote=F)


#CP:REACTOME
genesetsuse <- subset(genesets, gs_subcat=="CP:REACTOME", select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse$gs_name))) #查看有多少条通路（186个）
print(length(intersect(genelist,unique(genesetsuse$gene_symbol))))

res <- enricher(genelist, TERM2GENE = genesetsuse, pvalueCutoff=1, qvalueCutoff = 1)
saveRDS(res,file=paste0(args[3],".enricher_REACTOME.rds"))
res <- data.frame(res)
write.table(res, file=paste0(args[3],".enricher_REACTOME.xls"), row.names = F, sep="\t", quote=F)


#CP:BIOCARTA
genesetsuse <- subset(genesets, gs_subcat=="CP:BIOCARTA", select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse$gs_name))) #查看有多少条通路（186个）
print(length(intersect(genelist,unique(genesetsuse$gene_symbol))))

res <- enricher(genelist, TERM2GENE = genesetsuse, pvalueCutoff=1, qvalueCutoff = 1)
saveRDS(res,file=paste0(args[3],".enricher_BIOCARTA.rds"))
res <- data.frame(res)
write.table(res, file=paste0(args[3],".enricher_BIOCARTA.xls"), row.names = F, sep="\t", quote=F)


#H
genesetsH = msigdbr(species = "Homo sapiens", category = "H")
genesetsuse <- subset(genesetsH, select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse$gs_name))) #查看有多少条通路（186个）
print(length(intersect(genelist,unique(genesetsuse$gene_symbol))))

res <- enricher(genelist, TERM2GENE = genesetsuse, pvalueCutoff=1, qvalueCutoff = 1)
saveRDS(res,file=paste0(args[3],".enricher_Hallmark.rds"))
res <- data.frame(res)
write.table(res, file=paste0(args[3],".enricher_Hallmark.xls"), row.names = F, sep="\t", quote=F)


#GO
genesetsGO = msigdbr(species = "Homo sapiens", category = "C5")
print(unique(genesetsGO$gs_subcat))
#[1] "GO:BP" "GO:CC" "GO:MF" "HPO"

#GO:BP
genesetsuse <- subset(genesetsGO, gs_subcat %in% c("GO:BP","GO:CC","GO:MF"), select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse$gs_name))) #查看有多少条通路（186个）
print(length(intersect(genelist,unique(genesetsuse$gene_symbol))))

res <- enricher(genelist, TERM2GENE = genesetsuse, pvalueCutoff=1, qvalueCutoff = 1)
saveRDS(res,file=paste0(args[3],".enricher_GO.rds"))
res <- data.frame(res)
write.table(res, file=paste0(args[3],".enricher_GO.xls"), row.names = F, sep="\t", quote=F)


