args <- commandArgs(TRUE)
#if(length(args)!=7)
#{
#        print("<word_dir> <input rds> <output rds> <resolution such as \"RNA_snn_res.0.3\"> <genelist_file> <cluster id for the Monocle2 root> <genes_to_show_withpseudotime>")
#        print("For the first time, the root be set to FALSE. For the second time, the root be set to a cluster id.")
#        print("Author: Yang Junhui")
#        q()
#}


setwd(args[1])
library(Seurat)
library(cowplot)
library(harmony)
library(monocle)

lung_fib.my.clean=readRDS(file = args[2])
lung_fib.my.clean
unique(lung_fib.my.clean@active.ident)

lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value=args[4])
unique(lung_fib.my.clean@active.ident)

print("");print("DefaultAssay")
DefaultAssay(lung_fib.my.clean)
#[1] "RNA"

metadata <- lung_fib.my.clean@meta.data
write.table(metadata,paste0(args[3],".metadata_ori.txt"),quote=F,sep="\t", row.names=TRUE, col.names=NA)


newimport <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts

    #if(class(data) == "data.frame") {
      #data <- as(as.matrix(data), "sparseMatrix")
    #}

    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)

      message("This Seurat object doesn't provide any meta data");
      pd
    })

    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }

    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0

    #if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    #} else if(any(data < 0)){
      #expressionFamily <- uninormal()
    #} else {
     # expressionFamily <- tobit()
    #}

    valid_data <- data[, row.names(pd)]

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)

    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL

        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS

      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }

    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)

    }
    monocle_cds@auxClusteringData$seurat <- mist_list

  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")

    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset

    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }

    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 

    } else {
      mist_list <- list()
    }

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list

    monocle_cds@auxOrderingData$scran <- mist_list

  } else {
    stop('the object type you want to export to is not supported yet')
  }

  return(monocle_cds)
}

print("import into monocle")
lung_fib.my.clean.mono=newimport(lung_fib.my.clean)
class(lung_fib.my.clean.mono)
dim(lung_fib.my.clean.mono)
length(lung_fib.my.clean.mono$Size_Factor)
length(unique(lung_fib.my.clean.mono$Size_Factor))
unique(lung_fib.my.clean.mono$Size_Factor)

cat("dim for lung_fib.my.clean.mono@phenoData: ")
dim(lung_fib.my.clean.mono@phenoData)
cat("\ncolnames for lung_fib.my.clean.mono@phenoData: ")
colnames(lung_fib.my.clean.mono@phenoData)
cat("\n")

print("estimateSizeFactors")
lung_fib.my.clean.mono <- estimateSizeFactors(lung_fib.my.clean.mono) # Size factors help us normalize for differences in mRNA recovered across cells
head(sort(lung_fib.my.clean.mono$Size_Factor))
tail(sort(lung_fib.my.clean.mono$Size_Factor))

print("estimateDispersions")
lung_fib.my.clean.mono <- estimateDispersions(lung_fib.my.clean.mono) #"dispersion" values will help us perform differential expression analysis later. 
dim(lung_fib.my.clean.mono)

#CellDataSet objects provide a convenient place to store per-cell scoring data: the phenoData slot. 
#Once you've excluded cells that do not pass your quality control filters, you should verify that the expression values stored in your CellDataSet follow a distribution that is roughly lognormal: 
# Log-transform each value in the expression matrix.
print("For qplot:")
L <- log(exprs(lung_fib.my.clean.mono))
dim(L)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
library(reshape2)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
dim(melted_dens_df)

# Plot the distribution of the standardized gene expression values.
pdf(paste0(args[3],".qplot.pdf"))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")
dev.off()

#######################Constructing Single Cell Trajectories
#all(lung_fib.my.clean$RNA_snn_res.0.3==Idents(lung_fib.my.clean))
#[1] TRUE

lung_fib.my.clean.markers=read.table(args[5],header=T,row.names=1,stringsAsFactors=F,sep="\t")
#lung_fib.my.clean.markers.sig=subset(lung_fib.my.clean.markers, avg_log2FC>=0.5 & p_val_adj<0.01)
sel.gene <- rownames(lung_fib.my.clean.markers)

if(ncol(lung_fib.my.clean.mono)>40000){
print("There are more than 4w cells, random select 4w")
set.seed(1);cellids=sample(1:ncol(lung_fib.my.clean.mono),size=40000)

lung_fib.my.clean.mono.40k=lung_fib.my.clean.mono[,cellids]
print(dim(exprs(lung_fib.my.clean.mono.40k)))
print(dim(pData(lung_fib.my.clean.mono.40k)))
print(dim(fData(lung_fib.my.clean.mono.40k)))
}else{
print("There are less than 4w cells, use all cells")
lung_fib.my.clean.mono.40k=lung_fib.my.clean.mono
}
lung_fib.my.clean.mono.40k <- monocle::setOrderingFilter(lung_fib.my.clean.mono.40k, sel.gene)
pdf(paste0(args[3],".step1.plot_ordering_genes.pdf"))
plot_ordering_genes(lung_fib.my.clean.mono.40k)
dev.off()
## dimension reduciton
lung_fib.my.clean.mono.40k <- monocle::reduceDimension(lung_fib.my.clean.mono.40k, method = 'DDRTree')

## ordering cells
lung_fib.my.clean.mono.40k <- monocle::orderCells(lung_fib.my.clean.mono.40k)  #need to run 2h

pdf(paste0(args[3],".plot_cell_trajectory.pdf"))
plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = args[4])
plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "State")
plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "Pseudotime")
dev.off()
num=length(unique(lung_fib.my.clean.mono.40k@phenoData@data[,args[4]])) # lung_fib.my.clean.mono.40k@phenoData@data is same as pData(lung_fib.my.clean.mono.40k)
if(num<length(unique(lung_fib.my.clean.mono.40k@phenoData@data[,"State"])))
{
num=length(unique(lung_fib.my.clean.mono.40k@phenoData@data[,"State"]))
}
lung_fib.my.clean.mono.40k$forplot=lung_fib.my.clean.mono.40k@phenoData@data[,args[4]]
pdf(paste0(args[3],".plot_cell_trajectory.split.pdf"),width=5*num)
plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = args[4]) + facet_wrap(~forplot, nrow = 1)
plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "State") + facet_wrap(~State, nrow = 1)
dev.off()

saveRDS(lung_fib.my.clean.mono.40k,file=paste0(args[3],".mono.rds"))
#Visualize gene expression
if(length(args)==7){
human_genes=intersect(rownames(read.table(args[7],header=T,row.names=1,sep="\t")),rownames(lung_fib.my.clean))
num=length(human_genes)
x=floor(num/5)-1
y=num-(x+1)*5

pdf(paste0(args[3],".plot_genes_in_pseudotime.pdf"),width=20,height=5*(x+2))
if(x>=0){
for(i in c(0:x)){
print(human_genes[(i*5+1)])
lung_fib.my.clean.mono.40k_subset <- lung_fib.my.clean.mono.40k[human_genes[(i*5+1):(i*5+5)],]
print(plot_genes_in_pseudotime(lung_fib.my.clean.mono.40k_subset, color_by = args[4]))
}
}
if(y>0){
print(warnings())
lung_fib.my.clean.mono.40k_subset <- lung_fib.my.clean.mono.40k[human_genes[(i*5+6):(i*5+5+y)],]
print(plot_genes_in_pseudotime(lung_fib.my.clean.mono.40k_subset, color_by = args[4]))
dev.off()
print(warnings())
}
}

if(args[6] != "FALSE")
{
## ordering cells by assigning root nodes
GM_state <- function(cds){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds@phenoData@data[,args[4]])[,args[6]]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cat("The selectd root State is:")
print(GM_state(lung_fib.my.clean.mono.40k))
lung_fib.my.clean.mono.40k <- monocle::orderCells(lung_fib.my.clean.mono.40k, root_state =  GM_state(lung_fib.my.clean.mono.40k))

pdf(paste0(args[3],".plot_cell_trajectory.rooted.pdf"))
print(plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = args[4]))
print(plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "State"))
print(plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "Pseudotime"))
dev.off()
print(warnings())
pdf(paste0(args[3],".plot_cell_trajectory.rooted.split.pdf"),width=5*num)
print(plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = args[4]) + facet_wrap(~forplot, nrow = 1))
print(plot_cell_trajectory(lung_fib.my.clean.mono.40k, color_by = "State") + facet_wrap(~State, nrow = 1))
dev.off()
print(warnings())

#Visualize gene expression
if(length(args)==7){
human_genes=intersect(rownames(read.table(args[7],header=T,row.names=1,sep="\t")),rownames(lung_fib.my.clean))
num=length(human_genes)
x=floor(num/5)-1
y=num-(x+1)*5

pdf(paste0(args[3],".plot_genes_in_pseudotime.40k.rooted.pdf"),width=20,height=5*(x+2))
if(x>0){
for(i in c(0:x)){
lung_fib.my.clean.mono.40k_subset <- lung_fib.my.clean.mono.40k[human_genes[(i*5+1):(i*5+5)],]
print(plot_genes_in_pseudotime(lung_fib.my.clean.mono.40k_subset, color_by = args[4]))
}
}
if(y>0){
print(warnings())
lung_fib.my.clean.mono.40k_subset <- lung_fib.my.clean.mono.40k[human_genes[(i*5+6):(i*5+5+y)],]
print(plot_genes_in_pseudotime(lung_fib.my.clean.mono.40k_subset, color_by = args[4]))
dev.off()
print(warnings())
}
}
}

saveRDS(lung_fib.my.clean.mono.40k,file=paste0(args[3],".mono.rds"))
all(lung_fib.my.clean@meta.data[rownames(lung_fib.my.clean.mono.40k@phenoData@data),]==lung_fib.my.clean.mono.40k@phenoData@data[,1:ncol(lung_fib.my.clean@meta.data)],na.rm=T)

dim(lung_fib.my.clean@meta.data)
print("add Pseudotime and State to the metadata for the input rds")
lung_fib.my.clean$Pseudotime=rep("",nrow(lung_fib.my.clean@meta.data))
lung_fib.my.clean$State=rep("",nrow(lung_fib.my.clean@meta.data))
for(n in rownames(lung_fib.my.clean@meta.data))
{
 if(n %in% rownames(lung_fib.my.clean.mono.40k@phenoData@data))
 {
  lung_fib.my.clean@meta.data[n,"Pseudotime"]=lung_fib.my.clean.mono.40k@phenoData@data[n,"Pseudotime"]
  lung_fib.my.clean@meta.data[n,"State"]=lung_fib.my.clean.mono.40k@phenoData@data[n,"State"]
 }
}
dim(lung_fib.my.clean@meta.data)
metadata <- lung_fib.my.clean@meta.data
write.table(metadata,paste0(args[3],".metadata_add.txt"),quote=F,sep="\t", row.names=TRUE, col.names=NA)
saveRDS(lung_fib.my.clean, file = args[3])
