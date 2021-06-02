library(Seurat)
library(dplyr)
library(Matrix)
library(SingleCellExperiment)
library(magrittr)
library(tidyverse)
#library(scran)
library(scDblFinder)
#library(SingleR)

## Load Cellranger output files
data_dir <- '/home/sellgren/covid/main/feature_matrices/001/filtered_feature_bc_matrix'
list.files(data_dir)
expression_matrix <- Read10X(data.dir = data_dir)
data_dir <- '/home/sellgren/covid/main/feature_matrices/002/filtered_feature_bc_matrix'
expression_matrix2 <- Read10X(data.dir = data_dir)
data_dir <- '/home/sellgren/covid/main/feature_matrices/003/filtered_feature_bc_matrix'
expression_matrix3 <- Read10X(data.dir = data_dir)
data_dir <- '/home/sellgren/covid/main/feature_matrices/004/filtered_feature_bc_matrix'
expression_matrix4 <- Read10X(data.dir = data_dir)
data_dir <- '/home/sellgren/covid/main/feature_matrices/005/filtered_feature_bc_matrix'
expression_matrix5 <- Read10X(data.dir = data_dir)
data_dir <- '/home/sellgren/covid/main/feature_matrices/006/filtered_feature_bc_matrix'
expression_matrix6 <- Read10X(data.dir = data_dir)


cms <- list()
cms<- list(expression_matrix, expression_matrix2,expression_matrix3,expression_matrix4,expression_matrix5,expression_matrix6)
str(cms)
names(cms) <- c("10X_bsl3_001", "10X_bsl3_002", "10X_bsl3_003","10X_bsl3_004", "10X_bsl3_005","10X_bsl3_006")
str(cms)


sameGeneNames <- function(counts) {
  ## create place holder for gene names from each count matrix
  gene_names <- list()
  ## go through count matrices
  for ( i in names(counts) ) {
    ## extract gene names
    gene_names[[i]] <- rownames(counts[[i]])
  }
  ## check if all gene names are the same (using the first matrix as a reference
  ## for the others)
  return(all(sapply(gene_names, FUN = identical, gene_names[[1]])))
}

sameGeneNames(cms)

x=1
colnames(cms[[x]]) <- paste("10X_bsl3_001", colnames(cms[[x]]) , sep="_")
x=2
colnames(cms[[x]]) <- paste("10X_bsl3_002", colnames(cms[[x]]) , sep="_")
x=3
colnames(cms[[x]]) <- paste("10X_bsl3_003", colnames(cms[[x]]) , sep="_")
x=4
colnames(cms[[x]]) <- paste("10X_bsl3_004", colnames(cms[[x]]) , sep="_")
x=5
colnames(cms[[x]]) <- paste("10X_bsl3_005", colnames(cms[[x]]) , sep="_")
x=6
colnames(cms[[x]]) <- paste("10X_bsl3_006", colnames(cms[[x]]) , sep="_")

experiment.data <- do.call("cbind", cms)
dim(experiment.data)
experiment.aggregate <- CreateSeuratObject(experiment.data, project="Organoids Covid19", min.cells=3, min.features=200)
experiment.aggregate
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^MT-")
slotNames(experiment.aggregate)
experiment.aggregate$Chemistry <- "10X Chromium V3"
head(experiment.aggregate[[]])
experiment.aggregate$Species <- "Human"

## Set up meta data of all cells
meta<- rownames(experiment.aggregate@meta.data)
dim(meta)
class(meta)
meta<- as.data.frame(meta)
#head(meta)
control <- c(colnames(cms[[1]]), colnames(cms[[2]]))
length(control)
twentyfourhpi <- c(colnames(cms[[4]]), colnames(cms[[3]]))
seventytwohpi <- c(colnames(cms[[5]]), colnames(cms[[6]]))
meta$Sample <- ifelse(meta$meta%in%control, "Control", ifelse(meta$meta%in%twentyfourhpi, "24 HPI", ifelse(meta$meta%in%seventytwohpi, "72 HPI", "NA")))
table(meta$Sample)
#head(meta)
rownames(meta) <- meta$meta
meta<- select(meta, -meta)

## Identify doublets
experiment.aggregate$Condition <- meta
sce<- as.SingleCellExperiment(experiment.aggregate)
sce <- scDblFinder(sce)
experiment.aggregate$scDblFinder.class <- colData(sce)$scDblFinder.class

#Get cell IDs
lane1 <- colnames(cms[[1]])
lane2<- colnames(cms[[2]])
lane3 <- colnames(cms[[3]])
lane4 <- colnames(cms[[4]])
lane5 <- colnames(cms[[5]])
lane6<- colnames(cms[[6]])
meta$ControllerLane <- ifelse(rownames(meta)%in%lane1, "Lane 1", ifelse(rownames(meta)%in%lane2, "Lane 2", 
                                                                        ifelse(rownames(meta)%in%lane3, "Lane 3",
                                                                          ifelse(rownames(meta)%in%lane4, "Lane 4",  
                                                                                 ifelse(rownames(meta)%in%lane5, "Lane 5", ifelse(rownames(meta)%in%lane6, "Lane 6", "NA"))))))
meta<- select(meta, -Sample)
experiment.aggregate@meta.data <- meta
## Save the seurat object and count matrices
saveRDS(experiment.aggregate, "Covid_PreprocessedSeuratobj.rds")                                                                                                                             
saveRDS(cms, "Covid_CountMatrices.rds")

session_info()
                                                                                                                                