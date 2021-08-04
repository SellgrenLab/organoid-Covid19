library(GeneOverlap)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)

cov<- readRDS("./Covid_final_dataset.rds")
## Read in the microglial cluster signatures, Olah et al., 2020
dam<- readxl::read_xls("./Desktop/Microglia/MG_signatures.xls")


geneslist<- list()
dam.genes<- dam[which(dam$up_type==1),1]
geneslist[[1]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==2),1]
geneslist[[2]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==3),1]
geneslist[[3]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==4),1]
geneslist[[4]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==5),1]
geneslist[[5]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==6),1]
geneslist[[6]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==7),1]
geneslist[[7]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==8),1]
geneslist[[8]]<- as.data.frame(dam.genes)$gene
dam.genes<- dam[which(dam$up_type==9),1]
geneslist[[9]]<- as.data.frame(dam.genes)$gene

names(geneslist) <- c("Olah_MG_cluster1","Olah_MG_cluster2","Olah_MG_cluster3","Olah_MG_cluster4","Olah_MG_cluster5","Olah_MG_cluster6","Olah_MG_cluster7","Olah_MG_cluster8", "Olah_MG_cluster9")


microglia <- subset(cov, ident="Microglia")
Idents(microglia) <- "Condition"
DGE_cell_selection_72 <- FindMarkers(microglia, ident.1= "72 HPI", ident.2="Control" , logfc.threshold = -Inf,test.use = "wilcox", min.pct = 0,  min.diff.pct = 0, only.pos = FALSE,  max.cells.per.ident = 1000, return.thresh = 1,  assay = "RNA")
DGE_cell_selection_72 <- dplyr::filter(DGE_cell_selection_72, p_val_adj<0.05)
mg.genes <- dplyr::filter(DGE_cell_selection_72, avg_log2FC>0)

x=1
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=2
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=3
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=4
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=5
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=6
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=7
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=8
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])
x=9
goobj.list[[x]] <- newGeneOverlap(rownames(mg.genes), geneslist[[x]], genome.size = 8000)
goobj.list[[x]] <- testGeneOverlap(goobj.list[[x]])
print(goobj.list[[x]])






