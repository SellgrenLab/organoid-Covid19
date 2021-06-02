### Astrocyte Subclustering

library(Seurat)
library(dplyr)


cov<- readRDS("./Covid_final_dataset.rds")

Idents(cov) <- "celltype"
astro <- subset(cov, idents="Astrocytes")
astro<- CreateSeuratObject(counts = astro[["RNA"]]@counts, meta.data = astro@meta.data)
#astro<- NormalizeData(astro)
astro<- SCTransform(astro, vars.to.regress = c("CC.Difference"), verbose = T)
astro<- RunPCA(astro)
astro<- RunUMAP(astro, dims = 1:30)
astro <- FindNeighbors(astro, dims = 1:30, verbose = FALSE)
astro <- FindClusters(astro, verbose = FALSE)
DimPlot(astro, label = TRUE) + NoLegend()

# Astro phagocytic machinery
VlnPlot(astro, features = c("MEGF10", "MERTK", "TYRO3", "AXL", "GAS6", "PROS1", "GULP1", "ABCA1"))


### Comparison to reactive A1 astrocytes
## A1 astrocyte single-cell data extracted from Barbar et al., 2020 https://doi.org/10.1016/j.neuron.2020.05.014 

A1<- Read10X("./data/A1/")

