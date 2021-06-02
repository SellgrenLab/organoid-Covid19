## Second round of clustering

## Clusters for SARS-CoV-2 Analysis

library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Nebulosa)
library(clustifyr)

## Subset to include clusters that passed correlation and marker expression inspection (manual)
sars <- readRDS("./data/Covid_clustered_1.rds")
sars_subset <- SubsetData(sars, ident.use = c("Newborn Neurons (Inhibitory)","Excitatory Neurons","Astro-Radial Glia","Oligodendroglia", 
                                              "Proliferative Neuroepithelium",
                                              "Unknown 1 (Cell cycle arrest)",
                                              "Proliferative Neural Precursors",
                                              "Choroid",
                                              "Pericytes",
                                              "Microglia","Endothelial", "Epithelial ACE2-expressing",
                                              "Vascular Smooth Muscle Cells (vSMCs)"))

## Perform Integration and recluster upon SCtransformation
options(future.globals.maxSize = 4000 * 1024^2)

meta<- sars_subset@meta.data
meta<- dplyr::select(meta, c("Chemistry", "percent.mito", "Tissue", "Species", 
                             "Condition", "ControllerLane", "S.Score", "G2M.Score",
                             "Phase", "CC.Difference"))
covid<- CreateSeuratObject(counts = sars_subset[["RNA"]]@counts, meta.data = meta)

#Check number of cells per condition
meta$Condition <- factor(meta$Condition, levels = c("Control", "24 HPI", "72 HPI")) # arrange

meta %>% 
  ggplot(aes(x=Condition, fill=Condition)) + 
  geom_bar() +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Condition")


organoid.list <- SplitObject(covid, split.by = "Condition")
for (i in 1:length(organoid.list)) {
  organoid.list[[i]] <- SCTransform(organoid.list[[i]], verbose = FALSE, vars.to.regress = c("percent.mito", "CC.Difference"))
}
organoid.features <- SelectIntegrationFeatures(object.list = organoid.list, nfeatures = 3000)
organoid.list <- PrepSCTIntegration(object.list = organoid.list, anchor.features = organoid.features, 
                                    verbose = FALSE)

organoid.anchors <- FindIntegrationAnchors(object.list = organoid.list, normalization.method = "SCT", 
                                           anchor.features = organoid.features, verbose = FALSE)

## Above warning for future lapply can be ignored. Should not be in Seurat v4. 

organoid.integrated <- IntegrateData(anchorset = organoid.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

organoid.integrated <- RunPCA(organoid.integrated, verbose = FALSE)
ElbowPlot(organoid.integrated,ndims = 40) #check how many pcs to use
covid.integrated <- RunUMAP(organoid.integrated, dims = 1:40)

## Visualize the new clusters and reannotate if necessary

plots <- DimPlot(covid.integrated, group.by = c("Condition", "Celltype_temp"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 4, byrow = TRUE, 
                                                                     override.aes = list(size = 2.5)))
## Recluster at different resolutions
covid.integrated <- FindNeighbors(covid.integrated,  dims = 1:40)
covid.integrated <- FindClusters(covid.integrated, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))



## Save the seurat object
saveRDS(covid.integrated, "./data/Covid_Integrated.rds")

## proceeding with resolution 0.6
Idents(object = covid.integrated) <- "integrated_snn_res.0.6"
DimPlot(covid.integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
## Check for uninteresting sources of variation in our clustering
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mito")
FeaturePlot(covid.integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order=TRUE,
            min.cutoff = 'q10',
            label = TRUE)
## Find the conserved markers across conditions and rename idents

# function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(combined,
                       ident.1 = cluster,
                       grouping.var = "Condition",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
n= 20 # number of clusters found
conserved_markers <- map_dfr(c(0:n), get_conserved)


## to look at gene expression, change the default assay from "integrated" to "RNA

DefaultAssay(covid.integrated) <- "RNA"

covid.integrated <- NormalizeData(covid.integrated, verbose = FALSE)

plot_density(covid.integrated, "OLIG2")


## Cluster Annotation

# Create reference matrix directly from UCSC cell browser
get_ucsc_reference <- function(cb_url,
                               cluster_col,
                               ...){
  
  url <- httr::parse_url(cb_url)
  base_url <- url
  ds <- url$query$ds
  base_url$query <- ""
  
  mdata_url <- httr::modify_url(base_url,
                                path = file.path(ds, "meta.tsv"))
  if(!httr::http_error(mdata_url)){
    mdata <- readr::read_tsv(mdata_url)
  } else {
    stop("unable to find metadata at url: ", mdata_url)
  }
  
  mat_url <- httr::modify_url(base_url,
                              path = file.path(ds, "exprMatrix.tsv.gz"))
  if(!httr::http_error(mat_url)){
    mat <- readr::read_tsv(mat_url)
  } else {
    stop("unable to find matrix at url: ", mat_url)
  }
  
  mat <- tibble::column_to_rownames(mat, "gene")
  mat <- as.matrix(mat)
  average_clusters(mat, mdata, cluster_col = cluster_col, ...)
}

#Primary Fetal cortex, Nowakowski et al., 2016
new_ref_matrix<- get_ucsc_reference(cb_url = "https://cells.ucsc.edu/?ds=cortex-dev",
                                    cluster_col = "WGCNAcluster")
## Choroid plexus organoids, Pelligrini et al., 2020
#new_ref_matrix<- get_ucsc_reference(cb_url = "https://cells.ucsc.edu/?ds=chporg",
#                                    cluster_col = "Cluster")


## Get a correlation matrix
res <- clustify(
  input = covid.integrated[["RNA"]]@data, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = covid.integrated@meta.data, # meta.data table containing cell clusters
  cluster_col = "Celltype_temp", # name of column in meta.data containing cell clusters
  ref_mat = new_ref_matrix) # matrix of RNA-seq expression data for each cell type


# Create heatmap of correlation coefficients 
plot_cor_heatmap(cor_mat = res, col = rev(brewer.pal(11,"RdBu")))

## Annotate cluster identities
new.idents<- setNames(c("Neurons", 
                        "Immature Neurons",
                        "Neurons", 
                        "Perivascular",
                        "IPC", 
                        "Radial Glia-1", 
                        "Vascular",
                        "Radial Glia-2",
                        "Perivascular",
                        "Astrocytes",
                        "Radial Glia-2",
                        "Unknown-1",
                        "Perivascular",
                        "Pericytes",
                        "Unknown-2", 
                        "Choroid Plexus",
                        "Cycling Radial Glia",
                        "oRG",
                        "Cycling Radial Glia",
                        "Choroid Plexus",
                        "IPC",
                        "Cycling Progenitors",
                        "Cycling Progenitors",
                        "Microglia",
                        "Pinealocytes", ## filtered out
                        "Epithelial-like",## filtered out
                        "Perivascular",
                        "Neurons",
                        "Choroid Plexus",
                        "Endothelial Precursors",
                        "Schwann Precursors",
                        "Perivascular",
                        "Radial Glia-2",
                        "Ciliogenesis"), # filtered out
                      levels(covid.integrated))
cov <- RenameIdents(covid.integrated, new.idents)

saveRDS(cov, "./data/Covid_final_data.rds")

sessionInfo()
