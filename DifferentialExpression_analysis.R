## Differential Gene Expression Analysis



library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(ashr)
library(msigdbr)
library(fgsea)
library(DEGreport)

seurat <- readRDS(".data/Covid_final_dataset.rds")

## DE Analysis using Seurat

microglia<- subset(cov, cells = WhichCells(cov, idents = c("Microglia")))
Idents(microglia) <- "Condition"
DEGs_microglia<- FindAllMarkers(microglia, logfc.threshold = 0.2, test.use = "MAST", 
                                   min.pct = 0.1, min.diff.pct = 0.2,
                                   assay = "RNA")
## Save the DEGs in anexcel sheet
write.xlsx(DEGs_microglia, "Microglia_DEgenes across condition.xlsx")

### Gene set enrichment analysis (FGSEA)

#Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

#List available gene sets
unique(msigdbgmt$gs_subcat)

#Subset which gene set you want to use.
msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == "GO:BP",]
gmt <- lapply( unique(msigdbgmt_subset$gs_name),function(x){msigdbgmt_subset [msigdbgmt_subset$gs_name == x ,"gene_symbol"]} )
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name,"_",msigdbgmt_subset$gs_exact_source))

### Subset objects per condition
Idents(cov) <- "Condition"
cov72<- subset(cov, ident=c("72 HPI", "Control"))
cov24 <- subset(cov, ident=c("24 HPI", "Control"))
cov7224 <- subset(cov, ident=c("72 HPI", "24 HPI"))


Idents(cov) <- "celltype"
cell_selection <- subset(cov72, cells = WhichCells(cov72, idents = c("Microglia"))) ## change celltype here
cell_selection <- SetIdent(cell_selection, value = "Condition")
cell_selection$Condition <- factor(cell_selection$Condition, levels = c("72 HPI", "Control"))
DGE_cell_selection <- FindMarkers(cell_selection, ident.1= "72 HPI", ident.2="Control" , logfc.threshold = -Inf,test.use = "wilcox", min.pct = 0,  min.diff.pct = 0, only.pos = FALSE,  max.cells.per.ident = 1000, return.thresh = 1,  assay = "RNA")
gene_rank <- setNames( DGE_cell_selection$avg_logFC, casefold(rownames(DGE_cell_selection),upper=T))
fgseaRes <- fgsea( pathways=gmt, stats=gene_rank, minSize = 15, maxSize = 500,nperm = 10000)
fgseaRes <- fgseaRes[ order(fgseaRes$pval) ,]
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
#cols <- c("non-significant" = "grey", "significant" = "red")
fgseaResTidy <- dplyr::filter(fgseaResTidy, fgseaResTidy$adjPvalue=="significant")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = NES)) +
  geom_col() + theme_classic()+
  scale_fill_viridis(discrete = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO:BP Microglia 72 HPI vs Ctrl")
mg72<- fgseaResTidy

Idents(cov) <- "celltype"
cell_selection <- subset(cov24, cells = WhichCells(cov24, idents = c("Microglia"))) ## change celltype here
cell_selection <- SetIdent(cell_selection, value = "Condition")
cell_selection$Condition <- factor(cell_selection$Condition, levels = c("24 HPI", "Control"))
DGE_cell_selection <- FindMarkers(cell_selection, ident.1= "24 HPI", ident.2="Control" , logfc.threshold = -Inf,test.use = "wilcox", min.pct = 0,  min.diff.pct = 0, only.pos = FALSE,  max.cells.per.ident = 1000, return.thresh = 1,  assay = "RNA")
gene_rank <- setNames( DGE_cell_selection$avg_logFC, casefold(rownames(DGE_cell_selection),upper=T))
fgseaRes <- fgsea( pathways=gmt, stats=gene_rank, minSize = 15, maxSize = 500,nperm = 10000)
fgseaRes <- fgseaRes[ order(fgseaRes$pval) ,]
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
#cols <- c("non-significant" = "grey", "significant" = "red")
fgseaResTidy <- dplyr::filter(fgseaResTidy, fgseaResTidy$adjPvalue=="significant")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = NES)) +
  geom_col() + theme_classic()+
  scale_fill_viridis(discrete = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO:BP Microglia 24 HPI vs Ctrl")
mg24<- fgseaResTidy

##Save results
write.xlsx(list("MG 24hpi vs Ctrl"=mg24, "MG 72hpi vs Ctrl"=mg72), file = "FGSEA-Microglia Per condition.xlsx")


###Pseudobulk DE analysis with DESeq2

dir.create("DESeq2")
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$sample_id <- factor(seurat@sample_id)
metadata$cluster_id <- factor(seurat$seurat_clusters)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids
# Total number of clusters
nk <- length(kids)
nk
# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))
# Total number of samples 
ns <- length(sids)
ns
# Generate sample level metadata
## Determine the number of cells per sample
table(sce$sample_id)
## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
class(pb)
dim(pb)

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)
str(pb)
# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)

# Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples
samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()
# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 

metadata        

# Generate vector of cluster IDs
metadata$cluster_id <- factor(metadata$cluster_id)
clusters <- levels(metadata$cluster_id)
clusters

### Function adapted to run pseudobulk LRT test on all clusters

get_dds_LRTresults <- function(x){
  
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  #all(rownames(cluster_metadata) == colnames(cluster_counts))        
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  
  dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
  
  # Extract results
  res_LRT <- results(dds_lrt)
  
  # Create a tibble for LRT results
  res_LRT_tb <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Save all results
  write.csv(res_LRT_tb,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset to return genes with padj < 0.05
  sigLRT_genes <- res_LRT_tb %>% 
    filter(padj < 0.05)
  
  # Save sig results
  write.csv(sigLRT_genes,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Transform counts for data visualization
  rld <- rlog(dds_lrt, blind=TRUE)
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  
  # Obtain rlog values for those significant genes
  cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
  
  cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
  
  # # Remove samples without replicates
  # cluster_rlog <- cluster_rlog[, -1]
  # cluster_metadata <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
  
  
  # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
  cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig, time = "group_id", col=NULL)
  ggsave(paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.png"))
  
  # Let's see what is stored in the `df` component
  write.csv(cluster_groups$df,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  saveRDS(cluster_groups, paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.rds"))
  save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, file = paste0("DESeq2/lrt/", clusters[x], "_all_LRTresults.Rdata"))
  
}

map(1:length(clusters), get_dds_LRTresults)



