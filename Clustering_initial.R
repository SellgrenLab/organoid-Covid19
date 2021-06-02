library(tidyverse)
library(Seurat)
library(scran)
library(patchwork)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(scDblFinder)
library(SingleR)



custom_colors <- list()

## This script is customised from https://romanhaa.github.io/projects/scrnaseq_workflow/

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)

niceFormat <- function(number) {
  formatC(number, format = 'f', big.mark = ',', digits = 0)
}



## Read in  count matrices
transcript_counts<- readRDS("Covid_CountMatrices.rds")

transcripts <- list(
  raw = list()
)
## Save raw counts
transcripts$raw <- transcript_counts

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

averagenCountnFeature <- function(cells) {
  output <- tibble(
    'sample' = character(),
    'mean(nCount)' = numeric(),
    'median(nCount)' = numeric(),
    'mean(nFeature)' = numeric(),
    'median(nFeature)' = numeric()
  )
  for ( i in levels(cells$sample) ) {
    tmp <- tibble(
      'sample' = i,
      'mean(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% mean(),
      'median(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% median(),
      'mean(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% mean(),
      'median(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% median()
    )
    output <- bind_rows(output, tmp)
  }
  return(output)
}

## To save output
dir.create('data')
dir.create('plots')

sameGeneNames(transcripts$raw)

sample_names <- names(transcript_counts)

transcripts$raw$merged <- do.call("cbind", transcript_counts)

## Quality control

cells <- experiment.aggregate@meta.data

## Doublets plot
p <- cells %>%
  filter(scDblFinder.class != 'singlet') %>%
  group_by(ControllerLane) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = ControllerLane, y = count, fill = ControllerLane)) +
  geom_col(color = 'black') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$ControllerLane) +
                                  scale_y_continuous(name = 'Number of doublets', labels = scales::comma) +
                                  theme(
                                    axis.title.y = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                #ggsave('plots/qc_number_of_doublets_by_sample.png',p, height = 3, width = 6)
                                
                                
                                p <- ggplot(cells, aes(nCount_RNA, nFeature_RNA, color = percent.mito)) +
                                  geom_point(size = 0.5) +
                                  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
                                  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
                                  theme_bw() +
                                  scale_color_viridis(
                                    name = 'Percent MT\ntranscripts',
                                    limits = c(0,1),
                                    labels = scales::percent,
                                    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
                                  )
                                
                                #ggsave('plots/qc_nFeature_over_nCount_before_filtering.png', p, height = 4, width = 6)
                                
                                ## Calculate threshold values
                                
                                median_nCount <- median(cells$nCount_RNA)
                                # 4255
                                mad_nCount <- mad(cells$nCount_RNA)
                                # 2630.874
                                
                                median_nFeature <- median(cells$nFeature_RNA)
                                # 803
                                mad_nFeature <- mad(cells$nFeature_RNA)
                                # 544.1142
                                
                                median_percent_MT <- median(cells$percent.mito)
                                # 0.02760973
                                mad_percent_MT <- mad(cells$percent.mito)
                                # 0.02103674
                                
                                thresholds_nCount <- c(0, median_nCount + 5*mad_nCount)
                                # 0.00 17409.37
                                thresholds_nFeature <- c(0, median_nFeature + 5*mad_nFeature)
                                # 0.000 3523.571
                                thresholds_percent_MT <- c(0, median_percent_MT + 5*mad_percent_MT)
                                
                                ## Plots with thresholds
                                cells$sample <- cells$ControllerLane
                                
                                p1 <- ggplot(cells, aes(x = sample, y = nCount_RNA, fill = sample)) +
                                  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
                                  geom_hline(yintercept = median_nCount, color = 'black') +
                                  geom_hline(yintercept = thresholds_nCount, color = 'red') +
                                  theme_bw() +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_x_discrete(limits = rev(levels(cells$sample))) +
                                  scale_y_continuous(labels = scales::comma) +
                                  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
                                  theme(
                                    axis.title = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                p2 <- ggplot(cells, aes(x = sample, y = nCount_RNA, fill = sample)) +
                                  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
                                  geom_hline(yintercept = median_nCount, color = 'black') +
                                  geom_hline(yintercept = thresholds_nCount, color = 'red') +
                                  theme_bw() +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_x_discrete(limits = rev(levels(cells$sample))) +
                                  scale_y_log10(labels = scales::comma) +
                                  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
                                  theme(
                                    axis.title = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                p3 <- ggplot(cells, aes(x = sample, y = nFeature_RNA, fill = sample)) +
                                  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
                                  geom_hline(yintercept = median_nFeature, color = 'black') +
                                  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
                                  theme_bw() +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_x_discrete(limits = rev(levels(cells$sample))) +
                                  scale_y_continuous(labels = scales::comma) +
                                  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
                                  theme(
                                    axis.title = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                p4 <- ggplot(cells, aes(x = sample, y = nFeature_RNA, fill = sample)) +
                                  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
                                  geom_hline(yintercept = median_nFeature, color = 'black') +
                                  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
                                  theme_bw() +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_x_discrete(limits = rev(levels(cells$sample))) +
                                  scale_y_log10(labels = scales::comma) +
                                  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
                                  theme(
                                    axis.title = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                p5 <- ggplot(cells, aes(x = sample, y = percent.mito, fill = sample)) +
                                  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
                                  geom_hline(yintercept = median_percent_MT, color = 'black') +
                                  geom_hline(yintercept = thresholds_percent_MT, color = 'red') +
                                  theme_bw() +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_x_discrete(limits = rev(levels(cells$sample))) +
                                  scale_y_continuous(labels = scales::percent) +
                                  labs(title = 'Percent MT transcripts', subtitle = 'linear scale') +
                                  theme(
                                    axis.title = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    legend.position = 'none'
                                  ) +
                                  coord_flip()
                                
                                
                                #ggsave('plots/qc_histogram_nCount_nFeature_percentMT_thresholds.png',p1 + p3 + p5 +  p2 + p4 + plot_layout(ncol = 3),height = 7, width = 10  )
                                
                                
                                
                                p <- ggplot(cells, aes(nCount_RNA, nFeature_RNA, color = percent.mito)) +
                                  geom_point(size = 0.5) +
                                  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
                                  geom_vline(xintercept = thresholds_nCount, color = 'red') +
                                  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
                                  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
                                  theme_bw() +
                                  scale_color_viridis(
                                    name = 'Percent MT\ntranscripts',
                                    limits = c(0,1),
                                    labels = scales::percent,
                                    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
                                  )
                                
                                #ggsave('plots/qc_nFeature_over_nCount_thresholds.png', p, height = 4, width = 6)
                                
                                ## Filtering out cells and genes
                                
                                cells_filtered <- cells %>%
                                  dplyr::filter(
                                    nCount >= thresholds_nCount[1],
                                    nCount <= thresholds_nCount[2],
                                    nFeature >= thresholds_nFeature[1],
                                    nFeature <= thresholds_nFeature[2],
                                    percent_MT >= thresholds_percent_MT[1],
                                    percent_MT <= thresholds_percent_MT[2]
                                  )
                                
                                DataFrame(cells_filtered)
                                
                                cells_to_keep <- cells_filtered$cell
                                length(cells_to_keep)
                                
                                genes <- tibble(
                                  gene = rownames(transcripts$raw$merged),
                                  count = Matrix::rowSums(transcripts$raw$merged),
                                  cells = Matrix::rowSums(transcripts$raw$merged != 0)
                                )
                                
                                genes_to_keep <- genes %>% dplyr::filter(cells >= 5) %>% pull(gene)
                                length(genes_to_keep)
                                
                                
                                ## Post filtering
                                
                                transcripts$raw$filtered <- transcripts$raw$merged[genes_to_keep,cells_to_keep]
                                
                                cells_per_sample_after_filtering <- tibble(
                                  sample = character(),
                                  before = numeric(),
                                  after = numeric()
                                )
                                
                                for ( i in sample_names ) {
                                  tmp <- tibble(
                                    sample = i,
                                    before = grep(colnames(transcripts$raw$merged), pattern = paste0('-', i, '$')) %>% length(),
                                    after = grep(colnames(transcripts$raw$filtered), pattern = paste0('-', i, '$')) %>% length()
                                  )
                                  cells_per_sample_after_filtering <- bind_rows(cells_per_sample_after_filtering, tmp)
                                }
                                
                                knitr::kable(cells_per_sample_after_filtering)
                                
                                seurat <- CreateSeuratObject(
                                  counts = transcripts$raw$filtered,
                                  min.cells = 0,
                                  min.features = 0, meta.data= cells_filtered
                                )
                                
                                
                                
                                temp_labels <- seurat@meta.data %>%
                                  group_by(sample) %>%
                                  tally()
                                
                                p1 <- ggplot() +
                                  geom_half_violin(
                                    data = seurat@meta.data, aes(sample, nCount_RNA, fill = sample),
                                    side = 'l', show.legend = FALSE, trim = FALSE
                                  ) +
                                  geom_half_boxplot(
                                    data = seurat@meta.data, aes(sample, nCount_RNA, fill = sample),
                                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
                                  ) +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = sample, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_color_manual(values = custom_colors$discrete) +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_y_continuous(labels = scales::comma, expand = c(0.08,0)) +
                                  theme_bw() +
                                  labs(x = '', y = 'Number of transcripts') +
                                  theme(
                                    panel.grid.major.x = element_blank(),
                                    axis.title.x = element_blank()
                                  )
                                
                                p2 <- ggplot() +
                                  geom_half_violin(
                                    data = seurat@meta.data, aes(sample, nFeature_RNA, fill = sample),
                                    side = 'l', show.legend = FALSE, trim = FALSE
                                  ) +
                                  geom_half_boxplot(
                                    data = seurat@meta.data, aes(sample, nFeature_RNA, fill = sample),
                                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
                                  ) +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = sample, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_color_manual(values = custom_colors$discrete) +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_y_continuous(
                                    name = 'Number of expressed genes',
                                    labels = scales::comma, expand = c(0.08,0)
                                  ) +
                                  theme_bw() +
                                  theme(
                                    panel.grid.major.x = element_blank(),
                                    axis.title.x = element_blank()
                                  )
                                
                                #ggsave(
                                #  'plots/nCount_nFeature_by_sample.png',
                                #  p1 + p2 + plot_layout(ncol = 2), height = 4, width = 10
                                #)
                                
                                
                                ## Integrate data from all 6 samples together using Seurat
                                
                                #seurat <- SCTransform(seurat, assay = 'RNA')
                                
                                ## Data integration
                                
                                seurat_list <- SplitObject(seurat, split.by = 'sample')
                                
                                seurat_list <- lapply(
                                  X = seurat_list,
                                  FUN = function(x) {
                                    x <- SCTransform(x)
                                  }
                                )
                                
                                seurat_features <- SelectIntegrationFeatures(
                                  seurat_list,
                                  nfeatures = 3000
                                )
                                
                                seurat_list <- PrepSCTIntegration(
                                  seurat_list,
                                  anchor.features = seurat_features
                                )
                                
                                seurat_list <- lapply(
                                  X = seurat_list,
                                  FUN = RunPCA,
                                  verbose = FALSE,
                                  features = seurat_features
                                )
                                
                                seurat_anchors <- FindIntegrationAnchors(
                                  object.list = seurat_list,
                                  anchor.features = seurat_features,
                                  normalization.method = 'SCT'
                                )
                                
                                
                                seurat <- IntegrateData(
                                  anchorset = seurat_anchors,
                                  normalization.method = 'SCT'
                                )
                                
                                seurat@meta.data$sample <- Cells(seurat) %>%
                                  strsplit('-') %>%
                                  vapply(FUN.VALUE = character(1), `[`, 2)
                                seurat@meta.data$sample <- seurat@meta.data$sample %>%
                                  factor(levels = sample_names[which(sample_names %in% .)])
                                
                                seurat <- RunPCA(seurat, assay = 'SCT', npcs = 50)
                                
                                intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 10)
                                # 10.44694
                                
                                p <- tibble(
                                  PC = 1:50,
                                  stdev = seurat@reductions$pca@stdev
                                ) %>%
                                  ggplot(aes(PC, stdev)) +
                                  geom_point() +
                                  geom_vline(xintercept = 10, color = 'blue') +
                                  geom_vline(xintercept = 15, color = 'red') +
                                  theme_bw() +
                                  labs(x = 'Principal components', y = 'Standard deviation')
                                
                                #ggsave('plots/principal_components.png', p, height = 4, width = 5)
                                
                                
                                ### Clustering
                                
                                ## First round of clustering
                                
                                seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:30)
                                seurat <- FindClusters(seurat, resolution = 0.8)
                                
                                
                                # check no of transcripts
                                temp_labels <- seurat@meta.data %>%
                                  group_by(seurat_clusters) %>%
                                  tally()
                                
                                p1 <- ggplot() +
                                  geom_half_violin(
                                    data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
                                    side = 'l', show.legend = FALSE, trim = FALSE
                                  ) +
                                  geom_half_boxplot(
                                    data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
                                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
                                  ) +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_color_manual(values = custom_colors$discrete) +
                                  scale_fill_manual(values = custom_colors$discrete) +
                                  scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
                                  theme_bw() +
                                  theme(
                                    panel.grid.major.x = element_blank(),
                                    axis.title.x = element_blank()
                                  )
                                
                                ## Cluster stability
                                
                                sce <- as.SingleCellExperiment(seurat)
                                reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:15, drop = FALSE]
                                
                                ass_prob <- bootstrapCluster(sce, FUN = function(x) {
                                  g <- buildSNNGraph(x, use.dimred = 'PCA_sub')
                                  igraph::cluster_walktrap(g)$membership
                                },
                                clusters = sce$seurat_clusters
                                )
                                
                                p <- ass_prob %>%
                                  as_tibble() %>%
                                  rownames_to_column(var = 'cluster_1') %>%
                                  pivot_longer(
                                    cols = 2:ncol(.),
                                    names_to = 'cluster_2',
                                    values_to = 'probability'
                                  ) %>%
                                  mutate(
                                    cluster_1 = as.character(as.numeric(cluster_1) - 1),
                                    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
                                    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
                                  ) %>%
                                  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
                                  geom_tile(color = 'white') +
                                  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
                                  scale_x_discrete(name = 'Cluster', position = 'top') +
                                  scale_y_discrete(name = 'Cluster') +
                                  scale_fill_gradient(
                                    name = 'Probability', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
                                    limits = c(0,1),
                                    guide = guide_colorbar(
                                      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
                                      title.theme = element_text(hjust = 1, angle = 90),
                                      barwidth = 0.75, barheight = 10
                                    )
                                  ) +
                                  coord_fixed() +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'right',
                                    panel.grid.major = element_blank()
                                  )
                                
                                #ggsave('plots/cluster_stability.png', p, height = 6, width = 7)
                                
                                
                                ## cluster similarity- Supplementary
                                
                                sce <- as.SingleCellExperiment(seurat)
                                
                                reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:15, drop = FALSE]
                                
                                g <- scran::buildSNNGraph(sce, use.dimred = 'PCA_sub')
                                
                                ratio <- scran::clusterModularity(g, seurat@meta.data$seurat_clusters, as.ratio = TRUE)
                                
                                ratio_to_plot <- log10(ratio+1)
                                
                                p <- ratio_to_plot %>%
                                  as_tibble() %>%
                                  rownames_to_column(var = 'cluster_1') %>%
                                  pivot_longer(
                                    cols = 2:ncol(.),
                                    names_to = 'cluster_2',
                                    values_to = 'probability'
                                  ) %>%
                                  mutate(
                                    cluster_1 = as.character(as.numeric(cluster_1) - 1),
                                    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
                                    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
                                  ) %>%
                                  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
                                  geom_tile(color = 'white') +
                                  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
                                  scale_x_discrete(name = 'Cluster', position = 'top') +
                                  scale_y_discrete(name = 'Cluster') +
                                  scale_fill_gradient(
                                    name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
                                    guide = guide_colorbar(
                                      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
                                      title.theme = element_text(hjust = 1, angle = 90),
                                      barwidth = 0.75, barheight = 10
                                    )
                                  ) +
                                  coord_fixed() +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'right',
                                    panel.grid.major = element_blank()
                                  )
                                
                                #ggsave('plots/cluster_similarity.png', p, height = 6, width = 7)
                                
                                ## cluster tree
                                seurat <- BuildClusterTree(
                                  seurat,
                                  dims = 1:15,
                                  reorder = FALSE,
                                  reorder.numeric = FALSE
                                )
                                
                                tree <- seurat@tools$BuildClusterTree
                                tree$tip.label <- paste0("Cluster ", tree$tip.label)
                                
                                p <- ggtree::ggtree(tree, aes(x, y)) +
                                  scale_y_reverse() +
                                  ggtree::geom_tree() +
                                  ggtree::theme_tree() +
                                  ggtree::geom_tiplab(offset = 1) +
                                  ggtree::geom_tippoint(color = custom_colors$discrete[1:length(tree$tip.label)], shape = 16, size = 5) +
                                  coord_cartesian(clip = 'off') +
                                  theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))
                                
                                #ggsave('plots/cluster_tree.png', p, height = 4, width = 6)
                                
                                ## composition plots
                                table_samples_by_clusters <- seurat@meta.data %>%
                                  group_by(sample, seurat_clusters) %>%
                                  summarize(count = n()) %>%
                                  spread(seurat_clusters, count, fill = 0) %>%
                                  ungroup() %>%
                                  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
                                  dplyr::select(c('sample', 'total_cell_count', everything())) %>%
                                  arrange(factor(sample, levels = levels(seurat@meta.data$sample)))
                                
                                knitr::kable(table_samples_by_clusters)
                                
                                table_clusters_by_samples <- seurat@meta.data %>%
                                  dplyr::rename('cluster' = 'seurat_clusters') %>%
                                  group_by(cluster, sample) %>%
                                  summarize(count = n()) %>%
                                  spread(sample, count, fill = 0) %>%
                                  ungroup() %>%
                                  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
                                  select(c('cluster', 'total_cell_count', everything())) %>%
                                  arrange(factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)))
                                
                                knitr::kable(table_clusters_by_samples)
                                
                                
                                temp_labels <- seurat@meta.data %>%
                                  group_by(sample) %>%
                                  tally()
                                
                                p1 <- table_samples_by_clusters %>%
                                  select(-c('total_cell_count')) %>%
                                  reshape2::melt(id.vars = 'sample') %>%
                                  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
                                  ggplot(aes(sample, value)) +
                                  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
                                  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
                                  coord_cartesian(clip = 'off') +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'left',
                                    plot.title = element_text(hjust = 0.5),
                                    text = element_text(size = 16),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
                                  )
                                
                                temp_labels <- seurat@meta.data %>%
                                  group_by(seurat_clusters) %>%
                                  tally() %>%
                                  dplyr::rename('cluster' = seurat_clusters)
                                
                                p2 <- table_clusters_by_samples %>%
                                  select(-c('total_cell_count')) %>%
                                  reshape2::melt(id.vars = 'cluster') %>%
                                  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
                                  ggplot(aes(cluster, value)) +
                                  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
                                  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
                                  coord_cartesian(clip = 'off') +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'right',
                                    plot.title = element_text(hjust = 0.5),
                                    text = element_text(size = 16),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.title = element_blank(),
                                    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
                                  )
                                
                                #ggsave(
                                # 'plots/composition_samples_clusters_by_number.png',
                                #  p1 + p2 +
                                #    plot_layout(ncol = 2, widths = c(
                                #      seurat@meta.data$sample %>% unique() %>% length(),
                                #      seurat@meta.data$seurat_clusters %>% unique() %>% length()
                                #    )),
                                #  width = 18, height = 8
                                #)
                                
                                ## for percentage of composition
                                temp_labels <- seurat@meta.data %>%
                                  group_by(sample) %>%
                                  tally()
                                
                                p1 <- table_samples_by_clusters %>%
                                  select(-c('total_cell_count')) %>%
                                  reshape2::melt(id.vars = 'sample') %>%
                                  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
                                  ggplot(aes(sample, value)) +
                                  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
                                  geom_text(
                                    data = temp_labels,
                                    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
                                  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
                                  coord_cartesian(clip = 'off') +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'left',
                                    plot.title = element_text(hjust = 0.5),
                                    text = element_text(size = 16),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
                                  )
                                
                                temp_labels <- seurat@meta.data %>%
                                  group_by(seurat_clusters) %>%
                                  tally() %>%
                                  dplyr::rename('cluster' = seurat_clusters)
                                
                                p2 <- table_clusters_by_samples %>%
                                  select(-c('total_cell_count')) %>%
                                  reshape2::melt(id.vars = 'cluster') %>%
                                  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
                                  ggplot(aes(cluster, value)) +
                                  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
                                  geom_text(
                                    data = temp_labels, aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
                                    color = 'black', size = 2.8
                                  ) +
                                  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
                                  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
                                  coord_cartesian(clip = 'off') +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'right',
                                    plot.title = element_text(hjust = 0.5),
                                    text = element_text(size = 16),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.title = element_blank(),
                                    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
                                  )
                                
                                #ggsave(
                                #  'plots/composition_samples_clusters_by_percent.png',
                                #  p1 + p2 +
                                #    plot_layout(ncol = 2, widths = c(
                                #      seurat@meta.data$sample %>% unique() %>% length(),
                                #      seurat@meta.data$seurat_clusters %>% unique() %>% length()
                                #    )),
                                #  width = 18, height = 8
                                #)
                                
                                
                                ## Score all cells based on their cell cycle
                                
                                seurat <- CellCycleScoring(
                                  seurat,
                                  assay = 'SCT',
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes
                                )
                                
                                seurat@meta.data$cell_cycle_seurat <- seurat@meta.data$Phase
                                seurat@meta.data$Phase <- NULL
                                seurat@meta.data$cell_cycle_seurat <- factor(
                                  seurat@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M')
                                )
                                
                                
                                ## UMAP
                                seurat <- RunUMAP(
                                  seurat,
                                  reduction.name = 'UMAP',
                                  reduction = 'pca',
                                  dims = 1:30,
                                  seed.use = 100
                                )
                                
                                
                                
                                ## Alluvial Plots- supplementary
                                
                                ## get sample and cluster names
                                samples <- levels(seurat@meta.data$sample)
                                clusters <- levels(seurat@meta.data$seurat_clusters)
                                
                                ## create named vector holding the color assignments for both samples and
                                ## clusters
                                color_assignments <- setNames(
                                  c(custom_colors$discrete[1:length(samples)], custom_colors$discrete[1:length(clusters)]),
                                  c(samples,clusters)
                                )
                                
                                ## prepare data for the plot; factor() calls are necessary for the right order
                                ## of columns (first samples then clusters) and boxes within each column (
                                ## cluster 1, 2, 3, ..., not 1, 10, 11, ...)
                                data <- seurat@meta.data %>%
                                  group_by(sample,seurat_clusters) %>%
                                  tally() %>%
                                  ungroup() %>%
                                  gather_set_data(1:2) %>%
                                  dplyr::mutate(
                                    x = factor(x, levels = unique(x)),
                                    y = factor(y, levels = unique(y))
                                  )
                                
                                DataFrame(data)
                                
                                
                                ## create sample and cluster labels; hjust defines whether a label will be
                                ## aligned to the right (1) or to the left (0); the nudge_x parameter is used
                                ## to move the label outside of the boxes
                                data_labels <- tibble(
                                  group = c(
                                    rep('sample', length(samples)),
                                    rep('seurat_clusters', length(clusters))
                                  )
                                ) %>%
                                  mutate(
                                    hjust = ifelse(group == 'sample', 1, 0),
                                    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
                                  )
                                
                                DataFrame(data_labels)
                                
                                
                                ## create plot
                                p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
                                  geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
                                  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
                                  geom_text(
                                    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
                                    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
                                  ) +
                                  scale_x_discrete(labels = c('Sample','Cluster')) +
                                  scale_fill_manual(values = color_assignments) +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'none',
                                    axis.title = element_blank(),
                                    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
                                    axis.text.y = element_blank(),
                                    axis.ticks = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank()
                                  )
                                
                                clusters <- levels(seurat@meta.data$seurat_clusters)
                                cell_types <- sort(unique(seurat@meta.data$cell_type_singler_blueprintencode_main))
                                
                                color_assignments <- setNames(
                                  c(custom_colors$discrete[1:length(clusters)], custom_colors$discrete[1:length(cell_types)]),
                                  c(clusters,cell_types)
                                )
                                
                                data <- seurat@meta.data %>%
                                  dplyr::rename(cell_type = cell_type_singler_blueprintencode_main) %>%
                                  dplyr::mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
                                  group_by(seurat_clusters, cell_type) %>%
                                  tally() %>%
                                  ungroup() %>%
                                  gather_set_data(1:2) %>%
                                  dplyr::mutate(
                                    x = factor(x, levels = unique(x)),
                                    y = factor(y, levels = c(clusters,cell_types))
                                  )
                                
                                data_labels <- tibble(
                                  group = c(
                                    rep('seurat_clusters', length(clusters)),
                                    rep('cell_type', length(cell_types))
                                  )
                                ) %>%
                                  mutate(
                                    hjust = ifelse(group == 'seurat_clusters', 1, 0),
                                    nudge_x = ifelse(group == 'seurat_clusters', -0.1, 0.1)
                                  )
                                
                                p2 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
                                  geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
                                  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
                                  geom_text(
                                    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
                                    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
                                  ) +
                                  scale_x_discrete(labels = c('Cluster','Cell type')) +
                                  scale_fill_manual(values = color_assignments) +
                                  theme_bw() +
                                  theme(
                                    legend.position = 'none',
                                    axis.title = element_blank(),
                                    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
                                    axis.text.y = element_blank(),
                                    axis.ticks = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank()
                                  )
                                
                                ## Supplementary figure
                                #ggsave(
                                #  'plots/samples_clusters_cell_types_alluvial.png',
                                #  p1 + p2 + plot_layout(ncol = 2),
                                #  height = 6, width = 8
                                #)
                                
  #### Celltype Identification and Annotation
                                
    new.idents<- setNames(c("Vascular Meningeal Cells",
                                                        "Newborn Neurons (Inhibitory)",
                                                        "Excitatory Neurons", 
                                                        "Astro-Radial Glia",
                                                        "Ectomesenchyme 1",
                                                        "Ectomesenchyme 2",
                                                        "NG2-expressing Mural cells",
                                                        "Unknown 1 (Mesoderm)",
                                                        "Myoblast_Skeletel Differentiation",
                                                        "Oligodendroglia",
                                                        "Proliferative Neuroepithelium",
                                                        "Unknown 1 (Cell cycle arrest)",
                                                        "Epithelial ACE2-expressing",
                                                        "Proliferative Neural Precursors",
                                                        "Choroid",
                                                        "Unknown 2 (Mesoderm)",
                                                        "Myocytes",
                                                        "Pericytes",
                                                        "Microglia",
                                                        "Unknown 2 (Nervous System Regulation)",
                                                        "Endothelial"), levels(seurat))

sars <- RenameIdents(seurat, new.idents)
saveRDS(sars, "Covid_clustered_1.rds")