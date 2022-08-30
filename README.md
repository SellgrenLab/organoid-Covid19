# Organoid-Covid19 Neurotropism

Welcome to the Github repository associated to the biorxiv preprint: 

**SARS-CoV-2 promotes microglial synapse elimination in human brain organoids**, *Samudyata, Ana O. Oliveira, Susmita Malwade et al.* 


doi: https://doi.org/10.1101/2021.07.07.451463

*Corresponding Author: carl.sellgren@ki.se*

## Abstract
 
Neuropsychiatric manifestations are common in both the acute and post-acute phase of SARS-CoV-2 infection, but the mechanisms of these effects are unknown. In a newly established brain organoid model with innately developing microglia, we demonstrate that SARS-CoV-2 infection initiate neuronal cell death and cause a loss of post-synaptic termini. Despite limited neurotropism and a decelerating viral replication, we observe a threefold increase in microglial engulfment of postsynaptic termini after SARS-CoV-2 exposure. We define the microglial responses to SARS-CoV-2 infection by single cell transcriptomic profiling and observe an upregulation of interferon-responsive genes as well as genes promoting migration and synapse engulfment. To a large extent, SARS-CoV-2 exposed microglia adopt a transcriptomic profile overlapping with neurodegenerative disorders that display an early synapse loss as well as an increased incident risk after a SARS-CoV-2 infection. Our results reveal that brain organoids infected with SARS-CoV-2 display disruption in circuit integrity via microglia-mediated synapse elimination and identifies a potential novel mechanism contributing to cognitive impairments in patients recovering from COVID-19.
 
<p float="left">
  <img src="https://user-images.githubusercontent.com/56251389/187479169-e74aac1d-3d5b-4b6d-8b33-0b7e67805d2a.png" width="300" />
  <img src="https://user-images.githubusercontent.com/56251389/187480588-aa1af512-f1cd-4de5-ac20-6d7b1a54602d.png" width="300" />
  <img src="https://user-images.githubusercontent.com/56251389/187479107-21443e85-04bc-4d7d-8805-ccdf9d5f7e12.png" width="300" /> 
</p>



 ### You will find here
 
- oragnoid_env.yml: Conda environment file containing all softwares and packgaes along with their versions used and required for reproducibility.  
- Preprocessing.R: R Code to process 10X cellranger output files for all samples prior to downstream analysis. 
- Initial_clustering.R: R code to perform quality control on cells and genes and removal of doublets. Downstream processing to perform normalization, dimensionality reduction, scaling and first round of clustering to annotate celltypes was done using Seurat.
- Final_clustering.R: Second round of clustering upon cell cycle regression and integration of all samples and final celltype annotations used in this study using identified conserved cluster markers. Transcriptomic comparison to other datasets performed was also included.
- DifferentialExpression_analysis: Code to perform DE analysis on celltypes using Seurat as well as a pseudo-bulk approach in DESeq2. Functional analysis of pathwyas obtained via GSEA.
- Astrocyte_subclustering.R: Code applied to obtain subclusters of celltypes such as astrocyte or choroid plexus (as seen in the study) across all experimental condition. 
- Olag_MG_clusteroverlap.R: Comparison of organoid-grown microglial states upon infection to Microglia transcriptomic identitites found by Olah et al,2020.

*Note: This repository is only intended to provide insights into the analysis related to the key findings of the study and is not maintained or designed to be run directly. If you have any questions/doubts, please don't hesitate to contact the authors or raise an issue on the github repository.*


