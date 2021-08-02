# Organoid-Covid19 Infectivity

Welcome to the Github repository associated to the biorxiv preprint: 
**SARS-CoV-2 Neurotropism and Single Cell Responses in Brain Organoids Containing Innately Developing Microglia 

Samudyata1#, Ana Osório Oliveira1#, Susmita Malwade1, Nuno Rufino de Sousa2, Sravan Goparaju1, Funda Orhan1, Laura Steponaviciute2, Steven S Sheridan3, Roy H Perlis3, Antonio Rothfuchs2, Carl M. Sellgren1,4* 

1Department of Physiology and Pharmacology, Karolinska Institute, Stockholm, Sweden 

2Department of Microbiology, Tumor and Cell Biology, Karolinska Institute, Stockholm, Sweden 

3Center for Genomic Medicine and Department of Psychiatry, Massachusetts General Hospital, Boston, MA, USA. 

4Centre for Psychiatry Research, Department of Clinical Neuroscience, Karolinska Institutet & Stockholm Health Care Services, Stockholm County Council, Karolinska University Hospital, Stockholm, Sweden.* 

doi: https://doi.org/10.1101/2021.07.07.451463

Corresponding Author:  carl.sellgren@ki.se

## Abstract
 
Neuropsychiatric manifestations are common in both acute and post-acute phase of SARS-CoV-2 infection, but the mechanism of these effects is unknown. Here, we derive human brain organoids with innately developing microglia to investigate the cellular responses to SARS-CoV-2 infection on a single cell level. We find evidence of limited tropism to SARS-CoV-2 for all major cell types and observe extensive neuronal cell death that also include non-infected cells. Single cell transcriptome profiling reveals distinct responses in microglia and astrocytes that share features with cellular states observed in neurodegenerative diseases, includes upregulation of genes with relevance for synaptic stripping, and suggests altered blood brain barrier integrity. Across all cell types, we observe a global translational shut-down as well as altered carbohydrate metabolism and cellular respiration. Together, our findings provide insights into cellular responses of the resident brain immune cells to SARS-CoV-2 and pinpoint mechanisms that may be of relevance for the neuropathological changes observed in COVID-19 patients. 
 
 
 
 
 ### You will find here
- oragnoid_env.yml: Conda environment file containing all softwares and packgaes along with their versions used and required for reproducibility.  
- Preprocessing.R: R Code to process 10X cellranger output files for all samples prior to downstream analysis. 
- Initial_clustering.R: R code to perform quality control on cells and genes and removal of doublets. Downstream processing to perform normalization, dimensionality reduction, scaling and first round of clustering to annotate celltypes was done using Seurat.
- Final_clustering: Second round of clustering upon cell cycle regression and integration of all samples and final celltype annotations used in this study using identified conserved cluster markers. Transcriptomic comparison to other datasets performed was also included.
- DifferentialExpression_analysis: Code to perform DE analysis on celltypes using Seurat as well as a pseudo-bulk approach in DESeq2. Functional analysis of pathwyas obtained via GSEA.
- Astrocyte_subclustering: Code applied to obtain subclusters of celltypes such as astrocyte or choroid plexus (as seen in the study) across all experimental condition. 




