# Organoid-Covid19 Infectivity

Welcome to the Github repository associated to our preprint: 
**SARS-CoV-2 Neurotropism and Single Cell Responses in Brain Organoids Containing Innately Developing Microglia 

Samudyata1#, Ana Osório Oliveira1#, Susmita Malwade1, Nuno Rufino de Sousa2, Sravan Goparaju1, Funda Orhan1, Laura Steponaviciute2, Steven S Sheridan3, Roy H Perlis3, Antonio Rothfuchs2, Carl M. Sellgren1,4* 

1Department of Physiology and Pharmacology, Karolinska Institute, Stockholm, Sweden 

2Department of Microbiology, Tumor and Cell Biology, Karolinska Institute, Stockholm, Sweden 

3Center for Genomic Medicine and Department of Psychiatry, Massachusetts General Hospital, Boston, MA, USA. 

4Centre for Psychiatry Research, Department of Clinical Neuroscience, Karolinska Institutet & Stockholm Health Care Services, Stockholm County Council, Karolinska University Hospital, Stockholm, Sweden.* 



Corresponding Author:  carl.sellgren@ki.se

# Abstract
 
CNS manifestations are common both in the acute and post-acute phase of SARS-CoV-2 infection. Here, we derive human cerebral organoids with innately developing microglia to investigate the cellular responses to SARS-CoV-2 infection on a single cell level. We find the major cell types in the brain can harbor SARS-CoV-2 although the percentage of infected cells are relatively low. We observe an extensive neuronal cell death that also include non-infected cells. Single-cell RNA sequencing data from infected organoids reveals innate immune responses in microglia and astrocytes. Further, microglia acquire a phagocytic signature reminiscent of the disease associated phenotype, also observed in neurodegenerative diseases. Astrocytes transition into a proliferative and, reactive transcriptomic state with dysregulation of genes necessary for energizing neurons and maintaining the integrity of the blood-brain barrier. Across all cell types, we observe, dysregulated ions indicative of global translational shut-down, as well as altered carbohydrate metabolism and cellular respiration. Together, our findings provide insights into cellular responses of the resident brain immune cells to SARS-CoV-2 and pinpoint mechanisms that may be of relevance for the neuropathological changes observed in COVID-19 patients. 
 
 
 

 
 
 # You will find here
- library.R: 
- Covid_preprocessing.rmd: Summary of the processing steps performed to produce Figures 4,5 and 6. 
  - R Code used to integrate 10x datasets of 
  - Python Code used to visualize expression dynamics of 15q13.3 genes within these datasets. 
  
   
- organoids-infection_analysis.rmd: Summary of the processing steps performed to produce Figures 1 and 3. 
  - R Code used to caluclate maturations scores of single cells at E13.5 from MGE and E14.5 from CGE. (Adapted from Mayer et al.) 
  - R Code to visualize 15q13.3 gene expression in mitotic and postmitotic cells obtained from MGE and CGE. 
  
Mouse embryonic data used for this analysis was downloaded from [here](https://github.com/ChristophH/in-lineage/tree/master/data), and [here](https://github.com/jeremymsimon/MouseCortex).

Mouse Adult data used for this analysis was downloaded from [here](https://portal.brain-map.org/atlases-and-data/rnaseq). 

