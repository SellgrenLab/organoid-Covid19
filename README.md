# Organoid-Covid19 Infectivity

Welcome to the Github repository associated to our preprint: 
**Preprint name and authors* 

1 Department of Physiology and Phamacology, Karolinska Institutet, Stockholm, Sweden

Corresponding Author:  carl.sellgren@ki.se

# Abstract
 
 
 
 
 
 
 
 
 
 # You will find here
- library.R: Modifications made to Mayer et al.'s  lib.R file. Initially published [here](https://github.com/ChristophH/in-lineage).
- Covid_preprocessing.rmd: Summary of the processing steps performed to produce Figures 4,5 and 6. 
  - R Code used to integrate 10x datasets of developing mouse brain at E13.5, E18.5 and P10  and P56 to identify interneuron precursors during development. (Adapted from Mayer et al.) 
  - Python Code used to visualize expression dynamics of 15q13.3 genes within these datasets. 
  
   
- organoids-infection_analysis.rmd: Summary of the processing steps performed to produce Figures 1 and 3. 
  - R Code used to caluclate maturations scores of single cells at E13.5 from MGE and E14.5 from CGE. (Adapted from Mayer et al.) 
  - R Code to visualize 15q13.3 gene expression in mitotic and postmitotic cells obtained from MGE and CGE. 
  
Mouse embryonic data used for this analysis was downloaded from [here](https://github.com/ChristophH/in-lineage/tree/master/data), and [here](https://github.com/jeremymsimon/MouseCortex).

Mouse Adult data used for this analysis was downloaded from [here](https://portal.brain-map.org/atlases-and-data/rnaseq). 

