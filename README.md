# Human Fetal Heart Atlas

## Description
+ Forty first-trimester human hearts were studied to lay groundwork for further studies of principles underlying congenital heart defects.
+ We first sampled 49,227 cardiac nuclei from three fetuses at 8.6, 9.0, and 10.7 post-conceptional weeks (pcw) for single-nucleus RNA
+ sequencing, enabling distinction of six classes and 21 cell types. Improved resolution led to identification of novel cardiomyocytes
+ and minority autonomic and lymphatic endothelial transcriptomes, among others. After integration with 5-7 pcw heart single-cell RNAseq,
+ we identified a human cardiomyofibroblast progenitor preceding diversification of cardiomyocyte and stromal lineages. Analysis of six
+ Visium sections from two additional hearts was aided by deconvolution, and key spatial markers validated on sectioned and whole hearts
+ in two- and three-dimensional space and over time. Altogether, anatomical-positional features including innervation, conduction and
+ subdomains of the atrioventricular septum translate latent molecular identity into specialized cardiac functions. This atlas adds
+ unprecedented spatial and temporal resolution to the characterization of human-specific aspects of early heart formation.
+ 
+ The preprint can be downloaded from ([bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.21.624698v1)).


## Directory structure
* The code directory contains the R scripts used for:
   + snRNA data analysis,
   + visualization,
   + spatial transcriptomic analysis,
   + snRNA and spatial transcriptomic data integration procedures,
   + deconvolution, 
   + trajectory analysis, and
   + application of CellChat - https://github.com/jinworks/CellChat
   
  
* figures/ directory contains html and pdf files of the generated figures 

* results/ directory contains results of snRNA data analysis, its integration with spatial transcriptomics data and trajectory analysis  

* **docs/ 


## Data availability

* Sequencing data, images and snRNAseq matrices are available at ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283967))
* Accompanying Movies, large images and high-resolution gross anatomy images are available from [Figshare](https://figshare.com/projects/Multi-modal_refinement_of_the_human_heart_atlas_during_the_first_gestational_trimester/213151).


## Help

* Contact details: heather.etchevers@inserm.fr; stephane.zaffran@inserm.fr
