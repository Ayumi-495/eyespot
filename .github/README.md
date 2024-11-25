[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13147019.svg)](https://doi.org/10.5281/zenodo.13147019)

# A systematic review and meta-analysis of eyespot anti-predator mechanisms

This repository contains the all Rscript and data of **_A systematic review and meta-analysis of eyespot anti-predator mechanisms_**. Our meta-analytical study examined why butterfly eyespots intimidate birds and tested two hypotheses: the eye mimicry hypothesis and the conspicuousness hypothesis. Despite the existence of conflicting past empirical research, our results supported the idea of the conspicuousness hypothesis. We also found a key factor: larger eyespots were more effective in predator avoidance. Our findings suggest that birds avoid eyespots mainly due to their conspicuousness and emphasise the importance of such pattern characteristics in predator avoidance.

HTML file (https://ayumi-495.github.io/eyespot/) has code, metadata, and data within it. You can download them directly from there. 
Supplementary pdf is included supplimenal figures and tables.

## Description of the files in each folder:

**data**

- all.csv: Raw data including both predator and prey datasets. This dataset was used for calculating effect sizes.
 
- butterfly_eyespots.csv:	Raw data on butterfly eyespot patterns as measured by one of the authors using ImageJ v.1.53i.
 
- clearned_alldata.csv: Dataset including the calculated effect sizes and their variances. This dataset contains all the data used in the analysis.
 
- bird_phy.nex:	1000 phylogenetic tree data used in our study, obtained from birdtree.org (https://birdtree.org/).
 
**figs**

- fig1 - figS3.pdf are the final figures for our paper.
 
**R**

- analysis.R: R script for running meta-analysis / meta-regressions and creating figures.

- function.R:	R script for calculating effect sizes and their variances.
 
- phylogeny.R:	R script for checking whether bird phylogenetic relationships should be considered in our study.

 
**Rdata**		

- ma_50.RDS:	Results of 50 meta-analyses with 50 phylogenetic trees - calculated in "phylogeny.R".
 
**references**		

- Primary article data used in the initial and full-text screening.
