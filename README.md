# Systematic review and meta-analysis of anti-predator mechanism of eyespots: conspicuous pattern vs eye mimicry (version 1.0)

This repository contains the all code and data of **_Systematic review and meta-analysis of eyespot anti-predator mechanisms: conspicuous patterns vs. eye mimicry_**. Our meta-analytical study examined why butterfly eyespots intimidate birds and tested two hypotheses: the eye mimicry hypothesis and the conspicuousness hypothesis. Despite conflicting empirical studies, our meta-analysis supported the idea of the conspicuousness hypothesis. We also found a key factor—larger eyespots were more effective in predator avoidance. Our findings suggest that birds avoid eyespots mainly due to their conspicuousness and emphasise the importance of such pattern characteristics in predator avoidance.

HTML file (https://ayumi-495.github.io/eyespot/) contains the main supplementary material within it, including main CSV files and R scripts. You can download them directly from there. 

## Description of each folder and file:	

**data**

- all_21092023.csv: Raw data including both predator and prey datasets. This dataset was used for calculating effect sizes.
 
- butterfly_eyespots.csv:	Raw data on butterfly eyespot patterns as measured by one of the authors using ImageJ v.1.53i.
 
- clearned_alldata.csv: Dataset including the calculated effect sizes and their variances. This dataset contains all the data used in the analysis.
 
- bird_phy.nex:	1000 phylogenetic trees data used in our study, obtained from birdtree.org (https://birdtree.org/).
 
**figs**

- fig1 - figS3.pdf are the final figures for our paper.
 
**R**	

- function.R:	R script to calculate effect sizes and their variances.
 
- phylogeny.R:	R script to check whether bird phylogenetic relationships should be considered in our study.
 
**Rdata**		

- ma_50.RDS:	Results of 50 meta-analyses with 50 phylogenetic trees - calculated in "phylogeny.R", then used.
 
**references**		

- Primary article data used in the initial and full-text screening.
