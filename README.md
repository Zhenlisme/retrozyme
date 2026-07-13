RETROZYME readme

This repository contains data files and scripts for annotating retrozyme elements in a given genomic sequence.
The main analysis scripts are:
1.	RetrozymeSearch230207.py  #  Script for detection of retrozymes in any given genome. 
2.	Portable_PLESearch.py  #  Script for the detection of autonomous Penelope-like elements in any given genome.

When running these scripts to analyse several genomes, it is convenient to use a wrapper such as:
3.	run_plesearch.sh # Wrapper around RetrozymeSearch230207.py
4.	run_retrozyme.sh # Wrapper around Portable_PLESearch.py

The following scripts are useful to summarise and analyse various features on the genomic distribution of retrozymes.

1.	Retrozyme_feature_summary.R  #  R script to summarise retrozymes features 
2.	RetrozymeFeature.R  # R script to visualise the retrozyme distribution pattern in each species
3.	retrozyme_inchrm.R  #  R script to visualise the retrozyme distribution on each chromosome genome
4.	retrozyme_coord_correlation_with_genes.R  #  Overlaps between retrozymes and genomic features
5.	PLE_distribution_across_vert.R  #  Distribution of PLEs in each species’ genome
6.	tree.nwk  #  Time tree for sampled vertebrate species
7.	RNAfoldEnergy.R
8.	Retrozyme_expression.R
9.	Retrozyme_feature_taxonomics.R
10.	Retrozyme_tree_with_TE.R
11.	Rtz_353_331_Ax.R
12.	genomesize.sh
13.	retrozyme_coord_correlation_with_gene.R
