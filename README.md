# RETROZYME readme

This repository contains data files and scripts for annotating retrozyme elements in a given genomic sequence.

## Main analysis scripts 
The main analysis scripts are:

### 1.	RetrozymeSearch230207.py  #  Script for detection of retrozymes in any given genome.
The following programs are required to run this script:
- [bedtools = 2.30.0](https://bedtools.readthedocs.io/en/stable/)
- [blastn = 2.2.31](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [rnabob = 2.2.1](http://eddylab.org/software/rnabob/rnabob.tar.gz)
- [RNAfold = 2.6.4](https://www.tbi.univie.ac.at/RNA/)
- [RNAeval = 2.6.4](https://www.tbi.univie.ac.at/RNA/)
- [trf = 4.09](https://tandem.bu.edu/trf/trf.html)
- [vsearch = 2.11.0](https://github.com/torognes/vsearch)
        
### 2.	Portable_PLESearch.py  #  Script for the detection of autonomous Penelope-like elements in any given genome.
The following programs are required to run this script:
- [bedtools = 2.30.0](https://bedtools.readthedocs.io/en/stable/)
- [getorf from emboss = 6.6.0](http://emboss.open-bio.org/)
- [hmmer  = 3.3.2](http://hmmer.org/)

### General dependencies are: 
  - python = 3.9.0
  - r-base = 4.1
  - biopython
  - pybedtools = 0.9.0 =py39hd65a603_2
  - r-bedtoolsr
  - r-seqinr = 4.2_16 = r41h06615bd_0
  - genometools-genometools = 1.6.2 = py39h58cc16e_6

When running these scripts to analyse several genomes, it is convenient to use a wrapper such as:
- ```run_plesearch.sh``` # Wrapper around RetrozymeSearch230207.py
- ```run_retrozyme.sh``` # Wrapper around Portable_PLESearch.py


## Downstream analysis scripts

The following scripts are useful for summarising and analysing various features of the genomic distribution of retrozymes.

### Retrozyme_feature_summary.R
R script to summarise retrozyme features

### RetrozymeFeature.R
R script to visualise the retrozyme distribution pattern in each species

### retrozyme_inchrm.R
R script to visualise the retrozyme distribution on each chromosome genome

### retrozyme_coord_correlation_with_genes.R
Overlaps between retrozymes and genomic features

### PLE_distribution_across_vert.R
Distribution of PLEs in each species’ genome

### RNAfoldEnergy.R
Compute minimum energy of predicted folds and plot as a function of retrozyme family abundance (see rnafold_appearence.png).

### Retrozyme_expression.R
Compute and plot depth of coverage of either ESTs or piRNAs for selected retrozyme elements.

### Retrozyme_feature_taxonomics.R
Compute and plot several features of retrozyme elements across genomes following taxonomy: active counts, monomer size, Maximum of repeated monomers ... in retrozyme arrays (see RepeatLimits_range.png Feature_accross_trees.png).

### Retrozyme_tree_with_TE.R
Plot retrozyme and PLE tree.

### Rtz_353_331_Ax.R
Computes minimum energy of predicted folds for Rtz331_353_Ax sequences and plots an histogram (see Rtz331_353_Ax_rnafold.png).

### genomesize.sh
Compute and sort genomic sequence by size.

### retrozyme_coord_correlation_with_gene.R
Compute overlap between Xenopus laevis and Xenopus tropicalis retrozymes and gene features: exons, introns, intergenic (see Genomic_region_proportion.png).
