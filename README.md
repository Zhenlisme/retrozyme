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

### Usage of RetrozymeSearch
```usage: RetrozymeSearch230207.py [-h] -g GENOME [-gdb BLASTNDB] [-p PROCESS] [-t {0,1}] -o OPDIR```

#### Where: 
- ```g``` mandatory argument. This is the name of a FASTA-formatted genomic sequence for the target genome. Ideally, this file will be in a folder named after the species (e.g. ```Xenopus_tropicalis/GENOME.fasta```).
- ```gdb``` optional argument. This is the name of a fasta formatted sequence to run BLASTn; it is typically the same as ```g```.
- ```p``` optional argument. This is the number of threads to use.
- ```t``` optional argument. This is to select the use of trf (false by default).
- ```OPDIR``` mandatory argument. This is the name of the output directory.

#### Example: 
```python3 RetrozymeSearch230207.py -g genomedb/Xenopus_tropicalis/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna -gdb genomedb/Xenopus_tropicalis/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna -p 32 -o Retrozymes_detection/Xenopus_tropicalis```

#### Output: 
The output directory should contain the following files and directories: 
- ```retrozyme.tbl```  is the raw output of RetrozymeSearch230207.py providing the coordinates of detected retrozymes. Columns are: Chromosome, Start, End, Repeats, Consensus_size, Strand, TRF, Name, HHR_count, Active_count, Repeat_time (active retrozyme), Consensus, HHR type,  Family, Genome occurrence
- ```repeat_summary.tbl```  is an extended table for retrozyme.tbl. It contains the RNAfold Minimum Folding Energy (MFE) value.
- ```active_retrozyme.bed```  is a table for the coordinates of active retrozyme sequences in its genome. Can be obtained from repeat_summary.tbl. In this table, the MFE value of HHR within each retrozyme should be lower than -5. Columns are: Chromosome, Start, End, Name, MFE of HHR (rnafold), Strand, Family, Repeats, Consensus_size, Span
- ```retrozyme.trf.monomer.fa``` FASTA-formatted of predicted monomer sequences.
- ```retrozyme.trf.monomer.fa.clust.info``` Clustering information following the clustering of ```retrozyme.trf.monomer.fa```.
- ```retrozyme.trf.monomer.cons.fa``` Family-level monomer consensus sequences in FASTA-format.
- ```active_family.fa```  contains the FASTA-formatted file of active retrozyme sequences.
- ```trf/```   is the tandem repeats annotation obtained using the trf program when running RetrozymeSearch230207.py with the -t option
- ```HMM_cluster/```  contains the relative position and Minimum Folding Energy (MFE) of HHRs in each retrozyme.
- ```rtztbl/```  contains the genomic coordinates of each retrozyme sorted by chromosome.
- ```Clusters/```  contains bed files obtained after clustering based on genomic positions for each retrozyme family.
- ```GenomeDB.n*``` These files will be present in the output dir if not set in the invocation; they are the blastdb files.
- ```hh1m.descr``` and ```hh1.descr``` Descriptors of hammerhead motifs searched using rnabob.
- ```genomes/``` contains the FASTA-formatted files of genomic sequences retained for the analysis.

        
### 2.	Portable_PLESearch.py  #  Script for the detection of autonomous Penelope-like elements in any given genome.
The following programs are required to run this script:
- [bedtools = 2.30.0](https://bedtools.readthedocs.io/en/stable/)
- [getorf from emboss = 6.6.0](http://emboss.open-bio.org/)
- [hmmer  = 3.3.2](http://hmmer.org/)

### Usage of Portable_PLESearch
```usage: Portable_PLESearch.py [-h] -hmm HMM -g GENOME -o OPDIR [-p PROCESS]```

#### Where: 
- ```g``` is the name of a fasta-formatted genomic sequence for the target genome.
- ```h``` is the name of the HMM to use.
- ```p``` is the number of threads.
- ```OPDIR``` is the name of the output directory.

#### Example: 
```python3 Portable_PLESearch.py -hmm ./PLE.hmm -g /home/li/genomedb/Xenopus_tropicalis/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna -o PLE_identify/Xenopus_tropicalis -p 32 >Xenopus_tropicalis.log 2>&1```

#### Output: 
The ```PLE_identify``` directory should contain the output of Portable_PLESearch.py for annotation of PLEs for each species. 
- ```PLE.tbl```  #  Table records the coordinates and classification of PLE. The columns are: Chromosome, Start, End, Coord_RT, Coord_GIY, Strand, RT_classify, GIY_classify

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
