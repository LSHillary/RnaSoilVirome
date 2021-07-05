# RnaSoilVirome
This contains the R script used to analyse the read-mapping and HMMer results from the preprint "Diverse soil RNA viral communities have the potential to influence grassland ecosystems across multiple trophic levels" https://doi.org/10.1101/2021.06.28.448043. Raw sequencing reads are available from the European Nucleotide Archive (PRJEB45714). A description of the bioinformatics pipeline used in producing the input files for the data analysis hosted here is available in the methods section of the manuscript. This github repository contains the following files:

* Figure_generation.R - The R script used to analyse the data
* makeLongdb.R - A custom R function for combining multiple csv files into one data table
* mapping/ - A folder containing tsv results files from read mapping using bbmap
* RdRP_hits.txt - contigs identified as containing RdRP containing genes
* diamond_hits.txt - contigs identified as viral using Diamond (this analysis was not taken further due to using the RdRP hallmark gene to identify RNA viral contigs)
* RdRp-wolf-search-results.txt - HMMer search results
* wolf.csv - phylogenetic data from reference viruses used in Wolf et al. 2018
