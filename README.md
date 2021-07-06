# RnaSoilVirome
This contains the R script used to analyse the read-mapping and HMMer results from the preprint "Diverse soil RNA viral communities have the potential to influence grassland ecosystems across multiple trophic levels" https://doi.org/10.1101/2021.06.28.448043. Raw sequencing reads are available from the European Nucleotide Archive (PRJEB45714). A description of the bioinformatics pipeline used in producing the input files for the data analysis hosted here is available in the methods section of the manuscript. This github repository contains the following files:

* elevation profile.csv - elevation data for producing the elevation profile
* elevation_profile.R - R script for generating the elevation profile
* makeLongdb.R - A custom R function for combining multiple csv files into one data table
* RdRP_hits.txt - contigs identified as containing RdRP containing genes
* RdRp-wolf-search-results.txt - HMMer search results
* RnaVirome_FigureGeneration.R - The R script used to analyse the data
* wolf.csv - phylogenetic data from reference viruses used in Wolf et al. 2018

Seperate folders also contain additional files used in the analyses:
* faa_files/ - A folder containing faa files used to generate multiple and also necessary for generating itol annotations
* mapping/ - A folder containing tsv results files from read mapping using bbmap
* multiple_sequence_alignments - A folder containing multiple sequence alignments used for building phylogenetic trees
* phylogenetic_trees - A folder containing phylogenetic trees and their itol annotation files
