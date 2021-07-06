#### 0 - Clear workspace and load dependencies ####
# 0.1 - Clear workspace
rm(list=ls())
# 0.2 - Load packages
if(!require(ggsci)){install.packages('ggsci'); require(ggsci)}
if(!require(pheatmap)){install.packages('pheatmap'); require(pheatmap)}
if(!require(vegan)){install.packages('vegan'); require(vegan)}
if(!require(factoextra)){install.packages('factoextra'); require(factoextra)}
if(!require(doBy)){install.packages('doBy'); require(doBy)}
if(!require(rhmmer)){install.packages('rhmmer'); require(rhmmer)}
if(!require(tidyverse)){install.packages('tidyverse'); require(tidyverse)}
library(rhmmer)
library(dplyr)
library(tidyr)
library(stringr)
library(conflicted)
library(data.table)
library(UpSetR)
library(vegan)
library(RColorBrewer)
library(seqinr)

# Fix conflict
conflict_prefer("select", "dplyr")

# 0.3 - Set working directory - note, this needs altering if you wish to run the script yourself!
setwd("~/Dropbox/PhD/Chapter 4 - Soil RNA Virome/R/RdRp+Diamond/wolf/wolf300")
#setwd("C:/Users/bspa46/Dropbox/PhD/Chapter 4 - Soil RNA Virome/R/RdRp+Diamond/wolf/wolf300")

# 0.4 - Load/ define functions
# Change location to suit your own setup
source('../../../functions/makeLongDB.R')

# Comb determines if a vOTU is present in a sample, and then combines the results for each replicate within a sample type
Comb <- function(X, J){
  # J - data table
  # X - AF cutoff
  J$Presence <- ifelse(J$Covered_percent < X, 0, 1)
  K <- data.table(J$vOTU, J$Sample, J$Presence)
  names(K) <- c("vOTU", "Sample", "Presence")
  K <- spread(K, key = Sample, value = Presence)
  K$`Upland-peatland` <- ifelse(rowSums(cbind(K$`Upland-peatland 1`,
                                       K$`Upland-peatland 2`,
                                       K$`Upland-peatland 3`))
                         > 0, 1, 0)
  K$`Upland grassland` <- ifelse(rowSums(cbind(K$`Upland-grassland 1`,
                                             K$`Upland-grassland 2`,
                                             K$`Upland-grassland 3`))
                               > 0, 1, 0)
  K$`Semi-improved` <- ifelse(rowSums(cbind(K$`Semi-improved 1`,
                                            K$`Semi-improved 2`,
                                            K$`Semi-improved 3`))
                              > 0, 1, 0)
  K$`Lowland grassland` <- ifelse(rowSums(cbind(K$`Lowland-grassland 1`,
                                              K$`Lowland-grassland 2`,
                                              K$`Lowland-grassland 3`))
                                > 0, 1, 0)
  K$`Coastal-grassland` <- ifelse(rowSums(cbind(K$`Coastal-grassland 1`,
                                      K$`Coastal-grassland 2`,
                                      K$`Coastal-grassland 3`))
                        > 0, 1, 0)
  L <- data.table(K$vOTU, K$`Upland-peatland`, K$`Upland grassland`, K$`Semi-improved`, K$`Lowland grassland`, K$`Coastal-grassland`)
  names(L) <- c("vOTU", "Upland-peatland", "Upland-grassland", "Semi-improved", "Lowland-grassland", "Coastal-grassland")
  L <- as.matrix(L[,-1], rownames = L$vOTU)
  L <- L[rowSums(L)>0,]
  return(L)
}

# fpkm2tpm takes fpkm values and converts to tpm values
fpkm2tpm <- function(fpkm){
  tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
  tpm[which(is.na(tpm))] <- 0
  return(tpm)}

# presence takes a data table of contigs and mapping data, determines which contigs are present, and reshapes it for use in creating an UpSet plot
presence <- function(dt,X){
  dt$Presence <- ifelse(dt$Covered_percent < X, 0, 1)
  dt <- select(dt,c("vOTU", "Presence","Sample"))
  dt <- spread(dt, Sample, Presence)
  row.names(dt) <- dt$vOTU
  dt <-dt[,c(2,3,4)]
  dt <- dt[rowSums(dt)>0,]
  return(dt)
}

# RemoveContaminants removes all contigs that had any read map from either of the negative libraries
RemoveContaminants <- function(dt){
  dt_Negative_Extraction <- dt[dt$Sample == 'Negative-Extraction',]
  Negative_Extraction_Pass <- dt_Negative_Extraction[dt_Negative_Extraction$Reads == 0,]$vOTU
  dt2 <- dt[dt$vOTU %in% Negative_Extraction_Pass,]
  dt2 <- dt2[dt2$Sample != 'Negative-Extraction',]
  dt_Negagive_Library <- dt[dt$Sample == 'Negative-Library',]
  Negative_Library_Pass <- dt_Negagive_Library[dt_Negagive_Library$Reads == 0,]$vOTU
  dt3 <- dt2[dt2$vOTU %in% Negative_Library_Pass,]
  dt3 <- dt3[dt3$Sample != 'Negative-Library',]
  return(dt3)
}
#### 1 - Load coverage/ RPPKM data and create subsets ####

# 1.1 - Load coverage/ RPKM data for each sample into one large data table
CovRaw <- makeLongDB(path = "mapping/", pattern = "-coverage.tsv", chopFromFileName = "-coverage.tsv", header = TRUE)
RpkmRaw <- makeLongDB(path = "mapping/", pattern = "-rpkm.tsv", chopFromFileName = "-rpkm.tsv", header = TRUE)

# 1.2 - Change name of vOTU columns to "vOTU"

colnames(RpkmRaw)[colnames(RpkmRaw)=="X.Name"] <- "vOTU"
colnames(CovRaw)[colnames(CovRaw)=="X.ID"] <- "vOTU"
TableRaw <- full_join(CovRaw, RpkmRaw, by = c("Sample", "vOTU", "Length"))

# 1.2 - Create list of data tables for each vOTU where if a vOTU has no mapped reads, a null data table is returned
ListRaw <- by(TableRaw, TableRaw$vOTU, FUN = function(x){if(sum(x$Covered_percent) != 0){return(x)}})

# 1.3 - Create vector of vOTUs recording if the data table is empty
ListRawCheck <- unlist(lapply(ListRaw, FUN = function(x){ifelse(!is.null(x), yes = TRUE, no = FALSE)}))

# 1.4 - Filter list of data frames based on if the data table is null and convert to a single data table
TableMapped <- do.call(rbind,ListRaw[ListRawCheck == TRUE])
TableMapped$vOTU <- word(TableMapped$vOTU, 1)
TableMapped$Sample <- as.character(TableMapped$Sample)

# 1.5 - Rename samples as something meaningful and set factors
TableMapped[TableMapped == "A1"] <- "Upland-peatland 1"
TableMapped[TableMapped == "A2"] <- "Upland-peatland 2"
TableMapped[TableMapped == "A3"] <- "Upland-peatland 3"
TableMapped[TableMapped == "B1"] <- "Upland-grassland 1"
TableMapped[TableMapped == "B2"] <- "Upland-grassland 2"
TableMapped[TableMapped == "B3"] <- "Upland-grassland 3"
TableMapped[TableMapped == "C1"] <- "Semi-improved 1"
TableMapped[TableMapped == "C2"] <- "Semi-improved 2"
TableMapped[TableMapped == "C3"] <- "Semi-improved 3"
TableMapped[TableMapped == "D1"] <- "Lowland-grassland 1"
TableMapped[TableMapped == "D2"] <- "Lowland-grassland 2"
TableMapped[TableMapped == "D3"] <- "Lowland-grassland 3"
TableMapped[TableMapped == "E1"] <- "Coastal-grassland 1"
TableMapped[TableMapped == "E2"] <- "Coastal-grassland 2"
TableMapped[TableMapped == "E3"] <- "Coastal-grassland 3"
TableMapped[TableMapped == "N1"] <- "Negative-Extraction"
TableMapped[TableMapped == "N2"] <- "Negative-Library"

TableMapped$Sample <- factor(TableMapped$Sample,levels = 
                               c("Negative-Library", "Negative-Extraction", "Coastal-grassland 3", "Coastal-grassland 2", "Coastal-grassland 1",
                                 "Lowland-grassland 3", "Lowland-grassland 2", "Lowland-grassland 1",
                                 "Semi-improved 3", "Semi-improved 2", "Semi-improved 1",
                                 "Upland-grassland 3", "Upland-grassland 2", "Upland-grassland 1",
                                 "Upland-peatland 3", "Upland-peatland 2", "Upland-peatland 1"))

#### 2 - Identify Diamond and HMMer hits and add information to the data table ####

# 2.1 - Read in results files
rdrp_hits <- fread("RdRp_hits.txt", header = F)

# 2.2 - Add results to the data table
TableMapped$RdRp <- F
TableMapped$RdRp[TableMapped$vOTU %in% rdrp_hits$V1] <- T

# 2.3 - Run decon function
decon <- RemoveContaminants(TableMapped)

#### 3  - Identify contigs that are present in a sample and plot Upset figure ####

# 3.1 - Identify contigs that are present and subset data table to only contain those that are
# samples <- unique(decon$Sample)
decon$Presence <- TRUE
decon[decon$Covered_percent <50,]$Presence <- FALSE
RdRpTable <- decon[decon$RdRp == T,]

# 3.2 - Produce UpSet diagrams at different coverage thresholds
Comb25 <- Comb(25, RdRpTable)
Comb50 <- Comb(50, RdRpTable)
Comb75 <- Comb(75, RdRpTable)
Comb95 <- Comb(95, RdRpTable)

dtComb25 <- as.data.table(Comb25)
dtComb50 <- as.data.table(Comb50)
dtComb75 <- as.data.table(Comb75)
dtComb95 <- as.data.table(Comb95)

dtComb25$vOTU <- row.names(Comb25)
dtComb50$vOTU <- row.names(Comb50)
dtComb75$vOTU <- row.names(Comb75)
dtComb95$vOTU <- row.names(Comb95)

U25 <- upset(dtComb25, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

U50 <- upset(dtComb50, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

U75 <- upset(dtComb75, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

U95 <- upset(dtComb95, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

U50

# 3.3 - Produce UpSet diagrams of the triplicates from each sample type at 50% coverage
PL <- c("Upland-peatland 1", "Upland-peatland 2", "Upland-peatland 3")
UG <- c("Upland-grassland 1", "Upland-grassland 2", "Upland-grassland 3")
SI <- c("Semi-improved 1", "Semi-improved 2", "Semi-improved 3")
LG <- c("Lowland-grassland 1", "Lowland-grassland 2", "Lowland-grassland 3")
CS <- c("Coastal-grassland 1", "Coastal-grassland 2", "Coastal-grassland 3")

Upland <- RdRpTable[RdRpTable$Sample %in% UG & RdRpTable$Covered_percent,]
SemiImproved <- RdRpTable[RdRpTable$Sample %in% SI,]
Lowland <- RdRpTable[RdRpTable$Sample %in% LG,]
`Upland-peatland` <- RdRpTable[RdRpTable$Sample %in% PL,]
`Coastal-grassland` <- RdRpTable[RdRpTable$Sample %in% CS,]

U50 <- presence(Upland, 50)
P50 <- presence(`Upland-peatland`, 50)
L50 <- presence(Lowland, 50)
S50 <- presence(SemiImproved, 50)
C50 <- presence(`Coastal-grassland`, 50)

upset(U50, sets = c("Upland-grassland 1","Upland-grassland 2", "Upland-grassland 3"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

upset(P50, sets = c("Upland-peatland 1","Upland-peatland 2", "Upland-peatland 3"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

upset(L50, sets = c("Lowland-grassland 1","Lowland-grassland 2", "Lowland-grassland 3"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

upset(S50, sets = c("Semi-improved 1","Semi-improved 2", "Semi-improved 3"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

upset(C50, sets = c("Coastal-grassland 1","Coastal-grassland 2", "Coastal-grassland 3"), order.by = "freq", keep.order = TRUE,
      mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
      text.scale = c(1.5,1.5,1.5,1.5,2,1.5))

#### 4 - Add virus type to data table ####
# 4.1 - Read in HMMer results and filter for score (results already filtered for e value when HMMer was run)
HMMER_results <- read_tblout('RdRp-wolf-search-results.txt')
HMMER_filtered <- HMMER_results[HMMER_results$sequence_score>=50,]

# 4.2 - Extract virus types
retrieve_contig <- function(name){
  text <- as.character(name)
  group <- unlist(strsplit(text,'_'))
  contig_name <- paste(group[1], group[2], sep = "_")
  return(contig_name)
}

# 4.3 - Reformat contig names
HMMER_filtered$contig_names <- sapply(HMMER_filtered$domain_name, retrieve_contig)

# 4.4 - Create a data table of results which removes dublicate entries 
ht <- as.data.table(HMMER_filtered %>% select(contig_names, query_name, sequence_evalue))
ht2 <- ht[order(ht$contig_names, abs(ht$sequence_evalue) ), ]
ht3 <- ht2[ !duplicated(ht2$contig_names), ]

# 4.5 - creates lists of contig names by virus type
RdRpTable$type = FALSE
b1 <- ht3$contig_names[ht3$query_name == "RNAvirome.S4A"]
b2 <- ht3$contig_names[ht3$query_name == "RNAvirome.S4B"]
b3 <- ht3$contig_names[ht3$query_name == "RNAvirome.S4C"]
b4 <- ht3$contig_names[ht3$query_name == "RNAvirome.S4D"]
b5 <- ht3$contig_names[ht3$query_name == "RNAvirome.S4E"]

# 4.6 - Add information to data table

RdRpTable$type[RdRpTable$vOTU %in% b1] <- "Branch 1"
RdRpTable$type[RdRpTable$vOTU %in% b2] <- "Branch 2"
RdRpTable$type[RdRpTable$vOTU %in% b3] <- "Branch 3"
RdRpTable$type[RdRpTable$vOTU %in% b4] <- "Branch 4"
RdRpTable$type[RdRpTable$vOTU %in% b5] <- "Branch 5"

# 4.7 - Removes information from data table for viruses that are not present within that sample
RdRpTable$type[RdRpTable$Covered_percent < 50] <- NA

# #### 6 - Presence/ Absence analysis ####
# RdRp_summary <- summaryBy(Presence ~ Sample + type, data = RdRpTable, FUN = sum)
# RdRp_summary <- spread(RdRp_summary, type, Presence.sum)
# community <- RdRp_summary[,2:6]
# row.names(community) <- RdRp_summary[,1]
# community[is.na(community)] <- 0
# community_nmds <- metaMDS(community)
# stressplot(community_nmds)
# plot(community_nmds)
# 
# site=c(rep("Upland-peatland",3),rep("Upland grassland",3),rep("Semi-improved grassland",3),
#        rep("Lowland grassland",3), rep("Coastal-grassland",3))
# #ordiplot(community_nmds,type="n")
# ordihull(community_nmds,groups=site,draw="polygon",col=c("#451A54","#3B528B","#FEE833","#3C908C","#5FC863"),label=F)
# orditorp(community_nmds,display="species",col="red",air=0.1)
# orditorp(community_nmds,display="sites",col=c(rep("green",5),rep("blue",5)),
#          air=0.01,cex=1.25)
# 
# community4  <- data.frame(RdRp_summary[,2:5])
# row.names(community4) <- community[,1]
# community4_nmds <- metaMDS(community4)
# stressplot(community4_nmds)
# plot(community4_nmds)
# ordihull(community4_nmds,groups=site,draw="polygon",col=c("#451A54","#3B528B","#FEE833","#3C908C","#5FC863"),label=F)
# orditorp(community4_nmds,display="species",col="red",air=0.1)
# orditorp(community4_nmds,display="sites",col=c(rep("green",5),rep("blue",5)),
#          air=0.01,cex=1.25)

#### 5 - TPM analysis ####

# 5.1 - Calculates tpm values in each sample for each contig that is regarded as present
samples <- unique(decon$Sample)
RdRpTable$tpm <- NA
RdRpTable$FPKM[RdRpTable$Presence == F] <- NA
for (sample in samples){
  RdRpTable$tpm[RdRpTable$Sample == sample] <- fpkm2tpm(RdRpTable$FPKM[RdRpTable$Sample == sample])
}

#RdRpTable[RdRpTable$Sample == samples[1],]

#RdRpTable$tpm <- fpkm2tpm(RdRpTable$FPKM)
#RdRpTable$tpm[RdRpTable$Presence == F] <- 0
#RdRp_tpm_summary <- spread(RdRp_tpm_summary, type, tpm.sum)

# 5.2 - Creates a summary table for each virus type and each sample
RdRp_tpm_summary <- summaryBy(tpm ~ Sample + type, data = RdRpTable, FUN = sum)
RdRp_tpm_summary$Percentage <- RdRp_tpm_summary$tpm.sum/10000
RdRp_tpm_summary <- RdRp_tpm_summary[is.na(RdRp_tpm_summary$type) == F,]

# 5.3 - Renames samples and sets factors
RdRp_tpm_summary <- droplevels(RdRp_tpm_summary)
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Coastal-grassland 3"] <- "Coastal-grassland 3"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Coastal-grassland 2"] <- "Coastal-grassland 2"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Coastal-grassland 1"] <- "Coastal-grassland 1"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Lowland-grassland 3"] <- "Lowland grassland 3"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Lowland-grassland 2"] <- "Lowland grassland 2"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Lowland-grassland 1"] <- "Lowland grassland 1"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Semi-improved 3"] <- "Semi-improved grassland 3"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Semi-improved 2"] <- "Semi-improved grassland 2"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Semi-improved 1"] <- "Semi-improved grassland 1"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-grassland 3"] <- "Upland grassland 3"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-grassland 2"] <- "Upland grassland 2"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-grassland 1"] <- "Upland grassland 1"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-peatland 3"] <- "Upland peatland 3"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-peatland 2"] <- "Upland peatland 2"
RdRp_tpm_summary$name[RdRp_tpm_summary$Sample == "Upland-peatland 1"] <- "Upland peatland 1"

RdRp_tpm_summary$name <- factor(RdRp_tpm_summary$name, levels = c(
                                "Upland peatland 1", "Upland peatland 2", "Upland peatland 3",
                                "Upland grassland 1", "Upland grassland 2", "Upland grassland 3",
                                "Semi-improved grassland 1", "Semi-improved grassland 2", "Semi-improved grassland 3",
                                "Lowland grassland 1", "Lowland grassland 2", "Lowland grassland 3",
                                "Coastal-grassland 1", "Coastal-grassland 2", "Coastal-grassland 3"))

# 5.4 - plots figure
ggplot(RdRp_tpm_summary, aes(fill = type, y = Percentage, x = name)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance (%)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.position = "none")

#### 6 - Calculating diversity stats for TPM ####
# 6.1 - Creates a data table suitable for input for NMDS
CommunityTable <- select(RdRpTable, c(Sample, vOTU, tpm))
CommunityTable <- spread(CommunityTable, vOTU, tpm)
rownames(CommunityTable) <- CommunityTable$Sample
CommunityTable <- select(CommunityTable, -c(Sample))

# 6.2 - Generates NMDS object and plots using ggplot
CommunityTable_nmds <- metaMDS(CommunityTable, distance = "bray")

data.scores <- as.data.frame(scores(CommunityTable_nmds))
data.scores$site <- rownames(data.scores)
data.scores$grp <- c(rep("Coastal-grassland",3),rep("Lowland grassland",3),rep("Semi-improved grassland",3),
                     rep("Upland grassland",3), rep("Upland-peatland",3)) 
species.scores <- as.data.frame(scores(CommunityTable_nmds, "species"))
species.scores$species <- rownames(species.scores)

ggplot() + 
  geom_polygon(data=data.scores[data.scores$grp == "Coastal-grassland",],aes(x=NMDS1,y=NMDS2), fill = "#451A54", alpha = 0.1) +
  geom_polygon(data=data.scores[data.scores$grp == "Lowland grassland",],aes(x=NMDS1,y=NMDS2, fill = "#3B528B", alpha = 0.1)) +
  geom_polygon(data=data.scores[data.scores$grp == "Semi-improved grassland",],aes(x=NMDS1,y=NMDS2, fill = "#FEE833", alpha = 0.1)) +
  geom_polygon(data=data.scores[data.scores$grp == "Upland grassland",],aes(x=NMDS1,y=NMDS2, fill = "#3C908C", alpha = 0.1)) +
  geom_polygon(data=data.scores[data.scores$grp == "Upland-peatland",],aes(x=NMDS1,y=NMDS2, fill = "#5FC863", alpha = 0.1)) +
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site, colour=as.factor(site)), vjust=0) +  # add the site labels
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2),shape = 17, size = 1, fill = "black", alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, colour=grp)) + # add the point markers
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

#### 7 - Phylogeny ####
# retrieve_info 
retrieve_info <- function(branch, study){
  fasta <- read.fasta(paste(branch,"_",study,".faa", sep = ""))
  seqs <- names(fasta)
  dt <- data.table(seqs,branch,study)
  return(dt)
}

branches <- c("b1","b2","b3","b4","b5")

phylo_dt <- data.table(seqs = character(),branch = character(),study = character())

for(Branch in branches){
  starr_dt <- retrieve_info(Branch,"starr")
  wolf_dt <- retrieve_info(Branch,"wolf")
  bangor_dt <- retrieve_info(Branch, "seqs")
  phylo_dt <- rbind(phylo_dt, bangor_dt, starr_dt, wolf_dt)
}

names(phylo_dt) <- c("RdRp","branch","study")


# 8.2 - Import phylogenetic data on Wolf et al viruses
wolf_phylo <- data.frame(read.csv("wolf.csv"))

# source function for generating itol annotation tables
source("table2itol/table2itol.R")

# 8.3 - Create master annotation table
temp <- merge(phylo_dt, wolf_phylo, all.x = TRUE, by.x = "RdRp", by.y = "RdRp.GenBank.Acc")

# 8.4 - Generate branch specific annotation files (note that csv outputs were converted to .xlsx format by importing into Microsoft Excel)
### Branch 1 ###
b1 <- temp[temp$branch == "b1"]
write.csv(b1, "b1_dt.csv", row.names = FALSE)
create_itol_files(infiles = "b1_dt.xlsx",
                  identifier = "Label", na.strings = "X", width = 1)

### Branch 2 ###
b2 <- temp[temp$branch == "b2"]
write.csv(b2, "b2_dt.csv", row.names = FALSE)
create_itol_files(infiles = "b2_dt.xlsx",
                  identifier = "Label", na.strings = "X", width = 1)

### Branch 3 ###
b3 <- temp[temp$branch == "b3"]
write.csv(b3, "b3_dt.csv", row.names = FALSE)
create_itol_files(infiles = "b3_dt.xlsx",
                  identifier = "RdRp", na.strings = "X", width = 1)

### Branch 4 ###
b4 <- temp[temp$branch == "b4"]
write.csv(b4, "b4_dt.csv", row.names = FALSE)
create_itol_files(infiles = "b4_dt.xlsx",
                  identifier = "RdRp", na.strings = "X", width = 1)

### Branch 5 ###
b5 <- temp[temp$branch == "b5"]
write.csv(b5, "b5_dt.csv", row.names = FALSE)
create_itol_files(infiles = "b5_dt.xlsx",
                  identifier = "RdRp", na.strings = "X", width = 1)

##  8.5 - producing UpSet plots for supplementary information

noda_list <- c("k127_2771575", "k127_542449", "k127_3249532", "k127_2458157", "k127_931303", "k127_1820739", "k127_3393807", "k127_724126", "k127_1403543", "k127_1044899", "k127_2442395", "k127_1031125", "k127_2616941", "k127_1669650", "k127_3193470", "k127_2831791", "k127_2850814", "k127_2466107", "k127_384702", "k127_240129", "k127_3291927", "k127_1468481", "k127_1134143", "k127_1844969", "k127_2169659")
Noda <- RdRpTable[RdRpTable$vOTU %in% noda_list,]
CombNoda <- Comb(50, Noda)
dtCombNoda <- as.data.table(CombNoda)
dtCombNoda$vOTU <- row.names(CombNoda)

UNoda <- upset(dtCombNoda, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
             mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
             text.scale = c(1.5,1.5,1.5,1.5,2,1.5))
UNoda

dicistro_list <- c("k127_3303538", "k127_2168122", "k127_1529838", "k127_1404685", "k127_767399", "k127_2401058", "k127_2410764", "k127_2643016", "k127_2385373", "k127_1543493", "k127_898091", "k127_2754849", "k127_1846526", "k127_179655", "k127_726481", "k127_1475360", "k127_2249241", "k127_2613944", "k127_436183", "k127_2936501", "k127_1377121", "k127_148344", "k127_3191666", "k127_1452805", "k127_1069181", "k127_639836", "k127_2338401", "k127_2949364", "k127_2847127", "k127_1386140", "k127_2647977", "k127_1439664", "k127_2467920", "k127_3136826", "k127_2079389", "k127_2073162", "k127_3117831", "k127_3229607", "k127_3280409", "k127_3397315", "k127_2639070", "k127_2018362", "k127_2771575", "k127_542449", "k127_3249532", "k127_2458157", "k127_931303", "k127_1820739", "k127_3393807", "k127_724126", "k127_1403543", "k127_1044899", "k127_2442395", "k127_1031125", "k127_2616941", "k127_1669650", "k127_3193470", "k127_2831791", "k127_2850814", "k127_2466107", "k127_384702", "k127_240129", "k127_3291927", "k127_1468481", "k127_1134143", "k127_1844969", "k127_2169659")
dicistro <- RdRpTable[RdRpTable$vOTU %in% dicistro_list,]
Combdicistro <- Comb(50, dicistro)
dtCombdicistro <- as.data.table(Combdicistro)
dtCombdicistro$vOTU <- row.names(Combdicistro)

Udicistro <- upset(dtCombdicistro, sets = c("Coastal-grassland","Lowland-grassland", "Semi-improved", "Upland-grassland","Upland-peatland"), order.by = "freq", keep.order = TRUE,
               mainbar.y.label = "Sample Intersections", sets.x.label = "Viral contigs",
               text.scale = c(1.5,1.5,1.5,1.5,2,1.5))
Udicistro