################################
# map_IEDB_epitopes_to_hg19.R
################################

# Download IEDB epitopes fro synapse (https://www.synapse.org/; Zapata et al., id syn11935058)
#############################################################################################
epi<- readRDS("data/IEDB_epitopes.rds")

# Get transcript information and amino acid coordinates
########################################################
epi_ENST<- epi[,3]
epi_start<- as.numeric(gsub("\\-.*","",epi[,6]))
epi_end<- as.numeric(gsub(".*\\-","",epi[,6]))

# Map these coordinates to hg 19
################################

# Need latest version of iedb (2.4.1) to ise a nice function to do this (oroteinToGenome)!
# See https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html#5_mapping_protein_coordinates_to_the_genome
# Only worked after updating R, Bioconductor 3.7 is used. 
# Done seperately as all other analysis done easing earlier version
# source("https://bioconductor.org/biocLite.R")
# biocLite("ensembldb")
library(ensembldb) # This loads 2.4.1 which does seem to work!
library(EnsDb.Hsapiens.v75)
edb75<- EnsDb.Hsapiens.v75

# Put in IRanges object
epi_prt <- IRanges(start = epi_start, end = epi_end, names = epi_ENST)

# Map to GRanges: takes severall hours (+/-7hr)!
# Use loop as this avoids errors for unretrieved CDS etc
start_time<- Sys.time()
epi_grm <- suppressMessages(proteinToGenome(epi_prt[1], edb75, idType = "tx_id"))
for(i in (length(epi_grm)+1):(length(epi_prt))){
  cat(i," ")
  epi_grm_temp <- suppressMessages(proteinToGenome(epi_prt[i], edb75, idType = "tx_id"))
  epi_grm<- c(epi_grm,epi_grm_temp) 
  if(i%%1000==0) saveRDS(epi_grm,file="temp/epi_grm_temp.rds")
}
Sys.time()-start_time
saveRDS(epi_grm,file="temp/epi_grm.rds")

# Fuse data in one GRanges file
for(i in 1:length(epi_grm)){
  cat(i," ")
  if(i==1) epi_grm_fused<- epi_grm[[i]]
  else epi_grm_fused<- c(epi_grm_fused,epi_grm[[i]])
}
epi_gr<- epi_grm_fused
epi_gr<- sortSeqlevels(epi_gr) # Sort chromosomes
epi_gr<- sort(epi_gr) # Sort GRanges
saveRDS(epi_gr,file="data/IEDB_epitopes_gr.rds")
