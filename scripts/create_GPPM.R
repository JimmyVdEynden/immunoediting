##############################################################
# create_GPPM.R: create GPPM objects for all selected genes 
##############################################################

# Load data
###########
load("raw/TCGA_maf_mutect2_v7.RData") # TCGA_maf, needed to get ENSP information

# Load libraries
#################

# Internal Libraries
library(data.table)

# External bioconductor libraries (from user library, not packrat due to size)
library(ensembldb)
library(EnsDb.Hsapiens.v75) # hg19
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Load functions
#################
source("scripts/functions/create_GPosProt_from_gene.R")
source("scripts/functions/get_peptides_from_GPP.R")
source("scripts/functions/create_GPPM_from_GPP.R")

# # Select genes to use
# #####################
# 
# # All genes: Get HGNC protein coding genes
# # genes<- t(fread("../../downloads/data/HGNC/protein-coding_gene.txt",select="symbol"))
# 
# # Selected genes:
# genes<- c("PTEN","TP53")

# Create GPPM objects for selected genes
########################################

# Set location to save GPPM obejcts
GPPM_saveLocation<- "temp/GPPM/"
if(!file.exists(GPPM_saveLocation)) dir.create(GPPM_saveLocation)

# Clean errorfile
cat(file = "log/GPPM_errorfile.txt")

# Set progress bar
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
  
# Get GPPM  
for(i in 1:length(genes)){
  
  gene<- genes[i]
  
  # update progress bar
  setTxtProgressBar(pb, i)
  
  possibleError <- tryCatch({
    #   start_time<- Sys.time()
    
    # Get basic exomic information
    GPP_temp<- create_GPosProt_from_gene(gene,EnsDb.Hsapiens.v75)
    
    # Get nonapeptides
    GPP_temp$nonaPep<- get_peptides_from_GPP(GPP_temp,9)
    
    # Add gene
    GPP_temp$gene<- gene
    
    # Simulate all mutations
    GPPM_temp<- create_GPPM_from_GPP(GPP_temp,TCGA_maf)
    
    # Save GPPM 
    saveRDS(GPPM_temp,file = paste0(GPPM_saveLocation,gene,"_GPPM.rds"))
    
    # print(Sys.time()-start_time)
  }, error=function(e) e)
  if(inherits(possibleError, "error")){
    cat(gene,"\n",file = "log/GPPM_errorfile.txt",append = TRUE)
    next
  } 
}
