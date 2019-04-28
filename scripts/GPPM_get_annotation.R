###################################
# manuscript_HLA_annotation.R
###################################

# Load GPPM
############
{
  GPPM<- readRDS(file = "data/GPPM_inclHLAAlleles.rds")
}

# Filter GPPM
##############
{
  # No need for mutation info for this analysis!
  idx_mut<- which(!duplicated(paste(GPPM$gene,pos(GPPM),sep="_")))
  GPPM<- GPPM[idx_mut] 
}

# Create annotation bed file for overall binding (harmonic mean) and different alleles
#######################################################################################
{
  # Granges for annotated genome
  GPPM_annotated<- GPPM
  mcols(GPPM_annotated)<- NULL
  GPPM_annotated<- sortSeqlevels(GPPM_annotated) # Sort chromosomes
  GPPM_annotated<- sort(GPPM_annotated) # Sort
  GPPM_annotated<- union(GPPM_annotated,GPPM_annotated) # Convert in genomic ranges
  annotation_gr_all<- GPPM_annotated
  
  # GRanges for HLA binding regions
  for(i in 1:7){
    
    # cat(i," ")
    # First GRanges file with binding regions
    idx_allele<- grep("HLA_aff",names(mcols(GPPM)))[i] # A1, A2, B1, B2, C1, C2, mean
    allele_name<- c("HLA_A02_01","HLA_A01_01","HLA_B07_02","HLA_B08_01","HLA_C07_01","HLA_C07_02","HLA_Mean")[i] # Cfr read.table("temp/HLA_types_freq.txt",stringsAsFactors = F)
    GPPM_binding<- GPPM[mcols(GPPM)[,idx_allele]<500]
    
    # Turn into Granges
    GPPM_binding<- sortSeqlevels(GPPM_binding) # Sort chromosomes
    GPPM_binding<- sort(GPPM_binding) # Sort
    GPPM_binding<- union(GPPM_binding,GPPM_binding) # Convert in genomic ranges
    
    # Create bed
    HLA_bed <- data.frame(seqnames=gsub("chr","",as.character(seqnames(GPPM_binding))),
                          starts=format(start(GPPM_binding)-1,scientific = F),
                          ends=format(end(GPPM_binding),scientific = F),
                          names=c(rep(".", length(GPPM_binding))),
                          scores=c(rep(".", length(GPPM_binding))),
                          strands=strand(GPPM_binding),stringsAsFactors = F)
    # Format to avoid issues with scientific annotation in igv!!!
    
    # Save bed
    write.table(HLA_bed, file=paste0("results/bed/",allele_name,".bed"), quote=F, sep="\t", row.names=F, col.names=F)
    
    # Keep GRanges
    assign(paste("annotation_gr",allele_name,sep="_"),GPPM_binding)
  }
  # Save Granges
  save(list=ls()[grepl("annotation",ls())],file = "data/annotation_gr.RData")
}


