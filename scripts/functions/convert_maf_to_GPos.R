##########################################
# convert_maf_to_GPos.R
##########################################

# Function to convert mutation data from maf to GPos

Convert_maf_to_GPos.R<- function(maf,metadata_cols,metadata_names=NULL){
  
  # Check metadata info
  if(is.null(metadata_names)) metadata_names<- metadata_cols
  if(length(metadata_cols)!=length(metadata_names)) stop("Check length of metadata cols and names")
  
  # Turn data in GPos
  maf_GP <- GRanges(seqnames = gsub("chr","",maf[,"Chromosome"]), 
                         ranges = IRanges(as.numeric(maf[,"Start_Position"]), end = as.numeric(maf[,"End_Position"])))
  maf_GP <- GPos(maf_GP)
  
  # Add metadata
  for(i in 1:length(metadata_cols)){
    maf_GP$m<- maf[,metadata_cols[i]]
    colnames(mcols(maf_GP))[i]<- metadata_names[i]
  }
  
  # Sort
  maf_GP<- sortSeqlevels(maf_GP) # Sort chromosomes
  maf_GP<- sort(maf_GP) # Sort GRanges
  
  # Return
  return(maf_GP)
}
  
  


