create_GPPM_from_GPP<- function(GPP,TCGA_maf){
  
  # Remove alt_init_codons from GENETIC_CODE
  # Otherwise fro some reason CTG will ALWAYS be translated to methionine!!! e.g. pos 58858728 on A1BG
  # If alt_init_codons are the first ones, they will be used instead --> But apparently always done?
  # CHange this manually
  genetic_code_no_init<- GENETIC_CODE
  attributes(genetic_code_no_init)$alt_init_codons<- "ATG" # Easiest solution is to just change to normal initiation codon
  
  # Replicate all positions 4 times (4 different nt)
  GPPM_temp<- rep(GPP,each=4)
  GPPM_temp$nt_mut<- DNAStringSet(rep(c("A","C","G","T"),length(GPP)))
  
  # Remove substitutions by same nucleotide
  GPPM_temp<- GPPM_temp[!GPPM_temp$nt==GPPM_temp$nt_mut]
  
  # Convert genome build to hg19 (otherwise 37)
  genome(GPPM_temp)<- "hg19"
  
  # Get mono- and trinucleotide subsitution types
  subst_type_ref<- GPPM_temp$nt
  subst_type_alt<- GPPM_temp$nt_mut
  subst_type_triNT<- GPPM_temp$triNt
  for(nt in c("A","G")){
    idx_subst_type_temp<- which(subst_type_ref==nt)
    subst_type_ref[idx_subst_type_temp]<- reverseComplement(subst_type_ref[idx_subst_type_temp])
    subst_type_alt[idx_subst_type_temp]<- reverseComplement(subst_type_alt[idx_subst_type_temp])
    subst_type_triNT[idx_subst_type_temp]<- reverseComplement(subst_type_triNT[idx_subst_type_temp])
  }
  subst_type_triNT_alt<- paste(substr(subst_type_triNT,1,1),subst_type_alt,substr(subst_type_triNT,3,3),sep="")
  subst_type<- paste(subst_type_ref,subst_type_alt,sep=">")
  subst_type3<- paste(subst_type_triNT,subst_type_triNT_alt,sep=">")
  GPPM_temp$subst_type<- subst_type
  GPPM_temp$subst_type3<- subst_type3
  
  # Select data for 1 ENSP protein
  # Take ENSP information that is used in TCGA data. Most common if multiple! First if nothing in TCGA for the specific gene.
  ENSP<- names(sort(table(TCGA_maf[TCGA_maf[,1]==gene,"ENSP"]),decreasing = T)[1]) # Take most common if multiple  
  if(is.null(ENSP)||!ENSP%in%colnames(GPPM_temp$nonaPep)) ENSP<- colnames(GPPM_temp$nonaPep)[1]
  
  # Build final GPPM object
  GPPM<- GPos(GPPM_temp)
  GPPM$gene<- GPPM_temp$gene
  GPPM$ENSP<- ENSP
  GPPM$ref<- GPPM_temp$nt
  GPPM$alt<- GPPM_temp$nt_mut
  GPPM$subst_type<- GPPM_temp$subst_type
  GPPM$subst_type3<- GPPM_temp$subst_type3
  GPPM$nonaPep_wt<- GPPM_temp$nonaPep[,ENSP]
  GPPM$aaPos<- GPPM_temp$aaPos[,ENSP]
  GPPM<- GPPM[GPPM$aaPos!=0] # Only locations for nonapep
}
