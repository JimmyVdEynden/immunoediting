##############################
# create_subst_type_matrix.R
##############################

# Create matrices with total number of subst. types in N+, N-, S+ & S-, ... 
# 1) main subst. types 
# 2) trinucleotide subst. types
# 3) SSB7
# 4) pentaNT subst types

# Complete GPPM
################
{
  # Load
  GPPM<- readRDS("data/GPPM_inclHLAAlleles.rds")
  CpG_conversion_vector<- readRDS("data/SSB7_conversion_vector.rds")
  
  # Functions
  source("scripts/functions/get_subst_types5.R")
  
  # Process
  GPPM$isHA<- GPPM$HLA_aff_mean<500
  
  # Different subst models
  subst_type_t<- table(GPPM$subst_type,GPPM$variant,GPPM$isHA) # main subst types
  subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA) # triNT subst types
  SSB7_t<- table(CpG_conversion_vector[GPPM$subst_type3],GPPM$variant,GPPM$isHA) # SSB7
  
  # Kd cut-off of 50
  GPPM$isHA<- GPPM$HLA_aff_mean<50
  Kd50_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA) 

  # Best binder instead of harmonic mean
  GPPM$isHA<- rowMins(data.matrix(mcols(GPPM)[grep("HLA_aff_.[[:digit:]]",names(mcols(GPPM)))]))<500          
  bestBinder_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA) # triNT subst types
  
  # PentaNT
  # Get variables and rm rest to reduce ram
  chr<- as.character(seqnames(GPPM))
  pos<- pos(GPPM)
  alt_allele<- as.character(GPPM$alt)
  isHA<- GPPM$HLA_aff_mean<500
  variant<- GPPM$variant
  rm(GPPM)
  gc()
  GPPM_subst_type5<- get_subst_types5(chr = chr,pos = pos,alt_allele = alt_allele)
  saveRDS(GPPM_subst_type5,file="temp/GPPM_subst_type5.rds")
  # Create subst matrix
  subst_type5_t<- table(GPPM_subst_type5,variant,isHA) # pentaNT subst types

  # Put in 1 object
  subst_type_matrix_ls<- list()
  subst_type_matrix_ls$main<- subst_type_t
  subst_type_matrix_ls$triNT<- subst_type3_t
  subst_type_matrix_ls$SSB7<- SSB7_t
  subst_type_matrix_ls$triNT_Kd50<- Kd50_subst_type3_t
  subst_type_matrix_ls$triNT_bestBinder<- bestBinder_subst_type3_t
  subst_type_matrix_ls$pentaNT<- subst_type5_t
  
  # Save
  saveRDS(subst_type_matrix_ls,file="data/subst_type_matrix_ls.rds")
  
}

# Complete GPPM subset (with mutated peptides)
###################################################
{
  # Load
  GPPM<- readRDS("data/GPPM_subset_mut.rds")
  load("temp/CpG_conversion.RData")
  
  # Process: watch out for infinite for stopgains!
  GPPM$isHA_wt<- GPPM$HLA_aff_wt_mean<500
  GPPM$isHA_mut<- GPPM$HLA_aff_mut_mean<500
  
  # main subst types
  wt_subst_type_t<- table(GPPM$subst_type,GPPM$variant,GPPM$isHA_wt)
  mut_subst_type_t<- table(GPPM$subst_type,GPPM$variant,GPPM$isHA_mut)
  
  # triNT subst types
  wt_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_wt)
  mut_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_mut)
  
  # SSB7
  wt_SSB7_t<- table(CpG_conversion_vector[GPPM$subst_type3],GPPM$variant,GPPM$isHA_wt)
  mut_SSB7_t<- table(CpG_conversion_vector[GPPM$subst_type3],GPPM$variant,GPPM$isHA_mut)

  # Kd cut-off of 50
  GPPM$isHA_wt<- GPPM$HLA_aff_wt_mean<50
  GPPM$isHA_mut<- GPPM$HLA_aff_mut_mean<50
  wt_Kd50_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_wt) 
  mut_Kd50_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_mut) 
  
  # Best binder instead of harmonic mean
  GPPM$isHA_wt<- rowMins(data.matrix(GPPM$HLA_aff_wt))<500          
  GPPM$isHA_mut<- rowMins(data.matrix(GPPM$HLA_aff_mut))<500          
  wt_bestBinder_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_wt) # triNT subst types
  mut_bestBinder_subst_type3_t<- table(GPPM$subst_type3,GPPM$variant,GPPM$isHA_mut) # triNT subst types
  
  # PentaNT
  # Get variables and rm rest to reduce ram
  chr<- as.character(seqnames(GPPM))
  pos<- pos(GPPM)
  alt_allele<- as.character(GPPM$alt)
  isHA_wt<- GPPM$HLA_aff_wt_mean<500
  isHA_mut<- GPPM$HLA_aff_mut_mean<500
  variant<- GPPM$variant
  rm(GPPM)
  gc()
  subst_type5<- get_subst_types5(chr = chr,pos = pos,alt_allele = alt_allele)
  wt_subst_type5_t<- table(subst_type5,variant,isHA_wt)
  mut_subst_type5_t<- table(subst_type5,variant,isHA_mut) 
  
  # Put in 1 object
  subst_type_matrix_ls<- list()
  subst_type_matrix_ls$main_wt<- wt_subst_type_t
  subst_type_matrix_ls$triNT_wt<- wt_subst_type3_t
  subst_type_matrix_ls$SSB7_wt<- wt_SSB7_t
  subst_type_matrix_ls$main_mut<- mut_subst_type_t
  subst_type_matrix_ls$triNT_mut<- mut_subst_type3_t
  subst_type_matrix_ls$SSB7_mut<- mut_SSB7_t
  subst_type_matrix_ls$Kd50_wt<- wt_Kd50_subst_type3_t
  subst_type_matrix_ls$Kd50_mut<- mut_Kd50_subst_type3_t
  subst_type_matrix_ls$best_binder_wt<- wt_bestBinder_subst_type3_t
  subst_type_matrix_ls$best_binder_mut<- mut_bestBinder_subst_type3_t
  subst_type_matrix_ls$pentaNT_wt<- wt_subst_type5_t
  subst_type_matrix_ls$pentaNT_mut<- mut_subst_type5_t
  
  # Save
  saveRDS(subst_type_matrix_ls,file="data/GPPM_subset_subst_type_matrix_ls.rds")
}
