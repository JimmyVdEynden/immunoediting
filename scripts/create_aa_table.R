#############################################
# manuscript_create_aa_table.R
#############################################

# Create triNT x aa x variant & related tables from GPPM data  

# Load
#######
{
  GPPM<- readRDS(file="data/GPPM_inclHLAAlleles.rds") # Takes some time to load
}

# Create table
##############
{
  # triNT x aa x variant
  triNT_aa_var_t<- table(GPPM$subst_type3,GPPM$ref_aa,GPPM$variant)
  # gene x aa
  gene_aa_t<- table(GPPM$gene,GPPM$ref_aa)
  # gene x aa x HLA aff
  gene_aa_HA_t<- table(GPPM$gene,GPPM$ref_aa,GPPM$HLA_aff_mean<500)
  # gene x variant x HLA aff
  gene_var_HA_t<- table(GPPM$gene,GPPM$variant,GPPM$HLA_aff_mean<500)
  
}

# Define aa classes
###################
{
  aa_hp<- c("G","A","P","V","L","I","M","W","F")
  aa_polar<- c("S","T","Y","N","Q","C")
  aa_charged<- c("K","R","H","D","E")
  aa_all<- sort(c(aa_hp,aa_polar,aa_charged))
}

# Save
##############
{
  saveRDS(triNT_aa_var_t,"data/GPPM_triNT_aa_var_t.rds")  
  saveRDS(gene_aa_t,"data/GPPM_gene_aa_t.rds")  
  saveRDS(gene_aa_HA_t,"data/GPPM_gene_aa_HA_t.rds")  
  saveRDS(gene_var_HA_t,"data/GPPM_gene_var_HA_t.rds")  
  save(aa_hp,aa_hp,aa_polar,aa_charged,aa_all,file = "data/aa_classes.RData")
}




