###################################
# add_mainHLAType_to_maf_obs.R
###################################

# Add main HLA type (i.e. most freq alleles) to maf file. 

# Load maf file
###############
TCGA_maf<- readRDS("data/TCGA_maf_sim.rds")

# Load harmonic mean funcion
##############################
source("scripts/functions/harmonic_mean.R")

# Change HLA type to main type (most freq alleles: "HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B08:01", "HLA-C07:01", "HLA-C07:02")
######################################################################################################################################
HLA_alleles<- c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")
for(i in 1:6){TCGA_maf[,HLA_alleles[i]]<- c("HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B08:01", "HLA-C07:01", "HLA-C07:02")[i]}

# Erase affinities 
####################
TCGA_maf[,grepl("wt_HLA",colnames(TCGA_maf))]<- NA
TCGA_maf[,grepl("mut_HLA",colnames(TCGA_maf))]<- NA

# Calculate affinities for wt
##############################
gene_mhc_nM<- readRDS(file = "temp/wt_out_sim.rds")
gene_mhc_nM<- gene_mhc_nM[,c("HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B08:01", "HLA-C07:01", "HLA-C07:02")]
aff_wt_9_ls<- NULL
for(i in 1:6){ # 6 HLA alleles
  cat(i," ")
  col_idx<- match(TCGA_maf[,HLA_alleles[i]],colnames(gene_mhc_nM)) # Identify HLA allele
  aff_wt_9<- matrix(NA,nrow(TCGA_maf),9)
  for(j in 1:9){ # 9 peptides
    pep_temp<- sapply(TCGA_maf[,"wt_pep"],function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_wt_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  TCGA_maf[,grepl("wt_HLA",colnames(TCGA_maf))][,i]<- rowMins(aff_wt_9,na.rm=T)
  aff_wt_9_ls<- c(aff_wt_9_ls,list(aff_wt_9))
}
names(aff_wt_9_ls)<- HLA_alleles
saveRDS(aff_wt_9_ls,file="temp/TCGA_sim_aff_wt_9_mainHLAType.rds")
saveRDS(TCGA_maf,file="temp/TCGA_maf_sim_mainHLAType.rds")

# Calculate affinities for mut
##############################
gene_mhc_nM<- readRDS(file = "temp/mut_out_sim.rds")
gene_mhc_nM<- gene_mhc_nM[,c("HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B08:01", "HLA-C07:01", "HLA-C07:02")]
nonapep_mut_all<- na.omit(unique(unlist(TCGA_maf$mut_pep)))
gene_mhc_nM_wt<- readRDS(file = "temp/wt_out_sim.rds")
gene_mhc_nM_wt<- gene_mhc_nM_wt[,c("HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B08:01", "HLA-C07:01", "HLA-C07:02")]
gene_mhc_nM_wt<- gene_mhc_nM_wt[intersect(rownames(gene_mhc_nM_wt),nonapep_mut_all),]
gene_mhc_nM<- rbind(gene_mhc_nM,gene_mhc_nM_wt)

aff_mut_9_ls<- NULL
for(i in 1:6){ # 6 HLA alleles
  cat(i," ")
  col_idx<- match(TCGA_maf[,HLA_alleles[i]],colnames(gene_mhc_nM)) # Identify HLA allele
  aff_mut_9<- matrix(NA,nrow(TCGA_maf),9)
  for(j in 1:9){ # 9 peptides
    pep_temp<- sapply(TCGA_maf[,"mut_pep"],function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_mut_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  TCGA_maf[,grepl("mut_HLA",colnames(TCGA_maf))][,i]<- rowMins(aff_mut_9,na.rm=T)
  aff_mut_9_ls<- c(aff_mut_9_ls,list(aff_mut_9))
}
names(aff_mut_9_ls)<- HLA_alleles
saveRDS(aff_mut_9_ls,file="temp/TCGA_sim_aff_mut_9_mainHLAType.rds")

# Add harmonic means to obs and mut
###################################
TCGA_maf[,"wt_HLA_mean_aff"]<- apply(TCGA_maf[,grepl("wt_HLA_.._.*",colnames(TCGA_maf))],1,harmonic_mean)
TCGA_maf[,"mut_HLA_mean_aff"]<- apply(TCGA_maf[,grepl("mut_HLA_.._.*",colnames(TCGA_maf))],1,harmonic_mean)

# Save final TCGA_maf
#####################
saveRDS(TCGA_maf,file="data/TCGA_maf_sim_mainHLAType.rds")

# Add alternative ways to calculate aggregate HLA binding
##########################################################

# TCGA_maf<- readRDS(file="data/TCGA_maf_sim_mainHLAType.rds")

# 1) All 6 alleles seperate
n_mut<- nrow(TCGA_maf)
TCGA_maf_6<- TCGA_maf[rep(1:n_mut,each=6),]
wt_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_6[,grepl("wt_HLA_.._aff",colnames(TCGA_maf_6))])*diag(6)[rep(1:6,n_mut),])
mut_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_6[,grepl("mut_HLA_.._aff",colnames(TCGA_maf_6))])*diag(6)[rep(1:6,n_mut),])
TCGA_maf_6[,"wt_HLA_mean_aff"]<- wt_HLA_aff
TCGA_maf_6[,"mut_HLA_mean_aff"]<- mut_HLA_aff
saveRDS(TCGA_maf_6,file="data/TCGA_maf_sim_mainHLAType_6.rds")



