######################################
# create_GPPM_subset_mut.R
######################################

# Aim: calculate for each subst type (10000 each) the mutated peptide HLA affinities for prototypical HLA
# Note: Solution to avoid having to do this on complete GPPM files, which is computationally harder

GPPM<- readRDS("data/GPPM_inclHLAAlleles.rds")

# Take 10000 random muts from each subst type
subst_types3<- unique(GPPM$subst_type3)
for(i in 1:length(subst_types3)){
  cat(i," ")
  subst_type3<- subst_types3[i]
  idx_temp<- which(GPPM$subst_type3==subst_type3)
  idx_temp<- sample(idx_temp,10000,replace = F)
  if(i==1) GPPM_mut<- GPPM[idx_temp]
  else GPPM_mut<- c(GPPM_mut,GPPM[idx_temp])
}
rm(GPPM)
gc()

# Get mutated peptides
aa_ref_temp<- GPPM_mut$ref_aa
aa_alt_temp<- GPPM_mut$alt_aa
nonapep_mut<- NULL
for(i in 1:9){
  cat(i," ")
  nonapep_mut_temp<- sapply(GPPM_mut$nonaPep_wt,function(x) x[[i]])
  substr(nonapep_mut_temp,i,i)<- aa_alt_temp
  nonapep_mut<- cbind(nonapep_mut,nonapep_mut_temp)
}
nonapep_mut_unique<- setdiff(unique(nonapep_mut),NA)
length(nonapep_mut_unique) #8,237,909
write.table(nonapep_mut_unique,file = "temp/GPPM_subset_mut.pep",quote=F,row.names = F,col.names = F)

# Add mutated peptides to GPPM_mut (list, similar as wt)
GPPM_mut$nonaPep_mut<- as.list(data.frame(t(nonapep_mut),stringsAsFactors=F))

# Save GPPM_mut
saveRDS(GPPM_mut,file = "data/GPPM_subset_mut.rds")

# Run NetMHC Pan (1x): mostly around 12GB, total time +/- 6.5U --> On this subset could consider all HLA types and create sample-specific correction
# system("~/tools/netMHCpan-3.0/netMHCpan -a `head -1 temp/HLA_types_freq_oneLine.txt | tail -1` -p temp/GPPM_subset_mut.pep -l 9 -xls -xlsfile temp/GPPM_subset_mut_out.xls > /dev/null &") # Ignore warnings,works!

# Get affinities
gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/GPPM_subset_mut_out.xls"),n=1),"\t")),"")
gene_mhc_nM<- read.table(paste0("temp/GPPM_subset_mut_out.xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
colnames(gene_mhc_nM)<- gene_mhc_colnames
pep<- read.table(paste0("temp/GPPM_subset_mut_out.xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/GPPM_subset_mut_out.rds")

# Add affinities to GPPM
aff_mut<- matrix(NA,length(GPPM_mut),6,dimnames=list(NULL,c("A1","A2","B1","B2","C1","C2")))
for(i in 1:6){ # 6 HLA alleles
  cat("\n",i," ")
  col_idx<- i # Identify HLA allele
  aff_mut_9<- matrix(NA,length(GPPM_mut),9)
  for(j in 1:9){ # 9 peptides
    cat(j, " ")
    pep_temp<- sapply(GPPM_mut$nonaPep_mut,function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_mut_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  aff_mut[,i]<- rowMins(aff_mut_9,na.rm=T)
}
colnames(aff_mut)<- gene_mhc_colnames

# Add harmonic mean
source("scripts/functions/harmonic_mean.R")
mut_HLA_mean_aff<- apply(aff_mut,1,harmonic_mean)

# Add to GPPM
GPPM_mut$mean_HLA_aff_mut<- mut_HLA_mean_aff

# Add affinities to GPPM
GPPM_mut$HLA_aff_mut<- aff_mut

# Some general processing
# Put ind affinities of wt in one matrix
aff_wt<- data.matrix(mcols(GPPM_mut)[,9:14])
colnames(aff_wt)<- gene_mhc_colnames
GPPM_mut$HLA_aff_wt<- aff_wt

# Reoder and rename some metadata
mcols(GPPM_mut)<- mcols(GPPM_mut)[c('gene','ENSP','variant','ref','alt','subst_type','subst_type3','aaPos','ref_aa','alt_aa','PP2','nonaPep_wt','HLA_aff_wt','HLA_aff_mean','nonaPep_mut','HLA_aff_mut','mean_HLA_aff_mut')]
names(mcols(GPPM_mut))<- gsub("HLA_aff_mean","HLA_aff_wt_mean",names(mcols(GPPM_mut)))
names(mcols(GPPM_mut))<- gsub("mean_HLA_aff_mut","HLA_aff_mut_mean",names(mcols(GPPM_mut)))

# Save GPPM_mut
saveRDS(GPPM_mut,file = "data/GPPM_subset_mut.rds")

# Add all TCGA HLA alleles
###########################

# Run NetMHC Pan (10x)
# system("runNetMHCPan_GPPM_subset_mut.sh")

# Fuse results
for(i in 1:10){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/GPPM_subset_mut_out_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/GPPM_subset_mut_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/GPPM_subset_mut_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/GPPM_subset_mut_out_allHLA.rds")

# Add affinities to GPPM
GPPM_mut<- readRDS(file = "data/GPPM_subset_mut.rds")
aff_mut<- matrix(NA,length(GPPM_mut),ncol(gene_mhc_nM),dimnames=list(NULL,colnames(gene_mhc_nM)))
aff_mut_9_ls<- list()
for(j in 1:9){ # 9 peptides
  cat(j, " ")
  pep_temp<- sapply(GPPM_mut$nonaPep_mut,function(x) x[[j]])
  row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
  aff_mut_9_ls<- c(aff_mut_9_ls,list(gene_mhc_nM[row_idx,]))
}
# Best binders per HLA type
for(j in 1:ncol(aff_mut)){
  cat(j, " ")
  aff_mut_temp<- cbind(aff_mut_9_ls[[1]][,j],aff_mut_9_ls[[2]][,j],aff_mut_9_ls[[3]][,j],aff_mut_9_ls[[4]][,j],aff_mut_9_ls[[5]][,j],aff_mut_9_ls[[6]][,j],aff_mut_9_ls[[7]][,j],aff_mut_9_ls[[8]][,j],aff_mut_9_ls[[9]][,j])
  aff_mut[,j]<- rowMins(aff_mut_temp,na.rm=T)
}
GPPM_mut$HLA_aff_mut_allHLA<- aff_mut
# saveRDS(aff_mut,file="temp/GPPM_subset_aff_mut.rds")
saveRDS(GPPM_mut, file = "temp/GPPM_subset_mut.rds")

# Calculate harmonic mean for each TCGA sample!
#################################################
GPPM_mut<- readRDS(file = "data/GPPM_subset_mut.rds")
aff_mut<- readRDS(file = "temp/GPPM_subset_aff_mut.rds")
load("data/TCGA_HLA_types.RData")
source("scripts/functions/harmonic_mean.R")

subst_type_matrix_ls<- list()
for(i in 1:nrow(TCGA_HLA_types_netMHC)){
# for(i in 1:10){
  cat(i," ")
  HLA_temp<- TCGA_HLA_types_netMHC[i,]
  col_idx<- match(HLA_temp,colnames(aff_mut)) # Identify HLA type
  mean_aff_temp<- apply(aff_mut[,col_idx],1,harmonic_mean)
  # Subst_type_matrix
  subst_type_matrix<- table(GPPM_mut$subst_type3,GPPM_mut$variant,mean_aff_temp<500)
  subst_type_matrix_ls[[i]]<- subst_type_matrix
  # if (i==100) saveRDS(subst_type_matrix_ls,file="temp/subst_type_matrix_ls.rds")
}
names(subst_type_matrix_ls)<- rownames(TCGA_HLA_types_netMHC)
saveRDS(subst_type_matrix_ls,file="data/GPPM_subset_subst_type_matrix_mut_per_sample.rds")

# Repeat for simulated mutation data
####################################

# TCGA barcode had scrambled HLA alleles and hence different genotypes
TCGA_maf_sim<- readRDS("data/TCGA_maf_sim.rds")
TCGA_HLA_types_sim<- TCGA_maf_sim[,c("Tumor_Sample_Barcode","HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")]
TCGA_HLA_types_sim<- TCGA_HLA_types_sim[!duplicated(TCGA_HLA_types_sim$Tumor_Sample_Barcode),]
rownames(TCGA_HLA_types_sim)<- TCGA_HLA_types_sim[,1]
TCGA_HLA_types_sim<- TCGA_HLA_types_sim[,-1]

# Calculate
subst_type_matrix_sim_ls<- list()
for(i in 1:nrow(TCGA_HLA_types_sim)){
  cat(i," ")
  HLA_temp<- TCGA_HLA_types_sim[i,]
  col_idx<- match(HLA_temp,colnames(aff_mut)) # Identify HLA type
  mean_aff_temp<- apply(aff_mut[,col_idx],1,harmonic_mean)
  # Subst_type_matrix
  subst_type_matrix<- table(GPPM_mut$subst_type3,GPPM_mut$variant,mean_aff_temp<500)
  subst_type_matrix_sim_ls[[i]]<- subst_type_matrix
  if (i==100) saveRDS(subst_type_matrix_sim_ls,file="temp/subst_type_matrix_sim_ls.rds")
}
names(subst_type_matrix_sim_ls)<- substr(rownames(TCGA_HLA_types_sim),1,12)
saveRDS(subst_type_matrix_sim_ls,file="data/GPPM_subset_subst_type_matrix_sim_mut_per_sample.rds")

# Calculate expected wt/mut for prototypical HLA type 
#####################################################

aff_mut<- readRDS(file = "temp/GPPM_subset_aff_mut.rds")

# observed

# Simulated
TCGA_maf_sim<- readRDS("data/TCGA_maf_sim.rds")
main_HLA_type<- c("HLA-A02:01","HLA-A01:01","HLA-B07:02","HLA-B08:01","HLA-C07:01","HLA-C07:02") 

subst_type_matrix_main_sim_ls<- list()
HLA_temp<- main_HLA_type
col_idx<- match(HLA_temp,colnames(aff_mut)) 
# Harmonic mean
mean_aff_temp<- apply(aff_mut[,col_idx],1,harmonic_mean)
subst_type_matrix_main_sim_ls$Kd500<- table(GPPM_mut$subst_type3,GPPM_mut$variant,mean_aff_temp<500)
subst_type_matrix_main_sim_ls$Kd50<- table(GPPM_mut$subst_type3,GPPM_mut$variant,mean_aff_temp<50)
# Best binder
best_aff_temp<- apply(aff_mut[,col_idx],1,rowMins)
subst_type_matrix_main_sim_ls$best_binder<- table(GPPM_mut$subst_type3,GPPM_mut$variant,min_aff_temp<500)


