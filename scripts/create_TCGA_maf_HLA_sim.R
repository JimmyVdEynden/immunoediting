###################################
# create_TCGA_maf_HLA_obs.R
###################################

# Load data
###########

# TCGA maf observation data
TCGA_maf_sim<- readRDS(file="data/TCGA_maf.rds")

# TCGA HLA types
load("data/TCGA_HLA_types.RData")

# Remove all data except HLA types, subs types, sample & cancer
###############################################################
data_to_keep<- c(
  "Tumor_Sample_Barcode",
  "Cancer",
  "subst_type",
  "subst_type3"
)

TCGA_maf_sim[,!colnames(TCGA_maf_sim)%in%data_to_keep]<- NA

# Add random HLA types
######################

# Shuffle HLA alleles in matrix with HLA types 
HLA_alleles<- c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")
TCGA_HLA_types_netMHC_shuffled<- TCGA_HLA_types_netMHC
for (HLA_allele in HLA_alleles){
  TCGA_HLA_types_netMHC_shuffled[,HLA_allele]<- sample(TCGA_HLA_types_netMHC_shuffled[,HLA_allele])
}

# Add to TCGA_maf_sim
TCGA_barcodes<- substr(TCGA_maf_sim[,"Tumor_Sample_Barcode"],1,12)
TCGA_maf_sim<- TCGA_maf_sim[TCGA_barcodes%in%rownames(TCGA_HLA_types_netMHC_shuffled),] # Only keep when HLA type determined  
TCGA_maf_HLA_types<- TCGA_HLA_types_netMHC_shuffled[substr(TCGA_maf_sim[,"Tumor_Sample_Barcode"],1,12),]
TCGA_maf_sim[,HLA_alleles]<- TCGA_maf_HLA_types

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")

# Simulate mutations based on triNT substitution type
#####################################################

# Load GPP
GPPM_all<- readRDS("temp/GPPM_all.rds")

subst_type3_all<- unique(as.character(TCGA_maf_sim[,"subst_type3"]))
pb <- txtProgressBar(min = 0, max = length(subst_type3_all), style = 3)
for(i in 1:length(subst_type3_all)){
  setTxtProgressBar(pb, i)
  # cat(i," ")
  subst_type3_temp<- subst_type3_all[i]
  idx_subst_type<- which(as.character(TCGA_maf_sim[,"subst_type3"])==subst_type3_temp)
  GPPM_subst_type3<- GPPM_all[GPPM_all$subst_type3==subst_type3_temp]
  GPPM_sim<- GPPM_subst_type3[sample(1:length(GPPM_subst_type3),size=length(idx_subst_type),replace = T)] 
  TCGA_maf_sim[idx_subst_type,"Hugo_Symbol"]<- GPPM_sim$gene
  TCGA_maf_sim[idx_subst_type,"Chromosome"]<- as.character(seqnames(GPPM_sim))
  TCGA_maf_sim[idx_subst_type,"Start_Position"]<- pos(GPPM_sim)
  TCGA_maf_sim[idx_subst_type,"End_Position"]<- pos(GPPM_sim)
  TCGA_maf_sim[idx_subst_type,"Reference_Allele"]<- GPPM_sim$ref
  TCGA_maf_sim[idx_subst_type,"Tumor_Seq_Allele2"]<- GPPM_sim$alt
  TCGA_maf_sim[idx_subst_type,"Protein_position"]<- GPPM_sim$aaPos
  TCGA_maf_sim[idx_subst_type,"ENSP"]<- GPPM_sim$ENSP
  TCGA_maf_sim[idx_subst_type,]$wt_pep<- GPPM_sim$nonaPep_wt     
}

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")

# Run annovar to add variant, ref_aa, alt_aa, PP2, CDS, ...
###########################################################

# Get annovar input
Annovar_input<- cbind(Chromosome=gsub("chr","",TCGA_maf_sim[,"Chromosome"]),Start_Position=TCGA_maf_sim[,"Start_Position"],End_Position=TCGA_maf_sim[,"End_Position"],Reference_Allele=TCGA_maf_sim[,"Reference_Allele"],Tumor_Seq_Allele2=TCGA_maf_sim[,"Tumor_Seq_Allele2"])
write.table(Annovar_input,"temp/TCGA_maf_sim_Annovar_input.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

# Run annovar
system(paste("perl ~/tools/annovar_old_version/table_annovar.pl temp/TCGA_maf_sim_Annovar_input.txt ../../../tools/annovar_old_version/humandb/ -buildver hg19 -out temp/TCGA_maf_sim_myanno -remove -protocol refGene,ljb26_all -operation g,f -nastring .",sep=""))

# Process annovar output
sim_avo<- read.table("temp/TCGA_maf_sim_myanno.hg19_multianno.txt",header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")
avo_temp<- sim_avo[,c("Gene.refGene","ExonicFunc.refGene","Polyphen2_HVAR_score","AAChange.refGene")]
colnames(avo_temp)<-c("Hugo_Symbol","Variant_Classification","PP2","Amino_acids")
avo_temp[,"Variant_Classification"]<- sapply(avo_temp[,"Variant_Classification"], function(x) strsplit(x, ";")[[1]][1]) # sapply much faster than loop!!!
avo_sil_idx<-which(avo_temp[,"Variant_Classification"]=="synonymous SNV")
avo_temp[avo_sil_idx,"PP2"]<- 0
avo_trunc_idx<-which(avo_temp[,"Variant_Classification"]=="stopgain")
avo_temp[avo_trunc_idx,"PP2"]<- 1
avo_temp[,"Amino_acids"]<- gsub("[[:digit:]]+","/",gsub(".*p\\.","",avo_temp[,"Amino_acids"]))
avo_temp[,"Amino_acids"]<- gsub("X","\\*",avo_temp[,"Amino_acids"]) # same annotation as maf data

TCGA_maf_sim[,"Variant_Classification"]<- avo_temp[,"Variant_Classification"]
TCGA_maf_sim[,"ref_aa"]<- substr(avo_temp[,"Amino_acids"],1,1)
TCGA_maf_sim[,"alt_aa"]<- substr(avo_temp[,"Amino_acids"],3,3)
TCGA_maf_sim[,"PP2"]<- as.numeric(avo_temp[,"PP2"])

# Remove lines where amino acid doesn't match
idx_noMatch<- which(as.character(sapply(TCGA_maf_sim[,"wt_pep"],function(x) substr(x[[1]],1,1)))!=TCGA_maf_sim[,"ref_aa"])
length(idx_noMatch) # 13575 --> Probably because other isoform was used
TCGA_maf_sim<- TCGA_maf_sim[-idx_noMatch,]
idx_noMatch<- which(as.character(sapply(TCGA_maf_sim[,"wt_pep"],function(x) substr(x[[9]],9,9)))!=TCGA_maf_sim[,"ref_aa"]) # In case first are NA!
length(idx_noMatch) # 287
TCGA_maf_sim<- TCGA_maf_sim[-idx_noMatch,]

# Remove lines where amino acid not determined by annovar
idx_noAA<- which(is.na(TCGA_maf_sim[,"ref_aa"]))
length(idx_noAA) # 8777
TCGA_maf_sim<- TCGA_maf_sim[-idx_noAA,]

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")

# # Test based on dN/dS!!!
# # If correct: n/s should be lowest in SKCM --> Ok!
# sort(
#   table(TCGA_maf_sim[,"Cancer"],TCGA_maf_sim[,"Variant_Classification"]!="synonymous SNV")[,"TRUE"]/
#   table(TCGA_maf_sim[,"Cancer"],TCGA_maf_sim[,"Variant_Classification"]=="synonymous SNV")[,"TRUE"]
#   )

# Add mRNA
##########

load(file="raw/TCGA_mRNA.RData")

row_idx<- match(TCGA_maf_sim[,"Hugo_Symbol"],rownames(TCGA_mRNA))
col_idx<- match(TCGA_maf_sim[,"Tumor_Sample_Barcode"],colnames(TCGA_mRNA))
mRNA<- TCGA_mRNA[cbind(row_idx,col_idx)]
TCGA_maf_sim[,"mRNA"]<- mRNA

# Add CNV class
###############

load("../../core_data/data/TCGA_CNV_gistic.RData")

row_idx<- match(TCGA_maf_sim[,"Hugo_Symbol"],rownames(CNV_gistic))
col_idx<- match(TCGA_maf_sim[,"Tumor_Sample_Barcode"],colnames(CNV_gistic))
CNV_class<- CNV_gistic[cbind(row_idx,col_idx)]
TCGA_maf_sim[,"CNV_class"]<- CNV_class

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")

# Add mutated nonapeptides
###########################

# Add mutated pep to TCGA_maf!!!
nonapep_mut<- TCGA_maf_sim[,"wt_pep"]
pb <- txtProgressBar(min = 0, max = nrow(TCGA_maf_sim), style = 3)
for(i in 1:nrow(TCGA_maf_sim)){
  setTxtProgressBar(pb, i)
  aa_ref_temp<- as.character(TCGA_maf_sim[i,"ref_aa"])
  aa_alt_temp<- as.character(TCGA_maf_sim[i,"alt_aa"])
  nonapep_wt_temp<- unlist(TCGA_maf_sim[i,"wt_pep"])
  for(j in 1:9){ 
    pep_temp<- nonapep_wt_temp[j]
    if(is.na(pep_temp)) next
    if(substr(pep_temp,j,j)!=aa_ref_temp) stop("Check aa")
    else substr(nonapep_mut[[i]][j],j,j)<- aa_alt_temp
  }
}
TCGA_maf_sim$mut_pep<- nonapep_mut

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")

# Add HLA affinities: exactly the same as for obs!!!
###################################################

# Add variable
TCGA_maf_sim<- cbind(TCGA_maf_sim,
                      wt_HLA_A1_aff=NA,
                      wt_HLA_A2_aff=NA,
                      wt_HLA_B1_aff=NA,
                      wt_HLA_B2_aff=NA,
                      wt_HLA_C1_aff=NA,
                      wt_HLA_C2_aff=NA,
                      mut_HLA_A1_aff=NA,
                      mut_HLA_A2_aff=NA,
                      mut_HLA_B1_aff=NA,
                      mut_HLA_B2_aff=NA,
                      mut_HLA_C1_aff=NA,
                      mut_HLA_C2_aff=NA
)
HLA_alleles<- c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")

# 1. WT PEPTIDES

# Get all pep wt & save
nonapep_wt_all<- na.omit(unique(unlist(TCGA_maf_sim$wt_pep)))
length(nonapep_wt_all) #7.623.881
write.table(nonapep_wt_all,file = "temp/wt_sim.pep",quote=F,row.names = F,col.names = F)

# Run NetMHC Pan (10x)
system("scripts/other/runNetMHCPan_wt.sh") 

# Fuse results
for(i in 1:10){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/wt_out_sim_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/wt_out_sim_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/wt_out_sim_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/wt_out_sim.rds")

# Get affinities
aff_wt_9_ls<- NULL
for(i in 1:6){ # 6 HLA alleles
  cat(i," ")
  col_idx<- match(TCGA_maf_sim[,HLA_alleles[i]],colnames(gene_mhc_nM)) # Identify HLA allele
  aff_wt_9<- matrix(NA,nrow(TCGA_maf_sim),9)
  for(j in 1:9){ # 9 peptides
    pep_temp<- sapply(TCGA_maf_sim[,"wt_pep"],function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_wt_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  TCGA_maf_sim[,grepl("wt_HLA",colnames(TCGA_maf_sim))][,i]<- rowMins(aff_wt_9,na.rm=T)
  # Save all affinities as well
  aff_wt_9_ls<- c(aff_wt_9_ls,list(aff_wt_9))
}
names(aff_wt_9_ls)<- HLA_alleles
saveRDS(aff_wt_9_ls,file="temp/TCGA_sim_aff_wt_9.rds")

saveRDS(TCGA_maf_sim,file="temp/TCGA_maf_sim.rds")


# 2 MUTATED PEPTIDES

# Get all pep mut & save
nonapep_mut_all<- na.omit(unique(unlist(TCGA_maf_sim$mut_pep)))
nonapep_mut_all<- nonapep_mut_all[!grepl("\\*",nonapep_mut_all)] # Stopgains --> No peptides
nonapep_mut_all_excl<- setdiff(nonapep_mut_all,nonapep_wt_all) # Not needed to calculate if already done for wts
length(nonapep_mut_all_excl)
write.table(nonapep_mut_all_excl,file = "temp/mut_sim.pep",quote=F,row.names = F,col.names = F)

# Run NetMHC Pan (10x)
system("scripts/other/runNetMHCPan_mut.sh")

# Fuse results
for(i in 1:10){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/mut_out_sim_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/mut_out_sim_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/mut_out_sim_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/mut_out_sim.rds")

# Add affinities from wt 
nonapep_mut_all<- na.omit(unique(unlist(TCGA_maf_sim$mut_pep)))
gene_mhc_nM_wt<- readRDS(file = "temp/wt_out_sim.rds")
gene_mhc_nM_wt<- gene_mhc_nM_wt[intersect(rownames(gene_mhc_nM_wt),nonapep_mut_all),]
gene_mhc_nM<- rbind(gene_mhc_nM,gene_mhc_nM_wt)

# Get affinities
aff_mut_9_ls<- NULL
for(i in 1:6){ # 6 HLA alleles
  cat(i," ")
  col_idx<- match(TCGA_maf_sim[,HLA_alleles[i]],colnames(gene_mhc_nM)) # Identify HLA allele
  aff_mut_9<- matrix(NA,nrow(TCGA_maf_sim),9)
  for(j in 1:9){ # 9 peptides
    pep_temp<- sapply(TCGA_maf_sim[,"mut_pep"],function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_mut_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  TCGA_maf_sim[,grepl("mut_HLA",colnames(TCGA_maf_sim))][,i]<- rowMins(aff_mut_9,na.rm=T)
  # Save all affinities as well
  aff_mut_9_ls<- c(aff_mut_9_ls,list(aff_mut_9))
}
names(aff_mut_9_ls)<- HLA_alleles
saveRDS(aff_mut_9_ls,file="temp/TCGA_sim_aff_mut_9.rds")

# Add harmonic means
source("scripts/functions/harmonic_mean.R")
TCGA_maf_sim<- cbind(TCGA_maf_sim,wt_HLA_mean_aff=apply(TCGA_maf_sim[,grepl("wt_HLA",colnames(TCGA_maf_sim))],1,harmonic_mean))
TCGA_maf_sim<- cbind(TCGA_maf_sim,mut_HLA_mean_aff=apply(TCGA_maf_sim[,grepl("mut_HLA",colnames(TCGA_maf_sim))],1,harmonic_mean))

# Save final TCGA_maf
saveRDS(TCGA_maf_sim,file="data/TCGA_maf_sim.rds")


# Add alternative ways to calculate aggregate HLA binding
##########################################################

# # 1) All 6 alleles seperate
# n_mut<- nrow(TCGA_maf_sim)
# TCGA_maf_sim_6<- TCGA_maf_sim[rep(1:n_mut,each=6),]
# wt_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_sim_6[,grepl("wt_HLA_.._aff",colnames(TCGA_maf_sim_6))])*diag(6)[rep(1:6,n_mut),])
# mut_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_sim_6[,grepl("mut_HLA_.._aff",colnames(TCGA_maf_sim_6))])*diag(6)[rep(1:6,n_mut),])
# TCGA_maf_sim_6[,"wt_HLA_mean_aff"]<- wt_HLA_aff
# TCGA_maf_sim_6[,"mut_HLA_mean_aff"]<- mut_HLA_aff
# saveRDS(TCGA_maf_sim_6,file="data/TCGA_maf_sim_6.rds")
