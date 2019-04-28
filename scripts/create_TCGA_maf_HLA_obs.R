###################################
# create_TCGA_maf_HLA_obs.R
###################################

# Load data
###########

# TCGA maf data
load("raw/TCGA_maf_mutect2_v7.RData")

# TCGA HLA types
load("data/TCGA_HLA_types.RData")

# Filter TCGA_maf
#################

# Filter 1: restrict analysis to SNPs 
TCGA_maf<- TCGA_maf[TCGA_maf[,"Variant_Type"]=="SNP",]

# Filter 2: restrict analysis to known variants, also exclude few stoplosses (hard to deal with) 
TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"]),]
TCGA_maf<- TCGA_maf[TCGA_maf[,"Variant_Classification"]!="unknown",]
TCGA_maf<- TCGA_maf[TCGA_maf[,"Variant_Classification"]!="stoploss",]

# Filter 3: Change aa information needs to be present
TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Amino_acids"])&TCGA_maf[,"Amino_acids"]!="",]

# Add aa information to maf
###############################
aa_change<- TCGA_maf[,"Amino_acids"]
ref_aa<- gsub("/.*","",aa_change)
alt_aa<- gsub(".*/","",aa_change)
TCGA_maf<- cbind(TCGA_maf,ref_aa,alt_aa)

# Restrict variables to the ones necessary for study
####################################################
selected_variables<- c(
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode",
  "ref_aa",
  "alt_aa",
  "t_ref_count",
  "t_alt_count",
  "CDS_position",
  "Protein_position",
  "Codons",
  "TRANSCRIPT_STRAND",
  "ENSP",
  "Cancer",                      
  "vaf",
  "Variant_Classification",
  "PP2",
  "mRNA",
  "CNV_class",
  "Purity"
)

TCGA_maf<- TCGA_maf[,selected_variables]

# Add information on subst types
################################
require(Biostrings)
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)

# Get triNT ref information
mutTriNT_Nl<- getSeq(Hsapiens,TCGA_maf[,"Chromosome"],as.numeric(TCGA_maf[,"Start_Position"])-1,as.numeric(TCGA_maf[,"End_Position"])+1)

# Get triNT var information
mutTriNT_T<-mutTriNT_Nl
mutNT_T<-TCGA_maf[,"Tumor_Seq_Allele2"]
mutNT_T[mutNT_T=="0"]<-"" #Zero not recognized by DNAString
mutNT_T<-DNAStringSet(mutNT_T)
mutNT_T[width(mutNT_T)!=1]<- ""
subseq(mutTriNT_T,2,2)<-mutNT_T

# Take pyrimidine as ref
idx_compl<-which(substr(mutTriNT_Nl,2,2)%in%c("A","G"))
mutTriNT_Nl[idx_compl]<-reverseComplement(mutTriNT_Nl[idx_compl])
mutTriNT_T[idx_compl]<-reverseComplement(mutTriNT_T[idx_compl])
mutTriNT_Nl[width(mutTriNT_Nl)!=3]<- DNAString("N")
mutTriNT_T[width(mutTriNT_T)!=3]<- DNAString("N")

# Calculate substitution type
subst_type<- paste0(substr(mutTriNT_Nl,2,2),">",substr(mutTriNT_T,2,2))

# Calculate triNT substitution type
subst_type3<- paste0(mutTriNT_Nl,">",mutTriNT_T)

# Add information to TCGA_maf
TCGA_maf<- cbind(TCGA_maf,subst_type,subst_type3)

# Remove samples where given ref allele doesn't match hg19 allele
hg19_ref<- getSeq(Hsapiens,TCGA_maf[,"Chromosome"],as.numeric(TCGA_maf[,"Start_Position"]),as.numeric(TCGA_maf[,"End_Position"]))
idx_temp<- which(as.character(hg19_ref)!=TCGA_maf[,"Reference_Allele"]) 
length(idx_temp) # 6547
TCGA_maf<- TCGA_maf[-idx_temp,]

# Add HLA types
###############
TCGA_barcodes<- substr(TCGA_maf[,"Tumor_Sample_Barcode"],1,12)
TCGA_maf<- TCGA_maf[TCGA_barcodes%in%rownames(TCGA_HLA_types_netMHC),] # Only keep when HLA type determined  
TCGA_maf_HLA_types<- TCGA_HLA_types_netMHC[substr(TCGA_maf[,"Tumor_Sample_Barcode"],1,12),]
TCGA_maf<- cbind(TCGA_maf,TCGA_maf_HLA_types)
# 
# saveRDS(TCGA_maf,file="temp/TCGA_maf.rds")

# Add nonapeptides
###################

# Load GPP
GPPM_all<- readRDS("temp/GPPM_all.rds")

# Get nonapep_wt
GPPM_pep<- GPPM_all$nonaPep_wt
names(GPPM_pep)<- paste(GPPM_all$gene,pos(GPPM_all),GPPM_all$subst_type,sep="_")
TCGA_maf_ids<- paste(TCGA_maf[,"Hugo_Symbol"],TCGA_maf[,"Start_Position"],TCGA_maf[,"subst_type"],sep="_")
GPPM_pep<- GPPM_pep[names(GPPM_pep)%in%TCGA_maf_ids]
GPPM_pep_idx<- match(TCGA_maf_ids,names(GPPM_pep)) # >100x faster than GPPM_pep[TCGA_maf_ids]
nonapep_wt<- GPPM_pep[GPPM_pep_idx]

# Add to TCGA_maf
TCGA_maf<- cbind(TCGA_maf,wt_pep=NA)
TCGA_maf$wt_pep<- nonapep_wt

# Remove lines where no information on nonapep_wt
idx_noNonapep<- which(is.na(names(nonapep_wt)))
length(idx_noNonapep) # 101519
TCGA_maf<- TCGA_maf[-idx_noNonapep,]

# Remove lines where amino acid doesn't match
idx_noMatch<- which(as.character(sapply(TCGA_maf$wt_pep,function(x) substr(x[[1]],1,1)))!=TCGA_maf[,"ref_aa"])
length(idx_noMatch) # 245
TCGA_maf<- TCGA_maf[-idx_noMatch,]
idx_noMatch<- which(as.character(sapply(TCGA_maf$wt_pep,function(x) substr(x[[9]],9,9)))!=TCGA_maf[,"ref_aa"]) # In case first are NA!
length(idx_noMatch) # 27
TCGA_maf<- TCGA_maf[-idx_noMatch,]

# saveRDS(TCGA_maf,file="temp/TCGA_maf.rds")

# Add mutated pep to TCGA_maf!!!
TCGA_maf<- cbind(TCGA_maf,mut_pep=NA)
nonapep_mut<- TCGA_maf$wt_pep
pb <- txtProgressBar(min = 0, max = nrow(TCGA_maf), style = 3)
for(i in 1:nrow(TCGA_maf)){
  setTxtProgressBar(pb, i)
  aa_ref_temp<- as.character(TCGA_maf[i,"ref_aa"])
  aa_alt_temp<- as.character(TCGA_maf[i,"alt_aa"])
  nonapep_wt_temp<- unlist(TCGA_maf[i,"wt_pep"])
  for(j in 1:9){ 
    pep_temp<- nonapep_wt_temp[j]
    if(is.na(pep_temp)) next
    if(substr(pep_temp,j,j)!=aa_ref_temp) stop("Check aa")
    else substr(nonapep_mut[[i]][j],j,j)<- aa_alt_temp
  }
}
TCGA_maf$mut_pep<- nonapep_mut

saveRDS(TCGA_maf,file="temp/TCGA_maf.rds")

# Add HLA affinities
####################

# Add variable
TCGA_maf<- cbind(TCGA_maf,
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
nonapep_wt_all<- na.omit(unique(unlist(nonapep_wt)))
write.table(nonapep_wt_all,file = "temp/wt.pep",quote=F,row.names = F,col.names = F)

# Run NetMHC Pan (10x)
system("scripts/other/runNetMHCPan_wt.sh") # Ignore warnings,works!
# Computing time increases linearly with 1) n peptides or 2) n HLA alleles
# Total time duration +/- 54h!

# Fuse results
for(i in 1:10){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/wt_out_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/wt_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/wt_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/wt_out.rds")

# Get affinities
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
  # Save all affinities as well
  aff_wt_9_ls<- c(aff_wt_9_ls,list(aff_wt_9))
}
names(aff_wt_9_ls)<- HLA_alleles
saveRDS(aff_wt_9_ls,file="temp/TCGA_aff_wt_9.rds")

saveRDS(TCGA_maf,file="temp/TCGA_maf.rds")

# 2 MUTATED PEPTIDES

# Get all pep mut & save
nonapep_mut_all<- na.omit(unique(unlist(nonapep_mut)))
nonapep_mut_all<- nonapep_mut_all[!grepl("\\*",nonapep_mut_all)] # Stopgains --> No peptides
nonapep_mut_all_excl<- setdiff(nonapep_mut_all,nonapep_wt_all) # Not needed to calculate if already done for wts
write.table(nonapep_mut_all_excl,file = "temp/mut.pep",quote=F,row.names = F,col.names = F)

# Run NetMHC Pan (10x) --> +/- 81h in total?
system("scripts/other/runNetMHCPan_mut.sh")

# Fuse results
for(i in 1:10){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/mut_out_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/mut_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/mut_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/mut_out.rds")

# Add affinities from wt 
nonapep_mut_all<- na.omit(unique(unlist(TCGA_maf$mut_pep)))
gene_mhc_nM_wt<- readRDS(file = "temp/wt_out.rds")
gene_mhc_nM_wt<- gene_mhc_nM_wt[intersect(rownames(gene_mhc_nM_wt),nonapep_mut_all),]
gene_mhc_nM<- rbind(gene_mhc_nM,gene_mhc_nM_wt)
  
# Get affinities
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
  # Save all affinities as well
  aff_mut_9_ls<- c(aff_mut_9_ls,list(aff_mut_9))
}
names(aff_mut_9_ls)<- HLA_alleles
saveRDS(aff_mut_9_ls,file="temp/TCGA_aff_mut_9.rds")

# Add harmonic means
source("scripts/functions/harmonic_mean.R")
TCGA_maf<- cbind(TCGA_maf,wt_HLA_mean_aff=apply(TCGA_maf[,grepl("wt_HLA",colnames(TCGA_maf))],1,harmonic_mean))
TCGA_maf<- cbind(TCGA_maf,mut_HLA_mean_aff=apply(TCGA_maf[,grepl("mut_HLA",colnames(TCGA_maf))],1,harmonic_mean))

# Save final TCGA_maf
saveRDS(TCGA_maf,file="data/TCGA_maf.rds")

# # # Some checks -> ok
# TCGA_maf<- readRDS(file="data/TCGA_maf.rds")
# TCGA_maf_old<- readRDS("../../immunogenomics/data/TCGA_maf_HLA_mut.rds")
# TCGA_maf_temp<- TCGA_maf[TCGA_maf[,1]=="PTEN",]
# TCGA_maf_old_temp<- TCGA_maf_old[TCGA_maf_old[,1]=="PTEN",]
# View(cbind(TCGA_maf_old_temp[,"wt_HLA_A1_aff"],TCGA_maf_temp[,"wt_HLA_A1_aff"]))

# Add alternative ways to calculate aggregate HLA binding
##########################################################

# # 1) All 6 alleles seperate
# n_mut<- nrow(TCGA_maf)
# TCGA_maf_6<- TCGA_maf[rep(1:n_mut,each=6),]
# wt_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_6[,grepl("wt_HLA_.._aff",colnames(TCGA_maf_6))])*diag(6)[rep(1:6,n_mut),])
# mut_HLA_aff<- rowMaxs(as.matrix(TCGA_maf_6[,grepl("mut_HLA_.._aff",colnames(TCGA_maf_6))])*diag(6)[rep(1:6,n_mut),])
# TCGA_maf_6[,"wt_HLA_mean_aff"]<- wt_HLA_aff
# TCGA_maf_6[,"mut_HLA_mean_aff"]<- mut_HLA_aff
# saveRDS(TCGA_maf_6,file="data/TCGA_maf_6.rds")
# 
