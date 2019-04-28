###################################
# create_GPPM_wt.R
###################################

GPPM<- readRDS("temp/GPPM_all.rds")

# Get peptides
nonapep_GPPM_all<- na.omit(unique(unlist(GPPM$nonaPep_wt)))
length(nonapep_GPPM_all) 
write.table(nonapep_GPPM_all,file = "temp/GPPM_wt.pep",quote=F,row.names = F,col.names = F)

# Most frequent HLA alleles from TCGA
load("data/TCGA_HLA_types.RData")
HLA_A_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-A1"],TCGA_HLA_types_netMHC[,"HLA-A2"])),decreasing = T)
HLA_B_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-B1"],TCGA_HLA_types_netMHC[,"HLA-B2"])),decreasing = T)
HLA_C_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-C1"],TCGA_HLA_types_netMHC[,"HLA-C2"])),decreasing = T)
HLA_freq<- c(HLA_A_t[1:2],HLA_B_t[1:2],HLA_C_t[1:2])
HLA_freq<- HLA_freq/nrow(TCGA_HLA_types_netMHC)
for(i in 1:6){cat(names(HLA_freq)[i],HLA_freq[i],"\n")}
# HLA-A02:01 0.4351026 
# HLA-A01:01 0.2741971 
# HLA-B07:02 0.2255798 
# HLA-B08:01 0.1725022 
# HLA-C07:01 0.2614853 
# HLA-C07:02 0.2569135 
write.table(names(HLA_freq),"temp/HLA_types_freq.txt",sep=",",quote=FALSE,row.names = F,col.names = F)

# Run NetMHC Pan (10x)
system("scripts/other/runNetMHCPan_GPPM_wt.sh") # Ignore warnings,works!

# Fuse results
for(i in 1:6){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/GPPM_wt_out_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/GPPM_wt_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
pep<- read.table(paste0("temp/GPPM_wt_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),1),c("NULL","NULL")))[,1]
rownames(gene_mhc_nM)<- pep
saveRDS(gene_mhc_nM,file = "temp/GPPM_wt_out.rds")

# Get affinities
# aff_wt_9_ls<- NULL
aff_wt<- matrix(NA,length(GPPM),6,dimnames=list(NULL,c("A1","A2","B1","B2","C1","C2")))
for(i in 1:6){ # 6 HLA alleles
  cat("\n",i," ")
  col_idx<- i # Identify HLA allele
  aff_wt_9<- matrix(NA,length(GPPM),9)
  for(j in 1:9){ # 9 peptides
    cat(j, " ")
    pep_temp<- sapply(GPPM$nonaPep_wt,function(x) x[[j]])
    row_idx<- match(pep_temp,rownames(gene_mhc_nM)) # Identify pep
    aff_wt_9[,j]<- gene_mhc_nM[cbind(row_idx,col_idx)]
  }
  aff_wt[,i]<- rowMins(aff_wt_9,na.rm=T)

}
source("scripts/functions/harmonic_mean.R")
wt_HLA_mean_aff<- apply(aff_wt,1,harmonic_mean)
saveRDS(aff_wt,file = "temp/GPPM_all_wt_HLA_aff.rds")
saveRDS(wt_HLA_mean_aff,file = "temp/GPPM_all_wt_mean_HLA_aff.rds")

# Add to GPPM
GPPM$mean_HLA_aff_wt<- wt_HLA_mean_aff

# Add variants using annovar
############################

Annovar_input<- cbind(Chromosome=gsub("chr","",as.character(seqnames(GPPM))),Start_Position=pos(GPPM),End_Position=pos(GPPM),Reference_Allele=as.character(GPPM$ref),Tumor_Seq_Allele2=as.character(GPPM$alt))

# Split in 3 to avoid overload (otherwise error!)
write.table(Annovar_input[1:30000000,],"temp/GPPM_Annovar_input1.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(Annovar_input[30000001:60000000,],"temp/GPPM_Annovar_input2.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(Annovar_input[60000001:nrow(Annovar_input),],"temp/GPPM_Annovar_input3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
for(i in 1:3){
  system(paste0("perl ~/tools/annovar_old_version/table_annovar.pl temp/GPPM_Annovar_input",i,".txt ~/tools/annovar_old_version/humandb/ -buildver hg19 -out temp/GPPM_myanno",i," -remove -protocol refGene,ljb26_all -operation g,f -nastring ."))
}

# Process annovar output
sim_avo<- NULL
for(i in 1:3){
  cat(i," ")
  sim_avo_temp<- read.table(paste0("temp/GPPM_myanno",i,".hg19_multianno.txt"),header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")
  sim_avo_temp<- sim_avo_temp[,c("Gene.refGene","ExonicFunc.refGene","Polyphen2_HVAR_score","AAChange.refGene")]
  sim_avo<- rbind(sim_avo,sim_avo_temp)
}
colnames(sim_avo)<-c("Hugo_Symbol","Variant_Classification","PP2","Amino_acids")
sim_avo[,"Variant_Classification"]<- sapply(sim_avo[,"Variant_Classification"], function(x) strsplit(x, ";")[[1]][1]) # sapply much faster than loop!!!
avo_sil_idx<- which(sim_avo[,"Variant_Classification"]=="synonymous SNV")
sim_avo[avo_sil_idx,"PP2"]<- 0
avo_trunc_idx<- which(sim_avo[,"Variant_Classification"]=="stopgain")
sim_avo[avo_trunc_idx,"PP2"]<- 1
sim_avo[,"Amino_acids"]<- gsub("[[:digit:]]+","/",gsub(".*p\\.","",sim_avo[,"Amino_acids"]))
sim_avo[,"Amino_acids"]<- gsub("X","\\*",sim_avo[,"Amino_acids"]) # same annotation as maf data

GPPM$variant<- sim_avo[,"Variant_Classification"]
GPPM$ref_aa<- substr(sim_avo[,"Amino_acids"],1,1)
GPPM$alt_aa<- substr(sim_avo[,"Amino_acids"],3,3)
GPPM$PP2<- as.numeric(sim_avo[,"PP2"])

# Remove lines where amino acid doesn't match
idx_noMatch<- union(
  which(as.character(sapply(GPPM$nonaPep_wt,function(x) substr(x[[1]],1,1)))!=GPPM$ref_aa),
  which(as.character(sapply(GPPM$nonaPep_wt,function(x) substr(x[[9]],9,9)))!=GPPM$ref_aa) # In case first are NA!
)
length(idx_noMatch) # 632830
GPPM<- GPPM[-idx_noMatch]

# Remove lines where amino acid not determined by annovar
idx_noAA<- which(is.na(GPPM$ref_aa))
length(idx_noAA) # 432801
GPPM<- GPPM[-idx_noAA]

saveRDS(GPPM,file="data/GPPM_inclHLA.rds")

# # Also add information on mutated peptides?
# ##########################################
# 
# # Would take months to calculate with earlier approach, following only takes on hour!
# GPPM<- readRDS(file="data/GPPM_inclHLA.rds")
# aa_ref_temp<- GPPM$ref_aa
# aa_alt_temp<- GPPM$alt_aa
# start_time<- Sys.time()
# nonapep_mut_ls<- list()
# for(i in 1:9){
#   cat(i," ")
#   nonapep_mut_temp<- sapply(GPPM$nonaPep_wt,function(x) x[[i]])
#   substr(nonapep_mut_temp,i,i)<- aa_alt_temp
#   nonapep_mut_ls<- c(nonapep_mut_ls,list(nonapep_mut_temp))
# }
# Sys.time()-start_time
# saveRDS(nonapep_mut_ls,file = "temp/nonapep_mut_ls.rds")
# 
# # Synonymous should remain the same, check!
# # ok
# idx_syn<- which(GPPM$variant=="synonymous SNV")
# # i<- 1
# i<- 9
# nonapep_wt_temp<- sapply(GPPM$nonaPep_wt,function(x) x[[i]])
# nonapep_mut_temp<- nonapep_mut_ls[[i]]
# idx_nonsim<- which(nonapep_mut_temp==nonapep_wt_temp)
# cat(sum(!idx_nonsim%in%idx_syn)==0)
# 
# Put in same format as wt

# NOT A GOOD SOLUTION YET. CHECK LATER WHETHER NECESSARY TO CALCULATE?

# apply(cbind(nonapep_mut_temp[1:100],nonapep_mut_temp[1:100]),MARGIN = 1,function(x) list(x))
# nonapep_mut_ls_temp<- list(
#   nonapep_mut_ls[[1]][1:100],
#   nonapep_mut_ls[[2]][1:100],
#   nonapep_mut_ls[[3]][1:100],
#   nonapep_mut_ls[[4]][1:100],
#   nonapep_mut_ls[[5]][1:100],
#   nonapep_mut_ls[[6]][1:100],
#   nonapep_mut_ls[[7]][1:100],
#   nonapep_mut_ls[[8]][1:100],
#   nonapep_mut_ls[[9]][1:100]
# )
#
# sapply(nonapep_mut_ls_temp,function(x) list(x[1],x[2]))

# Get peptides
nonapep_mut_ls<- readRDS(file = "temp/nonapep_mut_ls.rds")
# object.size(nonapep_mut_ls) # +/- 41GB
nonapep_mut_unique<- setdiff(unique(unlist(nonapep_mut_ls)),NA)
length(nonapep_mut_unique) # 536,376,430 --> Too many to determine? Only if really necessary
# # Test duration and memory for 1 million: 45' on 6 alleles
# write.table(nonapep_mut_unique[1:1000000],file = "temp/GPPM_mut_temp1.pep",quote=F,row.names = F,col.names = F)
# SPlit in 10 seperate pep files (5 crashed?)
l<- length(nonapep_mut_unique)
stop_pos<- 0
for(i in 1:10){
  cat(i," ")
  start_pos<- stop_pos+1 
  stop_pos<- ceiling(i/10*l)
  # cat(start_pos," ",stop_pos,"-")
  write.table(nonapep_mut_unique[start_pos:stop_pos],file = paste0("temp/GPPM_mut_",i,".pep"),quote=F,row.names = F,col.names = F)
}

# write.table(nonapep_mut_unique,file = "temp/GPPM_mut.pep",quote=F,row.names = F,col.names = F)

# # # Run NetMHC Pan: 536,376,430 peptides --> Too many to determine? Only if really necessary
# # system("scripts/other/runNetMHCPan_GPPM_mut.sh") # Ignore warnings,works!

# # test
# i=""
# gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/GPPM_mut_temp1_out_",i,".xls"),n=1),"\t")),"")
# gene_mhc_nM_temp<- read.table(paste0("temp/GPPM_mut_temp1_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
# colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
# object.size(gene_mhc_nM_temp) #48MB

# Add HLA-specific information to GPPM (currently only contains harmonic mean!)
###############################################################################

GPPM<- readRDS("temp/GPPM_all.rds")
aff_wt<- readRDS(file = "temp/GPPM_all_wt_HLA_aff.rds")
wt_HLA_mean_aff<- readRDS(file = "temp/GPPM_all_wt_mean_HLA_aff.rds")
GPPM$HLA_aff_A1<- aff_wt[,1]
GPPM$HLA_aff_A2<- aff_wt[,2]
GPPM$HLA_aff_B1<- aff_wt[,3]
GPPM$HLA_aff_B2<- aff_wt[,4]
GPPM$HLA_aff_C1<- aff_wt[,5]
GPPM$HLA_aff_C2<- aff_wt[,6]
GPPM$HLA_aff_mean<- wt_HLA_mean_aff

# Add annovar annotations
sim_avo<- NULL
for(i in 1:3){
  cat(i," ")
  sim_avo_temp<- read.table(paste0("temp/GPPM_myanno",i,".hg19_multianno.txt"),header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")
  sim_avo_temp<- sim_avo_temp[,c("Gene.refGene","ExonicFunc.refGene","Polyphen2_HVAR_score","AAChange.refGene")]
  sim_avo<- rbind(sim_avo,sim_avo_temp)
}
colnames(sim_avo)<-c("Hugo_Symbol","Variant_Classification","PP2","Amino_acids")
sim_avo[,"Variant_Classification"]<- sapply(sim_avo[,"Variant_Classification"], function(x) strsplit(x, ";")[[1]][1]) # sapply much faster than loop!!!
avo_sil_idx<- which(sim_avo[,"Variant_Classification"]=="synonymous SNV")
sim_avo[avo_sil_idx,"PP2"]<- 0
avo_trunc_idx<- which(sim_avo[,"Variant_Classification"]=="stopgain")
sim_avo[avo_trunc_idx,"PP2"]<- 1
sim_avo[,"Amino_acids"]<- gsub("[[:digit:]]+","/",gsub(".*p\\.","",sim_avo[,"Amino_acids"]))
sim_avo[,"Amino_acids"]<- gsub("X","\\*",sim_avo[,"Amino_acids"]) # same annotation as maf data

GPPM$variant<- sim_avo[,"Variant_Classification"]
GPPM$ref_aa<- substr(sim_avo[,"Amino_acids"],1,1)
GPPM$alt_aa<- substr(sim_avo[,"Amino_acids"],3,3)
GPPM$PP2<- as.numeric(sim_avo[,"PP2"])

# # Remove lines where amino acid doesn't match
# idx_noMatch<- union(
#   which(as.character(sapply(GPPM$nonaPep_wt,function(x) substr(x[[1]],1,1)))!=GPPM$ref_aa),
#   which(as.character(sapply(GPPM$nonaPep_wt,function(x) substr(x[[9]],9,9)))!=GPPM$ref_aa) # In case first are NA!
# )
# length(idx_noMatch) # 632830
# GPPM<- GPPM[-idx_noMatch]
# 
# # Remove lines where amino acid not determined by annovar
# idx_noAA<- which(is.na(GPPM$ref_aa))
# length(idx_noAA) # 432801
# GPPM<- GPPM[-idx_noAA]

saveRDS(GPPM,file = "data/GPPM_inclHLAAlleles.rds")

