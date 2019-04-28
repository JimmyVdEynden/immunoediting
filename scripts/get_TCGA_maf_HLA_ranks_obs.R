##################################
# get_TCGA_maf_HLA_ranks.R"
##################################

# 1) For observed data
########################
{
  # Load TCGA maf file used in study (with main HLA types)
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  HLA_alleles<- c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")
  
  # Extract pct-based affinities
  ###############################
  
  # 1. WT
  #########
  
  # Fuse results
  for(i in 1:10){
    cat(i, " ")
    gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/wt_out_",i,".xls"),n=1),"\t")),"")
    gene_mhc_rank_temp<- read.table(paste0("temp/wt_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","NULL","numeric"),length(gene_mhc_colnames)),c("NULL","NULL")))
    colnames(gene_mhc_rank_temp)<- gene_mhc_colnames
    if(i==1) gene_mhc_rank<- gene_mhc_rank_temp
    else gene_mhc_rank<- cbind(gene_mhc_rank,gene_mhc_rank_temp)
  }
  pep<- read.table(paste0("temp/wt_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
  rownames(gene_mhc_rank)<- pep
  saveRDS(gene_mhc_rank,file = "temp/wt_out_rank.rds")
  
  # Get affinity ranks for maf file
  aff_wt_9_ls<- NULL
  for(i in 1:6){ # 6 HLA alleles
    cat(i," ")
    col_idx<- match(TCGA_maf[,HLA_alleles[i]],colnames(gene_mhc_rank)) # Identify HLA allele
    aff_wt_9<- matrix(NA,nrow(TCGA_maf),9)
    for(j in 1:9){ # 9 peptides
      pep_temp<- sapply(TCGA_maf[,"wt_pep"],function(x) x[[j]])
      row_idx<- match(pep_temp,rownames(gene_mhc_rank)) # Identify pep
      aff_wt_9[,j]<- gene_mhc_rank[cbind(row_idx,col_idx)]
    }
    aff_wt_9_ls<- c(aff_wt_9_ls,list(aff_wt_9))
  }
  names(aff_wt_9_ls)<- HLA_alleles
  saveRDS(aff_wt_9_ls,file="temp/TCGA_affRank_wt_9.rds")
  
  # Add informatio to TCGA_maf
  TCGA_maf[,grep("_aff",colnames(TCGA_maf))]<- NA
  TCGA_maf$wt_HLA_A1_aff<- rowMins(aff_wt_9_ls[["HLA-A1"]],na.rm=T)
  TCGA_maf$wt_HLA_A2_aff<- rowMins(aff_wt_9_ls[["HLA-A2"]],na.rm=T)
  TCGA_maf$wt_HLA_B1_aff<- rowMins(aff_wt_9_ls[["HLA-B1"]],na.rm=T)
  TCGA_maf$wt_HLA_B2_aff<- rowMins(aff_wt_9_ls[["HLA-B2"]],na.rm=T)
  TCGA_maf$wt_HLA_C1_aff<- rowMins(aff_wt_9_ls[["HLA-C1"]],na.rm=T)
  TCGA_maf$wt_HLA_C2_aff<- rowMins(aff_wt_9_ls[["HLA-C2"]],na.rm=T)
  
  # 2. MUT
  #########
  
  # Fuse results
  for(i in 1:10){
    cat(i, " ")
    gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/mut_out_",i,".xls"),n=1),"\t")),"")
    gene_mhc_rank_temp<- read.table(paste0("temp/mut_out_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","NULL","numeric"),length(gene_mhc_colnames)),c("NULL","NULL")))
    colnames(gene_mhc_rank_temp)<- gene_mhc_colnames
    if(i==1) gene_mhc_rank<- gene_mhc_rank_temp
    else gene_mhc_rank<- cbind(gene_mhc_rank,gene_mhc_rank_temp)
  }
  pep<- read.table(paste0("temp/mut_out_",1,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","character","NULL"),rep(c("NULL","NULL","NULL","NULL"),19),c("NULL","NULL")))[,1]
  rownames(gene_mhc_rank)<- pep
  saveRDS(gene_mhc_rank,file = "temp/mut_out_rank.rds")
  # Add affinities from wt 
  nonapep_mut_all<- na.omit(unique(unlist(TCGA_maf$mut_pep)))
  gene_mhc_rank_wt<- readRDS(file = "temp/wt_out_rank.rds")
  gene_mhc_rank_wt<- gene_mhc_rank_wt[intersect(rownames(gene_mhc_rank_wt),nonapep_mut_all),]
  gene_mhc_rank<- rbind(gene_mhc_rank,gene_mhc_rank_wt)
  
  # Get affinity ranks for maf file
  aff_mut_9_ls<- NULL
  for(i in 1:6){ # 6 HLA alleles
    cat(i," ")
    col_idx<- match(TCGA_maf[,HLA_alleles[i]],colnames(gene_mhc_rank)) # Identify HLA allele
    aff_mut_9<- matrix(NA,nrow(TCGA_maf),9)
    for(j in 1:9){ # 9 peptides
      pep_temp<- sapply(TCGA_maf[,"mut_pep"],function(x) x[[j]])
      row_idx<- match(pep_temp,rownames(gene_mhc_rank)) # Identify pep
      aff_mut_9[,j]<- gene_mhc_rank[cbind(row_idx,col_idx)]
    }
    aff_mut_9_ls<- c(aff_mut_9_ls,list(aff_mut_9))
  }
  names(aff_mut_9_ls)<- HLA_alleles
  saveRDS(aff_mut_9_ls,file="temp/TCGA_affRank_mut_9.rds")
  
  # Add informatio to TCGA_maf
  TCGA_maf$mut_HLA_A1_aff<- rowMins(aff_mut_9_ls[["HLA-A1"]],na.rm=T)
  TCGA_maf$mut_HLA_A2_aff<- rowMins(aff_mut_9_ls[["HLA-A2"]],na.rm=T)
  TCGA_maf$mut_HLA_B1_aff<- rowMins(aff_mut_9_ls[["HLA-B1"]],na.rm=T)
  TCGA_maf$mut_HLA_B2_aff<- rowMins(aff_mut_9_ls[["HLA-B2"]],na.rm=T)
  TCGA_maf$mut_HLA_C1_aff<- rowMins(aff_mut_9_ls[["HLA-C1"]],na.rm=T)
  TCGA_maf$mut_HLA_C2_aff<- rowMins(aff_mut_9_ls[["HLA-C2"]],na.rm=T)
  
  # Add harmonic means
  source("scripts/functions/harmonic_mean.R")
  TCGA_maf$wt_HLA_mean_aff<- apply(TCGA_maf[,grepl("wt_HLA",colnames(TCGA_maf))],1,harmonic_mean)
  TCGA_maf$mut_HLA_mean_aff<- apply(TCGA_maf[,grepl("mut_HLA",colnames(TCGA_maf))],1,harmonic_mean)
  
  # Save final TCGA_maf
  saveRDS(TCGA_maf,file="data/TCGA_maf_ranks.rds")
}
