#########################################
# manuscript_IEDB.R
#########################################

# Load
#######
{
  # Need GPPM data to calculate expected N/S rates 
  GPPM<- readRDS(file="data/GPPM_inclHLAAlleles.rds") # Takes some time to load
  # Need main and simulated mutation databases to determine whether any signals of neo-antigen depletion can be detected
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  # IEDB data
  epi<- readRDS("data/IEDB_epitopes.rds")
  epi_gr<- readRDS(file="data/IEDB_epitopes_gr.rds")
  # GPPM subst matrix
  subst_type_matrix_ls<- readRDS("data/subst_type_matrix_ls.rds")
}

# Libraries & functions
######################
{
  library(svglite)
  source("scripts/functions/convert_maf_to_GPos.R")
  source("scripts/functions/calculate_HBMR.R")
}

# Preprocess and filtering
###########################
{
  # Obs
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]
  TCGA_maf[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf[,"Tumor_Sample_Barcode"])
  isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  TCGA_maf<- cbind(TCGA_maf,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
  
  # GPPM
  idx_NS<- which(GPPM$variant%in%c("nonsynonymous SNV","synonymous SNV"))
  GPPM<- GPPM[idx_NS]
  GPPM$isHA<- GPPM$HLA_aff_mean<500
}

# Some data about epitopes
##########################
{
  # How many mapped?
  epi_ENST<- epi[,3]
  epi_start<- as.numeric(gsub("\\-.*","",epi[,6]))
  epi_end<- as.numeric(gsub(".*\\-","",epi[,6]))
  ENST_unmapped_unique<- unique(paste0(epi_ENST,"_",epi_start,"_",epi_end))
  cat("Total number of epitopes used:",length(ENST_unmapped_unique),"  \n") # 66698
  ENST_mapped_unique<- unique(paste0(epi_gr$tx_id,"_",epi_gr$protein_start,"_",epi_gr$protein_end))
  cat("Total number of epitopes mapped:",length(ENST_mapped_unique),"  \n") # 66536
  cat("Proportion of epitopes mapped:",length(ENST_mapped_unique)/length(ENST_unmapped_unique),"  \n") # 99.8ù

  # Length?
  epi_length<- epi_gr$protein_end-epi_gr$protein_start+1
  epi_gr$protein_length<- epi_length
  idx_above12<- which(epi_length>12)
  idx_below8<- which(epi_length<8)
  epi_length[idx_above12]<- ">12"
  epi_length[idx_below8]<- "<8"
  epi_length_t<- table(epi_length)
  epi_length_t<- epi_length_t[c("<8","8","9","10","11","12",">12")]
  print.data.frame(cbind(as.data.frame(epi_length_t),prop=as.numeric(prop.table(epi_length_t))))
  # epi_length  Freq        prop
  # 1         <8   347 0.004429128
  # 2          8  4329 0.055255600
  # 3          9 39245 0.500925394
  # 4         10 13393 0.170949008
  # 5         11  8639 0.110268683
  # 6         12  3317 0.042338375
  # 7        >12  9075 0.115833812
  
  # Create pie chart of length
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_length_pie.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_length_pie.pdf"),width = 10)
    pie(epi_length_t,labels = paste0(names(epi_length_t),"-mers"),main="IEDB epitopes")
    dev.off()
  }
}

# Map epitopes to mutation data & general description
######################################################
{
  # Convert maf to GPos format
  metadata_cols<- c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","subst_type","subst_type3","wt_HLA_mean_aff","mut_HLA_mean_aff","isHA_wt","isHA_mut","Cancer")
  metadata_names<- c("gene","sample","variant","subst_type","subst_type3","HLA_aff_wt","HLA_aff_mut","isHA_wt","isHA_mut","cancer")
  TCGA_maf_GP<- Convert_maf_to_GPos.R(maf = TCGA_maf,metadata_cols = metadata_cols, metadata_names = metadata_names)

  # Define epitope regions
  TCGA_maf_GP_epi_idx<- findOverlaps(TCGA_maf_GP,epi_gr,ignore.strand=T)
  TCGA_maf_GP$isEpi<- (1:length(TCGA_maf_GP))%in%queryHits(TCGA_maf_GP_epi_idx)

  # # Some numbers about epitope regions
  # mean(TCGA_maf_GP$isEpi) # 0.052
  # tapply(TCGA_maf_GP$HLA_aff_wt,TCGA_maf_GP$isEpi,"median") # 2631.641 1441.639 
  # tapply(TCGA_maf_GP$HLA_aff_wt,TCGA_maf_GP$isEpi,"median") # 2631.641 1441.639 
  # tapply(TCGA_maf_GP$isHA_wt,TCGA_maf_GP$isEpi,"mean") # 0.2201506 0.3198404 

  # Are these numbers different if only nonapeptides are considered? YES, very clear
  TCGA_maf_GP_epi_nona_idx<- findOverlaps(TCGA_maf_GP,epi_gr[epi_gr$protein_length=="9"],ignore.strand=T)
  TCGA_maf_GP$isEpiNona<- (1:length(TCGA_maf_GP))%in%queryHits(TCGA_maf_GP_epi_nona_idx)
  tapply(TCGA_maf_GP$HLA_aff_wt,TCGA_maf_GP$isEpiNona,"median") # 2622.2579  769.4802
  epi_t<- table(isHA=TCGA_maf_GP$isHA_wt,isEpi=TCGA_maf_GP$isEpi)
  epiNona_t<- table(isHA=TCGA_maf_GP[TCGA_maf_GP$isEpi]$isHA_wt,isEpiNona=TCGA_maf_GP[TCGA_maf_GP$isEpi]$isEpiNona)
  epi_t<- cbind(epi_t[,1],epiNona_t[,c("TRUE","FALSE")])
  colnames(epi_t)<- c("No","9","Other")
  p_nona<- fisher.test(epi_t[,c("No","9")])$p.value
  p_other<- fisher.test(epi_t[,c("No","Other")])$p.value
  epi_t_prop<- prop.table(epi_t,2)["TRUE",]
  print.data.frame(100*as.data.frame(epi_t_prop),digits = 3,row.names = c("No epitope","9-mer epitope","other epitope"))
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_HLABinding.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_HLABinding.pdf"),width = 10)
    bp<- barplot(epi_t_prop,names.arg = c("No epitope","9-mer epitope","other epitope"),ylab="% HLA binding")
    text(bp[2],0.35,labels = paste0("P=",format(p_nona,digits=3)),xpd=NA)
    text(bp[3],0.35,labels = paste0("P=",format(p_other,digits=3)),xpd=NA)
    dev.off()
  }
  
  # Doesn't this plot make more sense when mapped to GPPM? Than it is completely linked to annotation, now to mutations, which is confusing!
  epi_gr2<- epi_gr 
  seqlevels(epi_gr2)<- paste0("chr",seqlevels(epi_gr2)) # Match seqnames
  genome(epi_gr2)<- NA # Match genome build (hg19 GRCh37)
  GPPM_epi_idx<- findOverlaps(GPPM,epi_gr2,ignore.strand=T) # Porblems for noncommon chromosomes!
  GPPM$isEpi<- (1:length(GPPM))%in%queryHits(GPPM_epi_idx)
  # mean(GPPM$isEpi) # 0.058 Comparable
  # tapply(GPPM$isHA,GPPM$isEpi,"mean") # 32.5% vs 21.5% 
  GPPM_epi_nona_idx<- findOverlaps(GPPM,epi_gr2[epi_gr2$protein_length=="9"],ignore.strand=T)
  GPPM$isEpiNona<- (1:length(GPPM))%in%queryHits(GPPM_epi_nona_idx)
  epi_t<- table(isHA=GPPM$isHA,isEpi=GPPM$isEpi)
  epiNona_t<- table(isHA=GPPM[GPPM$isEpi]$isHA,isEpiNona=GPPM[GPPM$isEpi]$isEpiNona)
  epi_t<- cbind(epi_t[,1],epiNona_t[,c("TRUE","FALSE")])
  colnames(epi_t)<- c("No","9","Other")
  p_nona<- fisher.test(epi_t[,c("No","9")])$p.value
  p_other<- fisher.test(epi_t[,c("No","Other")])$p.value
  epi_t_prop<- prop.table(epi_t,2)["TRUE",]
  print.data.frame(100*as.data.frame(epi_t_prop),digits = 3,row.names = c("No epitope","9-mer epitope","other epitope"))
  
  # Plot proportion of genome thas has epitopes mapped
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_mapped.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_mapped.pdf"),width = 10)
    epitope_mapped_t<- table(GPPM$isEpi)
    pie(epitope_mapped_t,labels = c(NA,"Mapped epitopes"),main="IEDB")
    text(0,0,paste0(100*signif(prop.table(epitope_mapped_t)["TRUE"],2),"%"),adj = c(0,0))
    dev.off()
  }
    
  # Plot proportion of binder per epitope type
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_HLABinding_GPPM.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_HLABinding_GPPM.pdf"),width = 10)
    bp<- barplot(prop.table(epi_t,2)[c("TRUE","FALSE"),],names.arg = c("No epitope","9-mer epitope","other epitope"),axes=F)
    dev.off()
  }
}

# HBMR calculations based on epitopes: called EMR
##################################################
{
  # Calculate for observed data
  cancers<- sort(names(table(TCGA_maf[,"Cancer"])[table(TCGA_maf[,"Cancer"])>10000])) # Only the ones with a miniml number of mutations

  obs_EMR_triNT<- calculate_HBMR(TCGA_maf_GP$sample,TCGA_maf_GP$variant,TCGA_maf_GP$isEpi,TCGA_maf_GP$subst_type3,subst_type_matrix_ls$triNT,TCGA_maf_GP$cancer)
  obs_EMR_isEpi<- obs_EMR_triNT[cancers,c("HBMR_p","HBMR","HBMR_ciup")]
  colnames(obs_EMR_isEpi)<- c("p","EMR","CI")
  
  # sort
  obs_EMR_isEpi<- obs_EMR_isEpi[order(obs_EMR_isEpi[,"EMR"]),]
  
  # Print results
  print.data.frame(as.data.frame(format(obs_EMR_isEpi,digits=3)))

  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_EMR_obs.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_EMR_obs.pdf"),width = 10)
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(obs_EMR_isEpi[,"EMR"],las=2,ylab="EMR",beside=T,ylim=c(0,1.3))
    # error bars
    arrows(bp,obs_EMR_isEpi[,"EMR"],bp,obs_EMR_isEpi[,"CI"],lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)    
    abline(h=1,lty=2)
    dev.off()
  }
  
}

# EMR in binding vs non-binding epitopes?
##########################################
{
  # Calculate for HLA binding
  obs_EMR_B_triNT<- calculate_HBMR(TCGA_maf_GP$sample,TCGA_maf_GP$variant,TCGA_maf_GP$isEpi&TCGA_maf_GP$isHA_wt,TCGA_maf_GP$subst_type3,subst_type_matrix_ls$triNT,TCGA_maf_GP$cancer)
  obs_EMR_isEpi_B<- obs_EMR_B_triNT[cancers,c("HBMR_p","HBMR","HBMR_ciup")]
  colnames(obs_EMR_isEpi_B)<- c("p","EMR","CI")

  # sort
  obs_EMR_isEpi_B<- obs_EMR_isEpi_B[order(obs_EMR_isEpi_B[,"EMR"]),]
  
  # Print results
  print.data.frame(as.data.frame(format(obs_EMR_isEpi_B,digits=3)))
  
  # Calculate for HLA non-binding
  obs_EMR_nB_triNT<- calculate_HBMR(TCGA_maf_GP$sample,TCGA_maf_GP$variant,TCGA_maf_GP$isEpi&!TCGA_maf_GP$isHA_wt,TCGA_maf_GP$subst_type3,subst_type_matrix_ls$triNT,TCGA_maf_GP$cancer)
  obs_EMR_isEpi_nB<- obs_EMR_nB_triNT[cancers,c("HBMR_p","HBMR","HBMR_ciup")]
  colnames(obs_EMR_isEpi_nB)<- c("p","EMR","CI")

  # sort
  obs_EMR_isEpi_nB<- obs_EMR_isEpi_nB[order(obs_EMR_isEpi_nB[,"EMR"]),]
  
  # Print results
  print.data.frame(as.data.frame(format(obs_EMR_isEpi_nB,digits=3)))
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_epitope_EMR_isHA.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_epitope_EMR_isHA.pdf"),width = 10)
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(rbind(obs_EMR_isEpi_B[,"EMR"],obs_EMR_isEpi_nB[rownames(obs_EMR_isEpi_B),"EMR"]),las=2,ylab="EMR",beside=T,ylim=c(0,1.5))
    # error bars
    arrows(bp,rbind(obs_EMR_isEpi_B[,"EMR"],obs_EMR_isEpi_nB[rownames(obs_EMR_isEpi_B),"EMR"]),bp,rbind(obs_EMR_isEpi_B[,"CI"],obs_EMR_isEpi_nB[rownames(obs_EMR_isEpi_B),"CI"]),lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)    
    abline(h=1,lty=2)
    dev.off()
  }
}

# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  save(epi_t,obs_EMR_isEpi,obs_EMR_isEpi_B,obs_EMR_isEpi_nB,file=paste0("results/data/",FULL_FILENAME,"Data"))
}

