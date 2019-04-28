##############################
# manuscript_dNdS.R
##############################

# Calculate dN/dS values for different cancers, taking into account HLA binding properties using triNT and SSB7 model

# Load
{
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  subst_type_matrix_ls <- readRDS("data/subst_type_matrix_ls.rds")
}

# Libraries & functions
##########################
{
  library(svglite)
  source("scripts/functions/calculate_dNdS_HLA.R")
}

# 0) Preprocess & Filtering
#############################
{
  # Filter 1) Restrict analysis to synonymous and non-synonymous mutations (indels were removed previously, also remove nonsense)
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]
  
  # Turn barcodes in factors
  TCGA_maf[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf[,"Tumor_Sample_Barcode"])
  
  # Define high affinity HLA binders (default cut-off Kd<500nM)
  isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  TCGA_maf<- cbind(TCGA_maf,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
  
  # Add SSB7 subst types
  CpG_conversion_vector<- readRDS("data/SSB7_conversion_vector.rds")
  TCGA_maf$SSB7<- CpG_conversion_vector[as.character(TCGA_maf$subst_type3)]
}

# triNT subst model
##################
{
  cancers<- names(table(TCGA_maf$Cancer))[table(TCGA_maf$Cancer)>10000]
  
  # Calculate
  dNdS_matrix<- calculate_dNdS_HLA(TCGA_maf$Tumor_Sample_Barcode,variant=TCGA_maf$Variant_Classification,isHA=TCGA_maf$isHA_wt,subst_type3=TCGA_maf$subst_type3,exp_mut=subst_type_matrix_ls$triNT,cancer=TCGA_maf$Cancer)
  dNdS_matrix<- dNdS_matrix[cancers,]
  dNdS_matrix<- dNdS_matrix[order(dNdS_matrix[,"dNdS_HLA"]),]
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_dNdS_triNT.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_dNdS_triNT.pdf"),width = 10)
    par(oma=c(6,2,2,2))  
    bp<- barplot(rbind(dNdS_matrix[,"dNdS_HLA"],dNdS_matrix[,"dNdS_nonHLA"]),las=2,ylab="HLA-corrected dN/dS",beside=T,ylim=c(0,1.2),col=c("black","grey"))
    arrows(bp,rbind(dNdS_matrix[,"dNdS_HLA"],dNdS_matrix[,"dNdS_nonHLA"]),bp,rbind(dNdS_matrix[,"dNdS_HLA_ciup"],dNdS_matrix[,"dNdS_nonHLA_ciup"]),lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)  
    abline(h=1)
    dev.off()
  }
}


# SSB7 subst model
##################
{
  # Calculate
  dNdS_matrix_SSB7<- calculate_dNdS_HLA(TCGA_maf$Tumor_Sample_Barcode,variant=TCGA_maf$Variant_Classification,isHA=TCGA_maf$isHA_wt,subst_type3=TCGA_maf$SSB7,exp_mut=subst_type_matrix_ls$SSB7,cancer=TCGA_maf$Cancer)
  dNdS_matrix_SSB7<- dNdS_matrix_SSB7[cancers,]
  dNdS_matrix_SSB7<- dNdS_matrix_SSB7[order(dNdS_matrix_SSB7[,"dNdS_HLA"]),]
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_dNdS_SSB7.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_dNdS_SSB7.pdf"),width = 10)
    par(oma=c(6,2,2,2))  
    bp<- barplot(rbind(dNdS_matrix_SSB7[,"dNdS_HLA"],dNdS_matrix_SSB7[,"dNdS_nonHLA"]),las=2,ylab="HLA-corrected dN/dS",beside=T,ylim=c(0,1.2),col=c("black","grey"))
    arrows(bp,rbind(dNdS_matrix_SSB7[,"dNdS_HLA"],dNdS_matrix_SSB7[,"dNdS_nonHLA"]),bp,rbind(dNdS_matrix_SSB7[,"dNdS_HLA_ciup"],dNdS_matrix_SSB7[,"dNdS_nonHLA_ciup"]),lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)  
    abline(h=1)
    dev.off()
  }
}


# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  dNdS_matrix_triNT<- dNdS_matrix
  save(dNdS_matrix_triNT,dNdS_matrix_SSB7,file = paste0("results/data/",FULL_FILENAME,"Data"))
}
