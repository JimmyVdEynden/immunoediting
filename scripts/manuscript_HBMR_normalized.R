#######################################
# manuscript_HBMR_normalized.R
#######################################

# Load data
############

# Aim: analysis of ratio Obs/Exp % HLA binding non-syn mut (Rsim)

# Load & process data
######################
{
  # Observed Mutation data
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  TCGA_maf$isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  TCGA_maf$isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  # Add SSB7
  CpG_conversion_vector<- readRDS("data/SSB7_conversion_vector.rds")
  TCGA_maf$SSB7<- CpG_conversion_vector[as.character(TCGA_maf$subst_type3)]

  # Focus on cancers from study
  cancers<- sort(names(table(TCGA_maf[,"Cancer"]))[table(TCGA_maf[,"Cancer"])>10000])

  # Expectd mut per subst type
  subst_type_matrix_ls<- readRDS(file="data/subst_type_matrix_ls.rds")
  
  # Cancer colors
  source("scripts/functions/cancer_cols.R")
}

# Load functions
################
{
  source("scripts/functions/calculate_HBMR.R")
  library(svglite)
}

# Calculate HBMR using different normalization models
#######################################################
{
  # Calculate
  HBMR_main<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                           variant = TCGA_maf$Variant_Classification,
                           isHA = TCGA_maf$isHA_wt,
                           subst_type3 = TCGA_maf$subst_type,
                           exp_mut = subst_type_matrix_ls$main,
                           cancer =  TCGA_maf$Cancer)
  
  HBMR_SSB7<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                              variant = TCGA_maf$Variant_Classification,
                              isHA = TCGA_maf$isHA_wt,
                              subst_type3 = TCGA_maf$SSB7,
                              exp_mut = subst_type_matrix_ls$SSB7,
                              cancer =  TCGA_maf$Cancer)
  
  HBMR_triNT<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                              variant = TCGA_maf$Variant_Classification,
                              isHA = TCGA_maf$isHA_wt,
                              subst_type3 = TCGA_maf$subst_type3,
                              exp_mut = subst_type_matrix_ls$triNT,
                              cancer =  TCGA_maf$Cancer)
  
  # Combine
  HBMR_norm<- cbind(HBMR_main[,"HBMR"],HBMR_main[,"HBMR_norm"],HBMR_SSB7[,"HBMR_norm"],HBMR_triNT[,"HBMR_norm"])
  HBMR_norm_CI<- cbind(HBMR_main[,"HBMR_ciup"],HBMR_main[,"HBMR_norm_ciup"],HBMR_SSB7[,"HBMR_norm_ciup"],HBMR_triNT[,"HBMR_norm_ciup"])

  # Restrict to analyzed cancers, sortedon unnormalized HBMR
  cancers_sorted<- names(sort(HBMR_norm[cancers,1]))
  HBMR_norm<- HBMR_norm[cancers_sorted,]
  HBMR_norm_CI<- HBMR_norm_CI[cancers_sorted,]
  
}

# Plot
#######
{
  # Barplot
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_normalized.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_normalized.svg"))
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(t(HBMR_norm),las=2,ylab="HBMR",beside=T,ylim=c(0,1.3))
    # error bars
    arrows(bp,t(HBMR_norm),bp,t(HBMR_norm_CI),lwd = 1.5, angle = 90,code = 2, length = 0.03,xpd=NA)    
    abline(h=1,lty=2)
    legend("topleft",fill=grey.colors(4),border=NA,legend=c("None","6 main substitution classes","SSB7 model","96 trinucleotide substitution classes"),bty="n")
    dev.off()
  }
}

# Save results
###############
{
  HBMR_norm<- cbind(
    HBMR=HBMR_main[cancers,"HBMR"], HBMR_p=HBMR_main[cancers,"HBMR_p"], HBMR_q=p.adjust(HBMR_main[cancers,"HBMR_p"],"fdr"),
    HBMRnorm_main=HBMR_main[cancers,"HBMR_norm"], HBMRnorm_main_p=HBMR_main[cancers,"HBMR_norm_p"], HBMRnorm_main_q=p.adjust(HBMR_main[cancers,"HBMR_norm_p"],"fdr"),
    HBMRnorm_SSB7=HBMR_SSB7[cancers,"HBMR_norm"], HBMRnorm_SSB7_p=HBMR_SSB7[cancers,"HBMR_norm_p"], HBMRnorm_SSB7_q=p.adjust(HBMR_SSB7[cancers,"HBMR_norm_p"],"fdr"),
    HBMRnorm_triNT=HBMR_triNT[cancers,"HBMR_norm"], HBMRnorm_triNT_p=HBMR_triNT[cancers,"HBMR_norm_p"], HBMRnorm_triNT_q=p.adjust(HBMR_triNT[cancers,"HBMR_norm_p"],"fdr")
  )
  saveRDS(HBMR_norm,file = "results/data/HBMR_norm.rds")
}
