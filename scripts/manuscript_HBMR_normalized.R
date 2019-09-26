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
  source("scripts/functions/get_subst_types5.R")
  library(svglite)
}

# Add SSB7 & pentaNT information
################################
{
  # Add SSB7
  CpG_conversion_vector<- readRDS("data/SSB7_conversion_vector.rds")
  TCGA_maf$SSB7<- CpG_conversion_vector[as.character(TCGA_maf$subst_type3)]
  
  # Add pentaNT
  TCGA_maf$subst_type5<- get_subst_types5(chr = TCGA_maf$Chromosome,pos = TCGA_maf$Start_Position,alt_allele = TCGA_maf$Tumor_Seq_Allele2)
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

  HBMR_pentaNT<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                                variant = TCGA_maf$Variant_Classification,
                                isHA = TCGA_maf$isHA_wt,
                                subst_type3 = TCGA_maf$subst_type5,
                                exp_mut = subst_type_matrix_ls$pentaNT,
                                cancer =  TCGA_maf$Cancer)
  
  
  # Combine
  # HBMR_norm<- cbind(HBMR_main[,"HBMR"],HBMR_main[,"HBMR_norm"],HBMR_SSB7[,"HBMR_norm"],HBMR_triNT[,"HBMR_norm"])
  # HBMR_norm_CI<- cbind(HBMR_main[,"HBMR_ciup"],HBMR_main[,"HBMR_norm_ciup"],HBMR_SSB7[,"HBMR_norm_ciup"],HBMR_triNT[,"HBMR_norm_ciup"])
  HBMR_norm<- cbind(HBMR_main[,"HBMR"],HBMR_main[,"HBMR_norm"],HBMR_SSB7[,"HBMR_norm"],HBMR_triNT[,"HBMR_norm"],HBMR_pentaNT[,"HBMR_norm"])
  HBMR_norm_CI<- cbind(HBMR_main[,"HBMR_ciup"],HBMR_main[,"HBMR_norm_ciup"],HBMR_SSB7[,"HBMR_norm_ciup"],HBMR_triNT[,"HBMR_norm_ciup"],HBMR_pentaNT[,"HBMR_norm_ciup"])
  
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
    bp<- barplot(t(HBMR_norm),las=2,ylab="HBMR",beside=T,ylim=c(0,1.5))
    # error bars
    arrows(bp,t(HBMR_norm),bp,t(HBMR_norm_CI),lwd = 1.5, angle = 90,code = 2, length = 0.03,xpd=NA)    
    abline(h=1,lty=2)
    # legend("topleft",fill=grey.colors(4),border=NA,legend=c("None","6 main substitution classes","SSB7 model","96 trinucleotide substitution classes"),bty="n")
    legend("topleft",fill=grey.colors(5),border=NA,legend=c("None","main substitution classes","SSB7 model","trinucleotide substitution classes","pentanucleotide substitution classes"),bty="n")
    dev.off()
  }
}

# Compare deviations from 1
############################
{
  HBMR_dev<- abs(HBMR_norm-1)

  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_deviation.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_deviation.svg"))
    par(mfrow=c(1,2))
    vp<- viopoints::viopoints(HBMR_dev,pch=16,col="blue",group.names=c("None","main","SSB7","triNT","pentaNT"),ylab="HBMR deviation from 1",frame.plot=F,horizontal=F)
    for(l in 1:5){
      lines(c(vp[l]-0.2,vp[l]+0.2),rep(colMedians(HBMR_dev)[l],2),lwd=3)
    }
    # Add figure with cancer colors
    for(i in 1:5){
      if(i==1) plot(rep(i,nrow(HBMR_dev)),HBMR_dev[,i],col=c25[rownames(HBMR_dev)],pch=16,xlim=c(1,5),ylab="HBMR deviation from 1",frame.plot=F,xlab=NA)
      else points(rep(i,nrow(HBMR_dev)),HBMR_dev[,i],col=c25[rownames(HBMR_dev)],pch=16)
    }
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
    HBMRnorm_triNT=HBMR_triNT[cancers,"HBMR_norm"], HBMRnorm_triNT_p=HBMR_triNT[cancers,"HBMR_norm_p"], HBMRnorm_triNT_q=p.adjust(HBMR_triNT[cancers,"HBMR_norm_p"],"fdr"),
    HBMRnorm_pentaNT=HBMR_pentaNT[cancers,"HBMR_norm"], HBMRnorm_pentaNT_p=HBMR_pentaNT[cancers,"HBMR_norm_p"], HBMRnorm_pentaNT_q=p.adjust(HBMR_pentaNT[cancers,"HBMR_norm_p"],"fdr")
  )
  saveRDS(HBMR_norm,file = "results/data/HBMR_norm.rds")
}
