#################################################
# manuscript_Rsim_HLASpecific.R
#################################################

# Aim: analysis of ratio n+/n- Obs/Exp(Rsim) for HLA specific expectations

# Load & process data
######################
{
  # Observed Mutation data
  if(useSimulated_data) TCGA_maf<- readRDS("data/TCGA_maf_sim.rds")
  else TCGA_maf<- readRDS("data/TCGA_maf.rds")
  TCGA_maf$isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  
  # Focus on cancers from study
  cancers<- sort(names(table(TCGA_maf[,"Cancer"]))[table(TCGA_maf[,"Cancer"])>10000])
  
  # Exclude Non-expressed & driver genes
  TCGA_maf<- TCGA_maf[TCGA_maf[,"mRNA"]>0,] 
  driver_genes<- readRDS("data/CGC_v83.rds")
  TCGA_maf<- TCGA_maf[!TCGA_maf[,1]%in%driver_genes,]
  
  # Expected mut per subst type
  if(useSimulated_data) subst_type_matrix_ls_perSample<- readRDS(file="data/GPPM_subset_subst_type_matrix_sim_mut_per_sample.rds")
  else subst_type_matrix_ls_perSample<- readRDS(file="data/GPPM_subset_subst_type_matrix_mut_per_sample.rds")
  subst_type_matrix<- readRDS("data/subst_type_matrix_ls.rds")$triNT

  # Cancer colors
  source("scripts/functions/cancer_cols.R")
  # library(BSDA) # sign test
}

# Load functions
################
{
  source("scripts/functions/calculate_Rsim_HLASpecific.R")
  library(svglite)
}

# Get metrics per sample
########################
{
  # Rsim (calculated using simulated mutations - observed data)
  Rsim_ls<- calculate_Rsim_HLASpecific(
    sample = TCGA_maf$Tumor_Sample_Barcode,
    variant = TCGA_maf$Variant_Classification,
    isHA = TCGA_maf$isHA_mut,
    subst_type3 = TCGA_maf$subst_type3,
    exp_mut = subst_type_matrix,
    exp_mut_ls = subst_type_matrix_ls_perSample,
    cancer = TCGA_maf$Cancer
  )
  Rsim<- Rsim_ls[["Rsim"]]
  ratioHLA<- Rsim_ls[["ratioHLA"]]
  ratioHLAe<- Rsim_ls[["ratioHLAe"]]
}

# Plot Rsim
############
{
  # Get cancers for each sample 
  cancer_id<- TCGA_maf$Cancer
  names(cancer_id)<- substr(TCGA_maf$Tumor_Sample_Barcode,1,12)
  cancer_id<- cancer_id[!duplicated(names(cancer_id))]
  
  # Calculate median & p values per cancer
  Rsim_med<- sort(tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers],decreasing = T)
  Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[names(Rsim_med)]
  Rsim_q<- p.adjust(Rsim_p,"fdr")
  # Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) SIGN.test(x,md = 1)$p.value)[cancers]

  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Rsim_mut_HLASpec.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Rsim_mut_HLASpec.pdf"))
    par(mfrow=c(1,2))
    boxplot(Rsim~factor(cancer_id[names(Rsim)],levels = names(Rsim_med)),beside = T,las=2,outline=F,horizontal=T,frame.plot=F,col=c25[names(Rsim_med)],staplewex=0,notch=T,xlab="Observed/Expected")
    if(sum(Rsim_q<0.1)!=0) text(rep(2,sum(Rsim_q<0.1)),which(Rsim_q<0.1),paste0("P=",signif(Rsim_p[which(Rsim_q<0.1)],2)),xpd=NA)
    abline(v=1,lty=2)
    dev.off()
  }
  
  results_t<- cbind(median=Rsim_med[cancers],p=Rsim_p[cancers],q=p.adjust(Rsim_p[cancers],"fdr"))
  results_t<- results_t[order(results_t[,"median"]),]
  print.data.frame(as.data.frame(format(results_t,digits=3)))
}

# Compare observed to expected per samples
############
{
  ratioHLA_matrix<- cbind(obs=ratioHLA[names(ratioHLAe)],exp=ratioHLAe)
  ratioHLA_matrix<- ratioHLA_matrix[!is.infinite(ratioHLA_matrix[,"obs"]),]
  ratioHLA_matrix<- ratioHLA_matrix[as.character(cancer_id[rownames(ratioHLA_matrix)])%in%cancers,] # Only focus on cancers studied
  
  # Plot
  ratioHLA_med<- sort(tapply(ratioHLA_matrix[,"obs"],cancer_id[rownames(ratioHLA_matrix)],"median",na.rm=T)[cancers],decreasing = T)
  ratioHLAe_med<- sort(tapply(ratioHLA_matrix[,"exp"],cancer_id[rownames(ratioHLA_matrix)],"median",na.rm=T)[cancers],decreasing = T)

  # Compare
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_ratioHLA_obs_vs_exp_mut_HLASpec.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_ratioHLA_obs_vs_exp_mut_HLASpec.pdf"))
    par(mfrow=c(1,1))
    plot(ratioHLA_matrix,frame.plot=F,col=rgb(0.5,0.5,0.5,0.5),pch=16,xlab="Observed ratio",ylab="Expected ratio",xlim=c(0,2),ylim=c(0,1))
    cor_temp<- cor.test(ratioHLA_matrix[,"obs"],ratioHLA_matrix[,"exp"])
    # plot(log10(ratioHLA_matrix),frame.plot=F,col=rgb(0.5,0.5,0.5,0.5),pch=16,xlab="log10(observed ratio)",ylab="log10(expected ratio)")
    # cor_temp<- cor.test(ratioHLA_matrix[,"obs"],ratioHLA_matrix[,"exp"],method="spearman") 
    legend("topleft",legend = c(paste0("P=",signif(cor_temp$p.value,3)),paste0("r=",signif(cor_temp$estimate,3))),bty="n")
    dev.off()
  }

  # # Per cancer
  # ratioHLA_cancer_cor<- matrix(NA,length(cancers),2,dimnames = list(cancers,c("r","p")))
  # for(cancer in cancers){
  #   ratioHLA_matrix_temp<- ratioHLA_matrix[cancer_id[rownames(ratioHLA_matrix)]==cancer,]
  #   cor_temp<- cor.test(ratioHLA_matrix_temp[,"obs"],ratioHLA_matrix_temp[,"exp"],method="spearman")
  #   ratioHLA_cancer_cor[cancer,]<- c(cor_temp$estimate,cor_temp$p.value)
  # }
  
}


# Save results
###############
{
  Rsim_med_cancer<- results_t
  if(useSimulated_data) save(Rsim,Rsim_med_cancer,ratioHLA,ratioHLAe,ratioHLA_med,ratioHLAe_med,file = "results/data/Rsim_HLASpecific_simData.RData")
  else save(Rsim,Rsim_med_cancer,ratioHLA,ratioHLAe,ratioHLA_med,ratioHLAe_med,file = "results/data/Rsim_HLASpecific.RData")
}



