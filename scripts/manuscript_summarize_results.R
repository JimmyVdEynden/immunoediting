################################
# manuscript_summarize_results.R
#################################

# Aim: Get main results from metrics used in manuscript

# HBMR
#######
{
  # Load
  if(useSimulated_data) TCGA_maf<- readRDS("data/TCGA_maf_sim_mainHLAType.rds")
  else TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  
  # Get HLA binders
  TCGA_maf$isHA<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  
  # Get information on all sites & substitutions
  subst_type_matrix<- readRDS(file="data/subst_type_matrix_ls.rds")$triNT
  
  # Get cancers from study
  cancers<- sort(c('BRCA','COADREAD','LUAD','OV','HNSC','UCEC','LUSC','PRAD','LGG','STAD','SKCM','GBM','BLCA','LIHC','CESC','KIRP','SARC','PAAD','ESCA'))
  
  # Functions         
  source("scripts/functions/calculate_HBMR.R")

  # Calculate
  HBMR<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                        variant = TCGA_maf$Variant_Classification,
                        isHA = TCGA_maf$isHA,
                        subst_type3 = TCGA_maf$subst_type3,
                        exp_mut = subst_type_matrix,
                        cancer =  TCGA_maf$Cancer)[cancers,]
  
  # FDR correction
  HBMR<- cbind(HBMR,HBMR_q=p.adjust(HBMR[,"HBMR_p"],method = "fdr"),HBMR_norm_q=p.adjust(HBMR[,"HBMR_norm_p"],method = "fdr"))
  
  # summarize
  HBMR_sum<- as.data.frame(HBMR[,c("n_samples","n-","s-","n+","s+","n","s","HBMR","HBMR_p","HBMR_q","HBMR_exp","HBMR_norm","HBMR_norm_p","HBMR_norm_q")])
}

# dNHLA/dNnonHLA
##################
{  
  # Prototypical genotype - wt peptides (as HBMR)
  #################################################################
  {
    # Load data: same as above
    # Get HLA binders: same as above

    # Get information on all sites & substitutions
    subst_type_matrix_sites<- readRDS(file="data/subst_type_matrix_ls.rds")$triNT
    
    # Get information on genotypes-specific substitutions
    subst_type_matrix<- readRDS(file="data/GPPM_subset_subst_type_matrix_ls.rds")
    subst_type_matrix<- subst_type_matrix$triNT_wt

    # Functions         
    source("scripts/functions/calculate_Rsim3.R")
    
    # Exclude Non-expressed & driver genes
    TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"mRNA"])&TCGA_maf[,"mRNA"]>0,]
    driver_genes<- readRDS("data/CGC_v83.rds")
    TCGA_maf<- TCGA_maf[!TCGA_maf[,1]%in%driver_genes,]
    
    # Calculate
    Rsim_ls<- calculate_Rsim3(TCGA_maf$Tumor_Sample_Barcode,
                                         TCGA_maf$Variant_Classification,
                                         TCGA_maf$isHA,
                                         as.character(TCGA_maf$subst_type3),
                                         exp_sites = subst_type_matrix_sites,
                                         exp_mut = subst_type_matrix,
                                         cancer = as.character(TCGA_maf$Cancer))
    
    # Get cancers for each sample 
    cancer_id<- as.character(TCGA_maf$Cancer)
    names(cancer_id)<- TCGA_maf$Tumor_Sample_Barcode
    cancer_id<- cancer_id[!duplicated(names(cancer_id))]
    TCGA_cancer_id<- cancer_id
    
    # Median & p values
    Rsim<- Rsim_ls[["Rsim"]]
    Rsim_med<- tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers]
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[cancers]
    Rsim_q<- p.adjust(Rsim_p,method="fdr")
    
    # Get mutations per sample
    sample_t<- Rsim_ls$sample_t
    rownames(sample_t)<- substr(rownames(sample_t),1,15)
    sample_t<- sample_t[,c("n_nonHLA","s_nonHLA","n_HLA","s_HLA")]
    sample_t<- cbind(sample_t,n=rowSums(sample_t[,c("n_nonHLA","n_HLA")]),s=rowSums(sample_t[,c("s_nonHLA","s_HLA")]))
    colnames(sample_t)<- c("n-","s-","n+","s+","n","s")
    
    # Get observed and exp ratios
    obs<- Rsim_ls$ratioHLA
    exp<- Rsim_ls$ratioHLAe # Per cancer!
    exp<- exp[cancer_id[names(obs)]]
    names(exp)<- names(obs)
    
    # Summarize (Genotype-specific expectation model only calculated for mutated affinities!)
    results_t_sample_mainWt<- as.data.frame(cbind(sample_t,obs=obs[rownames(sample_t)],exp=exp[rownames(sample_t)],dNHLA_dNnonHLA=Rsim[rownames(sample_t)]))
    results_t_can_mainWt<- as.data.frame(cbind(dNHLA_dNnonHLA_med=Rsim_med[cancers],dNHLA_dNnonHLA_p=Rsim_p[cancers],dNHLA_dNnonHLA_q=Rsim_q[cancers]))
  }
  
  # Prototypical genotype - mut peptides
  ######################################
  {
    # Load data: same as above
    # Get HLA binders: 
    TCGA_maf$isHA<- TCGA_maf[,"mut_HLA_mean_aff"]<500
    
    # Get information on all sites & substitutions
    subst_type_matrix_sites<- readRDS(file="data/subst_type_matrix_ls.rds")$triNT
    
    # Get information on genotypes-specific substitutions
    subst_type_matrix<- readRDS(file="data/GPPM_subset_subst_type_matrix_ls.rds")
    subst_type_matrix<- subst_type_matrix$triNT_mut
    
    # Calculate
    Rsim_ls<- calculate_Rsim3(TCGA_maf$Tumor_Sample_Barcode,
                              TCGA_maf$Variant_Classification,
                              TCGA_maf$isHA,
                              as.character(TCGA_maf$subst_type3),
                              exp_sites = subst_type_matrix_sites,
                              exp_mut = subst_type_matrix,
                              cancer = as.character(TCGA_maf$Cancer))
    
    # Median & p values
    Rsim<- Rsim_ls[["Rsim"]]
    Rsim_med<- tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers]
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[cancers]
    Rsim_q<- p.adjust(Rsim_p,method="fdr")
    
    # Get mutations per sample
    sample_t<- Rsim_ls$sample_t
    rownames(sample_t)<- substr(rownames(sample_t),1,15)
    sample_t<- sample_t[,c("n_nonHLA","s_nonHLA","n_HLA","s_HLA")]
    sample_t<- cbind(sample_t,n=rowSums(sample_t[,c("n_nonHLA","n_HLA")]),s=rowSums(sample_t[,c("s_nonHLA","s_HLA")]))
    colnames(sample_t)<- c("n-","s-","n+","s+","n","s")
    
    # Get observed and exp ratios
    obs<- Rsim_ls$ratioHLA
    exp<- Rsim_ls$ratioHLAe # Per cancer!
    exp<- exp[cancer_id[names(obs)]]
    names(exp)<- names(obs)
    
    # Summarize (Genotype-specific expectation model only calculated for mutated affinities!)
    results_t_sample_mainMut<- as.data.frame(cbind(sample_t,obs=obs[rownames(sample_t)],exp=exp[rownames(sample_t)],dNHLA_dNnonHLA=Rsim[rownames(sample_t)]))
    results_t_can_mainMut<- as.data.frame(cbind(dNHLA_dNnonHLA_med=Rsim_med[cancers],dNHLA_dNnonHLA_p=Rsim_p[cancers],dNHLA_dNnonHLA_q=Rsim_q[cancers]))
  }
  
  
  # HLA genotype - mutated peptides (as in main analysis)
  #######################################################
  {
    # Load
    if(useSimulated_data) TCGA_maf<- readRDS("data/TCGA_maf_sim.rds")
    else TCGA_maf<- readRDS("data/TCGA_maf.rds")
    
    # Get HLA binders
    TCGA_maf$isHA<- TCGA_maf[,"mut_HLA_mean_aff"]<500
    
    # Get information on all sites & substitutions
    subst_type_matrix_sites<- readRDS(file="data/subst_type_matrix_ls.rds")$triNT
    
    # Get information on genotypes-specific substitutions
    if(useSimulated_data) subst_type_matrix_ls_perSample<- readRDS(file="data/GPPM_subset_subst_type_matrix_sim_mut_per_sample.rds")
    else subst_type_matrix_ls_perSample<- readRDS(file="data/GPPM_subset_subst_type_matrix_mut_per_sample.rds")
    
    # Exclude Non-expressed & driver genes
    TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"mRNA"])&TCGA_maf[,"mRNA"]>0,]
    driver_genes<- readRDS("data/CGC_v83.rds")
    TCGA_maf<- TCGA_maf[!TCGA_maf[,1]%in%driver_genes,]
    
    # Functions         
    source("scripts/functions/calculate_Rsim_HLASpecific.R")
    
    # Calculate
    Rsim_ls<- calculate_Rsim_HLASpecific(TCGA_maf$Tumor_Sample_Barcode,
                                         TCGA_maf$Variant_Classification,
                                         TCGA_maf$isHA,
                                         as.character(TCGA_maf$subst_type3),
                                         exp_mut = subst_type_matrix_sites,
                                         exp_mut_ls = subst_type_matrix_ls_perSample,
                                         cancer = as.character(TCGA_maf$Cancer))
    # Get cancers for each sample 
    cancer_id<- as.character(TCGA_maf$Cancer)
    names(cancer_id)<- substr(TCGA_maf$Tumor_Sample_Barcode,1,12)
    cancer_id<- cancer_id[!duplicated(names(cancer_id))]
    
    # Median & p values
    Rsim<- Rsim_ls[["Rsim"]]
    Rsim_med<- tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers]
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[cancers]
    Rsim_q<- p.adjust(Rsim_p,method="fdr")
    
    # Get mutations per sample
    sample_t<- Rsim_ls$sample_t
    rownames(sample_t)<- substr(rownames(sample_t),1,12)
    sample_t<- sample_t[,c("n_nonHLA","s_nonHLA","n_HLA","s_HLA")]
    sample_t<- cbind(sample_t,n=rowSums(sample_t[,c("n_nonHLA","n_HLA")]),s=rowSums(sample_t[,c("s_nonHLA","s_HLA")]))
    colnames(sample_t)<- c("n-","s-","n+","s+","n","s")
    
    # Get observed and exp ratios
    obs<- Rsim_ls$ratioHLA
    exp<- Rsim_ls$ratioHLAe # Per cancer!
    
    # Summarize (Genotype-specific expectation model only calculated for mutated affinities!)
    results_t_sample<- as.data.frame(cbind(sample_t,obs=obs[rownames(sample_t)],exp=exp[rownames(sample_t)],dNHLA_dNnonHLA=Rsim[rownames(sample_t)]))
    results_t_can<- as.data.frame(cbind(dNHLA_dNnonHLA_med=Rsim_med[cancers],dNHLA_dNnonHLA_p=Rsim_p[cancers],dNHLA_dNnonHLA_q=Rsim_q[cancers]))
  }
}

# Save results
##############
{
  if(useSimulated_data) save(TCGA_cancer_id,HBMR_sum,results_t_sample,results_t_can,results_t_sample_mainWt,results_t_can_mainWt,results_t_sample_mainMut,results_t_can_mainMut,file="results/data/summary_results_sim.RData")
  else save(TCGA_cancer_id,HBMR_sum,results_t_sample,results_t_can,results_t_sample_mainWt,results_t_can_mainWt,results_t_sample_mainMut,results_t_can_mainMut,file="results/data/summary_results.RData")
}
