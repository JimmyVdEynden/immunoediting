# TCGA_maf<- readRDS("data/TCGA_maf.rds")
# sample<- TCGA_maf$Tumor_Sample_Barcode
# variant<- TCGA_maf$Variant_Classification
# isHA<- TCGA_maf$mut_HLA_mean_aff<500
# subst_type3<- TCGA_maf$subst_type3
# exp_mut_ls<- readRDS("temp/subst_type_matrix_ls.rds")
# load("data/TCGA_HLA_types.RData")
# names(exp_mut_ls)<- rownames(TCGA_HLA_types_netMHC)[1:length(exp_mut_ls)]
# cancer<- TCGA_maf$Cancer
# # exp_mut<- readRDS("data/subst_type_matrix_ls.rds")$triNT
# 
# sample = TCGA_maf$Tumor_Sample_Barcode
# variant = TCGA_maf$Variant_Classification
# isHA = TCGA_maf$isHA_mut
# subst_type3 = TCGA_maf$subst_type3
# exp_mut = subst_type_matrix
# exp_mut_ls = subst_type_matrix_ls_perSample
# cancer = TCGA_maf$Cancer

calculate_Rsim_HLASpecific<- function(sample,variant,isHA,subst_type3,exp_mut,exp_mut_ls,cancer,isPentaNT=F){
  
  # 1. EXPECTED RATIO
  
  # Get cancer ids
  sample_cancer<- as.character(cancer[!duplicated(sample)])
  names(sample_cancer)<- substr(sample[!duplicated(sample)],1,12)
  
  # Exclude samples with no information
  idxExcl<- which(sapply(exp_mut_ls,function(x) length(x)==0))
  if(length(idxExcl)!=0) exp_mut_ls<- exp_mut_ls[-idxExcl]

  # Get ms per cancer & normalize by total n triNT sites
  subst_type_t<- table(cancer,subst_type3)
  n_triNT<- rowSums(exp_mut)
  ms_cancer<- prop.table(t(t(subst_type_t[,names(n_triNT)])/n_triNT))
  
  # Get expected n HLA and nonHLA per subst_type3
  exp_nHLA_triNT_matrix<- sapply(exp_mut_ls,function(x) x[,"nonsynonymous SNV","TRUE"])
  exp_nNonHLA_triNT_matrix<- sapply(exp_mut_ls,function(x) x[,"nonsynonymous SNV","FALSE"])

  # Correct for the fact that previous numbers are based on 10000 equal sites and probabilities are not!!!
  if(isPentaNT){
    n_pentaNT<- n_triNT
    n_triNT<- tapply(n_pentaNT,paste0(substr(names(n_pentaNT),2,4),">",substr(names(n_pentaNT),8,10)),"sum")
    n_triNT<- as.numeric(n_triNT[paste0(substr(names(n_pentaNT),2,4),">",substr(names(n_pentaNT),8,10))])
    names(n_triNT)<- names(n_pentaNT)
  }
  exp_nHLA_triNT_matrix<- exp_nHLA_triNT_matrix[names(n_triNT),]*n_triNT/10000   
  exp_nNonHLA_triNT_matrix<- exp_nNonHLA_triNT_matrix[names(n_triNT),]*n_triNT/10000

  # Get expected ratio per cancer
  exp_ratioHLA<- NULL
  common_samples<- intersect(names(exp_mut_ls),names(sample_cancer))
  for(sample_temp in common_samples){
    cancer_temp<- sample_cancer[sample_temp]
    exp_nHLA_triNT<- exp_nHLA_triNT_matrix[,sample_temp]
    exp_nNonHLA_triNT<- exp_nNonHLA_triNT_matrix[,sample_temp]
    exp_ratioHLA_temp<- sum(exp_nHLA_triNT*t(ms_cancer[cancer_temp,names(exp_nHLA_triNT)]))/sum(exp_nNonHLA_triNT*t(ms_cancer[cancer_temp,names(exp_nNonHLA_triNT)]))
    exp_ratioHLA<- c(exp_ratioHLA,exp_ratioHLA_temp)
  }
  names(exp_ratioHLA)<- common_samples

  # 2. OBSERVED RATIO 

  # Create table with mutation information for each sample
  sample<- substr(sample,1,12)
  sample_t<- table(sample,variant,isHA)
  sample_t<- cbind(sample_t[,c("nonsynonymous SNV","synonymous SNV"),"TRUE"],sample_t[,c("nonsynonymous SNV","synonymous SNV"),"FALSE"])
  colnames(sample_t)<- c("n_HLA","s_HLA","n_nonHLA","s_nonHLA")

  # Calculate proportion of mutations in HLA-bining regions 
  ratio_HLA<- sample_t[,"n_HLA"]/sample_t[,"n_nonHLA"]
  
  # # 3. Calculate R
  # names(ratio_HLA)<- substr(names(ratio_HLA),1,12)
  # common_samples<- intersect(names(ratio_HLA),names(exp_ratioHLA))
  # Rsim<- ratio_HLA[common_samples]/exp_ratioHLA[common_samples]
  
  # 3. EXPECTED & OBSERVED numbers
  # names(ratio_HLA)<- substr(names(ratio_HLA),1,12)
  common_samples<- intersect(names(ratio_HLA),names(exp_ratioHLA))
  obs<- sample_t[common_samples,"n_HLA"]
  exp<- sample_t[common_samples,"n_nonHLA"]*exp_ratioHLA[common_samples]
  
  # # 3. Calculate R
  Rsim<- obs/exp
  
  # Return
  Rsim_ls<- list(Rsim=Rsim,ratioHLA=ratio_HLA,ratioHLAe=exp_ratioHLA,obs=obs,exp=exp,sample_t=sample_t)
  return(Rsim_ls)
}

