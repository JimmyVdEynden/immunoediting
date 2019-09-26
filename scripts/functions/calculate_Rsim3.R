calculate_Rsim3<- function(sample,variant,isHA,subst_type3,exp_sites,exp_mut,cancer,isPentaNT=F){
  
  # 1. EXPECTED
  
  if(identical(exp_sites,exp_mut)) isSubset=F
  else isSubset=F
  
  # Get ms per cancer & normalize by total n triNT sites
  subst_type_t<- table(cancer,subst_type3)
  n_triNT<- rowSums(exp_sites)
  ms_cancer<- prop.table(t(t(subst_type_t[,names(n_triNT)])/n_triNT),1)
  
  # Get expected propHLA per subst_type3
  exp_nHLA_triNT<- exp_mut[,"nonsynonymous SNV","TRUE"]
  exp_nNonHLA_triNT<- exp_mut[,"nonsynonymous SNV","FALSE"]
  
  # Correct for the fact that previous numbers are based on 10000 equal sites and probabilities are not!!!
  if(isSubset){
    if(isPentaNT){
      n_pentaNT<- n_triNT
      n_triNT<- tapply(n_pentaNT,paste0(substr(names(n_pentaNT),2,4),">",substr(names(n_pentaNT),8,10)),"sum")
      n_triNT<- as.numeric(n_triNT[paste0(substr(names(n_pentaNT),2,4),">",substr(names(n_pentaNT),8,10))])
      names(n_triNT)<- names(n_pentaNT)
    }
    exp_nHLA_triNT<- exp_nHLA_triNT[names(n_triNT)]*n_triNT/10000   
    exp_nNonHLA_triNT<- exp_nNonHLA_triNT[names(n_triNT)]*n_triNT/10000
  }
  
  # Get expected propHLA per cancer
  exp_ratioHLA<- colSums(exp_nHLA_triNT*t(ms_cancer[,names(exp_nHLA_triNT)]))/colSums(exp_nNonHLA_triNT*t(ms_cancer[,names(exp_nNonHLA_triNT)]))
  
  # 2. OBSERVED  

  # Create table with mutation information for each sample
  sample_t<- table(sample,variant,isHA)
  sample_t<- cbind(sample_t[,c("nonsynonymous SNV","synonymous SNV"),"TRUE"],sample_t[,c("nonsynonymous SNV","synonymous SNV"),"FALSE"])
  colnames(sample_t)<- c("n_HLA","s_HLA","n_nonHLA","s_nonHLA")

  # Calculate proportion of mutations in HLA-bining regions 
  ratio_HLA<- sample_t[,"n_HLA"]/sample_t[,"n_nonHLA"]
  
  # 3. Calculate R
  sample_cancer<- as.character(cancer[!duplicated(sample)])
  names(sample_cancer)<- sample[!duplicated(sample)]
  Rsim<- ratio_HLA/exp_ratioHLA[sample_cancer[names(ratio_HLA)]]
  
  # Return
  Rsim_ls<- list(Rsim=Rsim,ratioHLA=ratio_HLA,ratioHLAe=exp_ratioHLA,sample_t=sample_t)
  return(Rsim_ls)
}
