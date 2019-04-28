calculate_HBMR<- function(sample,variant,isHA,subst_type3,exp_mut,cancer){
  
  cancers<- unique(cancer)
  
  # Proportions exp mut
  exp_mut<- prop.table(exp_mut,1)
  
  # Calculate expected HBMR
  mutProb<- prop.table(table(cancer,subst_type3),1)
  HBMR_exp<- NULL
  for(c in cancers){
    prob_cancer<- mutProb[c,]
    exp_muts<- exp_mut*prob_cancer[rownames(exp_mut)]
    isHA_t<- colSums(exp_muts[,,"TRUE"])
    isNonHA_t<- colSums(exp_muts[,,"FALSE"])
    HBMR_exp<- c(HBMR_exp,(isHA_t["nonsynonymous SNV"]/isHA_t["synonymous SNV"])/(isNonHA_t["nonsynonymous SNV"]/isNonHA_t["synonymous SNV"]))
  }
  names(HBMR_exp)<- cancers
  
  # Calculate HBMR
  HBMR_t<- table(cancer,variant,isHA)
  sample_t<- table(cancer[!duplicated(sample)])
  HBMR_matrix<- matrix(NA,length(cancers),16,dimnames = list(cancers,c("n_samples","n-","s-","n+","s+","n","s","HBMR","HBMR_p","HBMR_cilo","HBMR_ciup","HBMR_norm","HBMR_norm_p","HBMR_norm_cilo","HBMR_norm_ciup","HBMR_exp")))
  for(c in cancers){
    
    # General numbers
    HBMR_temp_t<- cbind(c(HBMR_t[c,"synonymous SNV","FALSE"],HBMR_t[c,"nonsynonymous SNV","FALSE"]),c(HBMR_t[c,"synonymous SNV","TRUE"],HBMR_t[c,"nonsynonymous SNV","TRUE"]))
    s_nonHLA<- HBMR_temp_t[1,1]
    n_nonHLA<- HBMR_temp_t[2,1]
    s_HLA<- HBMR_temp_t[1,2]
    n_HLA<- HBMR_temp_t[2,2]
    n<- n_nonHLA + n_HLA
    s<- s_nonHLA + s_HLA
    n_samples<- sample_t[c] 
    
    # HBMR
    HBMR_fish<- fisher.test(HBMR_temp_t)
    HBMR_temp<- HBMR_fish$estimate
    HBMR_temp_ci<- HBMR_fish$conf.int
    HBMR_temp_p<- HBMR_fish$p.value
    
    # Normalized HBMR
    HBMR_norm_fish<- fisher.test(HBMR_temp_t,or = HBMR_exp[c])
    HBMR_norm_temp_p<- HBMR_norm_fish$p.value
    HBMR_norm_temp<- HBMR_temp/HBMR_exp[c]
    HBMR_norm_ci<- HBMR_temp_ci/HBMR_exp[c]
    
    # Add to matrix
    HBMR_matrix[c,]<- c(n_samples,n_nonHLA, s_nonHLA, n_HLA,s_HLA,n,s,HBMR_temp,HBMR_temp_p,HBMR_temp_ci,HBMR_norm_temp,HBMR_norm_temp_p,HBMR_norm_ci,HBMR_exp[c])
  }
  
  # Return matrix
  return(HBMR_matrix)
}


