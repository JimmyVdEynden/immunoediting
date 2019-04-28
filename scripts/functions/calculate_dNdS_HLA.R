calculate_dNdS_HLA<- function(sample,variant,isHA,subst_type3,exp_mut,cancer){
  
  cancers<- unique(cancer)
  
  # Make sure same number (e.g. 10e6) of mut in each subst class!
  exp_mut<- 1000000*prop.table(exp_mut,1)
  
  # Calculate expected N/S
  mutProb<- prop.table(table(cancer,subst_type3),1)
  NS_HLA<- NULL
  NS_nonHLA<- NULL
  for(c in cancers){
    prob_cancer<- mutProb[c,]
    exp_muts<- exp_mut*prob_cancer[rownames(exp_mut)]
    isHA_t<- colSums(exp_muts[,,"TRUE"])
    isNonHA_t<- colSums(exp_muts[,,"FALSE"])
    NS_HLA<- c(NS_HLA,isHA_t["nonsynonymous SNV"]/isHA_t["synonymous SNV"])
    NS_nonHLA<- c(NS_nonHLA,isNonHA_t["nonsynonymous SNV"]/isNonHA_t["synonymous SNV"])
  }
  names(NS_HLA)<- cancers
  names(NS_nonHLA)<- cancers
  
  # Calculate observed n/s
  can_var_t<- table(cancer,variant,isHA)
  n_HLA<- can_var_t[,"nonsynonymous SNV","TRUE"]   
  n_nonHLA<- can_var_t[,"nonsynonymous SNV","FALSE"]   
  s_HLA<- can_var_t[,"synonymous SNV","TRUE"]   
  s_nonHLA<- can_var_t[,"synonymous SNV","FALSE"]   
  
  # Calculate dN/dS
  dNdS_matrix<- matrix(NA,length(cancers),8,dimnames = list(cancers,c("dNdS_nonHLA","dNdS_nonHLA_p","dNdS_nonHLA_cilo","dNdS_nonHLA_ciup","dNdS_HLA","dNdS_HLA_p","dNdS_HLA_cilo","dNdS_HLA_ciup")))
  for(c in cancers){
    
    # calculate nonHLA
    fish_temp<- fisher.test(cbind(c(n_nonHLA[c],s_nonHLA[c]),c(round(1000000*NS_nonHLA[c]),1000000))) 
    dNdS_temp<- fish_temp$estimate
    dNdS_p_temp<- fish_temp$p.value
    dNdS_ci<- fish_temp$conf.int
    dNdS_matrix[c,1:4]<- c(dNdS_temp,dNdS_p_temp,dNdS_ci[1],dNdS_ci[2])
      
    # calculate HLA
    fish_temp<- fisher.test(cbind(c(n_HLA[c],s_HLA[c]),c(round(1000000*NS_HLA[c]),1000000))) 
    dNdS_temp<- fish_temp$estimate
    dNdS_p_temp<- fish_temp$p.value
    dNdS_ci<- fish_temp$conf.int
    dNdS_matrix[c,5:8]<- c(dNdS_temp,dNdS_p_temp,dNdS_ci[1],dNdS_ci[2])
  }
  
  # Return matrix
  return(dNdS_matrix)
}

