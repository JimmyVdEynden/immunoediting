calculate_R<- function(sample,variant,isHA,subst_type3,isExpr=NULL){
  
  # 1. EXPECTED
  
  # Expected N/S: Npred
  NS_t<- table(subst_type3,variant)
  NS<- NS_t[,"nonsynonymous SNV"]/NS_t[,"synonymous SNV"]
  NS<- NS[!is.infinite(NS)] # 4 mutations can never result in silent mutation, remove!
  
  sample_s_t<- table(sample,subst_type3,variant)
  sample_s_t<- sample_s_t[,,"synonymous SNV"]
  # N_pred<- rowSums(sample_s_t[,names(NS)]*NS,na.rm=T) # MISTAKE!!!!
  N_pred<- rowSums(t(NS*t(sample_s_t[,names(NS)])),na.rm=T)
  
  # Expected prop binder: Bpred
  B_t<- table(subst_type3,variant,isHA)
  B<- B_t[,"nonsynonymous SNV","TRUE"]/NS_t[,"nonsynonymous SNV"]
  B_pred<- rowSums(t(NS*t(sample_s_t[,names(NS)])*B[names(NS)]),na.rm=T)
  
  # 2. OBSERVED  
  if(!is.null(isExpr)){
    sample<- sample[isExpr]
    variant<- variant[isExpr]
    isHA<- isHA[isExpr]
    subst_type3<- subst_type3[isExpr]
  }
  
  # Create table with mutation information for each sample
  sample_t<- table(sample,variant,isHA)
  sample_t<- cbind(sample_t[,c("nonsynonymous SNV","synonymous SNV"),"TRUE"],sample_t[,c("nonsynonymous SNV","synonymous SNV"),"FALSE"])
  colnames(sample_t)<- c("n_HLA","s_HLA","n_nonHLA","s_nonHLA")
  sample_t<- cbind(sample_t,n=(sample_t[,"n_HLA"]+sample_t[,"n_nonHLA"]),s=(sample_t[,"s_HLA"]+sample_t[,"s_nonHLA"]))
  
  # Calculate proportion of mutations in HLA-bining regions 
  prop_HLA<- prop.table(sample_t[,c("n_HLA","n_nonHLA")],1)[,"n_HLA"]
  
  # 3. Calculate R
  R<- prop_HLA/(B_pred/N_pred)[names(prop_HLA)]
  
  # Return
  R_ls<- list(R=R,propHLA=prop_HLA,Bpred=B_pred,Npred=N_pred)
  return(R_ls)
}
