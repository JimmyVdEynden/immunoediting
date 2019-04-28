calculate_dNdS<- function(gene,variant,subst_type3,exp_mut){
  
  # Get mutation probability (normalized for triNT sites)
  n_triNT<- rowSums(colSums(exp_mut))
  subst_type_t<- table(as.character(subst_type3))
  mutProb<- prop.table(subst_type_t/n_triNT[names(subst_type_t)])

  # Calculate expected N/S per gene
  genes<- sort(unique(gene))
  NS<- NULL
  for(i in 1:length(genes)){
    g<- genes[i]
    exp_mut_gene<- exp_mut[g,,]
    S<- sum(exp_mut_gene[,"synonymous SNV"]*mutProb[rownames(exp_mut_gene)])
    N<- sum((rowSums(exp_mut_gene)-S)*mutProb[rownames(exp_mut_gene)])
    NS<- c(NS,N/S)    
  }
  names(NS)<- genes
  NS<- NS[!is.na(NS)]
  genes<- names(NS)
  
  # Calculate observed n & s per gene
  gene_var_t<- table(gene,variant)
  s<- gene_var_t[,"synonymous SNV"]   
  n<- rowSums(gene_var_t) - s  
  
  # Calculate dN/dS
  dNdS_matrix<- matrix(NA,length(genes),5,dimnames = list(genes,c("dNdS","p","cilo","ciup","q")))
  for(i in 1:length(genes)){
    g<- genes[i]
    fish_temp<- fisher.test(cbind(c(n[g],s[g]),c(round(1000000*NS[g]),1000000))) 
    dNdS_temp<- fish_temp$estimate
    dNdS_p_temp<- fish_temp$p.value
    dNdS_ci<- fish_temp$conf.int
    dNdS_matrix[g,c("dNdS","p","cilo","ciup")]<- c(dNdS_temp,dNdS_p_temp,dNdS_ci[1],dNdS_ci[2])
  }
  dNdS_matrix[,"q"]<- p.adjust(dNdS_matrix[,"p"],"fdr") 
  
  # Return matrix
  return(dNdS_matrix)
}

